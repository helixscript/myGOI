library(dplyr)
library(readr)
library(GenomicRanges)
library(parallel)
library(data.table)
library(ggplot2)
library(ggrepel)
library(grDevices)
library(RColorBrewer)
library(kableExtra)
library(UpSetR)
source('lib.R')

# Read in the congifuration file.
config <- yaml::read_yaml('myGOI.yml')

# Read intSite data and limit to sites near genes.
d <- read_tsv(config$inputDataPath, progress = FALSE, col_types = cols())

# Limit sites to study those within config$maxDistNearestGene nt of a gene boundary
d <- subset(d, abs(nearestFeatureDist) <= config$maxDistNearestGene)

# Since we are studying genes, expand instances where a site intersects more than one gene.
# "chr7+543431  ABC1,XYZ1" will be expanded to "chr7+543431  ABC1" and "chr7+543431  XYZ1"
d$nearestFeature <- toupper(gsub('\\s', '', d$nearestFeature))
d <- mutate(d, nearestFeature = strsplit(nearestFeature, ",")) %>% tidyr::unnest(nearestFeature)


# Limit sites to those from samples with at least config$minSampleAbund estimated cells.
sampleTotalAbundance <- group_by(d, GTSP) %>%
                        summarise(totalAbund = sum(estAbund)) %>%
                        ungroup() %>%
                        filter(totalAbund >= config$minSampleAbund)

d <- subset(d, GTSP %in% sampleTotalAbundance$GTSP)


# Limit sites to those with nearest genes seen in at least config$minGeneSubjects patients.
genePatients <- group_by(d, nearestFeature) %>%
                summarise(nPatients = n_distinct(patient)) %>%
                ungroup() %>%
                filter(nPatients >= config$minGeneSubjects)

d <- subset(d, nearestFeature %in% genePatients$nearestFeature)


# Create a gene-centric table of maximal abundance and relative abundance values.
# This table will be joined to the main data table (d).

geneAbundances <- group_by(subset(d, d$timePointDays > config$earlyVsLateCutoffDays), nearestFeature) %>%
                  summarise(maxAbund = max(estAbund),
                            maxRelAbund = max(relAbund)) %>%
                  ungroup()

# Load external data files.
HGNC     <- loadHGNCfull()
HGNCprev <- loadHGNCprev()


# Use the retired HGNC gene ids to update gene symbols to current symbols.
d$nearestFeature <- updateGeneSymbols(d$nearestFeature)

# Separate the data into early and late data sets.
a <- d[d$timePointDays <= config$earlyVsLateCutoffDays,]
b <- d[d$timePointDays >  config$earlyVsLateCutoffDays,]

# Count unique sites in each time point grouping.
a_totalSites <- n_distinct(a$posid)
b_totalSites <- n_distinct(b$posid)

genesToSubjects <- group_by(d, nearestFeature) %>% 
                   summarise(subjects = n_distinct(patient),
                             totalSites = n_distinct(posid)) %>% 
                   ungroup()


# Build a table describing genes in later timepoints.
lateGeneData <- subset(b, timePointDays >= config$longitudinal_minTimeDays) %>%
                group_by(nearestFeature) %>%
                  summarise(latestTimePointDays = sort(unique(timePointDays), decreasing = TRUE)[1],
                  longitudinalSubjects = n_distinct(patient),
                  longitudinalSites = n_distinct(posid),
                  longitudinalTimePoints = n_distinct(timePointDays)) %>%
                ungroup()

# For early and late data groups (a and b), create tables of site counts.
ag <- group_by(a, nearestFeature) %>% summarise(nSites = n_distinct(posid)) %>% ungroup()
bg <- group_by(b, nearestFeature) %>% summarise(nSites = n_distinct(posid)) %>% ungroup()            
              
# Update table names to reflect their source since they will be jointed together.
names(ag) <- paste0(names(ag), '_a')
names(bg) <- paste0(names(bg), '_b')              

# Count TOTAL sites in early and late groupings since they will be constants used
# in the following operation.
a_totalSites <- n_distinct(a$posid)
b_totalSites <- n_distinct(b$posid)


# Join the early and late tables together by gene symbols.
# Then, for each gene, determine change in integration frequency and the significance
# of that change using Fisher's Exact tests. After each gene test is completed,
# join pre-built tables to the result.

k <- left_join(ag, bg, by = c('nearestFeature_a' = 'nearestFeature_b')) %>% 
     tidyr::drop_na() %>%
     group_by(nearestFeature_a) %>%
       mutate(gene = nearestFeature_a,
              pVal = fisher.test(matrix(c(a_totalSites - nSites_a, 
                                          b_totalSites - nSites_b, 
                                          nSites_a , 
                                          nSites_b), 
                                        ncol = 2, byrow = FALSE))$p.value,
              earlyCount = nSites_a,
              lateCount  = nSites_b,
              earlyFreq  = nSites_a / a_totalSites,
              lateFreq   = nSites_b / b_totalSites,
              percentChange = ((lateFreq - earlyFreq) / earlyFreq) * 100) %>%
     ungroup() %>% 
     select(-nearestFeature_a) %>%
     left_join(genesToSubjects, by = c('gene' = 'nearestFeature')) %>%
     left_join(geneAbundances, by = c('gene' = 'nearestFeature')) %>%
     left_join(lateGeneData, by = c('gene' = 'nearestFeature')) %>% 
     select(-nSites_a, -nSites_b) %>%
     as.data.table()
 

# We will study the result both with and without correction for multiple comparisons. 
# Here we save the uncorrected pValues (raw), correct for multiple comparisons (adj)
# Then set pVal to NA. pVal will be toggled latter between raw and adj values.
k$pVal.raw <- k$pVal
k$pVal.adj <- p.adjust(k$pVal, method = 'BH')
k$pVal     <- NA


# Read in oncogene lists.
oncoGeneLists <- list('cosmic'     = unique(updateGeneSymbols(unique(readLines(config$COSMIC_oncogene_table)))),
                      'cosmic_tsg' = unique(updateGeneSymbols(unique(readLines(config$COSMIC_tsg_table)))),
                      'allOnco'    = unique(updateGeneSymbols(unique(readLines(config$allOnco_oncogene_table)))))

# The abundance metric includes genes whoes abundances are within the top 1% of all abundances.
# This scalar stores the lowest abundance in the top 1% and is set by setCategories().
currentAbundantCategoryThreshold <- NA


# Test if gene symbols are in the included oncogene lists.
k$oncoGene <- ifelse(k$gene %in% unique(c(oncoGeneLists[['allOnco']], oncoGeneLists[['cosmic']])), 'yes', '')


# First, set the working pValues to the BH corrected pValues.
k$pVal <- k$pVal.adj      # Set working pValue
k <- setCategories(k)     # Set categories based on working pValue
pValAdjVolcanoPlot <- volcanoPlot(k)      # Build volvano plot.
correctedGeneListTests <- geneListTests() # Run a series of enrichment tests on 3rd party gene lists.


# Second, set the working pValues to the uncorrected pValues.
k$pVal <- k$pVal.raw
k <- setCategories(k)
pValVolcanoPlot <- volcanoPlot(k)
uncorrectedGeneListTests <- geneListTests()

geneListCompTab <- bind_rows(lapply(c('Depleated', 'Enriched', 'Abundant', 'Longitudinal'), function(criteria){
  crit <- substr(criteria, 1, 1)
  bind_cols(tibble('Criteria' = criteria), 
            tibble('Genes' = n_distinct(k[grepl(crit, k$flag),]$gene)),
            bind_cols(lapply(c('allOnco', 'cosmic', 'cosmic_tsg'), function(test){
              a <- subset(uncorrectedGeneListTests, Class == crit)
              b <- subset(correctedGeneListTests, Class == crit)
              o <- tibble(x = paste0(ifelse(a[paste0(test, '__pval')] <= 0.05, '*', '-'), '/', 
                                     ifelse(b[paste0(test, '__pval')] <= 0.05, '*', '-'), ' ', 
                                sprintf("%.1f%%", a[paste0(test, '__percentGenes')] *100)))
             names(o) <- test
             o
  })))     
}))

upSetPlot <- buildUpSetPlot()

rmarkdown::render('myGOI.Rmd',
                  output_file = 'myGOI.pdf',
                  params = list('date'  = format(Sys.Date(), format="%B %d, %Y"),
                                'title' = 'Gene of Interest Report (GOI)',
                                'author' = 'CART research group'))


# Checks...

# L_test <- bind_rows(lapply(k[grepl('L', k$flag),]$gene, function(x){
#             o <- dplyr::filter(subset(d, nearestFeature == x), timePointDays >= config$longitudinal_minTimeDays)
#             tibble(gene = x, numSubjects = n_distinct(o$patient), timePoints = n_distinct(o$timePointDays), sites = n_distinct(o$posid))
#           }))



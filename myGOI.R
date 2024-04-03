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

config <- yaml::read_yaml('myGOI.yml')

# Read intSite data and limit to sites near genes.
d <- read_tsv(config$inputDataPath, progress = FALSE, col_types = cols())

# Limit sites to study those within config$maxDistNearestGene nt of a gene boundary
d <- subset(d, abs(nearestFeatureDist) <= config$maxDistNearestGene)

# Since we are studying genes, expand instances where a site intersects more than one gene.
d$nearestFeature <- toupper(gsub('\\s', '', d$nearestFeature))
d <- mutate(d, nearestFeature = strsplit(nearestFeature, ",")) %>% tidyr::unnest(nearestFeature)


# Limit sites to study to those from samples with at least config$minSampleAbund estimated cells.
sampleTotalAbundance <- group_by(d, GTSP) %>%
                        summarise(totalAbund = sum(estAbund)) %>%
                        ungroup() %>%
                        filter(totalAbund >= config$minSampleAbund)

d <- subset(d, GTSP %in% sampleTotalAbundance$GTSP)

# Limit sites to study to those with nearest genes seen in at least config$minGeneSubjects patients.
genePatients <- group_by(d, nearestFeature) %>%
                summarise(nPatients = n_distinct(patient)) %>%
                ungroup() %>%
                filter(nPatients >= config$minGeneSubjects)

d <- subset(d, nearestFeature %in% genePatients$nearestFeature)

geneAbundances <- group_by(subset(d, d$timePointDays > config$earlyVsLateCutoffDays), nearestFeature) %>%
                  summarise(maxAbund = max(estAbund),
                            maxRelAbund = max(relAbund)) %>%
                  ungroup()

HGNC     <- loadHGNCfull()
HGNCprev <- loadHGNCprev()

d$nearestFeature <- updateGeneSymbols(d$nearestFeature)

# Separate the data into early and late data sets.
a <- d[d$timePointDays <= config$earlyVsLateCutoffDays,]
b <- d[d$timePointDays >  config$earlyVsLateCutoffDays,]

a_totalSites <- n_distinct(a$posid)
b_totalSites <- n_distinct(b$posid)

genesToSubjects <- group_by(d, nearestFeature) %>% 
                   summarise(subjects = n_distinct(patient),
                             totalSites = n_distinct(posid)) %>% 
                   ungroup()

earlyGenesToSubjects <- group_by(a, nearestFeature) %>% 
                        summarise(earlySubjects = n_distinct(patient)) %>% 
                        ungroup()

lateGeneData <- group_by(subset(b, timePointDays >= config$longitudinal_minTimeDays), nearestFeature) %>%
                summarise(latestTimePointDays = sort(unique(timePointDays), decreasing = TRUE)[1],
                longitudinalSubjects = n_distinct(patient),
                longitudinalSites = n_distinct(posid),
                longitudinalTimePoints = n_distinct(timePointDays)) %>%
                ungroup()

ag <- group_by(a, nearestFeature) %>% summarise(nSites = n_distinct(posid)) %>% ungroup()
bg <- group_by(b, nearestFeature) %>% summarise(nSites = n_distinct(posid)) %>% ungroup()            
              
names(ag) <- paste0(names(ag), '_a')
names(bg) <- paste0(names(bg), '_b')              

a_totalSites <- n_distinct(a$posid)
b_totalSites <- n_distinct(b$posid)

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
     left_join(earlyGenesToSubjects, by = c('gene' = 'nearestFeature')) %>%
     left_join(geneAbundances, by = c('gene' = 'nearestFeature')) %>%
     left_join(lateGeneData, by = c('gene' = 'nearestFeature')) %>% 
     select(-nSites_a, -nSites_b) %>%
     as.data.table()
 
k$pVal.raw <- k$pVal
k$pVal.adj <- p.adjust(k$pVal, method = 'BH')
k$pVal     <- NA

oncoGeneLists <- list('cosmic'     = unique(updateGeneSymbols(unique(readLines(config$COSMIC_oncogene_table)))),
                      'cosmic_tsg' = unique(updateGeneSymbols(unique(readLines(config$COSMIC_tsg_table)))),
                      'allOnco'    = unique(updateGeneSymbols(unique(readLines(config$allOnco_oncogene_table)))))

currentAbundantCategoryThreshold <- NA

k$oncoGene <- ifelse(k$gene %in% unique(c(oncoGeneLists[['allOnco']], oncoGeneLists[['cosmic']])), 'yes', '')

# First, set the working pValues to the BH corrected pValues.
k$pVal <- k$pVal.adj
k <- setCategories(k)
pValAdjVolcanoPlot <- volcanoPlot(k)
correctedGeneListTests <- geneListTests()


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



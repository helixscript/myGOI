library(dplyr)
library(readr)

HGNC <- readr::read_tsv('HGNC_symbols.tsv', col_types = cols())
HGNC <- select(HGNC, `Approved symbol`, `Previous symbols`)
names(HGNC) <- c('symbol', 'prev')
HGNC$symbol <- toupper(HGNC$symbol)
HGNC$prev <- toupper(gsub('\\s', '', HGNC$prev))
HGNC <- distinct(HGNC)
HGNC <- mutate(HGNC, prev = strsplit(prev, ",")) %>%
        tidyr::unnest(prev) %>%
        tidyr::drop_na() %>%
        distinct()

o <- read_tsv('COSMIC_oncoGenes_20240104_Tier_1.tsv')

cosmic <- tibble(gene = toupper(o$`Gene Symbol`), tsg = grepl('TSG', o$`Role in Cancer`))
cosmic <- left_join(cosmic, HGNC, by = c('gene' = 'prev'))
cosmic$gene <- ifelse(is.na(cosmic$symbol), cosmic$gene, cosmic$symbol)

write(cosmic$gene, 'COSMIC_oncogenes.txt', ncolumns = 1)
write(subset(cosmic, tsg == TRUE)$gene, 'COSMIC_oncogenes_tumor_suppressors.txt', ncolumns = 1)

o <- read_tsv('allOnco_June2021.tsv')
allOnco <- tibble(gene = unique(o$symbol))

allOnco <- bind_rows(lapply(split(allOnco, allOnco$gene), function(x){
          x$geneUpdated <- FALSE
          
          if(x$gene[1] %in% HGNC$prev & ! x$gene[1] %in% HGNC$symbol){
            o <- subset(HGNC, prev == x$gene[1])
    
            if(nrow(o) == 1){
              x$gene <- o$symbol
              x$geneUpdated <- TRUE
            }
          }
          x
         }))

allOnco <- allOnco[! grepl('^\\d', allOnco$gene),]
write(allOnco$gene, 'allOnco.txt', ncolumns = 1)



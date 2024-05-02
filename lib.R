
ppNum <- function (n) format(n, big.mark = ",", scientific = FALSE, trim = TRUE)

loadHGNCprev <- function(){
  HGNC <- readr::read_tsv('data/HGNC_symbols.tsv', col_types = cols())
  
  HGNC <- select(HGNC, `Approved symbol`, `Previous symbols`)
  names(HGNC) <- c('symbol', 'prev')
  
  HGNC$symbol <- toupper(HGNC$symbol)
  
  HGNC$prev <- toupper(gsub('\\s', '', HGNC$prev))
  
  HGNC <- distinct(HGNC)
  
  # Create a table that only has alternative symbols.
  mutate(HGNC, prev = strsplit(prev, ",")) %>%
  tidyr::unnest(prev) %>%
  tidyr::drop_na() %>%
  distinct() %>%
  as.data.table()
}


loadHGNCfull <- function(){
  HGNC <- readr::read_tsv('data/HGNC_symbols.tsv', col_types = cols())
  
  HGNC <- select(HGNC, `Approved symbol`, `Previous symbols`, Status, `Alias symbols`)
  names(HGNC) <- c('symbol', 'prevSymbols', 'status', 'aliasSymbols')
  
  HGNC$symbol <- toupper(HGNC$symbol)
  HGNC$prevSymbols <- toupper(gsub('\\s', '', HGNC$prevSymbols))
  HGNC$aliasSymbols <- toupper(gsub('\\s', '', HGNC$aliasSymbols))

  as.data.table(distinct(HGNC))
}


buildUpSetPlot <- function(){
  z <- bind_rows(lapply(split(k, 1:nrow(k)), function(x){
         if(x$flag == 'none') return(tibble())
         o <- tibble('gene' = x$gene, 'Depleated' = 0, 'Enriched' = 0, 'Abundant' = 0, 'Longintudinal' = 0)
         if(grepl('D', x$flag)) o$Depleated <- 1
         if(grepl('E', x$flag)) o$Enriched <- 1
         if(grepl('A', x$flag)) o$Abundant <- 1
         if(grepl('L', x$flag)) o$Longintudinal <- 1
         o
  })) %>% as.data.frame()
  rownames(z) <- z$gene
  z$gene <- NULL
        
  upset(z, order.by="freq", text.scale = 1.2, point.size = 3)
}


geneListTests <- function(){
  bind_rows(lapply(c('D', 'E', 'A', 'L'), function(x){
    a <- k[grepl(x, k$flag),]
    b <- k[! grepl(x, k$flag),]
    
    bind_cols(tibble('Class' = x), mapply(function(genes, label){
      o <- tibble(
                  enriched = (sum(a$gene %in% genes) / n_distinct(a$gene)) > (sum(b$gene %in% genes) / n_distinct(b$gene)),
                  percentGenes = sum(a$gene %in% genes) / n_distinct(a$gene),
                  pval = fisher.test(matrix(c(sum(a$gene %in% genes),
                                              sum(! a$gene %in% genes),
                                              sum(b$gene %in% genes),
                                              sum(! b$gene %in% genes)), ncol = 2, byrow = FALSE))$p.val)
                  
      names(o) <- c(paste0(label, '__enriched'), paste0(label, '__percentGenes'), paste0(label, '__pval'))
      o
      
    }, oncoGeneLists, names(oncoGeneLists), SIMPLIFY = FALSE))
  }))
}

logValue <- function(x) log2(abs(x)) * ifelse(x < 0, -1, 1)


updateGeneSymbols <- function(g){
  g <- toupper(gsub('\\s', '', g))
  g2 <- g
  for(gene in unique(g)){
    if(gene %in% HGNCprev$prev & ! gene %in% HGNCprev$symbol){
      o <- HGNCprev[prev == gene]
      if(nrow(o) == 1){
        message(gene, '-> ', o$symbol)
        g2[g2 == gene] <- o$symbol
      }
    }
  }
  
  message(sprintf("%.2f%%", (sum(g != g2) / length(g))*100), ' gene symbols updated.')
  g2
}

setCategories <- function(k){
  depletedGenes <- subset(k, percentChange <= 0 & pVal <= 0.05)
  enrichedGenes <- subset(k, percentChange > 0 & pVal <= 0.05)
  
  currentAbundantCategoryThreshold <<- arrange(k, desc(maxAbund))[ceiling(nrow(k)*0.01),]$maxAbund
  
  abundantGenes <- subset(k, maxAbund >= currentAbundantCategoryThreshold)

  k$flag <- ''
  k[k$percentChange <= 0 & 
    k$gene %in% depletedGenes$gene,]$flag <- paste0(k[k$percentChange <= 0 & 
                                                        k$gene %in% depletedGenes$gene,]$flag, 'D')
  
  k[k$percentChange > 0 & 
    k$gene %in% enrichedGenes$gene,]$flag <- paste0(k[k$percentChange > 0 & 
                                                      k$gene %in% enrichedGenes$gene,]$flag, 'E')
  
  k[k$gene %in% abundantGenes$gene,]$flag <- paste0(k[k$gene %in% abundantGenes$gene,]$flag, 'A')
  
  k[k$longitudinalSubjects >= config$longitudinal_minNumSubjects & 
    k$longitudinalSites >= config$longitudinal_minNumSites & 
 ###   k$longitudinalTimePoints >= config$longitudinal_minNumTimepoints &
    ! grepl('D', k$flag),]$flag <- paste0(k[k$longitudinalSubjects >= config$longitudinal_minNumSubjects & 
                                            k$longitudinalSites >= config$longitudinal_minNumSites & 
                                  ###          k$longitudinalTimePoints >= config$longitudinal_minNumTimepoints &
                                            ! grepl('D', k$flag),]$flag, 'L')
  k[k$flag == '',]$flag <- 'none'
  k
}

volcanoPlot <- function(k){
  k$volcanoPlot_x <- logValue(k$percentChange)
  k$volcanoPlot_y <- log(1 / k$pVal)
  
  k$geneLabel <- k$gene
  o <- split(k, k$percentChange >= 0)
  
  o[[1]] <- arrange(o[[1]], desc(volcanoPlot_y))
  o[[1]][(config$volcanoPlot_numTopGeneLabels+1):nrow(o[[1]]),]$geneLabel <- 'none'
  
  o[[2]] <- arrange(o[[2]], desc(volcanoPlot_y))
  o[[2]][(config$volcanoPlot_numTopGeneLabels+1):nrow(o[[2]]),]$geneLabel <- 'none'
  
  k <- bind_rows(o)
  
  # Update
  k[k$flag == 'EAL',]$geneLabel <- k[k$flag == 'EAL',]$gene
  
  flagLevels <- unique(k$flag)
  
  flagLevels <- flagLevels[flagLevels != 'none']
  flagLevels <- c('none', flagLevels[order(nchar(flagLevels))])
  
  k$flag <- factor(k$flag, levels = flagLevels)
  
  colors <- c('gray50', grDevices::colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(n_distinct(k$flag) - 1))
  
  jitter_pos <- position_jitter(width= 0.35, height = 0.35, seed = 1)
  
  k <- bind_rows(list(subset(k, flag == 'none'), subset(k, flag != 'none')))
  
  k$geneLabel <- ifelse(k$geneLabel == 'none', '', k$geneLabel)
  
  # Force label for current report -- move forced gene labels to config file.
  k[k$gene == 'FOXP1',]$geneLabel <- 'FOXP1'
  
  p <- volcanoPlot <- ggplot(subset(k, abs(percentChange) >= 2), aes(volcanoPlot_x, volcanoPlot_y, label = geneLabel, fill = flag)) + 
       scale_fill_manual(name = 'Class', values = colors, labels = flagLevels, drop = FALSE) +
       scale_x_continuous(limits = c(min(k$volcanoPlot_x - 3), max(k$volcanoPlot_x))) +
       scale_y_continuous(limits = c(-1, max(k$volcanoPlot_y+3))) +
       geom_hline(yintercept = log(1/0.05), color = 'black', linetype = 'dashed') +
       geom_jitter(position = jitter_pos, size = 2, shape = 21, alpha = 0.9) +  
       geom_text_repel(position = jitter_pos, size = 2, point.size = 2, direction = 'both', seed = 1,
                       max.overlaps = Inf,  max.time = 10, max.iter = 50000) +
       labs(x = 'log2(percent gene integration frequency change)', y = 'log(1/p-value)') +
       theme(legend.position = "bottom",
             legend.key=element_blank(),
             text = element_text(size = 12),
             panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
             panel.background = element_blank(), axis.line = element_line(colour = "black")) +
       guides(fill = guide_legend(override.aes = list(size=4)))
  
  p
}

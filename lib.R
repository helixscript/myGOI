logValue <- function(x) log2(abs(x)) * ifelse(x < 0, -1, 1)

setCategories <- function(k){
  depletedGenes <- subset(k, percentChange <= 0 & pVal <= 0.05)
  enrichedGenes <- subset(k, percentChange > 0 & pVal <= 0.05)
  abundantGenes <- subset(k, maxAbund >= arrange(k, desc(maxAbund))[ceiling(nrow(k)*0.01),]$maxAbund)

  k$flag <- ''
  k[k$percentChange <= 0 & k$gene %in% depletedGenes$gene,]$flag <- paste0(k[k$percentChange <= 0 & k$gene %in% depletedGenes$gene,]$flag, 'D')
  k[k$percentChange > 0  & k$gene %in% enrichedGenes$gene,]$flag <- paste0(k[k$percentChange > 0  & k$gene %in% enrichedGenes$gene,]$flag, 'E')
  k[k$gene %in% abundantGenes$gene,]$flag <- paste0(k[k$gene %in% abundantGenes$gene,]$flag, 'A')
  k[k$longitudinalSubjects >= 3 & k$longitudinalSites >= 3 & ! grepl('D', k$flag),]$flag <- paste0(k[k$longitudinalSubjects >= 3 & k$longitudinalSites >= 3 & ! grepl('D', k$flag),]$flag, 'L')

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
  flagLevels <- unique(k$flag)
  
  flagLevels <- flagLevels[flagLevels != 'none']
  flagLevels <- c('none', flagLevels[order(nchar(flagLevels))])
  
  k$flag <- factor(k$flag, levels = flagLevels)
  
  colors <- c('gray50', grDevices::colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(n_distinct(k$flag) - 1))
  
  jitter_pos <- position_jitter(width= 0.35, height = 0.35, seed = 1)
  
  k <- bind_rows(list(subset(k, flag == 'none'), subset(k, flag != 'none')))
  
  k$geneLabel <- ifelse(k$geneLabel == 'none', '', k$geneLabel)
  
  p <- volcanoPlot <- ggplot(subset(k, abs(percentChange) >= 2), aes(volcanoPlot_x, volcanoPlot_y, label = geneLabel, fill = flag)) + 
       scale_fill_manual(name = 'Class', values = colors, labels = flagLevels, drop = FALSE) +
       scale_x_continuous(limits = c(min(k$volcanoPlot_x - 3), max(k$volcanoPlot_x))) +
       scale_y_continuous(limits = c(-1, max(k$volcanoPlot_y+3))) +
       geom_hline(yintercept = log(1/0.05), color = 'black', linetype = 'dashed') +
       geom_jitter(position = jitter_pos, size = 2, shape = 21, alpha = 0.9) +  
       geom_text_repel(position = jitter_pos, size = 3, point.size = 6, direction = 'both', seed = 1,
                       max.overlaps = 15,  max.time = 5, max.iter = 5000) +
       labs(x = 'log2(percent gene integration frequency change)', y = 'log(1/p-value)') +
       theme(legend.position = "bottom",
             legend.key=element_blank(),
             text = element_text(size = 12),
             panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
             panel.background = element_blank(), axis.line = element_line(colour = "black")) +
       guides(fill = guide_legend(override.aes = list(size=4)))
  
  k$volcanoPlot_x <- NULL
  k$volcanoPlot_y <- NULL
  p
}
library(dplyr)
library(readr)
library(GenomicRanges)
library(parallel)
library(data.table)
library(ggrepel)

CPUs <- 24

earlyLateDaysCutoffs <- c(0, 14)
maxDistToNearestFeatures <- c(0, 25000, 50000)
volcanoPlot_numTopGeneLabels <- 10

cluster <- makeCluster(CPUs)

o <- rbindlist(lapply(earlyLateDaysCutoffs, function(daysCutOff){
       rbindlist(lapply(maxDistToNearestFeatures, function(maxDistToGene){
       message('Starting calc for early vs. late days cutoff: ', daysCutOff, ' and maxDistToGene: ', maxDistToGene, '.')
    
       # Read intSite data and limit to sites near genes.
       o <- read_tsv('expandedIntSiteData.tsv.gz', progress = FALSE, col_types = cols())
       o <- o[abs(o$nearestFeatureDist) <= maxDistToGene,]
    
       # Duplicate sites that interset more than one transcription unit.
       o <- data.table(tidyr::unnest(o, nearestFeature = strsplit(nearestFeature, ',')))
    
       # Separate the data into early and late data sets.
       a <- o[timePointDays <= daysCutOff]
       b <- o[timePointDays >  daysCutOff]
    
       a_totalSites <- n_distinct(a$posid)
       b_totalSites <- n_distinct(b$posid)
    
       # Limit data to sites near genes found in both early and late time point groups.
       i <- base::intersect(a$nearestFeature, b$nearestFeature)
       a <- a[nearestFeature %in% i]
       b <- b[nearestFeature %in% i]
       
       message('Exporting data frames to cluster.')
       clusterExport(cluster, c('a', 'b'), envir = environment())
       
       # Run Fisher's exact test on each nearest gene, early vs late.
       message('Starting Fishers exact tests for genes in the data set.')
       
       k <- rbindlist(parLapply(cluster, split(i, ntile(1:length(i), CPUs)), function(g){
              library(data.table)
      
              rbindlist(lapply(g, function(x){
                m <- matrix(c(length(unique(a[nearestFeature != x]$posid)),
                              length(unique(b[nearestFeature != x]$posid)),
                              length(unique(a[nearestFeature == x]$posid)),
                              length(unique(b[nearestFeature == x]$posid))), 
                            ncol = 2, byrow = FALSE)
        
               data.table(gene = x, earlyCount = m[1,2], lateCount = m[2,2], pVal = fisher.test(m)$p.value)            
             }))
           }))
       
       message('Fishers exact tests done.')
       
      # Correct pvalues for multiple comparisons.
      k$pVal.adj <- p.adjust(k$pVal, method = 'BH')
      
      k$daysCutOff <- daysCutOff
      k$maxDistToGene <- maxDistToGene
      k
    }))
  }))




# Standardize site counts and calcualte percent change.

k2 <- bind_rows(lapply(split(k, 1:nrow(k)), function(x){
        a <- (x$earlyCount / a_totalSites)
        b <- (x$lateCount  / b_totalSites)
        
        x$percentChange  <- ((b - a) / a) * 100
        x$plotScore      <- log(1 / x$pVal.adj)
        x
      }))



logValue <- function(x) log2(abs(x)) * ifelse(x < 0, -1, 1)

k2$geneLabel <- k2$gene
d <- split(k2, k2$percentChange >= 0)
d[[1]] <- arrange(d[[1]], desc(plotScore))
d[[1]][volcanoPlot_numTopGeneLabels:nrow(d[[1]]),]$geneLabel <- ''
d[[2]] <- arrange(d[[2]], desc(plotScore))
d[[2]][volcanoPlot_numTopGeneLabels:nrow(d[[2]]),]$geneLabel <- ''
d <- bind_rows(d)

jitter_pos <- position_jitter(width=.25, height = 0.25, seed = 1)

p <- ggplot(subset(d, abs(percentChange) >= 3), aes(logValue(percentChange), plotScore, label = geneLabel)) + 
     geom_hline(yintercept = log(1/0.05), color = 'red') +
     geom_jitter(position = jitter_pos, shape = 21, alpha = 0.8) +
     geom_text_repel(position = jitter_pos, size = 2.5, point.padding = 1) +
     scale_x_continuous(limits = c(-7.5, 13)) +
     labs(x = 'log2(percent gene integration frequency change)', y = 'log(1/p-value)') +
     theme(text = element_text(size = 10),
           panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
           panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggsave(p, file = 'volcanoPlot.pdf', units = 'in', width = 4, height = 8, dpi = 300)
readr::write_tsv(subset(d, abs(percentChange) >= 3), 'volcanoPlot.tsv')


library(dplyr)
library(GenomicRanges)
library(parallel)
library(ggrepel)

CPUs <- 20
earlyLateDaysCutoff <- 14
maxDistToNearestFeature <- 0 # Set to zero for in TU only.

cluster <- makeCluster(CPUs)

# Read intSite data and limit to sites near genes.
o <- readr::read_tsv('expandedIntSiteData.tsv.gz')
o <- subset(o, abs(o$nearestFeatureDist) <= maxDistToNearestFeature)

# Duplicate sites that interset more than one transcription unit.
o <- tidyr::unnest(o, nearestFeature = strsplit(nearestFeature, ','))

# Separate the data into early and late data sets.
a <- subset(o, timePointDays <  earlyLateDaysCutoff)
b <- subset(o, timePointDays >= earlyLateDaysCutoff)

a_totalSites <- n_distinct(a$posid)
b_totalSites <- n_distinct(b$posid)

# Limit data to sites near genes found in both early and late time point groups.
i <- intersect(a$nearestFeature, b$nearestFeature)
a <- subset(a, nearestFeature %in% i)
b <- subset(b, nearestFeature %in% i)
clusterExport(cluster, c('a', 'b'))


# Run Fisher's exact test on each nearest gene, early vs late.
k <- bind_rows(parLapply(cluster, split(i, ntile(1:length(i), CPUs)), function(g){
       library(dplyr)
       bind_rows(lapply(g, function(x){
         m <- matrix(c(n_distinct(subset(a, nearestFeature != x)$posid),
                       n_distinct(subset(b, nearestFeature != x)$posid),
                       n_distinct(subset(a, nearestFeature == x)$posid),
                       n_distinct(subset(b, nearestFeature == x)$posid)), 
                       ncol = 2, byrow = FALSE)
   
         tibble(gene = x, earlyCount = m[1,2], lateCount = m[2,2], pVal = fisher.test(m)$p.value)            
       }))
    }))


# Correct pvalues for multiple comparisons.
k$pVal.adj <- p.adjust(k$pVal, method = 'BH')


# Standardize site counts and calcualte percent change.
k2 <- bind_rows(lapply(split(k, 1:nrow(k)), function(x){
        # if(x$pVal.adj < 0.01) browser()
        a <- (x$earlyCount / a_totalSites)
        b <- (x$lateCount  / b_totalSites)
        
        x$percentChange  <- ((b - a) / a) * 100
        x$logChange      <- -log2(b/a)
        x$plotScore      <- log(1 / x$pVal.adj)
        x
      }))

logValue <- function(x) log2(abs(x)) * ifelse(x < 0, -1, 1)

k2$geneLabel <- k2$gene
d <- split(k2, k2$percentChange >= 0)
d[[1]] <- arrange(d[[1]], desc(plotScore))
d[[1]][9:nrow(d[[1]]),]$geneLabel <- ''
d[[2]] <- arrange(d[[2]], desc(plotScore))
d[[2]][9:nrow(d[[2]]),]$geneLabel <- ''
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

save.image('myImage.RData')

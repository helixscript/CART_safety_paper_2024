options(stringsAsFactors = FALSE, useFancyQuotes = FALSE)
library(gt23)    # https://github.com/everettJK/gt23
library(RMySQL)  # loads DBI
library(GenomicRanges)
library(dplyr)
library(parallel)

# Minimum fragment width allowed.
minRangeWidth <- 10

# Read in external sample id file
sampleIDs <- readLines('ZetaSampleList')

# Retrieve all sample data from the specimen database.
dbConn  <- dbConnect(MySQL(), group = 'specimen_management')
sampleData <- dbGetQuery(dbConn, 'select * from gtsp')
sampleData <- subset(sampleData, SpecimenAccNum %in% sampleIDs)
dbDisconnect(dbConn)

# Retrieve all integration site data from the production databse.
dbConn  <- dbConnect(MySQL(), group = 'intsites_miseq')
analyzedSamples <- dbGetQuery(dbConn, 'select * from samples')
analyzedSamples <- sub('\\-\\d+', '', analyzedSamples$sampleName)
dbDisconnect(dbConn)

# Retrieve fragment data from db and limit fragments to hg38 reference genome.
intSites <- gt23::getDBgenomicFragments(unique(sampleIDs), 'specimen_management', 'intsites_miseq')
intSites <- subset(intSites, refGenome == 'hg38')

# Remove known NTNG2 artifact.
a <- GenomicRanges::makeGRangesFromDataFrame(data.frame(seqnames = 'chr9', strand = '*', start = 132162058, end = 132244526))
o <- findOverlaps(a, intSites, ignore.strand = TRUE)
intSites <- intSites[-subjectHits(o)]

# Remove very short ranges because they are likely not real and may break downstream fragment standardization.
intSites <- intSites[width(intSites) >= minRangeWidth]

# Remove negative timepoints.
intSites <- intSites[! grepl('\\-', intSites$timePoint),]

# correct odd timepoints


intSites[intSites$timePoint == 'W2']$timePoint  <- 'D14'
intSites[intSites$timePoint == 'D14']$timePointDays   <- 14
intSites[intSites$timePoint == 'D14']$timePointMonths <- 0

intSites[intSites$timePoint == 'W8']$timePoint  <- 'M2'
intSites[intSites$timePoint == 'M2']$timePointDays   <- 60
intSites[intSites$timePoint == 'M2']$timePointMonths <- 2

# Limit sample data to samples with one more integration sites.
sampleData <- subset(sampleData, SpecimenAccNum %in% intSites$GTSP)


# Build a virtual CPU cluster for parallelization.
cluster <- makeCluster(8)
clusterExport(cluster, 'intSites')


# Standardize fragments, call sites, and annotate sites.
processingList <- split(sampleData, paste(sampleData$Trial, sampleData$Patient))

o <- unlist(GRangesList(parLapply(cluster, processingList, function(x){
       options(stringsAsFactors = FALSE)
       library(gt23)  
       library(GenomicRanges)
       library(dplyr)
  
       a <- gt23::stdIntSiteFragments(subset(intSites, GTSP %in% x$SpecimenAccNum)) %>% 
            gt23::collapseReplicatesCalcAbunds() %>%
            gt23::annotateIntSites(CPUs = 3)
     
       a$trial <- x$Trial[1]
       
       a
     })))


# Test control sample from database.
if(! subset(o, GTSP == 'GTSP9999')$estAbund == 4){
  stop('Error -- control sample failed to produce the expected number of fragments.')
} else {
  o <- subset(o, GTSP != 'GTSP9999')
}


# Build a sample level summary table
summary <- group_by(data.frame(o), trial, patient, GTSP, timePoint, cellType) %>%
           summarise(num_sites = n_distinct(posid),
                     num_reads = sum(reads),
                     total_inferred_cells = sum(estAbund),
                     max_inffered_cells = max(estAbund),
                     max_relative_abundance = max(relAbund),
                     nearest_gene = paste0(nearestFeature[relAbund == max(relAbund)], collapse = '|')) %>%
           ungroup()

openxlsx::write.xlsx(summary, 'zetaIntSiteSampleSummary.xlsx')

readr::write_tsv(data.frame(o), 'zetaIntSiteData.tsv.gz')

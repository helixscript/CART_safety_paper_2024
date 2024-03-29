options(stringsAsFactors = FALSE, useFancyQuotes = FALSE)
library(gt23)    # https://github.com/everettJK/gt23
library(RMySQL)  # loads DBI
library(GenomicRanges)
library(dplyr)
library(parallel)

# Minimum fragment width allowed.
minRangeWidth <- 10

# Read in external sample id file
sampleIDs <- readLines('masterSampleList')

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

# Remove very short ranges because they are likely not real and may break downstream fragment standardization.
intSites <- intSites[width(intSites) >= minRangeWidth]

# Remove negative timepoints.
intSites <- intSites[! grepl('\\-', intSites$timePoint),]

# correct odd timepoints
intSites[intSites$timePoint == '10']$timePoint   <- 'D10'
intSites[intSites$timePoint == '14']$timePoint   <- 'D14'
intSites[intSites$timePoint == '28']$timePoint   <- 'D28'
intSites[intSites$timePoint == 'D07']$timePoint  <- 'D7'
intSites[intSites$timePoint == 'W52']$timePoint  <- 'Y1'

intSites[intSites$timePoint == 'Y1']$timePointDays   <- 365
intSites[intSites$timePoint == 'Y1']$timePointMonths <- 12


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
       
       ### write(paste(x$Trial[1], x$Patient[1], 'done', date()), 'log', append = TRUE)
       
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

openxlsx::write.xlsx(summary, 'intSiteSampleSummary.xlsx')

readr::write_tsv(data.frame(o), 'intSiteData.tsv.gz')

# Expand current sample list to include ALL samples from the same patient.
dbConn  <- dbConnect(MySQL(), group = 'specimen_management')
sampleData <- dbGetQuery(dbConn, 'select * from gtsp')
dbDisconnect(dbConn)


# Create unique pairs of Trial and Patient ids rooted in the master sample list.
d <- distinct(select(subset(sampleData, SpecimenAccNum %in% sampleIDs), Trial, Patient))


# Find additional samples associated with these pairings not in the already compiled data set.
additional_ids <- unique(unlist(lapply(split(d, 1:nrow(d)), function(x){
                    ids <- subset(sampleData, Trial == x$Trial & Patient == x$Patient)$SpecimenAccNum
                    ids[! ids %in% o$GTSP]
                  })))


# Call and annotate expanded intSites.
intSites <- gt23::getDBgenomicFragments(additional_ids, 'specimen_management', 'intsites_miseq')
intSites <- subset(intSites, refGenome == 'hg38')
intSites <- intSites[width(intSites) >= minRangeWidth]

oe <- unlist(GRangesList(parLapply(cluster, split(intSites, intSites$patient), function(x){
         options(stringsAsFactors = FALSE)
         library(gt23)  
         library(GenomicRanges)
         library(dplyr)
  
         gt23::stdIntSiteFragments(x)  %>% 
         gt23::collapseReplicatesCalcAbunds() %>%
         gt23::annotateIntSites(CPUs = 3)
       })))

oe <- oe[! grepl('\\-', oe$timePoint),]


# Join trial data to expanded intSite data.
sampleData$trial <- sampleData$Trial
oe <-  makeGRangesFromDataFrame(left_join(data.frame(oe), select(sampleData, SpecimenAccNum, trial), by = c('GTSP' = 'SpecimenAccNum')), keep.extra.columns = TRUE)

o2 <- c(o, oe)

# Write out expanded intSite data.
readr::write_tsv(data.frame(o2), 'expandedIntSiteData.tsv.gz')

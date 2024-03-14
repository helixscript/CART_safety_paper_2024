library(gt23)
library(RMySQL)
library(tidyverse)
library(GenomicRanges)


numCores <- 4
if(!exists('TRIAL')) {TRIAL <- 'ALL'} # ALL, CLL, CALL= ALL+CLL
if(!exists('RESP')) {RESP <- 'na'} # whether filter by response class, RE= responders, nRE non responders, na= no filter



#import response data

library(readxl)
CARTSite_Response <- read_excel("CARTSite.Response.BOR.01142022.xlsx", 
                                col_types = c("text", "numeric", "skip", 
                                              "date", "skip", "skip", "skip", "date", 
                                              "numeric", "text", "numeric", "text"))
# BOR= best overal response
colnames(CARTSite_Response) <- c('Trial','ID','inf_date','response_date','response_timepoint','BOR','BORc','note')


#load intsites
if( !file.exists(here::here('all_intsites_20230113.rds'))) {
  num_cores <- 4
  trials <- c('UPENN_CART19_CLL', 'UPENN_CART19_ALL', 'UPENN_CART19_AML', 'Gill_CART19_18415')
  # Retrieve all relevant samples from sample database.
  invisible(sapply(dbListConnections(MySQL()), dbDisconnect))
  dbConn  <- dbConnect(MySQL(), group='specimen_management')
  arr <-paste0("select * from gtsp where Trial IN ('",paste(trials,collapse = "','"),"')")
  samples <- dbGetQuery(dbConn, arr)
  intSites <- getDBgenomicFragments(samples$SpecimenAccNum, 
                                    'specimen_management', 'intsites_miseq') %>%
    GenomicRanges::as.data.frame() %>% filter(refGenome == 'hg38') %>%
    makeGRangesFromDataFrame(keep.extra.columns=TRUE) %>%
    stdIntSiteFragments(CPUs = num_cores ) %>%
    collapseReplicatesCalcAbunds() %>%
    annotateIntSites(CPUs = num_cores)
  saveRDS(intSites, here::here('all_intsites_20230113.rds'))
}
intSites <- readRDS(here::here('all_intsites_20230113.rds'))



# match GTSP, PID and BOR
# responder = BOR 3-4
dbConn  <- dbConnect(MySQL(), group='specimen_management')
samples_all <- dbGetQuery(dbConn,'select * from gtsp')%>%
  filter(SpecimenAccNum %in% intSites$GTSP)

gtsp_to_pid <- samples_all %>% select(c("SpecimenAccNum","SamplePatientCode","Patient","Trial")) %>%
  mutate(SamplePatientCode = str_replace(SamplePatientCode,'^CHP','CHOP')) %>%
  separate(SamplePatientCode,c('Trial2',NA),'(?<=^(CHOP|UPCC))(?=\\d*(_|\\-))') %>%
  mutate(Trial2=ifelse(Trial2=='xxx','CHOP',Trial2)) %>%
  mutate(Patient=str_replace(Patient,'^\\D+','')) %>%
  mutate(Patient=str_replace(Patient,'\\-?[^0-9\\-].*','')) %>%
  mutate(PID=paste0(Trial2,Patient)) %>%
  mutate(GTSP=SpecimenAccNum) %>%
  select(c(GTSP,PID,Trial))

pid_to_bor <- CARTSite_Response %>%
  group_by(Trial,ID) %>%
  summarise(BORc=max(BORc), .groups = 'drop') %>%
  mutate(PID=paste0(Trial,'-',sprintf("%02d", ID))) %>%
  select(c(PID,BORc))

gtsp_to_bor <- left_join(gtsp_to_pid,pid_to_bor,by='PID') 



## transform the database to be compatible with chis' code ----------------
new_condensed <- left_join(intSites %>% as.data.frame(),gtsp_to_bor,by='GTSP') %>% 
  filter(seqnames %in% paste0('chr',c(1:22,'X','Y'))) %>% 
  dplyr::rename("specimen"=GTSP,
                'refgenome'=refGenome,
                'celltype'=cellType,
                'timepoint'=timePoint) %>% 
  mutate(timepoint=tolower(timepoint),
         in_gene=ifelse(nearestFeatureDist==0,nearestFeature,FALSE),
         in_geneOrt=ifelse(in_gene!='FALSE',nearestFeatureStrand,NA)) %>% 
  dplyr::rename('nearest_geneDist'=nearestFeatureDist,
                'nearest_gene'=nearestFeature,
                'nearest_geneOrt'= nearestFeatureStrand) %>% 
  mutate(gene_id_wo_annot=ifelse(in_gene == "FALSE",nearest_gene,in_gene)) #%>% 

new_condensed$gene_id_wo_annot <- sapply(
  strsplit(new_condensed$gene_id_wo_annot, ","), "[[", 1
)

new_condensed <- new_condensed %>% 
  mutate(gene_id=gene_id_wo_annot) %>% 
  mutate(gene_id=ifelse(in_gene=='FALSE',gene_id,paste0(gene_id,'*'))) %>% 
  mutate(gene_id=ifelse(abs(nearestOncoFeatureDist)<50000,paste0(gene_id,'~'),gene_id)) %>% 
  mutate(gene_id=ifelse(abs(nearestlymphomaFeatureDist)<50000,paste0(gene_id,'!'),gene_id)) %>% 
  mutate(patient =PID ) %>% 
  mutate(PID= NULL)


if(TRIAL=='ALL') {
  new_condensed <- new_condensed %>%  
    filter(Trial=="UPENN_CART19_ALL")
} else if (TRIAL=='CLL') {
  new_condensed <- new_condensed %>%  
    filter(Trial=="UPENN_CART19_CLL")
} else if (TRIAL=='CALL') {
  new_condensed <- new_condensed %>%  
    filter((Trial=="UPENN_CART19_CLL") | (Trial=="UPENN_CART19_ALL"))
} 

if(RESP=='RE') {
  new_condensed <- new_condensed %>%  
    filter(BORc>=3)
} else if (RESP=='nRE') {
  new_condensed <- new_condensed %>%  
    filter(BORc<=3)
}


new_cond_granges <- new_condensed %>% # removing patient 107 as they developed an additional unrelated cancer
  filter(patient!="CHOP959-107") %>% # add here ALL filter
  GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)

seqlevels(new_cond_granges) <- paste0("chr", c(1:22, "X", "Y", "M"))

genome_sequence <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
seqinfo(new_cond_granges) <- seqinfo(genome_sequence)[paste0("chr", c(1:22, "X", "Y", "M"))]

saveRDS(new_cond_granges,here::here('condensed_intsites',paste0(TRIAL,'_',RESP,'_conintsites.rds')))


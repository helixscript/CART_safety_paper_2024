library(gt23)
library(RMySQL)
library(tidyverse)
library(GenomicRanges)
dbConn  <- dbConnect(MySQL(), group='specimen_management')
ss <- c('UPENN_CART19_CLL','UPENN_CART19_ALL','Gill_CART19_18415')
arr <-paste0("SELECT * FROM gtsp WHERE Trial IN ('",paste(ss,collapse = "','"),"')")
samples <- dbGetQuery(dbConn,arr)

numCores <- 4

#import response data

library(readxl)
CARTSite_Response <- read_excel("CARTSite.Response.BOR.01142022.xlsx", 
                                col_types = c("text", "numeric", "skip", 
                                              "date", "skip", "skip", "skip", "date", 
                                              "numeric", "text", "numeric", "text"))
# BOR= best overal response
colnames(CARTSite_Response) <- c('Trial','ID','inf_date','response_date','response_timepoint','BOR','BORc','note')

if( !file.exists('ALL_intSites.rds') ) {
  intSites <- getDBgenomicFragments(samples$SpecimenAccNum, 'specimen_management', 'intsites_miseq') %>% 
    as.data.frame() %>% 
    filter(refGenome == 'hg38') %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE) %>% 
    stdIntSiteFragments(CPUs = numCores ) %>% 
    collapseReplicatesCalcAbunds() %>% 
    annotateIntSites(CPUs = numCores)
  saveRDS(intSites, 'ALL_intSites.rds')
} else {
  intSites <- readRDS('ALL_intSites.rds')
}

paste0('chr',c(1:22,'X','Y'))
condensed_intsites <- readRDS("~/data/CART19/CART19_from_git2/data/condensed_intsites.rds") %>%
  as.data.frame(row.names = NULL)

# match GTSP, PID and BOR
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

# ttt <- left_join(intSites %>% as.data.frame(),gtsp_to_bor,by='GTSP')



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

setdiff(colnames(condensed_intsites),colnames(new_condensed)) ## TODO migth still need to add relRank
head(condensed_intsites)
head(new_condensed)

new_cond_granges <- new_condensed %>% 
  GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)

seqlevels(new_cond_granges) <- paste0("chr", c(1:22, "X", "Y", "M"))

genome_sequence <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
seqinfo(new_cond_granges) <- seqinfo(genome_sequence)

saveRDS(new_cond_granges,'ALL_condensed_intsites.rds')

# condensed_intsites %>% 
#   filter(nearest_gene != gene_id_wo_annot) %>%
#   select(c(gene_id_wo_annot,in_gene,nearest_gene)) %>% head()
# new_condensed %>% filter(grepl(',', nearest_gene)) %>% select(c(gene_id_wo_annot,in_gene,nearest_gene))


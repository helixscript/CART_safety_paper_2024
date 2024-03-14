library(rlang)
library(tidyverse)
if(!exists('workingDir')) {workingDir <- getwd()}

# JKE
#--------------------
source('required.R')


# ALL, CLL, CALL= ALL+CLL
# whether filter by response class, RE= responders, nRE non responders, na= no filter


# generating a grid of all the cancer/response group to compute

# JKE 
# to_run <- expand_grid(x= c('ALL','CLL','CALL'),y= c('RE','na')) 
to_run <- tibble(x = 'ANY', y = 'na')

message('Starting 01_condensed_intsites.R')
tt<- to_run %>%   
purrr::pmap(~source(file.path(workingDir,'01_condensed_intsites.R'),local=env(TRIAL=.x,RESP=.y)))

message('Starting 02_gene_stats.R')
tt<- to_run %>% 
purrr::pmap(~source(file.path(workingDir,'02_gene_stats.R'),local=env(TRIAL=.x,RESP=.y)))

message('Starting 03_gene_impact.R')
tt<- to_run %>% 
  purrr::pmap(~source(file.path(workingDir,'03_gene_impact.R'),local=env(TRIAL=.x,RESP=.y)))

#get text for reports

get_group_title <- function(xx) {
  if(xx=='ALL') {
    retval<-'ALL'
  } else if (xx=='CLL'){
    retval <- 'CLL'
  } else if (xx=='CALL') {
    retval <- 'ALL or CLL'
  } else {
    retval <- '----'
  }
  return(retval)
}

get_resp_title <- function(xx) {
  if(xx=='RE') {
    retval<-'CR/PRtd'
  } else if (xx=='nRE'){
    retval <- 'PR/NR'
  } else if (xx=='na') {
    retval <- 'CR/PRtd & PR/NR'
  } else {
    retval <- '----'
  }
  return(retval)
}

message('Starting to build reports.')
# compile goi reports pdf
tt <- to_run %>% 
  purrr::pmap(~rmarkdown::render( 
  input       = '04_report.Rmd', 
  envir       = new.env(), 
  params      = list(pTRIAL=.x,
                     pRESP=.y,
                     tTRIAL=get_group_title(.x),
                     tRESP=get_resp_title(.y)),      
  output_file = file.path('goi_reports',paste0(.x,'_',.y,'_goi.pdf'))))

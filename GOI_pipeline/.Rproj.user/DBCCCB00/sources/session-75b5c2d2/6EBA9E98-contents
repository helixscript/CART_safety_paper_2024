
library(tidyverse)
library(GenomicRanges)

if(!exists('TRIAL')) {TRIAL <- 'ALL'} # ALL, CLL, CALL= ALL+CLL
if(!exists('RESP')) {RESP <- 'na'} # whether filter by response class, RE= responders, nRE non responders, na= no filter

# Compile Data for each gene -------------------------------------------------
refGenes <- readRDS(file.path('utils', "hg38.refSeq.rds"))

mod_refGenes <- unlist(GenomicRanges::reduce(split(
  refGenes, refGenes$name2), min.gapwidth = 100000L))

mod_refGenes$name <- names(mod_refGenes)


## Annotate Cluster FDR across gene set ======================================

cart19_clusters <- readRDS(file =file.path('condensed_intsites',paste0(TRIAL,'_',RESP,'_cart19_clusters.rds')))

mod_refXclusters <- suppressWarnings(
  findOverlaps(mod_refGenes, cart19_clusters))
mod_refGenes$Within_Cluster <- FALSE
mod_refGenes$Cluster_target.min <- as.numeric(NA)
mod_refGenes[queryHits(mod_refXclusters)]$Within_Cluster <- TRUE
mod_refGenes[queryHits(mod_refXclusters)]$Cluster_target.min <- 
  cart19_clusters$target.min[subjectHits(mod_refXclusters)]


## List of all genes =========================================================
gene_list <- mod_refGenes$name

## Format TP and TDN sites for analysis ======================================
cond_uniq_sites <- readRDS(file.path('condensed_intsites',paste0(TRIAL,'_',RESP,'_conintsites.rds')))
tdn_sites <- cond_uniq_sites[cond_uniq_sites$timepoint == "d0"]
timepoint_sites <- cond_uniq_sites[cond_uniq_sites$timepoint != "d0"]

format_sites <- function(sites){
  df <- GenomicRanges::as.data.frame(sites, row.names = NULL) %>%
    mutate(geneid = gene_id_wo_annot) %>%
    dplyr::select(
      seqnames, start, strand, patient, timepoint, 
      celltype, posid, estAbund, relAbund) %>%
    as.data.frame()
  ranges <- IRanges(start = df$start, width = rep(1, nrow(df)))
  gr <- GRanges(
    seqnames = df$seqnames,
    ranges = ranges,
    strand = df$strand,
    seqlengths = seqlengths(sites),
    seqinfo = seqinfo(sites)
  )
  mcols(gr) <- dplyr::select(df, -seqnames, -start, -strand)
  gintools::unique_granges(gr)
}

tdn_gr <- format_sites(tdn_sites)
tp_gr <- format_sites(timepoint_sites)


# Determine which integrations were within a transcriptional unit ------------
tdn_in_gene <- findOverlaps(
  mod_refGenes, tdn_gr, maxgap = 5000, ignore.strand = TRUE)

tp_in_gene <- findOverlaps(
  mod_refGenes, tp_gr, maxgap = 5000, ignore.strand = TRUE)



# Assemble stats for each transcriptional unit -------------------------------
## From all patients

get_loci_id <- function(gr){
  paste0(seqnames(gr), strand(gr), start(gr), ":", end(gr))
}

convert_time <- function(timecode){
  x <- tolower(timecode)
  scale <- ifelse(
    grepl("d", x), 1, 
    ifelse(grepl("m", x), 2, 
           ifelse(grepl("y", x), 3, 0)))
  timepoint <- as.numeric(stringr::str_extract(x, "[0-9.]+"))
  if(any(scale == 0)){
    stop("Not all timecodes in acceptable format.")
  }
  time_convert <- c("d" = 1, "m" = 30, "y" = 365)
  unname(timepoint * time_convert[scale])
}

CR_pats <- cond_uniq_sites %>%
  as.data.frame() %>%
  filter((BORc==3)|(BORc==4)) %>% 
  pull(var=patient) %>% 
  unique()

NR_pats <- cond_uniq_sites %>%
  as.data.frame() %>%
  filter((BORc==1)|(BORc==2)) %>% 
  pull(var=patient) %>% 
  unique()


stats_tdn <- as.data.frame(
  tdn_gr[subjectHits(tdn_in_gene)], row.names = NULL) %>%
  dplyr::mutate(
    loci = get_loci_id(mod_refGenes[queryHits(tdn_in_gene)]),
    gene_name = as.character(mod_refGenes$name[queryHits(tdn_in_gene)]),
    gene_ort = as.character(strand(mod_refGenes[queryHits(tdn_in_gene)])),
    strand = as.character(strand)) %>%
  dplyr::select(
    loci, gene_name, gene_ort, posid, strand, estAbund, relAbund, patient) %>%
  dplyr::group_by(loci, gene_name, gene_ort) %>%
  dplyr::summarise(
    "TDN_num_patients" = n_distinct(patient),
    "TDN_num_pats_CR" = n_distinct(patient[patient %in% CR_pats]),
    "TDN_num_pats_NR" = n_distinct(patient[patient %in% NR_pats]),
    "TDN_num_sites" = n_distinct(posid),
    "TDN_num_sites_CR" = n_distinct(posid[patient %in% CR_pats]),
    "TDN_num_sites_NR" = n_distinct(posid[patient %in% NR_pats]),
    "TDN_sum_abund" = sum(estAbund),
    "TDN_sum_abund_CR" = sum(estAbund[patient %in% CR_pats]),
    "TDN_sum_abund_NR" = sum(estAbund[patient %in% NR_pats]),
    "TDN_peak_abund" = max(estAbund),
    "TDN_peak_abund_CR" = max(c(estAbund[patient %in% CR_pats],0)),
    #      "TDN_peak_abund_CR" = ifelse(TDN_peak_abund_CR < 0, 0, TDN_peak_abund_CR),
    "TDN_peak_abund_NR" = max(c(estAbund[patient %in% NR_pats],0)),
    #      "TDN_peak_abund_NR" = ifelse(TDN_peak_abund_NR < 0, 0, TDN_peak_abund_NR),
    "TDN_peak_relAbund" = max(relAbund),
    "tdn_same" = sum(as.integer(strand == gene_ort)),
    "tdn_same_CR" = sum(as.integer(strand == gene_ort)[patient %in% CR_pats]),
    "tdn_same_NR" = sum(as.integer(strand == gene_ort)[patient %in% NR_pats]),
    "tdn_oppo" = sum(as.integer(strand != gene_ort)),
    "tdn_oppo_CR" = sum(as.integer(strand != gene_ort)[patient %in% CR_pats]),
    "tdn_oppo_NR" = sum(as.integer(strand != gene_ort)[patient %in% NR_pats]))

stats_tp <- as.data.frame(
  tp_gr[subjectHits(tp_in_gene)], row.names = NULL) %>%
  dplyr::mutate(
    loci = get_loci_id(mod_refGenes[queryHits(tp_in_gene)]),
    gene_name = as.character(mod_refGenes$name[queryHits(tp_in_gene)]),
    gene_ort = as.character(strand(mod_refGenes[queryHits(tp_in_gene)])),
    strand = as.character(strand),
    timepoint = convert_time(as.character(timepoint))) %>%
  dplyr::select(
    loci, gene_name, gene_ort, posid, strand, 
    estAbund, relAbund, patient, timepoint) %>%
  dplyr::group_by(loci, gene_name, gene_ort, patient, posid, strand) %>%
  dplyr::summarise(
    "long_count" = n_distinct(timepoint),
    "first_time" = min(timepoint),
    "last_time" = max(timepoint),
    "sum_abund" = sum(estAbund),
    "peak_abund" = max(estAbund),
    "peak_relAbund" = max(relAbund)) %>%
  dplyr::ungroup() %>% 
  dplyr::group_by(loci, gene_name, gene_ort) %>%
  dplyr::summarise(
    "TP_num_patients" = n_distinct(patient),
    "TP_num_pats_CR" = n_distinct(patient[patient %in% CR_pats]),
    "TP_num_pats_NR" = n_distinct(patient[patient %in% NR_pats]),
    "long_count" = max(long_count),
    "max_time" = max(last_time),
    "max_span" = max(last_time - first_time),
    "TP_num_sites" = n_distinct(posid),
    "TP_num_sites_CR" = n_distinct(posid[patient %in% CR_pats]),
    "TP_num_sites_NR" = n_distinct(posid[patient %in% NR_pats]),
    "TP_sum_abund" = sum(sum_abund),
    "TP_sum_abund_CR" = sum(sum_abund[patient %in% CR_pats]),
    "TP_sum_abund_NR" = sum(sum_abund[patient %in% NR_pats]),
    "TP_peak_abund" = max(peak_abund),
    "TP_peak_abund_CR" = max(c(0,peak_abund[patient %in% CR_pats])),
    #      "TP_peak_abund_CR" = ifelse(TP_peak_abund_CR < 0, 0, TP_peak_abund_CR),
    "TP_peak_abund_NR" = max(c(0,peak_abund[patient %in% NR_pats])),
    #      "TP_peak_abund_NR" = ifelse(TP_peak_abund_NR < 0, 0, TP_peak_abund_NR),
    "TP_peak_relAbund" = max(peak_relAbund),
    "abund_gini" = gintools::pop_calcs(peak_abund, calc = "gini"),
    "tp_same" = sum(as.integer(strand == gene_ort)),
    "tp_same_CR" = sum(as.integer(strand == gene_ort)[patient %in% CR_pats]),
    "tp_same_NR" = sum(as.integer(strand == gene_ort)[patient %in% NR_pats]),
    "tp_oppo" = sum(as.integer(strand != gene_ort)),
    "tp_oppo_CR" = sum(as.integer(strand != gene_ort)[patient %in% CR_pats]),
    "tp_oppo_NR" = sum(as.integer(strand != gene_ort)[patient %in% NR_pats]))

################ add section to get features for all times

tp_tdn_gr <- format_sites(cond_uniq_sites)
tp_tdn_in_gene <- findOverlaps(
  mod_refGenes, tp_tdn_gr, maxgap = 5000, ignore.strand = TRUE)

stats_tp_tdn <- as.data.frame(
  tp_tdn_gr[subjectHits(tp_tdn_in_gene)], row.names = NULL) %>%
  dplyr::mutate(
    loci = get_loci_id(mod_refGenes[queryHits(tp_tdn_in_gene)]),
    gene_name = as.character(mod_refGenes$name[queryHits(tp_tdn_in_gene)]),
    gene_ort = as.character(strand(mod_refGenes[queryHits(tp_tdn_in_gene)])),
    strand = as.character(strand),
    timepoint = convert_time(as.character(timepoint))) %>%
  dplyr::select(
    loci, gene_name, gene_ort, posid, strand, 
    estAbund, relAbund, patient, timepoint) %>% #TODO not sure it should be strand specific
  dplyr::group_by(loci, gene_name, gene_ort, patient, posid, strand) %>%
  dplyr::summarise(
    "TP_TDN_long_count" = n_distinct(timepoint),
    "TP_TDN_first_time" = min(timepoint),
    "TP_TDN_last_time" = max(timepoint),
    "TP_TDN_sum_abund" = sum(estAbund),
    "TP_TDN_peak_abund" = max(estAbund),
    "TP_TDN_peak_relAbund" = max(relAbund),
    "TP_TDN_span" = TP_TDN_last_time - TP_TDN_first_time ) %>%
  dplyr::ungroup() %>% 
  dplyr::group_by(loci, gene_name, gene_ort) %>%
  dplyr::summarise(
    "TP_TDN_num_patients" = n_distinct(patient),
    "TP_TDN_num_pats_CR" = n_distinct(patient[patient %in% CR_pats]),
    "TP_TDN_num_pats_NR" = n_distinct(patient[patient %in% NR_pats]),
    "TP_TDN_long_count" = max(TP_TDN_long_count),
    "TP_TDN_max_time" = max(TP_TDN_last_time),
    "TP_TDN_max_time_2" = ifelse(length(TP_TDN_last_time)>=2,
                             sort(TP_TDN_last_time, decreasing=TRUE)[2],
                             NA),
    "TP_TDN_max_time_3" = ifelse(length(TP_TDN_last_time)>=3,
                                 sort(TP_TDN_last_time, decreasing=TRUE)[3],
                                 NA),
    "TP_TDN_max_span" = max(TP_TDN_span),
    "TP_TDN_span_2" = ifelse(length(TP_TDN_span)>=2,
                             sort(TP_TDN_span, decreasing=TRUE)[2],
                             NA),
    "TP_TDN_span_3" = ifelse(length(TP_TDN_span)>=3,
                             sort(TP_TDN_span, decreasing=TRUE)[3],
                             NA),  
    "TP_TDN_num_sites" = n_distinct(posid),
    "TP_TDN_num_sites_CR" = n_distinct(posid[patient %in% CR_pats]),
    "TP_TDN_num_sites_NR" = n_distinct(posid[patient %in% NR_pats]),
    "TP_TDN_sum_abund" = sum(TP_TDN_sum_abund),
    "TP_TDN_sum_abund_CR" = sum(TP_TDN_sum_abund[patient %in% CR_pats]),
    "TP_TDN_sum_abund_NR" = sum(TP_TDN_sum_abund[patient %in% NR_pats]),
    "TP_TDN_peak_abund" = max(TP_TDN_peak_abund),
    "TP_TDN_peak_abund_CR" = max(c(0,TP_TDN_peak_abund[patient %in% CR_pats])),
    #      "TP_peak_abund_CR" = ifelse(TP_peak_abund_CR < 0, 0, TP_peak_abund_CR),
    "TP_TDN_peak_abund_NR" = max(c(0,TP_TDN_peak_abund[patient %in% NR_pats])),
    #      "TP_peak_abund_NR" = ifelse(TP_peak_abund_NR < 0, 0, TP_peak_abund_NR),
    "TP_TDN_peak_relAbund" = max(TP_TDN_peak_relAbund),
    "TP_TDN_abund_gini" = gintools::pop_calcs(TP_TDN_peak_abund, calc = "gini"),
    "TP_TDN_same" = sum(as.integer(strand == gene_ort)),
    "TP_TDN_same_CR" = sum(as.integer(strand == gene_ort)[patient %in% CR_pats]),
    "TP_TDN_same_NR" = sum(as.integer(strand == gene_ort)[patient %in% NR_pats]),
    "TP_TDN_oppo" = sum(as.integer(strand != gene_ort)),
    "TP_TDN_oppo_CR" = sum(as.integer(strand != gene_ort)[patient %in% CR_pats]),
    "TP_TDN_oppo_NR" = sum(as.integer(strand != gene_ort)[patient %in% NR_pats]),
    .groups = 'drop') %>% 
  rowwise() %>% 
  mutate(TP_TDN_top3avg=mean(c(TP_TDN_max_span,TP_TDN_span_2,TP_TDN_span_3),na.rm=TRUE)) %>% 
  mutate(TP_TDN_time_top3avg=mean(c(TP_TDN_max_time,TP_TDN_max_time_2,TP_TDN_max_time_3),na.rm=TRUE)) %>% 
  ungroup()


####################


sites_by_patient <- split(
  cond_uniq_sites[cond_uniq_sites$timepoint != "d0"], 
  cond_uniq_sites[cond_uniq_sites$timepoint != "d0"]$patient)

get_top_genes <- function(sites, genes, percent, rank_by){
  sites <- sites[order(mcols(sites)[,rank_by], decreasing = TRUE)]
  first_sites <- sites[!duplicated(sites$posid)]
  ranks <- rank(-(mcols(first_sites)[,rank_by]), ties.method = "min")
  cutoff <- round(length(first_sites)*percent/100)
  if(cutoff == 0){ cutoff <- 1 }
  posids <- first_sites[ranks <= cutoff]$posid
  top_sites <- sites[sites$posid %in% posids]
  unique(genes[subjectHits(findOverlaps(
    top_sites, genes, ignore.strand = TRUE))]$name2)
}

refGenes <- readRDS(file.path('utils', "hg38.refSeq.rds"))

top_one_pc_sonicAbund <- paste(
  unlist(lapply(
    sites_by_patient,
    get_top_genes,
    genes = refGenes,
    percent = 1,
    rank_by = "estAbund")),
  collapse = ", ")

top_ten_pc_sonicAbund <- paste(
  unlist(lapply(
    sites_by_patient,
    get_top_genes,
    genes = refGenes,
    percent = 10,
    rank_by = "estAbund")),
  collapse = ", ")

oncoGenesData <- read.delim(
  file.path('utils', "allOnco.human.v3.tsv"),
  header = TRUE, 
  sep = "\t",
  stringsAsFactors = FALSE
)

hgnc_complete <- data.table::fread(
  paste0("zcat ", file.path('utils', "hgnc_complete_set.180207.txt.gz")),
  sep = "\t", header = TRUE, 
  select = c(
    "HGNC ID", "Approved Symbol", "Approved Name", "Locus Group", "Locus Type", 
    "Synonyms", "Previous Symbols", "Entrez Gene ID", "Ensembl Gene ID", 
    "RefSeq (supplied by NCBI)", "UniProt ID (supplied by UniProt)"),
  data.table = FALSE
)

names(hgnc_complete) <- c(
  "hgnc_id", "symbol", "name", "locus_group", "locus_type", "alias_symbol", 
  "prev_symbol", "entrez_id", "ensembl_gene_id", "refseq_accession", 
  "uniprot_ids"
)

hgnc_complete <- dplyr::filter(hgnc_complete, !grepl("withdrawn", symbol)) %>%
  dplyr::mutate(
    kegg_id = paste0("hsa:", entrez_id),
    entrez_id = paste0(entrez_id, ":EZID"),
    alias_symbol = gsub(", ", "|", alias_symbol),
    prev_symbol = gsub(", ", "|", prev_symbol),
    extended_alias = paste0(
      alias_symbol, "|", prev_symbol, "|", ensembl_gene_id, "|", 
      refseq_accession, "|", uniprot_ids)
  )


oncoGenes <- unique(oncoGenesData[,"symbol"]) %>%
  spraphal::alias_arbiter(
    IDs = ., 
    RefIDs = hgnc_complete$symbol,
    aliasIDs = hgnc_complete$extended_alias,
    sep = "|", remove_absent_IDs = NULL, quiet = TRUE
  )


gene_stats <- data.frame(
  "loci" = get_loci_id(mod_refGenes),
  "gene_name" = as.character(mod_refGenes$name),
  "gene_ort" = as.character(strand(mod_refGenes)),
  "gene_width" = (width(mod_refGenes) + 10000) / 1000,
  "Within_Cluster" = mod_refGenes$Within_Cluster,
  "Cluster_target.min" = mod_refGenes$Cluster_target.min) %>%
  dplyr::full_join(., stats_tdn, by = c("loci", "gene_name", "gene_ort")) %>%
  dplyr::full_join(., stats_tp, by = c("loci", "gene_name", "gene_ort")) %>%
  dplyr::full_join(., stats_tp_tdn, by = c("loci", "gene_name", "gene_ort")) %>%
  dplyr::mutate(
    "On_Onco_List" = gene_name %in% oncoGenes,
    "Top_1pc_Abund" = gene_name %in% top_one_pc_sonicAbund,
    "Top_10pc_Abund" = gene_name %in% top_ten_pc_sonicAbund)


gene_stats[is.na(gene_stats)] <- 0
total_tdn_sites <- length(unique(paste(tdn_gr$patient, tdn_gr$posid)))
total_tdn_sites_CR <- length(unique(paste(
  tdn_gr$patient, tdn_gr$posid)[tdn_gr$patient %in% CR_pats]))
total_tdn_sites_NR <- length(unique(paste(
  tdn_gr$patient, tdn_gr$posid)[tdn_gr$patient %in% NR_pats]))
total_pat_sites <- length(unique(paste(tp_gr$patient, tp_gr$posid)))
total_pat_sites_CR <- length(unique(paste(
  tp_gr$patient, tp_gr$posid)[tp_gr$patient %in% CR_pats]))
total_pat_sites_NR <- length(unique(paste(
  tp_gr$patient, tp_gr$posid)[tp_gr$patient %in% NR_pats]))

gene_stats <- dplyr::filter(
  gene_stats, TDN_num_patients > 0 | TP_num_patients > 0
) %>%
  dplyr::group_by(loci, gene_name, gene_ort) %>%
  dplyr::mutate(
    "ort_fisher_test" = fisher.test(matrix(
      c(tdn_same, tdn_oppo, tp_same, tp_oppo), 
      nrow = 2, ncol = 2))$p.value,
    "ort_fisher_test_CR" = fisher.test(matrix(
      c(tdn_same_CR, tdn_oppo_CR, tp_same_CR, tp_oppo_CR), 
      nrow = 2, ncol = 2))$p.value,
    "ort_fisher_test_NR" = fisher.test(matrix(
      c(tdn_same_NR, tdn_oppo_NR, tp_same_NR, tp_oppo_NR), 
      nrow = 2, ncol = 2))$p.value) %>% 
  dplyr::ungroup() %>%
  dplyr::select(
    loci, gene_name, gene_ort, gene_width, TDN_num_patients, TDN_num_pats_CR, 
    TDN_num_pats_NR, TP_num_patients, TP_num_pats_CR, TP_num_pats_NR,
    TDN_num_sites, TDN_num_sites_CR, TDN_num_sites_NR, TP_num_sites, 
    TP_num_sites_CR, TP_num_sites_NR, TDN_peak_abund, TDN_peak_abund_CR, 
    TDN_peak_abund_NR, TP_peak_abund, TP_peak_abund_CR, TP_peak_abund_NR,
    TDN_peak_relAbund, TP_peak_relAbund, TDN_sum_abund, TDN_sum_abund_CR, 
    TDN_sum_abund_NR, TP_sum_abund, TP_sum_abund_CR, TP_sum_abund_NR, 
    long_count, max_time, max_span, abund_gini, ort_fisher_test, 
    ort_fisher_test_CR, ort_fisher_test_NR, On_Onco_List, Top_1pc_Abund, 
    Top_10pc_Abund, Within_Cluster, Cluster_target.min,
    TP_TDN_num_patients,
    TP_TDN_num_pats_CR,
    TP_TDN_num_pats_NR,
    TP_TDN_long_count,
    TP_TDN_max_time,
    TP_TDN_max_span,
    TP_TDN_num_sites,
    TP_TDN_num_sites_CR,
    TP_TDN_num_sites_NR,
    TP_TDN_sum_abund,
    TP_TDN_sum_abund_CR,
    TP_TDN_sum_abund_NR,
    TP_TDN_peak_abund,
    TP_TDN_peak_abund_CR,
    TP_TDN_peak_abund_NR,
    TP_TDN_peak_relAbund,
    TP_TDN_abund_gini,
    TP_TDN_top3avg,
    TP_TDN_time_top3avg
    ) %>%
  dplyr::mutate(
    TDN_freq = TDN_num_sites / (gene_width * total_tdn_sites),
    TDN_freq_CR = TDN_num_sites_CR / (gene_width * total_tdn_sites_CR),
    TDN_freq_NR = TDN_num_sites_NR / (gene_width * total_tdn_sites_NR),
    TP_freq = TP_num_sites/(gene_width * total_pat_sites),
    TP_freq_CR = TP_num_sites_CR/(gene_width * total_pat_sites_CR),
    TP_freq_NR = TP_num_sites_NR/(gene_width * total_pat_sites_NR),
    freq_diff = (TP_freq - TDN_freq),
    freq_diff_CR = (TP_freq_CR - TDN_freq_CR),
    freq_diff_NR = (TP_freq_NR - TDN_freq_NR),
    pct_chg = 100 * freq_diff / TDN_freq,
    pct_chg_CR = 100 * freq_diff_CR / TDN_freq_CR,
    pct_chg_NR = 100 * freq_diff_NR / TDN_freq_NR) %>%
  as.data.frame()



saveRDS(
  gene_stats,
  file=file.path('gene_impact',paste0(TRIAL,'_',RESP,'_gene_impact.rds')))
  


library(tidyverse)
library(GenomicRanges)
library(maditr)
library(gintools)

if(!exists('TRIAL')) {TRIAL <- 'ALL'} # ALL, CLL, CALL= ALL+CLL
if(!exists('RESP')) {RESP <- 'na'} # whether filter by response class, RE= responders, nRE non responders, na= no filter



cond_uniq_sites <- readRDS(here::here('condensed_intsites',paste0(TRIAL,'_',RESP,'_conintsites.rds')))


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

get_top_sites <- function(sites, percent, rank_by){
  sites <- sites[order(mcols(sites)[,rank_by], decreasing = TRUE)]
  first_sites <- sites[!duplicated(sites$posid)]
  ranks <- rank(-(mcols(first_sites)[,rank_by]), ties.method = "min")
  cutoff <- round(length(first_sites)*percent/100)
  if(cutoff == 0){ cutoff <- 1 }
  first_sites[ranks <= cutoff]$posid
}

annotate_scan_clusters <- function(scanned_ranges, tdn_sites, 
                                   timepoint_sites, refGenes){
  scanned_ranges$n_sites_tdn <- scan_sites_count(
    scanned_ranges, tdn_sites)
  scanned_ranges$n_sites_tp <- scan_sites_count(
    scanned_ranges, timepoint_sites)
  scanned_ranges$sum_abund_tdn <- scan_sites_abundance(
    scanned_ranges, tdn_sites)
  scanned_ranges$sum_abund_tp <- scan_sites_abundance(
    scanned_ranges, timepoint_sites)
  scanned_ranges$n_patients_tdn <- scan_patients(
    scanned_ranges, tdn_sites)
  scanned_ranges$n_patients_tp <- scan_patients(
    scanned_ranges, timepoint_sites)
  scanned_ranges$genes_in_cluster <- scan_genes(
    scanned_ranges, refGenes)
  scanned_ranges$in_gene_ort_tdn <- scan_orientation(
    scanned_ranges, tdn_sites)
  scanned_ranges$in_gene_ort_tp <- scan_orientation(
    scanned_ranges, timepoint_sites)
  scanned_ranges$ort_fisher_test <- scan_fisher_test_ort(
    scanned_ranges$in_gene_ort_tdn, scanned_ranges$in_gene_ort_tp)
  scanned_ranges
}

scan_sites_count <- function(clus_gr, sites){
  sapply(1:length(clus_gr), function(i){
    clus <- clus_gr[i]
    hits <- subjectHits(findOverlaps(clus, sites))
    length(unique(hits))
  })
}

scan_sites_abundance <- function(clus_gr, sites){
  sapply(1:length(clus_gr), function(i){
    clus <- clus_gr[i]
    hits <- subjectHits(findOverlaps(clus, sites))
    sum(sites[hits]$estAbund)
  })
}

scan_patients <- function(clus_gr, sites){
  sapply(1:length(clus_gr), function(i){
    clus <- clus_gr[i]
    hits <- subjectHits(findOverlaps(clus, sites))
    length(unique(sites[hits]$patient))
  })
}

scan_genes <- function(clus_gr, refGenes){
  sapply(1:length(clus_gr), function(i){
    clus <- clus_gr[i]
    hits <- subjectHits(findOverlaps(clus, refGenes))
    paste(unique(refGenes[hits]$name2), collapse = ", ")
  })
}

scan_orientation <- function(clus_gr, sites){
  sapply(1:length(clus_gr), function(i){
    clus <- clus_gr[i]
    site_hits <- findOverlaps(clus, sites)
    site_ort <- strand(sites[subjectHits(site_hits)])
    gene_ort <- sites[subjectHits(site_hits)]$in_geneOrt
    if(length(sites[subjectHits(site_hits)]) > 0){
      df <- data.frame(
        "posid" = generate_posid(sites[subjectHits(site_hits)]),
        "site_ort" = as.character(site_ort),
        "gene_ort" = as.character(gene_ort),
        "patient" = sites[subjectHits(site_hits)]$patient,
        stringsAsFactors = FALSE
      )
      df <- distinct(df)
      df$same_orientation <- df$site_ort == df$gene_ort
      score <- paste0(
        "T", length(grep("TRUE", df$same_orientation)), ":",
        "F", length(grep("FALSE", df$same_orientation)), ":",
        "N", length(grep("TRUE", is.na(df$same_orientation)))
      )
    }else{
      score <- "T0:F0:N0"
    }
    score
  })
}

scan_fisher_test_ort <- function(score1, score2){
  sapply(1:length(score1), function(i){
    grp1 <- unlist(strsplit(score1[i], ":"))
    grp2 <- unlist(strsplit(score2[i], ":"))
    x <- matrix(c(
      as.integer(substr(grp1[1], 2, 3)),
      as.integer(substr(grp1[2], 2, 3)),
      as.integer(substr(grp2[1], 2, 3)),
      as.integer(substr(grp2[2], 2, 3))),
      ncol = 2
    )
    fisher.test(x)$p.value
  })
}


####################

sites_by_patient <- split(
  cond_uniq_sites[cond_uniq_sites$timepoint != "d0"], 
  cond_uniq_sites[cond_uniq_sites$timepoint != "d0"]$patient)

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


refGenesOrt <- as.data.frame(refGenes, row.names = NULL) %>%
  dplyr::select(name2, strand) %>%
  dplyr::group_by(name2) %>%
  dplyr::summarise(
    gene_ort = names(sort(table(strand), decreasing = TRUE)[1])
  ) %>%
  dplyr::ungroup() %>% 
  dplyr::rename(gene = name2) %>%
  as.data.frame()


gene_stats <- as.data.frame(cond_uniq_sites, row.names = NULL) %>%
  dplyr::select(
    specimen, patient, celltype, timepoint, estAbund, relAbund, posid, 
    gene_id_wo_annot, in_gene, nearest_geneDist, start, strand, seqnames) %>%
  dplyr::mutate(
    type = ifelse(timepoint == "d0", "TDN", "TP"),
    gene_name = ifelse(is.na(gene_id_wo_annot), posid, gene_id_wo_annot),
    nearest_geneDist = ifelse(is.na(nearest_geneDist), 0, nearest_geneDist)) %>%
  dplyr::rename(ort = strand) %>%
  dplyr::group_by(type, patient, gene_name, posid, ort) %>%
  dplyr::summarise(
    long_count = n_distinct(timepoint),
    peak_abund = max(estAbund),
    peak_relAbund = max(relAbund),
    max_geneDist = max(ifelse(in_gene == FALSE, abs(nearest_geneDist), 0)),
    loci = unique(start)) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(., refGenesOrt, by = c("gene_name" = "gene")) %>%
  dplyr::mutate(gene_ort = ifelse(is.na(gene_ort), "*", gene_ort)) %>%
  dplyr::group_by(type, gene_name, gene_ort) %>%
  dplyr::summarise(
    num_patients = n_distinct(patient),
    num_sites = n_distinct(posid),
    long_count = max(long_count),
    abund_gini = pop_calcs(peak_abund, calc = "gini"),
    peak_abund = max(peak_abund),
    peak_relAbund = max(peak_relAbund),
    max_geneDist = max(max_geneDist),
    min_loci = min(loci),
    max_loci = max(loci),
    same_ort = length(which(ort == gene_ort)),
    oppo_ort = length(which(ort != gene_ort))) %>%
  dplyr::ungroup() %>% 
  as.data.frame() %>%
  melt(
    id.vars = c("type", "gene_name", "gene_ort"), 
    measure.vars = c(
      "num_patients", "num_sites", "long_count", 
      "peak_abund", "peak_relAbund", "abund_gini", "max_geneDist",
      "min_loci", "max_loci", "same_ort", "oppo_ort")) %>%
  dcast(gene_name + gene_ort ~ type + variable, fill = 0) %>%
  dplyr::group_by(gene_name, gene_ort) %>%
  dplyr::mutate(
    long_count = TP_long_count,
    abund_gini = TP_abund_gini,
    max_geneDist = max(TDN_max_geneDist, TP_max_geneDist),
    genomic_range = max(c(TDN_max_loci, TP_max_loci)) - 
      min(c(TDN_min_loci, TP_min_loci)),
    ort_fisher_test = fisher.test(
      matrix(
        c(TDN_same_ort, TDN_oppo_ort, TP_same_ort, TP_oppo_ort), 
        nrow = 2, ncol = 2))$p.value) %>%
  dplyr::select(
    gene_name, gene_ort, TDN_num_patients, TP_num_patients, TDN_num_sites, 
    TP_num_sites, TDN_peak_abund, TP_peak_abund, TDN_peak_relAbund, 
    TP_peak_relAbund, abund_gini, long_count, max_geneDist, genomic_range, 
    ort_fisher_test) %>%
  dplyr::ungroup() %>% 
  as.data.frame()

## Identify genes of interest from clusters and clonal expansions ------------

tdn_sites <- cond_uniq_sites[cond_uniq_sites$timepoint == "d0"]
timepoint_sites <- cond_uniq_sites[cond_uniq_sites$timepoint != "d0"]

### timepoint_clusters
scan_sites <- gintools:::scan_format(
  tdn_sites, timepoint_sites, grouping = "patient")

scanned_clus <- geneRxCluster::gRxCluster(
  object = scan_sites$chr, 
  starts = scan_sites$pos, 
  group = scan_sites$grp,
  kvals = c(10L:50L),
  nperm = 100L,
  cutpt.tail.expr = critVal.target(k, n, target = 7.5, posdiff = x),
  cutpt.filter.expr = apply(x, 2, quantile, probs = 0.35, na.rm = TRUE))


scanned_clus <- annotate_scan_clusters(
  scanned_clus, tdn_sites, timepoint_sites, refGenes)

scanned_clus <- scanned_clus[
  scanned_clus$n_patients_tdn > 1 | scanned_clus$n_patients_tp > 1]


enriched_timepoint_clus <- scanned_clus[
  (scanned_clus$n_sites_tdn/length(tdn_sites)) < 
    (scanned_clus$n_sites_tp/length(timepoint_sites))]

# Scans for clusters with orientation biases -----------------------------------
pos_strand_tdn_sites <- tdn_sites[strand(tdn_sites) == "+"]
pos_strand_timepoint_sites <- timepoint_sites[strand(timepoint_sites) == "+"]

neg_strand_tdn_sites <- tdn_sites[strand(tdn_sites) == "-"]
neg_strand_timepoint_sites <- timepoint_sites[strand(timepoint_sites) == "-"]

pos_strand_scan_sites <- gintools:::scan_format(
  pos_strand_tdn_sites, pos_strand_timepoint_sites, grouping = "patient")
neg_strand_scan_sites <- gintools:::scan_format(
  neg_strand_tdn_sites, neg_strand_timepoint_sites, grouping = "patient")

pos_scanned_clus <- geneRxCluster::gRxCluster(
  object = pos_strand_scan_sites$chr, 
  starts = pos_strand_scan_sites$pos, 
  group = pos_strand_scan_sites$grp,
  kvals = c(10L:50L),
  nperm = 100L,
  cutpt.tail.expr = critVal.target(k, n, target = 10, posdiff = x),
  cutpt.filter.expr = apply(x, 2, quantile, probs = 0.25, na.rm = TRUE))

#pos_scan_summary <- gRxSummary(pos_scanned_clus)

neg_scanned_clus <- geneRxCluster::gRxCluster(
  object = neg_strand_scan_sites$chr, 
  starts = neg_strand_scan_sites$pos, 
  group = neg_strand_scan_sites$grp,
  kvals = c(10L:50L),
  nperm = 100L,
  cutpt.tail.expr = critVal.target(k, n, target = 10, posdiff = x),
  cutpt.filter.expr = apply(x, 2, quantile, probs = 0.25, na.rm = TRUE))

#neg_scan_summary <- gRxSummary(neg_scanned_clus)

## Clusters are compared to all sites, not just the same strand
pos_scanned_clus <- annotate_scan_clusters(
  pos_scanned_clus, tdn_sites[strand(tdn_sites) == "+"], 
  timepoint_sites[strand(timepoint_sites) == "+"], refGenes)

neg_scanned_clus <- annotate_scan_clusters(
  neg_scanned_clus, tdn_sites[strand(tdn_sites) == "-"], 
  timepoint_sites[strand(timepoint_sites) == "-"], refGenes)

pos_scanned_clus <- pos_scanned_clus[
  pos_scanned_clus$n_patients_tdn > 1 | pos_scanned_clus$n_patients_tp > 1]

neg_scanned_clus <- neg_scanned_clus[
  neg_scanned_clus$n_patients_tdn > 1 | neg_scanned_clus$n_patients_tp > 1]

pos_strand_timepoint_enriched_clus <- pos_scanned_clus[
  (pos_scanned_clus$n_sites_tp/length(pos_strand_timepoint_sites)) > 
    (pos_scanned_clus$n_sites_tdn/length(pos_strand_tdn_sites))]

neg_strand_timepoint_enriched_clus <- neg_scanned_clus[
  (neg_scanned_clus$n_sites_tp/length(neg_strand_timepoint_sites)) > 
    (neg_scanned_clus$n_sites_tdn/length(neg_strand_tdn_sites))]

# Scans for clusters in higher abundance ---------------------------------------
# over lower abundance in other timepoints
cutoff <- 2.5 # percent

tp_high_abund_posids <- get_top_sites(timepoint_sites, cutoff, "estAbund")
timepoint_sites$abund_status <- ifelse(
  timepoint_sites$posid %in% tp_high_abund_posids,
  "High Abundance", "Low Abundance")

high_sites <- timepoint_sites[
  timepoint_sites$abund_status == "High Abundance"]
low_sites <- timepoint_sites[
  timepoint_sites$abund_status == "Low Abundance"]

high_low_sites <- gintools:::scan_format(
  low_sites, high_sites, grouping = "patient")

high_low_clus <- geneRxCluster::gRxCluster(
  object = high_low_sites$chr, 
  starts = high_low_sites$pos, 
  group = high_low_sites$grp,
  kvals = c(10L:35L),
  nperm = 100L,
  cutpt.tail.expr = critVal.target(k, n, target = 2, posdiff = x),
  cutpt.filter.expr = apply(x, 2, quantile, probs = 0.15, na.rm = TRUE))

#high_low_summary <- gRxSummary(high_low_clus)

high_low_clus <- high_low_clus[width(high_low_clus) > 1]
high_low_clus$n_sites_high <- scan_sites_count(high_low_clus, high_sites)
high_low_clus$n_sites_low <- scan_sites_count(high_low_clus, low_sites)
high_low_clus$n_patients_high <- scan_patients(high_low_clus, high_sites)
high_low_clus$n_patients_low <- scan_patients(high_low_clus, low_sites)
high_low_clus$genes_in_cluster <- scan_genes(high_low_clus, refGenes)
high_low_clus$in_gene_ort_high <- scan_orientation(high_low_clus, high_sites)
high_low_clus$in_gene_ort_low <- scan_orientation(high_low_clus, low_sites)
high_low_clus$ort_fisher_test <- scan_fisher_test_ort(
  high_low_clus$in_gene_ort_high, high_low_clus$in_gene_ort_low)
enriched_high_clus <- high_low_clus[
  (high_low_clus$n_sites_high/length(high_sites)) > 
    (high_low_clus$n_sites_low/length(low_sites))]
enriched_low_clus <- high_low_clus[
  (high_low_clus$n_sites_high/length(high_sites)) < 
    (high_low_clus$n_sites_low/length(low_sites))]


#### make cluster

clusters <- list(
  "timepoint_clusters" = enriched_timepoint_clus,
  "positive_strand_clusters" = pos_strand_timepoint_enriched_clus,
  "negative_strand_clusters" = neg_strand_timepoint_enriched_clus,
  "high_abund_clusters" = enriched_high_clus)

cart19_clusters <- unlist(GRangesList(lapply(
  1:4,
  function(i){
    clus_name <- c("Timepoint", "Pos-Strand", "Neg-Strand", "Abundance")[i]
    cluster_group <- granges(clusters[[i]])
    cluster_group$clus.origin <- clus_name
    cluster_group$target.min <- clusters[[i]]$target.min
    cluster_group
  }
)))

saveRDS(cart19_clusters,file =file.path('condensed_intsites',paste0(TRIAL,'_',RESP,'_cart19_clusters.rds')))


red_clusters <- GenomicRanges::reduce(cart19_clusters, with.revmap = TRUE)

red_clusters$cluster_origin <- sapply(red_clusters$revmap, function(x){
  paste(unique(cart19_clusters[x]$clus.origin), collapse = ", ")})

red_clusters$revmap <- NULL

red_clusters <- annotate_scan_clusters(
  red_clusters, tdn_sites, timepoint_sites, refGenes)




cluster_gene_list <- unique(unlist(sapply(
  red_clusters$genes_in_cluster, 
  function(x) unlist(strsplit(x, ", "))), use.names = FALSE))


## Annotate with which list the gene belongs =================================

oncoGenesData <- read.delim(
  file.path('utils', "allOnco.human.v3.tsv"),
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

oncoGenes <- unique(oncoGenesData[,"symbol"])

top_clonal_genes <- unique(unlist(strsplit(top_one_pc_sonicAbund, ", ")))
top_ten_pc_clonal_genes <- unique(unlist(strsplit(top_ten_pc_sonicAbund, ", ")))

gene_stats <- gene_stats %>%
  dplyr::mutate(
    "On_Onco_List" = gene_name %in% oncoGenes,
    "Top_1pc_Abund" = gene_name %in% top_clonal_genes,
    "Top_10pc_Abund" = gene_name %in% top_ten_pc_clonal_genes,
    "Within_Cluster" = gene_name %in% cluster_gene_list
  )

all_clusters <- cart19_clusters
all_clusters$genes_in_cluster <- scan_genes(all_clusters, refGenes)

gene_stats$Cluster_target.min <- sapply(
  gene_stats$gene_name, function(gene){
    clusters <- all_clusters[grep(gene, all_clusters$genes_in_cluster)]
    if(length(clusters) == 0){
      target <- NA
    }else{
      target <- min(clusters$target.min)
    }
    target
  })

write.csv(
  gene_stats,
  file = file.path('condensed_intsites',paste0(TRIAL,'_',RESP,'_cart19_gene_stats.csv')),
  quote = TRUE,
  row.names = FALSE
)

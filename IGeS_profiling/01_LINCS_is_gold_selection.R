print("start to prepare expression profile")
work_path <- "/picb/bigdata/project/FengFYM/drug_MoA_reanalysis/"
setwd(work_path)
result_path <- "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/03_drug_immuneSig_enrichment/data/"

dir.create(paste0(result_path, "drug_treated_expression_matrix/"))
dir.create(paste0(result_path, "drug_treated_design/"))
dir.create(paste0(result_path, "drug_treated_design/for_GSEA/"))
dir.create(paste0(result_path, "drug_treated_design/for_GSVA/"))
dir.create("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/06_results_summaryandplotting/data/LINCS_is_gold/")
library(prada)
library(rhdf5)
library(limma)

source("l1ktools-master/R/cmap/io.R")
options(stringsAsFactors = F)
options(quote = NULL)

dataset <- "70138"

inst_GSE70138 <- read.csv("LINCS_data/GSE70138/file/GSE70138_Broad_LINCS_inst_info.txt", header = T, sep = "\t")
sig_GSE70138 <-read.csv("LINCS_data/GSE70138/file/GSE70138_Broad_LINCS_sig_metrics.txt", header = T, sep = "\t")

# is_gold: distil_cc_q75 >= 0.2 and pct_self_rank_q25 <= 0.05.
sig_GSE70138_is_gold <- sig_GSE70138[sig_GSE70138$distil_cc_q75 >= 0.2 & sig_GSE70138$pct_self_rank_q25 <= 0.05, ]


inst_GSE70138_is_gold <- inst_GSE70138[is.element(inst_GSE70138$pert_id, unique(sig_GSE70138_is_gold$pert_id)), ]

sort(unique(inst_GSE70138_is_gold[inst_GSE70138_is_gold$cell_id == "A375", ]$pert_iname))

tmp = unique(inst_GSE70138_is_gold[c("cell_id", "pert_iname")])
table(tmp$cell_id)
write.table(inst_GSE70138_is_gold, "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/06_results_summaryandplotting/data/LINCS_is_gold/inst_GSE70138_is_gold.csv", quote = F, row.names = F, sep = '\t')

# GSE92742

inst_GSE92742 <- read.csv("LINCS_data/GSE92742/file/GSE92742_Broad_LINCS_inst_info.txt", header=T, sep="\t")
sig_GSE92742  <-read.csv("LINCS_data/GSE92742/file/GSE92742_Broad_LINCS_sig_metrics.txt", header = T, sep = "\t")

sig_GSE92742_is_gold <- sig_GSE92742[sig_GSE92742$distil_cc_q75 >= 0.2 & sig_GSE92742$pct_self_rank_q25 <= 0.05, ]


inst_GSE92742_is_gold <- inst_GSE92742[is.element(inst_GSE92742$pert_id, unique(sig_GSE92742_is_gold$pert_id)), ]

sort(unique(inst_GSE92742_is_gold[inst_GSE92742_is_gold$cell_id == "A375", ]$pert_iname))

tmp = unique(inst_GSE92742_is_gold[c("cell_id", "pert_iname")])
table(tmp$cell_id)
write.table(inst_GSE92742_is_gold, "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/06_results_summaryandplotting/data/LINCS_is_gold/inst_GSE92742_is_gold.csv", quote = F, row.names = F, sep = '\t')

# rm(list=ls())
work_path <- "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/"
setwd(work_path)
options(stringsAsFactors = F)

library(reshape2)
library(ggplot2)
library(ggplotify)
library(grid)
library(cowplot)
library(parallel)
library(patchwork)
library(dplyr)
library(stringr)
library(data.table)
require(GGally)
library(plot3D)
library(fmsb)
library(readxl)
library(dplyr)
library(UpSetR)
library(clusterProfiler)
library(enrichplot)
# source("06_results_summaryandplotting/scripts/plot_for_drug_function.R")
# source("06_results_summaryandplotting/scripts/ITS_analysis/ITS_gene_stat_function.R")
# source("06_results_summaryandplotting/scripts/immunesig_selection_function.R")
# source("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/scripts/drugITGSscore_lincs.R")
source("07_plotting_v2/scripts/functions/ITS_enrich.R")
source("07_plotting_v2/scripts/functions/ITS_function_detection.R")



# ----------------------------------------------------------------------------

dir.create("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/")
dir.create("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSp_sharedgene/drug_PMID35292771")
dir.create("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSp_sharedgene/drug_PMID35292771/GSEAenrich/")
dir.create("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSn_sharedgene/drug_PMID35292771")
dir.create("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSn_sharedgene/drug_PMID35292771/GSEAenrich/")

# LOAD GENES IN LINCS ------------------------------------------------------
lincs_gene = read.csv("/picb/bigdata/project/FengFYM/drug_MoA_reanalysis/LINCS_data/GSE70138/file/GSE92742_Broad_LINCS_gene_info.txt", sep = "\t")

# LOAD drug list from publication ------------------------------------------------------

drug_oncogenic <- as.data.frame(read_excel("07_plotting_v2/data/drug_target_lists_publications/PMID35292771_Tables.xlsx", sheet = 1))
drug_nononcogenic <- as.data.frame(read_excel("07_plotting_v2/data/drug_target_lists_publications/PMID35292771_Tables.xlsx", sheet = 2))


# ------------------------------------------------------------------------------
load("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITS_set.Rdata")

# ITS_p_s ------------------------------------------------------------------------------
ITS_p_s_notspecific_df <- read.table("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSp_sharedgene/ITS_p_s_notspecific_df.txt", header = T, sep = '\t')
# head(ITS_p_s_notspecific_df)



# enrich custom-GSEA -----------------------------------------------------------

load("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSp_sharedgene/enrich/custom_GSEA_s_scaled.Rdata")
res2 = data.frame(res)

# 1) drugTarget_drugbank-Drug -------------------------------------------------------------------
filepath <- ("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSp_sharedgene/drug_PMID35292771/GSEAenrich/ITSp_s_vs_drugtarget")

glist <- unlist(strsplit(res2[res2$ID == "ITSp_vs_drugtarget",]$core_enrichment, "/"))

druglist_OG <- unique(drug_oncogenic[is.element(drug_oncogenic$genename, glist),]$pert_iname)
drug_OG <- drug_oncogenic[is.element(drug_oncogenic$pert_iname, druglist_OG),]

druglist_nonOG <- unique(drug_nononcogenic[is.element(drug_nononcogenic$genename, glist),]$pert_iname)
drug_nonOG <- drug_nononcogenic[is.element(drug_nononcogenic$pert_iname, druglist_nonOG),]

write.table(drug_OG, paste0(filepath, "_OG.txt"), sep = '\t', quote = F, row.names = F)
write.table(drug_nonOG, paste0(filepath, "_nonOG.txt"), sep = '\t', quote = F, row.names = F)



# 2) OG-Drug -------------------------------------------------------------------
filepath <- ("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSp_sharedgene/drug_PMID35292771/GSEAenrich/ITSp_s_vs_OG_oncoKB")

glist <- unlist(strsplit(res2[res2$ID == "ITSp_vs_OG_oncoKB",]$core_enrichment, "/"))

druglist_OG <- unique(drug_oncogenic[is.element(drug_oncogenic$genename, glist),]$pert_iname)
drug_OG <- drug_oncogenic[is.element(drug_oncogenic$pert_iname, druglist_OG),]

druglist_nonOG <- unique(drug_nononcogenic[is.element(drug_nononcogenic$genename, glist),]$pert_iname)
drug_nonOG <- drug_nononcogenic[is.element(drug_nononcogenic$pert_iname, druglist_nonOG),]

write.table(drug_OG, paste0(filepath, "_OG.txt"), sep = '\t', quote = F, row.names = F)
write.table(drug_nonOG, paste0(filepath, "_nonOG.txt"), sep = '\t', quote = F, row.names = F)



# 3) TSG-Drug -------------------------------------------------------------------
filepath <- ("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSp_sharedgene/drug_PMID35292771/GSEAenrich/ITSp_s_vs_TSG_oncoKB")

glist <- unlist(strsplit(res2[res2$ID == "ITSp_vs_TSG_oncoKB",]$core_enrichment, "/"))


druglist_OG <- unique(drug_oncogenic[is.element(drug_oncogenic$genename, glist),]$pert_iname)
drug_OG <- drug_oncogenic[is.element(drug_oncogenic$pert_iname, druglist_OG),]

druglist_nonOG <- unique(drug_nononcogenic[is.element(drug_nononcogenic$genename, glist),]$pert_iname)
drug_nonOG <- drug_nononcogenic[is.element(drug_nononcogenic$pert_iname, druglist_nonOG),]

write.table(drug_OG, paste0(filepath, "_OG.txt"), sep = '\t', quote = F, row.names = F)
write.table(drug_nonOG, paste0(filepath, "_nonOG.txt"), sep = '\t', quote = F, row.names = F)





# enrich custom-ORA ---------------------------------------------------------------
dir.create("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSp_sharedgene/drug_PMID35292771/ORAenrich/")
filepath = "07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSp_sharedgene/enrich/"
load(paste0(filepath, "custom_ORA_s.Rdata"))

res2 = data.frame(res)

# 1) drugTarget_drugbank-Drug -------------------------------------------------------------------
filepath <- ("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSp_sharedgene/drug_PMID35292771/ORAenrich/ITSp_s_vs_drugtarget")

glist <- unlist(strsplit(res2[res2$ID == "ITSp_vs_drugtarget",]$geneID, "/"))
glist <- ITS_p_s_notspecific_df[is.element(ITS_p_s_notspecific_df$ENTREZID, glist), ]$gene_name


druglist_OG <- unique(drug_oncogenic[is.element(drug_oncogenic$genename, glist),]$pert_iname)
drug_OG <- drug_oncogenic[is.element(drug_oncogenic$pert_iname, druglist_OG),]

druglist_nonOG <- unique(drug_nononcogenic[is.element(drug_nononcogenic$genename, glist),]$pert_iname)
drug_nonOG <- drug_nononcogenic[is.element(drug_nononcogenic$pert_iname, druglist_nonOG),]

write.table(drug_OG, paste0(filepath, "_OG.txt"), sep = '\t', quote = F, row.names = F)
write.table(drug_nonOG, paste0(filepath, "_nonOG.txt"), sep = '\t', quote = F, row.names = F)





# 2) OG-Drug -------------------------------------------------------------------
filepath <- ("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSp_sharedgene/drug_PMID35292771/ORAenrich/ITSp_s_vs_OG_oncoKB")

glist <- unlist(strsplit(res2[res2$ID == "ITSp_vs_OG_oncoKB",]$geneID, "/"))
glist <- ITS_p_s_notspecific_df[is.element(ITS_p_s_notspecific_df$ENTREZID, glist), ]$gene_name


druglist_OG <- unique(drug_oncogenic[is.element(drug_oncogenic$genename, glist),]$pert_iname)
drug_OG <- drug_oncogenic[is.element(drug_oncogenic$pert_iname, druglist_OG),]

druglist_nonOG <- unique(drug_nononcogenic[is.element(drug_nononcogenic$genename, glist),]$pert_iname)
drug_nonOG <- drug_nononcogenic[is.element(drug_nononcogenic$pert_iname, druglist_nonOG),]

write.table(drug_OG, paste0(filepath, "_OG.txt"), sep = '\t', quote = F, row.names = F)
write.table(drug_nonOG, paste0(filepath, "_nonOG.txt"), sep = '\t', quote = F, row.names = F)


# 3) TSG-Drug -------------------------------------------------------------------
filepath <- ("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSp_sharedgene/drug_PMID35292771/ORAenrich/ITSp_s_vs_TSG_oncoKB")

glist <- unlist(strsplit(res2[res2$ID == "ITSp_vs_TSG_oncoKB",]$geneID, "/"))
glist <- ITS_p_s_notspecific_df[is.element(ITS_p_s_notspecific_df$ENTREZID, glist), ]$gene_name


druglist_OG <- unique(drug_oncogenic[is.element(drug_oncogenic$genename, glist),]$pert_iname)
drug_OG <- drug_oncogenic[is.element(drug_oncogenic$pert_iname, druglist_OG),]

druglist_nonOG <- unique(drug_nononcogenic[is.element(drug_nononcogenic$genename, glist),]$pert_iname)
drug_nonOG <- drug_nononcogenic[is.element(drug_nononcogenic$pert_iname, druglist_nonOG),]

write.table(drug_OG, paste0(filepath, "_OG.txt"), sep = '\t', quote = F, row.names = F)
write.table(drug_nonOG, paste0(filepath, "_nonOG.txt"), sep = '\t', quote = F, row.names = F)





# ITS_p_r ------------------------------------------------------------------------------
ITS_p_r_notspecific_df <- read.table("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSp_sharedgene/ITS_p_r_notspecific_df.txt", header = T, sep = '\t')
# head(ITS_p_r_notspecific_df)



# enrich custom-GSEA -----------------------------------------------------------

load("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSp_sharedgene/enrich/custom_GSEA_r_scaled.Rdata")
res2 = data.frame(res)

# 1) drugTarget_drugbank-Drug -------------------------------------------------------------------
filepath <- ("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSp_sharedgene/drug_PMID35292771/GSEAenrich/ITSp_r_vs_drugtarget")

glist <- unlist(strsplit(res2[res2$ID == "ITSp_vs_drugtarget",]$core_enrichment, "/"))


druglist_OG <- unique(drug_oncogenic[is.element(drug_oncogenic$genename, glist),]$pert_iname)
drug_OG <- drug_oncogenic[is.element(drug_oncogenic$pert_iname, druglist_OG),]

druglist_nonOG <- unique(drug_nononcogenic[is.element(drug_nononcogenic$genename, glist),]$pert_iname)
drug_nonOG <- drug_nononcogenic[is.element(drug_nononcogenic$pert_iname, druglist_nonOG),]

write.table(drug_OG, paste0(filepath, "_OG.txt"), sep = '\t', quote = F, row.names = F)
write.table(drug_nonOG, paste0(filepath, "_nonOG.txt"), sep = '\t', quote = F, row.names = F)


# 2) OG-Drug -------------------------------------------------------------------
filepath <- ("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSp_sharedgene/drug_PMID35292771/GSEAenrich/ITSp_r_vs_OG_oncoKB")

glist <- unlist(strsplit(res2[res2$ID == "ITSp_vs_OG_oncoKB",]$core_enrichment, "/"))

druglist_OG <- unique(drug_oncogenic[is.element(drug_oncogenic$genename, glist),]$pert_iname)
drug_OG <- drug_oncogenic[is.element(drug_oncogenic$pert_iname, druglist_OG),]

druglist_nonOG <- unique(drug_nononcogenic[is.element(drug_nononcogenic$genename, glist),]$pert_iname)
drug_nonOG <- drug_nononcogenic[is.element(drug_nononcogenic$pert_iname, druglist_nonOG),]

write.table(drug_OG, paste0(filepath, "_OG.txt"), sep = '\t', quote = F, row.names = F)
write.table(drug_nonOG, paste0(filepath, "_nonOG.txt"), sep = '\t', quote = F, row.names = F)




# 3) TSG-Drug -------------------------------------------------------------------
filepath <- ("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSp_sharedgene/drug_PMID35292771/GSEAenrich/ITSp_r_vs_TSG_oncoKB")

glist <- unlist(strsplit(res2[res2$ID == "ITSp_vs_TSG_oncoKB",]$core_enrichment, "/"))


druglist_OG <- unique(drug_oncogenic[is.element(drug_oncogenic$genename, glist),]$pert_iname)
drug_OG <- drug_oncogenic[is.element(drug_oncogenic$pert_iname, druglist_OG),]

druglist_nonOG <- unique(drug_nononcogenic[is.element(drug_nononcogenic$genename, glist),]$pert_iname)
drug_nonOG <- drug_nononcogenic[is.element(drug_nononcogenic$pert_iname, druglist_nonOG),]

write.table(drug_OG, paste0(filepath, "_OG.txt"), sep = '\t', quote = F, row.names = F)
write.table(drug_nonOG, paste0(filepath, "_nonOG.txt"), sep = '\t', quote = F, row.names = F)



# enrich custom-ORA ---------------------------------------------------------------
dir.create("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSp_sharedgene/drug_PMID35292771/ORAenrich/")
filepath = "07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSp_sharedgene/enrich/"
load(paste0(filepath, "custom_ORA_s.Rdata"))

res2 = data.frame(res)

# 1) drugTarget_drugbank-Drug -------------------------------------------------------------------
filepath <- ("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSp_sharedgene/drug_PMID35292771/ORAenrich/ITS_r_vs_drugtarget")

glist <- unlist(strsplit(res2[res2$ID == "ITSp_vs_drugtarget",]$geneID, "/"))
glist <- ITS_p_r_notspecific_df[is.element(ITS_p_r_notspecific_df$ENTREZID, glist), ]$gene_name


druglist_OG <- unique(drug_oncogenic[is.element(drug_oncogenic$genename, glist),]$pert_iname)
drug_OG <- drug_oncogenic[is.element(drug_oncogenic$pert_iname, druglist_OG),]

druglist_nonOG <- unique(drug_nononcogenic[is.element(drug_nononcogenic$genename, glist),]$pert_iname)
drug_nonOG <- drug_nononcogenic[is.element(drug_nononcogenic$pert_iname, druglist_nonOG),]

write.table(drug_OG, paste0(filepath, "_OG.txt"), sep = '\t', quote = F, row.names = F)
write.table(drug_nonOG, paste0(filepath, "_nonOG.txt"), sep = '\t', quote = F, row.names = F)




# 2) OG-Drug -------------------------------------------------------------------
filepath <- ("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSp_sharedgene/drug_PMID35292771/ORAenrich/ITSp_r_vs_OG_oncoKB")

glist <- unlist(strsplit(res2[res2$ID == "ITSp_vs_OG_oncoKB",]$geneID, "/"))
glist <- ITS_p_r_notspecific_df[is.element(ITS_p_r_notspecific_df$ENTREZID, glist), ]$gene_name


druglist_OG <- unique(drug_oncogenic[is.element(drug_oncogenic$genename, glist),]$pert_iname)
drug_OG <- drug_oncogenic[is.element(drug_oncogenic$pert_iname, druglist_OG),]

druglist_nonOG <- unique(drug_nononcogenic[is.element(drug_nononcogenic$genename, glist),]$pert_iname)
drug_nonOG <- drug_nononcogenic[is.element(drug_nononcogenic$pert_iname, druglist_nonOG),]

write.table(drug_OG, paste0(filepath, "_OG.txt"), sep = '\t', quote = F, row.names = F)
write.table(drug_nonOG, paste0(filepath, "_nonOG.txt"), sep = '\t', quote = F, row.names = F)




# 3) TSG-Drug -------------------------------------------------------------------
filepath <- ("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSp_sharedgene/drug_PMID35292771/ORAenrich/ITSp_r_vs_TSG_oncoKB")

glist <- unlist(strsplit(res2[res2$ID == "ITSp_vs_TSG_oncoKB",]$geneID, "/"))
glist <- ITS_p_r_notspecific_df[is.element(ITS_p_r_notspecific_df$ENTREZID, glist), ]$gene_name


druglist_OG <- unique(drug_oncogenic[is.element(drug_oncogenic$genename, glist),]$pert_iname)
drug_OG <- drug_oncogenic[is.element(drug_oncogenic$pert_iname, druglist_OG),]

druglist_nonOG <- unique(drug_nononcogenic[is.element(drug_nononcogenic$genename, glist),]$pert_iname)
drug_nonOG <- drug_nononcogenic[is.element(drug_nononcogenic$pert_iname, druglist_nonOG),]

write.table(drug_OG, paste0(filepath, "_OG.txt"), sep = '\t', quote = F, row.names = F)
write.table(drug_nonOG, paste0(filepath, "_nonOG.txt"), sep = '\t', quote = F, row.names = F)





# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ITS_n   ----------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
dir.create("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSn_sharedgene/drug_PMID35292771/")
dir.create("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSn_sharedgene/drug_PMID35292771/GSEAenrich/")


# ITS_n_s ------------------------------------------------------------------------------
ITS_n_s_notspecific_df <- read.table("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSn_sharedgene/ITS_n_s_notspecific_df.txt", header = T, sep = '\t')
# head(ITS_n_s_notspecific_df)



# enrich custom-GSEA -----------------------------------------------------------

load("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSn_sharedgene/enrich/custom_GSEA_s_scaled.Rdata")
res2 = data.frame(res)

# 1) drugTarget_drugbank-Drug -------------------------------------------------------------------
filepath <- ("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSn_sharedgene/drug_PMID35292771/GSEAenrich/ITSn_s_vs_drugtarget")

glist <- unlist(strsplit(res2[res2$ID == "ITSp_vs_drugtarget",]$core_enrichment, "/"))

druglist_OG <- unique(drug_oncogenic[is.element(drug_oncogenic$genename, glist),]$pert_iname)
drug_OG <- drug_oncogenic[is.element(drug_oncogenic$pert_iname, druglist_OG),]

druglist_nonOG <- unique(drug_nononcogenic[is.element(drug_nononcogenic$genename, glist),]$pert_iname)
drug_nonOG <- drug_nononcogenic[is.element(drug_nononcogenic$pert_iname, druglist_nonOG),]

write.table(drug_OG, paste0(filepath, "_OG.txt"), sep = '\t', quote = F, row.names = F)
write.table(drug_nonOG, paste0(filepath, "_nonOG.txt"), sep = '\t', quote = F, row.names = F)


# 2) OG-Drug -------------------------------------------------------------------
filepath <- ("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSn_sharedgene/drug_PMID35292771/GSEAenrich/ITSn_s_vs_OG_oncoKB")

glist <- unlist(strsplit(res2[res2$ID == "ITSp_vs_OG_oncoKB",]$core_enrichment, "/"))


druglist_OG <- unique(drug_oncogenic[is.element(drug_oncogenic$genename, glist),]$pert_iname)
drug_OG <- drug_oncogenic[is.element(drug_oncogenic$pert_iname, druglist_OG),]

druglist_nonOG <- unique(drug_nononcogenic[is.element(drug_nononcogenic$genename, glist),]$pert_iname)
drug_nonOG <- drug_nononcogenic[is.element(drug_nononcogenic$pert_iname, druglist_nonOG),]

write.table(drug_OG, paste0(filepath, "_OG.txt"), sep = '\t', quote = F, row.names = F)
write.table(drug_nonOG, paste0(filepath, "_nonOG.txt"), sep = '\t', quote = F, row.names = F)




# 3) TSG-Drug -------------------------------------------------------------------
filepath <- ("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSn_sharedgene/drug_PMID35292771/GSEAenrich/ITSn_s_vs_TSG_oncoKB")

glist <- unlist(strsplit(res2[res2$ID == "ITSp_vs_TSG_oncoKB",]$core_enrichment, "/"))


druglist_OG <- unique(drug_oncogenic[is.element(drug_oncogenic$genename, glist),]$pert_iname)
drug_OG <- drug_oncogenic[is.element(drug_oncogenic$pert_iname, druglist_OG),]

druglist_nonOG <- unique(drug_nononcogenic[is.element(drug_nononcogenic$genename, glist),]$pert_iname)
drug_nonOG <- drug_nononcogenic[is.element(drug_nononcogenic$pert_iname, druglist_nonOG),]

write.table(drug_OG, paste0(filepath, "_OG.txt"), sep = '\t', quote = F, row.names = F)
write.table(drug_nonOG, paste0(filepath, "_nonOG.txt"), sep = '\t', quote = F, row.names = F)




# enrich custom-ORA ---------------------------------------------------------------
dir.create("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSn_sharedgene/drug_PMID35292771/ORAenrich/")
filepath = "07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSn_sharedgene/enrich/"
load(paste0(filepath, "custom_ORA_s.Rdata"))

res2 = data.frame(res)

# 1) drugTarget_drugbank-Drug -------------------------------------------------------------------
filepath <- ("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSn_sharedgene/drug_PMID35292771/ORAenrich/ITSp_vs_drugtarget")

glist <- unlist(strsplit(res2[res2$ID == "ITSp_vs_drugtarget",]$geneID, "/"))
glist <- ITS_n_s_notspecific_df[is.element(ITS_n_s_notspecific_df$ENTREZID, glist), ]$gene_name


druglist_OG <- unique(drug_oncogenic[is.element(drug_oncogenic$genename, glist),]$pert_iname)
drug_OG <- drug_oncogenic[is.element(drug_oncogenic$pert_iname, druglist_OG),]

druglist_nonOG <- unique(drug_nononcogenic[is.element(drug_nononcogenic$genename, glist),]$pert_iname)
drug_nonOG <- drug_nononcogenic[is.element(drug_nononcogenic$pert_iname, druglist_nonOG),]

write.table(drug_OG, paste0(filepath, "_OG.txt"), sep = '\t', quote = F, row.names = F)
write.table(drug_nonOG, paste0(filepath, "_nonOG.txt"), sep = '\t', quote = F, row.names = F)


# 2) OG-Drug -------------------------------------------------------------------
filepath <- ("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSn_sharedgene/drug_PMID35292771/ORAenrich/ITSn_s_vs_OG_oncoKB")

glist <- unlist(strsplit(res2[res2$ID == "ITSp_vs_OG_oncoKB",]$geneID, "/"))
glist <- ITS_n_s_notspecific_df[is.element(ITS_n_s_notspecific_df$ENTREZID, glist), ]$gene_name


druglist_OG <- unique(drug_oncogenic[is.element(drug_oncogenic$genename, glist),]$pert_iname)
drug_OG <- drug_oncogenic[is.element(drug_oncogenic$pert_iname, druglist_OG),]

druglist_nonOG <- unique(drug_nononcogenic[is.element(drug_nononcogenic$genename, glist),]$pert_iname)
drug_nonOG <- drug_nononcogenic[is.element(drug_nononcogenic$pert_iname, druglist_nonOG),]

write.table(drug_OG, paste0(filepath, "_OG.txt"), sep = '\t', quote = F, row.names = F)
write.table(drug_nonOG, paste0(filepath, "_nonOG.txt"), sep = '\t', quote = F, row.names = F)



# 3) TSG-Drug -------------------------------------------------------------------
filepath <- ("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSn_sharedgene/drug_PMID35292771/ORAenrich/ITSn_s_vs_TSG_oncoKB")

glist <- unlist(strsplit(res2[res2$ID == "ITSp_vs_TSG_oncoKB",]$geneID, "/"))
glist <- ITS_n_s_notspecific_df[is.element(ITS_n_s_notspecific_df$ENTREZID, glist), ]$gene_name


druglist_OG <- unique(drug_oncogenic[is.element(drug_oncogenic$genename, glist),]$pert_iname)
drug_OG <- drug_oncogenic[is.element(drug_oncogenic$pert_iname, druglist_OG),]

druglist_nonOG <- unique(drug_nononcogenic[is.element(drug_nononcogenic$genename, glist),]$pert_iname)
drug_nonOG <- drug_nononcogenic[is.element(drug_nononcogenic$pert_iname, druglist_nonOG),]

write.table(drug_OG, paste0(filepath, "_OG.txt"), sep = '\t', quote = F, row.names = F)
write.table(drug_nonOG, paste0(filepath, "_nonOG.txt"), sep = '\t', quote = F, row.names = F)






# ITS_n_r ------------------------------------------------------------------------------
ITS_n_r_notspecific_df <- read.table("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSn_sharedgene/ITS_n_r_notspecific_df.txt", header = T, sep = '\t')
# head(ITS_n_r_notspecific_df)



# enrich custom-GSEA -----------------------------------------------------------

load("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSn_sharedgene/enrich/custom_GSEA_r_scaled.Rdata")
res2 = data.frame(res)

# 1) drugTarget_drugbank-Drug -------------------------------------------------------------------
filepath <- ("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSn_sharedgene/drug_PMID35292771/GSEAenrich/ITSn_r_vs_drugtarget")

glist <- unlist(strsplit(res2[res2$ID == "ITSp_vs_drugtarget",]$core_enrichment, "/"))

druglist_OG <- unique(drug_oncogenic[is.element(drug_oncogenic$genename, glist),]$pert_iname)
drug_OG <- drug_oncogenic[is.element(drug_oncogenic$pert_iname, druglist_OG),]

druglist_nonOG <- unique(drug_nononcogenic[is.element(drug_nononcogenic$genename, glist),]$pert_iname)
drug_nonOG <- drug_nononcogenic[is.element(drug_nononcogenic$pert_iname, druglist_nonOG),]

write.table(drug_OG, paste0(filepath, "_OG.txt"), sep = '\t', quote = F, row.names = F)
write.table(drug_nonOG, paste0(filepath, "_nonOG.txt"), sep = '\t', quote = F, row.names = F)

# 2) OG-Drug -------------------------------------------------------------------
filepath <- ("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSn_sharedgene/drug_PMID35292771/GSEAenrich/ITSn_r_vs_OG_oncoKB")

glist <- unlist(strsplit(res2[res2$ID == "ITSp_vs_OG_oncoKB",]$core_enrichment, "/"))

druglist_OG <- unique(drug_oncogenic[is.element(drug_oncogenic$genename, glist),]$pert_iname)
drug_OG <- drug_oncogenic[is.element(drug_oncogenic$pert_iname, druglist_OG),]

druglist_nonOG <- unique(drug_nononcogenic[is.element(drug_nononcogenic$genename, glist),]$pert_iname)
drug_nonOG <- drug_nononcogenic[is.element(drug_nononcogenic$pert_iname, druglist_nonOG),]

write.table(drug_OG, paste0(filepath, "_OG.txt"), sep = '\t', quote = F, row.names = F)
write.table(drug_nonOG, paste0(filepath, "_nonOG.txt"), sep = '\t', quote = F, row.names = F)



# 3) TSG-Drug -------------------------------------------------------------------
filepath <- ("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSn_sharedgene/drug_PMID35292771/GSEAenrich/ITSn_r_vs_TSG_oncoKB")

glist <- unlist(strsplit(res2[res2$ID == "ITSp_vs_TSG_oncoKB",]$core_enrichment, "/"))


druglist_OG <- unique(drug_oncogenic[is.element(drug_oncogenic$genename, glist),]$pert_iname)
drug_OG <- drug_oncogenic[is.element(drug_oncogenic$pert_iname, druglist_OG),]

druglist_nonOG <- unique(drug_nononcogenic[is.element(drug_nononcogenic$genename, glist),]$pert_iname)
drug_nonOG <- drug_nononcogenic[is.element(drug_nononcogenic$pert_iname, druglist_nonOG),]

write.table(drug_OG, paste0(filepath, "_OG.txt"), sep = '\t', quote = F, row.names = F)
write.table(drug_nonOG, paste0(filepath, "_nonOG.txt"), sep = '\t', quote = F, row.names = F)




# enrich custom-ORA ------------------------------------------
# ITSp_vs_drugtarget---------------------
dir.create("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSn_sharedgene/drug_PMID35292771/ORAenrich/")
filepath = "07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSn_sharedgene/enrich/"
load(paste0(filepath, "custom_ORA_r.Rdata"))

res2 = data.frame(res)

# 1) drugTarget_drugbank-Drug -------------------------------------------------------------------
filepath <- ("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSn_sharedgene/drug_PMID35292771/ORAenrich/ITSn_r_vs_drugtarget")

glist <- unlist(strsplit(res2[res2$ID == "ITSp_vs_drugtarget",]$geneID, "/"))
glist <- ITS_n_r_notspecific_df[is.element(ITS_n_r_notspecific_df$ENTREZID, glist), ]$gene_name

druglist_OG <- unique(drug_oncogenic[is.element(drug_oncogenic$genename, glist),]$pert_iname)
drug_OG <- drug_oncogenic[is.element(drug_oncogenic$pert_iname, druglist_OG),]

druglist_nonOG <- unique(drug_nononcogenic[is.element(drug_nononcogenic$genename, glist),]$pert_iname)
drug_nonOG <- drug_nononcogenic[is.element(drug_nononcogenic$pert_iname, druglist_nonOG),]

write.table(drug_OG, paste0(filepath, "_OG.txt"), sep = '\t', quote = F, row.names = F)
write.table(drug_nonOG, paste0(filepath, "_nonOG.txt"), sep = '\t', quote = F, row.names = F)


# 2) OG-Drug -------------------------------------------------------------------
filepath <- ("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSn_sharedgene/drug_PMID35292771/ORAenrich/ITSn_r_vs_OG_oncoKB")

glist <- unlist(strsplit(res2[res2$ID == "ITSp_vs_OG_oncoKB",]$geneID, "/"))
glist <- ITS_n_r_notspecific_df[is.element(ITS_n_r_notspecific_df$ENTREZID, glist), ]$gene_name

druglist_OG <- unique(drug_oncogenic[is.element(drug_oncogenic$genename, glist),]$pert_iname)
drug_OG <- drug_oncogenic[is.element(drug_oncogenic$pert_iname, druglist_OG),]

druglist_nonOG <- unique(drug_nononcogenic[is.element(drug_nononcogenic$genename, glist),]$pert_iname)
drug_nonOG <- drug_nononcogenic[is.element(drug_nononcogenic$pert_iname, druglist_nonOG),]

write.table(drug_OG, paste0(filepath, "_OG.txt"), sep = '\t', quote = F, row.names = F)
write.table(drug_nonOG, paste0(filepath, "_nonOG.txt"), sep = '\t', quote = F, row.names = F)

# 3) TSG-Drug -------------------------------------------------------------------
filepath <- ("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSn_sharedgene/drug_PMID35292771/ORAenrich/ITSn_r_vs_TSG_oncoKB")

glist <- unlist(strsplit(res2[res2$ID == "ITSp_vs_TSG_oncoKB",]$geneID, "/"))
glist <- ITS_n_r_notspecific_df[is.element(ITS_n_r_notspecific_df$ENTREZID, glist), ]$gene_name

druglist_OG <- unique(drug_oncogenic[is.element(drug_oncogenic$genename, glist),]$pert_iname)
drug_OG <- drug_oncogenic[is.element(drug_oncogenic$pert_iname, druglist_OG),]

druglist_nonOG <- unique(drug_nononcogenic[is.element(drug_nononcogenic$genename, glist),]$pert_iname)
drug_nonOG <- drug_nononcogenic[is.element(drug_nononcogenic$pert_iname, druglist_nonOG),]

write.table(drug_OG, paste0(filepath, "_OG.txt"), sep = '\t', quote = F, row.names = F)
write.table(drug_nonOG, paste0(filepath, "_nonOG.txt"), sep = '\t', quote = F, row.names = F)

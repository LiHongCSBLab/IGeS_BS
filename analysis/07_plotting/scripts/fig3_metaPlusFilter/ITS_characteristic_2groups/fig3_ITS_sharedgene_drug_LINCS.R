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
dir.create("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSp_sharedgene/")
dir.create("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSp_sharedgene/drug_LINCS")


# LOAD GENES IN LINCS ------------------------------------------------------
lincs_gene = read.csv("/picb/bigdata/project/FengFYM/drug_MoA_reanalysis/LINCS_data/GSE70138/file/GSE92742_Broad_LINCS_gene_info.txt", sep = "\t")

# LOAD drug-target-moa IN LINCS ------------------------------------------------------
drugTargetMoA <- read.table(paste0(work_path, "06_results_summaryandplotting/data/drugTargetMoA_merged.txt"), sep = '\t', header = T)

# ------------------------------------------------------------------------------
load( "07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITS_set.Rdata")

# ------------------------------------------------------------------------------
# stat: sensitive-related and resistant-related ITS ----------------------------
# ------------------------------------------------------------------------------
meta_result <- read.csv("05_immuneSig_ICBresponse/results/meta_results/res_metaRes.csv")
meta_selected <- meta_result[meta_result$Meta_Pval < 0.05, ]
meta_selected$immune_sig = meta_selected$X

load("05_immuneSig_ICBresponse/results_r0.4_ITSssgsea_CombatRef/finalmodel_17datasets_ssgseanorm/metaAnalysis_filter/en/elasticnet_model_pred.Rdata")

en_model_para = m_en$learner.model$glmnet.fit
m_en_imp <- data.frame(en_model_para$beta[,which(en_model_para$lambda == m_en$learner.model$lambda.min)])
names(m_en_imp) <- "weight"
m_en_imp$immune_sig = row.names(m_en_imp)
meta_selected = merge(meta_selected, m_en_imp, by = "immune_sig")
meta_selected = meta_selected[meta_selected$weight != 0, ]


meta_selected$flag = 'sensitive'
meta_selected[meta_selected$OR < 1, ]$flag = 'resistant'
ITS_S = meta_selected[meta_selected$flag == 'sensitive', ]$X
ITS_R = meta_selected[meta_selected$flag == 'resistant', ]$X

# ITS_S = ITS_S[ITS_S != "merck18_tide"]

# ------------------------------------------------------------------------------
# sensitive 

# ITS_p_s ------------------------------------------------------------------------------
ITS_p_s_notspecific_df <- read.table("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSp_sharedgene/ITS_p_s_notspecific_df.txt", header = T, sep = '\t')
# head(ITS_p_s_notspecific_df)



# enrich custom-GSEA -----------------------------------------------------------
load("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSp_sharedgene/enrich/custom_GSEA_s_scaled.Rdata")
res2 = data.frame(res)

# 1) drugtarget-Drug -------------------------------------------------------------------
filepath <- ("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSp_sharedgene/drug_LINCS/GSEAenrich/ITSp_s_vs_drugtarget")

glist <- unlist(strsplit(res2[res2$ID == "ITSp_vs_drugtarget",]$core_enrichment, "/"))
# DT_druglist <- unique(drugTargetMoA[is.element(drugTargetMoA$target, glist),]$drug_index)
# DT_drug <- drugTargetMoA[is.element(drugTargetMoA$drug_index, DT_druglist),]
DT_drug <- unique(drugTargetMoA[is.element(drugTargetMoA$target, glist),])

write.table(DT_drug, paste0(filepath, ".txt"), sep = '\t', quote = F, row.names = F)

# sort(table(unique(DT_drug[c(1,3)])$MoA))
# unique(DT_drug[DT_drug$MoA == "Arachidonate 5-lipoxygenase inhibitor",][c(1,2)])

# 2) OG-Drug -------------------------------------------------------------------
filepath <- ("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSp_sharedgene/drug_LINCS/GSEAenrich/ITSp_s_vs_OG_oncoKB")

glist <- unlist(strsplit(res2[res2$ID == "ITSp_vs_OG_oncoKB",]$core_enrichment, "/"))
# OG_druglist <- unique(drugTargetMoA[is.element(drugTargetMoA$target, glist),]$drug_index)
# OG_drug <- drugTargetMoA[is.element(drugTargetMoA$drug_index, OG_druglist),]
OG_drug <- unique(drugTargetMoA[is.element(drugTargetMoA$target, glist),])

write.table(OG_drug, paste0(filepath, ".txt"), sep = '\t', quote = F, row.names = F)
# sort(table(unique(OG_drug[c(1,3)])$MoA))
# unique(OG_drug[OG_drug$MoA == "Adenosine receptor antagonist",][c(1,2)])

# 3) TSG-Drug -------------------------------------------------------------------
filepath <- ("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSp_sharedgene/drug_LINCS/GSEAenrich/ITSp_s_vs_TSG_oncoKB")

glist <- unlist(strsplit(res2[res2$ID == "ITSp_vs_TSG_oncoKB",]$core_enrichment, "/"))
# TSG_druglist <- unique(drugTargetMoA[is.element(drugTargetMoA$target, glist),]$drug_index)
# TSG_drug <- drugTargetMoA[is.element(drugTargetMoA$drug_index, TSG_druglist),]
TSG_drug <- unique(drugTargetMoA[is.element(drugTargetMoA$target, glist),])
write.table(TSG_drug, paste0(filepath, ".txt"), sep = '\t', quote = F, row.names = F)

# sort(table(unique(TSG_drug[c(1,3)])$MoA))
# unique(TSG_drug[TSG_drug$MoA == "Histamine H1 receptor antagonist",][c(1,2)])


# enrich custom-ORA ---------------------------------------------------------------
dir.create("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSp_sharedgene/drug_LINCS/ORAenrich/")
filepath = "07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSp_sharedgene/enrich/"
load(paste0(filepath, "custom_ORA_s.Rdata"))

res2 = data.frame(res)

# 1) drugtarget-Drug -------------------------------------------------------------------
filepath <- ("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSp_sharedgene/drug_LINCS/ORAenrich/ITSp_vs_drugtarget")

glist <- unlist(strsplit(res2[res2$ID == "ITSp_vs_drugtarget",]$geneID, "/"))
glist <- ITS_p_s_notspecific_df[is.element(ITS_p_s_notspecific_df$ENTREZID, glist), ]

# DT_druglist <- unique(drugTargetMoA[is.element(drugTargetMoA$target, glist$gene_name),]$drug_index)
# DT_drug <- drugTargetMoA[is.element(drugTargetMoA$drug_index, DT_druglist),]
DT_drug <- unique(drugTargetMoA[is.element(drugTargetMoA$target, glist$gene_name),])

write.table(DT_drug, paste0(filepath, ".txt"), sep = '\t', quote = F, row.names = F)

names(DT_drug) = c(names(DT_drug)[1:3],"gene_name")
DT_drug2 <- merge(glist, DT_drug, by = 'gene_name',all=T)

# sort(table(unique(DT_drug[c(1,3)])$MoA))
# unique(DT_drug[DT_drug$MoA == "Histamine H1 receptor antagonist",][c(1,2)])

# head(DT_drug2[order(DT_drug2$weighted_rate,decreasing = T),])
# "Statin" https://bmccancer.biomedcentral.com/articles/10.1186/s12885-022-09385-8

# tmp=head(DT_drug2[order(DT_drug2$weighted_rate,decreasing = T),30)
# unique(tmp[order(tmp$drug_index),][c(1,7:9)])


# 2) OG-Drug -------------------------------------------------------------------
filepath <- ("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSp_sharedgene/drug_LINCS/ORAenrich/ITSp_s_vs_OG_oncoKB")

glist <- unlist(strsplit(res2[res2$ID == "ITSp_vs_OG_oncoKB",]$geneID, "/"))
glist <- ITS_p_s_notspecific_df[is.element(ITS_p_s_notspecific_df$ENTREZID, glist), ]

# OG_druglist <- unique(drugTargetMoA[is.element(drugTargetMoA$target, glist$gene_name),]$drug_index)
# OG_drug <- drugTargetMoA[is.element(drugTargetMoA$drug_index, OG_druglist),]
OG_drug <- unique(drugTargetMoA[is.element(drugTargetMoA$target, glist$gene_name),])

write.table(OG_drug, paste0(filepath, ".txt"), sep = '\t', quote = F, row.names = F)
# sort(table(unique(OG_drug[c(1,3)])$MoA))
# unique(OG_drug[OG_drug$MoA == "PI3-kinase class I inhibitor",][c(1,2)])



# 3) TSG-Drug -------------------------------------------------------------------
filepath <- ("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSp_sharedgene/drug_LINCS/ORAenrich/ITSp_s_vs_TSG_oncoKB")

glist <- unlist(strsplit(res2[res2$ID == "ITSp_vs_TSG_oncoKB",]$geneID, "/"))
glist <- ITS_p_s_notspecific_df[is.element(ITS_p_s_notspecific_df$ENTREZID, glist), ]

# TSG_druglist <- unique(drugTargetMoA[is.element(drugTargetMoA$target, glist$gene_name),]$drug_index)
# TSG_drug <- drugTargetMoA[is.element(drugTargetMoA$drug_index, TSG_druglist),]
TSG_drug <- unique(drugTargetMoA[is.element(drugTargetMoA$target, glist$gene_name),])

write.table(TSG_drug, paste0(filepath, ".txt"), sep = '\t', quote = F, row.names = F)
# sort(table(unique(TSG_drug[c(1,3)])$MoA))
# unique(TSG_drug[TSG_drug$MoA == "Ribonucleoside-diphosphate reductase RR1 inhibitor",][c(1,2)])




# ITS_p_s ------------------------------------------------------------------------------
ITS_p_r_notspecific_df <- read.table("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSp_sharedgene/ITS_p_r_notspecific_df.txt", header = T, sep = '\t')
# head(ITS_p_r_notspecific_df)



# enrich custom-GSEA -----------------------------------------------------------
load("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSp_sharedgene/enrich/custom_GSEA_r_scaled.Rdata")
res2 = data.frame(res)

# 1) drugtarget-Drug -------------------------------------------------------------------
filepath <- ("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSp_sharedgene/drug_LINCS/GSEAenrich/ITSp_r_vs_drugtarget")

glist <- unlist(strsplit(res2[res2$ID == "ITSp_vs_drugtarget",]$core_enrichment, "/"))
# DT_druglist <- unique(drugTargetMoA[is.element(drugTargetMoA$target, glist),]$drug_index)
# DT_drug <- drugTargetMoA[is.element(drugTargetMoA$drug_index, DT_druglist),]
DT_drug <- unique(drugTargetMoA[is.element(drugTargetMoA$target, glist),])
write.table(DT_drug, paste0(filepath, ".txt"), sep = '\t', quote = F, row.names = F)

# sort(table(unique(DT_drug[c(1,3)])$MoA))
# unique(DT_drug[DT_drug$MoA == "Tubulin inhibitor",][c(1,2)])

# 2) OG-Drug -------------------------------------------------------------------
filepath <- ("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSp_sharedgene/drug_LINCS/GSEAenrich/ITSp_r_vs_OG_oncoKB")

glist <- unlist(strsplit(res2[res2$ID == "ITSp_vs_OG_oncoKB",]$core_enrichment, "/"))
# OG_druglist <- unique(drugTargetMoA[is.element(drugTargetMoA$target, glist),]$drug_index)
# OG_drug <- drugTargetMoA[is.element(drugTargetMoA$drug_index, OG_druglist),]
OG_drug <- unique(drugTargetMoA[is.element(drugTargetMoA$target, glist),])
write.table(OG_drug, paste0(filepath, ".txt"), sep = '\t', quote = F, row.names = F)
# sort(table(unique(OG_drug[c(1,3)])$MoA))
# unique(OG_drug[OG_drug$MoA == "Adenosine receptor antagonist",][c(1,2)])

# 3) TSG-Drug -------------------------------------------------------------------
filepath <- ("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSp_sharedgene/drug_LINCS/GSEAenrich/ITSp_r_vs_TSG_oncoKB")

glist <- unlist(strsplit(res2[res2$ID == "ITSp_vs_TSG_oncoKB",]$core_enrichment, "/"))
# TSG_druglist <- unique(drugTargetMoA[is.element(drugTargetMoA$target, glist),]$drug_index)
# TSG_drug <- drugTargetMoA[is.element(drugTargetMoA$drug_index, TSG_druglist),]
TSG_drug <- unique(drugTargetMoA[is.element(drugTargetMoA$target, glist),])
write.table(TSG_drug, paste0(filepath, ".txt"), sep = '\t', quote = F, row.names = F)
# sort(table(unique(TSG_drug[c(1,3)])$MoA))
# unique(TSG_drug[TSG_drug$MoA == "Histamine H1 receptor antagonist",][c(1,2)])


# enrich custom-ORA ---------------------------------------------------------------
dir.create("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSp_sharedgene/drug_LINCS/ORAenrich/")
filepath = "07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSp_sharedgene/enrich/"
load(paste0(filepath, "custom_ORA_s.Rdata"))

res2 = data.frame(res)

# 1) drugtarget-Drug -------------------------------------------------------------------
filepath <- ("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSp_sharedgene/drug_LINCS/ORAenrich/ITSp_vs_drugtarget")

glist <- unlist(strsplit(res2[res2$ID == "ITSp_vs_drugtarget",]$geneID, "/"))
glist <- ITS_p_r_notspecific_df[is.element(ITS_p_r_notspecific_df$ENTREZID, glist), ]

# DT_druglist <- unique(drugTargetMoA[is.element(drugTargetMoA$target, glist$gene_name),]$drug_index)
# DT_drug <- drugTargetMoA[is.element(drugTargetMoA$drug_index, DT_druglist),]
DT_drug <- unique(drugTargetMoA[is.element(drugTargetMoA$target, glist$gene_name),])

write.table(DT_drug, paste0(filepath, ".txt"), sep = '\t', quote = F, row.names = F)

names(DT_drug) = c(names(DT_drug)[1:3],"gene_name")
DT_drug2 <- merge(glist, DT_drug, by = 'gene_name',all=T)

# sort(table(unique(DT_drug[c(1,3)])$MoA))
# unique(DT_drug[DT_drug$MoA == "Histamine H1 receptor antagonist",][c(1,2)])

# head(DT_drug2[order(DT_drug2$weighted_rate,decreasing = T),])
# "Statin" https://bmccancer.biomedcentral.com/articles/10.1186/s12885-022-09385-8

# tmp=head(DT_drug2[order(DT_drug2$weighted_rate,decreasing = T),30)
# unique(tmp[order(tmp$drug_index),][c(1,7:9)])


# 2) OG-Drug -------------------------------------------------------------------
filepath <- ("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSp_sharedgene/drug_LINCS/ORAenrich/ITSp_r_vs_OG_oncoKB")

glist <- unlist(strsplit(res2[res2$ID == "ITSp_vs_OG_oncoKB",]$geneID, "/"))
glist <- ITS_p_r_notspecific_df[is.element(ITS_p_r_notspecific_df$ENTREZID, glist), ]

# OG_druglist <- unique(drugTargetMoA[is.element(drugTargetMoA$target, glist$gene_name),]$drug_index)
# OG_drug <- drugTargetMoA[is.element(drugTargetMoA$drug_index, OG_druglist),]
OG_drug <- unique(drugTargetMoA[is.element(drugTargetMoA$target, glist$gene_name),])

write.table(OG_drug, paste0(filepath, ".txt"), sep = '\t', quote = F, row.names = F)
# sort(table(unique(OG_drug[c(1,3)])$MoA))
# unique(OG_drug[OG_drug$MoA == "PI3-kinase class I inhibitor",][c(1,2)])



# 3) TSG-Drug -------------------------------------------------------------------
filepath <- ("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSp_sharedgene/drug_LINCS/ORAenrich/ITSp_r_vs_TSG_oncoKB")

glist <- unlist(strsplit(res2[res2$ID == "ITSp_vs_TSG_oncoKB",]$geneID, "/"))
glist <- ITS_p_r_notspecific_df[is.element(ITS_p_r_notspecific_df$ENTREZID, glist), ]

# TSG_druglist <- unique(drugTargetMoA[is.element(drugTargetMoA$target, glist$gene_name),]$drug_index)
# TSG_drug <- drugTargetMoA[is.element(drugTargetMoA$drug_index, TSG_druglist),]
TSG_drug <- unique(drugTargetMoA[is.element(drugTargetMoA$target, glist$gene_name),])

write.table(TSG_drug, paste0(filepath, ".txt"), sep = '\t', quote = F, row.names = F)
# sort(table(unique(TSG_drug[c(1,3)])$MoA))
# unique(TSG_drug[TSG_drug$MoA == "Ribonucleoside-diphosphate reductase RR1 inhibitor",][c(1,2)])




# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ITS_n   ----------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
dir.create("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSn_sharedgene/")
dir.create("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSn_sharedgene/drug_LINCS/")
dir.create("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSn_sharedgene/drug_LINCS/GSEAenrich/")


# ITS_n_s ------------------------------------------------------------------------------
ITS_n_s_notspecific_df <- read.table("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSn_sharedgene/ITS_n_s_notspecific_df.txt", header = T, sep = '\t')
# head(ITS_n_s_notspecific_df)



# enrich custom-GSEA -----------------------------------------------------------
load("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSn_sharedgene/enrich/custom_GSEA_s_scaled.Rdata")
res2 = data.frame(res)

# 1) drugtarget-Drug -------------------------------------------------------------------
filepath <- ("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSn_sharedgene/drug_LINCS/GSEAenrich/ITSn_s_vs_drugtarget")

glist <- unlist(strsplit(res2[res2$ID == "ITSp_vs_drugtarget",]$core_enrichment, "/"))
# DT_druglist <- unique(drugTargetMoA[is.element(drugTargetMoA$target, glist),]$drug_index)
# DT_drug <- drugTargetMoA[is.element(drugTargetMoA$drug_index, DT_druglist),]
DT_drug <- unique(drugTargetMoA[is.element(drugTargetMoA$target, glist),])
write.table(DT_drug, paste0(filepath, ".txt"), sep = '\t', quote = F, row.names = F)

# sort(table(unique(DT_drug[c(1,3)])$MoA))
# unique(DT_drug[DT_drug$MoA == "Arachidonate 5-lipoxygenase inhibitor",][c(1,2)])

# 2) OG-Drug -------------------------------------------------------------------
filepath <- ("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSn_sharedgene/drug_LINCS/GSEAenrich/ITSn_s_vs_OG_oncoKB")

glist <- unlist(strsplit(res2[res2$ID == "ITSp_vs_OG_oncoKB",]$core_enrichment, "/"))
# OG_druglist <- unique(drugTargetMoA[is.element(drugTargetMoA$target, glist),]$drug_index)
# OG_drug <- drugTargetMoA[is.element(drugTargetMoA$drug_index, OG_druglist),]
OG_drug <- unique(drugTargetMoA[is.element(drugTargetMoA$target, glist),])
write.table(OG_drug, paste0(filepath, ".txt"), sep = '\t', quote = F, row.names = F)
# sort(table(unique(OG_drug[c(1,3)])$MoA))
# unique(OG_drug[OG_drug$MoA == "Vascular endothelial growth factor receptor inhibitor",][c(1,2)])

# 3) TSG-Drug -------------------------------------------------------------------
filepath <- ("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSn_sharedgene/drug_LINCS/GSEAenrich/ITSn_s_vs_TSG_oncoKB")

glist <- unlist(strsplit(res2[res2$ID == "ITSp_vs_TSG_oncoKB",]$core_enrichment, "/"))
# TSG_druglist <- unique(drugTargetMoA[is.element(drugTargetMoA$target, glist),]$drug_index)
# TSG_drug <- drugTargetMoA[is.element(drugTargetMoA$drug_index, TSG_druglist),]
TSG_drug <- unique(drugTargetMoA[is.element(drugTargetMoA$target, glist),])
write.table(TSG_drug, paste0(filepath, ".txt"), sep = '\t', quote = F, row.names = F)
# sort(table(unique(TSG_drug[c(1,3)])$MoA))
# unique(TSG_drug[TSG_drug$MoA == "PI3-kinase class I inhibitor",][c(1,2)])


# # enrich custom-ORA ---------------------------------------------------------------
# dir.create("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSn_sharedgene/drug_LINCS/ORAenrich/")
# filepath = "07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSn_sharedgene/enrich/"
# load(paste0(filepath, "custom_ORA_s.Rdata"))

# res2 = data.frame(res)

# # 1) drugtarget-Drug -------------------------------------------------------------------
# filepath <- ("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSn_sharedgene/drug_LINCS/ORAenrich/ITSp_vs_drugtarget")

# glist <- unlist(strsplit(res2[res2$ID == "ITSp_vs_drugtarget",]$geneID, "/"))
# glist <- ITS_n_s_notspecific_df[is.element(ITS_n_s_notspecific_df$ENTREZID, glist), ]

# DT_druglist <- unique(drugTargetMoA[is.element(drugTargetMoA$target, glist$gene_name),]$drug_index)
# DT_drug <- drugTargetMoA[is.element(drugTargetMoA$drug_index, DT_druglist),]
# write.table(DT_drug, paste0(filepath, ".txt"), sep = '\t', quote = F, row.names = F)

# names(DT_drug) = c(names(DT_drug)[1:3],"gene_name")
# DT_drug2 <- merge(glist, DT_drug, by = 'gene_name',all=T)

# # sort(table(unique(DT_drug[c(1,3)])$MoA))
# # unique(DT_drug[DT_drug$MoA == "Histamine H1 receptor antagonist",][c(1,2)])

# # head(DT_drug2[order(DT_drug2$weighted_rate,decreasing = T),])
# # "Statin" https://bmccancer.biomedcentral.com/articles/10.1186/s12885-022-09385-8

# # tmp=head(DT_drug2[order(DT_drug2$weighted_rate,decreasing = T),30)
# # unique(tmp[order(tmp$drug_index),][c(1,7:9)])


# # 2) OG-Drug -------------------------------------------------------------------
# filepath <- ("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSn_sharedgene/drug_LINCS/ORAenrich/ITSn_s_vs_OG_oncoKB")

# glist <- unlist(strsplit(res2[res2$ID == "ITSp_vs_OG_oncoKB",]$geneID, "/"))
# glist <- ITS_n_s_notspecific_df[is.element(ITS_n_s_notspecific_df$ENTREZID, glist), ]

# OG_druglist <- unique(drugTargetMoA[is.element(drugTargetMoA$target, glist$gene_name),]$drug_index)
# OG_drug <- drugTargetMoA[is.element(drugTargetMoA$drug_index, OG_druglist),]
# write.table(OG_drug, paste0(filepath, ".txt"), sep = '\t', quote = F, row.names = F)
# # sort(table(unique(OG_drug[c(1,3)])$MoA))
# # unique(OG_drug[OG_drug$MoA == "PI3-kinase class I inhibitor",][c(1,2)])



# # 3) TSG-Drug -------------------------------------------------------------------
# filepath <- ("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSn_sharedgene/drug_LINCS/ORAenrich/ITSn_s_vs_TSG_oncoKB")

# glist <- unlist(strsplit(res2[res2$ID == "ITSp_vs_TSG_oncoKB",]$geneID, "/"))
# glist <- ITS_n_s_notspecific_df[is.element(ITS_n_s_notspecific_df$ENTREZID, glist), ]

# TSG_druglist <- unique(drugTargetMoA[is.element(drugTargetMoA$target, glist$gene_name),]$drug_index)
# TSG_drug <- drugTargetMoA[is.element(drugTargetMoA$drug_index, TSG_druglist),]
# write.table(TSG_drug, paste0(filepath, ".txt"), sep = '\t', quote = F, row.names = F)
# # sort(table(unique(TSG_drug[c(1,3)])$MoA))
# # unique(TSG_drug[TSG_drug$MoA == "Ribonucleoside-diphosphate reductase RR1 inhibitor",][c(1,2)])




# ITS_n_s ------------------------------------------------------------------------------
ITS_n_r_notspecific_df <- read.table("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSn_sharedgene/ITS_n_r_notspecific_df.txt", header = T, sep = '\t')
# head(ITS_n_r_notspecific_df)



# enrich custom-GSEA -----------------------------------------------------------
load("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSn_sharedgene/enrich/custom_GSEA_r_scaled.Rdata")
res2 = data.frame(res)

# 1) drugtarget-Drug -------------------------------------------------------------------
filepath <- ("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSn_sharedgene/drug_LINCS/GSEAenrich/ITSn_r_vs_drugtarget")

glist <- unlist(strsplit(res2[res2$ID == "ITSp_vs_drugtarget",]$core_enrichment, "/"))
# DT_druglist <- unique(drugTargetMoA[is.element(drugTargetMoA$target, glist),]$drug_index)
# DT_drug <- drugTargetMoA[is.element(drugTargetMoA$drug_index, DT_druglist),]
DT_drug <- unique(drugTargetMoA[is.element(drugTargetMoA$target, glist),])
write.table(DT_drug, paste0(filepath, ".txt"), sep = '\t', quote = F, row.names = F)

# sort(table(unique(DT_drug[c(1,3)])$MoA))
# unique(DT_drug[DT_drug$MoA == "Tubulin inhibitor",][c(1,2)])

# 2) OG-Drug -------------------------------------------------------------------
filepath <- ("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSn_sharedgene/drug_LINCS/GSEAenrich/ITSn_r_vs_OG_oncoKB")

glist <- unlist(strsplit(res2[res2$ID == "ITSp_vs_OG_oncoKB",]$core_enrichment, "/"))
# OG_druglist <- unique(drugTargetMoA[is.element(drugTargetMoA$target, glist),]$drug_index)
# OG_drug <- drugTargetMoA[is.element(drugTargetMoA$drug_index, OG_druglist),]
OG_drug <- unique(drugTargetMoA[is.element(drugTargetMoA$target, glist),])
write.table(OG_drug, paste0(filepath, ".txt"), sep = '\t', quote = F, row.names = F)
# sort(table(unique(OG_drug[c(1,3)])$MoA))
# unique(OG_drug[OG_drug$MoA == "Adenosine receptor antagonist",][c(1,2)])

# 3) TSG-Drug -------------------------------------------------------------------
filepath <- ("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSn_sharedgene/drug_LINCS/GSEAenrich/ITSn_r_vs_TSG_oncoKB")

glist <- unlist(strsplit(res2[res2$ID == "ITSp_vs_TSG_oncoKB",]$core_enrichment, "/"))
# TSG_druglist <- unique(drugTargetMoA[is.element(drugTargetMoA$target, glist),]$drug_index)
# TSG_drug <- drugTargetMoA[is.element(drugTargetMoA$drug_index, TSG_druglist),]
TSG_drug <- unique(drugTargetMoA[is.element(drugTargetMoA$target, glist),])
write.table(TSG_drug, paste0(filepath, ".txt"), sep = '\t', quote = F, row.names = F)
# sort(table(unique(TSG_drug[c(1,3)])$MoA))
# unique(TSG_drug[TSG_drug$MoA == "Histamine H1 receptor antagonist",][c(1,2)])


# # enrich custom-ORA ---------------------------------------------------------------
# dir.create("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSn_sharedgene/drug_LINCS/ORAenrich/")
# filepath = "07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSn_sharedgene/enrich/"
# load(paste0(filepath, "custom_ORA_s.Rdata"))

# res2 = data.frame(res)

# # 1) drugtarget-Drug -------------------------------------------------------------------
# filepath <- ("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSn_sharedgene/drug_LINCS/ORAenrich/ITSp_vs_drugtarget")

# glist <- unlist(strsplit(res2[res2$ID == "ITSp_vs_drugtarget",]$geneID, "/"))
# glist <- ITS_n_r_notspecific_df[is.element(ITS_n_r_notspecific_df$ENTREZID, glist), ]

# DT_druglist <- unique(drugTargetMoA[is.element(drugTargetMoA$target, glist$gene_name),]$drug_index)
# DT_drug <- drugTargetMoA[is.element(drugTargetMoA$drug_index, DT_druglist),]
# write.table(DT_drug, paste0(filepath, ".txt"), sep = '\t', quote = F, row.names = F)

# names(DT_drug) = c(names(DT_drug)[1:3],"gene_name")
# DT_drug2 <- merge(glist, DT_drug, by = 'gene_name',all=T)

# # sort(table(unique(DT_drug[c(1,3)])$MoA))
# # unique(DT_drug[DT_drug$MoA == "Histamine H1 receptor antagonist",][c(1,2)])

# # head(DT_drug2[order(DT_drug2$weighted_rate,decreasing = T),])
# # "Statin" https://bmccancer.biomedcentral.com/articles/10.1186/s12885-022-09385-8

# # tmp=head(DT_drug2[order(DT_drug2$weighted_rate,decreasing = T),30)
# # unique(tmp[order(tmp$drug_index),][c(1,7:9)])


# # 2) OG-Drug -------------------------------------------------------------------
# filepath <- ("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSn_sharedgene/drug_LINCS/ORAenrich/ITSn_r_vs_OG_oncoKB")

# glist <- unlist(strsplit(res2[res2$ID == "ITSp_vs_OG_oncoKB",]$geneID, "/"))
# glist <- ITS_n_r_notspecific_df[is.element(ITS_n_r_notspecific_df$ENTREZID, glist), ]

# OG_druglist <- unique(drugTargetMoA[is.element(drugTargetMoA$target, glist$gene_name),]$drug_index)
# OG_drug <- drugTargetMoA[is.element(drugTargetMoA$drug_index, OG_druglist),]
# write.table(OG_drug, paste0(filepath, ".txt"), sep = '\t', quote = F, row.names = F)
# # sort(table(unique(OG_drug[c(1,3)])$MoA))
# # unique(OG_drug[OG_drug$MoA == "PI3-kinase class I inhibitor",][c(1,2)])



# # 3) TSG-Drug -------------------------------------------------------------------
# filepath <- ("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSn_sharedgene/drug_LINCS/ORAenrich/ITSn_r_vs_TSG_oncoKB")

# glist <- unlist(strsplit(res2[res2$ID == "ITSp_vs_TSG_oncoKB",]$geneID, "/"))
# glist <- ITS_n_r_notspecific_df[is.element(ITS_n_r_notspecific_df$ENTREZID, glist), ]

# TSG_druglist <- unique(drugTargetMoA[is.element(drugTargetMoA$target, glist$gene_name),]$drug_index)
# TSG_drug <- drugTargetMoA[is.element(drugTargetMoA$drug_index, TSG_druglist),]
# write.table(TSG_drug, paste0(filepath, ".txt"), sep = '\t', quote = F, row.names = F)
# # sort(table(unique(TSG_drug[c(1,3)])$MoA))
# # unique(TSG_drug[TSG_drug$MoA == "Ribonucleoside-diphosphate reductase RR1 inhibitor",][c(1,2)])

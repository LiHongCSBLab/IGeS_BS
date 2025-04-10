# rm(list=ls())
# Mapping drug DEGs to correlation results -------------------------------------
# work_path <- "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/"
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
 
source("06_results_summaryandplotting/scripts/immunesig_selection_function.R")
source("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/scripts/drugITGSscore_lincs.R")
source("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/scripts/drugITGSscore_GEOdataset.R")
source("06_results_summaryandplotting/scripts/plot_for_drug_function_new.R")

# ----------------------------------------------------------------------------
# drug-induced immune signatures NES
cancer = "SKCM"
tumor_subtype = "Metastatic"
purity_method = "TUMERIC"
dataset=70138
datatype="allgenes"
num_gene = 200


cor_method = "spearman"
cluster_method = "mclust"
similarity_method = "Otsuka_Ochiai" # "jaccard_index",
num_gene = 200
s1 = 0.4
ITSproof = 2
genevote = 0.5
model = "m2"


# ITSselected <- read.csv("06_results_summaryandplotting/data/immune_sig_selected_new/selected_ITS_0.4.txt", sep = '\t')
meta_res <- read.csv("05_immuneSig_ICBresponse/results/meta_results/res_metaRes.csv")
meta_res_selected <- meta_res[meta_res$Meta_Pval < 0.05, ]
sig_s = meta_res_selected[meta_res_selected$OR > 1,]
sig_r = meta_res_selected[meta_res_selected$OR < 1,]
sig_s = data.frame( immune_sig = sig_s$X, flag = "sensitive", setindex = '01_sensitive')
sig_r = data.frame( immune_sig = sig_r$X, flag = "resistant", setindex = '02_resistant')

load("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/05_immuneSig_ICBresponse/results_r0.4_ITSssgsea_CombatRef/finalmodel_17datasets_ssgseanorm/metaAnalysis_filter/en/elasticnet_model_pred.Rdata")

en_model_para = m_en$learner.model$glmnet.fit
m_en_imp <- data.frame(en_model_para$beta[,which(en_model_para$lambda == m_en$learner.model$lambda.min)])

names(m_en_imp) <- "weight"
m_en_imp$immune_sig = row.names(m_en_imp)


immunesig = rbind(sig_s, sig_r)
immunesig = merge(immunesig, m_en_imp, by = "immune_sig")
immunesig =immunesig[immunesig$weight != 0 ,]

# 1) if drug effect on cell lines and mice models are consistant -------------------------
# ------------------------------------------------------------------------------
setwd("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/")
options(stringsAsFactors = F)
GEO_ID = "GSE149825"
drugindex = "drug1487" 
s1=0.4
ITSproof=2
genevote = 0.5
cluster_method="mclust"
similarity_method="Otsuka_Ochiai" # jaccard_index

# gsea_resultpath = paste0("result_", GEO_ID, "/GSEA_result_",cluster_method,"_",similarity_method,"_", purity_method,"_r",s1,"_meta_oe_original/", 
#                          dataset, "/", cancer, "_", tumor_subtype, "/")
gsea_resultpath = paste0("result_", GEO_ID, "/drug_immunesig_combo_", purity_method,"_GSEA_result_r",s1,"_meta_oe_original/")


Pval=1
dir.create(paste0("result_", GEO_ID))
dir.create(paste0("result_", GEO_ID, "/drug_immunesig_combo_", purity_method,"_GSEA_result_r",s1,"_meta_oe_original_en_padj",Pval,"_final/"))
savepath = paste0("result_", GEO_ID, "/drug_immunesig_combo_", purity_method,"_GSEA_result_r",s1,"_meta_oe_original_en_padj",Pval,"_final/", drugindex)
dir.create(savepath)


drugITGSscore_GEO(dataset = dataset,
                  datatype = datatype,
                  cancer = cancer,
                  tumor_subtype = tumor_subtype,
                  purity_method = purity_method, 
                  immunesig = immunesig,
                  drugindex = drugindex,
                  type = "rawITS",
                  Pval = Pval,
                  gsea_resultpath = gsea_resultpath,
                  savepath = savepath)



drugITGSscore_GEO_tissueVScell(dataset = dataset,
                  datatype = datatype,
                  cancer = cancer,
                  tumor_subtype = tumor_subtype,
                  purity_method = purity_method, 
                  immunesig = immunesig,
                  drugindex = drugindex,
                  gsea_resultpath = gsea_resultpath,
                  compare_type ="rawITS",
                  Pval = Pval,
                  savepath = savepath)


gsea_resultpath = paste0("result_", GEO_ID, "/drug_immunesig_combo_", purity_method,"_GSEA_result_r",s1,"_meta_oe_original/")


Pval=0.1
dir.create(paste0("result_", GEO_ID, "/drug_immunesig_combo_", purity_method,"_GSEA_result_r",s1,"_meta_oe_original_en_padj",Pval,"_final/"))
savepath = paste0("result_", GEO_ID, "/drug_immunesig_combo_", purity_method,"_GSEA_result_r",s1,"_meta_oe_original_en_padj",Pval,"_final/", drugindex)
dir.create(savepath)

drugITGSscore_GEO(dataset = dataset,
                  datatype = datatype,
                  cancer = cancer,
                  tumor_subtype = tumor_subtype,
                  purity_method = purity_method, 
                  immunesig = immunesig,
                  drugindex = drugindex,
                  type = "rawITS",
                  Pval = Pval,
                  gsea_resultpath = gsea_resultpath,
                  savepath = savepath)



drugITGSscore_GEO_tissueVScell(dataset = dataset,
                  datatype = datatype,
                  cancer = cancer,
                  tumor_subtype = tumor_subtype,
                  purity_method = purity_method, 
                  immunesig = immunesig,
                  drugindex = drugindex,
                  gsea_resultpath = gsea_resultpath,
                  compare_type ="rawITS",
                  Pval = Pval,
                  savepath = savepath)


Pval=0.05
dir.create(paste0("result_", GEO_ID, "/drug_immunesig_combo_", purity_method,"_GSEA_result_r",s1,"_meta_oe_original_en_padj",Pval,"_final/"))
savepath = paste0("result_", GEO_ID, "/drug_immunesig_combo_", purity_method,"_GSEA_result_r",s1,"_meta_oe_original_en_padj",Pval,"_final/", drugindex)
dir.create(savepath)

drugITGSscore_GEO(dataset = dataset,
                  datatype = datatype,
                  cancer = cancer,
                  tumor_subtype = tumor_subtype,
                  purity_method = purity_method, 
                  immunesig = immunesig,
                  drugindex = drugindex,
                  type = "rawITS",
                  Pval = Pval,
                  gsea_resultpath = gsea_resultpath,
                  savepath = savepath)



drugITGSscore_GEO_tissueVScell(dataset = dataset,
                  datatype = datatype,
                  cancer = cancer,
                  tumor_subtype = tumor_subtype,
                  purity_method = purity_method, 
                  immunesig = immunesig,
                  drugindex = drugindex,
                  gsea_resultpath = gsea_resultpath,
                  compare_type ="rawITS",
                  Pval = Pval,
                  savepath = savepath)

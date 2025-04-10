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
cancer = "COAD"
tumor_subtype = "Primary"
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


savepath = paste0("06_results_summaryandplotting/results_meta_oe_original/", purity_method,"_", model,"_new/")
# savepath = paste0("06_results_summaryandplotting/results_meta_oe_original/", purity_method,"/")
dir.create(savepath)


immunesig <- read.csv("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/05_immuneSig_ICBresponse/results/meta_results/res_metaRes.csv")
immunesig <- immunesig[immunesig$Meta_Pval < 0.05, ]
names(immunesig) <- c("immune_sig", names(immunesig)[-1])
immunesig$flag <- "sensitive"
immunesig[immunesig$OR < 1,]$flag <- "resistant"
immunesig$setindex = "01_sensitive"
immunesig[immunesig$flag == "resistant", ]$setindex = "02_resistant"
# immunesig$weight = log2(immunesig$OR)


savepath = paste0("06_results_summaryandplotting/results_meta_oe_original/", purity_method,"_",model,"_meta_weight_ACATP0.05_en/") #
dir.create(savepath)

# immunesig_ITS_index <- read.table(paste0("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/03_drug_immuneSig_enrichment/data/ITS_mclust_Otsuka_Ochiai_TUMERIC_spearman_r_0.4_metaP_0.05_genevote_0.5/",cancer,"_immunesig_ITS_index.txt"), sep = "\t", header = T)
# immunesig <- inner_join(immunesig_ITS_index[-3], immunesig)
load("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/05_immuneSig_ICBresponse/results_r0.4_ITSssgsea_CombatRef/finalmodel_17datasets_ssgseanorm/metaAnalysis_filter/en/elasticnet_model_pred.Rdata")

en_model_para = m_en$learner.model$glmnet.fit
m_en_imp <- data.frame(en_model_para$beta[,which(en_model_para$lambda == m_en$learner.model$lambda.min)])
names(m_en_imp) <- "weight"
m_en_imp$immune_sig = row.names(m_en_imp)
immunesig = merge(immunesig, m_en_imp, by = "immune_sig")
immunesig =immunesig[immunesig$weight != 0 ,]

dir.create(paste0(savepath, "/"))
dir.create(paste0(savepath, "/",dataset,"/"))
dir.create(paste0(savepath, "/",dataset,"/", datatype,"/"))
dir.create(paste0(savepath, "/",dataset,"/", datatype,"/", cancer, "_", tumor_subtype))
  
  
df_immunesig <- drugITSscore(cancer = cancer,
                                tumor_subtype = tumor_subtype,
                                purity_method = purity_method,
                                dataset=dataset,
                                datatype=datatype,
                                num_gene = num_gene,
                                immunesig = immunesig, 
                                model = model,
                                ACAT_Pval = 0.05,
                                work_path = work_path,
                                savepath = savepath)


savepath2 = "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/data/drugITS_meta_oe_original/"
dir.create(savepath2)
drugITGSscore_preparation(df_immunesig = df_immunesig,
                          cancer = cancer,
                          tumor_subtype = tumor_subtype,
                          purity_method = purity_method,
                          dataset=dataset,
                          datatype="allgenes",
                          num_gene = 200,
                          work_path = work_path,
                          save_path = savepath2)

# ITSselected <- read.csv("06_results_summaryandplotting/data/immune_sig_selected_new/selected_ITS_0.4.txt", sep = '\t')
meta_res <- read.csv("05_immuneSig_ICBresponse/results/meta_results/res_metaRes.csv")
meta_res_selected <- meta_res[meta_res$Meta_Pval < 0.05, ]
sig_s = meta_res_selected[meta_res_selected$OR > 1,]
sig_r = meta_res_selected[meta_res_selected$OR < 1,]
sig_s = data.frame( immune_sig = sig_s$X, flag = "sensitive", setindex = '01_sensitive')
sig_r = data.frame( immune_sig = sig_r$X, flag = "resistant", setindex = '02_resistant')

# mergedP_res <- read.csv("05_immuneSig_ICBresponse/results/meta_results/res_merged_p.csv")
# sig_s = mergedP_res[mergedP_res$merged_RgreaterNR_p_ACAT < 0.05,]
# sig_r = mergedP_res[mergedP_res$merged_RlessNR_p_ACAT < 0.05,]
# sig_s = data.frame( immune_sig = intersect(sig_s$X, meta_res_selected$X), flag = "sensitive", setindex = '01_sensitive')
# sig_r = data.frame( immune_sig = intersect(sig_r$X, meta_res_selected$X), flag = "resistant", setindex = '02_resistant')


load("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/05_immuneSig_ICBresponse/results_r0.4_ITSssgsea_CombatRef/finalmodel_17datasets_ssgseanorm/metaAnalysis_filter/en/elasticnet_model_pred.Rdata")

en_model_para = m_en$learner.model$glmnet.fit
m_en_imp <- data.frame(en_model_para$beta[,which(en_model_para$lambda == m_en$learner.model$lambda.min)])

names(m_en_imp) <- "weight"
m_en_imp$immune_sig = row.names(m_en_imp)


immunesig = rbind(sig_s, sig_r)
immunesig = merge(immunesig, m_en_imp, by = "immune_sig")
immunesig =immunesig[immunesig$weight != 0 ,]
                          


# celecoxib: drug989 - with COAD drugITGS --------------------------------------
setwd("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/")
options(stringsAsFactors = F)
GEO_ID = "GSE160785"
drugindex = "drug989"


gsea_resultpath = paste0("result_", GEO_ID, "/drug_immunesig_", purity_method,"_GSEA_result_r",s1,"_meta_oe_original/", cancer, "_", tumor_subtype, "/")


# ITSproof=2
dir.create(paste0("result_", GEO_ID))
Pval = 1
dir.create(paste0("result_", GEO_ID, "/drug_immunesig_", purity_method,"_GSEA_result_r",s1,"_meta_oe_original_en_padj",Pval,"_final/"))
dir.create(paste0("result_", GEO_ID, "/drug_immunesig_", purity_method,"_GSEA_result_r",s1,"_meta_oe_original_en_padj",Pval,"_final/", cancer, "_", tumor_subtype, "/"))
savepath = paste0("result_", GEO_ID, "/drug_immunesig_", purity_method,"_GSEA_result_r",s1,"_meta_oe_original_en_padj",Pval,"_final/", cancer, "_", tumor_subtype, "/", drugindex)
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
                  compare_type ="rawITS",
                  Pval = Pval,
                  gsea_resultpath = gsea_resultpath,
                  savepath = savepath)


Pval = 0.1
dir.create(paste0("result_", GEO_ID, "/drug_immunesig_", purity_method,"_GSEA_result_r",s1,"_meta_oe_original_en_padj",Pval,"_final/"))
dir.create(paste0("result_", GEO_ID, "/drug_immunesig_", purity_method,"_GSEA_result_r",s1,"_meta_oe_original_en_padj",Pval,"_final/", cancer, "_", tumor_subtype, "/"))
savepath = paste0("result_", GEO_ID, "/drug_immunesig_", purity_method,"_GSEA_result_r",s1,"_meta_oe_original_en_padj",Pval,"_final/", cancer, "_", tumor_subtype, "/", drugindex)
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
                  compare_type ="rawITS",
                  Pval = Pval,
                  gsea_resultpath = gsea_resultpath,
                  savepath = savepath)

Pval = 0.05
dir.create(paste0("result_", GEO_ID, "/drug_immunesig_", purity_method,"_GSEA_result_r",s1,"_meta_oe_original_en_padj",Pval,"_final/"))
dir.create(paste0("result_", GEO_ID, "/drug_immunesig_", purity_method,"_GSEA_result_r",s1,"_meta_oe_original_en_padj",Pval,"_final/", cancer, "_", tumor_subtype, "/"))
savepath = paste0("result_", GEO_ID, "/drug_immunesig_", purity_method,"_GSEA_result_r",s1,"_meta_oe_original_en_padj",Pval,"_final/", cancer, "_", tumor_subtype, "/", drugindex)
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
                  compare_type ="rawITS",
                  Pval = Pval,
                  gsea_resultpath = gsea_resultpath,
                  savepath = savepath)



# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

# Mapping drug DEGs to correlation results -------------------------------------
# work_path <- "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/"
work_path <- "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/"

setwd(work_path)
options(stringsAsFactors = F)

cancer = "READ"
tumor_subtype = "Primary"
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


savepath = paste0("06_results_summaryandplotting/results_meta_oe_original/", purity_method,"_", model,"_new/")
# savepath = paste0("06_results_summaryandplotting/results_meta_oe_original/", purity_method,"/")
dir.create(savepath)



immunesig <- read.csv("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/05_immuneSig_ICBresponse/results/meta_results/res_metaRes.csv")
immunesig <- immunesig[immunesig$Meta_Pval < 0.05, ]
names(immunesig) <- c("immune_sig", names(immunesig)[-1])
immunesig$flag <- "sensitive"
immunesig[immunesig$OR < 1,]$flag <- "resistant"
immunesig$setindex = "01_sensitive"
immunesig[immunesig$flag == "resistant", ]$setindex = "02_resistant"
# immunesig$weight = log2(immunesig$OR)


savepath = paste0("06_results_summaryandplotting/results_meta_oe_original/", purity_method,"_",model,"_meta_weight_ACATP0.05_en/") #
dir.create(savepath)

# immunesig_ITS_index <- read.table(paste0("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/03_drug_immuneSig_enrichment/data/ITS_mclust_Otsuka_Ochiai_TUMERIC_spearman_r_0.4_metaP_0.05_genevote_0.5/",cancer,"_immunesig_ITS_index.txt"), sep = "\t", header = T)
# immunesig <- inner_join(immunesig_ITS_index[-3], immunesig)
load("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/05_immuneSig_ICBresponse/results_r0.4_ITSssgsea_CombatRef/finalmodel_17datasets_ssgseanorm/metaAnalysis_filter/en/elasticnet_model_pred.Rdata")

en_model_para = m_en$learner.model$glmnet.fit
m_en_imp <- data.frame(en_model_para$beta[,which(en_model_para$lambda == m_en$learner.model$lambda.min)])
names(m_en_imp) <- "weight"
m_en_imp$immune_sig = row.names(m_en_imp)
immunesig = merge(immunesig, m_en_imp, by = "immune_sig")
immunesig =immunesig[immunesig$weight != 0 ,]

dir.create(paste0(savepath, "/"))
dir.create(paste0(savepath, "/",dataset,"/"))
dir.create(paste0(savepath, "/",dataset,"/", datatype,"/"))
dir.create(paste0(savepath, "/",dataset,"/", datatype,"/", cancer, "_", tumor_subtype))
  
  
  
df_immunesig <- drugITSscore(cancer = cancer,
                                tumor_subtype = tumor_subtype,
                                purity_method = purity_method,
                                dataset=dataset,
                                datatype=datatype,
                                num_gene = num_gene,
                                immunesig = immunesig, 
                                model = model,
                                ACAT_Pval = 0.05,
                                work_path = work_path,
                                savepath = savepath)


savepath2 = "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/data/drugITS_meta_oe_original/"
dir.create(savepath2)
drugITGSscore_preparation(df_immunesig = df_immunesig,
                          cancer = cancer,
                          tumor_subtype = tumor_subtype,
                          purity_method = purity_method,
                          dataset=dataset,
                          datatype="allgenes",
                          num_gene = 200,
                          work_path = work_path,
                          save_path = savepath2)

# ITSselected <- read.csv("06_results_summaryandplotting/data/immune_sig_selected_new/selected_ITS_0.4.txt", sep = '\t')
meta_res <- read.csv("05_immuneSig_ICBresponse/results/meta_results/res_metaRes.csv")
meta_res_selected <- meta_res[meta_res$Meta_Pval < 0.05, ]
sig_s = meta_res_selected[meta_res_selected$OR > 1,]
sig_r = meta_res_selected[meta_res_selected$OR < 1,]
sig_s = data.frame( immune_sig = sig_s$X, flag = "sensitive", setindex = '01_sensitive')
sig_r = data.frame( immune_sig = sig_r$X, flag = "resistant", setindex = '02_resistant')

# mergedP_res <- read.csv("05_immuneSig_ICBresponse/results/meta_results/res_merged_p.csv")
# sig_s = mergedP_res[mergedP_res$merged_RgreaterNR_p_ACAT < 0.05,]
# sig_r = mergedP_res[mergedP_res$merged_RlessNR_p_ACAT < 0.05,]
# sig_s = data.frame( immune_sig = intersect(sig_s$X, meta_res_selected$X), flag = "sensitive", setindex = '01_sensitive')
# sig_r = data.frame( immune_sig = intersect(sig_r$X, meta_res_selected$X), flag = "resistant", setindex = '02_resistant')


load("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/05_immuneSig_ICBresponse/results_r0.4_ITSssgsea_CombatRef/finalmodel_17datasets_ssgseanorm/metaAnalysis_filter/en/elasticnet_model_pred.Rdata")

en_model_para = m_en$learner.model$glmnet.fit
m_en_imp <- data.frame(en_model_para$beta[,which(en_model_para$lambda == m_en$learner.model$lambda.min)])

names(m_en_imp) <- "weight"
m_en_imp$immune_sig = row.names(m_en_imp)


immunesig = rbind(sig_s, sig_r)
immunesig = merge(immunesig, m_en_imp, by = "immune_sig")
immunesig =immunesig[immunesig$weight != 0 ,]

# ------------------------------------------------------------------------------
# if drug effect on cell lines and mice models are consistant -------------------------
setwd("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/")
options(stringsAsFactors = F)
GEO_ID = "GSE160785"
drugindex = "drug989"


gsea_resultpath = paste0("result_", GEO_ID, "/drug_immunesig_", purity_method,"_GSEA_result_r",s1,"_meta_oe_original/", cancer, "_", tumor_subtype, "/")


# ITSproof=2

Pval = 1
dir.create(paste0("result_", GEO_ID, "/drug_immunesig_", purity_method,"_GSEA_result_r",s1,"_meta_oe_original_en_padj",Pval,"_final/"))
dir.create(paste0("result_", GEO_ID, "/drug_immunesig_", purity_method,"_GSEA_result_r",s1,"_meta_oe_original_en_padj",Pval,"_final/", cancer, "_", tumor_subtype, "/"))
savepath = paste0("result_", GEO_ID, "/drug_immunesig_", purity_method,"_GSEA_result_r",s1,"_meta_oe_original_en_padj",Pval,"_final/", cancer, "_", tumor_subtype, "/", drugindex)
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
                  compare_type ="rawITS",
                  Pval = Pval,
                  gsea_resultpath = gsea_resultpath,
                  savepath = savepath)


Pval = 0.1
dir.create(paste0("result_", GEO_ID, "/drug_immunesig_", purity_method,"_GSEA_result_r",s1,"_meta_oe_original_en_padj",Pval,"_final/"))
dir.create(paste0("result_", GEO_ID, "/drug_immunesig_", purity_method,"_GSEA_result_r",s1,"_meta_oe_original_en_padj",Pval,"_final/", cancer, "_", tumor_subtype, "/"))
savepath = paste0("result_", GEO_ID, "/drug_immunesig_", purity_method,"_GSEA_result_r",s1,"_meta_oe_original_en_padj",Pval,"_final/", cancer, "_", tumor_subtype, "/", drugindex)
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
                  compare_type ="rawITS",
                  Pval = Pval,
                  gsea_resultpath = gsea_resultpath,
                  savepath = savepath)

Pval = 0.05
dir.create(paste0("result_", GEO_ID, "/drug_immunesig_", purity_method,"_GSEA_result_r",s1,"_meta_oe_original_en_padj",Pval,"_final/"))
dir.create(paste0("result_", GEO_ID, "/drug_immunesig_", purity_method,"_GSEA_result_r",s1,"_meta_oe_original_en_padj",Pval,"_final/", cancer, "_", tumor_subtype, "/"))
savepath = paste0("result_", GEO_ID, "/drug_immunesig_", purity_method,"_GSEA_result_r",s1,"_meta_oe_original_en_padj",Pval,"_final/", cancer, "_", tumor_subtype, "/", drugindex)
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
                  compare_type ="rawITS",
                  Pval = Pval,
                  gsea_resultpath = gsea_resultpath,
                  savepath = savepath)

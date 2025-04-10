# rm(list=ls())
# ------------------------------------------------------------------------------
# source("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/06_results_summaryandplotting/scripts/immunesig_selection_function.R")
workpath <- "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/"
setwd(workpath)
options(stringsAsFactors = F)

# Checking and loading R packages -----------------------------------------
library(dplyr)
library(ggplot2)
library(pheatmap)
library(ROCR)
library(reshape2)

# Loading relied functions -----------------------------------------------------
# TIGS, AUC and Wilcox test
source("05_immuneSig_ICBresponse/scripts/functions/05_immunotherapy_immuneSig_analysis_function.R")
source("05_immuneSig_ICBresponse/scripts/functions/predictor_ITS.R")
source("06_results_summaryandplotting/scripts/immunesig_selection_function.R")

datasets <- list.files("05_immuneSig_ICBresponse/results_mergedITS_mclust_jaccard_r0.4_genevote_0.5/")
datasets <- datasets[-grep("ITS_wilcox_result", datasets)]
datasets <- datasets[-grep("res_r0.4_auc_2_mergedITS", datasets)]
datasets <- gsub(".Rdata", "", datasets)


# create savepath -----------------------------------------------------------------
dir.create("05_immuneSig_ICBresponse/results_r0.4_ITSssgsea_noCombat/")

# Hugo_dataset -----------------------------------------------------------------
cancer = "SKCM"
sample_type = "Metastatic"


# hugo --------------------------------------------------------------------------
# compute ITS and significant test between R and NR ----------------------------
# generating immune signatures -------------------------------------------------

processedDataPath = "05_immuneSig_ICBresponse/data/processedData/"
dir.create(processedDataPath)

# Hugo_dataset -----------------------------------------------------------------

dataset = "01_Hugo_dataset_outcome.Rdata"
load(paste0(processedDataPath, "/", dataset))
result_dir = "05_immuneSig_ICBresponse/results_r0.4_ITSssgsea_noCombat/01_Hugo_dataset_outcome/"
dir.create(result_dir)
dir.create("05_immuneSig_ICBresponse/results_r0.4_ITSssgsea_noCombat/ITS_wilcox_result_oe_original/")
dir.create("05_immuneSig_ICBresponse/results_r0.4_ITSssgsea_noCombat/immunesig_icbResponse_all_results_oe_original/")
cancer = "SKCM"
sample_type = "Metastatic"


immuneSig_res <- read.csv("05_immuneSig_ICBresponse/results/immunesig_icbResponse_test_results/01_Hugo_dataset_test_result.csv", row.names = 1)
names(immuneSig_res) <- c( "immuneSig_wilcox_P_value","immuneSig_R_greater_NR_P_value","immuneSig_R_less_NR_P_value", "immuneSig_auc","immune_sig"   )
immuneSig_res$immune_sig = tolower(immuneSig_res$immune_sig)
immuneSig_res = immunesig_name_unify(immuneSig_res)
savepath = paste0(result_dir, "/ITSasPredictor_oe_original/")
dir.create(savepath)

# savepath_allresult = "02_tumorGene_immuneSig_correlation/results_selected/immunecell_TCGAgenes_result_spearman_TUMERIC/"
savepath_allresult = "02_tumorGene_immuneSig_correlation/results_selected/ITS_oe_original_spearman_TUMERIC_r0.4/"
load(paste0(savepath_allresult, "for_GSVA/pcor/pcor_spearman_", cancer, "_",sample_type, "_positive_200.Rdata"))
load(paste0(savepath_allresult, "for_GSVA/pcor/pcor_spearman_", cancer, "_",sample_type, "_negative_200.Rdata"))
genelist_p_name = data.frame(immunesig = names(genelist_p), immune_sig = tolower(names(genelist_p)))
genelist_p_name = immunesig_name_unify(genelist_p_name)
names(genelist_p) = genelist_p_name$immune_sig

genelist_n_name = data.frame(immunesig = names(genelist_n), immune_sig = tolower(names(genelist_n)))
genelist_n_name = immunesig_name_unify(genelist_n_name)
names(genelist_n) = genelist_n_name$immune_sig


ori_TMEsig <- read.csv(paste0("05_immuneSig_ICBresponse/results/01_Hugo_dataset_outcome/SKCM_Metastatic_immuneSig_merged.csv"))
ori_TMEsig_name = data.frame(name = names(ori_TMEsig), immune_sig = tolower(names(ori_TMEsig)))
ori_TMEsig_name = immunesig_name_unify(ori_TMEsig_name)
names(ori_TMEsig) = ori_TMEsig_name$immune_sig

data =  as.matrix(log2(1+tpm_mat))
ITSgsva <- cor_immunesig_ITSgsva_ICBdata(data = data,
                              genelist_p = genelist_p,
                              genelist_n = genelist_n,
                              ori_TMEsig = ori_TMEsig,
                              type = "ITS", 
                              enrich_method = "ssgsea",
                              savepath = savepath,
                              savepath_dotplot = NULL)
                              

# riaz --------------------------------------------------------------------------

processedDataPath = "05_immuneSig_ICBresponse/data/processedData/"
dir.create(processedDataPath)

# Riaz_dataset -----------------------------------------------------------------

dataset = "02_Riaz_dataset_outcome.Rdata"
load(paste0(processedDataPath, "/", dataset))
result_dir = "05_immuneSig_ICBresponse/results_r0.4_ITSssgsea_noCombat/02_Riaz_dataset_outcome/"
dir.create(result_dir)
cancer = "SKCM"
sample_type = "Metastatic"


immuneSig_res <- read.csv("05_immuneSig_ICBresponse/results/immunesig_icbResponse_test_results/02_Riaz_dataset_test_result.csv", row.names = 1)
names(immuneSig_res) <- c( "immuneSig_wilcox_P_value","immuneSig_R_greater_NR_P_value","immuneSig_R_less_NR_P_value", "immuneSig_auc","immune_sig"   )
immuneSig_res$immune_sig = tolower(immuneSig_res$immune_sig)
immuneSig_res = immunesig_name_unify(immuneSig_res)
savepath = paste0(result_dir, "/ITSasPredictor_oe_original/")
dir.create(savepath)

savepath_allresult = "02_tumorGene_immuneSig_correlation/results_selected/ITS_oe_original_spearman_TUMERIC_r0.4/"
load(paste0(savepath_allresult, "for_GSVA/pcor/pcor_spearman_", cancer, "_",sample_type, "_positive_200.Rdata"))
load(paste0(savepath_allresult, "for_GSVA/pcor/pcor_spearman_", cancer, "_",sample_type, "_negative_200.Rdata"))
genelist_p_name = data.frame(immunesig = names(genelist_p), immune_sig = tolower(names(genelist_p)))
genelist_p_name = immunesig_name_unify(genelist_p_name)
names(genelist_p) = genelist_p_name$immune_sig

genelist_n_name = data.frame(immunesig = names(genelist_n), immune_sig = tolower(names(genelist_n)))
genelist_n_name = immunesig_name_unify(genelist_n_name)
names(genelist_n) = genelist_n_name$immune_sig



ori_TMEsig <- read.csv(paste0("05_immuneSig_ICBresponse/results/02_Riaz_dataset_outcome/SKCM_Metastatic_immuneSig_merged.csv"))
ori_TMEsig_name = data.frame(name = names(ori_TMEsig), immune_sig = tolower(names(ori_TMEsig)))
ori_TMEsig_name = immunesig_name_unify(ori_TMEsig_name)
names(ori_TMEsig) = ori_TMEsig_name$immune_sig

data =  as.matrix(log2(1+tpm_mat))

# gsva -----------------------------------------------------------------------
enrich_method = "ssgsea"
ITSgsva <- cor_immunesig_ITSgsva_ICBdata(data = data,
                              genelist_p = genelist_p,
                              genelist_n = genelist_n,
                              ori_TMEsig = ori_TMEsig,
                              type = "ITS", 
                              enrich_method = "ssgsea",
                              savepath = savepath,
                              savepath_dotplot = NULL)
                              

# Gide --------------------------------------------------------------------------
datasets[3]

processedDataPath = "05_immuneSig_ICBresponse/data/processedData/"
dir.create(processedDataPath)
# Gide_dataset -----------------------------------------------------------------

dataset = "05_Gide_dataset_outcome.Rdata"
load(paste0(processedDataPath, "/", dataset))
result_dir = "05_immuneSig_ICBresponse/results_r0.4_ITSssgsea_noCombat/05_Gide_dataset_outcome/"

dir.create(result_dir)
cancer = "SKCM"
sample_type = "Metastatic"

immuneSig_res <- read.csv("05_immuneSig_ICBresponse/results/immunesig_icbResponse_test_results/05_Gide_dataset_test_result.csv", row.names = 1)
names(immuneSig_res) <- c( "immuneSig_wilcox_P_value","immuneSig_R_greater_NR_P_value","immuneSig_R_less_NR_P_value", "immuneSig_auc","immune_sig"   )
immuneSig_res$immune_sig = tolower(immuneSig_res$immune_sig)
immuneSig_res = immunesig_name_unify(immuneSig_res)
savepath = paste0(result_dir, "/ITSasPredictor_oe_original/")
dir.create(savepath)

savepath_allresult = "02_tumorGene_immuneSig_correlation/results_selected/ITS_oe_original_spearman_TUMERIC_r0.4/"
load(paste0(savepath_allresult, "for_GSVA/pcor/pcor_spearman_", cancer, "_",sample_type, "_positive_200.Rdata"))
load(paste0(savepath_allresult, "for_GSVA/pcor/pcor_spearman_", cancer, "_",sample_type, "_negative_200.Rdata"))
genelist_p_name = data.frame(immunesig = names(genelist_p), immune_sig = tolower(names(genelist_p)))
genelist_p_name = immunesig_name_unify(genelist_p_name)
names(genelist_p) = genelist_p_name$immune_sig

genelist_n_name = data.frame(immunesig = names(genelist_n), immune_sig = tolower(names(genelist_n)))
genelist_n_name = immunesig_name_unify(genelist_n_name)
names(genelist_n) = genelist_n_name$immune_sig



ori_TMEsig <- read.csv(paste0("05_immuneSig_ICBresponse/results/05_Gide_dataset_outcome/SKCM_Metastatic_immuneSig_merged.csv"))
ori_TMEsig_name = data.frame(name = names(ori_TMEsig), immune_sig = tolower(names(ori_TMEsig)))
ori_TMEsig_name = immunesig_name_unify(ori_TMEsig_name)
names(ori_TMEsig) = ori_TMEsig_name$immune_sig


data =  as.matrix(log2(1+tpm_mat))

# gsva -----------------------------------------------------------------------

enrich_method = "ssgsea"
ITSgsva <- cor_immunesig_ITSgsva_ICBdata(data = data,
                              genelist_p = genelist_p,
                              genelist_n = genelist_n,
                              ori_TMEsig = ori_TMEsig,
                              type = "ITS", 
                              enrich_method = "ssgsea",
                              savepath = savepath,
                              savepath_dotplot = NULL)
                              


# Kim --------------------------------------------------------------------------
datasets[4]
processedDataPath = "05_immuneSig_ICBresponse/data/processedData/"
dir.create(processedDataPath)


# Kim_dataset ------------------------------------------------------------------
dataset = "06_Kim_dataset_outcome.Rdata"
load(paste0(processedDataPath, "/", dataset))
result_dir = "05_immuneSig_ICBresponse/results_r0.4_ITSssgsea_noCombat/06_Kim_dataset_outcome/"
dir.create(result_dir)
cancer = "STAD"
sample_type = "Primary"


immuneSig_res <- read.csv("05_immuneSig_ICBresponse/results/immunesig_icbResponse_test_results/06_Kim_dataset_test_result.csv", row.names = 1)
names(immuneSig_res) <- c( "immuneSig_wilcox_P_value","immuneSig_R_greater_NR_P_value","immuneSig_R_less_NR_P_value", "immuneSig_auc","immune_sig"   )
immuneSig_res$immune_sig = tolower(immuneSig_res$immune_sig)
immuneSig_res = immunesig_name_unify(immuneSig_res)
savepath = paste0(result_dir, "/ITSasPredictor_oe_original/")
dir.create(savepath)

savepath_allresult = "02_tumorGene_immuneSig_correlation/results_selected/ITS_oe_original_spearman_TUMERIC_r0.4/"
load(paste0(savepath_allresult, "for_GSVA/pcor/pcor_spearman_", cancer, "_",sample_type, "_positive_200.Rdata"))
load(paste0(savepath_allresult, "for_GSVA/pcor/pcor_spearman_", cancer, "_",sample_type, "_negative_200.Rdata"))
genelist_p_name = data.frame(immunesig = names(genelist_p), immune_sig = tolower(names(genelist_p)))
genelist_p_name = immunesig_name_unify(genelist_p_name)
names(genelist_p) = genelist_p_name$immune_sig

genelist_n_name = data.frame(immunesig = names(genelist_n), immune_sig = tolower(names(genelist_n)))
genelist_n_name = immunesig_name_unify(genelist_n_name)
names(genelist_n) = genelist_n_name$immune_sig


ori_TMEsig <- read.csv(paste0("05_immuneSig_ICBresponse/results/06_Kim_dataset_outcome/STAD_Metastatic_immuneSig_merged.csv"))
ori_TMEsig_name = data.frame(name = names(ori_TMEsig), immune_sig = tolower(names(ori_TMEsig)))
ori_TMEsig_name = immunesig_name_unify(ori_TMEsig_name)
names(ori_TMEsig) = ori_TMEsig_name$immune_sig


data =  as.matrix(log2(1+tpm_mat))

# gsva -----------------------------------------------------------------------
enrich_method = "ssgsea"
ITSgsva <- cor_immunesig_ITSgsva_ICBdata(data = data,
                              genelist_p = genelist_p,
                              genelist_n = genelist_n,
                              ori_TMEsig = ori_TMEsig,
                              type = "ITS", 
                              enrich_method = "ssgsea",
                              savepath = savepath,
                              savepath_dotplot = NULL)
                              

# IMvigor210 --------------------------------------------------------------------------
datasets[5]
processedDataPath = "/picb/bigdata/project/CancerSysBio/Immunotherapy_prediction/data/RNAseq/"
dir.create(processedDataPath)


# IMvigor210_dataset ------------------------------------------------------------------
dataset = "IMvigor210_bladder.Rdata"
load(paste0(processedDataPath, "/IMvigor210/processedData/", dataset))
result_dir = "05_immuneSig_ICBresponse/results_r0.4_ITSssgsea_noCombat/07_IMvigor210_bladder_immunotherapy/"
dir.create(result_dir)


tpm_mat <- tpm_mat[rowSums(tpm_mat) != 0, ]
tpm_mat <- tpm_mat[-which(rownames(tpm_mat)==""), ]


immuneSig_res <- read.csv("05_immuneSig_ICBresponse/results/immunesig_icbResponse_test_results/07_IMvigor210_dataset_test_result.csv", row.names = 1)
names(immuneSig_res) <- c( "immuneSig_wilcox_P_value","immuneSig_R_greater_NR_P_value","immuneSig_R_less_NR_P_value", "immuneSig_auc","immune_sig"   )
immuneSig_res$immune_sig = tolower(immuneSig_res$immune_sig)
immuneSig_res = immunesig_name_unify(immuneSig_res)
savepath = paste0(result_dir, "/ITSasPredictor_oe_original/")
dir.create(savepath)

# savepath_allresult = "02_tumorGene_immuneSig_correlation/results_selected/immunecell_TCGAgenes_result_spearman_TUMERIC/"
cancer = "BLCA"
sample_type = "Primary"
savepath_allresult = "02_tumorGene_immuneSig_correlation/results_selected/ITS_oe_original_spearman_TUMERIC_r0.4/"
load(paste0(savepath_allresult, "for_GSVA/pcor/pcor_spearman_", cancer, "_",sample_type, "_positive_200.Rdata"))
load(paste0(savepath_allresult, "for_GSVA/pcor/pcor_spearman_", cancer, "_",sample_type, "_negative_200.Rdata"))
genelist_p_name = data.frame(immunesig = names(genelist_p), immune_sig = tolower(names(genelist_p)))
genelist_p_name = immunesig_name_unify(genelist_p_name)
names(genelist_p) = genelist_p_name$immune_sig

genelist_n_name = data.frame(immunesig = names(genelist_n), immune_sig = tolower(names(genelist_n)))
genelist_n_name = immunesig_name_unify(genelist_n_name)
names(genelist_n) = genelist_n_name$immune_sig


ori_TMEsig <- read.csv(paste0("05_immuneSig_ICBresponse/results/07_IMvigor210_bladder_immunotherapy/BLCA_Metastatic_immuneSig_merged.csv"))
ori_TMEsig_name = data.frame(name = names(ori_TMEsig), immune_sig = tolower(names(ori_TMEsig)))
ori_TMEsig_name = immunesig_name_unify(ori_TMEsig_name)
names(ori_TMEsig) = ori_TMEsig_name$immune_sig


data =  as.matrix(log2(1+tpm_mat))

# gsva -----------------------------------------------------------------------
enrich_method = "ssgsea"
ITSgsva <- cor_immunesig_ITSgsva_ICBdata(data = data,
                              genelist_p = genelist_p,
                              genelist_n = genelist_n,
                              ori_TMEsig = ori_TMEsig,
                              type = "ITS", 
                              enrich_method = "ssgsea",
                              savepath = savepath,
                              savepath_dotplot = NULL)
                              

                      

# Braun  ----------------------------------------------------------
datasets[7]

processedDataPath = "/picb/bigdata/project/CancerSysBio/Immunotherapy_prediction/data/RNAseq/"
dir.create(processedDataPath)


# Braun_dataset_PD1 ------------------------------------------------------------------
dataset = "Braun_dataset_PD1.Rdata"
load(paste0(processedDataPath, "/Braun_dataset/processedData/", dataset))
result_dir = "05_immuneSig_ICBresponse/results_r0.4_ITSssgsea_noCombat/Braun_dataset/"
dir.create(result_dir)
dir.create("05_immuneSig_ICBresponse/results_r0.4_ITSssgsea_noCombat/ITS_wilcox_result_oe_original/")

dir.create(result_dir)
cancer = "KIRC"
sample_type = "Primary"

immuneSig_res <- read.csv("05_immuneSig_ICBresponse/results/immunesig_icbResponse_test_results/Braun_dataset_PD1_test_result.csv", row.names = 1)
names(immuneSig_res) <- c( "immuneSig_wilcox_P_value","immuneSig_R_greater_NR_P_value","immuneSig_R_less_NR_P_value", "immuneSig_auc","immune_sig"   )
immuneSig_res$immune_sig = tolower(immuneSig_res$immune_sig)
immuneSig_res = immunesig_name_unify(immuneSig_res)
savepath = paste0(result_dir, "/ITSasPredictor_oe_original/")
dir.create(savepath)

savepath_allresult = "02_tumorGene_immuneSig_correlation/results_selected/ITS_oe_original_spearman_TUMERIC_r0.4/"
load(paste0(savepath_allresult, "for_GSVA/pcor/pcor_spearman_", cancer, "_",sample_type, "_positive_200.Rdata"))
load(paste0(savepath_allresult, "for_GSVA/pcor/pcor_spearman_", cancer, "_",sample_type, "_negative_200.Rdata"))
genelist_p_name = data.frame(immunesig = names(genelist_p), immune_sig = tolower(names(genelist_p)))
genelist_p_name = immunesig_name_unify(genelist_p_name)
names(genelist_p) = genelist_p_name$immune_sig

genelist_n_name = data.frame(immunesig = names(genelist_n), immune_sig = tolower(names(genelist_n)))
genelist_n_name = immunesig_name_unify(genelist_n_name)
names(genelist_n) = genelist_n_name$immune_sig


ori_TMEsig <- read.csv(paste0("05_immuneSig_ICBresponse/results/Braun_dataset/KIRC_Primary_immuneSig_merged.csv"))
ori_TMEsig_name = data.frame(name = names(ori_TMEsig), immune_sig = tolower(names(ori_TMEsig)))
ori_TMEsig_name = immunesig_name_unify(ori_TMEsig_name)
names(ori_TMEsig) = ori_TMEsig_name$immune_sig


data =  as.matrix(log2(1+tpm_mat))

# gsva -----------------------------------------------------------------------

enrich_method = "ssgsea"
ITSgsva <- cor_immunesig_ITSgsva_ICBdata(data = data,
                              genelist_p = genelist_p,
                              genelist_n = genelist_n,
                              ori_TMEsig = ori_TMEsig,
                              type = "ITS", 
                              enrich_method = "ssgsea",
                              savepath = savepath,
                              savepath_dotplot = NULL)
                              



# GSE126044 ----------------------------------------------------------
datasets[8]
processedDataPath = "/picb/bigdata/project/CancerSysBio/Immunotherapy_prediction/data/RNAseq/"
dir.create(processedDataPath)


# GSE126044 ------------------------------------------------------------------
dataset = "GSE126044.Rdata"
load(paste0(processedDataPath, "/GSE126044/processedData/", dataset))
result_dir = "05_immuneSig_ICBresponse/results_r0.4_ITSssgsea_noCombat/GSE126044/"
dir.create(result_dir)
dir.create("05_immuneSig_ICBresponse/results_r0.4_ITSssgsea_noCombat/ITS_wilcox_result_oe_original/")
cancer = "LUAD"
sample_type = "Primary"

if(length(which(rowSums(tpm_mat) == 0)) > 0)tpm_mat <- tpm_mat[rowSums(tpm_mat) != 0, ]
if(length(which(rownames(tpm_mat)=="")) > 0)tpm_mat <- tpm_mat[-which(rownames(tpm_mat)==""), ]

# -----------------------------------------------------------------------

immuneSig_res <- read.csv("05_immuneSig_ICBresponse/results/immunesig_icbResponse_test_results/GSE126044_test_result.csv", row.names = 1)
names(immuneSig_res) <- c( "immuneSig_wilcox_P_value","immuneSig_R_greater_NR_P_value","immuneSig_R_less_NR_P_value", "immuneSig_auc","immune_sig"   )
immuneSig_res$immune_sig = tolower(immuneSig_res$immune_sig)
immuneSig_res = immunesig_name_unify(immuneSig_res)
savepath = paste0(result_dir, "/ITSasPredictor_oe_original/")
dir.create(savepath)

savepath_allresult = "02_tumorGene_immuneSig_correlation/results_selected/ITS_oe_original_spearman_TUMERIC_r0.4/"
load(paste0(savepath_allresult, "for_GSVA/pcor/pcor_spearman_", cancer, "_",sample_type, "_positive_200.Rdata"))
load(paste0(savepath_allresult, "for_GSVA/pcor/pcor_spearman_", cancer, "_",sample_type, "_negative_200.Rdata"))
genelist_p_name = data.frame(immunesig = names(genelist_p), immune_sig = tolower(names(genelist_p)))
genelist_p_name = immunesig_name_unify(genelist_p_name)
names(genelist_p) = genelist_p_name$immune_sig

genelist_n_name = data.frame(immunesig = names(genelist_n), immune_sig = tolower(names(genelist_n)))
genelist_n_name = immunesig_name_unify(genelist_n_name)
names(genelist_n) = genelist_n_name$immune_sig


ori_TMEsig <- read.csv(paste0("05_immuneSig_ICBresponse/results/GSE126044/LUAD_Primary_immuneSig_merged.csv"))
ori_TMEsig_name = data.frame(name = names(ori_TMEsig), immune_sig = tolower(names(ori_TMEsig)))
ori_TMEsig_name = immunesig_name_unify(ori_TMEsig_name)
names(ori_TMEsig) = ori_TMEsig_name$immune_sig



data =  as.matrix(log2(1+tpm_mat))

# gsva -----------------------------------------------------------------------
enrich_method = "ssgsea"
ITSgsva <- cor_immunesig_ITSgsva_ICBdata(data = data,
                              genelist_p = genelist_p,
                              genelist_n = genelist_n,
                              ori_TMEsig = ori_TMEsig,
                              type = "ITS", 
                              enrich_method = "ssgsea",
                              savepath = savepath,
                              savepath_dotplot = NULL)
                              

# GSE135222 ----------------------------------------------------------
datasets[9]


processedDataPath = "/picb/bigdata/project/CancerSysBio/Immunotherapy_prediction/data/RNAseq/"
dir.create(processedDataPath)


# GSE135222 ------------------------------------------------------------------
dataset = "GSE135222.Rdata"
load(paste0(processedDataPath, "/GSE135222/processedData/", dataset))
result_dir = "05_immuneSig_ICBresponse/results_r0.4_ITSssgsea_noCombat/GSE135222/"
dir.create(result_dir)
cancer = "LUAD"
sample_type = "Primary"

if(length(which(rowSums(tpm_mat) == 0)) > 0)tpm_mat <- tpm_mat[rowSums(tpm_mat) != 0, ]
if(length(which(rownames(tpm_mat)=="")) > 0)tpm_mat <- tpm_mat[-which(rownames(tpm_mat)==""), ]
#  --------------------------------------------------------------------------
immuneSig_res <- read.csv("05_immuneSig_ICBresponse/results/immunesig_icbResponse_test_results/GSE135222_test_result.csv", row.names = 1)
names(immuneSig_res) <- c( "immuneSig_wilcox_P_value","immuneSig_R_greater_NR_P_value","immuneSig_R_less_NR_P_value", "immuneSig_auc","immune_sig"   )
immuneSig_res$immune_sig = tolower(immuneSig_res$immune_sig)
immuneSig_res = immunesig_name_unify(immuneSig_res)
savepath = paste0(result_dir, "/ITSasPredictor_oe_original/")
dir.create(savepath)

savepath_allresult = "02_tumorGene_immuneSig_correlation/results_selected/ITS_oe_original_spearman_TUMERIC_r0.4/"
load(paste0(savepath_allresult, "for_GSVA/pcor/pcor_spearman_", cancer, "_",sample_type, "_positive_200.Rdata"))
load(paste0(savepath_allresult, "for_GSVA/pcor/pcor_spearman_", cancer, "_",sample_type, "_negative_200.Rdata"))

genelist_p_name = data.frame(immunesig = names(genelist_p), immune_sig = tolower(names(genelist_p)))
genelist_p_name = immunesig_name_unify(genelist_p_name)
names(genelist_p) = genelist_p_name$immune_sig

genelist_n_name = data.frame(immunesig = names(genelist_n), immune_sig = tolower(names(genelist_n)))
genelist_n_name = immunesig_name_unify(genelist_n_name)
names(genelist_n) = genelist_n_name$immune_sig


ori_TMEsig <- read.csv(paste0("05_immuneSig_ICBresponse/results/GSE135222/LUAD_Primary_immuneSig_merged.csv"))
ori_TMEsig_name = data.frame(name = names(ori_TMEsig), immune_sig = tolower(names(ori_TMEsig)))
ori_TMEsig_name = immunesig_name_unify(ori_TMEsig_name)
names(ori_TMEsig) = ori_TMEsig_name$immune_sig


data =  as.matrix(log2(1+tpm_mat))

# gsva -----------------------------------------------------------------------
enrich_method = "ssgsea"

ITSgsva <- cor_immunesig_ITSgsva_ICBdata(data = data,
                              genelist_p = genelist_p,
                              genelist_n = genelist_n,
                              ori_TMEsig = ori_TMEsig,
                              type = "ITS", 
                              enrich_method = "ssgsea",
                              savepath = savepath,
                              savepath_dotplot = NULL)
                              
# GSE145996 ----------------------------------------------------------
datasets[10]

processedDataPath = "/picb/bigdata/project/CancerSysBio/Immunotherapy_prediction/data/RNAseq/"
dir.create(processedDataPath)


# GSE145996 ------------------------------------------------------------------
dataset = "GSE145996.Rdata"
load(paste0(processedDataPath, "/GSE145996/processedData/", dataset))
result_dir = "05_immuneSig_ICBresponse/results_r0.4_ITSssgsea_noCombat/GSE145996/"
dir.create(result_dir)
cancer = "SKCM"
sample_type = "Metastatic"

if(length(which(rowSums(tpm_mat) == 0)) > 0)tpm_mat <- tpm_mat[rowSums(tpm_mat) != 0, ]
if(length(which(rownames(tpm_mat)=="")) > 0)tpm_mat <- tpm_mat[-which(rownames(tpm_mat)==""), ]
#  --------------------------------------------------------------------------
immuneSig_res <- read.csv("05_immuneSig_ICBresponse/results/immunesig_icbResponse_test_results/GSE145996_test_result.csv", row.names = 1)
names(immuneSig_res) <- c( "immuneSig_wilcox_P_value","immuneSig_R_greater_NR_P_value","immuneSig_R_less_NR_P_value", "immuneSig_auc","immune_sig"   )
immuneSig_res$immune_sig = tolower(immuneSig_res$immune_sig)
immuneSig_res = immunesig_name_unify(immuneSig_res)
savepath = paste0(result_dir, "/ITSasPredictor_oe_original/")
dir.create(savepath)

savepath_allresult = "02_tumorGene_immuneSig_correlation/results_selected/ITS_oe_original_spearman_TUMERIC_r0.4/"
load(paste0(savepath_allresult, "for_GSVA/pcor/pcor_spearman_", cancer, "_",sample_type, "_positive_200.Rdata"))
load(paste0(savepath_allresult, "for_GSVA/pcor/pcor_spearman_", cancer, "_",sample_type, "_negative_200.Rdata"))
genelist_p_name = data.frame(immunesig = names(genelist_p), immune_sig = tolower(names(genelist_p)))
genelist_p_name = immunesig_name_unify(genelist_p_name)
names(genelist_p) = genelist_p_name$immune_sig

genelist_n_name = data.frame(immunesig = names(genelist_n), immune_sig = tolower(names(genelist_n)))
genelist_n_name = immunesig_name_unify(genelist_n_name)
names(genelist_n) = genelist_n_name$immune_sig


ori_TMEsig <- read.csv(paste0("05_immuneSig_ICBresponse/results/GSE145996/SKCM_Metastatic_immuneSig_merged.csv"))
ori_TMEsig_name = data.frame(name = names(ori_TMEsig), immune_sig = tolower(names(ori_TMEsig)))
ori_TMEsig_name = immunesig_name_unify(ori_TMEsig_name)
names(ori_TMEsig) = ori_TMEsig_name$immune_sig



data =  as.matrix(log2(1+tpm_mat))

# gsva -----------------------------------------------------------------------
enrich_method = "ssgsea"

ITSgsva <- cor_immunesig_ITSgsva_ICBdata(data = data,
                              genelist_p = genelist_p,
                              genelist_n = genelist_n,
                              ori_TMEsig = ori_TMEsig,
                              type = "ITS", 
                              enrich_method = "ssgsea",
                              savepath = savepath,
                              savepath_dotplot = NULL)
                              
# ------------------------------------------------------------------------------


# Nathanson ----------------------------------------------------------
datasets[11]

processedDataPath = "/picb/bigdata/project/CancerSysBio/Immunotherapy_prediction/data/"

dataset = "Nathanson_dataset.Rdata"
load(paste0(processedDataPath, "RNAseq/Nathanson_dataset/processedData/", dataset))
result_dir = "05_immuneSig_ICBresponse/results_r0.4_ITSssgsea_noCombat/Nathanson_dataset_outcome/"
dir.create(result_dir)
cancer = "SKCM"
sample_type = "Metastatic"
if(length(which(rowSums(tpm_mat) == 0)) > 0)tpm_mat <- tpm_mat[rowSums(tpm_mat) != 0, ]
if(length(which(rownames(tpm_mat)=="")) > 0)tpm_mat <- tpm_mat[-which(rownames(tpm_mat)==""), ]

#  --------------------------------------------------------------------------
immuneSig_res <- read.csv("05_immuneSig_ICBresponse/results/immunesig_icbResponse_test_results/Nathanson_dataset_test_result.csv", row.names = 1)
names(immuneSig_res) <- c( "immuneSig_wilcox_P_value","immuneSig_R_greater_NR_P_value","immuneSig_R_less_NR_P_value", "immuneSig_auc","immune_sig"   )
immuneSig_res$immune_sig = tolower(immuneSig_res$immune_sig)
immuneSig_res = immunesig_name_unify(immuneSig_res)
savepath = paste0(result_dir, "/ITSasPredictor_oe_original/")
dir.create(savepath)

savepath_allresult = "02_tumorGene_immuneSig_correlation/results_selected/ITS_oe_original_spearman_TUMERIC_r0.4/"
load(paste0(savepath_allresult, "for_GSVA/pcor/pcor_spearman_", cancer, "_",sample_type, "_positive_200.Rdata"))
load(paste0(savepath_allresult, "for_GSVA/pcor/pcor_spearman_", cancer, "_",sample_type, "_negative_200.Rdata"))
genelist_p_name = data.frame(immunesig = names(genelist_p), immune_sig = tolower(names(genelist_p)))
genelist_p_name = immunesig_name_unify(genelist_p_name)
names(genelist_p) = genelist_p_name$immune_sig

genelist_n_name = data.frame(immunesig = names(genelist_n), immune_sig = tolower(names(genelist_n)))
genelist_n_name = immunesig_name_unify(genelist_n_name)
names(genelist_n) = genelist_n_name$immune_sig

ori_TMEsig <- read.csv(paste0("05_immuneSig_ICBresponse/results/Nathanson_dataset_outcome/SKCM_Metastatic_immuneSig_merged.csv"))
ori_TMEsig_name = data.frame(name = names(ori_TMEsig), immune_sig = tolower(names(ori_TMEsig)))
ori_TMEsig_name = immunesig_name_unify(ori_TMEsig_name)
names(ori_TMEsig) = ori_TMEsig_name$immune_sig



data =  as.matrix(log2(1+tpm_mat))

# gsva -----------------------------------------------------------------------
enrich_method = "ssgsea"
ITSgsva <- cor_immunesig_ITSgsva_ICBdata(data = data,
                              genelist_p = genelist_p,
                              genelist_n = genelist_n,
                              ori_TMEsig = ori_TMEsig,
                              type = "ITS", 
                              enrich_method = "ssgsea",
                              savepath = savepath,
                              savepath_dotplot = NULL)
                              

# ------------------------------------------------------------------------------

# PMID_29301960 ----------------------------------------------------------
datasets[12]

#  --------------------------------------------------------------------------
processedDataPath = "/picb/bigdata/project/CancerSysBio/Immunotherapy_prediction/data/"
# PMID_29301960 ------------------------------------------------------------------
dataset = "PMID_29301960.Rdata"
load(paste0(processedDataPath, "RNAseq/PMID_29301960/processedData/", dataset))
result_dir = "05_immuneSig_ICBresponse/results_r0.4_ITSssgsea_noCombat/PMID_29301960/"
dir.create(result_dir)
cancer = "KIRC"
sample_type = "Primary"

if(length(which(rowSums(tpm_mat) == 0)) > 0)tpm_mat <- tpm_mat[rowSums(tpm_mat) != 0, ]
if(length(which(rownames(tpm_mat)=="")) > 0)tpm_mat <- tpm_mat[-which(rownames(tpm_mat)==""), ]

#  --------------------------------------------------------------------------
immuneSig_res <- read.csv("05_immuneSig_ICBresponse/results/immunesig_icbResponse_test_results/PMID_29301960_test_result.csv", row.names = 1)
names(immuneSig_res) <- c( "immuneSig_wilcox_P_value","immuneSig_R_greater_NR_P_value","immuneSig_R_less_NR_P_value", "immuneSig_auc","immune_sig"   )
immuneSig_res$immune_sig = tolower(immuneSig_res$immune_sig)
immuneSig_res = immunesig_name_unify(immuneSig_res)
savepath = paste0(result_dir, "/ITSasPredictor_oe_original/")
dir.create(savepath)

savepath_allresult = "02_tumorGene_immuneSig_correlation/results_selected/ITS_oe_original_spearman_TUMERIC_r0.4/"
load(paste0(savepath_allresult, "for_GSVA/pcor/pcor_spearman_", cancer, "_",sample_type, "_positive_200.Rdata"))
load(paste0(savepath_allresult, "for_GSVA/pcor/pcor_spearman_", cancer, "_",sample_type, "_negative_200.Rdata"))
genelist_p_name = data.frame(immunesig = names(genelist_p), immune_sig = tolower(names(genelist_p)))
genelist_p_name = immunesig_name_unify(genelist_p_name)
names(genelist_p) = genelist_p_name$immune_sig

genelist_n_name = data.frame(immunesig = names(genelist_n), immune_sig = tolower(names(genelist_n)))
genelist_n_name = immunesig_name_unify(genelist_n_name)
names(genelist_n) = genelist_n_name$immune_sig


ori_TMEsig <- read.csv(paste0("05_immuneSig_ICBresponse/results/PMID_29301960/KIRC_Primary_immuneSig_merged.csv"))
ori_TMEsig_name = data.frame(name = names(ori_TMEsig), immune_sig = tolower(names(ori_TMEsig)))
ori_TMEsig_name = immunesig_name_unify(ori_TMEsig_name)
names(ori_TMEsig) = ori_TMEsig_name$immune_sig

#-------------------------------------------------------------------------------

data =  as.matrix(log2(1+tpm_mat))

# gsva -----------------------------------------------------------------------
enrich_method = "ssgsea"
ITSgsva <- cor_immunesig_ITSgsva_ICBdata(data = data,
                              genelist_p = genelist_p,
                              genelist_n = genelist_n,
                              ori_TMEsig = ori_TMEsig,
                              type = "ITS", 
                              enrich_method = "ssgsea",
                              savepath = savepath,
                              savepath_dotplot = NULL)
ITSp_res = ITSgsva[[1]]
ITSn_res = ITSgsva[[2]]


# zhao ----------------------------------------------------------
datasets[13]
# generating immune signatures -------------------------------------------------

processedDataPath = "/picb/bigdata/project/CancerSysBio/Immunotherapy_prediction/data/RNAseq/"
dir.create(processedDataPath)

# Zhao_dataset ------------------------------------------------------------------
dataset = "Zhao_dataset.Rdata"
load(paste0(processedDataPath, "/Zhao_dataset/processedData/", dataset))
result_dir = "05_immuneSig_ICBresponse/results_r0.4_ITSssgsea_noCombat/Zhao_dataset/"
dir.create(result_dir)
cancer = "GBM"
sample_type = "Primary"

if(length(which(rowSums(tpm_mat) == 0)) > 0)tpm_mat <- tpm_mat[rowSums(tpm_mat) != 0, ]
if(length(which(rownames(tpm_mat)=="")) > 0)tpm_mat <- tpm_mat[-which(rownames(tpm_mat)==""), ]
response <- response[!is.na(response$label), ]

#  --------------------------------------------------------------------------
#  --------------------------------------------------------------------------
immuneSig_res <- read.csv("05_immuneSig_ICBresponse/results/immunesig_icbResponse_test_results/Zhao_dataset_test_result.csv", row.names = 1)
names(immuneSig_res) <- c( "immuneSig_wilcox_P_value","immuneSig_R_greater_NR_P_value","immuneSig_R_less_NR_P_value", "immuneSig_auc","immune_sig"   )
immuneSig_res$immune_sig = tolower(immuneSig_res$immune_sig)
immuneSig_res = immunesig_name_unify(immuneSig_res)
savepath = paste0(result_dir, "/ITSasPredictor_oe_original/")
dir.create(savepath)

savepath_allresult = "02_tumorGene_immuneSig_correlation/results_selected/ITS_oe_original_spearman_TUMERIC_r0.4/"
load(paste0(savepath_allresult, "for_GSVA/pcor/pcor_spearman_", cancer, "_",sample_type, "_positive_200.Rdata"))
load(paste0(savepath_allresult, "for_GSVA/pcor/pcor_spearman_", cancer, "_",sample_type, "_negative_200.Rdata"))
genelist_p_name = data.frame(immunesig = names(genelist_p), immune_sig = tolower(names(genelist_p)))
genelist_p_name = immunesig_name_unify(genelist_p_name)
names(genelist_p) = genelist_p_name$immune_sig

genelist_n_name = data.frame(immunesig = names(genelist_n), immune_sig = tolower(names(genelist_n)))
genelist_n_name = immunesig_name_unify(genelist_n_name)
names(genelist_n) = genelist_n_name$immune_sig


ori_TMEsig <- read.csv(paste0("05_immuneSig_ICBresponse/results/Zhao_dataset/GBM_Primary_immuneSig_merged.csv"))
ori_TMEsig_name = data.frame(name = names(ori_TMEsig), immune_sig = tolower(names(ori_TMEsig)))
ori_TMEsig_name = immunesig_name_unify(ori_TMEsig_name)
names(ori_TMEsig) = ori_TMEsig_name$immune_sig


data =  as.matrix(log2(1+tpm_mat))

# gsva -----------------------------------------------------------------------
enrich_method = "ssgsea"
ITSgsva <- cor_immunesig_ITSgsva_ICBdata(data = data,
                              genelist_p = genelist_p,
                              genelist_n = genelist_n,
                              ori_TMEsig = ori_TMEsig,
                              type = "ITS", 
                              enrich_method = "ssgsea",
                              savepath = savepath,
                              savepath_dotplot = NULL)







# ----------------------------------------------------------
# ----------------------------------------------------------
# ----------------------------------------------------------
# ----------------------------------------------------------

# validation set--------------------------------------------------

# ----------------------------------------------------------
# ----------------------------------------------------------
# ----------------------------------------------------------
# ----------------------------------------------------------

# Lauss-ACT therapy  ----------------------------------------------------------
datasets[6]

# # GSE100797 ------------------------------------------------------------------
# dataset = "GSE100797.Rdata"
# load(paste0(processedDataPath, "/Lauss_dataset/processedData/", dataset))

processedDataPath = "05_immuneSig_ICBresponse/data/processedData/"
dir.create(processedDataPath)


dataset = "08_Lauss_melanoma.Rdata"
load(paste0(processedDataPath, "/", dataset))

result_dir = "05_immuneSig_ICBresponse/results_r0.4_ITSssgsea_noCombat/08_Lauss_dataset_outcome/"
dir.create(result_dir)
cancer = "SKCM"
sample_type = "Metastatic"


immuneSig_res <- read.csv("05_immuneSig_ICBresponse/results/immunesig_icbResponse_test_results/08_Lauss_dataset_test_result.csv", row.names = 1)
names(immuneSig_res) <- c( "immuneSig_wilcox_P_value","immuneSig_R_greater_NR_P_value","immuneSig_R_less_NR_P_value", "immuneSig_auc","immune_sig"   )
immuneSig_res$immune_sig = tolower(immuneSig_res$immune_sig)
immuneSig_res = immunesig_name_unify(immuneSig_res)
savepath = paste0(result_dir, "/ITSasPredictor_oe_original/")
dir.create(savepath)

savepath_allresult = "02_tumorGene_immuneSig_correlation/results_selected/ITS_oe_original_spearman_TUMERIC_r0.4/"
load(paste0(savepath_allresult, "for_GSVA/pcor/pcor_spearman_", cancer, "_",sample_type, "_positive_200.Rdata"))
load(paste0(savepath_allresult, "for_GSVA/pcor/pcor_spearman_", cancer, "_",sample_type, "_negative_200.Rdata"))
genelist_p_name = data.frame(immunesig = names(genelist_p), immune_sig = tolower(names(genelist_p)))
genelist_p_name = immunesig_name_unify(genelist_p_name)
names(genelist_p) = genelist_p_name$immune_sig

genelist_n_name = data.frame(immunesig = names(genelist_n), immune_sig = tolower(names(genelist_n)))
genelist_n_name = immunesig_name_unify(genelist_n_name)
names(genelist_n) = genelist_n_name$immune_sig


ori_TMEsig <- read.csv(paste0("05_immuneSig_ICBresponse/results/08_Lauss_dataset_outcome/SKCM_Metastatic_immuneSig_merged.csv"))
ori_TMEsig_name = data.frame(name = names(ori_TMEsig), immune_sig = tolower(names(ori_TMEsig)))
ori_TMEsig_name = immunesig_name_unify(ori_TMEsig_name)
names(ori_TMEsig) = ori_TMEsig_name$immune_sig



data =  as.matrix(log2(1+tpm_mat))

# gsva -----------------------------------------------------------------------
enrich_method = "ssgsea"
ITSgsva <- cor_immunesig_ITSgsva_ICBdata(data = data,
                              genelist_p = genelist_p,
                              genelist_n = genelist_n,
                              ori_TMEsig = ori_TMEsig,
                              type = "ITS", 
                              enrich_method = "ssgsea",
                              savepath = savepath,
                              savepath_dotplot = NULL)



xcell_remove <- xcell_sig_remover()


# Liu --------------------------------------------------------------------------
# compute ITS and significant test between R and NR ----------------------------
# generating immune signatures -------------------------------------------------
processedDataPath = "/picb/bigdata/project/CancerSysBio/Immunotherapy_prediction/data/RNAseq/"
dir.create(processedDataPath)
datasets <- list.files(processedDataPath)

# Liu_dataset_PD1 ------------------------------------------------------------------
dataset = "Liu_dataset.Rdata"
load(paste0(processedDataPath, "/Liu_dataset/processedData/", dataset))
response = response_CTLA4Naive_MAPKTxNaive_skin[c('patient', 'label')]
names(response) = c("patientID", "label")
tpm_mat <- tpm_mat[intersect(names(tpm_mat), response$patient)]
tpm_mat <- tpm_mat[-grep("Patient143", colnames(tpm_mat))]
response = response[is.element(response$patient, intersect(names(tpm_mat), response$patient)),]


result_dir = "05_immuneSig_ICBresponse/results_r0.4_ITSssgsea_noCombat/validation_results/Liu_dataset_CTLA4Naive_outcome/"
dir.create(result_dir)

cancer = "SKCM"
sample_type = "Metastatic"
# ------------------------------------------------------------------------------

savepath = result_dir
dir.create(paste0(savepath, "/ITSasPredictor_r0.4_oe_original/"))

savepath_allresult = "02_tumorGene_immuneSig_correlation/results_selected/ITS_oe_original_spearman_TUMERIC_r0.4/"
load(paste0(savepath_allresult, "for_GSVA/pcor/pcor_spearman_", cancer, "_",sample_type, "_positive_200.Rdata"))
load(paste0(savepath_allresult, "for_GSVA/pcor/pcor_spearman_", cancer, "_",sample_type, "_negative_200.Rdata"))
genelist_p_name = data.frame(immunesig = names(genelist_p), immune_sig = tolower(names(genelist_p)))
genelist_p_name = immunesig_name_unify(genelist_p_name)
names(genelist_p) = genelist_p_name$immune_sig

genelist_n_name = data.frame(immunesig = names(genelist_n), immune_sig = tolower(names(genelist_n)))
genelist_n_name = immunesig_name_unify(genelist_n_name)
names(genelist_n) = genelist_n_name$immune_sig



data =  as.matrix(log2(1+tpm_mat))

enrich_method = "ssgsea"
ITSp <- immune_sigAnalysis(data = data, 
                           gs = genelist_p,
                           method = enrich_method, 
                           kcdf = "Gaussian") 

ITSp_res <- as.data.frame(t(ITSp))
ITSp_res$patientID <- rownames(ITSp_res)
ITSp_res <- inner_join(response, ITSp_res, by = "patientID")
write.csv(ITSp_res, paste0(savepath, "/ITSasPredictor_r0.4_oe_original/ITSp_", enrich_method, ".csv"), 
          row.names =F, quote = F)


ITSn <- immune_sigAnalysis(data = data, 
                           gs = genelist_n,
                           method = enrich_method, 
                           kcdf = "Gaussian") 
ITSn_res <- as.data.frame(t(ITSn))
ITSn_res$patientID <- rownames(ITSn_res)
ITSn_res <- inner_join(response, ITSn_res)
write.csv(ITSn_res, paste0(savepath, "/ITSasPredictor_r0.4_oe_original/ITSn_", enrich_method, ".csv"), 
          row.names =F, quote = F)


# GSE115821 --------------------------------------------------------------------------
# compute ITS and significant test between R and NR ----------------------------
# generating immune signatures -------------------------------------------------
# dataset match : 05_immuneSig_ICBresponse/scripts/scripts_for_icb_datasets/dataset_reform.R
processedDataPath = "/picb/bigdata/project/CancerSysBio/Immunotherapy_prediction/data/"

dataset = "GSE115821_pre.Rdata"
load(paste0(processedDataPath, "RNAseq/GSE115821/processedData/",dataset))
tpm_mat <- tpm_mat_pre
response <- response_pre[c("patientID", "label")]
names(response) <- c("patientID", "label")
if(length(which(rowSums(tpm_mat) == 0)) > 0)tpm_mat <- tpm_mat[rowSums(tpm_mat) != 0, ]
if(length(which(rownames(tpm_mat)=="")) > 0)tpm_mat <- tpm_mat[-which(rownames(tpm_mat)==""), ]

result_dir = "05_immuneSig_ICBresponse/results_r0.4_ITSssgsea_noCombat/validation_results/GSE115821_pre_outcome/"
dir.create(result_dir)
cancer = "SKCM"
sample_type = "Metastatic" # NOT SURE

# ------------------------------------------------------------------------------

savepath = result_dir
dir.create(paste0(savepath, "/ITSasPredictor_r0.4_oe_original/"))

savepath_allresult = "02_tumorGene_immuneSig_correlation/results_selected/ITS_oe_original_spearman_TUMERIC_r0.4/"
load(paste0(savepath_allresult, "for_GSVA/pcor/pcor_spearman_", cancer, "_",sample_type, "_positive_200.Rdata"))
load(paste0(savepath_allresult, "for_GSVA/pcor/pcor_spearman_", cancer, "_",sample_type, "_negative_200.Rdata"))
genelist_p_name = data.frame(immunesig = names(genelist_p), immune_sig = tolower(names(genelist_p)))
genelist_p_name = immunesig_name_unify(genelist_p_name)
names(genelist_p) = genelist_p_name$immune_sig

genelist_n_name = data.frame(immunesig = names(genelist_n), immune_sig = tolower(names(genelist_n)))
genelist_n_name = immunesig_name_unify(genelist_n_name)
names(genelist_n) = genelist_n_name$immune_sig


data =  as.matrix(log2(1+tpm_mat))

enrich_method = "ssgsea"
ITSp <- immune_sigAnalysis(data = data, 
                           gs = genelist_p,
                           method = enrich_method, 
                           kcdf = "Gaussian") 

ITSp_res <- as.data.frame(t(ITSp))
ITSp_res$patientID <- rownames(ITSp_res)
ITSp_res <- inner_join(response, ITSp_res, by = "patientID")
write.csv(ITSp_res, paste0(savepath, "/ITSasPredictor_r0.4_oe_original/ITSp_", enrich_method, ".csv"), 
          row.names =F, quote = F)


ITSn <- immune_sigAnalysis(data = data, 
                           gs = genelist_n,
                           method = enrich_method, 
                           kcdf = "Gaussian") 
ITSn_res <- as.data.frame(t(ITSn))
ITSn_res$patientID <- rownames(ITSn_res)
ITSn_res <- inner_join(response, ITSn_res)
write.csv(ITSn_res, paste0(savepath, "/ITSasPredictor_r0.4_oe_original/ITSn_", enrich_method, ".csv"), 
          row.names =F, quote = F)



# GSE96619 --------------------------------------------------------------------------
# compute ITS and significant test between R and NR ----------------------------
# generating immune signatures -------------------------------------------------
# dataset match : 05_immuneSig_ICBresponse/scripts/scripts_for_icb_datasets/dataset_reform.R
processedDataPath = "/picb/bigdata/project/CancerSysBio/Immunotherapy_prediction/data/"

dataset = "GSE96619_pre.Rdata"
load(paste0(processedDataPath, "RNAseq/GSE96619/processedData/",dataset))
tpm_mat <- tpm_pre
response <- response_pre
names(response) <- c("patientID", "label")
if(length(which(rowSums(tpm_mat) == 0)) > 0)tpm_mat <- tpm_mat[rowSums(tpm_mat) != 0, ]
if(length(which(rownames(tpm_mat)=="")) > 0)tpm_mat <- tpm_mat[-which(rownames(tpm_mat)==""), ]

result_dir = "05_immuneSig_ICBresponse/results_r0.4_ITSssgsea_noCombat/validation_results/GSE96619_outcome/"
dir.create(result_dir)
cancer = "SKCM"
sample_type = "Metastatic" # NOT SURE

# ------------------------------------------------------------------------------

savepath = result_dir
dir.create(paste0(savepath, "/ITSasPredictor_r0.4_oe_original/"))

savepath_allresult = "02_tumorGene_immuneSig_correlation/results_selected/ITS_oe_original_spearman_TUMERIC_r0.4/"
load(paste0(savepath_allresult, "for_GSVA/pcor/pcor_spearman_", cancer, "_",sample_type, "_positive_200.Rdata"))
load(paste0(savepath_allresult, "for_GSVA/pcor/pcor_spearman_", cancer, "_",sample_type, "_negative_200.Rdata"))
genelist_p_name = data.frame(immunesig = names(genelist_p), immune_sig = tolower(names(genelist_p)))
genelist_p_name = immunesig_name_unify(genelist_p_name)
names(genelist_p) = genelist_p_name$immune_sig

genelist_n_name = data.frame(immunesig = names(genelist_n), immune_sig = tolower(names(genelist_n)))
genelist_n_name = immunesig_name_unify(genelist_n_name)
names(genelist_n) = genelist_n_name$immune_sig

data =  as.matrix(log2(1+tpm_mat))

enrich_method = "ssgsea"
ITSp <- immune_sigAnalysis(data = data, 
                           gs = genelist_p,
                           method = enrich_method, 
                           kcdf = "Gaussian") 

ITSp_res <- as.data.frame(t(ITSp))
ITSp_res$patientID <- gsub("-", "_", rownames(ITSp_res))
ITSp_res <- inner_join(response, ITSp_res, by = "patientID")
write.csv(ITSp_res, paste0(savepath, "/ITSasPredictor_r0.4_oe_original/ITSp_", enrich_method, ".csv"), 
          row.names =F, quote = F)


ITSn <- immune_sigAnalysis(data = data, 
                           gs = genelist_n,
                           method = enrich_method, 
                           kcdf = "Gaussian") 
ITSn_res <- as.data.frame(t(ITSn))
ITSn_res$patientID <- rownames(ITSn_res)
ITSn_res <- inner_join(response, ITSn_res)
write.csv(ITSn_res, paste0(savepath, "/ITSasPredictor_r0.4_oe_original/ITSn_", enrich_method, ".csv"), 
          row.names =F, quote = F)


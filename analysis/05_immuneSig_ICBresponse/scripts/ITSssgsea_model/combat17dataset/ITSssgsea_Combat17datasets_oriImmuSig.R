# load data --------------------------------------------------------------------
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
library(stringr)


# Loading relied functions -----------------------------------------------------
# TIGS, AUC and Wilcox test
source("05_immuneSig_ICBresponse/scripts/functions/05_immunotherapy_immuneSig_analysis_function.R")
source("05_immuneSig_ICBresponse/scripts/functions/predictor_ITS.R")
source("06_results_summaryandplotting/scripts/immunesig_selection_function.R")
# Loading relied functions -----------------------------------------------------
# TIL estimation, TIDE, TIP, IPS and Immune Resistance Pragram
source("02_tumorGene_immuneSig_correlation/scripts/functions/02_immuneSig_generation_function.R")
source("02_tumorGene_immuneSig_correlation/scripts/functions/TILestimator.R")
source("02_tumorGene_immuneSig_correlation/scripts/functions/predictor_TIDEandTIP.R")
source("02_tumorGene_immuneSig_correlation/scripts/functions/predictor_IPS.R")
source("02_tumorGene_immuneSig_correlation/scripts/functions/predictor_icb_ICI_resistance_program.R")
source("02_tumorGene_immuneSig_correlation/scripts/functions/immuneSig_sig160.R")
source("02_tumorGene_immuneSig_correlation/scripts/functions/immuneSig_TCGA_caf_immuneCLSssgsea.R")
# TIGS, AUC and Wilcox test
source("05_immuneSig_ICBresponse/scripts/functions/05_immunotherapy_immuneSig_analysis_function.R")
source("05_immuneSig_ICBresponse/scripts/functions/05_immunotherapy_immuneSig_analysis_function.R")
source("05_immuneSig_ICBresponse/scripts/functions/predictor_ITS.R")

source("06_results_summaryandplotting/scripts/immunesig_selection_function.R")

# Load data -----------------------------------------------------
load("05_immuneSig_ICBresponse/results_r0.4_ITSssgsea_Combat/ICBdataset_Exp17datasets.Rdata")
# data, dataAnnot

for(i in 1:length(data)){
  if(i == 1){
    cm = rownames(data[[1]])
  } else{
    cm = intersect(cm,rownames(data[[i]]))
  }
}

data_bc <- lapply(data,function(z){ return(t(z[cm,]))})
data_bc <- do.call(rbind,data_bc)
dataAnnot_all <- do.call(rbind, dataAnnot)
dataAnnot_all$dataset <- str_split_fixed(rownames(dataAnnot_all),'\\.',2)[,1]
data_bc <- as.data.frame(data_bc)
data_bc$SampleID <- row.names(data_bc)
data_bc <- merge(data_bc, dataAnnot_all, by = 'SampleID')
data_bc$batch <- data_bc$dataset


rownames(data_bc) <- data_bc$SampleID
edata <- t(data_bc[-which(is.element(names(data_bc), c("SampleID", "patientID","label","cancer","sample_type","dataset","batch")))])
edata_log <- edata
# ----------------------------------------------------------
# ----------------------------------------------------------
# batch correction with combat  ------------------------------------------------
# ----------------------------------------------------------
# ----------------------------------------------------------
library(sva)

data_ac <- ComBat(dat = edata_log, 
                  batch = data_bc$batch)
data_ac <- as.data.frame(t(data_ac))

dir.create("05_immuneSig_ICBresponse/results_r0.4_ITSssgsea_Combat/oriImmuSig/")
dir.create("05_immuneSig_ICBresponse/results_r0.4_ITSssgsea_Combat/oriImmuSig_alldataset_combat/")
# ----------------------------------------------------------
# whole dataset
dset = "alldatasets"
dsetAnnot <- dataAnnot

cancer = unique(dsetAnnot$cancer)
sample_type = unique(dsetAnnot$sample_type)
response = dsetAnnot[c('SampleID', 'label')]

tpm_mat <- t(2^data_ac - 1)

result_dir = paste0("05_immuneSig_ICBresponse/results_r0.4_ITSssgsea_Combat/oriImmuSig_alldataset_combat/")
dir.create(result_dir)

# generate expr matrix for cibersort webserver ---------------------------------
dir.create(paste0(result_dir, "/processedData_forCIBERSORT/"))
dir.create(paste0(result_dir, "/processedData_forImmunCellAI/"))
dir.create(paste0(result_dir, "/processedData_forTIDE/"))

expr4cibersort(expr = tpm_mat,
               dataset = dset,
               savepath = paste0(result_dir, "/processedData_forCIBERSORT/"))
expr4immuncelai(expr = tpm_mat,
                dataset = dset,
                savepath = paste0(result_dir, "/processedData_forImmunCellAI/"))

expr4tide(expr = tpm_mat,
          dataset = dset,
          savepath = paste0(result_dir, "/processedData_forTIDE/"))

# Immunesig --------------------------------------------------------------------
sourceCodePath <- '/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/02_tumorGene_immuneSig_correlation/scripts/softwares/'
immuneSig_sig160(data = tpm_mat,
                 cancer = cancer,
                 sample_type = sample_type,
                 sourceCodePath = sourceCodePath,
                 savePath = result_dir)

sig_dir <- "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/02_tumorGene_immuneSig_correlation/data/immune_signatures/"
immuneSig_sig160_othersigs(data = tpm_mat,
                           cancer = cancer,
                           sample_type = sample_type,
                           method = "ssgsea",
                           kcdf =  "Gaussian", 
                           min.sz=1, 
                           max.sz=1000,
                           sig_dir = sig_dir,
                           savePath = result_dir)

sig_dir <- "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/02_tumorGene_immuneSig_correlation/data/immune_signatures/"
TCGA_caf_immuneCLS(data = tpm_mat,
                   cancer = cancer,
                   sample_type = sample_type, 
                   method = "ssgsea",
                   kcdf =  "Gaussian", 
                   min.sz=0, 
                   max.sz=1000,
                   sig_dir = sig_dir,
                   savePath = result_dir)


# ICB Predictor ----------------------------------------------------------------
## TIDE and TIP
model_icb_tide_TIP(mat = tpm_mat,
                   cancername = cancer,
                   sample_type = sample_type,
                   sig_dir = "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r2_drug_sensitize_icb_prediction/",
                   save_path = result_dir)

res = read.csv(paste0(result_dir, "TIDE_TIP_result/TIDE_TIP_result_",cancer,"_",sample_type,".csv" ), row.names=1)
names(response) = c('patient', 'label')
res = inner_join(response, res)

## IPS
model_icb_IPS(gene_expression = tpm_mat,
              cancer_type = cancer,
              sample_type = sample_type,
              sig_dir = "02_tumorGene_immuneSig_correlation/scripts/softwares/Immunophenogram-master/",
              save_path = result_dir)
## ICI_resistance_program
model_icb_ICI_resistance_program(gene_expression = tpm_mat,
                                 cancer_type = cancer,
                                 sample_type = sample_type,
                                 sig_dir="02_tumorGene_immuneSig_correlation/scripts/softwares/",
                                 save_path = result_dir)

# TIL estimation ---------------------------------------------------------------
res_TIL_all1 <- immunecell_estiamtor(exprmat = tpm_mat, cancer = cancer,
                                     sample_type = sample_type, estimator_tool = "immunedeconv",  method = "all",
                                     savepath = result_dir)
res_TIL_all2 <- immunecell_estiamtor(exprmat = as.matrix(tpm_mat), cancer = 'SKCM',
                                     sample_type = sample_type, estimator_tool = "ConsensusTME", method = "gsva",
                                     savepath = result_dir)

rownames(res_TIL_all2) <- paste0(rownames(res_TIL_all2),"_ConsensusTMEgsva")
res_TIL_all = rbind(do.call(rbind, res_TIL_all1), res_TIL_all2)
dir.create(paste0(result_dir, "TIL_estimation/merged/"))
save(res_TIL_all, file = paste0(result_dir,"TIL_estimation/merged/",  
                                cancer, "_", sample_type, "_TILestimation.Rdata"))





result_dir = paste0("05_immuneSig_ICBresponse/results_r0.4_ITSssgsea_Combat/oriImmuSig/",dset,"/")
dir.create(result_dir)
# merge TME signature result files -----------------------------------------------------------
immuneSig_sig160_res <- read.csv(paste0(result_dir, "sig160_result/zscore_output/", cancer, "_", sample_type, ".csv"),row.names=1)
immuneSig_sig160_res <- as.data.frame(t(immuneSig_sig160_res))
colnames(immuneSig_sig160_res) <- paste0(colnames(immuneSig_sig160_res), "_sig160")
immuneSig_sig160_res$patient <- rownames(immuneSig_sig160_res)

load(paste0(result_dir, "sig160ssgsea_result/scaled_output/ssgsea_", cancer, "_", sample_type, "_immueSig.Rdata"))
immuneSig_sig160_othersigs_res <- esOut_scaled
immuneSig_sig160_othersigs_res <- as.data.frame(t(immuneSig_sig160_othersigs_res))
immuneSig_sig160_othersigs_res$patient <- rownames(immuneSig_sig160_othersigs_res)

load(paste0(result_dir, "TCGA_caf_immuneCLS_result/scaled_output/ssgsea_", cancer, "_", sample_type, "_immueSig.Rdata"))
TCGA_caf_immuneCLS_res <- esOut_scaled
TCGA_caf_immuneCLS_res <- as.data.frame(t(TCGA_caf_immuneCLS_res))
TCGA_caf_immuneCLS_res$patient <- rownames(TCGA_caf_immuneCLS_res)

TIDE_TIP_res <- read.csv(paste0(result_dir, "TIDE_TIP_result/TIDE_TIP_result_", cancer, "_", sample_type, ".csv"),row.names=1)
TIDE_TIP_res <- TIDE_TIP_res[-which(is.element(colnames(TIDE_TIP_res), c("class", "TIP_class")))]
# TIDE_TIP_res <- TIDE_TIP_res[c(1,3,5,6,9,11,18)]
TIDE_TIP_res <- TIDE_TIP_res[c("patient","Dysfunction","Exclusion",
                               "TIDE", "CD8","IFNG","TIP_signature")]
names(TIDE_TIP_res) <- c("patient","Dysfunction_TIDE","Exclusion_TIDE",
                         "TIDE_TIDE", "CD8_TIDE","IFNG_TIDE","TIP_signature_TIP")

IPS_res <- read.csv(paste0(result_dir, "IPS_result/IPS_result_", cancer, "_", sample_type, ".csv"),row.names=1)
names(IPS_res) <- c("patient","MHC_IPS","EC_IPS","SC_IPS","CP_IPS","AZ_IPS","IPS_IPS")

load(paste0(result_dir, "sig_OE_results/sig_OE_", cancer, "_", sample_type, ".Rdata"))
ICIRP_res <- df$res
names(ICIRP_res) <- paste0("OE_", names(ICIRP_res))
ICIRP_res$patient <- rownames(ICIRP_res)

load(paste0(result_dir,"TIL_estimation/merged/", cancer, "_", sample_type, "_TILestimation.Rdata"))
TIL_res <- as.data.frame(t(res_TIL_all))
TIL_res$patient <- rownames(TIL_res)
names(TIL_res) <- gsub("timer.", "", names(TIL_res))
names(TIL_res) <- gsub("quantiseq.", "", names(TIL_res))
names(TIL_res) <- gsub("mcp_counter.", "", names(TIL_res))
names(TIL_res) <- gsub("epic.", "", names(TIL_res))
names(TIL_res) <- gsub("xcell.", "", names(TIL_res))

TIDEweb_res <- read.csv(paste0(result_dir, "TIDEweb/export.csv"), row.names=1)
colnames(TIDEweb_res) <- paste0(colnames(TIDEweb_res),"_TIDEweb")
TIDEweb_res$patient = row.names(TIDEweb_res)
TIDEweb_res <- TIDEweb_res[c("patient", paste0(c("MSI.Expr.Sig","Merck18","MDSC", "CAF", "TAM.M2"),"_TIDEweb"))]

immuncellai_res <- read.csv(paste0(result_dir, "ImmuCellAI/ImmuCellAI_abundance_result.txt"), sep = '\t', row.names=1)
colnames(immuncellai_res) <- paste0(colnames(immuncellai_res),"_immuncellai")
immuncellai_res$patient = row.names(immuncellai_res)

cibersortx_res <- read.csv(paste0(result_dir, "/cibersortx/CIBERSORTx_Results.csv"), row.names=1)
cibersortx_res <- cibersortx_res[-c(23:25)]
colnames(cibersortx_res) <- paste0(colnames(cibersortx_res),"_CIBERSORT")
cibersortx_res$patient = row.names(cibersortx_res)

# cibersortxABS_res <- read.csv(paste0(result_dir, "/cibersortx/CIBERSORTxABS_Results.csv"), row.names=1)
# cibersortxABS_res <- cibersortxABS_res[-c(23:25)]
# colnames(cibersortxABS_res) <- paste0(colnames(cibersortxABS_res),"_cibersortxABS")
# cibersortxABS_res$patient = row.names(cibersortxABS_res)

merged_res <- merge(TIDE_TIP_res,
                    merge(IPS_res, 
                          merge(ICIRP_res, 
                                merge(TIL_res, 
                                      merge(TIDEweb_res,
                                            merge(immuncellai_res,
                                                  cibersortx_res,
                                                  by = "patient"), 
                                            by = "patient"), 
                                      by = "patient"), 
                                by = "patient"), 
                          by = "patient"), 
                    by = "patient")

merged_res <- merge(merged_res,
                    merge(immuneSig_sig160_res,
                          merge(immuneSig_sig160_othersigs_res,
                                TCGA_caf_immuneCLS_res,
                                by = "patient"),
                          by = "patient"), 
                    by = "patient")

write.csv(merged_res, paste0(result_dir, "/", cancer, "_", sample_type, "_immuneSig_merged.csv"), quote = F, row.names = F)


# ----------------------------------------------------------
# seperate datasets
for(dset in names(dataAnnot)){
  dsetAnnot <- dataAnnot[[dset]]
  cancer = unique(dsetAnnot$cancer)
  sample_type = unique(dsetAnnot$sample_type)
  response = dsetAnnot[c('SampleID', 'label')]
  tpm_mat_log2 <- data_ac[is.element(rownames(data_ac),dsetAnnot$SampleID),]
  tpm_mat <- t(2^tpm_mat_log2 - 1)
  
  result_dir = paste0("05_immuneSig_ICBresponse/results_r0.4_ITSssgsea_Combat/oriImmuSig/",dset,"/")
  dir.create(result_dir)
  
  # generate expr matrix for cibersort webserver ---------------------------------
  dir.create(paste0(result_dir, "/processedData_forCIBERSORT/"))
  dir.create(paste0(result_dir, "/processedData_forImmunCellAI/"))
  dir.create(paste0(result_dir, "/processedData_forTIDE/"))
  
  expr4cibersort(expr = tpm_mat,
                 dataset = dset,
                 savepath = paste0(result_dir, "/processedData_forCIBERSORT/"))
  expr4immuncelai(expr = tpm_mat,
                  dataset = dset,
                  savepath = paste0(result_dir, "/processedData_forImmunCellAI/"))
  
  expr4tide(expr = tpm_mat,
            dataset = dset,
            savepath = paste0(result_dir, "/processedData_forTIDE/"))
  
  # Immunesig --------------------------------------------------------------------
  sourceCodePath <- '/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/02_tumorGene_immuneSig_correlation/scripts/softwares/'
  immuneSig_sig160(data = tpm_mat,
                   cancer = cancer,
                   sample_type = sample_type,
                   sourceCodePath = sourceCodePath,
                   savePath = result_dir)
  
  sig_dir <- "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/02_tumorGene_immuneSig_correlation/data/immune_signatures/"
  immuneSig_sig160_othersigs(data = tpm_mat,
                             cancer = cancer,
                             sample_type = sample_type,
                             method = "ssgsea",
                             kcdf =  "Gaussian", 
                             min.sz=1, 
                             max.sz=1000,
                             sig_dir = sig_dir,
                             savePath = result_dir)
  
  sig_dir <- "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/02_tumorGene_immuneSig_correlation/data/immune_signatures/"
  TCGA_caf_immuneCLS(data = tpm_mat,
                     cancer = cancer,
                     sample_type = sample_type, 
                     method = "ssgsea",
                     kcdf =  "Gaussian", 
                     min.sz=0, 
                     max.sz=1000,
                     sig_dir = sig_dir,
                     savePath = result_dir)
  
  
  # ICB Predictor ----------------------------------------------------------------
  ## TIDE and TIP
  model_icb_tide_TIP(mat = tpm_mat,
                     cancername = cancer,
                     sample_type = sample_type,
                     sig_dir = "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r2_drug_sensitize_icb_prediction/",
                     save_path = result_dir)
  
  res = read.csv(paste0(result_dir, "TIDE_TIP_result/TIDE_TIP_result_",cancer,"_",sample_type,".csv" ), row.names=1)
  names(response) = c('patient', 'label')
  res = inner_join(response, res)
  
  ## IPS
  model_icb_IPS(gene_expression = tpm_mat,
                cancer_type = cancer,
                sample_type = sample_type,
                sig_dir = "02_tumorGene_immuneSig_correlation/scripts/softwares/Immunophenogram-master/",
                save_path = result_dir)
  ## ICI_resistance_program
  model_icb_ICI_resistance_program(gene_expression = tpm_mat,
                                   cancer_type = cancer,
                                   sample_type = sample_type,
                                   sig_dir="02_tumorGene_immuneSig_correlation/scripts/softwares/",
                                   save_path = result_dir)
  
  # TIL estimation ---------------------------------------------------------------
  res_TIL_all1 <- immunecell_estiamtor(exprmat = tpm_mat, cancer = cancer,
                                       sample_type = sample_type, estimator_tool = "immunedeconv",  method = "all",
                                       savepath = result_dir)
  res_TIL_all2 <- immunecell_estiamtor(exprmat = as.matrix(tpm_mat), cancer = cancer,
                                       sample_type = sample_type, estimator_tool = "ConsensusTME", method = "gsva",
                                       savepath = result_dir)
  
  rownames(res_TIL_all2) <- paste0(rownames(res_TIL_all2),"_ConsensusTMEgsva")
  res_TIL_all = rbind(do.call(rbind, res_TIL_all1), res_TIL_all2)
  dir.create(paste0(result_dir, "TIL_estimation/merged/"))
  save(res_TIL_all, file = paste0(result_dir,"TIL_estimation/merged/",  
                                  cancer, "_", sample_type, "_TILestimation.Rdata"))
  
}

for(dset in names(dataAnnot)){
  dsetAnnot <- dataAnnot[[dset]]
  cancer = unique(dsetAnnot$cancer)
  sample_type = unique(dsetAnnot$sample_type)
  response = dsetAnnot[c('SampleID', 'label')]
  tpm_mat_log2 <- data_ac[is.element(rownames(data_ac),dsetAnnot$SampleID),]
  tpm_mat <- t(2^tpm_mat_log2 - 1)
  
  result_dir = paste0("05_immuneSig_ICBresponse/results_r0.4_ITSssgsea_Combat/oriImmuSig/",dset,"/")
  dir.create(result_dir)
  # merge TME signature result files -----------------------------------------------------------
  immuneSig_sig160_res <- read.csv(paste0(result_dir, "sig160_result/zscore_output/", cancer, "_", sample_type, ".csv"),row.names=1)
  immuneSig_sig160_res <- as.data.frame(t(immuneSig_sig160_res))
  colnames(immuneSig_sig160_res) <- paste0(colnames(immuneSig_sig160_res), "_sig160")
  immuneSig_sig160_res$patient <- rownames(immuneSig_sig160_res)
  
  load(paste0(result_dir, "sig160ssgsea_result/scaled_output/ssgsea_", cancer, "_", sample_type, "_immueSig.Rdata"))
  immuneSig_sig160_othersigs_res <- esOut_scaled
  immuneSig_sig160_othersigs_res <- as.data.frame(t(immuneSig_sig160_othersigs_res))
  immuneSig_sig160_othersigs_res$patient <- rownames(immuneSig_sig160_othersigs_res)
  
  load(paste0(result_dir, "TCGA_caf_immuneCLS_result/scaled_output/ssgsea_", cancer, "_", sample_type, "_immueSig.Rdata"))
  TCGA_caf_immuneCLS_res <- esOut_scaled
  TCGA_caf_immuneCLS_res <- as.data.frame(t(TCGA_caf_immuneCLS_res))
  TCGA_caf_immuneCLS_res$patient <- rownames(TCGA_caf_immuneCLS_res)
  
  TIDE_TIP_res <- read.csv(paste0(result_dir, "TIDE_TIP_result/TIDE_TIP_result_", cancer, "_", sample_type, ".csv"),row.names=1)
  TIDE_TIP_res <- TIDE_TIP_res[-which(is.element(colnames(TIDE_TIP_res), c("class", "TIP_class")))]
  # TIDE_TIP_res <- TIDE_TIP_res[c(1,3,5,6,9,11,18)]
  TIDE_TIP_res <- TIDE_TIP_res[c("patient","Dysfunction","Exclusion",
                                 "TIDE", "CD8","IFNG","TIP_signature")]
  names(TIDE_TIP_res) <- c("patient","Dysfunction_TIDE","Exclusion_TIDE",
                           "TIDE_TIDE", "CD8_TIDE","IFNG_TIDE","TIP_signature_TIP")
  
  IPS_res <- read.csv(paste0(result_dir, "IPS_result/IPS_result_", cancer, "_", sample_type, ".csv"),row.names=1)
  names(IPS_res) <- c("patient","MHC_IPS","EC_IPS","SC_IPS","CP_IPS","AZ_IPS","IPS_IPS")
  
  load(paste0(result_dir, "sig_OE_results/sig_OE_", cancer, "_", sample_type, ".Rdata"))
  ICIRP_res <- df$res
  names(ICIRP_res) <- paste0("OE_", names(ICIRP_res))
  ICIRP_res$patient <- rownames(ICIRP_res)
  
  load(paste0(result_dir,"TIL_estimation/merged/", cancer, "_", sample_type, "_TILestimation.Rdata"))
  TIL_res <- as.data.frame(t(res_TIL_all))
  TIL_res$patient <- rownames(TIL_res)
  names(TIL_res) <- gsub("timer.", "", names(TIL_res))
  names(TIL_res) <- gsub("quantiseq.", "", names(TIL_res))
  names(TIL_res) <- gsub("mcp_counter.", "", names(TIL_res))
  names(TIL_res) <- gsub("epic.", "", names(TIL_res))
  names(TIL_res) <- gsub("xcell.", "", names(TIL_res))
  
  TIDEweb_res <- read.csv(paste0(result_dir, "TIDEweb/export.csv"), row.names=1)
  colnames(TIDEweb_res) <- paste0(colnames(TIDEweb_res),"_TIDEweb")
  TIDEweb_res$patient = row.names(TIDEweb_res)
  TIDEweb_res <- TIDEweb_res[c("patient", paste0(c("MSI.Expr.Sig","Merck18","MDSC", "CAF", "TAM.M2"),"_TIDEweb"))]
  
  immuncellai_res <- read.csv(paste0(result_dir, "ImmuCellAI/ImmuCellAI_abundance_result.txt"), sep = '\t', row.names=1)
  colnames(immuncellai_res) <- paste0(colnames(immuncellai_res),"_immuncellai")
  immuncellai_res$patient = row.names(immuncellai_res)
  
  cibersortx_res <- read.csv(paste0(result_dir, "/cibersortx/CIBERSORTx_Results.csv"), row.names=1)
  cibersortx_res <- cibersortx_res[-c(23:25)]
  colnames(cibersortx_res) <- paste0(colnames(cibersortx_res),"_CIBERSORT")
  cibersortx_res$patient = row.names(cibersortx_res)
  
  # cibersortxABS_res <- read.csv(paste0(result_dir, "/cibersortx/CIBERSORTxABS_Results.csv"), row.names=1)
  # cibersortxABS_res <- cibersortxABS_res[-c(23:25)]
  # colnames(cibersortxABS_res) <- paste0(colnames(cibersortxABS_res),"_cibersortxABS")
  # cibersortxABS_res$patient = row.names(cibersortxABS_res)
  
  merged_res <- merge(TIDE_TIP_res,
                      merge(IPS_res, 
                            merge(ICIRP_res, 
                                  merge(TIL_res, 
                                        merge(TIDEweb_res,
                                              merge(immuncellai_res,
                                                    cibersortx_res,
                                                    by = "patient"), 
                                              by = "patient"), 
                                        by = "patient"), 
                                  by = "patient"), 
                            by = "patient"), 
                      by = "patient")
  
  merged_res <- merge(merged_res,
                      merge(immuneSig_sig160_res,
                            merge(immuneSig_sig160_othersigs_res,
                                  TCGA_caf_immuneCLS_res,
                                  by = "patient"),
                            by = "patient"), 
                      by = "patient")
  
  write.csv(merged_res, paste0(result_dir, "/", cancer, "_", sample_type, "_immuneSig_merged.csv"), quote = F, row.names = F)
  
}
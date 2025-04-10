# Tumor genes - infiltrated immune cell fraction
workpath <- "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/"
setwd(workpath)

options(stringsAsFactors = F)

# Checking and loading R packages -----------------------------------------

if (!require(readxl)) {
  install.packages('readxl')
}
if (!require(data.table)) {
  install.packages('data.table')
}
if (!require(dplyr)) {
  install.packages('dplyr')
}

# if (!require(TCGAbiolinks)) {
#   BiocManager::install("TCGAbiolinks")
# }
if (!require(SummarizedExperiment)) {
  BiocManager::install("SummarizedExperiment")
}
if (!require(GSVA)) {
  BiocManager::install("GSVA")
}
if (!require(limma)) {
  BiocManager::install("limma")
}
if (!require(ConsensusClusterPlus)) {
  BiocManager::install("ConsensusClusterPlus")
}

if (!require(immunedeconv)) {
  devtools::install_github("grst/immunedeconv")
}

library(ROCR)
library(reshape2)

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

# set result saving path -------------------------------------------------------
result_path <- "05_immuneSig_ICBresponse/results/"


# generating immune signatures -------------------------------------------------

processedDataPath = "05_immuneSig_ICBresponse/data/processedData/"
dir.create(processedDataPath)
datasets <- list.files(processedDataPath)

# Hugo_dataset -----------------------------------------------------------------

dataset = "01_Hugo_dataset_outcome.Rdata"
load(paste0(processedDataPath, "/", dataset))
result_dir = "05_immuneSig_ICBresponse/results/01_Hugo_dataset_outcome/"
dir.create(result_dir)
cancer = "SKCM"
sample_type = "Metastatic"
merged_res <- read.csv(paste0(result_dir, "/", cancer, "_", sample_type, "_immuneSig_merged.csv"))
dim(merged_res)

icb_compr(result = merged_res,
          response = response,
          resultdir = "05_immuneSig_ICBresponse/results/immunesig_icbResponse_test_results/",
          filename = "01_Hugo_dataset_test_result.csv")


# Riaz_dataset -----------------------------------------------------------------
dataset = "02_Riaz_dataset_outcome.Rdata"
load(paste0(processedDataPath, "/", dataset))
result_dir = "05_immuneSig_ICBresponse/results/02_Riaz_dataset_outcome/"
dir.create(result_dir)
cancer = "SKCM"
sample_type = "Metastatic"

merged_res <- read.csv(paste0(result_dir, "/", cancer, "_", sample_type, "_immuneSig_merged.csv"))
dim(merged_res)

icb_compr(result = merged_res,
          response = response,
          resultdir = "05_immuneSig_ICBresponse/results/immunesig_icbResponse_test_results/",
          filename = "02_Riaz_dataset_test_result.csv")

# Kim_dataset ------------------------------------------------------------------
dataset = "06_Kim_dataset_outcome.Rdata"
load(paste0(processedDataPath, "/", dataset))
result_dir = "05_immuneSig_ICBresponse/results/06_Kim_dataset_outcome/"
dir.create(result_dir)
cancer = "STAD"
sample_type = "Metastatic"

merged_res <- read.csv(paste0(result_dir, "/", cancer, "_", sample_type, "_immuneSig_merged.csv"))
dim(merged_res)
icb_compr(result = merged_res,
          response = response,
          resultdir = "05_immuneSig_ICBresponse/results/immunesig_icbResponse_test_results/",
          filename = "06_Kim_dataset_test_result.csv")

# Gide_dataset -----------------------------------------------------------------

dataset = "05_Gide_dataset_outcome.Rdata"
load(paste0(processedDataPath, "/", dataset))
result_dir = "05_immuneSig_ICBresponse/results/05_Gide_dataset_outcome/"
dir.create(result_dir)
cancer = "SKCM"
sample_type = "Metastatic"

merged_res <- read.csv(paste0(result_dir, "/", cancer, "_", sample_type, "_immuneSig_merged.csv"))
dim(merged_res)
icb_compr(result = merged_res,
          response = response,
          resultdir = "05_immuneSig_ICBresponse/results/immunesig_icbResponse_test_results/",
          filename = "05_Gide_dataset_test_result.csv")


 #  Lauss_melanoma_immunotherapy -----------------------------------------------
dataset = "08_Lauss_melanoma.Rdata"
load(paste0(processedDataPath, "/", dataset))
result_dir = "05_immuneSig_ICBresponse/results/08_Lauss_dataset_outcome/"
dir.create(result_dir)
cancer = "SKCM"
sample_type = "Metastatic"

merged_res <- read.csv(paste0(result_dir, "/", cancer, "_", sample_type, "_immuneSig_merged.csv"))
dim(merged_res)
icb_compr(result = merged_res,
          response = response,
          resultdir = "05_immuneSig_ICBresponse/results/immunesig_icbResponse_test_results/",
          filename = "08_Lauss_dataset_test_result.csv")

# IMvigor210_dataset ------------------------------------------------------------------
dataset = "07_IMvigor210_bladder.Rdata"
load(paste0(processedDataPath, "/", dataset))
result_dir = "05_immuneSig_ICBresponse/results/07_IMvigor210_bladder_immunotherapy/"
dir.create(result_dir)
cancer = "BLCA"
sample_type = "Metastatic"
merged_res <- read.csv(paste0(result_dir, "/", cancer, "_", sample_type, "_immuneSig_merged.csv"))
dim(merged_res)

icb_compr(result = merged_res,
          response = response,
          resultdir = "05_immuneSig_ICBresponse/results/immunesig_icbResponse_test_results/",
          filename = "07_IMvigor210_dataset_test_result.csv")


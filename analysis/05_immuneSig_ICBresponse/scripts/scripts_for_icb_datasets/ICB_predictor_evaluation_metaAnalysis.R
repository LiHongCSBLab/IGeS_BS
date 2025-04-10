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

# TIGS, AUC and Wilcox test
source("05_immuneSig_ICBresponse/scripts/functions/05_immunotherapy_immuneSig_analysis_function.R")

# set result saving path -------------------------------------------------------
result_path <- "05_immuneSig_ICBresponse/results/"
dir.create("05_immuneSig_ICBresponse/results/meta_analysis_files/")


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
          resultdir = "05_immuneSig_ICBresponse/results/meta_analysis_files/",
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
dim(response)
length(intersect(merged_res$patient, response$patientID))

icb_compr(result = merged_res,
          response = response,
          resultdir = "05_immuneSig_ICBresponse/results/meta_analysis_files/",
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
dim(response)
length(intersect(merged_res$patient, response$patientID))
icb_compr(result = merged_res,
          response = response,
          resultdir = "05_immuneSig_ICBresponse/results/meta_analysis_files/",
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
dim(response)
length(intersect(merged_res$patient, response$patientID))
icb_compr(result = merged_res,
          response = response,
          resultdir = "05_immuneSig_ICBresponse/results/meta_analysis_files/",
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
dim(response)
length(intersect(merged_res$patient, response$patientID))
icb_compr(result = merged_res,
          response = response,
          resultdir = "05_immuneSig_ICBresponse/results/meta_analysis_files/",
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
dim(response)
length(intersect(merged_res$patient, response$patientID))

icb_compr(result = merged_res,
          response = response,
          resultdir = "05_immuneSig_ICBresponse/results/meta_analysis_files/",
          filename = "07_IMvigor210_dataset_test_result.csv")



# Braun_dataset_PD1 ------------------------------------------------------------------

processedDataPath = "/picb/bigdata/project/CancerSysBio/Immunotherapy_prediction/data/RNAseq/"
dir.create(processedDataPath)
# Braun_dataset_PD1 ------------------------------------------------------------------
dataset = "Braun_dataset_PD1.Rdata"
load(paste0(processedDataPath, "/Braun_dataset/processedData/", dataset))

result_dir = "05_immuneSig_ICBresponse/results/Braun_dataset/"
dir.create(result_dir)
cancer = "KIRC"
sample_type = "Primary"
dim(response)
length(intersect(merged_res$patient, response$patientID))

merged_res <- read.csv(paste0(result_dir, "/", cancer, "_", sample_type, "_immuneSig_merged.csv"))
dim(merged_res)
dim(response)
length(intersect(merged_res$patient, response$patientID))
if(nrow(response[is.na(response$label), ])>0)response <- response[!is.na(response$label), ]
icb_compr(result = merged_res,
          response = response,
          resultdir = "05_immuneSig_ICBresponse/results/meta_analysis_files/",
          filename = "Braun_dataset_PD1_test_result.csv")


# GSE126044 ------------------------------------------------------------------
processedDataPath = "/picb/bigdata/project/CancerSysBio/Immunotherapy_prediction/data/RNAseq/"
dir.create(processedDataPath)
datasets <- list.files(processedDataPath)
dataset = "GSE126044.Rdata"
load(paste0(processedDataPath, "/GSE126044/processedData/", dataset))
result_dir = "05_immuneSig_ICBresponse/results/GSE126044/"
dir.create(result_dir)
cancer = "LUAD"
sample_type = "Primary"

merged_res <- read.csv(paste0(result_dir, "/", cancer, "_", sample_type, "_immuneSig_merged.csv"))
dim(merged_res)
dim(response)
length(intersect(merged_res$patient, response$patientID))

icb_compr(result = merged_res,
          response = response,
          resultdir = "05_immuneSig_ICBresponse/results/meta_analysis_files/",
          filename = "GSE126044_test_result.csv")

# GSE135222 ------------------------------------------------------------------
# generating immune signatures -------------------------------------------------

processedDataPath = "/picb/bigdata/project/CancerSysBio/Immunotherapy_prediction/data/RNAseq/"
dir.create(processedDataPath)
datasets <- list.files(processedDataPath)

dataset = "GSE135222.Rdata"
load(paste0(processedDataPath, "/GSE135222/processedData/", dataset))
result_dir = "05_immuneSig_ICBresponse/results/GSE135222/"
dir.create(result_dir)
cancer = "LUAD"
sample_type = "Primary"

merged_res <- read.csv(paste0(result_dir, "/", cancer, "_", sample_type, "_immuneSig_merged.csv"))
dim(merged_res)
dim(response)
length(intersect(merged_res$patient, response$patientID))
icb_compr(result = merged_res,
          response = response,
          resultdir = "05_immuneSig_ICBresponse/results/meta_analysis_files/",
          filename = "GSE135222_test_result.csv")

# GSE145996 ------------------------------------------------------------------

processedDataPath = "/picb/bigdata/project/CancerSysBio/Immunotherapy_prediction/data/RNAseq/"
dir.create(processedDataPath)
datasets <- list.files(processedDataPath)
dataset = "GSE145996.Rdata"
load(paste0(processedDataPath, "/GSE145996/processedData/", dataset))
result_dir = "05_immuneSig_ICBresponse/results/GSE145996/"
dir.create(result_dir)
cancer = "SKCM"
sample_type = "Metastatic"

merged_res <- read.csv(paste0(result_dir, "/", cancer, "_", sample_type, "_immuneSig_merged.csv"))
dim(merged_res)
dim(response)
length(intersect(merged_res$patient, response$patientID))
icb_compr(result = merged_res,
          response = response,
          resultdir = "05_immuneSig_ICBresponse/results/meta_analysis_files/",
          filename = "GSE145996_test_result.csv")

# Nathanson_dataset ---------------------------------------------------------
# dataset match : 05_immuneSig_ICBresponse/scripts/scripts_for_icb_datasets/dataset_reform.R
processedDataPath = "/picb/bigdata/project/CancerSysBio/Immunotherapy_prediction/data/"

dataset = "Nathanson_dataset.Rdata"
load(paste0(processedDataPath, "RNAseq/Nathanson_dataset/processedData/", dataset))
result_dir = "05_immuneSig_ICBresponse/results/Nathanson_dataset_outcome/"
dir.create(result_dir)
cancer = "SKCM"
sample_type = "Metastatic"
merged_res <- read.csv(paste0(result_dir, "/", cancer, "_", sample_type, "_immuneSig_merged.csv"))
dim(merged_res)
dim(response)
length(intersect(merged_res$patient, response$patientID))
icb_compr(result = merged_res,
          response = response,
          resultdir = "05_immuneSig_ICBresponse/results/meta_analysis_files/",
          filename = "Nathanson_dataset_test_result.csv")

# PMID_29301960 ------------------------------------------------------------------
# generating immune signatures -------------------------------------------------

processedDataPath = "/picb/bigdata/project/CancerSysBio/Immunotherapy_prediction/data/RNAseq/"
dir.create(processedDataPath)
datasets <- list.files(processedDataPath)
dataset = "PMID_29301960.Rdata"
load(paste0(processedDataPath, "/PMID_29301960/processedData/", dataset))
result_dir = "05_immuneSig_ICBresponse/results/PMID_29301960/"
dir.create(result_dir)
cancer = "KIRC"
sample_type = "Primary"

merged_res <- read.csv(paste0(result_dir, "/", cancer, "_", sample_type, "_immuneSig_merged.csv"))
dim(merged_res)
dim(response)
length(intersect(merged_res$patient, response$patientID))
icb_compr(result = merged_res,
          response = response,
          resultdir = "05_immuneSig_ICBresponse/results/meta_analysis_files/",
          filename = "PMID_29301960_test_result.csv")

# Zhao_dataset ------------------------------------------------------------------
# generating immune signatures -------------------------------------------------

processedDataPath = "/picb/bigdata/project/CancerSysBio/Immunotherapy_prediction/data/RNAseq/"
dir.create(processedDataPath)
datasets <- list.files(processedDataPath)
dataset = "Zhao_dataset.Rdata"
load(paste0(processedDataPath, "/Zhao_dataset/processedData/", dataset))
result_dir = "05_immuneSig_ICBresponse/results/Zhao_dataset/"
dir.create(result_dir)
cancer = "GBM"
sample_type = "Primary"
merged_res <- read.csv(paste0(result_dir, "/", cancer, "_", sample_type, "_immuneSig_merged.csv"))

dim(merged_res)
dim(response)
length(intersect(merged_res$patient, response$patientID))

icb_compr(result = merged_res,
          response = response,
          resultdir = "05_immuneSig_ICBresponse/results/meta_analysis_files/",
          filename = "Zhao_dataset_test_result.csv")

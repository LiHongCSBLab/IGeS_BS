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

if (!require(TCGAbiolinks)) {
  BiocManager::install("TCGAbiolinks")
}
if (!require(SummarizedExperiment)) {
  BiocManager::install("SummarizedExperiment")
}

if (!require(limma)) {
  BiocManager::install("limma")
}


# Loading relied functions -----------------------------------------------------

source("02_tumorGene_immuneSig_correlation/scripts/functions/01_TCGA_RNAseq_dataTransformation_function.R")

rawDataPath <- "/picb/bigdata/project/FengFYM/s6_TCGAbiolinks_download_normalization/"
processedDataPath <- '02_tumorGene_immuneSig_correlation/data/TCGA_processedData/'

# Input cancer type ------------------------------------------------------------
project_ids <- c( "LAML","ACC","BLCA","LGG","BRCA","CESC","CHOL","COAD",
                  "ESCA","GBM","HNSC","KICH","KIRC","KIRP","LIHC",
                 "LUAD","LUSC",
                 "DLBC","MESO","OV","PAAD","PCPG","PRAD","READ","SARC","SKCM",
                 "STAD","TGCT","THYM","THCA","UCS","UCEC","UVM")

for(cancer in project_ids){
  if(cancer == 'SKCM'){
    TCGA_dataProcessing(cancer,
                        sample_type = 'Metastatic',
                        rawDataPath = rawDataPath,
                        processedDataPath = processedDataPath)
    TCGA_dataProcessing(cancer,
                        sample_type = 'Primary',
                        rawDataPath = rawDataPath,
                        processedDataPath = processedDataPath)
        
  }else{
    TCGA_dataProcessing(cancer,
                        sample_type = 'Primary',
                        rawDataPath = rawDataPath,
                        processedDataPath = processedDataPath)
  }
}

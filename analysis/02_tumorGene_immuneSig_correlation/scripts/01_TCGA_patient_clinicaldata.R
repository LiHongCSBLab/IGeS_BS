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

source("02_tumorGene_immuneSig_correlation/scripts/01_TCGA_RNAseq_dataTransformation_function.R")
 

# Input cancer type ------------------------------------------------------------
project_ids <- c("LGG","BRCA","CESC","CHOL","COAD", "ACC","BLCA", # "LAML",
                 "ESCA","GBM","HNSC","KICH","KIRC","KIRP","LIHC","LUAD","LUSC",
                 "SKCM", "PAAD","PCPG","PRAD","READ", "DLBC","MESO","STAD",
                 "TGCT","THYM","THCA","UCS","UCEC","OV","UVM","SARC")

for(cancer in project_ids){
  
  cancer_type = paste0('TCGA-', cancer)
  
  rawDataPath <- "/picb/bigdata/project/FengFYM/s6_TCGAbiolinks_download_normalization/"
  processedDataPath <- '02_tumorGene_immuneSig_correlation/data/TCGA_processedData/'
  processedDataPath <- paste0(processedDataPath, cancer,'/')
  dir.create(processedDataPath)

  sampleInfoPath = paste0(
    rawDataPath,
    "TCGA_sampleinfo/",
    cancer,
    ".Rdata"
  )
  
  if (!file.exists(sampleInfoPath)) {
    query <- GDCquery(
      project = cancer_type,
      data.category = "Transcriptome Profiling",
      data.type = "Gene Expression Quantification",
      workflow.type = "HTSeq - FPKM-UQ"
    )
    GDCdownload(query, method = "api", files.per.chunk = 200)
    data <- GDCprepare(query)
    sampleinfo <- colData(data)
    save(sampleinfo, file = sampleInfoPath)

  } 
}


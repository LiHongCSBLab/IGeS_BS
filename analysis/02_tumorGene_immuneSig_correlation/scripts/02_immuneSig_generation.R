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


# Loading relied functions -----------------------------------------------------

# source("02_tumorGene_immuneSig_correlation/scripts/03_cor_gene_immuneSig_function.R")
source("02_tumorGene_immuneSig_correlation/scripts/functions/02_immuneSig_generation_function.R")
source("02_tumorGene_immuneSig_correlation/scripts/functions/TILestimator.R")
source("02_tumorGene_immuneSig_correlation/scripts/functions/predictor_TIDEandTIP.R")
source("02_tumorGene_immuneSig_correlation/scripts/functions/predictor_IPS.R")
source("02_tumorGene_immuneSig_correlation/scripts/functions/predictor_icb_ICI_resistance_program.R")
source("02_tumorGene_immuneSig_correlation/scripts/functions/immuneSig_sig160.R")
source("02_tumorGene_immuneSig_correlation/scripts/functions/immuneSig_TCGA_caf_immuneCLSssgsea.R")


project_ids <- c("BRCA","COAD", "LIHC","LUAD","LUSC","DLBC",
                  "OV","PAAD","PRAD","READ",
                  "SKCM","UCEC","ACC","BLCA","ESCA","GBM","HNSC","KICH","KIRC","KIRP",
                  "LGG","MESO","PCPG","SARC",
                  "CESC","CHOL","STAD","TGCT","THYM","THCA","UCS","UVM", "LAML")

# Sig160 -----------------------------------------------------------------------
processedDataPath <- '02_tumorGene_immuneSig_correlation/data/TCGA_processedData/'
sourceCodePath <- '/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/02_tumorGene_immuneSig_correlation/scripts/softwares/'
originalSigPath <- '/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/02_tumorGene_immuneSig_correlation/data/immune_signatures/'
savePath <- "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/02_tumorGene_immuneSig_correlation/data/immune_signatures/"
sample_type = 'Primary'

for(cancer in project_ids){
  # cancer = 'LIHC'
  immuneSig_sig160_consistency(cancer = cancer,
                               sample_type = sample_type,
                               processedDataPath = processedDataPath,
                               sourceCodePath = sourceCodePath,
                               originalSigPath = originalSigPath, 
                               savePath = savePath)
}
cancer = 'SKCM'
sample_type = 'Metastatic'

immuneSig_sig160_consistency(cancer = cancer,
                             sample_type = sample_type,
                             processedDataPath = processedDataPath,
                             sourceCodePath = sourceCodePath,
                             originalSigPath = originalSigPath, 
                             savePath = savePath)

# sig160 - other immunesig using ssgsea ----------------------------------------
processedDataPath <- '02_tumorGene_immuneSig_correlation/data/TCGA_processedData/'
sig_dir <- "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/02_tumorGene_immuneSig_correlation/data/immune_signatures/"
originalSigPath <- '/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/02_tumorGene_immuneSig_correlation/data/immune_signatures/'
savePath <- "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/02_tumorGene_immuneSig_correlation/data/immune_signatures/"


sample_type = 'Primary'

for(cancer in project_ids){
  # cancer = 'LIHC'
  immuneSig_sig160_othersigs_consistency(cancer = cancer,
                               sample_type = sample_type,
                               method = "ssgsea",
                                 kcdf =  "Gaussian", 
                                 min.sz = 1, 
                                 max.sz = 1000,
                                 sig_dir = sig_dir,
                                 processedDataPath = processedDataPath,
                                 originalSigPath = originalSigPath, 
                                 savePath = savePath)
}
cancer = 'SKCM'
sample_type = 'Metastatic'

immuneSig_sig160_othersigs_consistency(cancer = cancer,
                             sample_type = sample_type,
                             method = "ssgsea",
                             kcdf =  "Gaussian", 
                             min.sz = 1, 
                             max.sz = 1000,
                             sig_dir = sig_dir,
                             processedDataPath = processedDataPath,
                             originalSigPath = originalSigPath, 
                             savePath = savePath)
                                 
# TCGA_caf_immuneCLS -----------------------------------------------------------
processedDataPath <- '02_tumorGene_immuneSig_correlation/data/TCGA_processedData/'
sig_dir <- "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/02_tumorGene_immuneSig_correlation/data/immune_signatures/"
originalSigPath <- "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/02_tumorGene_immuneSig_correlation/data/immune_signatures/"
savePath <- "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/02_tumorGene_immuneSig_correlation/data/immune_signatures/"


sample_type = 'Primary'

for (cancer in project_ids) {
  # cancer = 'LIHC'
  TCGA_caf_immuneCLS_consistency(
    cancer = cancer,
    sample_type = sample_type,
    method = "ssgsea",
    kcdf = "Gaussian",
    min.sz = 1,
    max.sz = 1000,
    sig_dir = sig_dir,
    processedDataPath = processedDataPath,
    originalSigPath = originalSigPath,
    savePath = savePath
  )
}

cancer = 'SKCM'
sample_type = 'Metastatic'
TCGA_caf_immuneCLS_consistency(cancer = cancer,
                               sample_type = sample_type, 
                               method = "ssgsea",
                               kcdf =  "Gaussian", 
                               min.sz = , 
                               max.sz = 1000,
                               sig_dir = sig_dir,
                               processedDataPath = processedDataPath,
                               originalSigPath = originalSigPath, 
                               savePath = savePath)



# immunotherapy predictor and TIL estimation -----------------------------------

# SKCM-Metastatic

cancer = 'SKCM'
sample_type = 'Metastatic'
cancer_type = paste0('TCGA-', cancer)

processedDataPath <- '02_tumorGene_immuneSig_correlation/data/TCGA_processedData/'
processedDataPath <- paste0(processedDataPath, cancer,'/')
dir.create(processedDataPath)


# Prepare TCGA expression matrix --------------------------------------------

# load in full expression matrix
# drop out low expression genes
# save to certain path

TPMMatPath = paste0(processedDataPath, "mixture_file_TCGA_", cancer,"_",sample_type, "_TPM.txt")
TPM_filtered <- read.csv(TPMMatPath, sep = '\t', row.names = 1)
colnames(TPM_filtered) <- gsub('\\.', '-',   colnames(TPM_filtered))

# 1) Calculate the rest immune cell and immune signatures

result_dir = "02_tumorGene_immuneSig_correlation/data/immune_signatures/"
dir.create(paste0(result_dir,"TIL_estimation/"))
dir.create(paste0(result_dir,"TIL_estimation/merged/"))
model_icb_tide_TIP(mat = TPM_filtered,
                   cancername = cancer,
                   sample_type = sample_type,
                   sig_dir = "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r2_drug_sensitize_icb_prediction/",
                   save_path = result_dir)

model_icb_IPS(gene_expression = TPM_filtered,
              cancer_type = cancer,
              sample_type = sample_type,
              sig_dir = "02_tumorGene_immuneSig_correlation/scripts/softwares/Immunophenogram-master/",
              save_path = result_dir)


# TIL estimation ---------------------------------------------------------------

res_TIL_all1 <- immunecell_estiamtor(
  exprmat = TPM_filtered,
  cancer = cancer,
  sample_type = sample_type,
  estimator_tool = "immunedeconv",
  method = "all",
  savepath = result_dir
)
res_TIL_all2 <- immunecell_estiamtor(
  exprmat = as.matrix(TPM_filtered),
  cancer = cancer,
  sample_type = sample_type,
  estimator_tool = "ConsensusTME",
  method = "gsva",
  savepath = result_dir
)

rownames(res_TIL_all2) <- paste0(rownames(res_TIL_all2),"_ConsensusTMEgsva")
res_TIL_all = rbind(do.call(rbind, res_TIL_all1),
                    res_TIL_all2)
save(res_TIL_all, file = paste0(result_dir,"TIL_estimation/merged/TCGA_", 
                                cancer, "_", sample_type, "_TILestimation.Rdata"))


# Input cancer type ------------------------------------------------------------

project_ids <- c("BRCA","COAD", "LIHC","LUAD","LUSC","DLBC","OV","PAAD","PRAD","READ",
                 "SKCM","UCEC","ACC","BLCA","ESCA","GBM","HNSC","KICH","KIRC","KIRP",
                 "LGG","MESO","PCPG","SARC",
                 "CESC","CHOL","STAD","TGCT","THYM","THCA","UCS", "LAML","UVM")

for(cancer in project_ids){
  # cancer = 'ACC'
  sample_type = 'Primary'
  cancer_type = paste0('TCGA-', cancer)
  
  processedDataPath <- '02_tumorGene_immuneSig_correlation/data/TCGA_processedData/'
  processedDataPath <- paste0(processedDataPath, cancer,'/')
  dir.create(processedDataPath)
  
  
  # Prepare TCGA expression matrix --------------------------------------------
  
  # load in full expression matrix
  # drop out low expression genes
  # save to certain path
  
  TPMMatPath = paste0(processedDataPath, "mixture_file_TCGA_", cancer,"_",sample_type, "_TPM.txt")
  TPM_filtered <- read.csv(TPMMatPath, sep = '\t', row.names = 1)
  colnames(TPM_filtered) <- gsub('\\.', '-',   colnames(TPM_filtered))
  
  # 1) Calculate the rest immune cell and immune signatures
  
  result_dir = "02_tumorGene_immuneSig_correlation/data/immune_signatures/"
  dir.create(paste0(result_dir,"TIDE_TIP_result/"))
  dir.create(paste0(result_dir,"TIL_estimation/"))
  dir.create(paste0(result_dir,"TIL_estimation/merged/"))
  
  
  
  model_icb_tide_TIP(mat = TPM_filtered,
                     cancername = cancer,
                     sample_type = sample_type,
                     sig_dir = "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r2_drug_sensitize_icb_prediction/",
                     save_path = result_dir)
  
  model_icb_IPS(gene_expression = TPM_filtered,
                cancer_type = cancer,
                sample_type = sample_type,
                sig_dir = "02_tumorGene_immuneSig_correlation/scripts/softwares/Immunophenogram-master/",
                save_path = result_dir)
  
  
  # TIL estimation ---------------------------------------------------------------
  
  res_TIL_all1 <- immunecell_estiamtor(
    exprmat = TPM_filtered,
    cancer = cancer,
    sample_type = sample_type,
    estimator_tool = "immunedeconv",
    method = "all",
    savepath = result_dir
  )
  
  if(cancer == 'LAML'){
    
    res_TIL_all2 <- matrix(NA,  19, ncol(exprmat))
    colnames(res_TIL_all2) <- colnames(exprmat)
    rownames(res_TIL_all2) <- c("B_cells","Cytotoxic_cells","Dendritic_cells",
                                "Eosinophils","Macrophages","Mast_cells",
                                "NK_cells","Neutrophils","T_cells_CD4",
                                "T_cells_CD8","T_cells_gamma_delta",
                                "T_regulatory_cells", "Macrophages_M1",
                                "Macrophages_M2","Endothelial", "Fibroblasts",
                                "Monocytes","Plasma_cells","Immune_Score" )
    
  }else{
    res_TIL_all2 <- immunecell_estiamtor(
      exprmat = as.matrix(TPM_filtered),
      cancer = cancer,
      sample_type = sample_type,
      estimator_tool = "ConsensusTME",
      method = "gsva",
      savepath = result_dir
    )
  }
  rownames(res_TIL_all2) <- paste0(rownames(res_TIL_all2),"_ConsensusTMEgsva")
  res_TIL_all = rbind(do.call(rbind, res_TIL_all1),
                      res_TIL_all2)
  
  
  save(res_TIL_all, file = paste0(result_dir,"TIL_estimation/merged/TCGA_", 
                                  cancer, "_", sample_type, "_TILestimation.Rdata"))
  
}

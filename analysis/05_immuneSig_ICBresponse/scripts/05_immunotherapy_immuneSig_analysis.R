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

# generate expr matrix for cibersort webserver ---------------------------------
expr4cibersort(expr = tpm_mat,
               dataset = "01_Hugo_dataset",
               savepath = "05_immuneSig_ICBresponse/data/processedData_forCIBERSORT/")
expr4tide(expr = tpm_mat,
          dataset = "01_Hugo_dataset",
          savepath = "05_immuneSig_ICBresponse/data/processedData_forTIDE/")


model_icb_tide_TIP(mat = tpm_mat,
                   cancername = cancer,
                   sample_type = sample_type,
                   sig_dir = "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r2_drug_sensitize_icb_prediction/",
                   save_path = result_dir)

res = read.csv(paste0(result_dir, "TIDE_TIP_result/TIDE_TIP_result_SKCM_Metastatic.csv" ), row.names=1)
names(response) = c('patient', 'label')
res = inner_join(response, res)

model_icb_IPS(gene_expression = tpm_mat,
              cancer_type = cancer,
              sample_type = sample_type,
              sig_dir = "02_tumorGene_immuneSig_correlation/scripts/softwares/Immunophenogram-master/",
              save_path = result_dir)

model_icb_ICI_resistance_program(gene_expression = tpm_mat,
                                 cancer_type = cancer,
                                 sample_type = sample_type,
                                 sig_dir="02_tumorGene_immuneSig_correlation/scripts/softwares/",
                                 save_path = result_dir)

# TIL estimation 
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

# merge result files
TIDE_TIP_res <- read.csv(paste0(result_dir, "TIDE_TIP_result/TIDE_TIP_result_", cancer, "_", sample_type, ".csv"),row.names=1)
TIDE_TIP_res <- TIDE_TIP_res[-which(is.element(colnames(TIDE_TIP_res), c("class", "TIP_class")))]
IPS_res <- read.csv(paste0(result_dir, "IPS_result/IPS_result_", cancer, "_", sample_type, ".csv"),row.names=1)
names(IPS_res) <- c("patient","MHC","EC","SC","CP","AZ","IPS")
load(paste0(result_dir, "sig_OE_results/sig_OE_", cancer, "_", sample_type, ".Rdata"))
ICIRP_res <- df$res
ICIRP_res$patient <- rownames(ICIRP_res)
load(paste0(result_dir,"TIL_estimation/merged/", cancer, "_", sample_type, "_TILestimation.Rdata"))
TIL_res <- as.data.frame(t(res_TIL_all))
TIL_res$patient <- rownames(TIL_res)

TIDEweb_res <- read.csv(paste0(result_dir, "TIDEweb/01_Hugo_dataset.csv"), row.names=1)
colnames(TIDEweb_res) <- paste0(colnames(TIDEweb_res),"_TIDEweb")
TIDEweb_res$patient = row.names(TIDEweb_res)
TIDEweb_res <- TIDEweb_res[c("patient", paste0(c("MSI.Expr.Sig","Merck18","MDSC", "CAF", "TAM.M2"),"_TIDEweb"))]

immuncellai_res <- read.csv(paste0(result_dir, "ImmuCellAI/ImmuCellAI_abundance_result.txt"), sep = '\t', row.names=1)
colnames(immuncellai_res) <- paste0(colnames(immuncellai_res),"_immuncellai")
immuncellai_res$patient = row.names(immuncellai_res)

cibersortx_res <- read.csv(paste0(result_dir, "/cibersortx/CIBERSORTx_Results.csv"), row.names=1)
cibersortx_res <- cibersortx_res[-c(23:25)]
colnames(cibersortx_res) <- paste0(colnames(cibersortx_res),"_cibersortx")
cibersortx_res$patient = row.names(cibersortx_res)

cibersortxABS_res <- read.csv(paste0(result_dir, "/cibersortx/CIBERSORTxABS_Results.csv"), row.names=1)
cibersortxABS_res <- cibersortxABS_res[-c(23:25)]
colnames(cibersortxABS_res) <- paste0(colnames(cibersortxABS_res),"_cibersortxABS")
cibersortxABS_res$patient = row.names(cibersortxABS_res)

merged_res <- merge(TIDE_TIP_res,
                merge(IPS_res, 
                merge(ICIRP_res, 
                merge(TIL_res, 
                merge(TIDEweb_res,
                merge(immuncellai_res,
                merge(cibersortx_res,
                      cibersortxABS_res,
                by = "patient"), 
                by = "patient"), 
                by = "patient"), 
                by = "patient"), 
                by = "patient"), 
                by = "patient"), 
                by = "patient")

write.csv(merged_res, paste0(result_dir, "/", cancer, "_", sample_type, "_immuneSig_merged.csv"), quote = F, row.names = F)


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

expr4cibersort(expr = tpm_mat,
               dataset = "02_Riaz_dataset",
               savepath = "05_immuneSig_ICBresponse/data/processedData_forCIBERSORT/")
expr4tide(expr = tpm_mat,
          dataset = "02_Riaz_dataset",
          savepath = "05_immuneSig_ICBresponse/data/processedData_forTIDE/")



model_icb_tide_TIP(mat = tpm_mat,
                   cancername = cancer,
                   sample_type = sample_type,
                   sig_dir = "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r2_drug_sensitize_icb_prediction/",
                   save_path = result_dir)


model_icb_IPS(gene_expression = tpm_mat,
              cancer_type = cancer,
              sample_type = sample_type,
              sig_dir = "02_tumorGene_immuneSig_correlation/scripts/softwares/Immunophenogram-master/",
              save_path = result_dir)

model_icb_ICI_resistance_program(gene_expression = tpm_mat,
                                 cancer_type = cancer,
                                 sample_type = sample_type,
                                 sig_dir="02_tumorGene_immuneSig_correlation/scripts/softwares/",
                                 save_path = result_dir)

# TIL estimation 
res_TIL_all1 <- immunecell_estiamtor(exprmat = tpm_mat, cancer = cancer,
  sample_type = sample_type, estimator_tool = "immunedeconv",  method = "all",
  savepath = result_dir)
res_TIL_all2 <- immunecell_estiamtor(exprmat = as.matrix(tpm_mat), cancer = cancer,
  sample_type = sample_type, estimator_tool = "ConsensusTME", method = "gsva",
  savepath = result_dir)

rownames(res_TIL_all2) <- paste0(rownames(res_TIL_all2),"_ConsensusTMEgsva")
res_TIL_all = rbind(do.call(rbind, res_TIL_all1),res_TIL_all2)
dir.create(paste0(result_dir, "TIL_estimation/merged/"))
save(res_TIL_all, file = paste0(result_dir,"TIL_estimation/merged/",  
                                cancer, "_", sample_type, "_TILestimation.Rdata"))


# merge result files
TIDE_TIP_res <- read.csv(paste0(result_dir, "TIDE_TIP_result/TIDE_TIP_result_", cancer, "_", sample_type, ".csv"),row.names=1)
TIDE_TIP_res <- TIDE_TIP_res[-which(is.element(colnames(TIDE_TIP_res), c("class", "TIP_class")))]
IPS_res <- read.csv(paste0(result_dir, "IPS_result/IPS_result_", cancer, "_", sample_type, ".csv"),row.names=1)
names(IPS_res) <- c("patient","MHC","EC","SC","CP","AZ","IPS")
load(paste0(result_dir, "sig_OE_results/sig_OE_", cancer, "_", sample_type, ".Rdata"))
ICIRP_res <- df$res
ICIRP_res$patient <- rownames(ICIRP_res)
load(paste0(result_dir,"TIL_estimation/merged/", cancer, "_", sample_type, "_TILestimation.Rdata"))
TIL_res <- as.data.frame(t(res_TIL_all))
TIL_res$patient <- rownames(TIL_res)

TIDEweb_res <- read.csv(paste0(result_dir, "TIDEweb/02_Riaz_dataset.csv"), row.names=1)
colnames(TIDEweb_res) <- paste0(colnames(TIDEweb_res),"_TIDEweb")
TIDEweb_res$patient = row.names(TIDEweb_res)
TIDEweb_res <- TIDEweb_res[c("patient", paste0(c("MSI.Expr.Sig","Merck18","MDSC", "CAF", "TAM.M2"),"_TIDEweb"))]

immuncellai_res <- read.csv(paste0(result_dir, "ImmuCellAI/ImmuCellAI_abundance_result.txt"), sep = '\t', row.names=1)
colnames(immuncellai_res) <- paste0(colnames(immuncellai_res),"_immuncellai")
immuncellai_res$patient = row.names(immuncellai_res)

cibersortx_res <- read.csv(paste0(result_dir, "/cibersortx/CIBERSORTx_Results.csv"), row.names=1)
cibersortx_res <- cibersortx_res[-c(23:25)]
colnames(cibersortx_res) <- paste0(colnames(cibersortx_res),"_cibersortx")
cibersortx_res$patient = row.names(cibersortx_res)

cibersortxABS_res <- read.csv(paste0(result_dir, "/cibersortx/CIBERSORTxABS_Results.csv"), row.names=1)
cibersortxABS_res <- cibersortxABS_res[-c(23:25)]
colnames(cibersortxABS_res) <- paste0(colnames(cibersortxABS_res),"_cibersortxABS")
cibersortxABS_res$patient = row.names(cibersortxABS_res)

merged_res <- merge(TIDE_TIP_res,
                merge(IPS_res, 
                merge(ICIRP_res, 
                merge(TIL_res, 
                merge(TIDEweb_res,
                merge(immuncellai_res,
                merge(cibersortx_res,
                cibersortxABS_res,
                by = "patient"), 
                by = "patient"), 
                by = "patient"), 
                by = "patient"), 
                by = "patient"), 
                by = "patient"), 
                by = "patient")
write.csv(merged_res, paste0(result_dir, "/", cancer, "_", sample_type, "_immuneSig_merged.csv"), quote = F, row.names = F)

merged_res <- read.csv(paste0(result_dir, "/", cancer, "_", sample_type, "_immuneSig_merged.csv"))
dim(merged_res)

icb_compr(result = merged_res,
          response = response,
          resultdir = "05_immuneSig_ICBresponse/results/immunesig_icbResponse_test_results/",
          filename = "02_Riaz_dataset_test_result.csv")

# Gide_dataset -----------------------------------------------------------------

dataset = "05_Gide_dataset_outcome.Rdata"
load(paste0(processedDataPath, "/", dataset))
result_dir = "05_immuneSig_ICBresponse/results/05_Gide_dataset_outcome/"
dir.create(result_dir)
cancer = "SKCM"
sample_type = "Metastatic"
expr4cibersort(expr = tpm_mat,
               dataset = "05_Gide_dataset",
               savepath = "05_immuneSig_ICBresponse/data/processedData_forCIBERSORT/")
expr4tide(expr = tpm_mat,
          dataset = "05_Gide_dataset",
          savepath = "05_immuneSig_ICBresponse/data/processedData_forTIDE/")


model_icb_tide_TIP(mat = tpm_mat,
                   cancername = cancer,
                   sample_type = sample_type,
                   sig_dir = "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r2_drug_sensitize_icb_prediction/",
                   save_path = result_dir)

model_icb_IPS(gene_expression = tpm_mat,
              cancer_type = cancer,
              sample_type = sample_type,
              sig_dir = "02_tumorGene_immuneSig_correlation/scripts/softwares/Immunophenogram-master/",
              save_path = result_dir)

model_icb_ICI_resistance_program(gene_expression = tpm_mat,
                                 cancer_type = cancer,
                                 sample_type = sample_type,
                                 sig_dir="02_tumorGene_immuneSig_correlation/scripts/softwares/",
                                 save_path = result_dir)

# TIL estimation 
res_TIL_all1 <- immunecell_estiamtor(exprmat = tpm_mat, cancer = cancer,
  sample_type = sample_type, estimator_tool = "immunedeconv",  method = "all",
  savepath = result_dir)
res_TIL_all2 <- immunecell_estiamtor(exprmat = as.matrix(tpm_mat), cancer = cancer,
  sample_type = sample_type, estimator_tool = "ConsensusTME", method = "gsva",
  savepath = result_dir)

rownames(res_TIL_all2) <- paste0(rownames(res_TIL_all2),"_ConsensusTMEgsva")
res_TIL_all = rbind(do.call(rbind, res_TIL_all1),res_TIL_all2)
dir.create(paste0(result_dir, "TIL_estimation/merged/"))
save(res_TIL_all, file = paste0(result_dir,"TIL_estimation/merged/",  
                                cancer, "_", sample_type, "_TILestimation.Rdata"))

# merge result files
TIDE_TIP_res <- read.csv(paste0(result_dir, "TIDE_TIP_result/TIDE_TIP_result_", cancer, "_", sample_type, ".csv"),row.names=1)
TIDE_TIP_res <- TIDE_TIP_res[-which(is.element(colnames(TIDE_TIP_res), c("class", "TIP_class")))]
IPS_res <- read.csv(paste0(result_dir, "IPS_result/IPS_result_", cancer, "_", sample_type, ".csv"),row.names=1)
names(IPS_res) <- c("patient","MHC","EC","SC","CP","AZ","IPS")
load(paste0(result_dir, "sig_OE_results/sig_OE_", cancer, "_", sample_type, ".Rdata"))
ICIRP_res <- df$res
ICIRP_res$patient <- rownames(ICIRP_res)
load(paste0(result_dir,"TIL_estimation/merged/", cancer, "_", sample_type, "_TILestimation.Rdata"))
TIL_res <- as.data.frame(t(res_TIL_all))
TIL_res$patient <- rownames(TIL_res)


TIDEweb_res <- read.csv(paste0(result_dir, "TIDEweb/05_Gide_dataset.csv"), row.names=1)
colnames(TIDEweb_res) <- paste0(colnames(TIDEweb_res),"_TIDEweb")
TIDEweb_res$patient = row.names(TIDEweb_res)
TIDEweb_res <- TIDEweb_res[c("patient", paste0(c("MSI.Expr.Sig","Merck18","MDSC", "CAF", "TAM.M2"),"_TIDEweb"))]

immuncellai_res <- read.csv(paste0(result_dir, "ImmuCellAI/ImmuCellAI_abundance_result.txt"), sep = '\t', row.names=1)
colnames(immuncellai_res) <- paste0(colnames(immuncellai_res),"_immuncellai")
immuncellai_res$patient = row.names(immuncellai_res)

cibersortx_res <- read.csv(paste0(result_dir, "/cibersortx/CIBERSORTx_Results.csv"), row.names=1)
cibersortx_res <- cibersortx_res[-c(23:25)]
colnames(cibersortx_res) <- paste0(colnames(cibersortx_res),"_cibersortx")
cibersortx_res$patient = row.names(cibersortx_res)

cibersortxABS_res <- read.csv(paste0(result_dir, "/cibersortx/CIBERSORTxABS_Results.csv"), row.names=1)
cibersortxABS_res <- cibersortxABS_res[-c(23:25)]
colnames(cibersortxABS_res) <- paste0(colnames(cibersortxABS_res),"_cibersortxABS")
cibersortxABS_res$patient = row.names(cibersortxABS_res)

merged_res <- merge(TIDE_TIP_res,
                merge(IPS_res, 
                merge(ICIRP_res, 
                merge(TIL_res, 
                merge(TIDEweb_res,
                merge(immuncellai_res,
                merge(cibersortx_res,
                cibersortxABS_res,
                by = "patient"), 
                by = "patient"), 
                by = "patient"), 
                by = "patient"), 
                by = "patient"), 
                by = "patient"), 
                by = "patient")

write.csv(merged_res, paste0(result_dir, "/", cancer, "_", sample_type, "_immuneSig_merged.csv"), quote = F, row.names = F)

merged_res <- read.csv(paste0(result_dir, "/", cancer, "_", sample_type, "_immuneSig_merged.csv"))
dim(merged_res)
icb_compr(result = merged_res,
          response = response,
          resultdir = "05_immuneSig_ICBresponse/results/immunesig_icbResponse_test_results/",
          filename = "05_Gide_dataset_test_result.csv")


# Kim_dataset ------------------------------------------------------------------
dataset = "06_Kim_dataset_outcome.Rdata"
load(paste0(processedDataPath, "/", dataset))
result_dir = "05_immuneSig_ICBresponse/results/06_Kim_dataset_outcome/"
dir.create(result_dir)
cancer = "STAD"
sample_type = "Metastatic"
expr4cibersort(expr = tpm_mat,
               dataset = "06_Kim_dataset",
               savepath = "05_immuneSig_ICBresponse/data/processedData_forCIBERSORT/")
expr4tide(expr = tpm_mat,
          dataset = "06_Kim_dataset",
          savepath = "05_immuneSig_ICBresponse/data/processedData_forTIDE/")

model_icb_tide_TIP(mat = tpm_mat,
                   cancername = cancer,
                   sample_type = sample_type,
                   sig_dir = "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r2_drug_sensitize_icb_prediction/",
                   save_path = result_dir)


model_icb_IPS(gene_expression = tpm_mat,
              cancer_type = cancer,
              sample_type = sample_type,
              sig_dir = "02_tumorGene_immuneSig_correlation/scripts/softwares/Immunophenogram-master/",
              save_path = result_dir)

model_icb_ICI_resistance_program(gene_expression = tpm_mat,
                                 cancer_type = cancer,
                                 sample_type = sample_type,
                                 sig_dir="02_tumorGene_immuneSig_correlation/scripts/softwares/",
                                 save_path = result_dir)

# TIL estimation 
res_TIL_all1 <- immunecell_estiamtor(exprmat = tpm_mat, cancer = cancer,
  sample_type = sample_type, estimator_tool = "immunedeconv",  method = "all",
  savepath = result_dir)
res_TIL_all2 <- immunecell_estiamtor(exprmat = as.matrix(tpm_mat), cancer = cancer,
  sample_type = sample_type, estimator_tool = "ConsensusTME", method = "gsva",
  savepath = result_dir)

rownames(res_TIL_all2) <- paste0(rownames(res_TIL_all2),"_ConsensusTMEgsva")
res_TIL_all = rbind(do.call(rbind, res_TIL_all1),res_TIL_all2)
dir.create(paste0(result_dir, "TIL_estimation/merged/"))
save(res_TIL_all, file = paste0(result_dir,"TIL_estimation/merged/",  
                                cancer, "_", sample_type, "_TILestimation.Rdata"))

# merge result files
TIDE_TIP_res <- read.csv(paste0(result_dir, "TIDE_TIP_result/TIDE_TIP_result_", cancer, "_", sample_type, ".csv"),row.names=1)
TIDE_TIP_res <- TIDE_TIP_res[-which(is.element(colnames(TIDE_TIP_res), c("class", "TIP_class")))]
IPS_res <- read.csv(paste0(result_dir, "IPS_result/IPS_result_", cancer, "_", sample_type, ".csv"),row.names=1)
names(IPS_res) <- c("patient","MHC","EC","SC","CP","AZ","IPS")
load(paste0(result_dir, "sig_OE_results/sig_OE_", cancer, "_", sample_type, ".Rdata"))
ICIRP_res <- df$res
ICIRP_res$patient <- rownames(ICIRP_res)
load(paste0(result_dir,"TIL_estimation/merged/", cancer, "_", sample_type, "_TILestimation.Rdata"))
TIL_res <- as.data.frame(t(res_TIL_all))
TIL_res$patient <- rownames(TIL_res)


TIDEweb_res <- read.csv(paste0(result_dir, "TIDEweb/06_Kim_dataset.csv"), row.names=1)
colnames(TIDEweb_res) <- paste0(colnames(TIDEweb_res),"_TIDEweb")
TIDEweb_res$patient = row.names(TIDEweb_res)
TIDEweb_res <- TIDEweb_res[c("patient", paste0(c("MSI.Expr.Sig","Merck18","MDSC", "CAF", "TAM.M2"),"_TIDEweb"))]

immuncellai_res <- read.csv(paste0(result_dir, "ImmuCellAI/ImmuCellAI_abundance_result.txt"), sep = '\t', row.names=1)
colnames(immuncellai_res) <- paste0(colnames(immuncellai_res),"_immuncellai")
immuncellai_res$patient = row.names(immuncellai_res)

cibersortx_res <- read.csv(paste0(result_dir, "/cibersortx/CIBERSORTx_Results.csv"), row.names=1)
cibersortx_res <- cibersortx_res[-c(23:25)]
colnames(cibersortx_res) <- paste0(colnames(cibersortx_res),"_cibersortx")
cibersortx_res$patient = row.names(cibersortx_res)

cibersortxABS_res <- read.csv(paste0(result_dir, "/cibersortx/CIBERSORTxABS_Results.csv"), row.names=1)
cibersortxABS_res <- cibersortxABS_res[-c(23:25)]
colnames(cibersortxABS_res) <- paste0(colnames(cibersortxABS_res),"_cibersortxABS")
cibersortxABS_res$patient = row.names(cibersortxABS_res)

merged_res <- merge(TIDE_TIP_res,
                merge(IPS_res, 
                merge(ICIRP_res, 
                merge(TIL_res, 
                merge(TIDEweb_res,
                merge(immuncellai_res,
                merge(cibersortx_res,
                cibersortxABS_res,
                by = "patient"), 
                by = "patient"), 
                by = "patient"), 
                by = "patient"), 
                by = "patient"), 
                by = "patient"), 
                by = "patient")

write.csv(merged_res, paste0(result_dir, "/", cancer, "_", sample_type, "_immuneSig_merged.csv"), quote = F, row.names = F)

merged_res <- read.csv(paste0(result_dir, "/", cancer, "_", sample_type, "_immuneSig_merged.csv"))
dim(merged_res)
icb_compr(result = merged_res,
          response = response,
          resultdir = "05_immuneSig_ICBresponse/results/immunesig_icbResponse_test_results/",
          filename = "06_Kim_dataset_test_result.csv")

 #  Lauss_melanoma_immunotherapy
dataset = "08_Lauss_melanoma.Rdata"
load(paste0(processedDataPath, "/", dataset))
result_dir = "05_immuneSig_ICBresponse/results/08_Lauss_dataset_outcome/"
dir.create(result_dir)
cancer = "SKCM"
sample_type = "Metastatic"

expr4cibersort(expr = tpm_mat,
               dataset = "08_Lauss_dataset",
               savepath = "05_immuneSig_ICBresponse/data/processedData_forCIBERSORT/")
expr4tide(expr = tpm_mat,
          dataset = "08_Lauss_dataset",
          savepath = "05_immuneSig_ICBresponse/data/processedData_forTIDE/")

model_icb_tide_TIP(mat = tpm_mat,
                   cancername = cancer,
                   sample_type = sample_type,
                   sig_dir = "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r2_drug_sensitize_icb_prediction/",
                   save_path = result_dir)


model_icb_IPS(gene_expression = tpm_mat,
              cancer_type = cancer,
              sample_type = sample_type,
              sig_dir = "02_tumorGene_immuneSig_correlation/scripts/softwares/Immunophenogram-master/",
              save_path = result_dir)

model_icb_ICI_resistance_program(gene_expression = tpm_mat,
                                 cancer_type = cancer,
                                 sample_type = sample_type,
                                 sig_dir="02_tumorGene_immuneSig_correlation/scripts/softwares/",
                                 save_path = result_dir)

# TIL estimation 
dir.create(paste0(result_dir, "TIL_estimation/"))
res_TIL_all1 <- immunecell_estiamtor(exprmat = tpm_mat, cancer = cancer,
  sample_type = sample_type, estimator_tool = "immunedeconv",  method = "all",
  savepath = result_dir)
res_TIL_all2 <- immunecell_estiamtor(exprmat = as.matrix(tpm_mat), cancer = cancer,
  sample_type = sample_type, estimator_tool = "ConsensusTME", method = "gsva",
  savepath = result_dir)

rownames(res_TIL_all2) <- paste0(rownames(res_TIL_all2),"_ConsensusTMEgsva")
res_TIL_all = rbind(do.call(rbind, res_TIL_all1),res_TIL_all2)
dir.create(paste0(result_dir, "TIL_estimation/merged/"))
save(res_TIL_all, file = paste0(result_dir,"TIL_estimation/merged/",  
                                cancer, "_", sample_type, "_TILestimation.Rdata"))

# merge result files
TIDE_TIP_res <- read.csv(paste0(result_dir, "TIDE_TIP_result/TIDE_TIP_result_", cancer, "_", sample_type, ".csv"),row.names=1)
TIDE_TIP_res <- TIDE_TIP_res[-which(is.element(colnames(TIDE_TIP_res), c("class", "TIP_class")))]
IPS_res <- read.csv(paste0(result_dir, "IPS_result/IPS_result_", cancer, "_", sample_type, ".csv"),row.names=1)
names(IPS_res) <- c("patient","MHC","EC","SC","CP","AZ","IPS")
load(paste0(result_dir, "sig_OE_results/sig_OE_", cancer, "_", sample_type, ".Rdata"))
ICIRP_res <- df$res
ICIRP_res$patient <- rownames(ICIRP_res)
load(paste0(result_dir,"TIL_estimation/merged/", cancer, "_", sample_type, "_TILestimation.Rdata"))
TIL_res <- as.data.frame(t(res_TIL_all))
TIL_res$patient <- rownames(TIL_res)

TIDEweb_res <- read.csv(paste0(result_dir, "TIDEweb/08_Lauss_dataset.csv"), row.names=1)
colnames(TIDEweb_res) <- paste0(colnames(TIDEweb_res),"_TIDEweb")
TIDEweb_res$patient = row.names(TIDEweb_res)
TIDEweb_res <- TIDEweb_res[c("patient", paste0(c("MSI.Expr.Sig","Merck18","MDSC", "CAF", "TAM.M2"),"_TIDEweb"))]

immuncellai_res <- read.csv(paste0(result_dir, "ImmuCellAI/ImmuCellAI_abundance_result.txt"), sep = '\t', row.names=1)
colnames(immuncellai_res) <- paste0(colnames(immuncellai_res),"_immuncellai")
immuncellai_res$patient = row.names(immuncellai_res)

cibersortx_res <- read.csv(paste0(result_dir, "/cibersortx/CIBERSORTx_Results.csv"), row.names=1)
cibersortx_res <- cibersortx_res[-c(23:25)]
colnames(cibersortx_res) <- paste0(colnames(cibersortx_res),"_cibersortx")
cibersortx_res$patient = row.names(cibersortx_res)

cibersortxABS_res <- read.csv(paste0(result_dir, "/cibersortx/CIBERSORTxABS_Results.csv"), row.names=1)
cibersortxABS_res <- cibersortxABS_res[-c(23:25)]
colnames(cibersortxABS_res) <- paste0(colnames(cibersortxABS_res),"_cibersortxABS")
cibersortxABS_res$patient = row.names(cibersortxABS_res)

merged_res <- merge(TIDE_TIP_res,
                merge(IPS_res, 
                merge(ICIRP_res, 
                merge(TIL_res, 
                merge(TIDEweb_res,
                merge(immuncellai_res,
                merge(cibersortx_res,
                cibersortxABS_res,
                by = "patient"), 
                by = "patient"), 
                by = "patient"), 
                by = "patient"), 
                by = "patient"), 
                by = "patient"), 
                by = "patient")


write.csv(merged_res, paste0(result_dir, "/", cancer, "_", sample_type, "_immuneSig_merged.csv"), quote = F, row.names = F)


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


expr4cibersort(expr = tpm_mat,
               dataset = "07_IMvigor210_dataset",
               savepath = "05_immuneSig_ICBresponse/data/processedData_forCIBERSORT/")
expr4tide(expr = tpm_mat,
          dataset = "07_IMvigor210_dataset",
          savepath = "05_immuneSig_ICBresponse/data/processedData_forTIDE/")

model_icb_tide_TIP(mat = tpm_mat,
                   cancername = cancer,
                   sample_type = sample_type,
                   sig_dir = "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r2_drug_sensitize_icb_prediction/",
                   save_path = result_dir)


model_icb_IPS(gene_expression = tpm_mat,
              cancer_type = cancer,
              sample_type = sample_type,
              sig_dir = "02_tumorGene_immuneSig_correlation/scripts/softwares/Immunophenogram-master/",
              save_path = result_dir)

model_icb_ICI_resistance_program(gene_expression = tpm_mat,
                                 cancer_type = cancer,
                                 sample_type = sample_type,
                                 sig_dir="02_tumorGene_immuneSig_correlation/scripts/softwares/",
                                 save_path = result_dir)

# TIL estimation 
dir.create(paste0(result_dir, "TIL_estimation/"))
res_TIL_all1 <- immunecell_estiamtor(exprmat = tpm_mat, 
                                     cancer = cancer, 
                                     sample_type = sample_type, 
                                     estimator_tool = "immunedeconv",  
                                     method = "all",
                                     savepath = result_dir)

res_TIL_all2 <- immunecell_estiamtor(exprmat = as.matrix(tpm_mat), 
                                     cancer = cancer,
                                     sample_type = sample_type, 
                                     estimator_tool = "ConsensusTME", 
                                     method = "gsva",
                                     savepath = result_dir)

rownames(res_TIL_all2) <- paste0(rownames(res_TIL_all2),"_ConsensusTMEgsva")
res_TIL_all = rbind(do.call(rbind, res_TIL_all1),res_TIL_all2)
dir.create(paste0(result_dir, "TIL_estimation/merged/"))
save(res_TIL_all, file = paste0(
  result_dir, "TIL_estimation/merged/",
  cancer, "_", sample_type, "_TILestimation.Rdata"
))

# merge result files
TIDE_TIP_res <- read.csv(paste0(result_dir, "TIDE_TIP_result/TIDE_TIP_result_", cancer, "_", sample_type, ".csv"),row.names=1)
TIDE_TIP_res <- TIDE_TIP_res[-which(is.element(colnames(TIDE_TIP_res), c("class", "TIP_class")))]
IPS_res <- read.csv(paste0(result_dir, "IPS_result/IPS_result_", cancer, "_", sample_type, ".csv"),row.names=1)
names(IPS_res) <- c("patient","MHC","EC","SC","CP","AZ","IPS")
load(paste0(result_dir, "sig_OE_results/sig_OE_", cancer, "_", sample_type, ".Rdata"))
ICIRP_res <- df$res
ICIRP_res$patient <- rownames(ICIRP_res)
load(paste0(result_dir,"TIL_estimation/merged/", cancer, "_", sample_type, "_TILestimation.Rdata"))
TIL_res <- as.data.frame(t(res_TIL_all))
TIL_res$patient <- rownames(TIL_res)


TIDEweb_res <- read.csv(paste0(result_dir, "TIDEweb/07_IMvigor210_dataset.csv"), row.names=1)
colnames(TIDEweb_res) <- paste0(colnames(TIDEweb_res),"_TIDEweb")
TIDEweb_res$patient = row.names(TIDEweb_res)
TIDEweb_res <- TIDEweb_res[c("patient", paste0(c("MSI.Expr.Sig","Merck18","MDSC", "CAF", "TAM.M2"),"_TIDEweb"))]

immuncellai_res <- read.csv(paste0(result_dir, "ImmuCellAI/ImmuCellAI_abundance_result.txt"), sep = '\t', row.names=1)
colnames(immuncellai_res) <- paste0(colnames(immuncellai_res),"_immuncellai")
immuncellai_res$patient = row.names(immuncellai_res)

cibersortx_res <- read.csv(paste0(result_dir, "/cibersortx/CIBERSORTx_Results.csv"), row.names=1)
cibersortx_res <- cibersortx_res[-c(23:25)]
colnames(cibersortx_res) <- paste0(colnames(cibersortx_res),"_cibersortx")
cibersortx_res$patient = row.names(cibersortx_res)

cibersortxABS_res <- read.csv(paste0(result_dir, "/cibersortx/CIBERSORTxABS_Results.csv"), row.names=1)
cibersortxABS_res <- cibersortxABS_res[-c(23:25)]
colnames(cibersortxABS_res) <- paste0(colnames(cibersortxABS_res),"_cibersortxABS")
cibersortxABS_res$patient = row.names(cibersortxABS_res)

merged_res <- merge(TIDE_TIP_res,
                merge(IPS_res, 
                merge(ICIRP_res, 
                merge(TIL_res, 
                merge(TIDEweb_res,
                merge(immuncellai_res,
                merge(cibersortx_res,
                cibersortxABS_res,
                by = "patient"), 
                by = "patient"), 
                by = "patient"), 
                by = "patient"), 
                by = "patient"), 
                by = "patient"), 
                by = "patient")

write.csv(merged_res, paste0(result_dir, "/", cancer, "_", sample_type, "_immuneSig_merged.csv"), quote = F, row.names = F)


merged_res <- read.csv(paste0(result_dir, "/", cancer, "_", sample_type, "_immuneSig_merged.csv"))
dim(merged_res)

icb_compr(result = merged_res,
          response = response,
          resultdir = "05_immuneSig_ICBresponse/results/immunesig_icbResponse_test_results/",
          filename = "07_IMvigor210_dataset_test_result.csv")


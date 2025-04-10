workpath <- "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/"
setwd(workpath)

options(stringsAsFactors = F)

# Checking and loading R packages -----------------------------------------

if (!require(readxl)) {  install.packages('readxl')}
if (!require(data.table)) {  install.packages('data.table')}
if (!require(dplyr)) {  install.packages('dplyr')}
# if (!require(TCGAbiolinks)) {  BiocManager::install("TCGAbiolinks") }
if (!require(SummarizedExperiment)) {  BiocManager::install("SummarizedExperiment")}
if (!require(GSVA)) {  BiocManager::install("GSVA")}
if (!require(limma)) {  BiocManager::install("limma")}
if (!require(ConsensusClusterPlus)) {  BiocManager::install("ConsensusClusterPlus")}
if (!require(immunedeconv)) {  devtools::install_github("grst/immunedeconv")}
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
source("05_immuneSig_ICBresponse/scripts/functions/predictor_ITS.R")

source("06_results_summaryandplotting/scripts/immunesig_selection_function.R")




library(mlr)
library(h2o)
set.seed(1234)

cores_method = "Slurm"
Ncores <- future::availableCores(method = cores_method)
options(mc.cores = Ncores)

load("05_immuneSig_ICBresponse/oriImmuneSig_res_scaled.Rdata")
mat <- do.call(rbind, merged_sig)
mat <- mat[!is.na(mat$label), ]
mat$tag <- 1
mat[mat$label == "NR",]$tag = 0
mat$tag <- factor(mat$tag, levels = c(0,1), order=TRUE)

mat <- mat[-c(1:5)]


task = makeClassifTask(id = "features", data = mat, target = "tag")


P_num = task$task.desc$class.distribution[[1]]
N_num = task$task.desc$class.distribution[[2]]

rdesc = makeResampleDesc("RepCV", folds = 10, reps = 10)

# elastic net from glmnet
elasticnet_opt_ps = elasticnet_tuning(task)

# set up model with optimal parameters
mod_elasticnet_opt = setHyperPars(learner = makeLearner("classif.cvglmnet", 
                                                        predict.type = "prob", 
                                                        fix.factors.prediction = TRUE,
                                                        nlambda = 1000L,
                                                        lambda.min.ratio = 1e-5,
                                                        nfolds = 5,
                                                        config = list(on.learner.error = "warn")) ,
                                  standardize =  elasticnet_opt_ps$standardize,
                                  s =  elasticnet_opt_ps$s,
                                  alpha =  elasticnet_opt_ps$alpha
)

mod_elasticnet_opt = makeUndersampleWrapper(mod_elasticnet_opt,  
                                            usw.rate = 1/(N_num/P_num))

# 10-fold cross validation to access the model performance
r_elasticnet = resample(mod_elasticnet_opt, 
                        task, 
                        rdesc, 
                        measures = list(mmce, tpr, fnr, fpr, tnr, acc, auc, f1, timetrain))

save(elasticnet_opt_ps,mod_elasticnet_opt,r_elasticnet,file = 'elasticnet_result_0.2.Rdata')


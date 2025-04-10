workpath <- "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/"
setwd(workpath)
set.seed(1234)

options(stringsAsFactors = F)

# Checking and loading R packages -----------------------------------------

if (!require(readxl)) {  install.packages('readxl')}
if (!require(data.table)) {  install.packages('data.table')}
if (!require(dplyr)) {  install.packages('dplyr')}
library(ROCR)
library(reshape2)
library(getopt)


command=matrix(c( 
  "featureFilterStrategy",  "s", 1, "character", "featureFilterStrategy in metaAnalysis, combinedP, voting, mlrFilter",
  "cancer",                 "c", 2, "character", "cancer",
  "mlrFilter",              "f", 2, "character", "mlrFilter in mlr",
  "mlrFilter_k",            "k", 2, "numeric",   "top k persent feature, k from 0-1",
  "dataset",                "d", 2, "character", "validation dataset",
  "help",                   "h", 0, "logical",   "help file"),
  byrow=T,ncol=5)

args=getopt(command)

featureFilterStrategy = args$featureFilterStrategy
cancer = args$cancer
mlrFilter = args$mlrFilter
mlrFilter_k = args$mlrFilter_k
dataset = args$dataset
# featureFilterStrategy = "mlrFilter" 
# mlrFilter =  'ranger_impurity'
# mlrFilter_k = 0.1
# dataset = "Liu_dataset_CTLA4Naive_outcome"

# Loading relied functions -----------------------------------------------------

# mlr --------------------------------------------------------------------------
source("05_immuneSig_ICBresponse/scripts/mlr_function/rf_learner.R")
source("05_immuneSig_ICBresponse/scripts/mlr_function/elasticnet_learner.R")
source("05_immuneSig_ICBresponse/scripts/mlr_function/ksvm_learner.R")
source("05_immuneSig_ICBresponse/scripts/mlr_function/ridge_learner.R")
source("05_immuneSig_ICBresponse/scripts/mlr_function/xgboost_learner.R")
source("05_immuneSig_ICBresponse/scripts/oriImmuSig_mlr_cancer/feature_filter_function.R")

source("06_results_summaryandplotting/scripts/immunesig_selection_function.R")

library(mlr)
library(mlrMBO)

# 
# ------------------------------------------------------------------------------
# feature selection ------———---------------------------------------------------
# and create directory according to strategy and model -------------------------
# ------------------------------------------------------------------------------

resultpath = "05_immuneSig_ICBresponse/validation_results/oriImmuSig_mlr_model/alldatasets_unscaled/"
dir.create(resultpath)

if(featureFilterStrategy == "mlrFilter"){
  if(mlrFilter %in% c('FSelector_gain.ratio',
                      'FSelector_chi.squared',
                      'FSelector_information.gain',
                      'FSelector_symmetrical.uncertainty')){
    dir.create(paste0(resultpath, "_", mlrFilter ,"/"))
    resultpath = paste0(resultpath, "_",mlrFilter ,"/")
    
  }else if(mlrFilter %in%  c("auc","variance","kruskal.test", # "univariate.model.score",
                             "FSelector_relief",
                             "praznik_CMIM","praznik_DISR","praznik_JMI","praznik_JMIM",
                             "praznik_MIM","praznik_MRMR","praznik_NJMIM",
                             # "party_cforest.importance",
                             "randomForest_importance", 
                             "randomForestSRC_importance", 
                             "randomForestSRC_var.select",
                             "ranger_impurity", "ranger_permutation")){
    dir.create(paste0(resultpath, "_",mlrFilter,"_", mlrFilter_k,"/"))
    resultpath = paste0(resultpath, "_",mlrFilter,"_", mlrFilter_k,"/")
    
  }
}else if(featureFilterStrategy == "mixed"){
  if(mlrFilter %in% c('FSelector_gain.ratio',
                      'FSelector_chi.squared',
                      'FSelector_information.gain',
                      'FSelector_symmetrical.uncertainty')){
    dir.create(paste0(resultpath, "_",featureFilterStrategy, "_", mlrFilter ,"/"))
    resultpath = paste0(resultpath, "_",featureFilterStrategy, "_",mlrFilter ,"/")
    
  }else if(mlrFilter %in%  c("auc","variance","kruskal.test", # "univariate.model.score",
                             "FSelector_relief",
                             "praznik_CMIM","praznik_DISR","praznik_JMI","praznik_JMIM",
                             "praznik_MIM","praznik_MRMR","praznik_NJMIM",
                             # "party_cforest.importance",
                             "randomForest_importance", 
                             "randomForestSRC_importance", 
                             "randomForestSRC_var.select",
                             "ranger_impurity", "ranger_permutation")){
    dir.create(paste0(resultpath, "_",featureFilterStrategy, "_",mlrFilter,"_", mlrFilter_k,"/"))
    resultpath = paste0(resultpath, "_",featureFilterStrategy, "_",mlrFilter,"_", mlrFilter_k,"/")

  }
}else if(featureFilterStrategy == "meta_mixed"){
  if(mlrFilter %in% c('FSelector_gain.ratio',
                      'FSelector_chi.squared',
                      'FSelector_information.gain',
                      'FSelector_symmetrical.uncertainty')){
    dir.create(paste0(resultpath, "_",featureFilterStrategy, "_", mlrFilter ,"/"))
    resultpath = paste0(resultpath, "_",featureFilterStrategy, "_",mlrFilter ,"/")

  }else if(mlrFilter %in%  c("auc","variance","kruskal.test", # "univariate.model.score",
                             "FSelector_relief",
                             "praznik_CMIM","praznik_DISR","praznik_JMI","praznik_JMIM",
                             "praznik_MIM","praznik_MRMR","praznik_NJMIM",
                             # "party_cforest.importance",
                             "randomForest_importance", 
                             "randomForestSRC_importance", 
                             "randomForestSRC_var.select",
                             "ranger_impurity", "ranger_permutation")){
    dir.create(paste0(resultpath, "_",featureFilterStrategy, "_",mlrFilter,"_", mlrFilter_k,"/"))
    resultpath = paste0(resultpath, "_",featureFilterStrategy, "_",mlrFilter,"_", mlrFilter_k,"/")

  }
}else{
  
  dir.create(paste0(resultpath, "/"))
  resultpath = paste0(resultpath, "/", featureFilterStrategy, "/")
  dir.create(paste0(resultpath, "/"))
  
}


# ------------------------------------------------------------------------------
# create savepath
# ------------------------------------------------------------------------------

savepath = paste0(resultpath, "/", dataset, "/")
dir.create(savepath)

# ------------------------------------------------------------------------------
# prediction on validation data ---------------------------------------------------------------------
# ------------------------------------------------------------------------------

# dataset = "Liu_dataset_CTLA4Naive_outcome"
cancer = "SKCM"
sample_type = "Metastatic"

result_dir = paste0("05_immuneSig_ICBresponse/validation_results/", dataset,"/")
data_test <- read.csv(paste0(result_dir, "/", cancer, "_", sample_type, "_immuneSig_merged_mlr.csv"))
data_test_name <- data.frame(name = names(data_test), immune_sig = tolower(names(data_test)))
data_test_name <- immunesig_name_unify(data_test_name)
names(data_test) <- data_test_name$immune_sig
data_test$tag <- 1
data_test[data_test$label == "NR",]$tag = 0
data_test$tag <- factor(data_test$tag,levels = c(1,0), order=TRUE)
data_test[is.na(data_test)] = 0

# ------------------------------------------------------------------------------
load(paste0(resultpath, "/rf/ranger_model_pred.Rdata"))
if(length(intersect(c(sort(m_rf$features), 'tag'), names(data_test))) == length(c(sort(m_rf$features), 'tag'))){
     data_test_rf <- data_test[c(sort(m_rf$features), 'tag')]
}else{
     data_test_rf <- data_test[intersect(c(sort(m_rf$features), 'tag'), names(data_test))]
     missing_feature <- setdiff(c(sort(m_rf$features), 'tag'), names(data_test))
     tmp = as.data.frame(matrix(0,nrow(data_test), length(missing_feature)))
     names(tmp) <- missing_feature
     data_test_rf <- cbind(data_test_rf, tmp)
}

ICBresponse_testtask = makeClassifTask(id = "ICBresponse", 
                                       data = data_test_rf, 
                                       target = "tag", 
                                       positive = 1)



pred_rf = predict(m_rf, task = ICBresponse_testtask)
calculateConfusionMatrix(pred_rf)
pred2 <- ROCR::prediction(pred_rf$data$prob.0 , pred_rf$data$truth)
aucPerf <- ROCR::performance( pred2, 'auc' )
AUCValue<-aucPerf@y.values[[1]]
AUCValue


modelroc <- pROC::roc( pred_rf$data$truth, pred_rf$data$prob.0 )
pdf(paste0(savepath, "/rf_ROC.pdf"))
plot(modelroc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),
     grid.col=c("green", "red"), max.auc.polygon=TRUE,
     auc.polygon.col="lightblue", print.thres=TRUE)
dev.off()

save(ranger_opt_ps_res, mod_ranger_opt, m_rf, pred_rf, 
     file = paste0(savepath, "/ranger_model_pred.Rdata"))


# ------------------------------------------------------------------------------
# elastic net from glmnet ------------------------------------------------------
# ------------------------------------------------------------------------------
load(paste0(resultpath, "/en/elasticnet_model_pred.Rdata"))
if(length(intersect(c(sort(m_en$features), 'tag'), names(data_test))) == length(c(sort(m_en$features), 'tag'))){
     data_test_en <- data_test[c(sort(m_en$features), 'tag')]
}else{
     data_test_en <- data_test[intersect(c(sort(m_en$features), 'tag'), names(data_test))]
     missing_feature <- setdiff(c(sort(m_en$features), 'tag'), names(data_test))
     tmp = as.data.frame(matrix(0,nrow(data_test), length(missing_feature)))
     names(tmp) <- missing_feature
     data_test_en <- cbind(data_test_en, tmp)
}

ICBresponse_testtask = makeClassifTask(id = "ICBresponse", 
                                       data = data_test_en, 
                                       target = "tag", 
                                       positive = 1)

pred_en = predict(m_en, task = ICBresponse_testtask)
calculateConfusionMatrix(pred_en)
#  mlr::performance(pred_en, mlr::auc)


pred2 <- ROCR::prediction(pred_en$data$prob.0 , pred_en$data$truth)
aucPerf <- ROCR::performance( pred2, 'auc' )
AUCValue<-aucPerf@y.values[[1]]
AUCValue


modelroc <- pROC::roc( pred_en$data$truth, pred_en$data$prob.0 )
pdf(paste0(savepath, "/en_ROC.pdf"))
plot(modelroc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),
     grid.col=c("green", "red"), max.auc.polygon=TRUE,
     auc.polygon.col="lightblue", print.thres=TRUE)
dev.off()

save(elasticnet_opt_ps_res, mod_elasticnet_opt, m_en, pred_en, 
     file = paste0(savepath, "/elasticnet_model_pred.Rdata"))

# ------------------------------------------------------------------------------
# ridge from glmnet ------------------------------------------------------
# ------------------------------------------------------------------------------
load(paste0(resultpath, "/ridge/ridge_model_pred.Rdata"))
if(length(intersect(c(sort(m_ridge$features), 'tag'), names(data_test))) == length(c(sort(m_ridge$features), 'tag'))){
     data_test_ridge <- data_test[c(sort(m_ridge$features), 'tag')]
}else{
     data_test_ridge <- data_test[intersect(c(sort(m_ridge$features), 'tag'), names(data_test))]
     missing_feature <- setdiff(c(sort(m_ridge$features), 'tag'), names(data_test))
     tmp = as.data.frame(matrix(0,nrow(data_test), length(missing_feature)))
     names(tmp) <- missing_feature
     data_test_ridge <- cbind(data_test_ridge, tmp)
}

ICBresponse_testtask = makeClassifTask(id = "ICBresponse", 
                                       data = data_test_ridge, 
                                       target = "tag", 
                                       positive = 0)

pred_ridge = predict(m_ridge, task = ICBresponse_testtask)
calculateConfusionMatrix(pred_ridge)
# mlr::performance(pred_ridge, mlr::auc)

# pred_ridge$data
pred2 <- ROCR::prediction(pred_ridge$data$prob.0 , pred_ridge$data$truth)
aucPerf <- ROCR::performance( pred2, 'auc' )
AUCValue<-aucPerf@y.values[[1]]
AUCValue


modelroc <- pROC::roc( pred_ridge$data$truth, pred_ridge$data$prob.0 )
pdf(paste0(savepath, "/ridge_ROC.pdf"))
plot(modelroc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),
     grid.col=c("green", "red"), max.auc.polygon=TRUE,
     auc.polygon.col="lightblue", print.thres=TRUE)
dev.off()

save(ridge_opt_ps_res, mod_ridge_opt, m_ridge, pred_ridge, 
     file = paste0(savepath, "/ridge_model_pred.Rdata"))


# ---------------------------------------------------------------------
# ksvm  ------------------------------------------------------
# ---------------------------------------------------------------------
load(paste0(resultpath, "/svm/svm_model_pred.Rdata"))
if(length(intersect(c(sort(m_svm$features), 'tag'), names(data_test))) == length(c(sort(m_svm$features), 'tag'))){
     data_test_svm <- data_test[c(sort(m_svm$features), 'tag')]
}else{
     data_test_svm <- data_test[intersect(c(sort(m_svm$features), 'tag'), names(data_test))]
     missing_feature <- setdiff(c(sort(m_svm$features), 'tag'), names(data_test))
     tmp = as.data.frame(matrix(0,nrow(data_test), length(missing_feature)))
     names(tmp) <- missing_feature
     data_test_svm <- cbind(data_test_svm, tmp)
}
ICBresponse_testtask = makeClassifTask(id = "ICBresponse", 
                                       data = data_test_svm, 
                                       target = "tag", 
                                       positive = 1)

pred_svm = predict(m_svm, task = ICBresponse_testtask)
calculateConfusionMatrix(pred_svm)
pred2 <- ROCR::prediction(pred_svm$data$prob.0 , pred_svm$data$truth)
aucPerf <- ROCR::performance( pred2, 'auc' )
AUCValue<-aucPerf@y.values[[1]]
AUCValue

modelroc <- pROC::roc( pred_svm$data$truth, pred_svm$data$prob.0 )
pdf(paste0(savepath, "/svm_ROC.pdf"))
plot(modelroc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),
     grid.col=c("green", "red"), max.auc.polygon=TRUE,
     auc.polygon.col="lightblue", print.thres=TRUE)
dev.off()

save(ksvm_opt_ps_res, mod_ksvm_opt, m_svm, pred_svm, 
     file = paste0(savepath, "/svm_model_pred.Rdata"))


# xgboost ----------------------------------------------------------------------

load(paste0(resultpath, "/xgboost/xgboost_model_pred.Rdata"))
if(length(intersect(c(sort(m_xgboost$features), 'tag'), names(data_test))) == length(c(sort(m_xgboost$features), 'tag'))){
     data_test_xgboost <- data_test[c(sort(m_xgboost$features), 'tag')]
}else{
     data_test_xgboost <- data_test[intersect(c(sort(m_xgboost$features), 'tag'), names(data_test))]
     missing_feature <- setdiff(c(sort(m_xgboost$features), 'tag'), names(data_test))
     tmp = as.data.frame(matrix(0,nrow(data_test), length(missing_feature)))
     names(tmp) <- missing_feature
     data_test_xgboost <- cbind(data_test_xgboost, tmp)
}
ICBresponse_testtask = makeClassifTask(id = "ICBresponse", 
                                       data = data_test_xgboost, 
                                       target = "tag", 
                                       positive = 1)

pred_xgboost = predict(m_xgboost, task = ICBresponse_testtask)
calculateConfusionMatrix(pred_xgboost)
# df = generateThreshVsPerfData(pred_xgboost, measures = list(fpr, tpr))
# plotROCCurves(df)
#     mlr::performance(pred_xgboost, mlr::auc)
pred_xgboost$data
pred2 <- ROCR::prediction(pred_xgboost$data$prob.0 , pred_xgboost$data$truth)
aucPerf <- ROCR::performance( pred2, 'auc' )
AUCValue<-aucPerf@y.values[[1]]
AUCValue

modelroc <- pROC::roc( pred_xgboost$data$truth, pred_xgboost$data$prob.0 )
pdf(paste0(savepath, "/xgboost_ROC.pdf"))
plot(modelroc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),
     grid.col=c("green", "red"), max.auc.polygon=TRUE,
     auc.polygon.col="lightblue", print.thres=TRUE)
dev.off()

save(xgboost_opt_ps_res, mod_xgboost_opt, m_xgboost, pred_xgboost, 
     file = paste0(savepath, "/xgboost_model_pred.Rdata"))

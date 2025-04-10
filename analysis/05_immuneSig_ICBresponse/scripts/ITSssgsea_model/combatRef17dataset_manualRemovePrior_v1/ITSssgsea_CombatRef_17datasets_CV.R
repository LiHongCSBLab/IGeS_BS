# rm(list=ls())
workpath <- "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/"
setwd(workpath)

options(stringsAsFactors = F)

# Checking and loading R packages -----------------------------------------

if (!require(readxl)) {  install.packages('readxl')}
if (!require(data.table)) {  install.packages('data.table')}
if (!require(dplyr)) {  install.packages('dplyr')}
library(ROCR)
library(reshape2)
library(getopt)


# command=matrix(c( 
#   "featureFilterStrategy",  "s", 1, "character", "featureFilterStrategy in metaAnalysis_removeprior_v2, combinedP, voting, mlrFilter",
#   "cancer",                 "c", 2, "character", "cancer",
#   "mlrFilter",              "f", 2, "character", "mlrFilter in mlr",
#   "mlrFilter_k",            "k", 2, "numeric", "top k persent feature, k from 0-1",
#   "help",                   "h", 0, "logical",   "help file"),
#   byrow=T,ncol=5)

# args=getopt(command)

# featureFilterStrategy = args$featureFilterStrategy
# cancer = args$cancer
# mlrFilter = args$mlrFilter
# mlrFilter_k = args$mlrFilter_k

# cancer = "SKCM"
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

# ------------------------------------------------------------------------------
# feature selection ------———---------------------------------------------------
# and create directory according to strategy and model -------------------------
# ------------------------------------------------------------------------------
result_dir <- "05_immuneSig_ICBresponse/results_r0.4_ITSssgsea_CombatRef/"


resultpath = paste0(result_dir, "modelCV_17datasets_ssgseanorm/")
dir.create(resultpath)

featureFilterStrategy = "metaAnalysis_removeprior_v1" 

  dir.create(paste0(resultpath, "/"))
  resultpath = paste0(resultpath, "/", featureFilterStrategy, "/")
  dir.create(paste0(resultpath, "/"))
  dir.create(paste0(resultpath, "/rf"))
  dir.create(paste0(resultpath, "/en"))
  dir.create(paste0(resultpath, "/ridge"))
  dir.create(paste0(resultpath, "/svm"))
  dir.create(paste0(resultpath, "/xgboost"))
  
  features_selected <- feature_filter(data,
                                      strategy = "metaAnalysis",
                                      p_value = 0.05,
                                      cancer = NULL, 
                                      mlrFilter = NULL,
                                      mlrFilter_k = NULL,
                                      workpath = workpath) 
  

print(length(features_selected))
# print((features_selected))

sig_filter_prior <- c("exhausted_immucellai", 
                      "merck18_tide",
                      "oe_tme.t.cd4.treg",
                      "mdsc_tcga_caf_immunecls",
                      "oe_co.culture.screen.hits10",
                      "oe_in.vivo.screen.gvax.vs.tcrako.enriched"
                      )

features_selected <- setdiff(features_selected, sig_filter_prior)


# load data --------------------------------------------------------------------------

meta_result <- read.csv("05_immuneSig_ICBresponse/results/meta_results/res_metaRes.csv")
meta_selected <- meta_result[meta_result$Meta_Pval < 0.05, ]
meta_selected$immune_sig = meta_selected$X

result_dir <- "05_immuneSig_ICBresponse/results_r0.4_ITSssgsea_CombatRef/"

load( paste0(result_dir, "ITS_afterCombat17datasets_ssgseanorm.Rdata"))
datasets = unique(ITSres_selected$dataset)


# pdf(paste0(result_dir, "ITS_afterCombat_ssgseanorm_heatmap.pdf"), width = 12, height = 8)
# for(i in 1:length(datasets)){
# ds <- ITSres_selected[ITSres_selected$dataset == datasets[i],]
# ds <- ds[order(ds$label),]
# ds_Annot <- ds[c('patientID','label','SampleID',
#                      'cancer','sample_type','dataset')]
# ds_selected2 <- ds[-which(is.element(names(ds), names(ds_Annot)))]
# pheatmap::pheatmap(ds_selected2, scale = 'column',
#                    cluster_row = F, cluster_col = F,
#                    annotation_row = ds_Annot['label'])
# }
# dev.off()

# data seperation --------------------------------------------------------------------------

vlSetName <- c("08_Lauss_dataset_outcome", "GSE115821_pre_outcome", "GSE96619_outcome","Liu_dataset_CTLA4Naive_outcome", "GSE176307")
dcSetName <- datasets[!is.element(datasets, vlSetName)]

dcSet <- ITSres_selected[is.element(ITSres_selected$dataset, dcSetName),]
vlSet <- ITSres_selected[is.element(ITSres_selected$dataset, vlSetName),]

dcSet_Annot <- dcSet[c('patientID','SampleID', 'cancer','sample_type','dataset')]
dcSet_ITS <- dcSet[-which(is.element(names(dcSet), names(dcSet_Annot)))]

# ------------------------------------------------------------------------------
# seperate dataset
# ------------------------------------------------------------------------------
set.seed(1234)
library(sampling)
# library(caret)  
res_df <- dcSet_ITS[!is.na(dcSet_ITS$label), ]

res_df = res_df[order(res_df$label), ]
res_df$label <- factor(res_df$label,levels = c("NR","R"), order=TRUE)
res_df$tag <- 1
res_df[res_df$label == "NR",]$tag = 0
res_df$tag <- factor(res_df$tag,levels = c(0, 1), order=TRUE)

folds<-caret::createFolds(y=res_df$label,k=5) 

modeling_task <- lapply(as.list(names(folds)), function(i){
  
  # i = 1
  ID_unit=folds[[i]]
  # sub_train = strata(res_df, strataname = ("label"), size= n, method = "srswor")
  
  data_train = res_df[-ID_unit, ]
  
  #   rownames(data_train) <- paste0(data_train$dataset,"_",data_train$patient)
  data_test = res_df[ID_unit, ]
  #   rownames(data_test) <- paste0(data_test$dataset,"_",data_test$patient)
  
  # ------------------------------------------------------------------------------
  # create task
  # ------------------------------------------------------------------------------
  
  data_train_filter <- data_train[c(intersect(features_selected, names(data_train)), 'tag')]
  
  filtered.task = makeClassifTask(id = "ICBresponse", 
                                  data = data_train_filter, 
                                  target = "tag", 
                                  positive = 1)
  
  print(filtered.task)
  # ------------------------------------------------------------------------------
  # feature selection for test set 
  # ------------------------------------------------------------------------------
  
  data_test <- data_test[c(getTaskFeatureNames(filtered.task), 'tag')]
  ICBresponse_testtask = makeClassifTask(id = "ICBresponse", 
                                         data = data_test, 
                                         target = "tag", 
                                         positive = 1)
  print(ICBresponse_testtask)
  
  # ------------------------------------------------------------------------------
  # tuning parameter
  # ------------------------------------------------------------------------------
  print("start modeling")
  
  # ------------------------------------------------------------------------------
  # ridge from glmnet ------------------------------------------------------
  # ------------------------------------------------------------------------------
  rdesc =  makeResampleDesc("CV", stratify = T,  iters = 10L)
  
  ridge_opt_ps_res = ridge_tuning(traindata = filtered.task)
  ridge_opt_ps = ridge_opt_ps_res$x
  # set up model with optimal parameters
  mod_ridge_opt = setHyperPars(learner = makeLearner("classif.cvglmnet", 
                                                     predict.type = "prob", 
                                                     fix.factors.prediction = TRUE,
                                                     nlambda = 1000L,
                                                     lambda.min.ratio = 1e-5,
                                                     nfolds = 5,
                                                     config = list(on.learner.error = "warn")) ,
                               standardize =  ridge_opt_ps$standardize,
                               s =  ridge_opt_ps$s,
                               alpha =  ridge_opt_ps$alpha
  )
  
  
  m_ridge = train(mod_ridge_opt, filtered.task)
  pred_ridge = predict(m_ridge, task = ICBresponse_testtask)
  calculateConfusionMatrix(pred_ridge)
  # mlr::performance(pred_ridge, mlr::auc)
  
  pred_ridge$data
  pred2 <- ROCR::prediction(pred_ridge$data$prob.1 , pred_ridge$data$truth)
  aucPerf <- ROCR::performance( pred2, 'auc' )
  AUCValue<-aucPerf@y.values[[1]]
  AUCValue
  
  
  modelroc <- pROC::roc( pred_ridge$data$truth, pred_ridge$data$prob.1 )
  pdf(paste0(resultpath, "/ridge/ROC_",i,".pdf"))
  plot(modelroc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),
       grid.col=c("green", "red"), max.auc.polygon=TRUE,
       auc.polygon.col="lightblue", print.thres=TRUE)
  dev.off()
  
  save(ridge_opt_ps_res, mod_ridge_opt, m_ridge, pred_ridge, 
       file = paste0(resultpath, "/ridge/ridge_model_pred_",i,".Rdata"))
})


# ------------------------------------------------------------------------------
# elastic net from glmnet ------------------------------------------------------
# ------------------------------------------------------------------------------
modeling_task <- lapply(as.list(names(folds)), function(i){
  # i = 1
  ID_unit=folds[[i]]
  # sub_train = strata(res_df, strataname = ("label"), size= n, method = "srswor")
  
  data_train = res_df[-ID_unit, ]
  #   rownames(data_train) <- paste0(data_train$dataset,"_",data_train$patient)
  data_test = res_df[ID_unit, ]
  #   rownames(data_test) <- paste0(data_test$dataset,"_",data_test$patient)
  
  # ------------------------------------------------------------------------------
  # create task
  # ------------------------------------------------------------------------------
  
  data_train_filter <- data_train[c(intersect(features_selected, names(data_train)), 'tag')]
  
  filtered.task = makeClassifTask(id = "ICBresponse", 
                                  data = data_train_filter, 
                                  target = "tag", 
                                  positive = 1)
  # ------------------------------------------------------------------------------
  # feature selection for test set 
  # ------------------------------------------------------------------------------
  
  data_test <- data_test[c(getTaskFeatureNames(filtered.task), 'tag')]
  ICBresponse_testtask = makeClassifTask(id = "ICBresponse", 
                                         data = data_test, 
                                         target = "tag", 
                                         positive = 1)
  
  # ------------------------------------------------------------------------------
  # tuning parameter
  # ------------------------------------------------------------------------------
  print("start modeling")
  
  rdesc =  makeResampleDesc("CV", stratify = T,  iters = 10L)
  
  elasticnet_opt_ps_res = elasticnet_tuning(traindata = filtered.task)
  elasticnet_opt_ps = elasticnet_opt_ps_res$x
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
  
  
  m_en = train(mod_elasticnet_opt, filtered.task)
  pred_en = predict(m_en, task = ICBresponse_testtask)
  calculateConfusionMatrix(pred_en)
  #  mlr::performance(pred_en, mlr::auc)
  
  pred_en$data
  pred2 <- ROCR::prediction(pred_en$data$prob.1 , pred_en$data$truth)
  aucPerf <- ROCR::performance( pred2, 'auc' )
  AUCValue<-aucPerf@y.values[[1]]
  AUCValue
  
  
  modelroc <- pROC::roc( pred_en$data$truth, pred_en$data$prob.1 )
  pdf(paste0(resultpath, "/en/ROC_",i,".pdf"))
  plot(modelroc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),
       grid.col=c("green", "red"), max.auc.polygon=TRUE,
       auc.polygon.col="lightblue", print.thres=TRUE)
  dev.off()
  
  save(elasticnet_opt_ps_res, mod_elasticnet_opt, m_en, pred_en, 
       file = paste0(resultpath, "/en/elasticnet_model_pred_",i,".Rdata"))
})



modeling_task <- lapply(as.list(names(folds)), function(i){
  # i = 1
  ID_unit=folds[[i]]
  # sub_train = strata(res_df, strataname = ("label"), size= n, method = "srswor")
  
  data_train = res_df[-ID_unit, ]
  #   rownames(data_train) <- paste0(data_train$dataset,"_",data_train$patient)
  data_test = res_df[ID_unit, ]
  #   rownames(data_test) <- paste0(data_test$dataset,"_",data_test$patient)
  
  # ------------------------------------------------------------------------------
  # create task
  # ------------------------------------------------------------------------------
  
  data_train_filter <- data_train[c(intersect(features_selected, names(data_train)), 'tag')]
  
  filtered.task = makeClassifTask(id = "ICBresponse", 
                                  data = data_train_filter, 
                                  target = "tag", 
                                  positive = 1)
  # ------------------------------------------------------------------------------
  # feature selection for test set 
  # ------------------------------------------------------------------------------
  
  data_test <- data_test[c(getTaskFeatureNames(filtered.task), 'tag')]
  ICBresponse_testtask = makeClassifTask(id = "ICBresponse", 
                                         data = data_test, 
                                         target = "tag", 
                                         positive = 1)
  
  # ------------------------------------------------------------------------------
  # tuning parameter
  # ------------------------------------------------------------------------------
  print("start modeling")
  
  # ranger: Random Forest from ranger
  # tuning parameters, return optimal hyperparameters
  ranger_opt_ps_res = rf_tuning(traindata = filtered.task)
  ranger_opt_ps = ranger_opt_ps_res$x
  # set up model with optimal parameters
  mod_ranger_opt = setHyperPars(learner = makeLearner("classif.ranger", 
                                                      predict.type = "prob", 
                                                      fix.factors.prediction = TRUE) ,
                                num.trees =ranger_opt_ps$num.trees, 
                                mtry = ranger_opt_ps$mtry,
                                min.node.size = ranger_opt_ps$min.node.size,
                                replace = ranger_opt_ps$replace)
  
  rdesc =  makeResampleDesc("CV", stratify = T,  iters = 10L)
  
  
  m_rf = train(mod_ranger_opt, filtered.task)
  pred_rf = predict(m_rf, task = ICBresponse_testtask)
  calculateConfusionMatrix(pred_rf)
  # df = generateThreshVsPerfData(pred, measures = list(fpr, tpr))
  # plotROCCurves(df)
  #  mlr::performance(pred_rf, mlr::auc)
  # pred$data
  pred2 <- ROCR::prediction(pred_rf$data$prob.1 , pred_rf$data$truth)
  aucPerf <- ROCR::performance( pred2, 'auc' )
  AUCValue<-aucPerf@y.values[[1]]
  AUCValue
  
  
  modelroc <- pROC::roc( pred_rf$data$truth, pred_rf$data$prob.1 )
  pdf(paste0(resultpath, "/rf/ROC_",i,".pdf"))
  plot(modelroc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),
       grid.col=c("green", "red"), max.auc.polygon=TRUE,
       auc.polygon.col="lightblue", print.thres=TRUE)
  dev.off()
  
  save(ranger_opt_ps_res, mod_ranger_opt, m_rf, pred_rf, 
       file = paste0(resultpath, "/rf/ranger_model_pred_",i,".Rdata"))
  
})




# ksvm  ------------------------------------------------------
modeling_task <- lapply(as.list(names(folds)), function(i){
  # i = 1
  ID_unit=folds[[i]]
  # sub_train = strata(res_df, strataname = ("label"), size= n, method = "srswor")
  
  data_train = res_df[-ID_unit, ]
  #   rownames(data_train) <- paste0(data_train$dataset,"_",data_train$patient)
  data_test = res_df[ID_unit, ]
  #   rownames(data_test) <- paste0(data_test$dataset,"_",data_test$patient)
  
  # ------------------------------------------------------------------------------
  # create task
  # ------------------------------------------------------------------------------
  
  data_train_filter <- data_train[c(intersect(features_selected, names(data_train)), 'tag')]
  
  filtered.task = makeClassifTask(id = "ICBresponse", 
                                  data = data_train_filter, 
                                  target = "tag", 
                                  positive = 1)
  # ------------------------------------------------------------------------------
  # feature selection for test set 
  # ------------------------------------------------------------------------------
  
  data_test <- data_test[c(getTaskFeatureNames(filtered.task), 'tag')]
  ICBresponse_testtask = makeClassifTask(id = "ICBresponse", 
                                         data = data_test, 
                                         target = "tag", 
                                         positive = 1)
    # new data ---------------------------------------------------------------------
  
  
  # rdesc = makeResampleDesc("RepCV", folds = 10, reps = 10)
  rdesc =  makeResampleDesc("CV", stratify = T,  iters = 10L)
  
  # ksvm: Support Vector Machines from kernlab
  # tuning parameters, return optimal hyperparameters
  ksvm_opt_ps_res = csvm_tuning(traindata = filtered.task)
  ksvm_opt_ps = ksvm_opt_ps_res$x
  
  if(ksvm_opt_ps$kernel == 'rbfdot'){
    # set up model with optimal parameters
    mod_ksvm_opt = setHyperPars(learner = makeLearner("classif.ksvm", 
                                                      predict.type = "prob", 
                                                      fix.factors.prediction = TRUE) ,
                                # type=ksvm_opt_ps$type,
                                C=ksvm_opt_ps$C,
                                kernel = ksvm_opt_ps$kernel,
                                sigma = ksvm_opt_ps$sigma)
    
  }else if(ksvm_opt_ps$kernel == 'vanilladot'){
    
    mod_ksvm_opt = setHyperPars(learner = makeLearner("classif.ksvm", 
                                                      predict.type = "prob", 
                                                      fix.factors.prediction = TRUE) ,
                                # type=ksvm_opt_ps$type,
                                C=ksvm_opt_ps$C,
                                kernel = ksvm_opt_ps$kernel
    )
    
  }else if(ksvm_opt_ps$kernel == 'polydot'){
    mod_ksvm_opt = setHyperPars(learner = makeLearner("classif.ksvm", 
                                                      predict.type = "prob", 
                                                      fix.factors.prediction = TRUE) ,
                                # type=ksvm_opt_ps$type,
                                C=ksvm_opt_ps$C,
                                kernel = ksvm_opt_ps$kernel,
                                degree = ksvm_opt_ps$degree,
                                scale = ksvm_opt_ps$scale,
                                offset = ksvm_opt_ps$offset,
                                sigma = ksvm_opt_ps$sigma)
    
  }else if(ksvm_opt_ps$kernel == 'laplacedot'){
    mod_ksvm_opt = setHyperPars(learner = makeLearner("classif.ksvm", 
                                                      predict.type = "prob", 
                                                      fix.factors.prediction = TRUE) ,
                                # type=ksvm_opt_ps$type,
                                C=ksvm_opt_ps$C,
                                kernel = ksvm_opt_ps$kernel,
                                sigma = ksvm_opt_ps$sigma)
    
  }else if(ksvm_opt_ps$kernel == 'besseldot'){
    mod_ksvm_opt = setHyperPars(learner = makeLearner("classif.ksvm", 
                                                      predict.type = "prob", 
                                                      fix.factors.prediction = TRUE) ,
                                # type=ksvm_opt_ps$type,
                                C=ksvm_opt_ps$C,
                                kernel = ksvm_opt_ps$kernel,
                                sigma = ksvm_opt_ps$sigma,
                                order = ksvm_opt_ps$order)
  }
  
  
  
  m_svm = train(mod_ksvm_opt, filtered.task)
  pred_svm = predict(m_svm, task = ICBresponse_testtask)
  calculateConfusionMatrix(pred_svm)
  # df = generateThreshVsPerfData(pred_svm, measures = list(fpr, tpr))
  # plotROCCurves(df)
  # mlr::performance(pred_svm, mlr::auc)
  pred_svm$data
  pred2 <- ROCR::prediction(pred_svm$data$prob.1 , pred_svm$data$truth)
  aucPerf <- ROCR::performance( pred2, 'auc' )
  AUCValue<-aucPerf@y.values[[1]]
  AUCValue
  
  modelroc <- pROC::roc( pred_svm$data$truth, pred_svm$data$prob.1 )
  pdf(paste0(resultpath, "/svm/ROC_",i,".pdf"))
  plot(modelroc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),
       grid.col=c("green", "red"), max.auc.polygon=TRUE,
       auc.polygon.col="lightblue", print.thres=TRUE)
  dev.off()
  
  save(ksvm_opt_ps_res, mod_ksvm_opt, m_svm, pred_svm, 
       file = paste0(resultpath, "/svm/svm_model_pred_",i,".Rdata"))
})






# xgboost ----------------------------------------------------------------------
modeling_task <- lapply(as.list(names(folds)), function(i){
  # i = 1
  ID_unit=folds[[i]]
  # sub_train = strata(res_df, strataname = ("label"), size= n, method = "srswor")
  
  data_train = res_df[-ID_unit, ]
  #   rownames(data_train) <- paste0(data_train$dataset,"_",data_train$patient)
  data_test = res_df[ID_unit, ]
  #   rownames(data_test) <- paste0(data_test$dataset,"_",data_test$patient)
  
  # ------------------------------------------------------------------------------
  # create task
  # ------------------------------------------------------------------------------
  
  data_train_filter <- data_train[c(intersect(features_selected, names(data_train)), 'tag')]
  
  filtered.task = makeClassifTask(id = "ICBresponse", 
                                  data = data_train_filter, 
                                  target = "tag", 
                                  positive = 1)
  # ------------------------------------------------------------------------------
  # feature selection for test set 
  # ------------------------------------------------------------------------------
  
  data_test <- data_test[c(getTaskFeatureNames(filtered.task), 'tag')]
  ICBresponse_testtask = makeClassifTask(id = "ICBresponse", 
                                         data = data_test, 
                                         target = "tag", 
                                         positive = 1)
    # new data ---------------------------------------------------------------------
  
  xgboost_opt_ps_res = xgboost_tuning(filtered.task)
  xgboost_opt_ps = xgboost_opt_ps_res$x
  if( xgboost_opt_ps$booster == 'gbtree'){
    # set up model with optimal parameters
    mod_xgboost_opt = setHyperPars(learner = makeLearner("classif.xgboost", 
                                                         predict.type = "prob", 
                                                         fix.factors.prediction = TRUE,
                                                         print_every_n = 1000L,
                                                         nthread = 1L,
                                                         eval_metric='logloss',                                                       
                                                         config = list(on.learner.error = "warn")) ,
                                   # optimal parameters
                                   booster = xgboost_opt_ps$booster,
                                   eta = xgboost_opt_ps$eta,
                                   gamma = xgboost_opt_ps$gamma,
                                   max_depth = xgboost_opt_ps$max_depth,
                                   min_child_weight = xgboost_opt_ps$min_child_weight,
                                   subsample = xgboost_opt_ps$subsample,
                                   colsample_bytree = xgboost_opt_ps$colsample_bytree,
                                   colsample_bylevel = xgboost_opt_ps$colsample_bylevel,
                                   base_score = xgboost_opt_ps$base_score,
                                   nrounds = xgboost_opt_ps$nrounds
    )
  }else if(xgboost_opt_ps$booster == 'gblinear'){
    
    mod_xgboost_opt = setHyperPars(learner = makeLearner("classif.xgboost", 
                                                         predict.type = "prob", 
                                                         fix.factors.prediction = TRUE,
                                                         print_every_n = 1000L,
                                                         nthread = 1L,
                                                         config = list(on.learner.error = "warn")) ,
                                   # optimal parameters
                                   booster = xgboost_opt_ps$booster,
                                   lambda = xgboost_opt_ps$lambda,
                                   lambda_bias = xgboost_opt_ps$lambda_bias,
                                   alpha = xgboost_opt_ps$alpha,
                                   base_score = xgboost_opt_ps$base_score,
                                   nrounds = xgboost_opt_ps$nrounds
    )
    
  }
  
  
  m_xgboost = train(mod_xgboost_opt, filtered.task)
  pred_xgboost = predict(m_xgboost, task = ICBresponse_testtask)
  calculateConfusionMatrix(pred_xgboost)
  # df = generateThreshVsPerfData(pred_xgboost, measures = list(fpr, tpr))
  # plotROCCurves(df)
  #     mlr::performance(pred_xgboost, mlr::auc)
  pred_xgboost$data
  pred2 <- ROCR::prediction(pred_xgboost$data$prob.1 , pred_xgboost$data$truth)
  aucPerf <- ROCR::performance( pred2, 'auc' )
  AUCValue<-aucPerf@y.values[[1]]
  AUCValue
  
  modelroc <- pROC::roc( pred_xgboost$data$truth, pred_xgboost$data$prob.1 )
  pdf(paste0(resultpath, "/xgboost/ROC_",i,".pdf"))
  plot(modelroc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),
       grid.col=c("green", "red"), max.auc.polygon=TRUE,
       auc.polygon.col="lightblue", print.thres=TRUE)
  dev.off()
  
  save(xgboost_opt_ps_res, mod_xgboost_opt, m_xgboost, pred_xgboost, 
       file = paste0(resultpath, "/xgboost/xgboost_model_pred_",i,".Rdata"))
  
})

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


command=matrix(c( 
  "featureFilterStrategy",  "s", 1, "character", "featureFilterStrategy in metaAnalysis, combinedP, voting, mlrFilter",
  "cancer",                 "c", 2, "character", "cancer",
  "mlrFilter",              "f", 2, "character", "mlrFilter in mlr",
  "mlrFilter_k",            "k", 2, "numeric", "top k persent feature, k from 0-1",
  "dataset",                "d", 1, "character", "validation dataset",
  "help",                   "h", 0, "logical",   "help file"),
  byrow=T,ncol=5)

args=getopt(command)

featureFilterStrategy = args$featureFilterStrategy
cancer = args$cancer
mlrFilter = args$mlrFilter
mlrFilter_k = args$mlrFilter_k
dataset = args$dataset

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
load("05_immuneSig_ICBresponse/results/oriImmuneSig_res_all.Rdata")
merged_sig1 <- lapply(as.list(names(merged_sig)), function(x){
  df = merged_sig[[x]]
  print(x)
  df <- df[!is.na(df$label), ]
  df$tag <- 1
  df[df$label == "NR",]$tag = 0
  df$tag <- factor(df$tag,levels = c(0,1), order=TRUE)
  return(df)
})
names(merged_sig1) <- names(merged_sig)
merged_sig1 <- do.call(rbind,merged_sig1)

rownames(merged_sig1) <- paste0(merged_sig1$dataset, merged_sig1$patient)
res_df <-merged_sig1[!is.na(merged_sig1$label), ]
res_df = res_df[order(res_df$label), ]

res_df = res_df[-c(grep("mesa", names(res_df)),
                   grep("melanocytes", names(res_df)),
                   grep("uncharacter", names(res_df)),
                   grep("astro", names(res_df)),
                   grep("hsc", names(res_df)),
                   grep("mep", names(res_df)),
                   grep("mus", names(res_df)), 
                   grep("mapki", names(res_df)), 
                   grep("melanoma.cell", names(res_df)), 
                   grep("az_", names(res_df)),
                   grep("ips_ips", names(res_df)),
                   grep("tide_tide", names(res_df)))]
res_df[is.na(res_df)] = 0

res_df <- res_df[res_df$cancer == cancer, ]

# ------------------------------------------------------------------------------
# new data ---------------------------------------------------------------------
# ------------------------------------------------------------------------------

# dataset = "Liu_dataset_CTLA4Naive"
# cancer = "SKCM"
sample_type = "Metastatic"

result_dir = paste0("05_immuneSig_ICBresponse/validation_results/", dataset,"/")
data_test <- read.csv(paste0(result_dir, "/", cancer, "_", sample_type, "_immuneSig_merged_mlr.csv"))
data_test_name <- data.frame(name = names(data_test), immune_sig = tolower(names(data_test)))
data_test_name <- immunesig_name_unify(data_test_name)
names(data_test) <- data_test_name$immune_sig
data_test$tag <- 1
data_test[data_test$label == "NR",]$tag = 0
data_test$tag <- factor(data_test$tag,levels = c(1,0), order=TRUE)

# ------------------------------------------------------------------------------
# feature selection ------———---------------------------------------------------
# and create directory according to strategy and model -------------------------
# ------------------------------------------------------------------------------
resultpath = paste0("05_immuneSig_ICBresponse/validation_results/oriImmuSig_mlr_model/",cancer,"_unscaled/")
dir.create(resultpath)

if(featureFilterStrategy == "mlrFilter"){
  if(mlrFilter %in% c('FSelector_gain.ratio',
                      'FSelector_chi.squared',
                      'FSelector_information.gain',
                      'FSelector_symmetrical.uncertainty')){
    dir.create(paste0(resultpath, "_", mlrFilter ,"/"))
    resultpath = paste0(resultpath, "_",mlrFilter ,"/")
    dir.create(paste0(resultpath, "/rf"))
    dir.create(paste0(resultpath, "/en"))
    dir.create(paste0(resultpath, "/ridge"))
    dir.create(paste0(resultpath, "/svm"))
    dir.create(paste0(resultpath, "/xgboost"))
    
    features_selected <- feature_filter(data = res_df[-c(1:5)],
                                        strategy = featureFilterStrategy,
                                        #   p_value = 0.05,
                                        cancer = cancer, 
                                        mlrFilter = mlrFilter,
                                        mlrFilter_k = NULL,
                                        workpath = workpath) 
    
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
    dir.create(paste0(resultpath, "/rf"))
    dir.create(paste0(resultpath, "/en"))
    dir.create(paste0(resultpath, "/ridge"))
    dir.create(paste0(resultpath, "/svm"))
    dir.create(paste0(resultpath, "/xgboost"))
    
    features_selected <- feature_filter(data = res_df[-c(1:5)],
                                        strategy = featureFilterStrategy,
                                        #  p_value = 0.05,
                                        cancer = cancer, 
                                        mlrFilter = mlrFilter,
                                        mlrFilter_k = mlrFilter_k,
                                        workpath = workpath) 
    
  }
}else if(featureFilterStrategy == "mixed"){
  if(mlrFilter %in% c('FSelector_gain.ratio',
                      'FSelector_chi.squared',
                      'FSelector_information.gain',
                      'FSelector_symmetrical.uncertainty')){
    dir.create(paste0(resultpath, "_",featureFilterStrategy,"_", mlrFilter ,"/"))
    resultpath = paste0(resultpath, "_",featureFilterStrategy, "_",mlrFilter ,"/")
    dir.create(paste0(resultpath, "/rf"))
    dir.create(paste0(resultpath, "/en"))
    dir.create(paste0(resultpath, "/ridge"))
    dir.create(paste0(resultpath, "/svm"))
    dir.create(paste0(resultpath, "/xgboost"))
    
    features_selected <- feature_filter(data = res_df[-c(1:5)],
                                        strategy = featureFilterStrategy,
                                        #   p_value = 0.05,
                                        cancer = cancer, 
                                        mlrFilter = mlrFilter,
                                        mlrFilter_k = NULL,
                                        workpath = workpath) 
    
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
    dir.create(paste0(resultpath, "/rf"))
    dir.create(paste0(resultpath, "/en"))
    dir.create(paste0(resultpath, "/ridge"))
    dir.create(paste0(resultpath, "/svm"))
    dir.create(paste0(resultpath, "/xgboost"))
    
    features_selected <- feature_filter(data = res_df[-c(1:5)],
                                        strategy = featureFilterStrategy,
                                        #  p_value = 0.05,
                                        cancer = cancer, 
                                        mlrFilter = mlrFilter,
                                        mlrFilter_k = mlrFilter_k,
                                        workpath = workpath) 
    
  }
}else{
  dir.create(paste0(resultpath, "/"))
  # resultpath = paste0(resultpath, "/")
  resultpath = paste0(resultpath, "/", featureFilterStrategy, "/")
  dir.create(paste0(resultpath, "/"))
  
  dir.create(paste0(resultpath, "/rf"))
  dir.create(paste0(resultpath, "/en"))
  dir.create(paste0(resultpath, "/ridge"))
  dir.create(paste0(resultpath, "/svm"))
  dir.create(paste0(resultpath, "/xgboost"))
  
  features_selected <- feature_filter(data,
                                      strategy = featureFilterStrategy,
                                      p_value = 0.05,
                                      cancer = cancer, 
                                      mlrFilter = NULL,
                                      mlrFilter_k = NULL,
                                      workpath = workpath) 
  
}

# ------------------------------------------------------------------------------
# seperate dataset
# ------------------------------------------------------------------------------
data_train = res_df
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


# ------------------------------------------------------------------------------
# tuning parameter
# ------------------------------------------------------------------------------

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
mlr::performance(pred_rf, mlr::auc)
# pred$data
pred2 <- ROCR::prediction(pred_rf$data$prob.0 , pred_rf$data$truth)
aucPerf <- ROCR::performance( pred2, 'auc' )
AUCValue<-aucPerf@y.values[[1]]
AUCValue


modelroc <- pROC::roc( pred_rf$data$truth, pred_rf$data$prob.0 )
pdf(paste0(resultpath, "/rf/ROC.pdf"))
plot(modelroc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),
     grid.col=c("green", "red"), max.auc.polygon=TRUE,
     auc.polygon.col="lightblue", print.thres=TRUE)
dev.off()

save(ranger_opt_ps_res, mod_ranger_opt, m_rf, pred_rf, 
     file = paste0(resultpath, "/rf/ranger_model_pred.Rdata"))


# ------------------------------------------------------------------------------
# elastic net from glmnet ------------------------------------------------------
# ------------------------------------------------------------------------------

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
mlr::performance(pred_en, mlr::auc)

pred_en$data
pred2 <- ROCR::prediction(pred_en$data$prob.0 , pred_en$data$truth)
aucPerf <- ROCR::performance( pred2, 'auc' )
AUCValue<-aucPerf@y.values[[1]]
AUCValue


modelroc <- pROC::roc( pred_en$data$truth, pred_en$data$prob.0 )
pdf(paste0(resultpath, "/en/ROC.pdf"))
plot(modelroc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),
     grid.col=c("green", "red"), max.auc.polygon=TRUE,
     auc.polygon.col="lightblue", print.thres=TRUE)
dev.off()

save(elasticnet_opt_ps_res, mod_elasticnet_opt, m_en, pred_en, 
     file = paste0(resultpath, "/en/elasticnet_model_pred.Rdata"))

# ------------------------------------------------------------------------------
# elastic net from glmnet ------------------------------------------------------
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
mlr::performance(pred_ridge, mlr::auc)

pred_ridge$data
pred2 <- ROCR::prediction(pred_ridge$data$prob.0 , pred_ridge$data$truth)
aucPerf <- ROCR::performance( pred2, 'auc' )
AUCValue<-aucPerf@y.values[[1]]
AUCValue


modelroc <- pROC::roc( pred_ridge$data$truth, pred_ridge$data$prob.0 )
pdf(paste0(resultpath, "/ridge/ROC.pdf"))
plot(modelroc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),
     grid.col=c("green", "red"), max.auc.polygon=TRUE,
     auc.polygon.col="lightblue", print.thres=TRUE)
dev.off()

save(ridge_opt_ps_res, mod_ridge_opt, m_ridge, pred_ridge, 
     file = paste0(resultpath, "/ridge/ridge_model_pred.Rdata"))


# ksvm  ------------------------------------------------------

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
mlr::performance(pred_svm, mlr::auc)
pred_svm$data
pred2 <- ROCR::prediction(pred_svm$data$prob.0 , pred_svm$data$truth)
aucPerf <- ROCR::performance( pred2, 'auc' )
AUCValue<-aucPerf@y.values[[1]]
AUCValue

modelroc <- pROC::roc( pred_svm$data$truth, pred_svm$data$prob.0 )
pdf(paste0(resultpath, "/svm/ROC.pdf"))
plot(modelroc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),
     grid.col=c("green", "red"), max.auc.polygon=TRUE,
     auc.polygon.col="lightblue", print.thres=TRUE)
dev.off()

save(ksvm_opt_ps_res, mod_ksvm_opt, m_svm, pred_svm, 
     file = paste0(resultpath, "/svm/svm_model_pred.Rdata"))


# xgboost ----------------------------------------------------------------------

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
mlr::performance(pred_xgboost, mlr::auc)
pred_xgboost$data
pred2 <- ROCR::prediction(pred_xgboost$data$prob.0 , pred_xgboost$data$truth)
aucPerf <- ROCR::performance( pred2, 'auc' )
AUCValue<-aucPerf@y.values[[1]]
AUCValue

modelroc <- pROC::roc( pred_xgboost$data$truth, pred_xgboost$data$prob.0 )
pdf(paste0(resultpath, "/xgboost/ROC.pdf"))
plot(modelroc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),
     grid.col=c("green", "red"), max.auc.polygon=TRUE,
     auc.polygon.col="lightblue", print.thres=TRUE)
dev.off()

save(xgboost_opt_ps_res, mod_xgboost_opt, m_xgboost, pred_xgboost, 
     file = paste0(resultpath, "/xgboost/xgboost_model_pred.Rdata"))



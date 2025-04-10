# rm(list=ls())
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
  "mlrFilter_k",            "k", 2, "numeric", "top k persent feature, k from 0-1",
  "dataset",                "d", 2, "character", "validation dataset",
  "help",                   "h", 0, "logical",   "help file"),
  byrow=T,ncol=5)

args=getopt(command)

featureFilterStrategy = args$featureFilterStrategy
cancer = args$cancer
mlrFilter = args$mlrFilter
mlrFilter_k = args$mlrFilter_k
dataset = args$dataset

# featureFilterStrategy = "metaAnalysis" # c("metaAnalysis", "combinedP", "voting", "mlrFilter"
# mlrFilter = 'auc'
# mlrFilter_k = 0.1
# cancer

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
merged_sig_scaled <- lapply(as.list(names(merged_sig)), function(x){
  df = merged_sig[[x]]
  print(x)
  df_scaled <- apply(df[-c(1:5)], 2, scale)
  df_scaled[is.nan(df_scaled)] <- 0
  df_scaled <- cbind(df[1:5], df_scaled)
  df_scaled <- df_scaled[!is.na(df_scaled$label), ]
  df_scaled$tag <- 1
  df_scaled[df_scaled$label == "NR",]$tag = 0
  df_scaled$tag <- factor(df_scaled$tag,levels = c(0,1), order=TRUE)
  return(df_scaled)
})
names(merged_sig_scaled) <- names(merged_sig)
merged_sig_scaled <- do.call(rbind,merged_sig_scaled)

rownames(merged_sig_scaled) <- paste0(merged_sig_scaled$dataset, merged_sig_scaled$patient)
res_df <-merged_sig_scaled[!is.na(merged_sig_scaled$label), ]
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
xcell_remove <- xcell_sig_remover()
res_df <- res_df[-which(is.element(names(res_df), tolower(xcell_remove)))]

res_df[is.na(res_df)] = 0

res_df <- res_df[res_df$cancer == cancer, ]

# ------------------------------------------------------------------------------
# feature selection ------———---------------------------------------------------
# and create directory according to strategy and model -------------------------
# ------------------------------------------------------------------------------

resultpath = paste0("05_immuneSig_ICBresponse/validation_results/oriImmuSig_mlr_model/",cancer,"_zscaled/")
dir.create(resultpath)

if(featureFilterStrategy == "mlrFilter"){
  if(mlrFilter %in% c('FSelector_gain.ratio',
                      'FSelector_chi.squared',
                      'FSelector_information.gain',
                      'FSelector_symmetrical.uncertainty')){
    dir.create(paste0(resultpath, "_",mlrFilter ,"/"))
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
                                        #   p_value = 0.05,
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
    dir.create(paste0(resultpath, "_",featureFilterStrategy, "_", mlrFilter ,"/"))
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
}else if(featureFilterStrategy == "meta_mixed"){
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
  #  resultpath = paste0(resultpath, "/")
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

print(length(features_selected))
print((features_selected))
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
save(ranger_opt_ps_res, mod_ranger_opt, m_rf, 
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

save(elasticnet_opt_ps_res, mod_elasticnet_opt, m_en, 
     file = paste0(resultpath, "/en/en_model_pred.Rdata"))

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

save(ridge_opt_ps_res, mod_ridge_opt, m_ridge, 
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
save(ksvm_opt_ps_res, mod_ksvm_opt, m_svm, 
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

save(xgboost_opt_ps_res, mod_xgboost_opt, m_xgboost, 
     file = paste0(resultpath, "/xgboost/xgboost_model_pred.Rdata"))



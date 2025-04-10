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
library(glmnet)

# Loading relied functions -----------------------------------------------------

# mlr --------------------------------------------------------------------------
source("05_immuneSig_ICBresponse/scripts/mlr_function/rf_learner.R")
source("05_immuneSig_ICBresponse/scripts/mlr_function/elasticnet_learner.R")
source("05_immuneSig_ICBresponse/scripts/mlr_function/ksvm_learner.R")
source("05_immuneSig_ICBresponse/scripts/mlr_function/lasso_learner.R")
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

featureFilterStrategy = "metaAnalysis" 

dir.create(paste0(resultpath, "/"))
resultpath = paste0(resultpath, "/", featureFilterStrategy, "/")
dir.create(paste0(resultpath, "/"))
dir.create(paste0(resultpath, "/glmnet_lasso"))

features_selected <- feature_filter(data,
                                    strategy = featureFilterStrategy,
                                    p_value = 0.05,
                                    cancer = NULL, 
                                    mlrFilter = NULL,
                                    mlrFilter_k = NULL,
                                    workpath = workpath) 


print(length(features_selected))
# print((features_selected))


# load data --------------------------------------------------------------------------

meta_result <- read.csv("05_immuneSig_ICBresponse/results/meta_results/res_metaRes.csv")
meta_selected <- meta_result[meta_result$Meta_Pval < 0.05, ]
meta_selected$immune_sig = meta_selected$X

result_dir <- "05_immuneSig_ICBresponse/results_r0.4_ITSssgsea_CombatRef/"

load( paste0(result_dir, "ITS_afterCombat17datasets_ssgseanorm.Rdata"))
datasets = unique(ITSres_selected$dataset)


# data seperation --------------------------------------------------------------------------

vlSetName <- c("08_Lauss_dataset_outcome", "GSE115821_pre_outcome", "GSE96619_outcome","Liu_dataset_CTLA4Naive_outcome", "GSE176307")
dcSetName <- datasets[!is.element(datasets, vlSetName)]

dcSet <- ITSres_selected[is.element(ITSres_selected$dataset, dcSetName),]
vlSet <- ITSres_selected[is.element(ITSres_selected$dataset, vlSetName),]

dcSet_Annot <- dcSet[c('patientID','SampleID', 'cancer','sample_type','dataset')]
dcSet_ITS <- dcSet[-which(is.element(names(dcSet), names(dcSet_Annot)))]

vlSet_Annot <- vlSet[c('patientID','SampleID', 'cancer','sample_type','dataset')]
vlSet_ITS <- vlSet[-which(is.element(names(vlSet), names(vlSet_Annot)))]
# ------------------------------------------------------------------------------
# 5-CV--------------------------------------------------------------------------

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

  lasso.fit <- glmnet(data_train_filter[-which(names(data_train_filter) == 'tag')], as.numeric(data_train_filter$tag), alpha = 0)
  cv.lasso <- cv.glmnet(as.matrix(data_train_filter[-which(names(data_train_filter) == 'tag')]), as.numeric(data_train_filter$tag), alpha =0)
  
  lasso_model_para = cv.lasso$glmnet.fit
  m_lasso_imp <- data.frame(lasso_model_para$beta[,which(lasso_model_para$lambda == cv.lasso$lambda.min)])
  names(m_lasso_imp) <- 'feature_weight'
  m_lasso_imp$feature_name <- rownames(m_lasso_imp)
  
  bestlam <- cv.lasso$lambda.min
#   cat("\nBestlam (with grid):",bestlam)
  
  pred <- predict(lasso.fit, s = bestlam, newx= as.matrix(data_test[-which(names(data_train_filter) == 'tag')]))
  pred <- cbind(pred, data_test[which(names(data_train_filter) == 'tag')])
  names(pred) <- c('prob.1','tag')
  
  pred2 <- ROCR::prediction(pred$prob.1 , pred$tag)
  aucPerf <- ROCR::performance( pred2, 'auc' )
  AUCValue<-aucPerf@y.values[[1]]
  AUCValue
  
  
  pdf(paste0(resultpath, "/glmnet_lasso/ROC_",i,".pdf"))
  modelroc <- pROC::roc(pred$tag, pred$prob.1)
  plot(modelroc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),
       grid.col=c("green", "red"), max.auc.polygon=TRUE,
       auc.polygon.col="lightblue", print.thres=TRUE)
  dev.off()
  
  save(cv.lasso, lasso.fit, pred, 
     file = paste0(resultpath, "/glmnet_lasso/glmnet_lasso_model_pred_",i,".Rdata"))

})




# ------------------------------------------------------------------------------
# FINAL MODEL ------------------------------------------------------------------
# ------------------------------------------------------------------------------
# seperate dataset
# ------------------------------------------------------------------------------
set.seed(1234)
library(sampling)
# library(caret)  
result_dir <- "05_immuneSig_ICBresponse/results_r0.4_ITSssgsea_CombatRef/"


resultpath = paste0(result_dir, "finalmodel_17datasets_ssgseanorm/")
dir.create(resultpath)

featureFilterStrategy = "metaAnalysis" 

dir.create(paste0(resultpath, "/"))
resultpath = paste0(resultpath, "/", featureFilterStrategy, "/")
dir.create(paste0(resultpath, "/"))
dir.create(paste0(resultpath, "/glmnet_lasso"))

res_df <- dcSet_ITS[!is.na(dcSet_ITS$label), ]
res_df = res_df[order(res_df$label), ]
res_df$label <- factor(res_df$label,levels = c("NR","R"), order=TRUE)
res_df$tag <- 1
res_df[res_df$label == "NR",]$tag = 0
res_df$tag <- factor(res_df$tag,levels = c(0, 1), order=TRUE)

# vlSet_ITS
vl_df <- vlSet_ITS[!is.na(vlSet_ITS$label), ]
vl_df = vl_df[order(vl_df$label), ]
vl_df$label <- factor(vl_df$label,levels = c("NR","R"), order=TRUE)
vl_df$tag <- 1
vl_df[vl_df$label == "NR",]$tag = 0
vl_df$tag <- factor(vl_df$tag,levels = c(0, 1), order=TRUE)


data_train = res_df
data_test = vl_df
#   rownames(data_test) <- paste0(data_test$dataset,"_",data_test$patient)

# ------------------------------------------------------------------------------
# create task
# ------------------------------------------------------------------------------
data_train_filter <- data_train[c(intersect(features_selected, names(data_train)), 'tag')]
dim(data_train_filter)
# ------------------------------------------------------------------------------
data_test <- data_test[names(data_train_filter)]



# Again but without the grid (allowing R to figure lambda out
lasso.fit <- glmnet(data_train_filter[-which(names(data_train_filter) == 'tag')], as.numeric(data_train_filter$tag), alpha = 1)
cv.lasso <- cv.glmnet(as.matrix(data_train_filter[-which(names(data_train_filter) == 'tag')]), as.numeric(data_train_filter$tag), alpha = 1)

lasso_model_para = cv.lasso$glmnet.fit
m_lasso_imp <- data.frame(lasso_model_para$beta[,which(lasso_model_para$lambda == cv.lasso$lambda.min)])
names(m_lasso_imp) <- 'feature_weight'
m_lasso_imp$feature_name <- rownames(m_lasso_imp)

bestlam <- cv.lasso$lambda.min
cat("\nBestlam (with grid):",bestlam)

pred <- predict(lasso.fit, s = bestlam, newx= as.matrix(data_test[-which(names(data_test) == 'tag')]))
pred <- cbind(pred, data_test[which(names(data_train_filter) == 'tag')])
names(pred) <- c('prob.1','tag')

pred2 <- ROCR::prediction(pred$prob.1 , pred$tag)
aucPerf <- ROCR::performance( pred2, 'auc' )
AUCValue<-aucPerf@y.values[[1]]
AUCValue


modelroc <- pROC::roc(pred$tag, pred$prob.1)
pdf(paste0(resultpath, "/glmnet_lasso/ROC.pdf"))
plot(modelroc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),
     grid.col=c("green", "red"), max.auc.polygon=TRUE,
     auc.polygon.col="lightblue", print.thres=TRUE)
dev.off()

save(lasso.fit, cv.lasso, pred, 
     file = paste0(resultpath, "/glmnet_lasso/glmnet_lasso_model_pred.Rdata"))

load(paste0(resultpath, "/glmnet_lasso/glmnet_lasso_model_pred.Rdata"))
lasso_model_para = cv.lasso$glmnet.fit
m_lasso_imp <- data.frame(lasso_model_para$beta[,which(lasso_model_para$lambda == cv.lasso$lambda.min)])

names(m_lasso_imp) <- 'feature_weight'
m_lasso_imp$feature_name <- rownames(m_lasso_imp)

m_lasso_imp[order(m_lasso_imp$feature_weight),]
m_lasso_imp[m_lasso_imp$feature_weight != 0, ]


res_auc <- lapply(as.list(vlSetName), function(dset){
    # dset = "08_Lauss_dataset_outcome"
    print(dset)
    tmp = pred[grep(dset, rownames(pred)),]
    pred2 <- ROCR::prediction(tmp$prob.1 , tmp$tag)
    aucPerf <- ROCR::performance( pred2, 'auc' )
    AUCValue<-aucPerf@y.values[[1]]
    res = data.frame(datasets = dset, auc = AUCValue)
    return(res)

    # modelroc <- pROC::roc(tmp$truth, tmp$prob.1)
    # res = data.frame(datasets = dset, auc = modelroc$auc[1])
    # return(res)

})
do.call(rbind, res_auc)

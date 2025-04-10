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

featureFilterStrategy = "metaAnalysis_filter_removeprior_v2" 

dir.create(paste0(resultpath, "/"))
resultpath = paste0(resultpath, "/", featureFilterStrategy, "/")
dir.create(paste0(resultpath, "/"))
dir.create(paste0(resultpath, "/glmnet_ridge"))

sig_select <- read.csv("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/07_plotting_v2/01_metaAnalysis/table/metaSelected_immunesigAnnot_subtypedetail_filter.txt", sep = '\t')
sig_filter_prior <- c("exhausted_immucellai", 
                      "merck18_tide",
                      "oe_tme.t.cd4.treg",
                      "mdsc_tcga_caf_immunecls"
                    #   "oe_co.culture.screen.hits10",
                    #   "oe_in.vivo.screen.gvax.vs.tcrako.enriched"
                      )

features_selected <- sig_select[!is.element(sig_select$immune_sig, sig_filter_prior), ]$immune_sig


print(length(features_selected))
# print((features_selected))


# load data --------------------------------------------------------------------------

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

  ridge.fit <- glmnet(data_train_filter[-which(names(data_train_filter) == 'tag')], as.numeric(data_train_filter$tag), alpha = 0)
  cv.ridge <- cv.glmnet(as.matrix(data_train_filter[-which(names(data_train_filter) == 'tag')]), as.numeric(data_train_filter$tag), alpha =0)
  
  ridge_model_para = cv.ridge$glmnet.fit
  m_ridge_imp <- data.frame(ridge_model_para$beta[,which(ridge_model_para$lambda == cv.ridge$lambda.min)])
  names(m_ridge_imp) <- 'feature_weight'
  m_ridge_imp$feature_name <- rownames(m_ridge_imp)
  
  bestlam <- cv.ridge$lambda.min
#   cat("\nBestlam (with grid):",bestlam)
  
  pred <- predict(ridge.fit, s = bestlam, newx= as.matrix(data_test[-which(names(data_train_filter) == 'tag')]))
  pred <- cbind(pred, data_test[which(names(data_train_filter) == 'tag')])
  names(pred) <- c('prob.1','tag')
  
  pred2 <- ROCR::prediction(pred$prob.1 , pred$tag)
  aucPerf <- ROCR::performance( pred2, 'auc' )
  AUCValue<-aucPerf@y.values[[1]]
  AUCValue
  
  
  pdf(paste0(resultpath, "/glmnet_ridge/ROC_",i,".pdf"))
  modelroc <- pROC::roc(pred$tag, pred$prob.1)
  plot(modelroc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),
       grid.col=c("green", "red"), max.auc.polygon=TRUE,
       auc.polygon.col="lightblue", print.thres=TRUE)
  dev.off()
  
  save(cv.ridge, ridge.fit, pred, 
     file = paste0(resultpath, "/glmnet_ridge/glmnet_ridge_model_pred_",i,".Rdata"))

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

featureFilterStrategy = "metaAnalysis_filter_removeprior_v2" 

dir.create(paste0(resultpath, "/"))
resultpath = paste0(resultpath, "/", featureFilterStrategy, "/")
dir.create(paste0(resultpath, "/"))
dir.create(paste0(resultpath, "/glmnet_ridge"))

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
ridge.fit <- glmnet(data_train_filter[-which(names(data_train_filter) == 'tag')], as.numeric(data_train_filter$tag), alpha = 0)
cv.ridge <- cv.glmnet(as.matrix(data_train_filter[-which(names(data_train_filter) == 'tag')]), as.numeric(data_train_filter$tag), alpha =0)

ridge_model_para = cv.ridge$glmnet.fit
m_ridge_imp <- data.frame(ridge_model_para$beta[,which(ridge_model_para$lambda == cv.ridge$lambda.min)])
names(m_ridge_imp) <- 'feature_weight'
m_ridge_imp$feature_name <- rownames(m_ridge_imp)

bestlam <- cv.ridge$lambda.min
cat("\nBestlam (with grid):",bestlam)

pred <- predict(ridge.fit, s = bestlam, newx= as.matrix(data_test[-which(names(data_test) == 'tag')]))
pred <- cbind(pred, data_test[which(names(data_train_filter) == 'tag')])
names(pred) <- c('prob.1','tag')

pred2 <- ROCR::prediction(pred$prob.1 , pred$tag)
aucPerf <- ROCR::performance( pred2, 'auc' )
AUCValue<-aucPerf@y.values[[1]]
AUCValue


modelroc <- pROC::roc(pred$tag, pred$prob.1)
pdf(paste0(resultpath, "/glmnet_ridge/ROC.pdf"))
plot(modelroc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),
     grid.col=c("green", "red"), max.auc.polygon=TRUE,
     auc.polygon.col="lightblue", print.thres=TRUE)
dev.off()

save(ridge.fit, cv.ridge, pred, 
     file = paste0(resultpath, "/glmnet_ridge/glmnet_ridge_model_pred.Rdata"))

load(paste0(resultpath, "/glmnet_ridge/glmnet_ridge_model_pred.Rdata"))
ridge_model_para = cv.ridge$glmnet.fit
m_ridge_imp <- data.frame(ridge_model_para$beta[,which(ridge_model_para$lambda == cv.ridge$lambda.min)])
names(m_ridge_imp) <- 'feature_weight'
m_ridge_imp$feature_name <- rownames(m_ridge_imp)

m_ridge_imp[order(m_ridge_imp$feature_weight),]



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

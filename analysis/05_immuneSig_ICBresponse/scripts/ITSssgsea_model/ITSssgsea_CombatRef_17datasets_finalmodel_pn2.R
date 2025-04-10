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

featureFilterStrategy = "metaAnalysis" 
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
source("05_immuneSig_ICBresponse/scripts/mlr_function/naive_bayes_learner.R")
source("05_immuneSig_ICBresponse/scripts/oriImmuSig_mlr_cancer/feature_filter_function.R")
source("06_results_summaryandplotting/scripts/immunesig_selection_function.R")

library(mlr)
library(mlrMBO)

# ------------------------------------------------------------------------------
# feature selection ------———---------------------------------------------------
# and create directory according to strategy and model -------------------------
# ------------------------------------------------------------------------------
result_dir <- "05_immuneSig_ICBresponse/results_r0.4_ITSssgsea_CombatRef/"


resultpath = paste0(result_dir, "finalmodel_17datasets_ssgseanorm/")
dir.create(resultpath)


  
  dir.create(paste0(resultpath, "/"))
  resultpath = paste0(resultpath, "/", featureFilterStrategy, "/")
  dir.create(paste0(resultpath, "/"))
  dir.create(paste0(resultpath, "/en_pn2"))


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

load(paste0(result_dir, "ITS_afterCombat17datasets_ssgseanorm_pn.Rdata"))


datasets = unique(ITSres_selected$dataset)

# pdf(paste0(result_dir, "ITS_afterCombat_ssgseanorm_p_heatmap.pdf"), width = 12, height = 8)
# for(i in 1:length(datasets)){
# ds <- ITSres_selected_p[ITSres_selected_p$dataset == datasets[i],]
# ds <- ds[order(ds$label),]
# ds_Annot <- ds[c('patientID','label','SampleID',
#                      'cancer','sample_type','dataset')]
# ds_selected2 <- ds[-which(is.element(names(ds), names(ds_Annot)))]
# pheatmap::pheatmap(ds_selected2, scale = 'column',
#                    cluster_row = F, cluster_col = F,
#                    annotation_row = ds_Annot['label'])
# }
# dev.off()

# pdf(paste0(result_dir, "ITS_afterCombat_ssgseanorm_n_heatmap.pdf"), width = 12, height = 8)
# for(i in 1:length(datasets)){
# ds <- ITSres_selected_n[ITSres_selected_n$dataset == datasets[i],]
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

names(ITSres_selected_p) <- c(names(ITSres_selected_p)[c(1:6)], paste0(names(ITSres_selected_p)[-c(1:6)], "_p"))
names(ITSres_selected_n) <- c(names(ITSres_selected_n)[c(1:6)],  paste0(names(ITSres_selected_n)[-c(1:6)], "_n"))
ITSres_selected <- inner_join(ITSres_selected_p,ITSres_selected_n)
rownames(ITSres_selected) <- ITSres_selected$SampleID

vlSetName <- c("08_Lauss_dataset_outcome", "GSE115821_pre_outcome", "GSE96619_outcome","Liu_dataset_CTLA4Naive_outcome","GSE176307")
dcSetName <- datasets[!is.element(datasets, vlSetName)]

dcSet <- ITSres_selected[is.element(ITSres_selected$dataset, dcSetName),]
vlSet <- ITSres_selected[is.element(ITSres_selected$dataset, vlSetName),]

dcSet_Annot <- dcSet[c('patientID','SampleID', 'cancer','sample_type','dataset')]
dcSet_ITS <- dcSet[-which(is.element(names(dcSet), names(dcSet_Annot)))]

vlSet_Annot <- vlSet[c('patientID','SampleID', 'cancer','sample_type','dataset')]
vlSet_ITS <- vlSet[-which(is.element(names(vlSet), names(vlSet_Annot)))]
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

data_train_filter <- data_train[-1]# [c(intersect(features_selected, names(data_train)), 'tag')]

filtered.task = makeClassifTask(id = "ICBresponse", 
                                data = data_train_filter, 
                                target = "tag", 
                                positive = 1)

print(filtered.task)
# ------------------------------------------------------------------------------
# feature selection for test set 
# ------------------------------------------------------------------------------

data_test <- data_test[c(getTaskFeatureNames(filtered.task), 'tag')]
rownames(data_test) <- rownames(vl_df)
ICBresponse_testtask = makeClassifTask(id = "ICBresponse", 
                                       data = data_test, 
                                       target = "tag", 
                                       positive = 1)
print(ICBresponse_testtask)

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

# pred_en$data
pred2 <- ROCR::prediction(pred_en$data$prob.1 , pred_en$data$truth)
aucPerf <- ROCR::performance( pred2, 'auc' )
AUCValue<-aucPerf@y.values[[1]]
AUCValue


modelroc <- pROC::roc( pred_en$data$truth, pred_en$data$prob.1 )
pdf(paste0(resultpath, "/en_pn2/ROC.pdf"))
plot(modelroc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),
     grid.col=c("green", "red"), max.auc.polygon=TRUE,
     auc.polygon.col="lightblue", print.thres=TRUE)
dev.off()

save(elasticnet_opt_ps_res, mod_elasticnet_opt, m_en, pred_en, 
     file = paste0(resultpath, "/en_pn2/elasticnet_model_pred.Rdata"))

load( paste0(resultpath, "/en_pn2/elasticnet_model_pred.Rdata"))

respred = pred_en$data
res_auc <- lapply(as.list(vlSetName), function(dset){
    # dset = "GSE115821_pre_outcome"
    print(dset)
    tmp = respred[grep(dset, rownames(respred)),]
    pred2 <- ROCR::prediction(tmp$prob.1 , tmp$truth)
    aucPerf <- ROCR::performance( pred2, 'auc' )
    AUCValue<-aucPerf@y.values[[1]]
    res = data.frame(datasets = dset, auc = AUCValue)
    return(res)

    # modelroc <- pROC::roc(tmp$tag, tmp$prob.1)
    # res = data.frame(datasets = dset, auc = modelroc$auc[1])
    # return(res)

})

do.call(rbind, res_auc)


library(ggplot2)

ds = pred_en$data[grep("Liu", rownames(pred_en$data)),]
wilcox.test(ds[ds$truth == 1,]$prob.1, ds[ds$truth == 0,]$prob.1, alternative = 'greater')

ggplot(data=ds, mapping=aes(x=truth, y=prob.1))+
geom_boxplot()



en_model_para = m_en$learner.model$glmnet.fit
m_en_imp <- data.frame(en_model_para$beta[,which(en_model_para$lambda == m_en$learner.model$lambda.min)])
names(m_en_imp) <- 'en_weight'
m_en_imp$immune_sig <- rownames(m_en_imp)
sigAnnot2$immune_sig <- paste0(sigAnnot2$X, '_p')
m_en_imp1 <- inner_join(m_en_imp, sigAnnot2)
sigAnnot2$immune_sig <- paste0(sigAnnot2$X, '_n')
m_en_imp2 <- inner_join(m_en_imp, sigAnnot2)
m_en_imp <- rbind(m_en_imp1,m_en_imp2)

plot(m_en_imp$en_weight,m_en_imp$OR)
abline(v=0,h=1)

m_en_imp <- m_en_imp[m_en_imp$en_weight !=0,]
m_en_imp <- m_en_imp[order(m_en_imp$immune_sig),]
m_en_imp[grep('cd8',m_en_imp$immune_sig),]
View(m_en_imp[order(m_en_imp$en_weight),])


m_en_imp_index <- aggregate(m_en_imp[c('en_weight')], by = list((m_en_imp$Index)), mean)
m_en_imp_type <- aggregate(m_en_imp[c('en_weight')], by = list((m_en_imp$Type)), mean)
m_en_imp_index[order(m_en_imp_index$en_weight),]
m_en_imp_type[order(m_en_imp_type$en_weight),]
m_en_imp[grep('mdsc',m_en_imp$immune_sig),]

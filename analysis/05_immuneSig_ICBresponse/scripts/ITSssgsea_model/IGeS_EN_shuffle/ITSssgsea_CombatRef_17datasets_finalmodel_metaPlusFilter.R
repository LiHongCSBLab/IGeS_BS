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
lapply(as.list(1:10), function(rand_index){

     result_dir <- paste0("07_plotting_v2/IGeS_EN_shuffle/rand", rand_index, "/")

     resultpath = paste0(result_dir, "finalmodel_17datasets_ssgseanorm/")
     dir.create(resultpath)
     featureFilterStrategy = "metaAnalysis_filter" 

     dir.create(paste0(resultpath, "/"))
     resultpath = paste0(resultpath, "/", featureFilterStrategy, "/")
     dir.create(paste0(resultpath, "/"))
     #   dir.create(paste0(resultpath, "/rf"))
     dir.create(paste0(resultpath, "/en"))
     #   dir.create(paste0(resultpath, "/ridge"))
     #   dir.create(paste0(resultpath, "/svm"))
     #   dir.create(paste0(resultpath, "/xgboost"))
     #   dir.create(paste0(resultpath, "/naiveBayes/"))


     sig_select <- read.csv("07_plotting_v2/01_metaAnalysis/table/metaSelected_immunesigAnnot_subtypedetail_filter.txt", sep = '\t')

     features_selected <- sig_select$immune_sig

     print(length(features_selected))
     # print((features_selected))


     # load data --------------------------------------------------------------------------
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
     load(paste0(result_dir, "Exp_afterCombatRef.Rdata"))
     load(paste0(result_dir, "rand_cohort_index.Rdata"))

     dcSetName<- trainset1_name# c("08_Lauss_dataset_outcome", "GSE115821_pre_outcome", "GSE96619_outcome","Liu_dataset_CTLA4Naive_outcome","GSE176307")
     vlSetName <- valset1_name

     dcSet <- ITSres_selected[is.element(ITSres_selected$dataset, dcSetName),]
     vlSet <- ITSres_selected[is.element(ITSres_selected$dataset, vlSetName),]

     dcSet_Annot <- dcSet[c('patientID','SampleID', 'cancer','sample_type','dataset')]
     dcSet_ITS <- dcSet[-which(is.element(names(dcSet), names(dcSet_Annot)))]

     vlSet_Annot <- vlSet[c('patientID','SampleID', 'cancer','sample_type','dataset')]
     vlSet_ITS <- vlSet[-which(is.element(names(vlSet), names(vlSet_Annot)))]
     write.csv(dcSet, '07_plotting_v2/02_ITSprediction_metaPlusFilter/dcSet.csv',row.names=F, quote = F)
     write.csv(vlSet, '07_plotting_v2/02_ITSprediction_metaPlusFilter/vlSet.csv',row.names=F, quote = F)
     write.csv(as.data.frame.array(table(dcSet[c('dataset','label')])), '07_plotting_v2/02_ITSprediction_metaPlusFilter/dcSet_summary.csv', quote = F)
     write.csv(as.data.frame.array(table(vlSet[c('dataset','label')])), '07_plotting_v2/02_ITSprediction_metaPlusFilter/vlSet_summary.csv', quote = F)

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
     # elastic net from glmnet ------------------------------------------------------
     # ------------------------------------------------------------------------------

     print("start modeling")
     # tuning parameter

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
                                   s = elasticnet_opt_ps$s,
                                   alpha = elasticnet_opt_ps$alpha
     )


     m_en = train(mod_elasticnet_opt, filtered.task)
     pred_en = predict(m_en, task = ICBresponse_testtask)
     pred_en$data[(pred_en$data$prob.1 >0.2),]$response=1
     calculateConfusionMatrix(pred_en)
     #  mlr::performance(pred_en, mlr::auc)

     # pred_en$data
     pred2 <- ROCR::prediction(pred_en$data$prob.1 , pred_en$data$truth)
     aucPerf <- ROCR::performance( pred2, 'auc' )
     AUCValue<-aucPerf@y.values[[1]]
     AUCValue


     modelroc <- pROC::roc( pred_en$data$truth, pred_en$data$prob.1 )
     pdf(paste0(resultpath, "/en/ROC.pdf"))
     plot(modelroc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),
          grid.col=c("green", "red"), max.auc.polygon=TRUE,
          auc.polygon.col="lightblue", print.thres=TRUE)
     dev.off()

     save(elasticnet_opt_ps_res, mod_elasticnet_opt, m_en, pred_en, 
          file = paste0(resultpath, "/en/elasticnet_model_pred.Rdata"))

     # load( paste0(resultpath, "/en/elasticnet_model_pred.Rdata"))

     respred = pred_en$data
     res_auc <- lapply(as.list(vlSetName), function(dset){
     # dset = "Liu_dataset_CTLA4Naive_outcome"
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

     res_auc_all = do.call(rbind, res_auc)
     write.csv(res_auc_all, paste0(resultpath, "/en/res_auc_all.csv"), quote=F)



     immunesig <- read.csv("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/05_immuneSig_ICBresponse/results/meta_results/res_metaRes.csv")
     immunesig <- immunesig[immunesig$Meta_Pval < 0.05, ]
     names(immunesig) <- c("immune_sig", names(immunesig)[-1])
     immunesig$flag <- "sensitive"
     immunesig[immunesig$OR < 1,]$flag <- "resistant"
     immunesig$setindex = "01_sensitive"
     immunesig[immunesig$flag == "resistant", ]$setindex = "02_resistant"

     # load(paste0(resultpath, "/elasticnet_model_pred.Rdata"))
     en_model_para = m_en$learner.model$glmnet.fit
     m_en_imp <- data.frame(en_model_para$beta[,which(en_model_para$lambda == m_en$learner.model$lambda.min)])
     names(m_en_imp) <- "weight"
     m_en_imp$immune_sig = row.names(m_en_imp)
     immunesig = merge(immunesig, m_en_imp, by = "immune_sig")
     immunesig =immunesig[immunesig$weight != 0 ,]
     write.csv(immunesig, paste0(resultpath, "/en/immunesig_", rand_index, ".csv"), quote=F)

})



# tmp = respred[-grep('96619', rownames(respred)),]
# pred2 <- ROCR::prediction(tmp$prob.1 , tmp$truth)
# aucPerf <- ROCR::performance( pred2, 'auc' )
# AUCValue<-aucPerf@y.values[[1]]
# AUCValue


# library(ggplot2)

# ds = pred_en$data[grep("115", rownames(pred_en$data)),]
# wilcox.test(ds[ds$truth == 1,]$prob.1, ds[ds$truth == 0,]$prob.1, alternative = 'greater')

# ggplot(data=ds, mapping=aes(x=truth, y=prob.1))+
# geom_boxplot()


# en_model_para = m_en$learner.model$glmnet.fit
# m_en_imp <- data.frame(en_model_para$beta[,which(en_model_para$lambda == m_en$learner.model$lambda.min)])
# names(m_en_imp) <- 'en_weight'
# m_en_imp$immune_sig <- rownames(m_en_imp)
# m_en_imp <- inner_join(m_en_imp, immunesig)
# plot(m_en_imp$en_weight, m_en_imp$OR)
# abline(v=0,h=1)

# m_en_imp <- m_en_imp[m_en_imp$en_weight !=0,]
# m_en_imp[m_en_imp$en_weight >0 & m_en_imp$OR >1,]
# m_en_imp[m_en_imp$en_weight <0 & m_en_imp$OR <1,]


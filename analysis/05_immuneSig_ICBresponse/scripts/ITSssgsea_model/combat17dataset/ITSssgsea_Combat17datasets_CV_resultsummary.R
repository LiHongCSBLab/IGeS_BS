# rm(list=ls())
workpath <- "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/"
setwd(workpath)

options(stringsAsFactors = F)
library(glmnet)
library(mlr)

result_dir <- "05_immuneSig_ICBresponse/results_r0.4_ITSssgsea_Combat/"
resultpath = paste0(result_dir, "modelCV_17datasets_ssgseanorm/metaAnalysis/")


load( paste0(result_dir, "ITS_afterCombat_ssgseanorm.Rdata"))
datasets = unique(ITSres_selected$dataset)
vlSetName <- c("08_Lauss_dataset_outcome", "GSE115821_pre_outcome", "GSE96619_outcome","Liu_dataset_CTLA4Naive_outcome")
dcSetName <- datasets[!is.element(datasets, vlSetName)]

discoverset_auc=list()
folds = paste0('Fold', 1:5)

# glmnet_ridge -----------------------------------------------------------------

res_pred <- lapply(as.list(folds), function(i){
    # i = 'Fold1'
    load(paste0(resultpath, "/glmnet_ridge/glmnet_ridge_model_pred_",i,".Rdata"))
    return(pred)
})

res_pred <- do.call(rbind, res_pred)
res_pred$SampleID <- row.names(res_pred)
res_pred <- merge(ITSres_selected[c('SampleID', 'label')], res_pred)

pred2 <- ROCR::prediction(res_pred$prob.1 , res_pred$tag)
aucPerf <- ROCR::performance( pred2, 'auc' )
AUCValue<-aucPerf@y.values[[1]]
AUCValue
summary(res_pred$prob.1)
res_pred$predtag = 0
res_pred[res_pred$prob.1 > 1.286,]$predtag = 1
table(res_pred[c('tag', 'predtag')])


modelroc <- pROC::roc(res_pred$tag, res_pred$prob.1)
plot(modelroc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),
     grid.col=c("green", "red"), max.auc.polygon=TRUE,
     auc.polygon.col="lightblue", print.thres=TRUE)


res_auc <- lapply(as.list(dcSetName), function(dset){
    # dset = "Zhao_dataset"
    print(dset)
    tmp = res_pred[grep(dset, res_pred$SampleID),]
    pred2 <- ROCR::prediction(tmp$prob.1 , tmp$label)
    aucPerf <- ROCR::performance( pred2, 'auc' )
    AUCValue<-aucPerf@y.values[[1]]
    res = data.frame(datasets = dset, auc = AUCValue)
    return(res)

    # modelroc <- pROC::roc(tmp$tag, tmp$prob.1)
    # res = data.frame(datasets = dset, auc = modelroc$auc[1])
    # return(res)

})

discoverset_auc[[1]] <- do.call(rbind, res_auc)
discoverset_auc

# elasticnet -----------------------------------------------------------------

res_pred <- lapply(as.list(folds), function(i){
    # i = 'Fold1'
    load(paste0(resultpath, "/en/elasticnet_model_pred_",i,".Rdata"))
    return(pred_en$data)
})

res_pred <- do.call(rbind, res_pred)
res_pred$SampleID <- row.names(res_pred)
res_pred <- merge(ITSres_selected[c('SampleID', 'label')], res_pred)

pred2 <- ROCR::prediction(res_pred$prob.1 , res_pred$truth)
aucPerf <- ROCR::performance( pred2, 'auc' )
AUCValue<-aucPerf@y.values[[1]]
AUCValue

table(res_pred[c('truth', 'response')])

pdf(paste0(resultpath, "/en/allset_ROC.pdf"))
modelroc <- pROC::roc(res_pred$truth, res_pred$prob.1)
plot(modelroc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),
     grid.col=c("green", "red"), max.auc.polygon=TRUE,
     auc.polygon.col="lightblue", print.thres=TRUE)
dev.off()

res_auc <- lapply(as.list(dcSetName), function(dset){
    # dset = "Zhao_dataset"
    print(dset)
    tmp = res_pred[grep(dset, res_pred$SampleID),]
    pred2 <- ROCR::prediction(tmp$prob.1 , tmp$label)
    aucPerf <- ROCR::performance( pred2, 'auc' )
    AUCValue<-aucPerf@y.values[[1]]
    res = data.frame(datasets = dset, auc = AUCValue)
    return(res)

    # modelroc <- pROC::roc(tmp$truth, tmp$prob.1)
    # res = data.frame(datasets = dset, auc = modelroc$auc[1])
    # return(res)

})

discoverset_auc[[2]] <- do.call(rbind, res_auc)
mean(discoverset_auc[[2]]$auc)

# rf -----------------------------------------------------------------

res_pred <- lapply(as.list(folds), function(i){
    # i = 'Fold1'
    load(paste0(resultpath, "/rf/ranger_model_pred_",i,".Rdata"))
    return(pred_rf$data)
})

res_pred <- do.call(rbind, res_pred)
res_pred$SampleID <- row.names(res_pred)
res_pred <- merge(ITSres_selected[c('SampleID', 'label')], res_pred)

pred2 <- ROCR::prediction(res_pred$prob.1 , res_pred$truth)
aucPerf <- ROCR::performance( pred2, 'auc' )
AUCValue<-aucPerf@y.values[[1]]
AUCValue


modelroc <- pROC::roc(res_pred$truth, res_pred$prob.1)
plot(modelroc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),
     grid.col=c("green", "red"), max.auc.polygon=TRUE,
     auc.polygon.col="lightblue", print.thres=TRUE)


res_auc <- lapply(as.list(dcSetName), function(dset){
    # dset = "Zhao_dataset"
    print(dset)
    tmp = res_pred[grep(dset, res_pred$SampleID),]
    pred2 <- ROCR::prediction(tmp$prob.1 , tmp$label)
    aucPerf <- ROCR::performance( pred2, 'auc' )
    AUCValue<-aucPerf@y.values[[1]]
    res = data.frame(datasets = dset, auc = AUCValue)
    return(res)

    # modelroc <- pROC::roc(tmp$truth, tmp$prob.1)
    # res = data.frame(datasets = dset, auc = modelroc$auc[1])
    # return(res)

})

discoverset_auc[[3]] <- do.call(rbind, res_auc)
mean(discoverset_auc[[3]]$auc)
discoverset_auc[[3]]
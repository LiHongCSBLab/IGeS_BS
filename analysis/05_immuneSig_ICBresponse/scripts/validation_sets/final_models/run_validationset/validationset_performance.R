# rm(list=ls())
setwd("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/05_immuneSig_ICBresponse/")
options(stringsAsFactor = F)

result_paths <- list.files("validation_results/oriImmuSig_mlr_model/")
result_paths <- result_paths[-grep(".csv", result_paths)]
model_types <- c("metaAnalysis","_mixed_ranger_impurity_0.1", "_mixed_ranger_impurity_0.2", "_ranger_impurity_0.1", "_ranger_impurity_0.2",'_meta_mixed_ranger_impurity_0.1','_meta_mixed_ranger_impurity_0.2','_meta_mixed_ranger_impurity_0.3','_meta_mixed_ranger_impurity_0.5')
datasets <- c("Liu_dataset_CTLA4Naive_outcome", "GSE96619_outcome", "GSE115821_pre_outcome", "08_Lauss_dataset_outcome")


res <- lapply(as.list(result_paths), function(resPath){
# for(resPath in result_paths){
    # resPath = "alldatasets_unscaled"
    m_res <- lapply(as.list(model_types), function(m){
        # m = 'metaAnalysis'
        # for(m in model_types){
        d_res <- lapply(as.list(datasets), function(d){
            resPath2 <- paste0("validation_results/oriImmuSig_mlr_model/", resPath, "/", m, "/", d)
            files <- list.files(resPath2)
            files <- files[grep("_model_pred.Rdata", files)]
            # for(f in files){
            f = "ranger_model_pred.Rdata"
            load(paste0(resPath2, "/", f))

            # threshold = c(`1` = 0.4,`0`=0.6)
            # pred = setThreshold(pred_rf, threshold = threshold)
            # mlr::performance(pred, acc)
            # df = generateThreshVsPerfData(list(rf=pred), measures = list(fpr, tpr))
            # generateThreshVsPerfData(pred, measures = list(fpr, tpr))
            # calculateConfusionMatrix(pred_rf)


            modelroc <- pROC::roc( pred_rf$data$truth, pred_rf$data$prob.1 )
            AUCvalue = modelroc$auc[1]
            return(AUCvalue)
            })

        names(d_res) <- datasets
        d_res <- do.call(cbind, d_res)
        return(d_res)
    })
    names(m_res) <- model_types
    m_res <- as.data.frame(do.call(rbind, m_res))
    rownames(m_res) <- model_types
    m_res$model_types <- model_types
    return(m_res)
})
names(res) = result_paths
res <- do.call(rbind, res)
res$model_types <- rownames(res) 
write.csv(res, "validation_results/oriImmuSig_mlr_model/validationset_performance_rf.csv", quote = F, row.names = F)
res[grep("metaAnalysis", res$model_types), ]



res <- lapply(as.list(result_paths), function(resPath){
# for(resPath in result_paths){
    m_res <- lapply(as.list(model_types), function(m){
        # for(m in model_types){
        d_res <- lapply(as.list(datasets), function(d){
            resPath2 <- paste0("validation_results/oriImmuSig_mlr_model/", resPath, "/", m, "/", d)
            files <- list.files(resPath2)
            files <- files[grep("_model_pred.Rdata", files)]
            # for(f in files){
            f = "ridge_model_pred.Rdata"
            load(paste0(resPath2, "/", f))
            modelroc <- pROC::roc( pred_ridge$data$truth, pred_ridge$data$prob.1 )
            AUCvalue = modelroc$auc[1]
            return(AUCvalue)
            })

        names(d_res) <- datasets
        d_res <- do.call(cbind, d_res)
        return(d_res)
    })
    names(m_res) <- model_types
    m_res <- as.data.frame(do.call(rbind, m_res))
    rownames(m_res) <- model_types
    m_res$model_types <- model_types
    return(m_res)
})
names(res) = result_paths
res <- do.call(rbind, res)
res$model_types <- rownames(res) 
write.csv(res, "validation_results/oriImmuSig_mlr_model/validationset_performance_ridge.csv", quote = F, row.names = F)





res <- lapply(as.list(result_paths), function(resPath){
# for(resPath in result_paths){
    m_res <- lapply(as.list(model_types), function(m){
        # for(m in model_types){
        d_res <- lapply(as.list(datasets), function(d){
            resPath2 <- paste0("validation_results/oriImmuSig_mlr_model/", resPath, "/", m, "/", d)
            files <- list.files(resPath2)
            files <- files[grep("_model_pred.Rdata", files)]
            # for(f in files){
            f = "svm_model_pred.Rdata"
            load(paste0(resPath2, "/", f))
            modelroc <- pROC::roc( pred_svm$data$truth, pred_svm$data$prob.1 )
            AUCvalue = modelroc$auc[1]
            return(AUCvalue)
            })

        names(d_res) <- datasets
        d_res <- do.call(cbind, d_res)
        return(d_res)
    })
    names(m_res) <- model_types
    m_res <- as.data.frame(do.call(rbind, m_res))
    rownames(m_res) <- model_types
    m_res$model_types <- model_types
    return(m_res)
})
names(res) = result_paths
res <- do.call(rbind, res)
res$model_types <- rownames(res) 
write.csv(res, "validation_results/oriImmuSig_mlr_model/validationset_performance_svm.csv", quote = F, row.names = F)



res <- lapply(as.list(result_paths), function(resPath){
# for(resPath in result_paths){
    m_res <- lapply(as.list(model_types), function(m){
        # for(m in model_types){
        d_res <- lapply(as.list(datasets), function(d){
            resPath2 <- paste0("validation_results/oriImmuSig_mlr_model/", resPath, "/", m, "/", d)
            files <- list.files(resPath2)
            files <- files[grep("_model_pred.Rdata", files)]
            # for(f in files){
            f = "elasticnet_model_pred.Rdata"
            load(paste0(resPath2, "/", f))
            modelroc <- pROC::roc( pred_en$data$truth, pred_en$data$prob.1 )
            AUCvalue = modelroc$auc[1]
            return(AUCvalue)
            })

        names(d_res) <- datasets
        d_res <- do.call(cbind, d_res)
        return(d_res)
    })
    names(m_res) <- model_types
    m_res <- as.data.frame(do.call(rbind, m_res))
    rownames(m_res) <- model_types
    m_res$model_types <- model_types
    return(m_res)
})
names(res) = result_paths
res <- do.call(rbind, res)
res$model_types <- rownames(res) 
write.csv(res, "validation_results/oriImmuSig_mlr_model/validationset_performance_en.csv", quote = F, row.names = F)

# ------------------------------------------------------------------------------

resfiles <- list.files("validation_results/oriImmuSig_mlr_model/")
resfiles <- resfiles[grep("validationset_performance_", resfiles)]

res_all <- lapply(as.list(resfiles), function(f){
  df = read.csv(paste0("validation_results/oriImmuSig_mlr_model/",f))
  df$files = gsub(".csv", "", gsub("validationset_performance_","", f))
  return(df)
})
res_all = do.call(rbind, res_all)

write.csv(res_all, "validation_results/oriImmuSig_mlr_model/validationset_performance_all.csv", quote = F, row.names = F)



# ------------------------------------------------------------------------------
datasets <- c("Liu_dataset_CTLA4Naive_outcome", "GSE96619_outcome", "GSE115821_pre_outcome", "08_Lauss_dataset_outcome")


res_othermodel <- lapply(as.list(datasets), function(f){
  df = read.csv(paste0("validation_results/",f,"/SKCM_Metastatic_otherModel_AUC.csv"))
  names(df) = c('immune_sig', f)
  return(df)
})
res = merge(merge(merge(res_othermodel[[1]],res_othermodel[[2]], by = 'immune_sig'),res_othermodel[[3]], by = 'immune_sig'),res_othermodel[[4]], by = 'immune_sig')
write.csv(res, "validation_results/oriImmuSig_mlr_model/validationset_othermodel_performance.csv", quote = F, row.names = F)

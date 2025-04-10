# rm(list=ls())
workpath <- "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/"
setwd(workpath)

options(stringsAsFactors = F)
library(glmnet)
library(mlr)
library(ggplot2)
library(ggsci)
library(colorspace)
library(patchwork)
library(tidyverse)
library(aplot)
library(reshape2)
# devtools::install_github("ricardo-bion/ggradar", dependencies=TRUE)
library(ggradar)

source("06_results_summaryandplotting/scripts/immunesig_selection_function.R")

dir.create("07_plotting_v2/02_ITSprediction_metaPlusFilter_removeprior_v1/")

result_dir <- "05_immuneSig_ICBresponse/results_r0.4_ITSssgsea_CombatRef/"
resultpath = paste0(result_dir, "finalmodel_17datasets_ssgseanorm/metaAnalysis_filter_removeprior_v1/")


load(paste0(result_dir, "ITS_afterCombat17datasets_ssgsea.Rdata"))
datasets = unique(ITSres_selected$dataset)

vlSetName <- c("08_Lauss_dataset_outcome", 
               "GSE115821_pre_outcome", 
               "GSE96619_outcome",
               "Liu_dataset_CTLA4Naive_outcome", 
               "GSE176307")


# ITS-based models -------------------------------------------------------------

models = c("glmnet_ridge", "elastic_net", "rf", "svm", "xgboost")
cv_auc = list()
vlSetName_auc = list()
vlSetName_auc_all = list()
vlSetName_roc_all = list()

for(i in seq(length(models))){
    m = models[i]
    folds = paste0('Fold', 1:5)

    if(m == 'elastic_net'){
            load(paste0(resultpath, "/en/elasticnet_model_pred.Rdata"))
            pred = pred_en$data
            names(pred)[grep('truth',names(pred))] = 'tag'
    }else if(m == 'rf'){
            load(paste0(resultpath, "/rf/ranger_model_pred.Rdata"))
            pred = pred_rf$data
            names(pred)[grep('truth',names(pred))] = 'tag'
    }else if(m == "svm"){
            load(paste0(resultpath, "/", m, "/", m ,"_model_pred.Rdata"))
            pred = pred_svm$data
            names(pred)[grep('truth',names(pred))] = 'tag'
    }else if(m == "xgboost"){
            load(paste0(resultpath, "/", m, "/", m ,"_model_pred.Rdata"))
            pred = pred_xgboost$data
            names(pred)[grep('truth',names(pred))] = 'tag'
    }else if(m == "glmnet_ridge"){
            load(paste0(resultpath, "/", m, "/", m ,"_model_pred.Rdata"))

    }

    res_pred <- pred
    res_pred$SampleID <- row.names(res_pred)
    res_pred <- merge(ITSres_selected[c('SampleID', 'label')], res_pred)

    res_pred$tag <- factor(res_pred$tag,levels = c(1,0), order=TRUE)

    modelroc <- pROC::roc(res_pred$tag, res_pred$prob.1)
    vlSetName_roc_all[[i]] <- modelroc

    vlSetName_auc_all[[i]] <- modelroc$auc[1]

    res_auc <- lapply(as.list(vlSetName), function(dset){
        # dset = "GSE126044"
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

    vlSetName_auc[[i]] <- do.call(rbind, res_auc)
    
}


names(vlSetName_roc_all) <- c("ridge", "elastic net", "random forest", "svm", "xgboost")

names(vlSetName_auc_all) <- c("ridge", "elastic net", "random forest", "svm", "xgboost")
vlSetName_auc_all <- as.data.frame(do.call(rbind, vlSetName_auc_all))
names(vlSetName_auc_all) <- 'AUC'
vlSetName_auc_all$model <- rownames(vlSetName_auc_all)

names(vlSetName_auc) <- c("ridge", "elastic net", "random forest", "svm", "xgboost")
vlSetName_auc <- do.call(rbind, vlSetName_auc)
vlSetName_auc$model <- do.call(c, lapply(strsplit(rownames(vlSetName_auc), split="\\."), function(x)x[1]))

write.csv(vlSetName_auc_all, "07_plotting_v2/02_ITSprediction_metaPlusFilter_removeprior_v1/vlSetName_auc_all.csv", quote = F, row.names = F)
write.csv(vlSetName_auc, "07_plotting_v2/02_ITSprediction_metaPlusFilter_removeprior_v1/vlSetName_auc.csv", quote = F, row.names = F)


# ------------------------------------------------------------------------------
# fig2a 
# ------------------------------------------------------------------------------
p_s2b = ggplot(vlSetName_auc, aes(factor(datasets), factor(model, levels = c("elastic net", "ridge", "random forest", "svm", "xgboost")),  fill = auc)) + 
        geom_tile(aes(fill = auc), colour = "white", lwd = 2, linetype = 1)+
        scale_fill_gradient2(low = "#228383",mid = "white",high = "#B33333", midpoint = 0.5) + 
        geom_text(aes(label=round(auc, 3)),col ="black",size = 4) +
        theme_bw() + 
        theme(axis.title.x=element_blank(), # remove title
                axis.ticks.x=element_blank(), # remove x axis
                axis.title.y=element_blank(), # remove y axis
                axis.text.x = element_text(angle = 90, hjust = 1, size = 10), 
                axis.text.y = element_text(hjust = 1, size = 10,)) 

ggsave("07_plotting_v2/02_ITSprediction_metaPlusFilter_removeprior_v1/fig_s2_vlSetName.pdf", p_s2b, width = 8, height = 6)


p1 = ggplot(data = vlSetName_auc_all, mapping = aes(x = factor(model, levels = rev(c("elastic net", "ridge", "random forest", "svm", "xgboost"))), y = AUC)) + 
        geom_bar(stat = 'identity', fill = '#b2b2d3', color = 'black') +
        geom_text(aes(label = round(AUC, 3)),col ="black",size = 4) +
        geom_hline(yintercept=0.5, lty=2) +
        coord_flip() + 
        theme_bw() + 
        theme(# axis.title.x=element_blank(), # remove title
            axis.ticks.x=element_blank(), # remove x axis
            axis.title.y=element_blank(), # remove y axis
            axis.text.x = element_text(angle = 90, hjust = 1, size = 14), 
            axis.text.y = element_text(hjust = 1, size = 14))  # remove text on y axis
# p1


# other models -------------------------------------------------------------
# ------------------------------------------------------------------------------
validationSet <- c("Liu_dataset_CTLA4Naive_outcome",
                   "GSE96619_outcome",
                   "GSE115821_pre_outcome",
                   "08_Lauss_dataset_outcome")

data_test_all <- lapply(as.list(vlSetName), function(dataset){
    if(dataset == 'GSE176307'){
        cancer = "BLCA"
        sample_type = "Primary"
        result_dir = paste0("05_immuneSig_ICBresponse/validation_results/", dataset,"/")
        data_test <- read.csv(paste0(result_dir, "/", cancer, "_", sample_type, "_immuneSig_merged.csv"))

        processedDataPath = "/picb/bigdata/project/FengFYM/Immunotherapy_prediction/data/"
        dataset = "GSE176307.Rdata"
        load(paste0(processedDataPath, "RNAseq/GSE176307/processedData/",dataset))
        
        names(response) = c('patient', 'label')
        response <- response[!is.na(response$label), ]
        data_test <- merge(response, data_test, by = 'patient')

    }else{
        cancer = "SKCM"
        sample_type = "Metastatic"
        result_dir = paste0("05_immuneSig_ICBresponse/validation_results/", dataset,"/")
        data_test <- read.csv(paste0(result_dir, "/", cancer, "_", sample_type, "_immuneSig_merged_mlr.csv"))

    }

    data_test_name <- data.frame(name = names(data_test), immune_sig = tolower(names(data_test)))
    data_test_name <- immunesig_name_unify(data_test_name)
    names(data_test) <- data_test_name$immune_sig
    data_test$tag <- 1
    data_test[data_test$label == "NR",]$tag = 0
    data_test$tag <- factor(data_test$tag,levels = c(1,0), order=TRUE)
    data_test[is.na(data_test)] = 0
    names(data_test) <- c("patientID", names(data_test)[-1])
    data_test$dataset <- dataset
    data_test$cancer <- cancer
    return(data_test)
})
names(data_test_all) <- vlSetName

otherModel = c('tide_tide', 'cd8_tide', 'ifng_tide', 'tip_signature_tip','ips_ips', 'merck18_tide', 'oe_res')
data_test_selected = lapply(data_test_all, function(x) x[c("patientID", "label","tag", "dataset","cancer", otherModel)])
data_test_selected <- do.call(rbind, data_test_selected)

rownames(data_test_selected) <- paste0(data_test_selected$dataset, data_test_selected$patient)
res_df <-data_test_selected[!is.na(data_test_selected$label), ]
res_df = res_df[order(res_df$label), ]
res_df = res_df[order(res_df$label), ]
res_df$label <- factor(res_df$label,levels = c("NR","R"), order=TRUE)
res_df$tag <- 1
res_df[res_df$label == "NR",]$tag = 0
res_df$tag <- factor(res_df$tag,levels = c(1,0), order=TRUE)


# ------------------------------------------------------------------------------
# all datasets
# ------------------------------------------------------------------------------
AUC_VAL <- apply(res_df[-c(1:5)], 2, function(x){
    df <- data.frame(x = x, tag= res_df$tag)
    if(length(which(is.na(df$x)) > 0 ))df[is.na(df$x),]$x = 0
    modelroc <- pROC::roc( df$tag, df$x)
    AUCvalue = modelroc$auc[1]
    # print(AUCvalue)
    return(AUCvalue)
})
AUC_VAL <- as.data.frame(AUC_VAL)
AUC_VAL$model <- rownames(AUC_VAL)
names(AUC_VAL) <- c('AUC', 'model')


# ------------------------------------------------------------------------------
othermodels = c('tide_tide', 'cd8_tide', 'ifng_tide', 'tip_signature_tip','ips_ips', 'merck18_tide', 'oe_res')
ROC_VAL <- lapply(as.list(othermodels), function(m){

    x= res_df[,m ]
    df <- data.frame(x = x, tag= res_df$tag)
    if(m == 'tide_tide'){ df$x = df$x*(-1)    }
    if(length(which(is.na(df$x)) > 0 ))df[is.na(df$x),]$x = 0
    df$tag <- factor(df$tag,levels = c(1,0), order=TRUE)
    modelroc <- pROC::roc( df$tag, df$x)
    return(modelroc)
})
names(ROC_VAL) = othermodels

# ------------------------------------------------------------------------------
# single dataset
# ------------------------------------------------------------------------------

VALname = names(data_test_all)
AnnotColumn <- c("patientID","label","tag","dataset", "cancer")
singledataset_AUC <- list()
for(i in 1:length(VALname)){
    tmp <- data_test_all[[VALname[i]]]
    AUC_df <- apply(tmp[-which(is.element(names(tmp), AnnotColumn))], 2, function(x){
        df <- data.frame(x = x, tag= tmp$tag)
        if(length(which(is.na(df$x)) >0 ))df[is.na(df$x),]$x = 0
        modelroc <- pROC::roc( df$tag, df$x)
        AUCvalue = modelroc$auc[1]
        return(AUCvalue)
    })
    singledataset_AUC[[i]] = AUC_df
}

names(singledataset_AUC) <- VALname
singledataset_AUC_selected = lapply(singledataset_AUC, function(x) x[otherModel])
singledataset_AUC_selected <- as.data.frame(t(do.call(rbind, singledataset_AUC_selected)))
singledataset_AUC_selected$model = row.names(singledataset_AUC_selected)
singledataset_AUC_selected <- melt(singledataset_AUC_selected)
singledataset_AUC_selected <- singledataset_AUC_selected[c(2,3,1)]
names(singledataset_AUC_selected) <- names(vlSetName_auc)
write.csv(singledataset_AUC_selected, "07_plotting_v2/02_ITSprediction_metaPlusFilter_removeprior_v1/vlSetName_auc_otherModel.csv", quote = F)


res = rbind(vlSetName_auc[vlSetName_auc$model == 'elastic net',],
            singledataset_AUC_selected)
res$model = factor(res$model, levels = rev(c('elastic net', 'cd8_tide', 'ifng_tide', 'tide_tide', 'merck18_tide', 'tip_signature_tip','ips_ips', 'oe_res')))

AUC_VAL_selected = rbind(vlSetName_auc_all[is.element(vlSetName_auc_all$model, "elastic net"), ], AUC_VAL)
AUC_VAL_selected$model = factor(AUC_VAL_selected$model, levels = rev(c('elastic net', 'cd8_tide', 'ifng_tide', 'tide_tide', 'merck18_tide', 'tip_signature_tip','ips_ips', 'oe_res')))

p2_1 = ggplot(data = AUC_VAL_selected, mapping = aes(x = factor(model), y = AUC)) + 
        geom_bar(stat = 'identity', fill = '#b2b2d3', color = 'black') +
        geom_text(aes(label = round(AUC, 3)),col ="black",size = 4) +
        geom_hline(yintercept=0.5, lty=2) +
        coord_flip() + 
        theme_bw() + 
        theme(# axis.title.x=element_blank(), # remove title
              axis.ticks.x=element_blank(), # remove x axis
              axis.title.y=element_blank(), # remove y axis
              axis.text.x = element_text(angle = 90, hjust = 1, size = 14),
              axis.text.y = element_blank()) 

p2_2 = ggplot(res, aes(factor(datasets), factor(model),  fill = auc)) + 
        geom_tile(aes(fill = auc), colour = "white", lwd = 2, linetype = 1)+
        scale_fill_gradient2(low = "#228383",mid = "white",high = "#B33333", midpoint = 0.5) + 
        geom_text(aes(label=round(auc, 3)),col ="black",size = 4) +
        theme_bw() + 
        theme(axis.title.x=element_blank(), # remove title
                axis.ticks.x=element_blank(), # remove x axis
                axis.title.y=element_blank(), # remove y axis
                axis.text.x = element_text(angle = 90, hjust = 1, size = 14), 
                axis.text.y = element_text(hjust = 1, size = 14,)) + 
        guides(fill='none')


res_mat = spread(res, model, auc)
res_mat = res_mat[c('datasets', 'elastic net', 'cd8_tide', 'ifng_tide', 'tide_tide', 'merck18_tide', 'tip_signature_tip','ips_ips', 'oe_res')]
p3 = ggradar(res_mat, legend.position = 'bottom')

p_fig2cd1 = (p1 + p3) + plot_layout(widths = c(1, 3))
# p_fig2cd2 = (p2_2 %>% insert_right(p2_1, width = .1))
p_fig2cd2 = (p2_2 | p2_1 )+ plot_layout(widths = c(3, 1))

p_fig2cd = (p_fig2cd1 / p_fig2cd2 ) + plot_layout(widths = c(2, 1))

ggsave("07_plotting_v2/02_ITSprediction_metaPlusFilter_removeprior_v1/fig_2cd_vlSetName.pdf", p_fig2cd, width = 15, height = 18)
 
library(pROC)
pROC::roc.test(vlSetName_roc_all[[2]], ROC_VAL[[1]], method = 'venkatraman')
pROC::roc.test(vlSetName_roc_all[[2]], ROC_VAL[[2]], method = 'venkatraman')
pROC::roc.test(vlSetName_roc_all[[2]], ROC_VAL[[3]], method = 'venkatraman')
pROC::roc.test(vlSetName_roc_all[[2]], ROC_VAL[[4]], method = 'venkatraman')
pROC::roc.test(vlSetName_roc_all[[2]], ROC_VAL[[5]], method = 'venkatraman')
pROC::roc.test(vlSetName_roc_all[[2]], ROC_VAL[[6]], method = 'venkatraman')
pROC::roc.test(vlSetName_roc_all[[2]], ROC_VAL[[7]], method = 'venkatraman')

# pROC::roc.test(vlSetName_roc_all[[2]], ROC_VAL[[1]], method = 'delong')
# pROC::roc.test(vlSetName_roc_all[[2]], ROC_VAL[[2]], method = 'delong')
# pROC::roc.test(vlSetName_roc_all[[2]], ROC_VAL[[3]], method = 'delong')
# pROC::roc.test(vlSetName_roc_all[[2]], ROC_VAL[[4]], method = 'delong')
# pROC::roc.test(vlSetName_roc_all[[2]], ROC_VAL[[5]], method = 'delong')
# pROC::roc.test(vlSetName_roc_all[[2]], ROC_VAL[[6]], method = 'delong')
# pROC::roc.test(vlSetName_roc_all[[2]], ROC_VAL[[7]], method = 'delong')

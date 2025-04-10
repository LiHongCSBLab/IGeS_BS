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

result_dir <- "05_immuneSig_ICBresponse/results_r0.4_ITSssgsea_CombatRef/"
resultpath = paste0(result_dir, "modelCV_17datasets_ssgseanorm/metaAnalysis_removeprior_v3/")

load(paste0(result_dir, "ITS_afterCombat17datasets_ssgsea.Rdata"))
datasets = unique(ITSres_selected$dataset)
vlSetName <- c("08_Lauss_dataset_outcome", "GSE115821_pre_outcome", "GSE96619_outcome","Liu_dataset_CTLA4Naive_outcome", "GSE176307")
dcSetName <- datasets[!is.element(datasets, vlSetName)]


# ITS-based models -------------------------------------------------------------

models = c("glmnet_ridge", "elastic_net", "rf", "svm", "xgboost")
cv_auc = list()
discoverset_auc = list()
discoverset_auc_all = list()
discoverset_roc_all = list()


for(i in seq(length(models))){
    # i= 2
    m = models[i]
    folds = paste0('Fold', 1:5)

    if(m == 'elastic_net'){
        res_pred <- lapply(as.list(folds), function(f){
            # f = 'Fold1'
            load(paste0(resultpath, "/en/elasticnet_model_pred_",f,".Rdata"))
            pred = pred_en$data
            names(pred)[grep('truth',names(pred))] = 'tag'
            return(pred)
        })
    }else if(m == 'rf'){
        res_pred <- lapply(as.list(folds), function(f){
            # f = 'Fold1'
            load(paste0(resultpath, "/rf/ranger_model_pred_",f,".Rdata"))
            pred = pred_rf$data
            names(pred)[grep('truth',names(pred))] = 'tag'
            return(pred)
        })
    }else if(m == "svm"){
        res_pred <- lapply(as.list(folds), function(f){
            # f = 'Fold1'
            load(paste0(resultpath, "/", m, "/", m ,"_model_pred_",f,".Rdata"))
            pred = pred_svm$data
            names(pred)[grep('truth',names(pred))] = 'tag'
            return(pred)
        })
    }else if(m == "xgboost"){
        res_pred <- lapply(as.list(folds), function(f){
            # f = 'Fold1'
            load(paste0(resultpath, "/", m, "/", m ,"_model_pred_",f,".Rdata"))
            pred = pred_xgboost$data
            names(pred)[grep('truth',names(pred))] = 'tag'
            return(pred)
        })
    }else if(m == "glmnet_ridge"){
        res_pred <- lapply(as.list(folds), function(f){
            # f = 'Fold1'
            load(paste0(resultpath, "/", m, "/", m ,"_model_pred_",f,".Rdata"))
            return(pred)
        })
    }

    cv_auc[[i]] <- data.frame(folds=folds, 
                              auc = do.call(c,lapply(res_pred, function(x){pROC::roc(x$tag, x$prob.1)$auc[1]})))


    res_pred <- do.call(rbind, res_pred)
    res_pred$SampleID <- row.names(res_pred)
    res_pred <- merge(ITSres_selected[c('SampleID', 'dataset','label')], res_pred)
    print(dim(res_pred))
    # pred2 <- ROCR::prediction(res_pred$prob.1 , res_pred$label)
    # aucPerf <- ROCR::performance( pred2, 'auc' )
    # AUCValue<-aucPerf@y.values[[1]]
    # AUCValue
    # summary(res_pred$prob.1)
    # res_pred$predtag = 0
    # res_pred[res_pred$prob.1 > 1.286,]$predtag = 1
    # table(res_pred[c('tag', 'predtag')])
    res_pred$tag <- factor(res_pred$tag, levels = c(1,0), order=TRUE)
    # boxplot(res_pred$prob.1 ~ res_pred$tag)
    # boxplot(res_pred$prob.0 ~ res_pred$tag)
    
    modelroc <- pROC::roc(res_pred$tag, res_pred$prob.1)
    discoverset_roc_all[[i]] = modelroc

    # plot(modelroc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),
    #     grid.col=c("green", "red"), max.auc.polygon=TRUE,
    #     auc.polygon.col="lightblue", print.thres=TRUE)
    discoverset_auc_all[[i]] <- modelroc$auc[1]

    res_auc <- lapply(as.list(dcSetName), function(dset){
        # dset = "Zhao_dataset"
        print(dset)
        tmp = res_pred[grep(dset, res_pred$SampleID),]
        pred2 <- ROCR::prediction(tmp$prob.1 , tmp$label)
        aucPerf <- ROCR::performance( pred2, 'auc' )
        AUCValue<-aucPerf@y.values[[1]]
        res = data.frame(datasets = dset, auc = AUCValue)
        return(res)

        # modelroc <- pROC::roc(tmp$tag, tmp$prob.1, levels = c(0,1), direction = '<')
        # # plot(modelroc)
        # res = data.frame(datasets = dset, auc = modelroc$auc[1])
        # return(res)

    })

    discoverset_auc[[i]] <- do.call(rbind, res_auc)
    
}

# cv_auc
names(cv_auc) <- c("ridge", "elastic net", "random forest", "svm", "xgboost")
cv_auc <- do.call(rbind, cv_auc)
cv_auc$model <- do.call(c, lapply(strsplit(rownames(cv_auc), split="\\."), function(x)x[1]))


names(discoverset_auc_all) <- c("ridge", "elastic net", "random forest", "svm", "xgboost")
discoverset_auc_all <- as.data.frame(do.call(rbind, discoverset_auc_all))
names(discoverset_auc_all) <- 'AUC'
discoverset_auc_all$model <- rownames(discoverset_auc_all)

names(discoverset_roc_all) <- c("ridge", "elastic net", "random forest", "svm", "xgboost")


names(discoverset_auc) <- c("ridge", "elastic net", "random forest", "svm", "xgboost")
discoverset_auc <- do.call(rbind, discoverset_auc)
discoverset_auc$model <- do.call(c, lapply(strsplit(rownames(discoverset_auc), split="\\."), function(x)x[1]))

write.csv(cv_auc, "07_plotting_v2/02_ITSprediction_removeprior_v3/cv_auc.csv", quote = F, row.names = F)
write.csv(discoverset_auc_all, "07_plotting_v2/02_ITSprediction_removeprior_v3/discoverset_auc_all.csv", quote = F, row.names = F)
write.csv(discoverset_auc, "07_plotting_v2/02_ITSprediction_removeprior_v3/discoverset_auc.csv", quote = F, row.names = F)


# ------------------------------------------------------------------------------
# fig2a 
# ------------------------------------------------------------------------------


p1 = ggplot(data = discoverset_auc_all, mapping = aes(x = factor(model, levels = rev(c("elastic net", "ridge", "random forest", "svm", "xgboost"))), y = AUC)) + 
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

p_s2a = ggplot(data = cv_auc, mapping = aes(x = factor(model, levels = c("elastic net", "ridge", "random forest", "svm", "xgboost")), y = auc, fill = folds)) + 
        geom_bar(stat = 'identity', position = 'dodge', color = 'black', alpha = 0.5) +
        geom_text(mapping = aes(x = factor(model), y = auc, label = round(auc, 3)),col ="black", size = 4, position = position_dodge(0.9)) +
        coord_flip() + 
        theme_bw() + 
        theme(axis.title.x=element_blank(), # remove title
            axis.ticks.x=element_blank(), # remove x axis
            axis.title.y=element_blank(), # remove y axis
            axis.text.x = element_text(angle = 90, hjust = 1, size = 10), 
            axis.text.y = element_text(hjust = 1, size = 10))  # remove text on y axis


p_s2b = ggplot(discoverset_auc, aes(factor(datasets), factor(model, levels = c("elastic net", "ridge", "random forest", "svm", "xgboost")),  fill = auc)) + 
        geom_tile(aes(fill = auc), colour = "white", lwd = 2, linetype = 1)+
        scale_fill_gradient2(low = "#228383",mid = "white",high = "#B33333", midpoint = 0.5) + 
        geom_text(aes(label=round(auc, 3)),col ="black",size = 4) +
        theme_bw() + 
        theme(axis.title.x=element_blank(), # remove title
                axis.ticks.x=element_blank(), # remove x axis
                axis.title.y=element_blank(), # remove y axis
                axis.text.x = element_text(angle = 90, hjust = 1, size = 10), 
                axis.text.y = element_text(hjust = 1, size = 10,)) 
p_s2 = (p_s2a + p_s2b )+ 
  plot_layout(widths = c(1, 2))

ggsave("07_plotting_v2/02_ITSprediction_removeprior_v3/fig_s2_discoverySet.pdf", p_s2, width = 15, height = 8)


# other models -------------------------------------------------------------
# ------------------------------------------------------------------------------
load("05_immuneSig_ICBresponse/results/oriImmuneSig_res_all.Rdata")
merged_sig1 <- lapply(as.list(names(merged_sig)), function(x){
  df = merged_sig[[x]]
  print(x)
  df <- df[!is.na(df$label), ]
  df$tag <- 1
  df[df$label == "NR",]$tag = 0
  df$tag <- factor(df$tag,levels = c(1,0), order=TRUE)
  return(df)
})

names(merged_sig1) <- names(merged_sig)
merged_sig1 <- do.call(rbind,merged_sig1)

rownames(merged_sig1) <- paste0(merged_sig1$dataset, merged_sig1$patient)
res_df <-merged_sig1[!is.na(merged_sig1$label), ]
res_df = res_df[order(res_df$label), ]
# res_df_otherModel <- res_df[c(grep("", names(res_df)),
#                               grep("ips_ips", names(res_df)),
#                               grep("tide_tide", names(res_df)))]
# res_df


res_df = res_df[order(res_df$label), ]
res_df$label <- factor(res_df$label,levels = c("NR","R"), order=TRUE)
res_df$tag <- 1
res_df[res_df$label == "NR",]$tag = 0
res_df$tag <- factor(res_df$tag,levels = c(1,0), order=TRUE)

# ------------------------------------------------------------------------------
# all datasets
# ------------------------------------------------------------------------------

AUC_allDS <- apply(res_df[-c(1:5,397)], 2, function(x){
    df <- data.frame(x = x, tag= res_df$tag)
    if(length(which(is.na(df$x)) > 0 ))df[is.na(df$x),]$x = 0
    modelroc <- pROC::roc( df$tag, df$x)
    AUCvalue = modelroc$auc[1]
    # print(AUCvalue)
    return(AUCvalue)
})
AUC_allDS <- as.data.frame(AUC_allDS)
AUC_allDS$model <- rownames(AUC_allDS)
names(AUC_allDS) <- c('AUC', 'model')
AUC_allDS_selected <- AUC_allDS[is.element(AUC_allDS$model, c('tide_tide', 'cd8_tide', 'ifng_tide', 'tip_signature_tip','ips_ips', 'merck18_tide', 'oe_res')), ]


# ------------------------------------------------------------------------------
othermodels = c('tide_tide', 'cd8_tide', 'ifng_tide', 'tip_signature_tip','ips_ips', 'merck18_tide', 'oe_res')
ROC_allDS <- lapply(as.list(othermodels), function(m){
    x= res_df[,m ]
    df <- data.frame(x = x, tag= res_df$tag)
    if(length(which(is.na(df$x)) > 0 ))df[is.na(df$x),]$x = 0
    modelroc <- pROC::roc( df$tag, df$x)
    return(modelroc)
})
names(ROC_allDS) = othermodels

# ------------------------------------------------------------------------------
# single dataset
# ------------------------------------------------------------------------------

dsname = unique(res_df$dataset)
singledataset_AUC <- list()
for(i in 1:length(dsname)){
  tmp <- res_df[res_df$dataset == dsname[i], ]
  AUC_df <- apply(tmp[-c(1:5,397)], 2, function(x){
    df <- data.frame(x = x, tag= tmp$tag)
    if(length(which(is.na(df$x)) > 0 ))df[is.na(df$x),]$x = 0
    #   pred <- prediction(df$x, df$tag)
    #   aucPerf <- ROCR::performance( pred, 'auc' )
    #   AUCValue<-aucPerf@y.values[[1]]
    modelroc <- pROC::roc( df$tag, df$x)
    AUCvalue = modelroc$auc[1]
    # print(AUCvalue)
    return(AUCvalue)
     
  })
  singledataset_AUC[[i]] = AUC_df
}
names(singledataset_AUC) <- dsname

singledataset_AUC <- as.data.frame(t(do.call(rbind, singledataset_AUC)))

# names(alldataset_AUCValue) <- "alldataset_AUCValue"
# alldataset_AUCValue <- cbind(alldataset_AUCValue, singledataset_AUC)

write.csv(singledataset_AUC, "07_plotting_v2/02_ITSprediction_removeprior_v3/discoverySet_otherModel_all.csv", quote = F)

singledataset_AUC_selected <- singledataset_AUC[is.element(rownames(singledataset_AUC), 
                                                         c('tide_tide', 'cd8_tide', 'ifng_tide', 'tip_signature_tip','ips_ips', 'merck18_tide', 'oe_res')),]

singledataset_AUC_selected$model = row.names(singledataset_AUC_selected)
singledataset_AUC_selected <- melt(singledataset_AUC_selected)
singledataset_AUC_selected <- singledataset_AUC_selected[c(2,3,1)]
names(singledataset_AUC_selected) <- names(discoverset_auc)
write.csv(singledataset_AUC_selected, "07_plotting_v2/02_ITSprediction_removeprior_v3/discoverset_auc_otherModel.csv", quote = F)


res = rbind(discoverset_auc[discoverset_auc$model == 'elastic net',],
            singledataset_AUC_selected)
res$model = factor(res$model, levels = rev(c('elastic net', 'cd8_tide', 'ifng_tide', 'tide_tide', 'merck18_tide', 'tip_signature_tip','ips_ips', 'oe_res')))

AUC_allDS_selected = rbind(discoverset_auc_all[is.element(discoverset_auc_all$model, "elastic net"), ], AUC_allDS_selected)
AUC_allDS_selected$model = factor(AUC_allDS_selected$model, levels = rev(c('elastic net', 'cd8_tide', 'ifng_tide', 'tide_tide', 'merck18_tide', 'tip_signature_tip','ips_ips', 'oe_res')))

p2_1 = ggplot(data = AUC_allDS_selected, mapping = aes(x = factor(model), y = AUC)) + 
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

p_fig2ab1 = (p1 + p3) + plot_layout(widths = c(1, 3))
# p_fig2ab2 = (p2_2 %>% insert_right(p2_1, width = .1))
p_fig2ab2 = (p2_2 | p2_1 )+ plot_layout(widths = c(3, 1))

p_fig2ab = (p_fig2ab1 / p_fig2ab2 ) + plot_layout(widths = c(2, 1))
ggsave("07_plotting_v2/02_ITSprediction_removeprior_v3/fig_2ab_discoverySet.pdf", p_fig2ab, width = 15, height = 18)



# ROC test ---------------------------------------------------------------------
discoverset_roc_all
ROC_allDS

pROC::roc.test(discoverset_roc_all[[2]], ROC_allDS[[1]], method = 'venkatraman')
pROC::roc.test(discoverset_roc_all[[2]], ROC_allDS[[2]], method = 'venkatraman')
pROC::roc.test(discoverset_roc_all[[2]], ROC_allDS[[3]], method = 'venkatraman')
pROC::roc.test(discoverset_roc_all[[2]], ROC_allDS[[4]], method = 'venkatraman')
pROC::roc.test(discoverset_roc_all[[2]], ROC_allDS[[5]], method = 'venkatraman')
pROC::roc.test(discoverset_roc_all[[2]], ROC_allDS[[6]], method = 'venkatraman')
pROC::roc.test(discoverset_roc_all[[2]], ROC_allDS[[7]], method = 'venkatraman')


pROC::roc.test(discoverset_roc_all[[2]], ROC_allDS[[1]], method = 'delong')
pROC::roc.test(discoverset_roc_all[[2]], ROC_allDS[[2]], method = 'delong')
pROC::roc.test(discoverset_roc_all[[2]], ROC_allDS[[3]], method = 'delong')
pROC::roc.test(discoverset_roc_all[[2]], ROC_allDS[[4]], method = 'delong')
pROC::roc.test(discoverset_roc_all[[2]], ROC_allDS[[5]], method = 'delong')
pROC::roc.test(discoverset_roc_all[[2]], ROC_allDS[[6]], method = 'delong')
pROC::roc.test(discoverset_roc_all[[2]], ROC_allDS[[7]], method = 'delong')
# load data --------------------------------------------------------------------
# rm(list=ls())
# ------------------------------------------------------------------------------
# source("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/06_results_summaryandplotting/scripts/immunesig_selection_function.R")
workpath <- "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/"
setwd(workpath)
options(stringsAsFactors = F)

# Checking and loading R packages -----------------------------------------
library(dplyr)
library(ggplot2)
library(pheatmap)
library(ROCR)
library(reshape2)
library(stringr)
library(sva)

library(FactoMineR)
library(ggplot2)
library(patchwork)
# Loading relied functions -----------------------------------------------------
# TIGS, AUC and Wilcox test
source("05_immuneSig_ICBresponse/scripts/functions/05_immunotherapy_immuneSig_analysis_function.R")
source("05_immuneSig_ICBresponse/scripts/functions/predictor_ITS.R")
source("06_results_summaryandplotting/scripts/immunesig_selection_function.R")

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# batch correction with combat  ------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

result_dir <- "05_immuneSig_ICBresponse/results_mouse_ITSssgsea_CombatRef/"
dir.create(result_dir)

# ------------------------------------------------------------------------------
load(file = "05_immuneSig_ICBresponse/results_r0.4_ITSssgsea_Combat/ICBdataset_Exp17datasets.Rdata")


for(i in 1:length(data)){
  if(i == 1){
    cm = rownames(data[[1]])
  } else{
    cm = intersect(cm,rownames(data[[i]]))
  }
}

data_bc_all <- lapply(data,function(z){ return(t(z[cm,]))})

data_bc <- do.call(rbind,data_bc_all[c(1:12)])
dataAnnot_all <- do.call(rbind, dataAnnot)
dataAnnot_all$dataset <- str_split_fixed(rownames(dataAnnot_all),'\\.',2)[,1]
data_bc <- as.data.frame(data_bc)
data_bc$SampleID <- row.names(data_bc)
data_bc <- merge(data_bc, dataAnnot_all, by = 'SampleID')
data_bc$batch <- data_bc$dataset

rownames(data_bc) <- data_bc$SampleID
edata_log <- t(data_bc[-which(is.element(names(data_bc), c("SampleID", "patientID","label","cancer","sample_type","dataset","batch")))])

data_bc2 <- data_bc[!is.na(data_bc$label),]
edata_log2 <- edata_log[,data_bc2$SampleID]

data_train <- ComBat(dat = edata_log2, 
                 batch = data_bc2$batch)
data_train <- as.data.frame(t(data_train))


mod = model.matrix(~as.factor(label), data=data_bc2)
data_train2 <- ComBat(dat = edata_log2, 
                      # mod = mod,
                      batch = data_bc2$batch)
data_train2 <- as.data.frame(t(data_train2))

rownames(dataAnnot_all) <- dataAnnot_all$SampleID



data_mouse <- read.csv("05_immuneSig_ICBresponse/data/mouse_zhao/noWT_addNewSet/TPMmat_TPM.txt", sep = '\t', header = T, row.names = 1) 

data_train2 = t(data_train2)
edata_s <- data_mouse[intersect(rownames(data_mouse), rownames(data_train2)), ]
data_train2 = data_train2[intersect(rownames(data_mouse), rownames(data_train2)), ]

mydata = cbind(data_train2, log2(1+edata_s))

dataBatch = c(rep('train', ncol(data_train2)), rep('test', ncol(edata_s)))
data_val2 <- ComBat(dat = mydata, 
                    ref.batch = 'train', 
                    batch = dataBatch)

                
data_val2 <- as.data.frame(t(data_val2))
pca_fmr = PCA(data_val2, 
                scale.unit = T, 
                ncp = 6, graph = F)  ## !! just one step !! ##

## generate PCA using ggplot2
PCs.2dim = pca_fmr$ind$coord[, 1:2]    # 读取前两维的PC的二维轴平面的坐标
colnames(PCs.2dim) = c("Dim1", "Dim2")
xlab = paste("Dim1 (", round(pca_fmr$eig[1,2], 2), "%)", sep = ""); xlab   # 设置x轴标题
ylab = paste("Dim2 (", round(pca_fmr$eig[2,2], 2), "%)", sep = ""); ylab   # 设置y轴标题

## plot results of PCA with ggplot2
PCs.2dim = as.data.frame(PCs.2dim)
PCs.2dim$SampleID = rownames(PCs.2dim)

data_val2$MLbatch = dataBatch
data_val2$SampleID = rownames(data_val2)
Fig1.PCs = merge(data_val2[c("SampleID", "MLbatch")], PCs.2dim, by = "SampleID") 

p_MLbatch <- 
    ggplot(Fig1.PCs, aes(x=Dim1, y=Dim2, color= MLbatch )) + 
    geom_point(size = 1) +
    theme_bw() +
    xlab(xlab) + ylab(ylab) + 
    scale_alpha_discrete(range = c(0.5,1)) +
    coord_fixed() +  
    ggtitle("MLbatch")
ggsave(paste0(result_dir, "Exp_valds_mouse_afterCombat_pca.pdf"), p_MLbatch, width = 10, height = 5)

data_ac = data_val2
save(data_ac, file = paste0(result_dir, "ExpMouse_afterCombatRef.Rdata"))

# ----------------------------------------------------------
# ----------------------------------------------------------
# ITSssgsea  -------------------------------------------------------------------
# ----------------------------------------------------------
# ----------------------------------------------------------
source("05_immuneSig_ICBresponse/scripts/functions/05_immunotherapy_immuneSig_analysis_function.R")
source("05_immuneSig_ICBresponse/scripts/functions/predictor_ITS.R")
source("06_results_summaryandplotting/scripts/immunesig_selection_function.R")

load(paste0(result_dir, "ExpMouse_afterCombatRef.Rdata"))

savepath_allresult = "02_tumorGene_immuneSig_correlation/results_selected/ITS_oe_original_spearman_TUMERIC_r0.4/"
# after combat: data_ac
savepath = paste0(result_dir, "ITSssgsea_afterCombat17mouse/")
dir.create(savepath)
dataAnnot_all <- dataAnnot_all[!is.na(dataAnnot_all$label), ]
datasets = unique(dataAnnot_all$dataset)
df <- t(data_ac)

cancer = 'LIHC'
sample_type = 'Primary'
load(paste0(savepath_allresult, "for_GSVA/pcor/pcor_spearman_", cancer, "_",sample_type, "_positive_200.Rdata"))
load(paste0(savepath_allresult, "for_GSVA/pcor/pcor_spearman_", cancer, "_",sample_type, "_negative_200.Rdata"))
genelist_p_name = data.frame(immunesig = names(genelist_p), immune_sig = tolower(names(genelist_p)))
genelist_p_name = immunesig_name_unify(genelist_p_name)
names(genelist_p) = genelist_p_name$immune_sig

genelist_n_name = data.frame(immunesig = names(genelist_n), immune_sig = tolower(names(genelist_n)))
genelist_n_name = immunesig_name_unify(genelist_n_name)
names(genelist_n) = genelist_n_name$immune_sig

enrich_method = "ssgsea"
ITSp <- immune_sigAnalysis(data = df, 
                            gs = genelist_p,
                            method = enrich_method, 
                            kcdf = "Gaussian",
                            ssgsea.norm = TRUE) 

ITSp_res <- as.data.frame(t(ITSp))
ITSp_res$SampleID <- rownames(ITSp_res)
ds = data_val2[c("SampleID", "MLbatch")]
ITSp_res <- inner_join(ds, ITSp_res, by = "SampleID")
write.csv(ITSp_res, paste0(savepath, "/mouse_ITSp_", enrich_method, "_ssgseanorm.csv"), 
            row.names =F, quote = F)


ITSn <- immune_sigAnalysis(data = df, 
                            gs = genelist_n,
                            method = enrich_method, 
                            kcdf = "Gaussian",
                            ssgsea.norm = TRUE) 
ITSn_res <- as.data.frame(t(ITSn))
ITSn_res$SampleID <- rownames(ITSn_res)
ITSn_res <- inner_join(ds, ITSn_res,by='SampleID')
write.csv(ITSn_res, paste0(savepath, "/mouseITSn_", enrich_method, "_ssgseanorm.csv"), 
            row.names =F, quote = F)


# prediction
library(mlr)
enrich_method = "ssgsea"
ITSp_res <- read.csv(paste0(savepath, "/mouse_ITSp_", enrich_method, "_ssgseanorm.csv"), row.names = 1)

mouseMetaInfo <- read.csv("05_immuneSig_ICBresponse/data/mouse_zhao/noWT_addNewSet/mMetadataMerged.csv")
rownames(mouseMetaInfo) = paste0(mouseMetaInfo$Genotype, "_", mouseMetaInfo$SampleID)
mouseMetaInfo$Genotype_SampleID = paste0(mouseMetaInfo$Genotype, "_", mouseMetaInfo$SampleID)
data_test = ITSp_res[ITSp_res$MLbatch == 'test', ]
data_test$SampleID = row.names(data_test)
data_test <- merge(mouseMetaInfo, data_test, by = 'SampleID')
row.names(data_test) =  paste0(data_test$Genotype, "_", data_test$SampleID)
data_test <- data_test[-is.element(names(data_test), c("X",
                                                      "Genotype",
                                                      "batch",
                                                      "subtype",
                                                      "design",
                                                      "Group"))]

dir.create(paste0(savepath, "/metaAnalysis_filter/"))
load("05_immuneSig_ICBresponse/results_r0.4_ITSssgsea_CombatRef/finalmodel_17datasets_ssgseanorm/metaAnalysis_filter/en/elasticnet_model_pred.Rdata")
sig_select <- read.csv("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/07_plotting_v2/01_metaAnalysis/table/metaSelected_immunesigAnnot_subtypedetail_filter.txt", sep = '\t')
rownames(sig_select) = sig_select$immune_sig
features_selected <- sig_select$immune_sig

data_test_selected <- data_test[c(intersect(features_selected, names(data_test)))]
pdf(paste0(savepath, "/metaAnalysis_filter/mouseITS_", enrich_method, "_IGeS.pdf"), 
    width = 15, height = 15)
pheatmap(data_test_selected, scale = 'column', 
         cluster_col = F,
         annotation_row = mouseMetaInfo["Genotype"],
         annotation_col = sig_select[c("flag","Type")])
dev.off()

data_test_selected$tag = "0"
mouse_testtask = makeClassifTask(id = "mouse", 
                                 data = data_test_selected,
                                 target = "tag")
pred_en = predict(m_en, task = mouse_testtask)
pred_mouse = pred_en$data
pred_mouse$Genotype_SampleID = rownames(pred_mouse)
pred_mouse <- merge(mouseMetaInfo, pred_mouse, by = 'Genotype_SampleID' )
pred_mouse = pred_mouse[order(pred_mouse$prob.1, decreasing = T), ]
gt_score = aggregate(pred_mouse['prob.1'], by = list(pred_mouse$Genotype), mean)
pred_mouse$Genotype = factor(pred_mouse$Genotype, levels = gt_score[order(gt_score$prob.1, decreasing = T), ]$Group.1)
p = ggplot(pred_mouse, aes(x=Genotype, y=prob.1)) +
      geom_boxplot() + 
      # geom_jitter(width = 0.1,alpha = 0.2)+
      geom_point(size=1) +
      theme_bw()+   
      theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(paste0(savepath, "/metaAnalysis_filter/mouseITS_", enrich_method, "_pred.pdf"), p, width = 10, height = 5)
write.csv(pred_mouse, paste0(savepath, "/metaAnalysis_filter/mouseITS_", enrich_method, "_enPred.csv"),
             row.names =F, quote = F)


dir.create(paste0(savepath, "/metaAnalysis_filter_removeprior_v2/"))
load("05_immuneSig_ICBresponse/results_r0.4_ITSssgsea_CombatRef/finalmodel_17datasets_ssgseanorm/metaAnalysis_filter_removeprior_v2/en/elasticnet_model_pred.Rdata")
sig_filter_prior <- c("exhausted_immucellai", 
                      "merck18_tide",
                      "oe_tme.t.cd4.treg",
                      "mdsc_tcga_caf_immunecls")
features_selected <- sig_select[!is.element(sig_select$immune_sig, sig_filter_prior), ]$immune_sig
data_test_selected <- data_test[c(intersect(features_selected, names(data_test)))]
pdf(paste0(savepath, "/metaAnalysis_filter_removeprior_v2/mouseITS_", enrich_method, "_IGeS.pdf"), 
    width = 15, height = 15)
pheatmap(data_test_selected, scale = 'column', 
         cluster_col = F,
         annotation_row = mouseMetaInfo["Genotype"],
         annotation_col = sig_select[c("flag","Type")])
dev.off()
data_test_selected$tag = "0"
mouse_testtask = makeClassifTask(id = "mouse", 
                                 data = data_test_selected,
                                 target = "tag")
pred_en = predict(m_en, task = mouse_testtask)
pred_mouse = pred_en$data
pred_mouse$Genotype_SampleID = rownames(pred_mouse)
pred_mouse <- merge(mouseMetaInfo, pred_mouse, by = 'Genotype_SampleID' )
pred_mouse = pred_mouse[order(pred_mouse$prob.1, decreasing = T), ]
gt_score = aggregate(pred_mouse['prob.1'], by = list(pred_mouse$Genotype), mean)
pred_mouse$Genotype = factor(pred_mouse$Genotype, levels = gt_score[order(gt_score$prob.1, decreasing = T), ]$Group.1)
p = ggplot(pred_mouse, aes(x=Genotype, y=prob.1)) +
      geom_boxplot() + 
      # geom_jitter(width = 0.1,alpha = 0.2)+
      geom_point(size=1) +
      theme_bw()+   
      theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(paste0(savepath, "/metaAnalysis_filter_removeprior_v2/mouseITS_", enrich_method, "_pred.pdf"), p, width = 10, height = 5)
write.csv(pred_mouse, paste0(savepath, "/metaAnalysis_filter_removeprior_v2/mouseITS_", enrich_method, "_enPred.csv"),
             row.names =F, quote = F)


dir.create(paste0(savepath, "/metaAnalysis_filter_removeprior_v3/"))
load("05_immuneSig_ICBresponse/results_r0.4_ITSssgsea_CombatRef/finalmodel_17datasets_ssgseanorm/metaAnalysis_filter_removeprior_v3/en/elasticnet_model_pred.Rdata")
sig_filter_prior <- c("exhausted_immucellai", 
                      "oe_tme.t.cd4.treg",
                      "mdsc_tcga_caf_immunecls")
features_selected <- sig_select[!is.element(sig_select$immune_sig, sig_filter_prior), ]$immune_sig
data_test_selected <- data_test[c(intersect(features_selected, names(data_test)))]
pdf(paste0(savepath, "/metaAnalysis_filter_removeprior_v3/mouseITS_", enrich_method, "_IGeS.pdf"), 
    width = 15, height = 15)
pheatmap(data_test_selected, scale = 'column', 
         cluster_col = F,
         annotation_row = mouseMetaInfo["Genotype"],
         annotation_col = sig_select[c("flag","Type")])
dev.off()
data_test_selected$tag = "0"
mouse_testtask = makeClassifTask(id = "mouse", 
                                 data = data_test_selected,
                                 target = "tag")
pred_en = predict(m_en, task = mouse_testtask)
pred_mouse = pred_en$data
pred_mouse$Genotype_SampleID = rownames(pred_mouse)
pred_mouse <- merge(mouseMetaInfo, pred_mouse, by = 'Genotype_SampleID' )
pred_mouse = pred_mouse[order(pred_mouse$prob.1, decreasing = T), ]
gt_score = aggregate(pred_mouse['prob.1'], by = list(pred_mouse$Genotype), mean)
pred_mouse$Genotype = factor(pred_mouse$Genotype, levels = gt_score[order(gt_score$prob.1, decreasing = T), ]$Group.1)
p = ggplot(pred_mouse, aes(x=Genotype, y=prob.1)) +
      geom_boxplot() + 
      # geom_jitter(width = 0.1,alpha = 0.2)+
      geom_point(size=1) +
      theme_bw()+   
      theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(paste0(savepath, "/metaAnalysis_filter_removeprior_v3/mouseITS_", enrich_method, "_pred.pdf"), p, width = 10, height = 5)
write.csv(pred_mouse, paste0(savepath, "/metaAnalysis_filter_removeprior_v3/mouseITS_", enrich_method, "_enPred.csv"),
             row.names =F, quote = F)


dir.create(paste0(savepath, "/metaAnalysis_filter_removeprior_v4/"))
load("05_immuneSig_ICBresponse/results_r0.4_ITSssgsea_CombatRef/finalmodel_17datasets_ssgseanorm/metaAnalysis_filter_removeprior_v4/en/elasticnet_model_pred.Rdata")
sig_filter_prior <- c("merck18_tide")
features_selected <- sig_select[!is.element(sig_select$immune_sig, sig_filter_prior), ]$immune_sig
data_test_selected <- data_test[c(intersect(features_selected, names(data_test)))]
pdf(paste0(savepath, "/metaAnalysis_filter_removeprior_v4/mouseITS_", enrich_method, "_IGeS.pdf"), 
    width = 15, height = 15)
pheatmap(data_test_selected, scale = 'column', 
         cluster_col = F,
         annotation_row = mouseMetaInfo["Genotype"],
         annotation_col = sig_select[c("flag","Type")])
dev.off()
data_test_selected$tag = "0"
mouse_testtask = makeClassifTask(id = "mouse", 
                                 data = data_test_selected,
                                 target = "tag")
pred_en = predict(m_en, task = mouse_testtask)
pred_mouse = pred_en$data
pred_mouse$Genotype_SampleID = rownames(pred_mouse)
pred_mouse <- merge(mouseMetaInfo, pred_mouse, by = 'Genotype_SampleID' )
pred_mouse = pred_mouse[order(pred_mouse$prob.1, decreasing = T), ]
gt_score = aggregate(pred_mouse['prob.1'], by = list(pred_mouse$Genotype), mean)
pred_mouse$Genotype = factor(pred_mouse$Genotype, levels = gt_score[order(gt_score$prob.1, decreasing = T), ]$Group.1)
p = ggplot(pred_mouse, aes(x=Genotype, y=prob.1)) +
      geom_boxplot() + 
      # geom_jitter(width = 0.1,alpha = 0.2)+
      geom_point(size=1) +
      theme_bw()+   
      theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(paste0(savepath, "/metaAnalysis_filter_removeprior_v4/mouseITS_", enrich_method, "_pred.pdf"), p, width = 10, height = 5)
write.csv(pred_mouse, paste0(savepath, "/metaAnalysis_filter_removeprior_v4/mouseITS_", enrich_method, "_enPred.csv"),
             row.names =F, quote = F)
             
dir.create(paste0(savepath, "/metaAnalysis/"))
load("05_immuneSig_ICBresponse/results_r0.4_ITSssgsea_CombatRef/finalmodel_17datasets_ssgseanorm/metaAnalysis/en/elasticnet_model_pred.Rdata")
en_model_para = m_en$learner.model$glmnet.fit
m_en_imp <- data.frame(en_model_para$beta[,which(en_model_para$lambda == m_en$learner.model$lambda.min)])
rownames(m_en_imp)

data_test$pdc_xcell=0
data_test$monocyte_immucellai=0
data_test$macrophage.m1_cibersort=0

data_test_selected <- data_test[c(intersect(rownames(m_en_imp), names(data_test)))]
pdf(paste0(savepath, "/metaAnalysis/mouseITS_", enrich_method, "_IGeS.pdf"), 
    width = 15, height = 15)
pheatmap(data_test_selected, scale = 'column', 
         cluster_col = F,
         annotation_row = mouseMetaInfo["Genotype"])
dev.off()
data_test_selected$tag = "0"
mouse_testtask = makeClassifTask(id = "mouse", 
                                 data = data_test_selected,
                                 target = "tag")
pred_en = predict(m_en, task = mouse_testtask)
pred_mouse = pred_en$data
pred_mouse$Genotype_SampleID = rownames(pred_mouse)
pred_mouse <- merge(mouseMetaInfo, pred_mouse, by = 'Genotype_SampleID' )
pred_mouse = pred_mouse[order(pred_mouse$prob.1, decreasing = T), ]
gt_score = aggregate(pred_mouse['prob.1'], by = list(pred_mouse$Genotype), mean)
pred_mouse$Genotype = factor(pred_mouse$Genotype, levels = gt_score[order(gt_score$prob.1, decreasing = T), ]$Group.1)
p = ggplot(pred_mouse, aes(x=Genotype, y=prob.1)) +
      geom_boxplot() + 
      # geom_jitter(width = 0.1,alpha = 0.2)+
      geom_point(size=1) +
      theme_bw()+   
      theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(paste0(savepath, "/metaAnalysis/mouseITS_", enrich_method, "_pred.pdf"), p, width = 10, height = 5)
write.csv(pred_mouse, paste0(savepath, "/metaAnalysis/mouseITS_", enrich_method, "_enPred.csv"),
             row.names =F, quote = F)
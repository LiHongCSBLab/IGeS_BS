# rm(list=ls())
# Tumor genes - infiltrated immune cell fraction
workpath <- "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/"
setwd(workpath)

options(stringsAsFactors = F)

# Checking and loading R packages -----------------------------------------

if (!require(readxl)) {  install.packages('readxl')}
if (!require(data.table)) {  install.packages('data.table')}
if (!require(dplyr)) {  install.packages('dplyr')}
# if (!require(TCGAbiolinks)) {  BiocManager::install("TCGAbiolinks") }
if (!require(SummarizedExperiment)) {  BiocManager::install("SummarizedExperiment")}
if (!require(GSVA)) {  BiocManager::install("GSVA")}
if (!require(limma)) {  BiocManager::install("limma")}
if (!require(ConsensusClusterPlus)) {  BiocManager::install("ConsensusClusterPlus")}
if (!require(immunedeconv)) {  devtools::install_github("grst/immunedeconv")}
library(ROCR)
library(reshape2)

# Loading relied functions -----------------------------------------------------
# TIL estimation, TIDE, TIP, IPS and Immune Resistance Pragram
source("02_tumorGene_immuneSig_correlation/scripts/functions/02_immuneSig_generation_function.R")
source("02_tumorGene_immuneSig_correlation/scripts/functions/TILestimator.R")
source("02_tumorGene_immuneSig_correlation/scripts/functions/predictor_TIDEandTIP.R")
source("02_tumorGene_immuneSig_correlation/scripts/functions/predictor_IPS.R")
source("02_tumorGene_immuneSig_correlation/scripts/functions/predictor_icb_ICI_resistance_program.R")
source("02_tumorGene_immuneSig_correlation/scripts/functions/immuneSig_sig160.R")
source("02_tumorGene_immuneSig_correlation/scripts/functions/immuneSig_TCGA_caf_immuneCLSssgsea.R")
# TIGS, AUC and Wilcox test
source("05_immuneSig_ICBresponse/scripts/functions/05_immunotherapy_immuneSig_analysis_function.R")
source("05_immuneSig_ICBresponse/scripts/functions/predictor_ITS.R")

source("06_results_summaryandplotting/scripts/immunesig_selection_function.R")


ITSselected2 <- read.table(paste0("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/06_results_summaryandplotting/data/immune_sig_selected_new/selected_ITS_0.4.txt"),
                           sep = "\t", header = T)


result_dir = "05_immuneSig_ICBresponse/results_r0.4_ITSssgsea_noCombat/"

datasets <- list.files(result_dir)
if(length(grep(".csv", datasets))>0) datasets <- datasets[-grep(".csv", datasets)]
if(length(grep(".Rdata", datasets))>0) datasets <- datasets[-grep(".Rdata", datasets)]
if(length(grep(".pdf", datasets))>0) datasets <- datasets[-grep(".pdf", datasets)]
if(length(grep("ITS_wilcox_result", datasets))>0) datasets <- datasets[-grep("ITS_wilcox_result", datasets)]
if(length(grep("immunesig_icbResponse_all_results", datasets))>0) datasets <- datasets[-grep("immunesig_icbResponse_all_results", datasets)]
if(length(grep("08_Lauss_dataset_outcome", datasets))>0) datasets <- datasets[-grep("08_Lauss_dataset_outcome", datasets)]
if(length(grep("validation_results", datasets))>0) datasets <- datasets[-grep("validation_results", datasets)]
if(length(grep("metaITS_pca", datasets))>0) datasets <- datasets[-grep("metaITS_pca", datasets)]
if(length(grep("ITS_pca", datasets))>0) datasets <- datasets[-grep("ITS_pca", datasets)]
# datasets <- datasets[-grep("ridge_gmm_model", datasets)]

merged_ITSsig <- lapply(as.list(seq(length(datasets))), function(i){
  
  dataset = datasets[i]
  
  merged_res_selected <- read.csv(paste0(result_dir, "/", dataset,"/ITSasPredictor_oe_original/ITSp_ssgsea.csv"))
  dim(merged_res_selected)
  merged_res_scaled <- merged_res_selected #as.data.frame(scale(merged_res_selected[-c(1:2)]))
  merged_res_scaled <- cbind(merged_res_selected[c(1,2)], merged_res_scaled)
  merged_res_scaled$dataset = dataset

return(merged_res_scaled)
})

names(merged_ITSsig) = datasets
lapply(merged_ITSsig, dim)
merged_ITSsig_name <- lapply(merged_ITSsig, names)

# ------------------------------------------------------------------------------
# # validation
# ------------------------------------------------------------------------------
validatesets <- list.files(paste0(result_dir, "/validation_results/"))


merged_ITSsig_validate <- lapply(as.list(seq(length(validatesets))), function(i){
  
  dataset = validatesets[i]
  
  merged_res_selected <- read.csv(paste0(result_dir, "validation_results/", dataset,"/ITSasPredictor_oe_original/ITSp_ssgsea.csv"))
  dim(merged_res_selected)
  merged_res_scaled <- merged_res_selected #as.data.frame(scale(merged_res_selected[-c(1:2)]))
  merged_res_scaled <- cbind(merged_res_selected[c(1,2)], merged_res_scaled)
  merged_res_scaled$dataset = dataset

return(merged_res_scaled)
})

names(merged_ITSsig_validate) = validatesets
lapply(merged_ITSsig_validate, dim)
merged_ITSsig_validate_name <- lapply(merged_ITSsig_validate, names)
# # merged_ITSsig

merged_ITSsig <- c(merged_ITSsig, merged_ITSsig_validate)




sharedITS <- intersect(intersect(intersect(intersect(intersect(
    intersect(
      intersect(
        intersect(
          intersect(
            intersect(
              intersect(
                intersect(
                  intersect(
                    intersect(names(merged_ITSsig[[1]]), names(merged_ITSsig[[2]])),
                    names(merged_ITSsig[[3]])),
                  names(merged_ITSsig[[4]])),
                names(merged_ITSsig[[5]])),
              names(merged_ITSsig[[6]])),
            names(merged_ITSsig[[7]])),
          names(merged_ITSsig[[8]])),
        names(merged_ITSsig[[9]])),
      names(merged_ITSsig[[10]])),
    names(merged_ITSsig[[11]])),
    names(merged_ITSsig[[12]])),
    names(merged_ITSsig[[13]])),
    names(merged_ITSsig[[14]])),
    names(merged_ITSsig[[15]]))

# sharedITS <- c("patientID", "label", intersect(sharedITS,ITSselected2$immune_sig),"dataset")
merged_ITSsig <- lapply(merged_ITSsig, function(x){x[sharedITS]})
lapply(merged_ITSsig, dim)

save(merged_ITSsig, file = paste0(result_dir, "merged_ITSsig_res_oe_original.Rdata"))


# ----------------------------------------------------------
# ----------------------------------------------------------
# feature pca
# ----------------------------------------------------------
# ----------------------------------------------------------


meta_result <- read.csv("05_immuneSig_ICBresponse/results/meta_results/res_metaRes.csv")
meta_selected <- meta_result[meta_result$Meta_Pval < 0.05, ]
meta_selected$immune_sig = meta_selected$X


data_dir = "05_immuneSig_ICBresponse/results_r0.4_ITSssgsea_noCombat/"
load(paste0(data_dir, "merged_ITSsig_res_oe_original.Rdata"))

merged_sig = do.call(rbind, merged_ITSsig)
# merged_sig = merged_sig[!is.element(merged_sig$dataset, c("07_IMvigor210_bladder_immunotherapy", "Braun_dataset")), ]

result_dir = "05_immuneSig_ICBresponse/results_r0.4_ITSssgsea_noCombat/ITS_pca/"
dir.create(result_dir)


mMetadataMerged <- read.csv("05_immuneSig_ICBresponse/dataset_annot.csv")
datasetAnnot <- unique(mMetadataMerged[c(2:4)])
validatesetAnnot <- data.frame(dataset = c("08_Lauss_dataset_outcome", "GSE115821_pre_outcome","GSE96619_outcome","Liu_dataset_CTLA4Naive_outcome" ),
                               cancer = rep("SKCM",4 ),
                               sample_type = rep("Metastatic",4 ))
                               
datasetAnnot<- rbind(datasetAnnot, validatesetAnnot)                           
merged_sig_annot = inner_join(datasetAnnot, merged_sig)
# library(fixbatch)
library(sva)
# merged_sig <- merged_sig[merged_sig$dataset != "Braun_dataset", ]
rownames(merged_sig_annot) <- paste0(merged_sig_annot$dataset, merged_sig_annot$patient)
merged_sig_annot <- merged_sig_annot[!is.na(merged_sig_annot$label), ]
merged_sig_scaled <- merged_sig_annot[-c(1:2)]
merged_sig_scaled <- merged_sig_scaled[intersect(names(merged_sig_scaled), meta_selected$immune_sig)]
# merged_sig_scaled <- as.data.frame(scale(merged_sig_annot[-c(1:5)]))

mMetadataMerged <- merged_sig_annot[c("label","dataset","cancer","sample_type")]
mMetadataMerged$SampleID <- rownames(merged_sig_annot)


library(pheatmap)
mMetadataMerged <- mMetadataMerged[order(mMetadataMerged$label), ]
pdf(paste0(result_dir, "meta_selected_ITS_heatmap.pdf"))
pheatmap(merged_sig_scaled[rownames(mMetadataMerged),], 
          cluster_rows = F,
         annotation_row = mMetadataMerged['label'])
dev.off()


# # analysis
# pvcaObj <- pvca(merged_sig_scaled, 
#                 merged_sig[c("label","dataset","cancer","sample_type")], 
#                 threshold = 0.7,inter = FALSE)
# # plot
# # pdf("result_batchEffect_rerun_20220117/batcheffect/count_pvca_factor.pdf")
# pvca_plot(pvcaObj)
# # dev.off()

# simple data consistency ------------------------------------------------------
# with Count
library(FactoMineR)
library(ggplot2)
library(patchwork)

pca_fmr = PCA(merged_sig_scaled, 
              scale.unit = T, ncp = 6, graph = F)  ## !! just one step !! ##

## generate PCA using ggplot2
PCs.2dim = pca_fmr$ind$coord[, 1:2]    # 读取前两维的PC的二维轴平面的坐标
colnames(PCs.2dim) = c("Dim1", "Dim2")
xlab = paste("Dim1 (", round(pca_fmr$eig[1,2], 2), "%)", sep = ""); xlab   # 设置x轴标题
ylab = paste("Dim2 (", round(pca_fmr$eig[2,2], 2), "%)", sep = ""); ylab   # 设置y轴标题

## plot results of PCA with ggplot2
PCs.2dim = as.data.frame(PCs.2dim)
PCs.2dim$SampleID = rownames(PCs.2dim)
Fig1.PCs = merge(mMetadataMerged, PCs.2dim, by = "SampleID") 

p_dataset <- 
  ggplot(Fig1.PCs, aes(x=Dim1, y=Dim2, color= dataset )) + 
  geom_point(size = 1) +
  theme_bw() +
  xlab(xlab) + ylab(ylab) + 
  scale_alpha_discrete(range = c(0.5,1)) +
  coord_fixed() +  
  ggtitle("dataset")

p_cancer <- 
  ggplot(Fig1.PCs, aes(x=Dim1, y=Dim2, color= cancer )) + 
  geom_point(size = 1) +
  theme_bw() +
  xlab(xlab) + ylab(ylab) + 
  scale_alpha_discrete(range = c(0.5,1)) +
  coord_fixed() +  
  ggtitle("cancer")

p_label <- 
  ggplot(Fig1.PCs, aes(x=Dim1, y=Dim2, color= label )) + 
  geom_point(size = 1) +
  theme_bw() +
  xlab(xlab) + ylab(ylab) + 
  scale_alpha_discrete(range = c(0.5,1)) +
  coord_fixed() +  
  ggtitle("label")


p = p_dataset | p_cancer | p_label
ggsave(paste0(result_dir, "meta_selected_ITS_pca.pdf"), p, width = 15, height = 5)


library(Rtsne)
set.seed(42) # Set a seed if you want reproducible results
tsne_out <- Rtsne(merged_sig_scaled) # Run TSNE
# mMetadataMerged 


# Show the objects in the 2D tsne representation
# plot(tsne_out$Y,col=mMetadataMerged$cancer)

tsne_plot <- data.frame(x = tsne_out$Y[,1], 
                        y = tsne_out$Y[,2], 
                        col = mMetadataMerged$cancer)
p_tsne <- ggplot(tsne_plot) + 
          geom_point(aes(x=x, y=y, color=col))+
          theme_bw()
ggsave(paste0(result_dir, "meta_selected_ITStsne.pdf"), p_tsne, width = 5, height = 5)


  # single dataset pca

library(FactoMineR)
library(ggplot2)
library(patchwork)
data_dir = "05_immuneSig_ICBresponse/results_r0.4_ITSssgsea_noCombat/"
load(paste0(data_dir, "merged_ITSsig_res_oe_original.Rdata"))
p_label <- list()

for(i in seq(length(merged_ITSsig))){
  res = merged_ITSsig[[i]]
  df = res[intersect(names(res), meta_selected$immune_sig)]
  rownames(df) = res$patientID
  df_annot = res[c("patientID", "label" ,"dataset")]

  pca_fmr = PCA(df, 
                scale.unit = T, ncp = 6, graph = F)  ## !! just one step !! ##

  ## generate PCA using ggplot2
  PCs.2dim = pca_fmr$ind$coord[, c(1,2)]    # 读取前两维的PC的二维轴平面的坐标
  colnames(PCs.2dim) = c("Dim1", "Dim2")
  xlab = paste("Dim1 (", round(pca_fmr$eig[1,2], 2), "%)", sep = ""); xlab   # 设置x轴标题
  ylab = paste("Dim2 (", round(pca_fmr$eig[2,2], 2), "%)", sep = ""); ylab   # 设置y轴标题

  ## plot results of PCA with ggplot2
  PCs.2dim = as.data.frame(PCs.2dim)
  PCs.2dim$patientID = rownames(PCs.2dim)
  Fig1.PCs = merge(df_annot, PCs.2dim, by = "patientID") 


  p_label[[i]] <- 
    ggplot(Fig1.PCs, aes(x=Dim1, y=Dim2, color= label )) + 
    geom_point(size = 1) +
    theme_bw() +
    xlab(xlab) + ylab(ylab) + 
    scale_alpha_discrete(range = c(0.5,1)) +
    coord_fixed() +  
    ggtitle(paste0(names(merged_ITSsig)[i], " - label"))
}

pdf(paste0(result_dir, "singledataset_metaITSsig_pca.pdf"),
            width = 5,
            height = 5, 
            onefile = TRUE)
for(x in seq(length(p_label))){
    print(p_label[[x]])
}
dev.off()
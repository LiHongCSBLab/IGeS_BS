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
library(Combat)

# Loading relied functions -----------------------------------------------------
# TIGS, AUC and Wilcox test
source("05_immuneSig_ICBresponse/scripts/functions/05_immunotherapy_immuneSig_analysis_function.R")
source("05_immuneSig_ICBresponse/scripts/functions/predictor_ITS.R")
source("06_results_summaryandplotting/scripts/immunesig_selection_function.R")


# ----------------------------------------------------------
# ----------------------------------------------------------
# merge expression matrix (log2)   -------------------------
# ----------------------------------------------------------
# ----------------------------------------------------------
data <- list()
datasetname <- list()
dataAnnot <- list()
# create savepath -----------------------------------------------------------------
dir.create("05_immuneSig_ICBresponse/results_r0.4_ITSssgsea_Combat/")

# Hugo_dataset -----------------------------------------------------------------
processedDataPath = "05_immuneSig_ICBresponse/data/processedData/"
dataset = "01_Hugo_dataset_outcome.Rdata"
load(paste0(processedDataPath, "/", dataset))
result_dir = "05_immuneSig_ICBresponse/results_r0.4_ITSssgsea_Combat/01_Hugo_dataset_outcome/"
cancer = "SKCM"
sample_type = "Metastatic"
data[[1]] =  as.matrix(log2(1+tpm_mat))
datasetname[[1]] = '01_Hugo_dataset_outcome'
colnames(data[[1]]) <- paste0(datasetname[[1]],  colnames(data[[1]]))
dataAnnot[[1]] <- response
dataAnnot[[1]]$SampleID <- paste0(datasetname[[1]],  dataAnnot[[1]]$patientID)
dataAnnot[[1]]$cancer <- cancer
dataAnnot[[1]]$sample_type <- sample_type

# Riaz_dataset -----------------------------------------------------------------
processedDataPath = "05_immuneSig_ICBresponse/data/processedData/"
dataset = "02_Riaz_dataset_outcome.Rdata"
load(paste0(processedDataPath, "/", dataset))
result_dir = "05_immuneSig_ICBresponse/results_r0.4_ITSssgsea_Combat/02_Riaz_dataset_outcome/"
cancer = "SKCM"
sample_type = "Metastatic"
data[[2]] =  as.matrix(log2(1+tpm_mat))
datasetname[[2]] = '02_Riaz_dataset_outcome'
colnames(data[[2]]) <- paste0(datasetname[[2]],  colnames(data[[2]]))
dataAnnot[[2]] <- response
dataAnnot[[2]]$SampleID <- paste0(datasetname[[2]],  dataAnnot[[2]]$patientID)
dataAnnot[[2]]$cancer <- cancer
dataAnnot[[2]]$sample_type <- sample_type

# Gide --------------------------------------------------------------------------
processedDataPath = "05_immuneSig_ICBresponse/data/processedData/"
dataset = "05_Gide_dataset_outcome.Rdata"
load(paste0(processedDataPath, "/", dataset))
result_dir = "05_immuneSig_ICBresponse/results_r0.4_ITSssgsea_Combat/05_Gide_dataset_outcome/"
cancer = "SKCM"
sample_type = "Metastatic"

data[[3]] =  as.matrix(log2(1+tpm_mat))
datasetname[[3]] = '05_Gide_dataset_outcome'
colnames(data[[3]]) <- paste0(datasetname[[3]],  colnames(data[[3]]))
dataAnnot[[3]] <- response
dataAnnot[[3]]$SampleID <- paste0(datasetname[[3]],  dataAnnot[[3]]$patientID)
dataAnnot[[3]]$cancer <- cancer
dataAnnot[[3]]$sample_type <- sample_type

# Kim_dataset ------------------------------------------------------------------
processedDataPath = "05_immuneSig_ICBresponse/data/processedData/"
dataset = "06_Kim_dataset_outcome.Rdata"
load(paste0(processedDataPath, "/", dataset))
result_dir = "05_immuneSig_ICBresponse/results_r0.4_ITSssgsea_Combat/06_Kim_dataset_outcome/"
cancer = "STAD"
sample_type = "Primary"
data[[4]] =  as.matrix(log2(1+tpm_mat))
datasetname[[4]] = '06_Kim_dataset_outcome'
colnames(data[[4]]) <- paste0(datasetname[[4]],  colnames(data[[4]]))
dataAnnot[[4]] <- response
dataAnnot[[4]]$SampleID <- paste0(datasetname[[4]],  dataAnnot[[4]]$patientID)
dataAnnot[[4]]$cancer <- cancer
dataAnnot[[4]]$sample_type <- sample_type

# IMvigor210_dataset ------------------------------------------------------------------
dataset = "IMvigor210_bladder.Rdata"
processedDataPath = "/picb/bigdata/project/CancerSysBio/Immunotherapy_prediction/data/RNAseq/"
load(paste0(processedDataPath, "/IMvigor210/processedData/", dataset))
result_dir = "05_immuneSig_ICBresponse/results_r0.4_ITSssgsea_Combat/07_IMvigor210_bladder_immunotherapy/"
tpm_mat <- tpm_mat[rowSums(tpm_mat) != 0, ]
tpm_mat <- tpm_mat[-which(rownames(tpm_mat)==""), ]
data[[5]] =  as.matrix(log2(1+tpm_mat))
datasetname[[5]] = '07_IMvigor210_bladder_immunotherapy'
colnames(data[[5]]) <- paste0(datasetname[[5]],  colnames(data[[5]]))
dataAnnot[[5]] <- response
dataAnnot[[5]]$SampleID <- paste0(datasetname[[5]],  dataAnnot[[5]]$patientID)
dataAnnot[[5]]$cancer <- cancer
dataAnnot[[5]]$sample_type <- sample_type

# Braun_dataset_PD1 ------------------------------------------------------------------
processedDataPath = "/picb/bigdata/project/CancerSysBio/Immunotherapy_prediction/data/RNAseq/"
dataset = "Braun_dataset_PD1.Rdata"
load(paste0(processedDataPath, "/Braun_dataset/processedData/", dataset))
result_dir = "05_immuneSig_ICBresponse/results_r0.4_ITSssgsea_Combat/Braun_dataset/"
cancer = "KIRC"
sample_type = "Primary"
data[[6]] =  as.matrix(log2(1+tpm_mat))
datasetname[[6]] = 'Braun_dataset'
colnames(data[[6]]) <- paste0(datasetname[[6]],  colnames(data[[6]]))
dataAnnot[[6]] <- response
dataAnnot[[6]]$SampleID <- paste0(datasetname[[6]],  dataAnnot[[6]]$patientID)
dataAnnot[[6]]$cancer <- cancer
dataAnnot[[6]]$sample_type <- sample_type

# GSE126044 ----------------------------------------------------------
processedDataPath = "/picb/bigdata/project/CancerSysBio/Immunotherapy_prediction/data/RNAseq/"
dataset = "GSE126044.Rdata"
load(paste0(processedDataPath, "/GSE126044/processedData/", dataset))
result_dir = "05_immuneSig_ICBresponse/results_r0.4_ITSssgsea_Combat/GSE126044/"
cancer = "LUAD"
sample_type = "Primary"
if(length(which(rowSums(tpm_mat) == 0)) > 0)tpm_mat <- tpm_mat[rowSums(tpm_mat) != 0, ]
if(length(which(rownames(tpm_mat)=="")) > 0)tpm_mat <- tpm_mat[-which(rownames(tpm_mat)==""), ]
data[[7]] =  as.matrix(log2(1+tpm_mat))
datasetname[[7]] = 'GSE126044'
colnames(data[[7]]) <- paste0(datasetname[[7]],  colnames(data[[7]]))
dataAnnot[[7]] <- response
dataAnnot[[7]]$SampleID <- paste0(datasetname[[7]],  dataAnnot[[7]]$patientID)
dataAnnot[[7]]$cancer <- cancer
dataAnnot[[7]]$sample_type <- sample_type

# GSE135222 ----------------------------------------------------------
processedDataPath = "/picb/bigdata/project/CancerSysBio/Immunotherapy_prediction/data/RNAseq/"
dataset = "GSE135222.Rdata"
load(paste0(processedDataPath, "/GSE135222/processedData/", dataset))
result_dir = "05_immuneSig_ICBresponse/results_r0.4_ITSssgsea_Combat/GSE135222/"
cancer = "LUAD"
sample_type = "Primary"

if(length(which(rowSums(tpm_mat) == 0)) > 0)tpm_mat <- tpm_mat[rowSums(tpm_mat) != 0, ]
if(length(which(rownames(tpm_mat)=="")) > 0)tpm_mat <- tpm_mat[-which(rownames(tpm_mat)==""), ]
data[[8]] =  as.matrix(log2(1+tpm_mat))
datasetname[[8]] = 'GSE135222'
colnames(data[[8]]) <- paste0(datasetname[[8]],  colnames(data[[8]]))
dataAnnot[[8]] <- response
dataAnnot[[8]]$SampleID <- paste0(datasetname[[8]],  dataAnnot[[8]]$patientID)
dataAnnot[[8]]$cancer <- cancer
dataAnnot[[8]]$sample_type <- sample_type


                              
# GSE145996 ----------------------------------------------------------
processedDataPath = "/picb/bigdata/project/CancerSysBio/Immunotherapy_prediction/data/RNAseq/"
dataset = "GSE145996.Rdata"
load(paste0(processedDataPath, "/GSE145996/processedData/", dataset))
result_dir = "05_immuneSig_ICBresponse/results_r0.4_ITSssgsea_Combat/GSE145996/"
cancer = "SKCM"
sample_type = "Metastatic"
if(length(which(rowSums(tpm_mat) == 0)) > 0)tpm_mat <- tpm_mat[rowSums(tpm_mat) != 0, ]
if(length(which(rownames(tpm_mat)=="")) > 0)tpm_mat <- tpm_mat[-which(rownames(tpm_mat)==""), ]
data[[9]] =  as.matrix(log2(1+tpm_mat))
datasetname[[9]] = 'GSE145996'
colnames(data[[9]]) <- paste0(datasetname[[9]],  colnames(data[[9]]))
dataAnnot[[9]] <- response
dataAnnot[[9]]$SampleID <- paste0(datasetname[[9]],  dataAnnot[[9]]$patientID)
dataAnnot[[9]]$cancer <- cancer
dataAnnot[[9]]$sample_type <- sample_type


# Nathanson ----------------------------------------------------------
processedDataPath = "/picb/bigdata/project/CancerSysBio/Immunotherapy_prediction/data/"
dataset = "Nathanson_dataset.Rdata"
load(paste0(processedDataPath, "RNAseq/Nathanson_dataset/processedData/", dataset))
result_dir = "05_immuneSig_ICBresponse/results_r0.4_ITSssgsea_Combat/Nathanson_dataset_outcome/"
cancer = "SKCM"
sample_type = "Metastatic"
if(length(which(rowSums(tpm_mat) == 0)) > 0)tpm_mat <- tpm_mat[rowSums(tpm_mat) != 0, ]
if(length(which(rownames(tpm_mat)=="")) > 0)tpm_mat <- tpm_mat[-which(rownames(tpm_mat)==""), ]
data[[10]] =  as.matrix(log2(1+tpm_mat))
datasetname[[10]] = 'Nathanson_dataset_outcome'
colnames(data[[10]]) <- paste0(datasetname[[10]],  colnames(data[[10]]))
dataAnnot[[10]] <- response
dataAnnot[[10]]$SampleID <- paste0(datasetname[[10]],  dataAnnot[[10]]$patientID)
dataAnnot[[10]]$cancer <- cancer
dataAnnot[[10]]$sample_type <- sample_type

# PMID_29301960 ----------------------------------------------------------
processedDataPath = "/picb/bigdata/project/CancerSysBio/Immunotherapy_prediction/data/"
dataset = "PMID_29301960.Rdata"
load(paste0(processedDataPath, "RNAseq/PMID_29301960/processedData/", dataset))
result_dir = "05_immuneSig_ICBresponse/results_r0.4_ITSssgsea_Combat/PMID_29301960/"
cancer = "KIRC"
sample_type = "Primary"
if(length(which(rowSums(tpm_mat) == 0)) > 0)tpm_mat <- tpm_mat[rowSums(tpm_mat) != 0, ]
if(length(which(rownames(tpm_mat)=="")) > 0)tpm_mat <- tpm_mat[-which(rownames(tpm_mat)==""), ]

data[[11]] =  as.matrix(log2(1+tpm_mat))
datasetname[[11]] = 'PMID_29301960'
colnames(data[[11]]) <- paste0(datasetname[[11]],  colnames(data[[11]]))
dataAnnot[[11]] <- response
dataAnnot[[11]]$SampleID <- paste0(datasetname[[11]],  dataAnnot[[11]]$patientID)
dataAnnot[[11]]$cancer <- cancer
dataAnnot[[11]]$sample_type <- sample_type

# Zhao_dataset ------------------------------------------------------------------
processedDataPath = "/picb/bigdata/project/CancerSysBio/Immunotherapy_prediction/data/RNAseq/"

dataset = "Zhao_dataset.Rdata"
load(paste0(processedDataPath, "/Zhao_dataset/processedData/", dataset))
result_dir = "05_immuneSig_ICBresponse/results_r0.4_ITSssgsea_Combat/Zhao_dataset/"
cancer = "GBM"
sample_type = "Primary"
data[[12]] =  as.matrix(log2(1+tpm_mat))
datasetname[[12]] = 'Zhao_dataset'
colnames(data[[12]]) <- paste0(datasetname[[12]],  colnames(data[[12]]))
dataAnnot[[12]] <- response
dataAnnot[[12]]$SampleID <- paste0(datasetname[[12]],  dataAnnot[[12]]$patientID)
dataAnnot[[12]]$cancer <- cancer
dataAnnot[[12]]$sample_type <- sample_type


# ------------------------------------------------------------------------------
# validation sets             --------------------------------------------------
# ------------------------------------------------------------------------------

# Lauss-ACT therapy  ----------------------------------------------------------
# # GSE100797 ------------------------------------------------------------------
# dataset = "GSE100797.Rdata"
# load(paste0(processedDataPath, "/Lauss_dataset/processedData/", dataset))

processedDataPath = "05_immuneSig_ICBresponse/data/processedData/"
dataset = "08_Lauss_melanoma.Rdata"
load(paste0(processedDataPath, "/", dataset))
result_dir = "05_immuneSig_ICBresponse/results_r0.4_ITSssgsea_Combat/08_Lauss_dataset_outcome/"
cancer = "SKCM"
sample_type = "Metastatic"
data[[13]] =  as.matrix(log2(1+tpm_mat))
datasetname[[13]] = '08_Lauss_dataset_outcome'
colnames(data[[13]]) <- paste0(datasetname[[13]],  colnames(data[[13]]))
dataAnnot[[13]] <- response
dataAnnot[[13]]$SampleID <- paste0(datasetname[[13]],  dataAnnot[[13]]$patientID)
dataAnnot[[13]]$cancer <- cancer
dataAnnot[[13]]$sample_type <- sample_type

# Liu_dataset_PD1 ------------------------------------------------------------------
processedDataPath = "/picb/bigdata/project/CancerSysBio/Immunotherapy_prediction/data/RNAseq/"
dataset = "Liu_dataset.Rdata"
load(paste0(processedDataPath, "/Liu_dataset/processedData/", dataset))
response = response_CTLA4Naive_MAPKTxNaive_skin[c('patient', 'label')]
names(response) = c("patientID", "label")
tpm_mat <- tpm_mat[intersect(names(tpm_mat), response$patient)]
tpm_mat <- tpm_mat[-grep("Patient143", colnames(tpm_mat))]
response = response[is.element(response$patient, intersect(names(tpm_mat), response$patient)),]
cancer = "SKCM"
sample_type = "Metastatic"
data[[14]] =  as.matrix(log2(1+tpm_mat))
datasetname[[14]] = 'Liu_dataset_CTLA4Naive_outcome'
colnames(data[[14]]) <- paste0(datasetname[[14]],  colnames(data[[14]]))
dataAnnot[[14]] <- response
dataAnnot[[14]]$SampleID <- paste0(datasetname[[14]],  dataAnnot[[14]]$patientID)
dataAnnot[[14]]$cancer <- cancer
dataAnnot[[14]]$sample_type <- sample_type


# GSE115821 --------------------------------------------------------------------------
# dataset match : 05_immuneSig_ICBresponse/scripts/scripts_for_icb_datasets/dataset_reform.R
processedDataPath = "/picb/bigdata/project/CancerSysBio/Immunotherapy_prediction/data/"

dataset = "GSE115821_pre.Rdata"
load(paste0(processedDataPath, "RNAseq/GSE115821/processedData/",dataset))
tpm_mat <- tpm_mat_pre
response <- response_pre[c("patientID", "label")]
names(response) <- c("patientID", "label")
if(length(which(rowSums(tpm_mat) == 0)) > 0)tpm_mat <- tpm_mat[rowSums(tpm_mat) != 0, ]
if(length(which(rownames(tpm_mat)=="")) > 0)tpm_mat <- tpm_mat[-which(rownames(tpm_mat)==""), ]

result_dir = "05_immuneSig_ICBresponse/results_r0.4_ITSssgsea_Combat/validation_results/GSE115821_pre_outcome/"
dir.create(result_dir)
cancer = "SKCM"
sample_type = "Metastatic" # NOT SURE
data[[15]] =  as.matrix(log2(1+tpm_mat))
datasetname[[15]] = 'GSE115821_pre_outcome'
colnames(data[[15]]) <- paste0(datasetname[[15]],  colnames(data[[15]]))
dataAnnot[[15]] <- response
dataAnnot[[15]]$SampleID <- paste0(datasetname[[15]],  dataAnnot[[15]]$patientID)
dataAnnot[[15]]$cancer <- cancer
dataAnnot[[15]]$sample_type <- sample_type


# GSE96619 --------------------------------------------------------------------------
# dataset match : 05_immuneSig_ICBresponse/scripts/scripts_for_icb_datasets/dataset_reform.R
processedDataPath = "/picb/bigdata/project/CancerSysBio/Immunotherapy_prediction/data/"

dataset = "GSE96619_pre.Rdata"
load(paste0(processedDataPath, "RNAseq/GSE96619/processedData/",dataset))
tpm_mat <- tpm_pre
response <- response_pre
names(response) <- c("patientID", "label")
if(length(which(rowSums(tpm_mat) == 0)) > 0)tpm_mat <- tpm_mat[rowSums(tpm_mat) != 0, ]
if(length(which(rownames(tpm_mat)=="")) > 0)tpm_mat <- tpm_mat[-which(rownames(tpm_mat)==""), ]

colnames(tpm_mat) <- gsub('-','_', colnames(tpm_mat))
cancer = "SKCM"
sample_type = "Metastatic" # NOT SURE
data[[16]] =  as.matrix(log2(1+tpm_mat))
datasetname[[16]] = 'GSE96619_outcome'
colnames(data[[16]]) <- paste0(datasetname[[16]],  colnames(data[[16]]))
dataAnnot[[16]] <- response
dataAnnot[[16]]$SampleID <- paste0(datasetname[[16]],  dataAnnot[[16]]$patientID)
dataAnnot[[16]]$cancer <- cancer
dataAnnot[[16]]$sample_type <- sample_type

# GSE176307  --------------------------------------------------------------------------
processedDataPath = "/picb/bigdata/project/CancerSysBio/Immunotherapy_prediction/data/"

dataset = "GSE176307.Rdata"
load(paste0(processedDataPath, "RNAseq/GSE176307/processedData/",dataset))
tpm_mat <- exp_mat

names(response) <- c("patientID", "label")
if(length(which(rowSums(tpm_mat) == 0)) > 0)tpm_mat <- tpm_mat[rowSums(tpm_mat) != 0, ]
if(length(which(rownames(tpm_mat)=="")) > 0)tpm_mat <- tpm_mat[-which(rownames(tpm_mat)==""), ]

colnames(tpm_mat) <- gsub('-','_', colnames(tpm_mat))
cancer = "BLCA"
sample_type = "Primary" 
data[[17]] =  as.matrix(log2(1+tpm_mat))
datasetname[[17]] = 'GSE176307'
colnames(data[[17]]) <- paste0(datasetname[[17]],  colnames(data[[17]]))
dataAnnot[[17]] <- response
dataAnnot[[17]]$SampleID <- paste0(datasetname[[17]],  dataAnnot[[17]]$patientID)
dataAnnot[[17]]$cancer <- cancer
dataAnnot[[17]]$sample_type <- sample_type


# ----------------------------------------------------------
# merge gene expression matrix  ----------------------------
# ----------------------------------------------------------
names(dataAnnot) <- datasetname

names(data) <- datasetname

save(data, dataAnnot, file = "05_immuneSig_ICBresponse/results_r0.4_ITSssgsea_Combat/ICBdataset_Exp17datasets.Rdata")


for(i in 1:length(data)){
  if(i == 1){
    cm = rownames(data[[1]])
  } else{
    cm = intersect(cm,rownames(data[[i]]))
  }
}

data_bc <- lapply(data,function(z){ return(t(z[cm,]))})
data_bc <- do.call(rbind,data_bc)
dataAnnot_all <- do.call(rbind, dataAnnot)
dataAnnot_all$dataset <- str_split_fixed(rownames(dataAnnot_all),'\\.',2)[,1]
data_bc <- as.data.frame(data_bc)
data_bc$SampleID <- row.names(data_bc)
data_bc <- merge(data_bc, dataAnnot_all, by = 'SampleID')
data_bc$batch <- data_bc$dataset


rownames(data_bc) <- data_bc$SampleID
edata <- t(data_bc[-which(is.element(names(data_bc), c("SampleID", "patientID","label","cancer","sample_type","dataset","batch")))])
edata_log <- edata

# ----------------------------------------------------------
# ----------------------------------------------------------
# batch correction with combat  ------------------------------------------------
# ----------------------------------------------------------
# ----------------------------------------------------------
library(sva)

data_ac <- ComBat(dat = edata_log, 
                 batch = data_bc$batch)
data_ac <- as.data.frame(t(data_ac))
# data_ac <- cbind(data_bc[,1:8],data_ac)
# data_ac <- data_ac[!is.na(data_ac$response),] ##remove patient with unknown response status


# ----------------------------------------------------------
# ----------------------------------------------------------
# Expression-level pca  --------------------------------------------------------
# ----------------------------------------------------------
# ----------------------------------------------------------
result_dir <- "05_immuneSig_ICBresponse/results_r0.4_ITSssgsea_Combat/"
rownames(dataAnnot_all) <- dataAnnot_all$SampleID

library(FactoMineR)
library(ggplot2)
library(patchwork)
# before Combat: edata_log

pca_fmr = PCA(t(edata_log), 
              scale.unit = T, ncp = 6, graph = F)  ## !! just one step !! ##

## generate PCA using ggplot2
PCs.2dim = pca_fmr$ind$coord[, 1:2]    # 读取前两维的PC的二维轴平面的坐标
colnames(PCs.2dim) = c("Dim1", "Dim2")
xlab = paste("Dim1 (", round(pca_fmr$eig[1,2], 2), "%)", sep = ""); xlab   # 设置x轴标题
ylab = paste("Dim2 (", round(pca_fmr$eig[2,2], 2), "%)", sep = ""); ylab   # 设置y轴标题

## plot results of PCA with ggplot2
PCs.2dim = as.data.frame(PCs.2dim)
PCs.2dim$SampleID = rownames(PCs.2dim)
Fig1.PCs = merge(dataAnnot_all, PCs.2dim, by = "SampleID") 

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


p_beforeCombat = p_dataset | p_cancer | p_label
ggsave(paste0(result_dir, "Exp_beforeCombat17datasets_pca.pdf"), p_beforeCombat, width = 15, height = 5)


# after combat: data_ac

pca_fmr = PCA(data_ac, 
              scale.unit = T, ncp = 6, graph = F)  ## !! just one step !! ##

## generate PCA using ggplot2
PCs.2dim = pca_fmr$ind$coord[, 1:2]    # 读取前两维的PC的二维轴平面的坐标
colnames(PCs.2dim) = c("Dim1", "Dim2")
xlab = paste("Dim1 (", round(pca_fmr$eig[1,2], 2), "%)", sep = ""); xlab   # 设置x轴标题
ylab = paste("Dim2 (", round(pca_fmr$eig[2,2], 2), "%)", sep = ""); ylab   # 设置y轴标题

## plot results of PCA with ggplot2
PCs.2dim = as.data.frame(PCs.2dim)
PCs.2dim$SampleID = rownames(PCs.2dim)
Fig1.PCs = merge(dataAnnot_all, PCs.2dim, by = "SampleID") 

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


p_afterCombat = p_dataset | p_cancer | p_label
ggsave(paste0(result_dir, "Exp_afterCombat17datasets_pca.pdf"), p_afterCombat, width = 15, height = 5)

# ----------------------------------------------------------
# ----------------------------------------------------------
# ITSssgsea  -------------------------------------------------------------------
# ----------------------------------------------------------
# ----------------------------------------------------------
# before Combat: edata_log
source("05_immuneSig_ICBresponse/scripts/functions/05_immunotherapy_immuneSig_analysis_function.R")
source("05_immuneSig_ICBresponse/scripts/functions/predictor_ITS.R")
source("06_results_summaryandplotting/scripts/immunesig_selection_function.R")

savepath_allresult = "02_tumorGene_immuneSig_correlation/results_selected/ITS_oe_original_spearman_TUMERIC_r0.4/"

savepath = paste0(result_dir, "ITSssgsea_beforeCombat17datasets/")
dir.create(savepath)
datasets = unique(dataAnnot_all$dataset)
for(dataset in datasets){
  # dataset = datasets[1]
  ds = dataAnnot_all[dataAnnot_all$dataset == dataset, ]
  df <- edata_log[, which(is.element(colnames(edata_log),ds$SampleID))]
  cancer = unique(ds$cancer)
  sample_type = unique(ds$sample_type)
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
                            ssgsea.norm = FALSE) 

  ITSp_res <- as.data.frame(t(ITSp))
  ITSp_res$SampleID <- rownames(ITSp_res)
  ITSp_res <- inner_join(ds, ITSp_res, by = "SampleID")
  write.csv(ITSp_res, paste0(savepath, "/", dataset, "_ITSp_", enrich_method, ".csv"), 
            row.names =F, quote = F)


  ITSn <- immune_sigAnalysis(data = df, 
                            gs = genelist_n,
                            method = enrich_method, 
                            kcdf = "Gaussian",
                            ssgsea.norm = FALSE) 
  ITSn_res <- as.data.frame(t(ITSn))
  ITSn_res$SampleID <- rownames(ITSn_res)
  ITSn_res <- inner_join(ds, ITSn_res,by='SampleID')
  write.csv(ITSn_res, paste0(savepath, "/", dataset, "ITSn_", enrich_method, ".csv"), 
            row.names =F, quote = F)


  enrich_method = "ssgsea"
  ITSp <- immune_sigAnalysis(data = df, 
                            gs = genelist_p,
                            method = enrich_method, 
                            kcdf = "Gaussian",
                            ssgsea.norm = TRUE) 

  ITSp_res <- as.data.frame(t(ITSp))
  ITSp_res$SampleID <- rownames(ITSp_res)
  ITSp_res <- inner_join(ds, ITSp_res, by = "SampleID")
  write.csv(ITSp_res, paste0(savepath, "/", dataset, "_ITSp_", enrich_method, "_ssgseanorm.csv"), 
            row.names =F, quote = F)


  ITSn <- immune_sigAnalysis(data = df, 
                            gs = genelist_n,
                            method = enrich_method, 
                            kcdf = "Gaussian",
                            ssgsea.norm = TRUE) 
  ITSn_res <- as.data.frame(t(ITSn))
  ITSn_res$SampleID <- rownames(ITSn_res)
  ITSn_res <- inner_join(ds, ITSn_res,by='SampleID')
  write.csv(ITSn_res, paste0(savepath, "/", dataset, "ITSn_", enrich_method, "_ssgseanorm.csv"), 
            row.names =F, quote = F)

}


# after combat: data_ac
savepath = paste0(result_dir, "ITSssgsea_afterCombat17datasets/")
dir.create(savepath)
datasets = unique(dataAnnot_all$dataset)
data_ac <- t(data_ac)
for(dataset in datasets){
  # dataset = datasets[1]
  ds = dataAnnot_all[dataAnnot_all$dataset == dataset, ]
  df <- data_ac[, which(is.element(colnames(data_ac),ds$SampleID))]
  cancer = unique(ds$cancer)
  sample_type = unique(ds$sample_type)
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
                            ssgsea.norm = FALSE) 

  ITSp_res <- as.data.frame(t(ITSp))
  ITSp_res$SampleID <- rownames(ITSp_res)
  ITSp_res <- inner_join(ds, ITSp_res, by = "SampleID")
  write.csv(ITSp_res, paste0(savepath, "/", dataset, "_ITSp_", enrich_method, ".csv"), 
            row.names =F, quote = F)


  ITSn <- immune_sigAnalysis(data = df, 
                            gs = genelist_n,
                            method = enrich_method, 
                            kcdf = "Gaussian",
                            ssgsea.norm = FALSE) 
  ITSn_res <- as.data.frame(t(ITSn))
  ITSn_res$SampleID <- rownames(ITSn_res)
  ITSn_res <- inner_join(ds, ITSn_res,by='SampleID')
  write.csv(ITSn_res, paste0(savepath, "/", dataset, "ITSn_", enrich_method, ".csv"), 
            row.names =F, quote = F)


  enrich_method = "ssgsea"
  ITSp <- immune_sigAnalysis(data = df, 
                            gs = genelist_p,
                            method = enrich_method, 
                            kcdf = "Gaussian",
                            ssgsea.norm = TRUE) 

  ITSp_res <- as.data.frame(t(ITSp))
  ITSp_res$SampleID <- rownames(ITSp_res)
  ITSp_res <- inner_join(ds, ITSp_res, by = "SampleID")
  write.csv(ITSp_res, paste0(savepath, "/", dataset, "_ITSp_", enrich_method, "_ssgseanorm.csv"), 
            row.names =F, quote = F)


  ITSn <- immune_sigAnalysis(data = df, 
                            gs = genelist_n,
                            method = enrich_method, 
                            kcdf = "Gaussian",
                            ssgsea.norm = TRUE) 
  ITSn_res <- as.data.frame(t(ITSn))
  ITSn_res$SampleID <- rownames(ITSn_res)
  ITSn_res <- inner_join(ds, ITSn_res,by='SampleID')
  write.csv(ITSn_res, paste0(savepath, "/", dataset, "ITSn_", enrich_method, "_ssgseanorm.csv"), 
            row.names =F, quote = F)

}


# ----------------------------------------------------------
# ----------------------------------------------------------
# ITS-level pca  ---------------------------------------------------------------
# ----------------------------------------------------------
# ----------------------------------------------------------
library(FactoMineR)
library(ggplot2)
library(patchwork)

meta_result <- read.csv("05_immuneSig_ICBresponse/results/meta_results/res_metaRes.csv")
meta_selected <- meta_result[meta_result$Meta_Pval < 0.05, ]
meta_selected$immune_sig = meta_selected$X

result_dir <- "05_immuneSig_ICBresponse/results_r0.4_ITSssgsea_Combat/"

# before Combat: edata_log -------------------------
savepath = paste0(result_dir, "ITSssgsea_beforeCombat17datasets/")
files = list.files(savepath)[-grep('ITSn', list.files(savepath))]
# ssgsea  ----------------------------
ITSfile <- files[-grep('_ssgseanorm',files)]
ITSres <- lapply(as.list(ITSfile), function(filename){
  res <- read.csv(paste0(savepath, filename))
  res_selected <- res[intersect(names(res), c('patientID','label','SampleID',
                     'cancer','sample_type','dataset', meta_selected$immune_sig))]
  return(res_selected)
})


for(i in 1:length(ITSres)){
  if(i == 1){
    its_shared = colnames(ITSres[[1]])
  } else{
    its_shared = intersect(its_shared,colnames(ITSres[[i]]))
  }
}

ITSres_selected <- lapply(ITSres,function(z){ return((z[its_shared]))})
ITSres_selected <- do.call(rbind,ITSres_selected)
rownames(ITSres_selected) <- ITSres_selected$SampleID
ITSres_Annot <- ITSres_selected[c('patientID','label','SampleID',
                     'cancer','sample_type','dataset')]
ITSres_selected2 <- ITSres_selected[-which(is.element(names(ITSres_selected), names(ITSres_Annot)))]                     
pca_fmr = PCA(ITSres_selected2, 
              scale.unit = T, ncp = 6, graph = F)  

## generate PCA using ggplot2
PCs.2dim = pca_fmr$ind$coord[, 1:2]    
colnames(PCs.2dim) = c("Dim1", "Dim2")
xlab = paste("Dim1 (", round(pca_fmr$eig[1,2], 2), "%)", sep = "")
ylab = paste("Dim2 (", round(pca_fmr$eig[2,2], 2), "%)", sep = "")

## plot results of PCA with ggplot2
PCs.2dim = as.data.frame(PCs.2dim)
PCs.2dim$SampleID = rownames(PCs.2dim)
Fig1.PCs = merge(ITSres_Annot, PCs.2dim, by = "SampleID") 

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


p_beforeCombat = p_dataset | p_cancer | p_label
ggsave(paste0(result_dir, "ITS_beforeCombat17datasets_ssgsea_pca.pdf"), p_beforeCombat, width = 15, height = 5)


# ssgsea.norm ----------------------------
ITSfile <- files[grep('_ssgseanorm',files)]
ITSres <- lapply(as.list(ITSfile), function(filename){
  res <- read.csv(paste0(savepath, filename))
  res_selected <- res[intersect(names(res), c('patientID','label','SampleID',
                     'cancer','sample_type','dataset', meta_selected$immune_sig))]
  return(res_selected)
})


for(i in 1:length(ITSres)){
  if(i == 1){
    its_shared = colnames(ITSres[[1]])
  } else{
    its_shared = intersect(its_shared,colnames(ITSres[[i]]))
  }
}

ITSres_selected <- lapply(ITSres,function(z){ return((z[its_shared]))})
ITSres_selected <- do.call(rbind,ITSres_selected)
rownames(ITSres_selected) <- ITSres_selected$SampleID
ITSres_Annot <- ITSres_selected[c('patientID','label','SampleID',
                     'cancer','sample_type','dataset')]
ITSres_selected2 <- ITSres_selected[-which(is.element(names(ITSres_selected), names(ITSres_Annot)))]                     
pca_fmr = PCA(ITSres_selected2, 
              scale.unit = T, ncp = 6, graph = F) 

## generate PCA using ggplot2
PCs.2dim = pca_fmr$ind$coord[, 1:2]  
colnames(PCs.2dim) = c("Dim1", "Dim2")
xlab = paste("Dim1 (", round(pca_fmr$eig[1,2], 2), "%)", sep = "")
ylab = paste("Dim2 (", round(pca_fmr$eig[2,2], 2), "%)", sep = "")

## plot results of PCA with ggplot2
PCs.2dim = as.data.frame(PCs.2dim)
PCs.2dim$SampleID = rownames(PCs.2dim)
Fig1.PCs = merge(ITSres_Annot, PCs.2dim, by = "SampleID") 

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


p_beforeCombat = p_dataset | p_cancer | p_label
ggsave(paste0(result_dir, "ITS_beforeCombat17datasets_ssgseanorm_pca.pdf"), p_beforeCombat, width = 15, height = 5)


# after combat: data_ac -------------------------
savepath = paste0(result_dir, "ITSssgsea_afterCombat17datasets/")
files = list.files(savepath)[-grep('ITSn', list.files(savepath))]
# ssgsea  ----------------------------
ITSfile <- files[-grep('_ssgseanorm',files)]
ITSres <- lapply(as.list(ITSfile), function(filename){
  res <- read.csv(paste0(savepath, filename))
  res_selected <- res[intersect(names(res), c('patientID','label','SampleID',
                     'cancer','sample_type','dataset', meta_selected$immune_sig))]
  return(res_selected)
})


for(i in 1:length(ITSres)){
  if(i == 1){
    its_shared = colnames(ITSres[[1]])
  } else{
    its_shared = intersect(its_shared,colnames(ITSres[[i]]))
  }
}

ITSres_selected <- lapply(ITSres,function(z){ return((z[its_shared]))})
ITSres_selected <- do.call(rbind,ITSres_selected)
rownames(ITSres_selected) <- ITSres_selected$SampleID
ITSres_Annot <- ITSres_selected[c('patientID','label','SampleID',
                     'cancer','sample_type','dataset')]
ITSres_selected2 <- ITSres_selected[-which(is.element(names(ITSres_selected), names(ITSres_Annot)))]                     

save(ITSres_selected, file = paste0(result_dir, "ITS_afterCombat17datasets_ssgsea.Rdata"))

pca_fmr = PCA(ITSres_selected2, 
              scale.unit = T, ncp = 6, graph = F) 
              
## generate PCA using ggplot2
PCs.2dim = pca_fmr$ind$coord[, 1:2]   
colnames(PCs.2dim) = c("Dim1", "Dim2")
xlab = paste("Dim1 (", round(pca_fmr$eig[1,2], 2), "%)", sep = "")
ylab = paste("Dim2 (", round(pca_fmr$eig[2,2], 2), "%)", sep = "")

## plot results of PCA with ggplot2
PCs.2dim = as.data.frame(PCs.2dim)
PCs.2dim$SampleID = rownames(PCs.2dim)
Fig1.PCs = merge(ITSres_Annot, PCs.2dim, by = "SampleID") 

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


p_afterCombat = p_dataset | p_cancer | p_label
ggsave(paste0(result_dir, "ITS_afterCombat17datasets_ssgsea_pca.pdf"), p_afterCombat, width = 15, height = 5)


# ssgsea.norm ----------------------------

ITSfile <- files[grep('_ssgseanorm',files)]

ITSres <- lapply(as.list(ITSfile), function(filename){
  res <- read.csv(paste0(savepath, filename))
  res_selected <- res[intersect(names(res), c('patientID','label','SampleID',
                     'cancer','sample_type','dataset', meta_selected$immune_sig))]
  return(res_selected)
})


for(i in 1:length(ITSres)){
  if(i == 1){
    its_shared = colnames(ITSres[[1]])
  } else{
    its_shared = intersect(its_shared,colnames(ITSres[[i]]))
  }
}


ITSres_selected <- lapply(ITSres,function(z){ return((z[its_shared]))})
ITSres_selected <- do.call(rbind,ITSres_selected)
rownames(ITSres_selected) <- ITSres_selected$SampleID

ITSres_Annot <- ITSres_selected[c('patientID','label','SampleID',
                     'cancer','sample_type','dataset')]
ITSres_selected2 <- ITSres_selected[-which(is.element(names(ITSres_selected), names(ITSres_Annot)))]

save(ITSres_selected, file = paste0(result_dir, "ITS_afterCombat17datasets_ssgseanorm.Rdata"))


# ssgsea.norm ----------------------------

pca_fmr = PCA(ITSres_selected2, 
              scale.unit = T, ncp = 6, graph = F) 

## generate PCA using ggplot2
PCs.2dim = pca_fmr$ind$coord[, 1:2]   
colnames(PCs.2dim) = c("Dim1", "Dim2")
xlab = paste("Dim1 (", round(pca_fmr$eig[1,2], 2), "%)", sep = "")
ylab = paste("Dim2 (", round(pca_fmr$eig[2,2], 2), "%)", sep = "")

## plot results of PCA with ggplot2
PCs.2dim = as.data.frame(PCs.2dim)
PCs.2dim$SampleID = rownames(PCs.2dim)
Fig1.PCs = merge(ITSres_Annot, PCs.2dim, by = "SampleID") 

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


p_afterCombat = p_dataset | p_cancer | p_label
ggsave(paste0(result_dir, "ITS_afterCombat17datasets_ssgseanorm_pca.pdf"), p_afterCombat, width = 15, height = 5)

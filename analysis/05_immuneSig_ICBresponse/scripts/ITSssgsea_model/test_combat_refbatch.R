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
edata <- t(data_bc[-which(is.element(names(data_bc), c("SampleID", "patientID","label","cancer","sample_type","dataset","batch")))])
edata_log <- edata

# ----------------------------------------------------------
# ----------------------------------------------------------
# batch correction with combat  ------------------------------------------------
# ----------------------------------------------------------
# ----------------------------------------------------------

# data_train <- ComBat(dat = edata_log, 
#                  batch = data_bc$batch)
# data_train <- as.data.frame(t(data_train))

data_bc2 <- data_bc[!is.na(data_bc$label),]
edata_log2 <- edata_log[,data_bc2$SampleID]
mod = model.matrix(~as.factor(label), data=data_bc2)
data_train2 <- ComBat(dat = edata_log2, 
                      mod = mod,
                      batch = data_bc2$batch)
data_train2 <- as.data.frame(t(data_train2))

rownames(dataAnnot_all) <- dataAnnot_all$SampleID

library(FactoMineR)
library(ggplot2)
library(patchwork)
# before Combat: edata_log


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


p_label <- 
  ggplot(Fig1.PCs, aes(x=Dim1, y=Dim2, color= label )) + 
  geom_point(size = 1) +
  theme_bw() +
  xlab(xlab) + ylab(ylab) + 
  scale_alpha_discrete(range = c(0.5,1)) +
  coord_fixed() +  
  ggtitle("label")

p_dataset
p_cancer
p_label




data_bc_s <- do.call(rbind,data_bc_all[c(14)])
dataAnnot_all <- do.call(rbind, dataAnnot)
dataAnnot_all$dataset <- str_split_fixed(rownames(dataAnnot_all),'\\.',2)[,1]
data_bc_s <- as.data.frame(data_bc_s)
data_bc_s$SampleID <- row.names(data_bc_s)
data_bc_s <- merge(data_bc_s, dataAnnot_all, by = 'SampleID')
data_bc_s$batch <- data_bc_s$dataset


rownames(data_bc_s) <- data_bc_s$SampleID
edata_s <- t(data_bc_s[-which(is.element(names(data_bc_s), c("SampleID", "patientID","label","cancer","sample_type","dataset","batch")))])

mydata = cbind(t(data_train2),edata_s)

# mydata = cbind(edata_log,edata_s)


data_bc <- do.call(rbind,data_bc_all[c(1:12,14)])
dataAnnot_all <- do.call(rbind, dataAnnot)
dataAnnot_all$dataset <- str_split_fixed(rownames(dataAnnot_all),'\\.',2)[,1]
data_bc <- as.data.frame(data_bc)
data_bc$SampleID <- row.names(data_bc)
data_bc <- merge(data_bc, dataAnnot_all, by = 'SampleID')
data_bc$batch <- data_bc$dataset


rownames(data_bc) <- data_bc$SampleID
data_bc$MLbatch = 'train'
data_bc[data_bc$batch == 'Liu_dataset_CTLA4Naive_outcome', ]$MLbatch = 'test'
data_bc$MLbatch  = factor(data_bc$MLbatch, levels = c('train','test'))
data_bc <- data_bc[colnames(mydata),]

# data_test1 <- ComBat(dat = mydata, 
#                     ref.batch = 'train', 
#                     par.prior=TRUE,
#                     batch = data_bc$MLbatch)
# data_test1[c(1:5), c( "Liu_dataset_CTLA4Naive_outcomePatient112", "Liu_dataset_CTLA4Naive_outcomePatient116")]
# mydata[c(1:5),c( "Liu_dataset_CTLA4Naive_outcomePatient112", "Liu_dataset_CTLA4Naive_outcomePatient116")]
# data_test1 <- as.data.frame(t(data_test1))

data_bc2 <- data_bc[!is.na(data_bc$label),]
mydata2 <- mydata[,data_bc2$SampleID]

mod = model.matrix(~as.factor(label), data=data_bc2)
data_test2 <- ComBat(dat = mydata2, 
                   mod = mod,
                   ref.batch = 'train',
                   batch = data_bc2$MLbatch)
data_test2[c(1:5), c( "Liu_dataset_CTLA4Naive_outcomePatient112", "Liu_dataset_CTLA4Naive_outcomePatient116")]
mydata[c(1:5),c( "Liu_dataset_CTLA4Naive_outcomePatient112", "Liu_dataset_CTLA4Naive_outcomePatient116")]
             
data_test2 <- as.data.frame(t(data_test2))


pca_fmr = PCA(data_test2, 
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
Fig1.PCs = merge(data_bc[c("SampleID", "MLbatch",'dataset','cancer','label')], PCs.2dim, by = "SampleID") 

p_MLbatch <- 
  ggplot(Fig1.PCs, aes(x=Dim1, y=Dim2, color= MLbatch )) + 
  geom_point(size = 1) +
  theme_bw() +
  xlab(xlab) + ylab(ylab) + 
  scale_alpha_discrete(range = c(0.5,1)) +
  coord_fixed() +  
  ggtitle("MLbatch")
p_MLbatch

p_dataset <- 
  ggplot(Fig1.PCs, aes(x=Dim1, y=Dim2, color= dataset )) + 
  geom_point(size = 1) +
  theme_bw() +
  xlab(xlab) + ylab(ylab) + 
  scale_alpha_discrete(range = c(0.5,1)) +
  coord_fixed() +  
  ggtitle("dataset")
p_dataset

p_cancer <- 
  ggplot(Fig1.PCs, aes(x=Dim1, y=Dim2, color= cancer )) + 
  geom_point(size = 1) +
  theme_bw() +
  xlab(xlab) + ylab(ylab) + 
  scale_alpha_discrete(range = c(0.5,1)) +
  coord_fixed() +  
  ggtitle("cancer")
p_cancer

p_label <- 
  ggplot(Fig1.PCs, aes(x=Dim1, y=Dim2, color= label )) + 
  geom_point(size = 1) +
  theme_bw() +
  xlab(xlab) + ylab(ylab) + 
  scale_alpha_discrete(range = c(0.5,1)) +
  coord_fixed() +  
  ggtitle("label")
p_label




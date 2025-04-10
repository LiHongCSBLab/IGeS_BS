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
dir.create("07_plotting_v2/IGeS_EN_shuffle/")

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


# randomize cohorts ------------------------------------------------------------
set.seed(42)

lapply(as.list(1:10), function(rand_index){

  print(paste0("------------------------------------------------------------"))
  print(paste0("------------------------------------------------------------"))
  print(paste0("Randomize for ", rand_index))
  print(paste0("------------------------------------------------------------"))
  print(paste0("------------------------------------------------------------"))
  # set save path ----------------------------------------------------------------
  result_dir <- paste0("07_plotting_v2/IGeS_EN_shuffle/rand", rand_index, "/")
  dir.create(result_dir)

  trainset1=sample(length(data_bc_all), 12, replace = FALSE, prob = NULL)
  valset1=seq(length(data_bc_all))[-trainset1]
  trainset1_name = names(data_bc_all)[trainset1]
  valset1_name = names(data_bc_all)[valset1]
  save(trainset1, trainset1_name, valset1, valset1_name,file = paste0(result_dir, "rand_cohort_index.Rdata"))

  # ------------------------------------------------------------------------------

  data_bc <- do.call(rbind,data_bc_all[c(trainset1)])
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

  # After Combat: data_train2
  pca_fmr = PCA(data_train, 
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



  p_afterCombat = p_dataset |p_cancer| p_label
  ggsave(paste0(result_dir, "Exp_train_afterCombat_pca.pdf"), p_afterCombat, width = 15, height = 5)


  data_val_all = list()

  for(i in valset1){
      # i = 3
      data_bc_s <- do.call(rbind,data_bc_all[i])
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


      data_bc <- do.call(rbind,data_bc_all[c(trainset1,i)])
      dataAnnot_all <- do.call(rbind, dataAnnot)
      dataAnnot_all$dataset <- str_split_fixed(rownames(dataAnnot_all),'\\.',2)[,1]
      data_bc <- as.data.frame(data_bc)
      data_bc$SampleID <- row.names(data_bc)
      data_bc <- merge(data_bc, dataAnnot_all, by = 'SampleID')
      data_bc$batch <- data_bc$dataset


      rownames(data_bc) <- data_bc$SampleID
      data_bc$MLbatch = 'train'
      data_bc[data_bc$batch == unique(data_bc_s$dataset), ]$MLbatch = 'test'
      data_bc$MLbatch  = factor(data_bc$MLbatch, levels = c('train','test'))
      data_bc <- data_bc[colnames(mydata),]

      # data_val1 <- ComBat(dat = mydata, 
      #                     ref.batch = 'train', 
      #                     par.prior=TRUE,
      #                     batch = data_bc$MLbatch)
      # data_val1[c(1:5), c( "Liu_dataset_CTLA4Naive_outcomePatient112", "Liu_dataset_CTLA4Naive_outcomePatient116")]
      # mydata[c(1:5),c( "Liu_dataset_CTLA4Naive_outcomePatient112", "Liu_dataset_CTLA4Naive_outcomePatient116")]
      # data_val1 <- as.data.frame(t(data_val1))

      data_bc2 <- data_bc[!is.na(data_bc$label),]
      mydata2 <- mydata[,data_bc2$SampleID]

      mod = model.matrix(~as.factor(label), data=data_bc2)
      data_val2 <- ComBat(dat = mydata2, 
                      # mod = mod,
                      ref.batch = 'train', 
                      batch = data_bc2$MLbatch)
      # data_val2[c(1:5), head(data_bc_s$SampleID)]
      # mydata[c(1:5),head(data_bc_s$SampleID)]
                  
      data_val2 <- as.data.frame(t(data_val2))
  
      data_val_all[[i]] <- data_val2[data_bc_s$SampleID,]
      # save(data_val2, file = paste0(result_dir, "Exp_valds_",i,"_afterCombat.Rdata"))

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
      Fig1.PCs = merge(data_bc[c("SampleID", "MLbatch",'dataset','cancer','label')], PCs.2dim, by = "SampleID") 

      p_MLbatch <- 
      ggplot(Fig1.PCs, aes(x=Dim1, y=Dim2, color= MLbatch )) + 
      geom_point(size = 1) +
      theme_bw() +
      xlab(xlab) + ylab(ylab) + 
      scale_alpha_discrete(range = c(0.5,1)) +
      coord_fixed() +  
      ggtitle("MLbatch")
      # p_MLbatch

      p_dataset <- 
      ggplot(Fig1.PCs, aes(x=Dim1, y=Dim2, color= dataset )) + 
      geom_point(size = 1) +
      theme_bw() +
      xlab(xlab) + ylab(ylab) + 
      scale_alpha_discrete(range = c(0.5,1)) +
      coord_fixed() +  
      ggtitle("dataset")
      # p_dataset

      p_cancer <- 
      ggplot(Fig1.PCs, aes(x=Dim1, y=Dim2, color= cancer )) + 
      geom_point(size = 1) +
      theme_bw() +
      xlab(xlab) + ylab(ylab) + 
      scale_alpha_discrete(range = c(0.5,1)) +
      coord_fixed() +  
      ggtitle("cancer")
      # p_cancer

      p_label <- 
      ggplot(Fig1.PCs, aes(x=Dim1, y=Dim2, color= label )) + 
      geom_point(size = 1) +
      theme_bw() +
      xlab(xlab) + ylab(ylab) + 
      scale_alpha_discrete(range = c(0.5,1)) +
      coord_fixed() +  
      ggtitle("label")
      # p_label
      p_afterCombat =(p_MLbatch|p_label )/ (p_dataset | p_cancer)

      ggsave(paste0(result_dir, "Exp_valds_",i,"_afterCombat_pca.pdf"), p_afterCombat, width = 10, height = 5)
  }

  data_val_all_df <- do.call(rbind,data_val_all[valset1])
  data_ac <- rbind(data_train2, data_val_all_df)
  save(data_ac, file = paste0(result_dir, "Exp_afterCombatRef.Rdata"))



  # ----------------------------------------------------------
  # ----------------------------------------------------------
  # ITSssgsea  -------------------------------------------------------------------
  # ----------------------------------------------------------
  # ----------------------------------------------------------
  source("05_immuneSig_ICBresponse/scripts/functions/05_immunotherapy_immuneSig_analysis_function.R")
  source("05_immuneSig_ICBresponse/scripts/functions/predictor_ITS.R")
  source("06_results_summaryandplotting/scripts/immunesig_selection_function.R")

  savepath_allresult = "02_tumorGene_immuneSig_correlation/results_selected/ITS_oe_original_spearman_TUMERIC_r0.4/"
  # after combat: data_ac
  savepath = paste0(result_dir, "ITSssgsea_afterCombat17datasets/")
  dir.create(savepath)
  dataAnnot_all <- dataAnnot_all[!is.na(dataAnnot_all$label), ]
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


    # ITSn <- immune_sigAnalysis(data = df, 
    #                           gs = genelist_n,
    #                           method = enrich_method, 
    #                           kcdf = "Gaussian",
    #                           ssgsea.norm = FALSE) 
    # ITSn_res <- as.data.frame(t(ITSn))
    # ITSn_res$SampleID <- rownames(ITSn_res)
    # ITSn_res <- inner_join(ds, ITSn_res,by='SampleID')
    # write.csv(ITSn_res, paste0(savepath, "/", dataset, "ITSn_", enrich_method, ".csv"), 
    #           row.names =F, quote = F)


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


    # ITSn <- immune_sigAnalysis(data = df, 
    #                           gs = genelist_n,
    #                           method = enrich_method, 
    #                           kcdf = "Gaussian",
    #                           ssgsea.norm = TRUE) 
    # ITSn_res <- as.data.frame(t(ITSn))
    # ITSn_res$SampleID <- rownames(ITSn_res)
    # ITSn_res <- inner_join(ds, ITSn_res,by='SampleID')
    # write.csv(ITSn_res, paste0(savepath, "/", dataset, "ITSn_", enrich_method, "_ssgseanorm.csv"), 
    #           row.names =F, quote = F)

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


  # after combat: data_ac -------------------------
  savepath = paste0(result_dir, "ITSssgsea_afterCombat17datasets/")
  files = list.files(savepath)
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

  ITSres_selected3 = ITSres_selected2[rownames(data_train2),]
  pca_fmr = PCA(ITSres_selected3, 
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

})
# rm(list=ls())
work_path <- "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/"

setwd(work_path)
options(stringsAsFactors = F)

library(reshape2)
library(parallel)
library(patchwork)
library(dplyr)
library(stringr)
library(data.table)

library(ggplot2)
library(ggplotify)
library(grid)
library(cowplot)
require(GGally)
library(plot3D)
library(fmsb)
library(readxl)

library(dplyr)
library(UpSetR)
library(vcd)
library(fmsb)
library(pheatmap)
# source("06_results_summaryandplotting/scripts/OriImmuSig_analysis/OriImmuSig_gene_stat_function.R")
# source("06_results_summaryandplotting/scripts/immunesig_selection_function.R")

source("03_drug_immuneSig_enrichment/scripts/02_GSEA_geneset_merge_metaAnalysis_function.R")
# source("07_plotting_v2/scripts/functions/OriImmuSig_function_detection.R")
source("03_drug_immuneSig_enrichment/scripts/02_GSEA_geneset_preparation_function.R")
source("03_drug_immuneSig_enrichment/scripts/02_GSEA_geneset_merge_metaAnalysis_function.R")
source("06_results_summaryandplotting/scripts/immunesig_selection_function.R")
source("06_results_summaryandplotting/scripts/OriImmuSig_analysis/OriImmuSig_function_detection.R")
source("06_results_summaryandplotting/scripts/OriImmuSig_analysis/OriImmuSigexpr_purity_detection.R")
source("06_results_summaryandplotting/scripts/OriImmuSig_analysis/OriImmuSigexpr_purity_detection.R")
source("06_results_summaryandplotting/scripts/OriImmuSig_analysis/immuneModulatorDetection.R")


# ----------------------------------------------------------------------------
dir.create("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/OriImmuSig_similarity/")
# ------------------------------------------------------------------------------
immuneSigInfoPath = "06_results_summaryandplotting/data/immene_cell_hieracy/"

immunesigInfo <- as.data.frame(readxl::read_excel(paste0(immuneSigInfoPath, "immune_sig_hieracy_new.xlsx"), sheet = 1))
immunesigInfo <- immunesigInfo[c(1,4,5:9)]
immunesigInfo$immune_sig = tolower(immunesigInfo$immune_sig)
immunesigInfo = immunesig_name_unify(immunesigInfo)
rownames(immunesigInfo) <- immunesigInfo$immune_sig

meta_result <- read.csv("05_immuneSig_ICBresponse/results/meta_results/res_metaRes.csv")
meta_selected <- meta_result[meta_result$Meta_Pval < 0.05, ]
meta_selected$immune_sig = meta_selected$X

load("05_immuneSig_ICBresponse/results_r0.4_ITSssgsea_CombatRef/finalmodel_17datasets_ssgseanorm/metaAnalysis_filter/en/elasticnet_model_pred.Rdata")

en_model_para = m_en$learner.model$glmnet.fit
m_en_imp <- data.frame(en_model_para$beta[,which(en_model_para$lambda == m_en$learner.model$lambda.min)])
names(m_en_imp) <- "weight"
m_en_imp$immune_sig = row.names(m_en_imp)
meta_selected = merge(meta_selected, m_en_imp, by = "immune_sig")
meta_selected = meta_selected[meta_selected$weight != 0, ]


# LOAD original immune signature gene sets ---------------------------------
cancerlist <- c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA", 
                "GBM", "HNSC", "KIHC", "KIRC", "KIRP", "LGG", "LIHC", "LUAD", 
                "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", 
                "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM")

OriImmuSig_set_all <- lapply(as.list(cancerlist), function(cancer){
  # cancer = cancerlist
  immuneSigGSpath = "02_tumorGene_immuneSig_correlation/data/immune_signatures/originalSignaturesGenes/"
  immunesigGS1 <- as.data.frame(read_excel(paste0(immuneSigGSpath, "immunesig_original_used.xlsx"), sheet = 1))
  immunesigGS1 <- immunesigGS1[-grep("various", immunesigGS1$Genes), ]
  immunesigGS1 <- immunesigGS1[-grep("SignatureMatrix", immunesigGS1$Genes), ]
  immunesigGS1 <- immunesigGS1[-grep("ReferenceCellExpression", immunesigGS1$Genes), ]
  immunesigGS1 <- immunesigGS1[-grep("mergedSignatures", immunesigGS1$Genes), ]
  immunesigGS1 <- immunesigGS1[-grep("L22Matrix", immunesigGS1$Genes), ]
  
  immunesigGS1 <- immunesigGS1[-c(2:5)]
  immunesigGS1Name <- as.list(immunesigGS1$immune_sig)
  immunesigGS1 <- lapply(immunesigGS1Name, function(x){
    y=immunesigGS1[immunesigGS1$immune_sig == x,][-1]
    unlist(y[-which(is.na(y))])
  })
  names(immunesigGS1) = unlist(immunesigGS1Name)
  
  immunesigGS2 <- read_excel(paste0(immuneSigGSpath, "immunesig_original_used.xlsx"), sheet = 2)
  immunesigGS2 <- as.data.frame(immunesigGS2[immunesigGS2$CancerType == cancer, ])
  immunesigGS2Name <- as.list(immunesigGS2$Signatures)
  immunesigGS2 <- lapply(immunesigGS2Name, function(x){
    y=immunesigGS2[immunesigGS2$Signatures == x,][-1]
    unlist(y[-which(is.na(y))])
  })
  
  names(immunesigGS2) = paste0(unlist(immunesigGS2Name),"_ConsensusTMEgsva")
  immunesigGS = c(immunesigGS1, immunesigGS2)
  immunesigGSname = data.frame(signame = names(immunesigGS), immune_sig = tolower(names(immunesigGS)))
  immunesigGSname = immunesig_name_unify(immunesigGSname)
  names(immunesigGS) = immunesigGSname$immune_sig
  return(immunesigGS)
})

names(OriImmuSig_set_all) <- cancerlist
save(OriImmuSig_set_all, file = "07_plotting_v2/03_ITScharacteristic_metaPlusFilter/OriImmuSig_similarity/OriImmuSig_set.Rdata")

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# 1) gene set similarity
# ------------------------------------------------------------------------------
#   jaccard_index() "jaccard_index"
dir.create("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/OriImmuSig_similarity/jcindex/")

OriImmuSig_set_all_jcindex <- lapply(as.list(names(OriImmuSig_set_all)), function(cancer){
  # cancer = names(OriImmuSig_set_all)[[1]]
  OriImmuSig_set <- OriImmuSig_set_all[[cancer]]
  OriImmuSig_set <- OriImmuSig_set[intersect(names(OriImmuSig_set),meta_selected$X)]
  OriImmuSig_jcindex <- list()
  for(i in 1:length(OriImmuSig_set)){
    x = OriImmuSig_set[[i]]
    OriImmuSig_jcindex[[i]] = sapply(OriImmuSig_set, function(y){
      # jaccard_index(x,y)
      gs_similarity(x, y, similarity_method = "jaccard_index")
      # length(intersect(x, y))
    })
  }
  names(OriImmuSig_jcindex) = names(OriImmuSig_set)
  OriImmuSig_jcindex = do.call(cbind, OriImmuSig_jcindex)

  pdf(paste0("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/OriImmuSig_similarity/jcindex/OriImmuSig_",cancer,".pdf"), width = 20, height = 20)
  tmp = pheatmap(OriImmuSig_jcindex,
                 annotation_col = immunesigInfo['Type'],
                 annotation_row = immunesigInfo['Type'])
  dev.off()
  
  tree = tmp$tree_col
  treedf = as.data.frame(as.matrix(do.call(cbind, tree[c(4,2:3)])))
  treedf$height = as.numeric(treedf$height)
  treedf$order = as.numeric(treedf$order)
  names(treedf) = c("immune_sig", "height","order")
  treedf$immune_sig = as.character(treedf$immune_sig)
  v = cutree(tree, k = 15)[tree$order] # , h = 1.5
  # gaps = which((v[-1] - v[-length(v)]) != 0)
  OriImmuSig_cls <- as.data.frame(as.matrix(v))
  # OriImmuSig_cls$immune_sig <- str_split(rownames(OriImmuSig_cls),'[.]',simplify = T)[,1]
  OriImmuSig_cls$immune_sig <- rownames(OriImmuSig_cls)
  OriImmuSig_cls = inner_join(treedf, OriImmuSig_cls)
  OriImmuSig_cls = OriImmuSig_cls[order(OriImmuSig_cls$V1),]
  OriImmuSig_cls$V1 = as.character(OriImmuSig_cls$V1)
  
  immunesigInfo2 <- merge(immunesigInfo, OriImmuSig_cls, by = 'immune_sig', all.y = T)
  rownames(immunesigInfo2) <- immunesigInfo2$immune_sig
  immunesigInfo2 <- immunesigInfo2[order(immunesigInfo2$V1),]
  OriImmuSig_jcindex3 <- OriImmuSig_jcindex
  OriImmuSig_jcindex3[OriImmuSig_jcindex3>0.8] <- 1
  OriImmuSig_jcindex3[OriImmuSig_jcindex3>0.6 & OriImmuSig_jcindex3<=0.8] <- 0.8
  OriImmuSig_jcindex3[OriImmuSig_jcindex3>0.4 & OriImmuSig_jcindex3<=0.6] <- 0.6
  OriImmuSig_jcindex3[OriImmuSig_jcindex3>0.2 & OriImmuSig_jcindex3<=0.4] <- 0.4
  OriImmuSig_jcindex3[OriImmuSig_jcindex3>0 & OriImmuSig_jcindex3<=0.2] <- 0.2
  OriImmuSig_jcindex3[OriImmuSig_jcindex3<=0] <- 0
  OriImmuSig_jcindex3 <- OriImmuSig_jcindex3[rev(OriImmuSig_cls$immune_sig), rev(OriImmuSig_cls$immune_sig)]
  
  # immunesigInfo <- immunesigInfo[order(immunesigInfo$Type), ]
  # OriImmuSig_jcindex3 <- OriImmuSig_jcindex3[intersect(immunesigInfo$immune_sig, rownames(OriImmuSig_jcindex3)),
  #                          intersect(immunesigInfo$immune_sig, colnames(OriImmuSig_jcindex3))]
  
  pdf(paste0("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/OriImmuSig_similarity/jcindex/OriImmuSig_",cancer,"_bins.pdf"), width = 20, height = 20)
  
  pheatmap(OriImmuSig_jcindex3,
           cluster_cols = F, 
           cluster_rows = F, 
           scale = "none",
           color = colorRampPalette(colors = c("blue","white","yellow","red"))(6),
           annotation_col = immunesigInfo2[c('Type','V1')], # 'SubType',
           annotation_row = immunesigInfo2[c('Type','V1')]) # 'SubType',
  dev.off()
  
  return(OriImmuSig_jcindex)
})

names(OriImmuSig_set_all_jcindex) <- names(OriImmuSig_set_all)





#   Otsuka_Ochiai_coefficient() "Otsuka_Ochiai"
dir.create("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/OriImmuSig_similarity/Ochiai/")
OriImmuSig_set_all_Ochiai <-  lapply(as.list(names(OriImmuSig_set_all)), function(cancer){
  # cancer = names(OriImmuSig_set_all)[[1]]
  OriImmuSig_set <- OriImmuSig_set_all[[cancer]]
  OriImmuSig_set <- OriImmuSig_set[intersect(names(OriImmuSig_set),meta_selected$X)]
  OriImmuSig_Ochiai <- list()

  for(i in 1:length(OriImmuSig_set)){
    x = OriImmuSig_set[[i]]
    OriImmuSig_Ochiai[[i]] = sapply(OriImmuSig_set, function(y){
      # jaccard_index(x,y)
      gs_similarity(x, y, similarity_method = "Otsuka_Ochiai")
      # length(intersect(x, y))
    })
  }
  names(OriImmuSig_Ochiai) = names(OriImmuSig_set)
  OriImmuSig_Ochiai = do.call(cbind, OriImmuSig_Ochiai)
  OriImmuSig_Ochiai[is.na(OriImmuSig_Ochiai)] = 0
  
  pdf(paste0("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/OriImmuSig_similarity/Ochiai/OriImmuSig_",cancer,".pdf"), width = 20, height = 20)
  tmp = pheatmap(OriImmuSig_Ochiai,
                 annotation_col = immunesigInfo['Type'],
                 annotation_row = immunesigInfo['Type'])
  dev.off()
  
  
  tree = tmp$tree_col
  treedf = as.data.frame(as.matrix(do.call(cbind, tree[c(4,2:3)])))
  treedf$height = as.numeric(treedf$height)
  treedf$order = as.numeric(treedf$order)
  names(treedf) = c("immune_sig", "height","order")
  treedf$immune_sig = as.character(treedf$immune_sig)
  v = cutree(tree, k = 15)[tree$order] # , h = 1.5
  # gaps = which((v[-1] - v[-length(v)]) != 0)
  OriImmuSig_cls <- as.data.frame(as.matrix(v))
  # OriImmuSig_cls$immune_sig <- str_split(rownames(OriImmuSig_cls),'[.]',simplify = T)[,1]
  OriImmuSig_cls$immune_sig <- rownames(OriImmuSig_cls)
  OriImmuSig_cls = inner_join(treedf, OriImmuSig_cls)
  OriImmuSig_cls = OriImmuSig_cls[order(OriImmuSig_cls$V1),]
  OriImmuSig_cls$V1 = as.character(OriImmuSig_cls$V1)
  
  immunesigInfo2 <- merge(immunesigInfo, OriImmuSig_cls, by = 'immune_sig', all.y = T)
  rownames(immunesigInfo2) <- immunesigInfo2$immune_sig
  immunesigInfo2 <- immunesigInfo2[order(immunesigInfo2$V1),]
  
  OriImmuSig_Ochiai3 <- OriImmuSig_Ochiai
  OriImmuSig_Ochiai3[OriImmuSig_Ochiai3>0.8] <- 1
  OriImmuSig_Ochiai3[OriImmuSig_Ochiai3>0.6 & OriImmuSig_Ochiai3<=0.8] <- 0.8
  OriImmuSig_Ochiai3[OriImmuSig_Ochiai3>0.4 & OriImmuSig_Ochiai3<=0.6] <- 0.6
  OriImmuSig_Ochiai3[OriImmuSig_Ochiai3>0.2 & OriImmuSig_Ochiai3<=0.4] <- 0.4
  OriImmuSig_Ochiai3[OriImmuSig_Ochiai3>0 & OriImmuSig_Ochiai3<=0.2] <- 0.2
  OriImmuSig_Ochiai3[OriImmuSig_Ochiai3 <= 0 ] <- 0
  OriImmuSig_Ochiai3 <- OriImmuSig_Ochiai3[OriImmuSig_cls$immune_sig, OriImmuSig_cls$immune_sig]
  
  # immunesigInfo <- immunesigInfo[order(immunesigInfo$Type), ]
  # OriImmuSig_Ochiai3 <- OriImmuSig_Ochiai3[intersect(immunesigInfo$immune_sig, rownames(OriImmuSig_Ochiai3)),
  #                          intersect(immunesigInfo$immune_sig, colnames(OriImmuSig_Ochiai3))]
  
  pdf(paste0("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/OriImmuSig_similarity/Ochiai/OriImmuSig_",cancer,"_bins.pdf"), width = 20, height = 20)
  
  pheatmap(OriImmuSig_Ochiai3,
           cluster_cols = F, 
           cluster_rows = F, 
           scale = "none",
           color = colorRampPalette(colors = c("blue","white","yellow","red"))(6),
           annotation_col = immunesigInfo2[c('Type', 'V1')], # 'SubType',
           annotation_row = immunesigInfo2[c('Type','V1')])  # 'SubType',
  dev.off()
  
  
  return(OriImmuSig_Ochiai)
})
names(OriImmuSig_set_all_Ochiai) <- names(OriImmuSig_set_all)

save(OriImmuSig_set_all_jcindex, OriImmuSig_set_all_Ochiai, file = "07_plotting_v2/03_ITScharacteristic_metaPlusFilter/OriImmuSig_similarity/OriImmuSig_gs_similarity.Rdata")

# ------------------------------------------------------------------------------
# 2) kappa statistics 
# ------------------------------------------------------------------------------
dir.create("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/OriImmuSig_similarity/kappa/")
OriImmuSig_set_all_kappa <- lapply(as.list(names(OriImmuSig_set_all)), function(cancer){
  # cancer = names(OriImmuSig_set_all)[[1]]
  OriImmuSig_set <- OriImmuSig_set_all[[cancer]]
  OriImmuSig_set <- OriImmuSig_set[intersect(names(OriImmuSig_set),meta_selected$X)]

  OriImmuSig_mat <- matrix(0, length(OriImmuSig_set), length(unique(unlist(OriImmuSig_set))))
  rownames(OriImmuSig_mat) <- names(OriImmuSig_set)
  colnames(OriImmuSig_mat) <- unique(unlist(OriImmuSig_set))
  
  for(sig in names(OriImmuSig_set)){
    # sig=names(OriImmuSig_set)[1]
    x = OriImmuSig_set[[sig]]
    for(g in x){
      OriImmuSig_mat[which(rownames(OriImmuSig_mat) == sig), which(colnames(OriImmuSig_mat) == g)] = 1
    }
  }
  
  OriImmuSig_kappa <- list()
  
  for(i in 1:nrow(OriImmuSig_mat)){
    # i = 1
    OriImmuSig_kappa_res <- list()
    x = OriImmuSig_mat[i,]
    for(j in 1:nrow(OriImmuSig_mat)){
      
      y = OriImmuSig_mat[j,]
      df = as.data.frame(cbind(x,y))
      df_a = df[df$x == 1 & df$y == 1, ]
      df_b = df[df$x == 0 & df$y == 1, ]
      df_c = df[df$x == 1 & df$y == 0, ]
      df_d = df[df$x == 0 & df$y == 0, ]
      
      a = nrow(df_a)
      b = nrow(df_b)
      c = nrow(df_c)
      d = nrow(df_d)
      n <- a + b + c + d
      p0 <- (a + d) / n
      pe <- ((a + b) * (a + c) + (c + d) * (b + d)) / n^2
      kappa <- (p0 - pe) / (1 - pe)
      
      sd <- sqrt(p0 * (1 - p0) / (1 - pe)^2)
      se <- sd / sqrt(n)
      z_alpha <- qnorm(0.025, lower.tail = F)
      CI <- c(kappa - z_alpha * se, kappa + z_alpha * se)
      kappa0 <- 0.7
      z <- (kappa - kappa0) / se
      pvalue <- 2 * pnorm(abs(z), lower.tail = F)
      
      df <- data.frame(i = rownames(OriImmuSig_mat)[i],
                       j = rownames(OriImmuSig_mat)[j], 
                       kappa = kappa, 
                       sd = sd, 
                       se = se,
                       CI_lower = CI[1],
                       CI_upper = CI[2],
                       z = z,
                       pvalue = pvalue)
      # K <- vcd::Kappa(matrix(c(a, b, c, d), byrow = T, nrow = 2))
      
      OriImmuSig_kappa_res[[j]] = df
    }
    OriImmuSig_kappa[[i]] = do.call(rbind, OriImmuSig_kappa_res)
  }
  
  OriImmuSig_kappa <- do.call(rbind, OriImmuSig_kappa)

  OriImmuSig_kappa2 <- dcast(OriImmuSig_kappa[c('i','j','kappa')],i~j)
  rownames(OriImmuSig_kappa2) <- OriImmuSig_kappa2$i
  OriImmuSig_kappa2 <- OriImmuSig_kappa2[-1]
  
  
 
  pdf(paste0("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/OriImmuSig_similarity/kappa/OriImmuSig_",cancer,".pdf"), width = 20, height = 20)
  tmp = pheatmap(OriImmuSig_kappa2,
                 annotation_col = immunesigInfo['Type'],
                 annotation_row = immunesigInfo['Type'])
  dev.off()
  
  
  tree = tmp$tree_col
  treedf = as.data.frame(as.matrix(do.call(cbind, tree[c(4,2:3)])))
  treedf$height = as.numeric(treedf$height)
  treedf$order = as.numeric(treedf$order)
  names(treedf) = c("immune_sig", "height","order")
  treedf$immune_sig = as.character(treedf$immune_sig)
  v = cutree(tree, k = 15)[tree$order] # , h = 1.5
  # gaps = which((v[-1] - v[-length(v)]) != 0)
  OriImmuSig_cls <- as.data.frame(as.matrix(v))
  # OriImmuSig_cls$immune_sig <- str_split(rownames(OriImmuSig_cls),'[.]',simplify = T)[,1]
  OriImmuSig_cls$immune_sig <- rownames(OriImmuSig_cls)
  OriImmuSig_cls = inner_join(treedf, OriImmuSig_cls)
  OriImmuSig_cls = OriImmuSig_cls[order(OriImmuSig_cls$V1),]
  OriImmuSig_cls$V1 = as.character(OriImmuSig_cls$V1)
  
  immunesigInfo2 <- merge(immunesigInfo, OriImmuSig_cls, by = 'immune_sig', all.y = T)
  rownames(immunesigInfo2) <- immunesigInfo2$immune_sig
  immunesigInfo2 <- immunesigInfo2[order(immunesigInfo2$V1),]
  
  OriImmuSig_kappa3 <- OriImmuSig_kappa2
  OriImmuSig_kappa3[OriImmuSig_kappa3>0.8] <- 1
  OriImmuSig_kappa3[OriImmuSig_kappa3>0.6 & OriImmuSig_kappa3<=0.8] <- 0.8
  OriImmuSig_kappa3[OriImmuSig_kappa3>0.4 & OriImmuSig_kappa3<=0.6] <- 0.6
  OriImmuSig_kappa3[OriImmuSig_kappa3>0.2 & OriImmuSig_kappa3<=0.4] <- 0.4
  OriImmuSig_kappa3[OriImmuSig_kappa3>0 & OriImmuSig_kappa3<=0.2] <- 0.2
  OriImmuSig_kappa3[OriImmuSig_kappa3<=0] <- 0
  OriImmuSig_kappa3 <- OriImmuSig_kappa3[OriImmuSig_cls$immune_sig, OriImmuSig_cls$immune_sig]
  
  # immunesigInfo <- immunesigInfo[order(immunesigInfo$Type), ]
  # OriImmuSig_kappa3 <- OriImmuSig_kappa3[intersect(immunesigInfo$immune_sig, rownames(OriImmuSig_kappa3)),
  #                          intersect(immunesigInfo$immune_sig, colnames(OriImmuSig_kappa3))]
  
  pdf(paste0("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/OriImmuSig_similarity/kappa/OriImmuSig_",cancer,"_bins.pdf"), width = 20, height = 20)
  
  pheatmap(OriImmuSig_kappa3,
           cluster_cols = F, 
           cluster_rows = F, 
           scale = "none",
           color = colorRampPalette(colors = c("blue","white","yellow","red"))(6),
           annotation_col = immunesigInfo2[c('Type','V1')], # 'SubType',
           annotation_row = immunesigInfo2[c('Type','V1')]) # 'SubType',
  dev.off()
  
  
  OriImmuSig_kappa$cancer = cancer
  return(OriImmuSig_kappa)
})

names(OriImmuSig_set_all_kappa) <- names(OriImmuSig_set_all)

save(OriImmuSig_set_all_kappa, file = "07_plotting_v2/03_ITScharacteristic_metaPlusFilter/OriImmuSig_similarity/OriImmuSig_gs_kappa.Rdata")



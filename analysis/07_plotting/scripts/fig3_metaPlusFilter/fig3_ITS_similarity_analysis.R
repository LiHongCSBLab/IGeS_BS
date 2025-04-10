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
# source("06_results_summaryandplotting/scripts/ITS_analysis/ITS_gene_stat_function.R")
# source("06_results_summaryandplotting/scripts/immunesig_selection_function.R")

source("03_drug_immuneSig_enrichment/scripts/02_GSEA_geneset_merge_metaAnalysis_function.R")
source("07_plotting_v2/scripts/functions/ITS_function_detection.R")

source("03_drug_immuneSig_enrichment/scripts/02_GSEA_geneset_preparation_function.R")
source("03_drug_immuneSig_enrichment/scripts/02_GSEA_geneset_merge_metaAnalysis_function.R")

source("06_results_summaryandplotting/scripts/immunesig_selection_function.R")
source("06_results_summaryandplotting/scripts/ITS_analysis/ITS_function_detection.R")
source("06_results_summaryandplotting/scripts/ITS_analysis/ITSexpr_purity_detection.R")
source("06_results_summaryandplotting/scripts/ITS_analysis/ITSexpr_purity_detection.R")
source("06_results_summaryandplotting/scripts/ITS_analysis/immuneModulatorDetection.R")


# ----------------------------------------------------------------------------

dir.create("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSsimilarity/")


# ------------------------------------------------------------------------------
immuneSigInfoPath = "06_results_summaryandplotting/data/immene_cell_hieracy/"

immunesigInfo <- as.data.frame(readxl::read_excel(paste0(immuneSigInfoPath, "immune_sig_hieracy_new.xlsx"), sheet = 1))
immunesigInfo <- immunesigInfo[c(1,4,5:9)]
immunesigInfo$immune_sig = tolower(immunesigInfo$immune_sig)
immunesigInfo = immunesig_name_unify(immunesigInfo)
rownames(immunesigInfo) <- immunesigInfo$immune_sig

# ------------------------------------------------------------------------------
load("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITS_set.Rdata")

# ------------------------------------------------------------------------------
# 1) gene set similarity
#   jaccard_index() "jaccard_index"
#   Otsuka_Ochiai_coefficient() "Otsuka_Ochiai"
# 
# 2) kappa statistics 
#
# 3) function similarity 
#   GOSemSim R package
#       
#
# 
# ------------------------------------------------------------------------------
# 1) gene set similarity
# ------------------------------------------------------------------------------
#   jaccard_index() "jaccard_index"
dir.create("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSsimilarity/jcindex/")

ITS_p_set_jcindex <- lapply(as.list(names(ITS_p_set)), function(cancer_type){
  # cancer_type = names(ITS_p_set)[1]
  ITS_set = ITS_p_set[[cancer_type]]
  # ITS_set = ITS_p_set[[1]]
  ITS_jcindex <- list()
  for(i in 1:length(ITS_set)){
    x = ITS_set[[i]]
    ITS_jcindex[[i]] = sapply(ITS_set, function(y){
      # jaccard_index(x,y)
      gs_similarity(x, y, similarity_method = "jaccard_index")
      # length(intersect(x, y))
    })
  }
  names(ITS_jcindex) = names(ITS_set)
  ITS_jcindex = do.call(cbind, ITS_jcindex)
  
  
  pdf(paste0("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSsimilarity/jcindex/ITS_p_",cancer_type,".pdf"), width = 20, height = 20)
  tmp = pheatmap(ITS_jcindex,
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
  its_cls <- as.data.frame(as.matrix(v))
  # its_cls$immune_sig <- str_split(rownames(its_cls),'[.]',simplify = T)[,1]
  its_cls$immune_sig <- rownames(its_cls)
  its_cls = inner_join(treedf, its_cls)
  its_cls = its_cls[order(its_cls$V1),]
  its_cls$V1 = as.character(its_cls$V1)
  
  immunesigInfo2 <- merge(immunesigInfo, its_cls, by = 'immune_sig', all.y = T)
  rownames(immunesigInfo2) <- immunesigInfo2$immune_sig
  immunesigInfo2 <- immunesigInfo2[order(immunesigInfo2$V1),]
  ITS_jcindex3 <- ITS_jcindex
  ITS_jcindex3[ITS_jcindex3>0.8] <- 1
  ITS_jcindex3[ITS_jcindex3>0.6 & ITS_jcindex3<=0.8] <- 0.8
  ITS_jcindex3[ITS_jcindex3>0.4 & ITS_jcindex3<=0.6] <- 0.6
  ITS_jcindex3[ITS_jcindex3<=0.4] <- 0.4
  ITS_jcindex3 <- ITS_jcindex3[rev(its_cls$immune_sig), rev(its_cls$immune_sig)]
  
  # immunesigInfo <- immunesigInfo[order(immunesigInfo$Type), ]
  # ITS_jcindex3 <- ITS_jcindex3[intersect(immunesigInfo$immune_sig, rownames(ITS_jcindex3)),
  #                          intersect(immunesigInfo$immune_sig, colnames(ITS_jcindex3))]
  
  pdf(paste0("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSsimilarity/jcindex/ITS_p_",cancer_type,"_bins.pdf"), width = 20, height = 20)
  
  pheatmap(ITS_jcindex3,
           cluster_cols = F, 
           cluster_rows = F, 
           scale = "none",
           color = colorRampPalette(colors = c("blue","white","yellow","red"))(5),
           annotation_col = immunesigInfo2[c('Type','V1')], # 'SubType',
           annotation_row = immunesigInfo2[c('Type','V1')]) # 'SubType',
  dev.off()
  
  return(ITS_jcindex)
})
names(ITS_p_set_jcindex) <- names(ITS_p_set)


# ------------------------------------------------------------------------------
#   Otsuka_Ochiai_coefficient() "Otsuka_Ochiai"
dir.create("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSsimilarity/Ochiai/")

ITS_p_set_Ochiai <- lapply(as.list(names(ITS_p_set)), function(cancer_type){
  # cancer_type = names(ITS_p_set)[1]
  ITS_set = ITS_p_set[[cancer_type]]
  ITS_Ochiai <- list()
  for(i in 1:length(ITS_set)){
    x = ITS_set[[i]]
    ITS_Ochiai[[i]] = sapply(ITS_set, function(y){
      # jaccard_index(x,y)
      gs_similarity(x, y, similarity_method = "Otsuka_Ochiai")
      # length(intersect(x, y))
    })
  }
  names(ITS_Ochiai) = names(ITS_set)
  ITS_Ochiai = do.call(cbind, ITS_Ochiai)
  
  ITS_Ochiai[is.na(ITS_Ochiai)] = 0
  
  pdf(paste0("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSsimilarity/Ochiai/ITS_p_",cancer_type,".pdf"), width = 20, height = 20)
  tmp = pheatmap(ITS_Ochiai,
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
  its_cls <- as.data.frame(as.matrix(v))
  # its_cls$immune_sig <- str_split(rownames(its_cls),'[.]',simplify = T)[,1]
  its_cls$immune_sig <- rownames(its_cls)
  its_cls = inner_join(treedf, its_cls)
  its_cls = its_cls[order(its_cls$V1),]
  its_cls$V1 = as.character(its_cls$V1)
  
  immunesigInfo2 <- merge(immunesigInfo, its_cls, by = 'immune_sig', all.y = T)
  rownames(immunesigInfo2) <- immunesigInfo2$immune_sig
  immunesigInfo2 <- immunesigInfo2[order(immunesigInfo2$V1),]

  ITS_Ochiai3 <- ITS_Ochiai
  ITS_Ochiai3[ITS_Ochiai3>0.8] <- 1
  ITS_Ochiai3[ITS_Ochiai3>0.6 & ITS_Ochiai3<=0.8] <- 0.8
  ITS_Ochiai3[ITS_Ochiai3>0.4 & ITS_Ochiai3<=0.6] <- 0.6
  ITS_Ochiai3[ITS_Ochiai3<=0.4] <- 0.4
  ITS_Ochiai3 <- ITS_Ochiai3[its_cls$immune_sig, its_cls$immune_sig]
  
  # immunesigInfo <- immunesigInfo[order(immunesigInfo$Type), ]
  # ITS_Ochiai3 <- ITS_Ochiai3[intersect(immunesigInfo$immune_sig, rownames(ITS_Ochiai3)),
  #                          intersect(immunesigInfo$immune_sig, colnames(ITS_Ochiai3))]
  
  pdf(paste0("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSsimilarity/Ochiai/ITS_p_",cancer_type,"_bins.pdf"), width = 20, height = 20)
  
  pheatmap(ITS_Ochiai3,
           cluster_cols = F, 
           cluster_rows = F, 
           scale = "none",
           color = colorRampPalette(colors = c("blue","white","yellow","red"))(5),
           annotation_col = immunesigInfo2[c('Type', 'V1')], # 'SubType',
           annotation_row = immunesigInfo2[c('Type','V1')])  # 'SubType',
  dev.off()
  
  
  
  return(ITS_Ochiai)
})
names(ITS_p_set_Ochiai) <- names(ITS_p_set)

save(ITS_p_set_jcindex, ITS_p_set_Ochiai, file = "07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSsimilarity/ITS_gs_similarity.Rdata")

# ------------------------------------------------------------------------------
# 2) kappa statistics 
# ------------------------------------------------------------------------------
dir.create("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSsimilarity/kappa/")
ITS_p_set_kappa <- lapply(as.list(names(ITS_p_set)), function(cancer_type){
  # cancer_type = names(ITS_p_set)[1]
  ITS_set = ITS_p_set[[cancer_type]]
  
  ITS_mat <- matrix(0, length(ITS_set), length(unique(unlist(ITS_set))))
  rownames(ITS_mat) <- names(ITS_set)
  colnames(ITS_mat) <- unique(unlist(ITS_set))
  
  for(sig in names(ITS_set)){
    # sig=names(ITS_set)[1]
    x = ITS_set[[sig]]
    for(g in x){
      ITS_mat[which(rownames(ITS_mat) == sig), which(colnames(ITS_mat) == g)] = 1
    }
  }
  
  ITS_kappa <- list()
  
  for(i in 1:nrow(ITS_mat)){
    # i = 1
    ITS_kappa_res <- list()
    x = ITS_mat[i,]
    for(j in 1:nrow(ITS_mat)){
      
      y = ITS_mat[j,]
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
      
      df <- data.frame(i = rownames(ITS_mat)[i],
                       j = rownames(ITS_mat)[j], 
                       kappa = kappa, 
                       sd = sd, 
                       se = se,
                       CI_lower = CI[1],
                       CI_upper = CI[2],
                       z = z,
                       pvalue = pvalue)
      # K <- vcd::Kappa(matrix(c(a, b, c, d), byrow = T, nrow = 2))
      
      ITS_kappa_res[[j]] = df
    }
    ITS_kappa[[i]] = do.call(rbind, ITS_kappa_res)
  }
  
  ITS_kappa <- do.call(rbind, ITS_kappa)
  ITS_kappa2 <- dcast(ITS_kappa[c('i','j','kappa')],i~j)
  rownames(ITS_kappa2) <- ITS_kappa2$i
  ITS_kappa2 <- ITS_kappa2[-1]
  
  
  pdf(paste0("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSsimilarity/kappa/ITS_p_",cancer_type,".pdf"), width = 20, height = 20)
  tmp = pheatmap(ITS_kappa2,
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
  its_cls <- as.data.frame(as.matrix(v))
  # its_cls$immune_sig <- str_split(rownames(its_cls),'[.]',simplify = T)[,1]
  its_cls$immune_sig <- rownames(its_cls)
  its_cls = inner_join(treedf, its_cls)
  its_cls = its_cls[order(its_cls$V1),]
  its_cls$V1 = as.character(its_cls$V1)
  
  immunesigInfo2 <- merge(immunesigInfo, its_cls, by = 'immune_sig', all.y = T)
  rownames(immunesigInfo2) <- immunesigInfo2$immune_sig
  immunesigInfo2 <- immunesigInfo2[order(immunesigInfo2$V1),]

  ITS_kappa3 <- ITS_kappa2
  ITS_kappa3[ITS_kappa3>0.8] <- 1
  ITS_kappa3[ITS_kappa3>0.6 & ITS_kappa3<=0.8] <- 0.8
  ITS_kappa3[ITS_kappa3>0.4 & ITS_kappa3<=0.6] <- 0.6
  ITS_kappa3[ITS_kappa3<=0.4] <- 0.4
  ITS_kappa3 <- ITS_kappa3[its_cls$immune_sig, its_cls$immune_sig]
  
  # immunesigInfo <- immunesigInfo[order(immunesigInfo$Type), ]
  # ITS_kappa3 <- ITS_kappa3[intersect(immunesigInfo$immune_sig, rownames(ITS_kappa3)),
  #                          intersect(immunesigInfo$immune_sig, colnames(ITS_kappa3))]
  
  pdf(paste0("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSsimilarity/kappa/ITS_p_",cancer_type,"_bins.pdf"), width = 20, height = 20)
  
  pheatmap(ITS_kappa3,
           cluster_cols = F, 
           cluster_rows = F, 
           scale = "none",
           color = colorRampPalette(colors = c("blue","white","yellow","red"))(5),
           annotation_col = immunesigInfo2[c('Type','V1')], # 'SubType',
           annotation_row = immunesigInfo2[c('Type','V1')]) # 'SubType',
  dev.off()
  
  
  ITS_kappa$cancer_type = cancer_type
  return(ITS_kappa)
  
  
  # ITS_kappa_fmsb <- list()
  
  # for(i in 1:length(ITS_set)){
  #   # i = 1
  #   x = ITS_set[[i]]
  #   ITS_kappa_fmsb[[i]] = sapply(ITS_set, function(y){
  #     # y = ITS_set[[2]]
  #     K <- fmsb::Kappa.test(x, y)
  #     data.frame(estimate = K$Result$estimate,
  #                z = K$Result$statistic,
  #                p = K$Result$p.value)
  #     # length(intersect(x, y))
  #   })
  # }
  
  # names(ITS_kappa) = names(ITS_set)
  # ITS_kappa = do.call(cbind, ITS_kappa)
  
})

names(ITS_p_set_kappa) <- names(ITS_p_set)

save(ITS_p_set_kappa, file = "07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSsimilarity/ITS_gs_kappa.Rdata")




# ------------------------------------------------------------------------------
# 3) function similarity
# ------------------------------------------------------------------------------
dir.create("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSsimilarity/ITS_bioFunction/")
savepath = "07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSsimilarity/ITS_bioFunction/"


ITS_p_biofunction <- lapply(as.list(names(ITS_p_set)), function(cancer_type){
  # cancer_type = names(ITS_p_set)[1]
  ITS_set = ITS_p_set[[cancer_type]]
  ITS_function_detection(genelist = ITS_set,
                         sign = "p",
                         cancer_type = cancer_type,
                         savepath = savepath) 
})


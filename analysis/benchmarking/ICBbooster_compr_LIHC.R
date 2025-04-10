# rm(list=ls())
work_path <- "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/"
setwd(work_path)
options(stringsAsFactors = F)

library(getopt)

library(dplyr)
library(reshape2)
library(parallel)
library(patchwork)
library(openxlsx)
library(stringr)
library(ggplot2)
# save path setting ------------------------------------------------------------
savepath = "r1_drug_immuneSig_correlation_analysis/07_plotting_v2/ICBbooster_compr/"
dir.create(savepath)


# 70138 ------------------------------------------------------------------------
cancer='LIHC'
sampleType = 'Primary'
cancer_type = paste0(cancer,"_", sampleType)

savePath <- paste0("r1_drug_immuneSig_correlation_analysis/07_plotting_v2/ICBpredictor_booster_pred/")
saveAllPath_mergedresult = paste0(savePath, "/result_summarized/")
    
# TIDE base predictor
res_TIDEdrug1 <- read.csv( paste0(saveAllPath_mergedresult,  "/70138/", cancer_type, "_TIDE_merged.csv"))
res_TIDEdrug1 <- res_TIDEdrug1[res_TIDEdrug1$sig == 'TIDE', ]
res_TIDEdrug1 <- res_TIDEdrug1[c('drug_index', 'mNES')]
names(res_TIDEdrug1)=c('drug_index', 'TIDE')
res_TIDEdrug1 = res_TIDEdrug1[order(res_TIDEdrug1$TIDE, decreasing=T), ]
res_TIDEdrug1$TIDE_rank <- rank(-res_TIDEdrug1$TIDE, ties.method = "min")

# IPS base predictor
res_IPSdrug1 <- read.csv(paste0(saveAllPath_mergedresult,  "/70138/", cancer_type, "_IPS_merged.csv"))
res_IPSdrug1 <- res_IPSdrug1[res_IPSdrug1$sig == 'IPS', ]
res_IPSdrug1 <- res_IPSdrug1[c('drug_index', 'mNES')]
names(res_IPSdrug1)=c('drug_index', 'IPS')
res_IPSdrug1 = res_IPSdrug1[order(res_IPSdrug1$IPS, decreasing=T), ]
res_IPSdrug1$IPS_rank <- rank(-res_IPSdrug1$IPS, ties.method = "min")

# TIP base predictor
res_TIPdrug1 <- read.csv(paste0(saveAllPath_mergedresult,  "/70138/", cancer_type, "_TIP_merged.csv"))
res_TIPdrug1=res_TIPdrug1[res_TIPdrug1$sig == 'TIP_signature', ]
res_TIPdrug1 <- res_TIPdrug1[c('drug_index', 'mNES')]
names(res_TIPdrug1)=c('drug_index', 'TIP_signature')
res_TIPdrug1 = res_TIPdrug1[order(res_TIPdrug1$TIP_signature, decreasing=T), ]
res_TIPdrug1$TIP_rank <- rank(-res_TIPdrug1$TIP_signature, ties.method = "min")

# # OE base predictor
# res_OEdrug1 <- read.csv(paste0(saveAllPath_mergedresult,  "/70138/", cancer_type, "_OE_resu_merged.csv"))
# res_OEdrug1 <- res_OEdrug1[c('drug_index', 'sig', 'mNES')]
# res_OEdrug1=res_OEdrug1[res_OEdrug1$sig == 'resu', ]
# res_OEdrug1=res_OEdrug1[res_OEdrug1$sig == 'resu', ]
# res_OEdrug1 <- res_OEdrug1[c('drug_index', 'mNES')]
# names(res_OEdrug1)=c('drug_index', 'OE_resu')


# CM-Drug result
res_cmdrug_70138_good <- read.xlsx("r1_drug_immuneSig_correlation_analysis/07_plotting_v2/CM_drug_res/result/cp/70138_results_with_CM_Drug.xlsx")
res_cmdrug_70138_good$pert_iname_lowercase = tolower(res_cmdrug_70138_good$Compound)
res_cmdrug_70138_good = res_cmdrug_70138_good[order(res_cmdrug_70138_good$pert_score_max, decreasing=T), ]
res_cmdrug_70138_good$cmdrug.good_rank <- rank(-res_cmdrug_70138_good$pert_score_max, ties.method = "min")

res_cmdrug_70138_all <- read.xlsx("r1_drug_immuneSig_correlation_analysis/07_plotting_v2/CM_drug_res/result/cp/70138_all_results_with_CM_Drug.xlsx")
res_cmdrug_70138_all = res_cmdrug_70138_all[order(res_cmdrug_70138_all$pert_score_max, decreasing=T), ]
res_cmdrug_70138_all$cmdrug.all_rank <- rank(-res_cmdrug_70138_all$pert_score_max, ties.method = "min")
names(res_cmdrug_70138_all)[2:4] = paste0(names(res_cmdrug_70138_all)[2:4] , '_nofilter')
res_cmdrug_70138_all$pert_iname_lowercase = tolower(res_cmdrug_70138_all$Compound)

res_cmdrug = merge(res_cmdrug_70138_all, res_cmdrug_70138_good, by = c('Compound', 'pert_iname_lowercase'), all=T)

# IGeS-BS: simple sum
IGeS_BS_simpleSum <- read.csv("r1_drug_immuneSig_correlation_analysis/07_plotting_v2/IGeS_BS_simpleSum/TUMERIC_m2_meta_en_ACATP0.05/70138/allgenes/LIHC_Primary/drugITSRank/model_weight_unweighted_simple_res.txt", sep='\t', header=T)
names(IGeS_BS_simpleSum)[c(3:6)] <- paste0(names(IGeS_BS_simpleSum)[c(3:6)], '_simpleSum')
# IGeS-BS: weighted sum
IGeS_BS <- read.csv("r1_drug_immuneSig_correlation_analysis/06_results_summaryandplotting/results_meta_filter_ITSssgsea_CombatRef_17datasets/TUMERIC_m2_meta_en_ACATP0.05_LIHC/70138/allgenes/LIHC_Primary/drugITSRank/model_weight_weighted_glmnet_weight_res.txt", sep='\t', header=T)
IGeS_BS$pert_iname_lowercase = tolower(IGeS_BS$pert_iname)
names(IGeS_BS)[c(3:6)] <- paste0(names(IGeS_BS)[c(3:6)], '_IGeS_BS')

# merge results 
res_compr <- merge(merge(IGeS_BS, IGeS_BS_simpleSum, by = c('pert_iname', 'drug_index')),  res_cmdrug, by =  'pert_iname_lowercase', all.x=T )
res_compr <- merge(merge(merge(res_compr, res_TIDEdrug1, by =  'drug_index', all.x=T ), res_IPSdrug1, by =  'drug_index', all.x=T ), res_TIPdrug1, by =  'drug_index', all.x=T )
dim(res_compr)
res_compr$rank_IGeS_BS = rank(-res_compr$score_IGeS_BS, ties.method = "min")
res_compr$rank_simpleSum = rank(-res_compr$score_simpleSum, ties.method = "min")
res_compr$cmdrug.all_rank = rank(-res_compr$pert_score_max_nofilter, ties.method = "min")
res_compr$cmdrug.good_rank = rank(-res_compr$pert_score_max, ties.method = "min")
res_compr[is.na(res_compr$pert_score_max),]$cmdrug.good_rank=NA
res_compr$TIDE_rank = rank(-res_compr$TIDE, ties.method = "min")
res_compr$IPS_rank = rank(-res_compr$IPS, ties.method = "min")
res_compr$TIP_rank = rank(-res_compr$TIP_signature, ties.method = "min")
write.xlsx(res_compr, paste0(savepath,  "res_compr_70138.xlsx"))


res_compr
cor(res_compr$pert_score_max_nofilter, res_compr$score_IGeS_BS )
cor(res_compr$score_simpleSum, res_compr$score_IGeS_BS )
cor(res_compr$score_IGeS_BS, res_compr$TIDE )
cor(res_compr$score_IGeS_BS, res_compr$IPS )
cor(res_compr$score_IGeS_BS, res_compr$TIP_signature )
cor(res_compr$score_IGeS_BS, res_compr$pert_score_max, use="pairwise.complete.obs")

# annotation of drug evidence
top20 = res_compr[res_compr$rank_IGeS_BS <= 20, ]
write.xlsx(top20, paste0(savepath,  "res_compr_70138_top20.xlsx"))

drug_proof <- read.table("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/data/FDAapproved_p3Trial_IOtherapy/drug_list_proof_cancertype.txt", header = T, sep = '\t')
res_compr <- merge(res_compr, drug_proof, by=c('drug_index', 'pert_iname'), all.x=T)
write.xlsx(res_compr, paste0(savepath,  "res_compr_70138.xlsx"))
tmp=unique(res_compr[which(res_compr$if_with_ICB == 'yes'), ][-c(19:27)])
tmp=tmp[order(tmp$score_IGeS_BS,decreasing=T), ]
nrow(tmp[tmp$score_IGeS_BS > 0,])
nrow(tmp[tmp$score_simpleSum > 0,])
nrow(tmp[tmp$pert_score_detail > 0,])
nrow(tmp[tmp$TIDE > 0,])
nrow(tmp[tmp$IPS > 0,])
nrow(tmp[tmp$TIP_signature > 0,])


# 92742 ------------------------------------------------------------------------
savePath <- paste0("r1_drug_immuneSig_correlation_analysis/07_plotting_v2/ICBpredictor_booster_pred/")
saveAllPath_mergedresult = paste0(savePath, "/result_summarized/")
    
# TIDE base predictor
res_TIDEdrug1 <- read.csv( paste0(saveAllPath_mergedresult,  "/92742/", cancer_type, "_TIDE_merged.csv"))
res_TIDEdrug1 <- res_TIDEdrug1[res_TIDEdrug1$sig == 'TIDE', ]
res_TIDEdrug1 <- res_TIDEdrug1[c('drug_index', 'mNES')]
names(res_TIDEdrug1)=c('drug_index', 'TIDE')
res_TIDEdrug1 = res_TIDEdrug1[order(res_TIDEdrug1$TIDE, decreasing=T), ]
res_TIDEdrug1$TIDE_rank <- rank(-res_TIDEdrug1$TIDE, ties.method = "min")

# IPS base predictor
res_IPSdrug1 <- read.csv(paste0(saveAllPath_mergedresult,  "/92742/", cancer_type, "_IPS_merged.csv"))
res_IPSdrug1 <- res_IPSdrug1[res_IPSdrug1$sig == 'IPS', ]
res_IPSdrug1 <- res_IPSdrug1[c('drug_index', 'mNES')]
names(res_IPSdrug1)=c('drug_index', 'IPS')
res_IPSdrug1 = res_IPSdrug1[order(res_IPSdrug1$IPS, decreasing=T), ]
res_IPSdrug1$IPS_rank <- rank(-res_IPSdrug1$IPS, ties.method = "min")

# TIP base predictor
res_TIPdrug1 <- read.csv(paste0(saveAllPath_mergedresult,  "/92742/", cancer_type, "_TIP_merged.csv"))
res_TIPdrug1=res_TIPdrug1[res_TIPdrug1$sig == 'TIP_signature', ]
res_TIPdrug1 <- res_TIPdrug1[c('drug_index', 'mNES')]
names(res_TIPdrug1)=c('drug_index', 'TIP_signature')
res_TIPdrug1 = res_TIPdrug1[order(res_TIPdrug1$TIP_signature, decreasing=T), ]
res_TIPdrug1$TIP_rank <- rank(-res_TIPdrug1$TIP_signature, ties.method = "min")

# # OE base predictor
# res_OEdrug1 <- read.csv(paste0(saveAllPath_mergedresult,  "/92742/", cancer_type, "_OE_resu_merged.csv"))
# res_OEdrug1 <- res_OEdrug1[c('drug_index', 'sig', 'mNES')]
# res_OEdrug1=res_OEdrug1[res_OEdrug1$sig == 'resu', ]
# res_OEdrug1=res_OEdrug1[res_OEdrug1$sig == 'resu', ]
# res_OEdrug1 <- res_OEdrug1[c('drug_index', 'mNES')]
# names(res_OEdrug1)=c('drug_index', 'OE_resu')


# CM-Drug result
res_cmdrug_92742_good <- read.xlsx("r1_drug_immuneSig_correlation_analysis/07_plotting_v2/CM_drug_res/result_92742/cp/92742_results_with_CM_Drug.xlsx")
res_cmdrug_92742_good$pert_iname_lowercase = tolower(res_cmdrug_92742_good$Compound)
res_cmdrug_92742_good = res_cmdrug_92742_good[order(res_cmdrug_92742_good$pert_score_max, decreasing=T), ]
res_cmdrug_92742_good$cmdrug.good_rank <- rank(-res_cmdrug_92742_good$pert_score_max, ties.method = "min")

res_cmdrug_92742_all <- read.xlsx("r1_drug_immuneSig_correlation_analysis/07_plotting_v2/CM_drug_res/result_92742/cp/92742_all_results_with_CM_Drug.xlsx")
res_cmdrug_92742_all = res_cmdrug_92742_all[order(res_cmdrug_92742_all$pert_score_max, decreasing=T), ]
res_cmdrug_92742_all$cmdrug.all_rank <- rank(-res_cmdrug_92742_all$pert_score_max, ties.method = "min")
names(res_cmdrug_92742_all)[2:4] = paste0(names(res_cmdrug_92742_all)[2:4] , '_nofilter')
res_cmdrug_92742_all$pert_iname_lowercase = tolower(res_cmdrug_92742_all$Compound)

res_cmdrug = merge(res_cmdrug_92742_all, res_cmdrug_92742_good, by = c('Compound', 'pert_iname_lowercase'), all=T)


# IGeS-BS: simple sum
IGeS_BS_simpleSum <- read.csv("r1_drug_immuneSig_correlation_analysis/07_plotting_v2/IGeS_BS_simpleSum/TUMERIC_m2_meta_en_ACATP0.05/92742/allgenes/LIHC_Primary/drugITSRank/model_weight_unweighted_simple_res.txt", sep='\t', header=T)
names(IGeS_BS_simpleSum)[c(3:6)] <- paste0(names(IGeS_BS_simpleSum)[c(3:6)], '_simpleSum')
# IGeS-BS: weighted sum
IGeS_BS <- read.csv("r1_drug_immuneSig_correlation_analysis/06_results_summaryandplotting/results_meta_filter_ITSssgsea_CombatRef_17datasets/TUMERIC_m2_meta_en_ACATP0.05/92742/allgenes/LIHC_Primary/drugITSRank/model_weight_weighted_glmnet_weight_res.txt", sep='\t', header=T)
IGeS_BS$pert_iname_lowercase = tolower(IGeS_BS$pert_iname)
names(IGeS_BS)[c(3:6)] <- paste0(names(IGeS_BS)[c(3:6)], '_IGeS_BS')

# merge results 
res_compr <- merge(merge(IGeS_BS, IGeS_BS_simpleSum, by = c('pert_iname', 'drug_index')),  res_cmdrug, by =  'pert_iname_lowercase', all.x=T )
res_compr <- merge(merge(merge(res_compr, res_TIDEdrug1, by =  'drug_index', all.x=T ), res_IPSdrug1, by =  'drug_index', all.x=T ), res_TIPdrug1, by =  'drug_index', all.x=T )

top20 = res_compr[res_compr$rank_IGeS_BS <= 20, ]
write.xlsx(top20, paste0(savepath,  "res_compr_92742_top20.xlsx"))

cor(res_compr$score_simpleSum, res_compr$score_IGeS_BS )
cor(res_compr$score_IGeS_BS, res_compr$TIDE )
cor(res_compr$score_IGeS_BS, res_compr$IPS )
cor(res_compr$score_IGeS_BS, res_compr$TIP_signature )
cor(res_compr$score_IGeS_BS, res_compr$pert_score_max, use="pairwise.complete.obs")

top20 = res_compr[res_compr$rank_IGeS_BS <= 20, ]
write.xlsx(top20, paste0(savepath,  "res_compr_92742_top20.xlsx"))

res_compr2 = res_compr[-grep('BRD-', res_compr$pert_iname),]
nrow(res_compr2)
res_compr2$rank_IGeS_BS = rank(-res_compr2$score_IGeS_BS, ties.method = "min")
res_compr2$rank_simpleSum = rank(-res_compr2$score_simpleSum, ties.method = "min")
res_compr2$cmdrug.all_rank = rank(-res_compr2$pert_score_max_nofilter, ties.method = "min")
res_compr2$cmdrug.good_rank = rank(-res_compr2$pert_score_max, ties.method = "min")
res_compr2[is.na(res_compr2$pert_score_max),]$cmdrug.good_rank=NA
res_compr2$TIDE_rank = rank(-res_compr2$TIDE, ties.method = "min")
res_compr2$IPS_rank = rank(-res_compr2$IPS, ties.method = "min")
res_compr2$TIP_rank = rank(-res_compr2$TIP_signature, ties.method = "min")
write.xlsx(res_compr2, paste0(savepath,  "res_compr_92742_removeBRD.xlsx"))
top20 = res_compr2[res_compr2$rank_IGeS_BS <= 20, ]
write.xlsx(top20, paste0(savepath,  "res_compr_92742_removeBRD_top20.xlsx"))


# annotation of drug evidence
drug_proof <- read.table("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/data/FDAapproved_p3Trial_IOtherapy/drug_list_proof_cancertype.txt", header = T, sep = '\t')
res_compr <- merge(res_compr, drug_proof, by=c('drug_index', 'pert_iname'), all.x=T)
write.xlsx(res_compr, paste0(savepath,  "res_compr_92742.xlsx"))
tmp=unique(res_compr[which(res_compr$if_with_ICB == 'yes'), ][-c(19:27)])
nrow(tmp[tmp$score_IGeS_BS > 0,])
nrow(tmp[tmp$score_simpleSum > 0,])
nrow(tmp[tmp$pert_score_detail > 0,])
nrow(tmp[tmp$TIDE > 0,])
nrow(tmp[tmp$IPS > 0,])
nrow(tmp[tmp$TIP_signature > 0,])

res_compr1=read.xlsx(paste0(savepath,  "res_compr_70138.xlsx"))
res_compr2=read.xlsx(paste0(savepath,  "res_compr_92742_removeBRD.xlsx"))
res_compr=rbind(res_compr1,res_compr2)
cor(res_compr$TIP_signature, res_compr$score_IGeS_BS, method='spearman')
# plot -------------------------------------------------------------------------
res_compr=read.xlsx(paste0(savepath,  "res_compr_70138.xlsx"))

R3R1=c('score_IGeS_BS', 'TIDE','IPS','TIP_signature')
df=unique(res_compr[c('drug_index', R3R1)])
rownames(df) = df$drug_index
cor_matrix <- cor(df[-1], method = "spearman")
melted_cor <- melt(cor_matrix)
p = ggplot(data = melted_cor, aes(x=Var1, y=Var2, fill=value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(value,3)), color = "black", size = 4) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(-1,1), space = "Lab",
                       name="Spearman\nCorrelation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1)) +
  coord_fixed()
p


library(psych)
# corPlot(cor_matrix, numbers = TRUE, upper = FALSE, diag = FALSE, main = "Spearman Correlation Heatmap")
pdf(file = paste0(savepath,  "res_compr_R3R1.pdf"), width = 10, height =10)
res_compr=read.xlsx(paste0(savepath,  "res_compr_70138.xlsx"))
df=unique(res_compr[c('drug_index', R3R1)])
rownames(df) = df$drug_index

pairs.panels(
  df[-1],
  method = "spearman",      # Use Spearman correlation
  hist.col = "#a6cee3",     # Color of histograms on the diagonal
  density = TRUE,           # Include density plots on the diagonal
  ellipses = F,          # Add correlation ellipses
  pch = 19,                 # Point character (use 21:25 for background color options)
  bg = "#1f78b4",           # Background color of the data points
  col = "#33a02c",          # Color of the data points' borders
  cex=0.5,
  cex.axis = 1.5,  
  cex.cor = 1.2,            # Size of the correlation text
  lm = TRUE,
  star=TRUE)                # the significance levels of the correlation coefficients.
    df=unique(res_compr[c(R3R1, 'if_with_ICB')])
    df[which(df$if_with_ICB == 'no'), ]$if_with_ICB = 'unknown'
    df[is.na(df$if_with_ICB), ]$if_with_ICB = 'unknown'
    # Define colors for groups

    custom_colors <- c("yes" = "red", "unknown" = "grey")

    library(GGally)
    # Create the scatter plot matrix
    ggpairs(
    df, 
    columns = 1:5,  # Columns to compare against score_IGeS_BS
    aes(color = if_with_ICB),  # Color by group
    upper = list(continuous = wrap("cor", method = 'spearman', size = 4)), # Scatter plots in the upper triangle
    lower = list(continuous = "points"),  # Scatter plots in the lower triangle
    diag = list(continuous = "barDiag")  # Diagonal: histograms
    )+theme_bw()+  # Apply black-and-white theme
    scale_color_manual(values = custom_colors)  # Use custom colors for the groups


res_compr=read.xlsx(paste0(savepath,  "res_compr_92742.xlsx"))
df=unique(res_compr[c('drug_index', R3R1)])
rownames(df) = df$drug_index
pairs.panels(
  df[-1],
  method = "spearman",      # Use Spearman correlation
  hist.col = "#a6cee3",     # Color of histograms on the diagonal
  density = TRUE,           # Include density plots on the diagonal
  ellipses = F,          # Add correlation ellipses
  pch = 19,                 # Point character (use 21:25 for background color options)
  bg = "#1f78b4",           # Background color of the data points
  col = "#33a02c",          # Color of the data points' borders
  cex=0.2,
  cex.axis = 1.5,  
  cex.cor = 1.2,            # Size of the correlation text
  lm = TRUE,
  star=TRUE)  
  df=unique(res_compr[c(R3R1, 'if_with_ICB')])
    df[which(df$if_with_ICB == 'no'), ]$if_with_ICB = 'unknown'
    df[is.na(df$if_with_ICB), ]$if_with_ICB = 'unknown'
    # Define colors for groups

    custom_colors <- c("yes" = "red", "unknown" = "grey")

    library(GGally)
    # Create the scatter plot matrix
    ggpairs(
    df, 
    columns = 1:5,  # Columns to compare against score_IGeS_BS
    aes(color = if_with_ICB),  # Color by group
    upper = list(continuous = wrap("cor", method = 'spearman', size = 4)), # Scatter plots in the upper triangle
    lower = list(continuous = "points"),  # Scatter plots in the lower triangle
    diag = list(continuous = "barDiag")  # Diagonal: histograms
    )+theme_bw()+  # Apply black-and-white theme
    scale_color_manual(values = custom_colors)  # Use custom colors for the groups

    
dev.off()

R4=c('score_IGeS_BS', 'score_simpleSum')
df=unique(res_compr[c('drug_index', R4)])
rownames(df) = df$drug_index
cor_matrix <- cor(df[-1], method = "spearman")
melted_cor <- melt(cor_matrix)
p = ggplot(data = melted_cor, aes(x=Var1, y=Var2, fill=value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(value,3)), color = "black", size = 4) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(-1,1), space = "Lab",
                       name="Spearman\nCorrelation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1)) +
  coord_fixed()
p


library(psych)
# corPlot(cor_matrix, numbers = TRUE, upper = FALSE, diag = FALSE, main = "Spearman Correlation Heatmap")
pdf(file = paste0(savepath,  "res_compr_R4.pdf"), width = 10, height =10)
res_compr=read.xlsx(paste0(savepath,  "res_compr_70138.xlsx"))
df=unique(res_compr[c('drug_index', R4)])
rownames(df) = df$drug_index

pairs.panels(
  df[-1],
  method = "spearman",      # Use Spearman correlation
  hist.col = "#a6cee3",     # Color of histograms on the diagonal
  density = TRUE,           # Include density plots on the diagonal
  ellipses = F,          # Add correlation ellipses
  pch = 19,                 # Point character (use 21:25 for background color options)
  bg = "#1f78b4",           # Background color of the data points
  col = "#33a02c",          # Color of the data points' borders
  cex=0.5,
  cex.axis = 1.5,  
  cex.cor = 1.2,            # Size of the correlation text
  lm = TRUE,
  star=TRUE)                # the significance levels of the correlation coefficients.
    df=unique(res_compr[c(R4, 'if_with_ICB')])
    df[which(df$if_with_ICB == 'no'), ]$if_with_ICB = 'unknown'
    df[is.na(df$if_with_ICB), ]$if_with_ICB = 'unknown'
    # Define colors for groups

    custom_colors <- c("yes" = "red", "unknown" = "grey")

    library(GGally)
    # Create the scatter plot matrix
    ggpairs(
    df, 
    columns = 1:2,  # Columns to compare against score_IGeS_BS
    aes(color = if_with_ICB),  # Color by group
    upper = list(continuous = wrap("cor", method = 'spearman', size = 4)), # Scatter plots in the upper triangle
    lower = list(continuous = "points"),  # Scatter plots in the lower triangle
    diag = list(continuous = "barDiag")  # Diagonal: histograms
    )+theme_bw()+  # Apply black-and-white theme
    scale_color_manual(values = custom_colors)  # Use custom colors for the groups


res_compr=read.xlsx(paste0(savepath,  "res_compr_92742.xlsx"))
df=unique(res_compr[c('drug_index', R4)])
rownames(df) = df$drug_index
pairs.panels(
  df[-1],
  method = "spearman",      # Use Spearman correlation
  hist.col = "#a6cee3",     # Color of histograms on the diagonal
  density = TRUE,           # Include density plots on the diagonal
  ellipses = F,          # Add correlation ellipses
  pch = 19,                 # Point character (use 21:25 for background color options)
  bg = "#1f78b4",           # Background color of the data points
  col = "#33a02c",          # Color of the data points' borders
  cex=0.2,
  cex.axis = 1.5,  
  cex.cor = 1.2,            # Size of the correlation text
  lm = TRUE,
  star=TRUE)  
  df=unique(res_compr[c(R4, 'if_with_ICB')])
    df[which(df$if_with_ICB == 'no'), ]$if_with_ICB = 'unknown'
    df[is.na(df$if_with_ICB), ]$if_with_ICB = 'unknown'
    # Define colors for groups

    custom_colors <- c("yes" = "red", "unknown" = "grey")

    library(GGally)
    # Create the scatter plot matrix
    ggpairs(
    df, 
    columns = 1:2,  # Columns to compare against score_IGeS_BS
    aes(color = if_with_ICB),  # Color by group
    upper = list(continuous = wrap("cor", method = 'spearman', size = 4)), # Scatter plots in the upper triangle
    lower = list(continuous = "points"),  # Scatter plots in the lower triangle
    diag = list(continuous = "barDiag")  # Diagonal: histograms
    )+theme_bw()+  # Apply black-and-white theme
    scale_color_manual(values = custom_colors)  # Use custom colors for the groups

    
dev.off()

library(ggpubr)
library(tidyverse)

pdf(file = paste0(savepath,  "res_compr_70138.pdf"), width = 7, height =7)

for(model in c('score_simpleSum', 'TIDE','IPS','TIP_signature')){
  # model='score_simpleSum'
  # Add the regression line
  res_compr=read.xlsx(paste0(savepath,  "res_compr_70138.xlsx"))

  df=unique(res_compr[c('score_IGeS_BS', model)])
  names(df) = c('score_IGeS_BS', 'method')
  p=ggplot(df, aes(x=score_IGeS_BS, y=method)) + 
    geom_point(size = 1, shape = 19) + 
    stat_cor(p.accuracy = 0.001, r.accuracy = 0.01)+
  geom_smooth(method=lm)+
    labs(y = model, x = "IGeS_BS") + 
        theme_bw()+
        theme(plot.title = element_text(hjust = 0.5, size = 18), 
                            axis.title = element_text(size = 15), 
                            axis.title.y = element_text(vjust = 0.5), 
                            axis.text = element_text(size = 11))
      scale_color_manual(values = custom_colors) 
  print(p)

  df=unique(res_compr[c('score_IGeS_BS', model, 'if_with_ICB')])
  df[which(df$if_with_ICB == 'no'), ]$if_with_ICB = 'unknown'
  df[is.na(df$if_with_ICB), ]$if_with_ICB = 'unknown'
  names(df) = c('score_IGeS_BS', 'method','if_with_ICB')
  p=ggplot(df, aes(x=score_IGeS_BS, y=method, color=if_with_ICB)) + 
    geom_point(size = 1, shape = 19) + 
    stat_cor(p.accuracy = 0.001, r.accuracy = 0.01)+
    labs(y = model, x = "IGeS_BS") + 
        theme_bw()+
        theme(plot.title = element_text(hjust = 0.5, size = 18), 
                            axis.title = element_text(size = 15), 
                            axis.title.y = element_text(vjust = 0.5), 
                            axis.text = element_text(size = 11))
      scale_color_manual(values = custom_colors) 
  print(p)
}

dev.off()



pdf(file = paste0(savepath,  "res_compr_92742.pdf"), width = 7, height =7)

for(model in c('score_simpleSum', 'TIDE','IPS','TIP_signature')){
  # model='score_simpleSum'
  # Add the regression line
  res_compr=read.xlsx(paste0(savepath,  "res_compr_92742.xlsx"))


  df=unique(res_compr[c('score_IGeS_BS', model)])
  names(df) = c('score_IGeS_BS', 'method')
  p=ggplot(df, aes(x=score_IGeS_BS, y=method)) + 
    geom_point(size = 1, shape = 19) + 
    stat_cor(p.accuracy = 0.001, r.accuracy = 0.01)+
  geom_smooth(method=lm)+
    labs(y = model, x = "IGeS_BS") + 
        theme_bw()+
        theme(plot.title = element_text(hjust = 0.5, size = 18), 
                            axis.title = element_text(size = 15), 
                            axis.title.y = element_text(vjust = 0.5), 
                            axis.text = element_text(size = 11))
      scale_color_manual(values = custom_colors) 
  print(p)

  df=unique(res_compr[c('score_IGeS_BS', model, 'if_with_ICB')])
  df[which(df$if_with_ICB == 'no'), ]$if_with_ICB = 'unknown'
  df[is.na(df$if_with_ICB), ]$if_with_ICB = 'unknown'
  names(df) = c('score_IGeS_BS', 'method','if_with_ICB')
  p=ggplot(df, aes(x=score_IGeS_BS, y=method, color=if_with_ICB)) + 
    geom_point(size = 1, shape = 19) + 
    stat_cor(p.accuracy = 0.001, r.accuracy = 0.01)+
    labs(y = model, x = "IGeS_BS") + 
        theme_bw()+
        theme(plot.title = element_text(hjust = 0.5, size = 18), 
                            axis.title = element_text(size = 15), 
                            axis.title.y = element_text(vjust = 0.5), 
                            axis.text = element_text(size = 11))
      scale_color_manual(values = custom_colors) 
  print(p)
}
dev.off()

library(tidyr)             # Load the tidyr package
res_compr_70138=read.xlsx(paste0(savepath,  "res_compr_70138.xlsx"))
res_compr_70138$dataset='70138'
res_compr_92742=read.xlsx(paste0(savepath,  "res_compr_92742.xlsx"))
res_compr_92742$dataset='92742'
res_compr = rbind(res_compr_70138, res_compr_92742)
df=unique(res_compr[c('score_IGeS_BS','rank_IGeS_BS', 'score_simpleSum', 'TIDE','IPS','TIP_signature', 'if_with_ICB')])
df[which(df$if_with_ICB == 'no'), ]$if_with_ICB = 'unknown'
df[is.na(df$if_with_ICB), ]$if_with_ICB = 'unknown'
# Define colors for groups
nrow(df[df$score_IGeS_BS > 0 & df$if_with_ICB  == 'yes', ])
nrow(df[df$score_simpleSum > 0 & df$if_with_ICB  == 'yes', ])
nrow(df[df$TIDE > 0 & df$if_with_ICB  == 'yes', ])
nrow(df[df$IPS > 0 & df$if_with_ICB  == 'yes', ])
nrow(df[df$TIP_signature > 0 & df$if_with_ICB  == 'yes', ])

nrow(df[df$score_IGeS_BS < 0 & df$if_with_ICB  == 'yes', ])
nrow(df[df$score_simpleSum  < 0 & df$if_with_ICB  == 'yes', ])
nrow(df[df$TIDE  < 0 & df$if_with_ICB  == 'yes', ])
nrow(df[df$IPS  < 0 & df$if_with_ICB  == 'yes', ])
nrow(df[df$TIP_signature <  0 & df$if_with_ICB  == 'yes', ])


# Calculate counts for each condition
counts <- data.frame(
  variable = rep(c("IGeS_BS", "IGeS_simpleSum", "TIDE", 'IPS', "TIP_sig"), each = 2),
  score = rep(c("> 0", "<= 0"), times = 5),
  count = c(
    nrow(df[df$score_IGeS_BS > 0 & df$if_with_ICB == 'yes', ]),
    nrow(df[df$score_IGeS_BS <= 0 & df$if_with_ICB == 'yes', ]),
    nrow(df[df$score_simpleSum > 0 & df$if_with_ICB == 'yes', ]),
    nrow(df[df$score_simpleSum <= 0 & df$if_with_ICB == 'yes', ]),
    nrow(df[df$TIDE > 0 & df$if_with_ICB == 'yes', ]),
    nrow(df[df$TIDE <= 0 & df$if_with_ICB == 'yes', ]),
    nrow(df[df$IPS > 0 & df$if_with_ICB == 'yes', ]),
    nrow(df[df$IPS <= 0 & df$if_with_ICB == 'yes', ]),
    nrow(df[df$TIP_signature > 0 & df$if_with_ICB == 'yes', ]),
    nrow(df[df$TIP_signature <= 0 & df$if_with_ICB == 'yes', ])
  )
)

counts$variable=factor(counts$variable, levels=c("IGeS_BS", "IGeS_simpleSum", "TIDE", 'IPS', "TIP_sig"))
p=ggplot(counts, aes(x = variable, y = count, fill = score)) +
  geom_bar(stat = "identity", position = "stack") +  # Change position to stack
  labs(x = "Variable", y = "Count", title = "positive predictor Value (with ICB combo evidence)") +
  scale_fill_manual(values = c("> 0" = "lightblue", "<= 0" = "grey")) +
  theme_bw() +  # Use black and white theme
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability

ggsave(paste0(savepath,  "res_compr_positive_score.pdf"),p, width=5, height=5)



# Calculate counts for each condition
counts <- data.frame(
  variable = rep(c("IGeS_BS", "TIDE", 'IPS', "TIP_sig"), each = 2),
  score = rep(c("> 0", "<= 0"), times = 4),
  count = c(
    nrow(df[df$score_IGeS_BS > 0 & df$if_with_ICB == 'yes', ]),
    nrow(df[df$score_IGeS_BS <= 0 & df$if_with_ICB == 'yes', ]),
    nrow(df[df$TIDE > 0 & df$if_with_ICB == 'yes', ]),
    nrow(df[df$TIDE <= 0 & df$if_with_ICB == 'yes', ]),
    nrow(df[df$IPS > 0 & df$if_with_ICB == 'yes', ]),
    nrow(df[df$IPS <= 0 & df$if_with_ICB == 'yes', ]),
    nrow(df[df$TIP_signature > 0 & df$if_with_ICB == 'yes', ]),
    nrow(df[df$TIP_signature <= 0 & df$if_with_ICB == 'yes', ])
  )
)

counts$variable=factor(counts$variable, levels=c("IGeS_BS", "TIDE", 'IPS', "TIP_sig"))
p=ggplot(counts, aes(x = variable, y = count, fill = score)) +
  geom_bar(stat = "identity", position = "stack") +  # Change position to stack
  labs(x = "Variable", y = "Count", title = "positive predictor Value (with ICB combo evidence)") +
  scale_fill_manual(values = c("> 0" = "lightblue", "<= 0" = "grey")) +
  theme_bw() +  # Use black and white theme
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability

ggsave(paste0(savepath,  "res_compr_positive_score_R3R1.pdf"),p, width=5, height=5)


# rm(list=ls())
work_path <- "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/"
setwd(work_path)
options(stringsAsFactors = F)

library(getopt)

library(dplyr)
library(reshape2)
library(parallel)
library(patchwork)
library(openxlsx)
library(stringr)

# save path setting ------------------------------------------------------------
savepath = "r1_drug_immuneSig_correlation_analysis/07_plotting_v2/ICBbooster_compr/"
dir.create(savepath)


# annotation of drug evidence
drug_proof <- read.table("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/data/FDAapproved_p3Trial_IOtherapy/drug_list_proof_cancertype.txt", header = T, sep = '\t')


# CM-Drug result------------------------------------------------------------------------
res_cmdrug <- read.xlsx("r1_drug_immuneSig_correlation_analysis/07_plotting_v2/CM_drug_res/result/cp/results_with_CM_Drug.xlsx")

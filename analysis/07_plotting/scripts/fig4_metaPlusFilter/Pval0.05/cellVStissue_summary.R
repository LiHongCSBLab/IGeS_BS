# rm(list=ls())
setwd("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/")
options(stringsAsFactors = F)
library(ggplot2)
library(reshape2)
library(patchwork)
datasets <- list.files()
datasets <- datasets[grep("result_", datasets)]

dir.create("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/07_plotting_v2/04_cellVSmouse_metaPlusFilter/")
dir.create("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/07_plotting_v2/04_cellVSmouse_metaPlusFilter/padj0.05/")
plot_savepath = "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/07_plotting_v2/04_cellVSmouse_metaPlusFilter/padj0.05/cellVStissue_consistency/"
dir.create(plot_savepath)


immunesig <- read.csv("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/05_immuneSig_ICBresponse/results/meta_results/res_metaRes.csv")
immunesig <- immunesig[immunesig$Meta_Pval < 0.05, ]
names(immunesig) <- c("immune_sig", names(immunesig)[-1])
immunesig$flag <- "sensitive"
immunesig[immunesig$OR < 1,]$flag <- "resistant"
immunesig$setindex = "01_sensitive"
immunesig[immunesig$flag == "resistant", ]$setindex = "02_resistant"

load("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/05_immuneSig_ICBresponse/results_r0.4_ITSssgsea_CombatRef/finalmodel_17datasets_ssgseanorm/metaAnalysis_filter/en/elasticnet_model_pred.Rdata")

en_model_para = m_en$learner.model$glmnet.fit
m_en_imp <- data.frame(en_model_para$beta[,which(en_model_para$lambda == m_en$learner.model$lambda.min)])

names(m_en_imp) <- "weight"
m_en_imp$immune_sig = row.names(m_en_imp)
immunesig = merge(immunesig, m_en_imp, by = "immune_sig")
immunesig =immunesig[which(immunesig$weight != 0),]


# ------------------------------------------------------------------------------

valDS <- list.files()
valDS <- valDS[grep('result_GSE',valDS)]
valDS_immu <- c("GSE149825_meSKCM","GSE149825_priSKCM", 
                "GSE152925",
                "GSE120500", 
                "GSE160785_COAD", "GSE160785_READ",
                "GSE114601_jq+", "GSE114601_jq-")

res_immu <- list()

# GSE149825 --------------------------------------------------------------------
file2 = "result_GSE149825/drug_immunesig_TUMERIC_GSEA_result_r0.4_meta_oe_original_en_padj0.05/drug1487/CellvstreatedMouse.csv"
if(!file.exists(file2)){
  res_immu[[1]] <- data.frame(Var1 = NA, Freq = NA)
}else{
  df = read.csv(file2)
  names(df)[1] <- 'immune_sig'
  df = merge(df, immunesig, by = 'immune_sig')
  df$treated_mouse_binary = sign(df$treated_mouse)
  df$treated_cell_binary = sign(df$drugindex)
  df$consistency <- df$treated_mouse_binary * df$treated_cell_binary
  res_immu[[1]] <- as.data.frame(table(df$consistency))
}


file2 = "result_GSE149825/drug_immunesig_CPE_GSEA_result_r0.4_meta_oe_original_en_padj0.05/drug1487/CellvstreatedMouse.csv"
if(!file.exists(file2)){
  res_immu[[2]] <- data.frame(Var1 = NA, Freq = NA)
}else{
  df = read.csv(file2)
  names(df)[1] <- 'immune_sig'
  df = merge(df, immunesig, by = 'immune_sig')
  df$treated_mouse_binary = sign(df$treated_mouse)
  df$treated_cell_binary = sign(df$drugindex)
  df$consistency <- df$treated_mouse_binary * df$treated_cell_binary
  res_immu[[2]] <- as.data.frame(table(df$consistency))
}


# GSE152925 --------------------------------------------------------------------
file2 = "result_GSE152925/drug_immunesig_TUMERIC_GSEA_result_r0.4_meta_oe_original_en_padj0.05/drug888/CellvstreatedMouse.csv"
if(!file.exists(file2)){
  res_immu[[3]] <- data.frame(Var1 = NA, Freq = NA)
}else{
  df = read.csv(file2)
  names(df)[1] <- 'immune_sig'
  df = merge(df, immunesig, by = 'immune_sig')
  df$treated_mouse_binary = sign(df$treated_mouse)
  df$treated_cell_binary = sign(df$drugindex)
  df$consistency <- df$treated_mouse_binary * df$treated_cell_binary
  res_immu[[3]] <- as.data.frame(table(df$consistency))
}

# GSE120500 --------------------------------------------------------------------
file2 = "result_GSE120500/drug_immunesig_TUMERIC_GSEA_result_r0.4_meta_oe_original_en_padj0.05/drug53/CellvstreatedMouse.csv"
if(!file.exists(file2)){
  res_immu[[4]] <- data.frame(Var1 = NA, Freq = NA)
}else{
  df = read.csv(file2)
  names(df)[1] <- 'immune_sig'
  df = merge(df, immunesig, by = 'immune_sig')
  df$treated_mouse_binary = sign(df$treated_mouse)
  df$treated_cell_binary = sign(df$drugindex)
  df$consistency <- df$treated_mouse_binary * df$treated_cell_binary
  res_immu[[4]] <- as.data.frame(table(df$consistency))
}


# GSE160785 --------------------------------------------------------------------
file2 = "result_GSE160785/drug_immunesig_TUMERIC_GSEA_result_r0.4_meta_oe_original_en_padj0.05/COAD_Primary/drug989/CellvstreatedMouse.csv"
if(!file.exists(file2)){
  res_immu[[5]] <- data.frame(Var1 = NA, Freq = NA)
}else{
  df = read.csv(file2)
  names(df)[1] <- 'immune_sig'
  df = merge(df, immunesig, by = 'immune_sig')
  df$treated_mouse_binary = sign(df$treated_mouse)
  df$treated_cell_binary = sign(df$drugindex)
  df$consistency <- df$treated_mouse_binary * df$treated_cell_binary
  res_immu[[5]] <- as.data.frame(table(df$consistency))
}


file2 = "result_GSE160785/drug_immunesig_TUMERIC_GSEA_result_r0.4_meta_oe_original_en_padj0.05/READ_Primary/drug989/CellvstreatedMouse.csv"
if(!file.exists(file2)){
  res_immu[[6]] <- data.frame(Var1 = NA, Freq = NA)
}else{
  df = read.csv(file2)
  names(df)[1] <- 'immune_sig'
  df = merge(df, immunesig, by = 'immune_sig')
  df$treated_mouse_binary = sign(df$treated_mouse)
  df$treated_cell_binary = sign(df$drugindex)
  df$consistency <- df$treated_mouse_binary * df$treated_cell_binary
  res_immu[[6]] <- as.data.frame(table(df$consistency))
}


# GSE114601 --------------------------------------------------------------------
file2 = "result_GSE114601/drug_immunesig_TUMERIC_GSEA_result_r0.4_meta_oe_original_en_padj0.05/drug19/CellvstreatedMouse.csv"
if(!file.exists(file2)){
  res_immu[[5]] <- data.frame(Var1 = NA, Freq = NA)
}else{
  df = read.csv(file2)
  names(df)[1] <- 'immune_sig'
  df = merge(df, immunesig, by = 'immune_sig')
  df$treated_mouse_binary = sign(df$treated_mouse)
  df$treated_cell_binary = sign(df$drugindex)
  df$consistency <- df$treated_mouse_binary * df$treated_cell_binary
  res_immu[[7]] <- as.data.frame(table(df$consistency))
}


file2 = "result_GSE114601/drug_immunesig_TUMERIC_GSEA_result_r0.4_meta_oe_original_en_padj0.05/drug222/CellvstreatedMouse.csv"
if(!file.exists(file2)){
  res_immu[[6]] <- data.frame(Var1 = NA, Freq = NA)
}else{
  df = read.csv(file2)
  names(df)[1] <- 'immune_sig'
  df = merge(df, immunesig, by = 'immune_sig')
  df$treated_mouse_binary = sign(df$treated_mouse)
  df$treated_cell_binary = sign(df$drugindex)
  df$consistency <- df$treated_mouse_binary * df$treated_cell_binary
  res_immu[[8]] <- as.data.frame(table(df$consistency))
}


# ------------------------------------------------------------------------------
names(res_immu) <- valDS_immu
res_immu = res_immu[-3]
res_immu <- lapply(as.list(names(res_immu)), function(x){
    tmp = res_immu[[x]]
    names(tmp) <- c("Type","Consistency")
    tmp$all <- sum(tmp$Consistency)
    tmp$Consistency_ratio <- tmp$Consistency/sum(tmp$Consistency)
    tmp$dataset <- x
    return(tmp)
})
res_immu <- do.call(rbind, res_immu)
res_immu[res_immu$Type==1,]
write.csv(res_immu, paste0(plot_savepath, "cellVStissue_immu_summary.csv"), row.names = F, quote = F)


write.csv(res_immu, paste0(plot_savepath, "cellVStissue_summary.csv"), row.names = F, quote = F)


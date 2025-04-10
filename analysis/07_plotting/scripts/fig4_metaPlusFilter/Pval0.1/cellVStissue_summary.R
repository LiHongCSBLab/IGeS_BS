# rm(list=ls())
setwd("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/")
options(stringsAsFactors = F)
library(ggplot2)
library(reshape2)
library(patchwork)
datasets <- list.files()
datasets <- datasets[grep("result_", datasets)]

dir.create("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/07_plotting_v2/04_cellVSmouse_metaPlusFilter/")
dir.create("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/07_plotting_v2/04_cellVSmouse_metaPlusFilter/padj0.1/")
plot_savepath = "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/07_plotting_v2/04_cellVSmouse_metaPlusFilter/padj0.1/cellVStissue_consistency/"
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
file2 = "result_GSE149825/drug_immunesig_TUMERIC_GSEA_result_r0.4_meta_oe_original_en_padj0.1/drug1487/CellvstreatedMouse.csv"
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


file2 = "result_GSE149825/drug_immunesig_CPE_GSEA_result_r0.4_meta_oe_original_en_padj0.1/drug1487/CellvstreatedMouse.csv"
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
file2 = "result_GSE152925/drug_immunesig_TUMERIC_GSEA_result_r0.4_meta_oe_original_en_padj0.1/drug888/CellvstreatedMouse.csv"
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
file2 = "result_GSE120500/drug_immunesig_TUMERIC_GSEA_result_r0.4_meta_oe_original_en_padj0.1/drug53/CellvstreatedMouse.csv"
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
file2 = "result_GSE160785/drug_immunesig_TUMERIC_GSEA_result_r0.4_meta_oe_original_en_padj0.1/COAD_Primary/drug989/CellvstreatedMouse.csv"
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


file2 = "result_GSE160785/drug_immunesig_TUMERIC_GSEA_result_r0.4_meta_oe_original_en_padj0.1/READ_Primary/drug989/CellvstreatedMouse.csv"
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
file2 = "result_GSE114601/drug_immunesig_TUMERIC_GSEA_result_r0.4_meta_oe_original_en_padj0.1/drug19/CellvstreatedMouse.csv"
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


file2 = "result_GSE114601/drug_immunesig_TUMERIC_GSEA_result_r0.4_meta_oe_original_en_padj0.1/drug222/CellvstreatedMouse.csv"
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

p_valDS_immu = ggplot( res_immu, aes( x = factor(dataset), weight = Consistency_ratio, fill = factor(Type, levels = c(-1,0,1))))+
  geom_bar( position = "stack") +
  geom_text(aes(label = paste0(round((Consistency_ratio*100),2), "%"), y=Consistency_ratio), 
            position = position_stack(vjust = 0.5), size = 5) +
  coord_flip() +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14)) 

# p_valDS_immu

ggsave(paste0(plot_savepath, "cellVStissue_immu_summary.pdf"), p_valDS_immu , width = 5, height = 6)

# ------------------------------------------------------------------------------
valDS_PDX <- c("GSE60939_COAD_all", "GSE60939_COAD_PDX020", "GSE60939_COAD_PDX098", "GSE60939_COAD_PDX102", 
               "GSE60939_READ_all", "GSE60939_READ_PDX020", "GSE60939_READ_PDX098", "GSE60939_READ_PDX102", 
               "GSE148538", 
               "GSE155923", 
               "GSE33366_letrozole", "GSE33366_tamoxifen", "GSE33366_BMS754807")
res_PDX <- list()

# PDX: GSE60939 --------------------------------------------------------------------
# COAD Cabozantinib all pdx ---------------------------------------------------------
file2 = "result_GSE60939/drug_immunesig_TUMERIC_GSEA_result_r0.4_meta_oe_original_en_padj0.1/COAD_drug170/CellvstreatedMouse.csv"
if(!file.exists(file2)){
  res_PDX[[1]] <- data.frame(Var1 = NA, Freq = NA)
}else{
  df = read.csv(file2)
  names(df)[1] <- 'immune_sig'
  df = merge(df, immunesig, by = 'immune_sig')
  df$treated_mouse_binary = sign(df$treated_mouse)
  df$treated_cell_binary = sign(df$drugindex)
  df$consistency <- df$treated_mouse_binary * df$treated_cell_binary
  res_PDX[[1]] <- as.data.frame(table(df$consistency))
}

# COAD Cabozantinib PDX020 ---------------------------------------------------------
file2 = "result_GSE60939/drug_immunesig_TUMERIC_GSEA_result_r0.4_meta_oe_original_PDX020_en_padj0.1/COAD_drug170/CellvstreatedMouse.csv"
if(!file.exists(file2)){
  res_PDX[[2]] <- data.frame(Var1 = NA, Freq = NA)
}else{
  df = read.csv(file2)
  names(df)[1] <- 'immune_sig'
  df = merge(df, immunesig, by = 'immune_sig')
  df$treated_mouse_binary = sign(df$treated_mouse)
  df$treated_cell_binary = sign(df$drugindex)
  df$consistency <- df$treated_mouse_binary * df$treated_cell_binary
  res_PDX[[2]] <- as.data.frame(table(df$consistency))
}

# COAD Cabozantinib PDX098 ---------------------------------------------------------
file2 = "result_GSE60939/drug_immunesig_TUMERIC_GSEA_result_r0.4_meta_oe_original_PDX098_en_padj0.1/COAD_drug170/CellvstreatedMouse.csv"
if(!file.exists(file2)){
  res_PDX[[3]] <- data.frame(Var1 = NA, Freq = NA)
}else{
  df = read.csv(file2)
  names(df)[1] <- 'immune_sig'
  df = merge(df, immunesig, by = 'immune_sig')
  df$treated_mouse_binary = sign(df$treated_mouse)
  df$treated_cell_binary = sign(df$drugindex)
  df$consistency <- df$treated_mouse_binary * df$treated_cell_binary
  res_PDX[[3]] <- as.data.frame(table(df$consistency))
}

# COAD Cabozantinib PDX102 ---------------------------------------------------------
file2 = "result_GSE60939/drug_immunesig_TUMERIC_GSEA_result_r0.4_meta_oe_original_PDX102_en_padj0.1/COAD_drug170/CellvstreatedMouse.csv"
if(!file.exists(file2)){
  res_PDX[[4]] <- data.frame(Var1 = NA, Freq = NA)
}else{
  df = read.csv(file2)
  names(df)[1] <- 'immune_sig'
  df = merge(df, immunesig, by = 'immune_sig')
  df$treated_mouse_binary = sign(df$treated_mouse)
  df$treated_cell_binary = sign(df$drugindex)
  df$consistency <- df$treated_mouse_binary * df$treated_cell_binary
  res_PDX[[4]] <- as.data.frame(table(df$consistency))
}

# READ Cabozantinib all pdx ---------------------------------------------------------
file2 = "result_GSE60939/drug_immunesig_TUMERIC_GSEA_result_r0.4_meta_oe_original_en_padj0.1/READ_drug170/CellvstreatedMouse.csv"
if(!file.exists(file2)){
  res_PDX[[5]] <- data.frame(Var1 = NA, Freq = NA)
}else{
df = read.csv(file2)
names(df)[1] <- 'immune_sig'
df = merge(df, immunesig, by = 'immune_sig')
df$treated_mouse_binary = sign(df$treated_mouse)
df$treated_cell_binary = sign(df$drugindex)
df$consistency <- df$treated_mouse_binary * df$treated_cell_binary
res_PDX[[5]] <- as.data.frame(table(df$consistency))
}

# READ Cabozantinib PDX020 ---------------------------------------------------------
file2 = "result_GSE60939/drug_immunesig_TUMERIC_GSEA_result_r0.4_meta_oe_original_PDX020_en_padj0.1/READ_drug170/CellvstreatedMouse.csv"
if(!file.exists(file2)){
  res_PDX[[6]] <- data.frame(Var1 = NA, Freq = NA)
}else{
  df = read.csv(file2)
  names(df)[1] <- 'immune_sig'
  df = merge(df, immunesig, by = 'immune_sig')
  df$treated_mouse_binary = sign(df$treated_mouse)
  df$treated_cell_binary = sign(df$drugindex)
  df$consistency <- df$treated_mouse_binary * df$treated_cell_binary
  res_PDX[[6]] <- as.data.frame(table(df$consistency))
}


# READ Cabozantinib PDX098 ---------------------------------------------------------
file2 = "result_GSE60939/drug_immunesig_TUMERIC_GSEA_result_r0.4_meta_oe_original_PDX098_en_padj0.1/READ_drug170/CellvstreatedMouse.csv"
if(!file.exists(file2)){
  res_PDX[[7]] <- data.frame(Var1 = NA, Freq = NA)
}else{
  df = read.csv(file2)
  names(df)[1] <- 'immune_sig'
  df = merge(df, immunesig, by = 'immune_sig')
  df$treated_mouse_binary = sign(df$treated_mouse)
  df$treated_cell_binary = sign(df$drugindex)
  df$consistency <- df$treated_mouse_binary * df$treated_cell_binary
  res_PDX[[7]] <- as.data.frame(table(df$consistency))
}

# READ Cabozantinib PDX102 ---------------------------------------------------------
file2 = "result_GSE60939/drug_immunesig_TUMERIC_GSEA_result_r0.4_meta_oe_original_PDX102_en_padj0.1/READ_drug170/CellvstreatedMouse.csv"
if(!file.exists(file2)){
  res_PDX[[8]] <- data.frame(Var1 = NA, Freq = NA)
}else{
  df = read.csv(file2)
  names(df)[1] <- 'immune_sig'
  df = merge(df, immunesig, by = 'immune_sig')
  df$treated_mouse_binary = sign(df$treated_mouse)
  df$treated_cell_binary = sign(df$drugindex)
  df$consistency <- df$treated_mouse_binary * df$treated_cell_binary
  res_PDX[[8]] <- as.data.frame(table(df$consistency))
}

# PDX: GSE148538 --------------------------------------------------------------------
# v02 GSE148538 Cabozantinib all pdx ---------------------------------------------------------
file2 = "result_GSE148538/drug_immunesig_TUMERIC_GSEA_result_r0.4_meta_oe_original_en_padj0.1/drug170/CellvstreatedMouse.csv"
if(!file.exists(file2)){
  res_PDX[[9]] <- data.frame(Var1 = NA, Freq = NA)
}else{
  df = read.csv(file2)
  names(df)[1] <- 'immune_sig'
  df = merge(df, immunesig, by = 'immune_sig')
  df$treated_mouse_binary = sign(df$treated_mouse)
  df$treated_cell_binary = sign(df$drugindex)
  df$consistency <- df$treated_mouse_binary * df$treated_cell_binary
  res_PDX[[9]] <- as.data.frame(table(df$consistency))
}

# PDX: GSE155923 --------------------------------------------------------------------
# v02	GSE155923 Everolimus---------------------------------------------------------
file2 = "result_GSE155923/drug_immunesig_TUMERIC_GSEA_result_r0.4_meta_oe_original_en_padj0.1/drug196/CellvstreatedMouse.csv"
if(!file.exists(file2)){
  res_PDX[[10]] <- data.frame(Var1 = NA, Freq = NA)
}else{
  df = read.csv(file2)
  names(df)[1] <- 'immune_sig'
  df = merge(df, immunesig, by = 'immune_sig')
  df$treated_mouse_binary = sign(df$treated_mouse)
  df$treated_cell_binary = sign(df$drugindex)
  df$consistency <- df$treated_mouse_binary * df$treated_cell_binary
  res_PDX[[10]] <- as.data.frame(table(df$consistency))
}

# PDX: GSE33366 --------------------------------------------------------------------
# v02 GSE33366 letrozole---------------------------------------------------------
file2 = "result_GSE33366/letrozole_immunesig_TUMERIC_GSEA_result_r0.4_meta_oe_original_en_padj0.1/drug789/CellvstreatedMouse.csv"
if(!file.exists(file2)){
  res_PDX[[11]] <- data.frame(Var1 = NA, Freq = NA)
}else{
  df = read.csv(file2)
  names(df)[1] <- 'immune_sig'
  df = merge(df, immunesig, by = 'immune_sig')
  df$treated_mouse_binary = sign(df$treated_mouse)
  df$treated_cell_binary = sign(df$drugindex)
  df$consistency <- df$treated_mouse_binary * df$treated_cell_binary
  res_PDX[[11]] <- as.data.frame(table(df$consistency))
}

# v02 GSE33366 tamoxifen---------------------------------------------------------
file2 = "result_GSE33366/tamoxifen_immunesig_TUMERIC_GSEA_result_r0.4_meta_oe_original_en_padj0.1/drug1792/CellvstreatedMouse.csv"
if(!file.exists(file2)){
  res_PDX[[12]] <- data.frame(Var1 = NA, Freq = NA)
}else{
  df = read.csv(file2)
  names(df)[1] <- 'immune_sig'
  df = merge(df, immunesig, by = 'immune_sig')
  df$treated_mouse_binary = sign(df$treated_mouse)
  df$treated_cell_binary = sign(df$drugindex)
  df$consistency <- df$treated_mouse_binary * df$treated_cell_binary
  res_PDX[[12]] <- as.data.frame(table(df$consistency))
}

# v02 GSE33366 BMS-754807---------------------------------------------------------
file2 = "result_GSE33366/bms754807_immunesig_TUMERIC_GSEA_result_r0.4_meta_oe_original_en_padj0.1/drug3045/CellvstreatedMouse.csv"
if(!file.exists(file2)){
  res_PDX[[13]] <- data.frame(Var1 = NA, Freq = NA)
}else{
  df = read.csv(file2)
  names(df)[1] <- 'immune_sig'
  df = merge(df, immunesig, by = 'immune_sig')
  df$treated_mouse_binary = sign(df$treated_mouse)
  df$treated_cell_binary = sign(df$drugindex)
  df$consistency <- df$treated_mouse_binary * df$treated_cell_binary
  res_PDX[[13]] <- as.data.frame(table(df$consistency))
}

# ------------------------------------------------------------------------------
names(res_PDX) <- valDS_PDX

res_PDX <- lapply(as.list(names(res_PDX)), function(x){
    tmp = res_PDX[[x]]
    names(tmp) <- c("Type","Consistency")
    tmp$all <- sum(tmp$Consistency)
    tmp$Consistency_ratio <- tmp$Consistency/sum(tmp$Consistency)
    tmp$dataset <- x
    return(tmp)
})
res_PDX <- do.call(rbind, res_PDX)
res_PDX[res_PDX$Type ==1,]
write.csv(res_PDX, paste0(plot_savepath, "cellVStissue_PDX_summary.csv"), row.names = F, quote = F)


p_valDS_PDX = ggplot( res_PDX, aes( x = factor(dataset), weight = Consistency_ratio, fill = factor(Type, levels = c(-1,0,1))))+
  geom_bar( position = "stack")+
  geom_text(aes(label = paste0(round((Consistency_ratio*100),2), "%"), y=Consistency_ratio), 
            position = position_stack(vjust = 0.5), size = 5)+
  coord_flip() +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14)) 


# p_valDS_PDX
ggsave(paste0(plot_savepath, "cellVStissue_PDX_summary.pdf"), p_valDS_PDX , width = 5, height = 8)

p = p_valDS_immu / p_valDS_PDX + 
  plot_layout(height = c(1, 2))
p
ggsave(paste0(plot_savepath, "cellVStissue_summary.pdf"), p , width = 7, height = 12)



# ------------------------------------------------------------------------------
res_immu2 = res_immu[which(res_immu$all > 5), ]
p_valDS_immu2 = ggplot(res_immu, aes( x = factor(dataset), weight = Consistency_ratio, fill = factor(Type, levels = c(-1,0,1))))+
  geom_bar( position = "stack") +
  geom_text(aes(label = paste0(round((Consistency_ratio*100),2), "%"), y=Consistency_ratio), 
            position = position_stack(vjust = 0.5), size = 5) +
  coord_flip() +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14)) 

# p_valDS_immu2
res_PDX2 = res_PDX[which(res_PDX$all > 5), ]

p_valDS_PDX2 = ggplot(res_PDX2, aes( x = factor(dataset), weight = Consistency_ratio, fill = factor(Type, levels = c(-1,0,1))))+
  geom_bar( position = "stack")+
  geom_text(aes(label = paste0(round((Consistency_ratio*100),2), "%"), y=Consistency_ratio), 
            position = position_stack(vjust = 0.5), size = 5)+
  coord_flip() +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14)) 


# p_valDS_PDX2

p2 = p_valDS_immu2 / p_valDS_PDX2 + 
  plot_layout(height = c(1, 2))
p2

ggsave(paste0(plot_savepath, "cellVStissue_summary_allgreater5.pdf"), p2 , width = 7, height = 12)


# ------------------------------------------------------------------------------
res_immu$mousetype='withImmuneSystem'
res_PDX$mousetype='PDX'
res = rbind(res_immu, res_PDX)

write.csv(res, paste0(plot_savepath, "cellVStissue_summary.csv"), row.names = F, quote = F)


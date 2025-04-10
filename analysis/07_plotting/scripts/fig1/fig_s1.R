# rm(list=ls())
work_path <- "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/"
setwd(work_path)
options(stringsAsFactors = F)
library(ggplot2)
library(pheatmap)
library(aplot)
library(reshape2)
library(ggsci)
library(see)
library(patchwork)


source("03_drug_immuneSig_enrichment/scripts/02_GSEA_geneset_preparation_function.R")
source("03_drug_immuneSig_enrichment/scripts/02_GSEA_geneset_merge_metaAnalysis_function.R")

source("06_results_summaryandplotting/scripts/immunesig_selection_function.R")
source("06_results_summaryandplotting/scripts/ITS_analysis/ITS_function_detection.R")
source("06_results_summaryandplotting/scripts/ITS_analysis/ITSexpr_purity_detection.R")
source("06_results_summaryandplotting/scripts/ITS_analysis/ITSexpr_purity_detection.R")
source("06_results_summaryandplotting/scripts/ITS_analysis/immuneModulatorDetection.R")


dir.create("07_plotting_v2/01_metaAnalysis/")
dir.create("07_plotting_v2/01_metaAnalysis/fig_s1/")

immuneSigInfoPath = "06_results_summaryandplotting/data/immene_cell_hieracy/"

immunesigInfo <- as.data.frame(readxl::read_excel(paste0(immuneSigInfoPath, "immune_sig_hieracy_new.xlsx"), sheet = 1))
immunesigInfo <- immunesigInfo[c(1,4,5:9)]
immunesigInfo$immune_sig = tolower(immunesigInfo$immune_sig)
immunesigInfo = immunesig_name_unify(immunesigInfo)


library(tidyverse)

meta_result <- read.csv("05_immuneSig_ICBresponse/results/meta_results/res_metaRes.csv")
meta_selected <- meta_result[meta_result$Meta_Pval < 0.05, ]
meta_selected$immune_sig = meta_selected$X
# meta_selected <- meta_selected[-c(grep("caf_tide", meta_selected$immune_sig),
#                                   grep("m2_tide", meta_selected$immune_sig)), ]

# ------------------------------------------------------------------------------
# 1. overall plot -----------------------------------------------------------------
# ------------------------------------------------------------------------------
# 1) density plot: OR of all immune signatures, marked P-values

# plot(density(meta_result$OR))
# plot(density(meta_result$Meta_Pval))
dt_dot = meta_result
dt_dot$col = "black"
dt_dot[dt_dot$Meta_Pval < 0.05 & dt_dot$OR > 1, ]$col = "#e79999"
# dt_dot[dt_dot$Meta_Pval < 0.05 & dt_dot$OR > 1.5 , ]$col = '#B33333'
dt_dot[dt_dot$Meta_Pval < 0.05 & dt_dot$OR < 1 , ]$col = "#338080"
# dt_dot[dt_dot$Meta_Pval < 0.05 & dt_dot$OR < 0.5 , ]$col =  "#adfdfd73"


p_dot <- ggplot(data = dt_dot, mapping = aes(x = log10(Meta_Pval), y = OR)) + 
  geom_point(colour = dt_dot$col, size = 2)+
  geom_vline(aes(xintercept=log10(0.05)),linetype=3,col="red")+
  geom_hline(aes(yintercept=1),linetype=3,col="red")+
  # geom_hline(aes(yintercept=1.5),linetype=3,col="red")+
  geom_hline(aes(yintercept=0.5),linetype=3,col="red") +
  scale_x_reverse() + # x轴翻转
  theme_bw() # 不要背景

p_dot


# 2) pie chart: persentage of selected immune signatures, and their classification
library(scales)
blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )

df_pie = data.frame(group = c("sensitive", "resistant", "filtered"),
                    value = c(nrow(meta_selected[meta_selected$OR > 1, ]), 
                              nrow(meta_selected[meta_selected$OR < 1, ]), 
                              nrow(meta_result[meta_result$Meta_Pval > 0.05, ])),
                    col = c("#e79999", "#338080", "gray"))
df_pie$percent = df_pie$value / sum(df_pie$value) * 100
pie <- ggplot(df_pie, aes(x="", y=percent, fill=group))+
  geom_bar( width = 1, stat = "identity", alpha = 0.5) + 
  geom_col(color='black')+
  coord_polar("y", start=0) +
  scale_fill_manual(values=c("gray70", "#338080","#e79999"))+
  blank_theme +
  theme(axis.text.x=element_blank())+
  geom_label(aes(label = paste0(round(percent,2), "%")), size = 5, position = position_stack(vjust = 0.5))
pie

# 
immunesigAnnot <- merge(unique(meta_selected[c("immune_sig")]), immunesigInfo, by='immune_sig')
dfpie_immunesigIndex = as.data.frame(table(immunesigAnnot$Classification))
names(dfpie_immunesigIndex) = c("group","value")
dfpie_immunesigIndex$percent = dfpie_immunesigIndex$value / sum(dfpie_immunesigIndex$value) * 100

pie_immunesigIndex <- ggplot(dfpie_immunesigIndex, aes(x="", y=percent, fill=group))+
  geom_bar( width = 1, stat = "identity", alpha = 0.5) + 
  geom_col(color='black')+
  coord_polar("y", start=0) +
  # scale_fill_brewer(palette="Dark2") +
  blank_theme +
  theme(axis.text.x=element_blank())+
  geom_label(aes(label = paste0(round(percent,2), "%")), size = 5, position = position_stack(vjust = 0.5))
pie_immunesigIndex

immunesigAnnot <- merge(unique(meta_selected[c("immune_sig")]), immunesigInfo, by='immune_sig')
dfpie_immunesigCL = as.data.frame(table(immunesigAnnot$Type))
names(dfpie_immunesigCL) = c("group","value")
dfpie_immunesigCL = dfpie_immunesigCL[order(dfpie_immunesigCL$value, decreasing = T), ]
dfpie_immunesigCL$percent = dfpie_immunesigCL$value / sum(dfpie_immunesigCL$value) * 100
pie_immunesigCL <- ggplot(dfpie_immunesigCL, aes(x="", y=percent, fill=group))+
#   geom_bar(stat = "identity", alpha = 0.5) + 
  geom_col(color='black')+
  coord_polar("y", start=0) +
  # scale_fill_brewer(palette="Dark2") + 
  blank_theme +
  theme(axis.text.x=element_blank())+
  geom_label(aes(label = paste0(round(percent,2), "%")), size = 5, position = position_stack(vjust = 0.5))

pie_immunesigCL


ggsave("07_plotting_v2/01_metaAnalysis/fig_s1/meta_dot_plot.pdf", p_dot, width = 5, height = 5)
ggsave("07_plotting_v2/01_metaAnalysis/fig_s1/pie_meta_selected_plot.pdf", pie, width = 5, height = 4)
ggsave("07_plotting_v2/01_metaAnalysis/fig_s1/pie_immunesigIndex.pdf", pie_immunesigIndex, width = 5, height = 4)
ggsave("07_plotting_v2/01_metaAnalysis/fig_s1/pie_immunesigCL.pdf", pie_immunesigCL, width = 10, height = 10)



p = (p_dot / pie )/ pie_immunesigIndex | pie_immunesigCL
p
ggsave("07_plotting_v2/01_metaAnalysis/fig_s1.pdf", p, width = 20, height = 20)

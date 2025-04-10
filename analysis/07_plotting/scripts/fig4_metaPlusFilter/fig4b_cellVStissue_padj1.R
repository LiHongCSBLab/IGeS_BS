# rm(list=ls())
setwd("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/")
options(stringsAsFactors = F)
library(ggplot2)
library(reshape2)
library(patchwork)

plot_savepath= "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/07_plotting_v2/04_cellVSmouse_metaPlusFilter/fig4b_padj1_"

res_savepath  = "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/07_plotting_v2/04_cellVSmouse_metaPlusFilter/padj1/cellVStissue_consistency/"
res <- read.csv(paste0(res_savepath, "cellVStissue_summary.csv"))

res_immu = res[res$mousetype == 'withImmuneSystem', ]
res_immu = res_immu[is.element(res_immu$dataset, c('GSE149825_meSKCM', 'GSE152925', 'GSE120500', 'GSE160785_COAD', 'GSE114601_jq+')), ]
res_immu$Type = factor(res_immu$Type, levels = c(-1,0,1))
res_immu$dataset = factor(res_immu$dataset)
p_valDS_immu = ggplot( res_immu, aes(x = dataset, weight = Consistency_ratio, fill = factor(Type, levels = c(-1,0,1))))+
  geom_bar( position = "stack") +
  geom_text(aes(label = paste0(round((Consistency_ratio*100),2), "%"), y=Consistency_ratio), 
            position = position_stack(vjust = 0.5), size = 5) +
  coord_flip() +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14)) 


res_PDX = res[res$mousetype == 'PDX', ]
res_PDX = res_PDX[is.element(res_PDX$dataset, 
                            c('GSE60939_COAD_all', 'GSE60939_COAD_PDX020', 'GSE60939_COAD_PDX098', 'GSE60939_COAD_PDX102', 
                              'GSE148538', 'GSE155923', 'GSE33366_letrozole', 'GSE33366_tamoxifen')), ]
res_PDX$Type = factor(res_PDX$Type, levels = c(-1,0,1))
res_PDX$dataset = factor(res_PDX$dataset)
p_valDS_PDX = ggplot(res_PDX, aes( x = dataset, weight = Consistency_ratio, fill = Type))+
  geom_bar( position = "stack")+
  geom_text(aes(label = paste0(round((Consistency_ratio*100),2), "%"), y=Consistency_ratio), 
            position = position_stack(vjust = 0.5), size = 5)+
  coord_flip() +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14)) 

p = p_valDS_immu / p_valDS_PDX + 
  plot_layout(height = c(2, 3))
p
ggsave(paste0(plot_savepath, "cellVStissue_summary.pdf"), p , width = 10, height = 12)



p_valDS_immu2 = ggplot(res_immu, aes( x = dataset, 
                           y = Consistency,
                           fill = factor(Type, levels = c(-1,0,1))))+
  geom_bar( position = "stack", stat = "identity", width = 0.6 ) +
  geom_text(aes(label = Consistency, y=Consistency), 
            position = position_stack(vjust = 0.5), size = 5)+
  coord_flip() +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14)) 


p_valDS_PDX2 = ggplot(res_PDX, aes( x = dataset, 
                           y = Consistency,
                           fill = factor(Type, levels = c(-1,0,1))))+
  geom_bar( position = "stack", stat = "identity", width = 0.6 ) +
  geom_text(aes(label = Consistency, y=Consistency), 
            position = position_stack(vjust = 0.5), size = 5)+
  coord_flip() +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14)) 

p2 = p_valDS_immu2 / p_valDS_PDX2 + 
  plot_layout(height = c(1, 2))
p2
ggsave(paste0(plot_savepath, "cellVStissue_summary2.pdf"), p2 , width = 10, height = 12)


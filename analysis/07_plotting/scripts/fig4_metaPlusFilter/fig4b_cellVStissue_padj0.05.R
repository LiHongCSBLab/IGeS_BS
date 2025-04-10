# rm(list=ls())
setwd("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/")
options(stringsAsFactors = F)
library(ggplot2)
library(reshape2)
library(patchwork)

plot_savepath= "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/07_plotting_v2/04_cellVSmouse_metaPlusFilter/fig4b_padj0.05_"

res_savepath  = "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/07_plotting_v2/04_cellVSmouse_metaPlusFilter/padj0.05/cellVStissue_consistency/"
res <- read.csv(paste0(res_savepath, "cellVStissue_summary.csv"))

res_immu = res
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

p_valDS_immu 
ggsave(paste0(plot_savepath, "cellVStissue_summary.pdf"), p_valDS_immu , width = 10, height = 6)



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

p_valDS_immu2
ggsave(paste0(plot_savepath, "cellVStissue_summary2.pdf"), p_valDS_immu2 , width = 10, height = 6)


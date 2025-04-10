# rm(list=ls())
setwd("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/")
options(stringsAsFactors = F)
library(ggplot2)
library(reshape2)
library(patchwork)
datasets <- list.files()
datasets <- datasets[grep("result_", datasets)]

dir.create("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/07_plotting_v2/04_cellVSmouse_metaPlusFilter/padj1/")

# COAD GSE160785 celecoxib---------------------------------------------------------
file2 = "result_GSE160785/drug_immunesig_TUMERIC_GSEA_result_r0.4_meta_oe_original_en_padj1/COAD_Primary/drug989/CellvstreatedMouse.csv"
df = read.csv(file2)
df$treated_mouse_binary = sign(df$treated_mouse)
df$treated_cell_binary = sign(df$drugindex)
dfs = df[df$sign == 's', ]
dfs <- dfs[order(dfs$treated_mouse, decreasing = T),]
dfs <- dfs[order(dfs$treated_cell_binary, decreasing = T),]
dfs$immunesig <- factor(dfs$immunesig, levels = rev(dfs$immunesig ))
data_s <- melt(dfs[c("immunesig", "treated_mouse_binary", "treated_cell_binary")])
names(data_s) <-  c("immunesig_s", "variable", "value")
p2_s=ggplot(data_s, aes(x=variable,y=immunesig_s))+ 
  geom_tile(aes(fill=value), color = 'white', alpha=0.8)+ 
  scale_fill_gradient2(low = "#6A82FB",mid="grey", high = "#FC5C7D")+
  theme_bw()+
  ggtitle("GSE160785-COAD-drug989-celecoxib")

dfr = df[df$sign == 'r', ]
dfr <- dfr[order(dfr$treated_mouse, decreasing = T),]
dfr <- dfr[order(dfr$treated_cell_binary, decreasing = T),]
dfr$immunesig <- factor(dfr$immunesig, levels = dfr$immunesig )
data_r <- melt(dfr[c("immunesig", "treated_mouse_binary", "treated_cell_binary")])
names(data_r) <-  c("immunesig_r", "variable", "value")
p2_r=ggplot(data_r, aes(x=variable,y=immunesig_r))+ 
  geom_tile(aes(fill=value), color = 'white', alpha=0.8)+ 
  scale_fill_gradient2(low = "#6A82FB",mid="grey", high = "#FC5C7D")+
  theme_bw()+
  ggtitle("GSE160785-COAD-drug989-celecoxib")
p2 = p2_s/p2_r

# p1 + p2 + plot_layout(guides = 'auto')  + plot_layout(guides = 'collect')
ggsave("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/07_plotting_v2/04_cellVSmouse_metaPlusFilter/padj1/GSE160785_COAD.pdf", p2, height = 10, width = 5)

# v02	GSE160785 celecoxib---------------------------------------------------------
file2 = "result_GSE160785/drug_immunesig_TUMERIC_GSEA_result_r0.4_meta_oe_original_en_padj1/READ_Primary/drug989/CellvstreatedMouse.csv"
df = read.csv(file2)
df$treated_mouse_binary = sign(df$treated_mouse)
df$treated_cell_binary = sign(df$drugindex)
dfs = df[df$sign == 's', ]
dfs <- dfs[order(dfs$treated_mouse, decreasing = T),]
dfs <- dfs[order(dfs$treated_cell_binary, decreasing = T),]
dfs$immunesig <- factor(dfs$immunesig, levels = rev(dfs$immunesig ))
data_s <- melt(dfs[c("immunesig", "treated_mouse_binary", "treated_cell_binary")])
names(data_s) <-  c("immunesig_s", "variable", "value")
p2_s=ggplot(data_s, aes(x=variable,y=immunesig_s))+ 
  geom_tile(aes(fill=value), color = 'white', alpha=0.8)+ 
  scale_fill_gradient2(low = "#6A82FB",mid="grey", high = "#FC5C7D")+
  theme_bw()+
  ggtitle("GSE160785-READ-drug989-celecoxib")

dfr = df[df$sign == 'r', ]
dfr <- dfr[order(dfr$treated_mouse, decreasing = T),]
dfr <- dfr[order(dfr$treated_cell_binary, decreasing = T),]
dfr$immunesig <- factor(dfr$immunesig, levels = dfr$immunesig )
data_r <- melt(dfr[c("immunesig", "treated_mouse_binary", "treated_cell_binary")])
names(data_r) <-  c("immunesig_r", "variable", "value")
p2_r=ggplot(data_r, aes(x=variable,y=immunesig_r))+ 
  geom_tile(aes(fill=value), color = 'white', alpha=0.8)+ 
  scale_fill_gradient2(low = "#6A82FB",mid="grey", high = "#FC5C7D")+
  theme_bw()+
  ggtitle("GSE160785-READ-drug989-celecoxib")
p2 = p2_s/p2_r

# p1 + p2 + plot_layout(guides = 'auto')  + plot_layout(guides = 'collect')
ggsave("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/07_plotting_v2/04_cellVSmouse_metaPlusFilter/padj1/GSE160785_READ.pdf", p2, height = 10, width = 5)

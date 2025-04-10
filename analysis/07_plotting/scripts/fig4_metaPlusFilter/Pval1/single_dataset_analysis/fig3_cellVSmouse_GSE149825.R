# rm(list=ls())
setwd("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/")
options(stringsAsFactors = F)
library(ggplot2)
library(reshape2)
library(patchwork)
datasets <- list.files()
datasets <- datasets[grep("result_", datasets)]

dir.create("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/07_plotting_v2/04_cellVSmouse_metaPlusFilter/padj1/")


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
immunesig =immunesig[immunesig$weight != 0 ,]


# v02	GSE149825 birinapant---------------------------------------------------------
file2 = "result_GSE149825/drug_immunesig_TUMERIC_GSEA_result_r0.4_meta_oe_original_en_padj1/drug1487/CellvstreatedMouse.csv"
df = read.csv(file2)
names(df)[1] <- 'immune_sig'
df = merge(df, immunesig, by = 'immune_sig')

df$treated_mouse_binary = sign(df$treated_mouse)
df$treated_cell_binary = sign(df$drugindex)
dfs = df[df$sign == 's', ]
dfs <- dfs[order(dfs$treated_mouse, decreasing = T),]
dfs <- dfs[order(dfs$treated_cell_binary, decreasing = T),]
dfs$immune_sig <- factor(dfs$immune_sig, levels = rev(dfs$immune_sig ))
data_s <- melt(dfs[c('immune_sig', "treated_mouse_binary", "treated_cell_binary")])
names(data_s) <-  c("immunesig_s", "variable", "value")
p2_s=ggplot(data_s, aes(x=variable,y=immunesig_s))+ 
  geom_tile(aes(fill=value), color = 'white', alpha=0.8)+ 
  scale_fill_gradient2(low = "#6A82FB",mid="grey", high = "#FC5C7D")+
  theme_bw()+
  ggtitle("GSE149825-SKCM-drug1487-birinapant")

dfr = df[df$sign == 'r', ]
dfr <- dfr[order(dfr$treated_mouse, decreasing = T),]
dfr <- dfr[order(dfr$treated_cell_binary, decreasing = T),]
dfr$immune_sig <- factor(dfr$immune_sig, levels = dfr$immune_sig )
data_r <- melt(dfr[c('immune_sig', "treated_mouse_binary", "treated_cell_binary")])
names(data_r) <-  c("immunesig_r", "variable", "value")
p2_r=ggplot(data_r, aes(x=variable,y=immunesig_r))+ 
  geom_tile(aes(fill=value), color = 'white', alpha=0.8)+ 
  scale_fill_gradient2(low = "#6A82FB",mid="grey", high = "#FC5C7D")+
  theme_bw()+
  ggtitle("GSE149825-SKCM-drug1487-birinapant")
p2 = p2_s/p2_r

# p1 + p2 + plot_layout(guides = 'auto')  + plot_layout(guides = 'collect')
ggsave("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/07_plotting_v2/04_cellVSmouse_metaPlusFilter/padj1/GSE149825.pdf", p2, height = 10, width = 5)

names(data_s) <-  c("immunesig", "variable", "value")
names(data_r) <-  c("immunesig", "variable", "value")
data <- rbind(data_s, data_r)
p2_all=ggplot(data, aes(x=variable,y=immunesig))+ 
  geom_tile(aes(fill=value), color = 'white', alpha=0.8)+ 
  scale_fill_gradient2(low = "#6A82FB",mid="grey", high = "#FC5C7D")+
  theme_bw()+
  ggtitle("GSE149825-SKCM-drug1487-birinapant")
ggsave("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/07_plotting_v2/04_cellVSmouse_metaPlusFilter/padj1/GSE149825_2.pdf", p2_all, height = 5, width = 5)


# # v02	GSE149825 birinapant---------------------------------------------------------
# file2 = "result_GSE149825/drug_immunesig_TUMERIC_GSEA_result_r0.4_meta_oe_original/drug1487_SKCM/CellvstreatedMouse.csv"
# df = read.csv(file2)
# df$treated_mouse_binary = sign(df$treated_mouse)
# df$treated_cell_binary = sign(df$drugindex)
# dfs = df[df$sign == 's', ]
# dfs <- dfs[order(dfs$treated_mouse, decreasing = T),]
# dfs <- dfs[order(dfs$treated_cell_binary, decreasing = T),]
# dfs$immunesig <- factor(dfs$immunesig, levels = rev(dfs$immunesig ))
# data_s <- melt(dfs[c('immune_sig', "treated_mouse_binary", "treated_cell_binary")])
# names(data_s) <-  c("immunesig_s", "variable", "value")
# p2_s=ggplot(data_s, aes(x=variable,y=immunesig_s))+ 
#   geom_tile(aes(fill=value), color = 'white', alpha=0.8)+ 
#   scale_fill_gradient2(low = "#6A82FB",mid="grey", high = "#FC5C7D")+
#   theme_bw()+
#   ggtitle("GSE149825-SKCM-drug1487-birinapant")

# dfr = df[df$sign == 'r', ]
# dfr <- dfr[order(dfr$treated_mouse, decreasing = T),]
# dfr <- dfr[order(dfr$treated_cell_binary, decreasing = T),]
# dfr$immunesig <- factor(dfr$immunesig, levels = dfr$immunesig )
# data_r <- melt(dfr[c('immune_sig', "treated_mouse_binary", "treated_cell_binary")])
# names(data_r) <-  c("immunesig_r", "variable", "value")
# p2_r=ggplot(data_r, aes(x=variable,y=immunesig_r))+ 
#   geom_tile(aes(fill=value), color = 'white', alpha=0.8)+ 
#   scale_fill_gradient2(low = "#6A82FB",mid="grey", high = "#FC5C7D")+
#   theme_bw()+
#   ggtitle("GSE149825-SKCM-drug1487-birinapant")
# p2 = p2_s/p2_r

# # p1 + p2 + plot_layout(guides = 'auto')  + plot_layout(guides = 'collect')
# ggsave("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/07_plotting_v2/04_cellVSmouse_metaPlusFilter/padj1/GSE149825_SKCM.pdf", p2, height = 10, width = 5)


# v02	GSE149825 birinapant---------------------------------------------------------
file2 = "result_GSE149825/drug_immunesig_CPE_GSEA_result_r0.4_meta_oe_original_en_padj1/drug1487/CellvstreatedMouse.csv"
df = read.csv(file2)
df$treated_mouse_binary = sign(df$treated_mouse)
df$treated_cell_binary = sign(df$drugindex)
dfs = df[df$sign == 's', ]
dfs <- dfs[order(dfs$treated_mouse, decreasing = T),]
dfs <- dfs[order(dfs$treated_cell_binary, decreasing = T),]
dfs$immunesig <- factor(dfs$immunesig, levels = rev(dfs$immunesig ))
data_s <- melt(dfs[c('immunesig', "treated_mouse_binary", "treated_cell_binary")])
names(data_s) <-  c("immunesig_s", "variable", "value")
p2_s=ggplot(data_s, aes(x=variable,y=immunesig_s))+ 
  geom_tile(aes(fill=value), color = 'white', alpha=0.8)+ 
  scale_fill_gradient2(low = "#6A82FB",mid="grey", high = "#FC5C7D")+
  theme_bw()+
  ggtitle("GSE149825-pSKCM-drug1487-birinapant")

dfr = df[df$sign == 'r', ]
dfr <- dfr[order(dfr$treated_mouse, decreasing = T),]
dfr <- dfr[order(dfr$treated_cell_binary, decreasing = T),]
dfr$immunesig <- factor(dfr$immunesig, levels = dfr$immunesig )
data_r <- melt(dfr[c('immunesig', "treated_mouse_binary", "treated_cell_binary")])
names(data_r) <-  c("immunesig_r", "variable", "value")
p2_r=ggplot(data_r, aes(x=variable,y=immunesig_r))+ 
  geom_tile(aes(fill=value), color = 'white', alpha=0.8)+ 
  scale_fill_gradient2(low = "#6A82FB",mid="grey", high = "#FC5C7D")+
  theme_bw()+
  ggtitle("GSE149825-pSKCM-drug1487-birinapant")
p2 = p2_s/p2_r

# p1 + p2 + plot_layout(guides = 'auto')  + plot_layout(guides = 'collect')
ggsave("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/07_plotting_v2/04_cellVSmouse_metaPlusFilter/padj1/GSE149825_primarySKCM.pdf", p2, height = 10, width = 5)

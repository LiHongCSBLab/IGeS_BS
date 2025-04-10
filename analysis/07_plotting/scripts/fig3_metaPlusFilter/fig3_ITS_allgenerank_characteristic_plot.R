# rm(list=ls())
work_path <- "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/"
setwd(work_path)
options(stringsAsFactors = F)

library(patchwork)
library(ggplot2)
library(pheatmap)
library(aplot)
library(reshape2)
library(ggsci)
library(see)

library(dplyr)
library(UpSetR)
library(clusterProfiler)
library(enrichplot)
# ------------------------------------------------------------------------------
load( "07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITS_set.Rdata")

# ------------------------------------------------------------------------------

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

meta_selected$flag = 'sensitive'
meta_selected[meta_selected$OR < 1, ]$flag = 'resistant'
ITS_S = meta_selected[meta_selected$flag == 'sensitive', ]$X
ITS_R = meta_selected[meta_selected$flag == 'resistant', ]$X

# ITS_S = ITS_S[ITS_S != "merck18_tide"]


# ------------------------------------------------------------------------------

load("07_plotting_v2/data/immunesig_original_geneset.Rdata")
# immunesigGS, immunesigGS_df


immunesigGS_S <- immunesigGS_df[is.element(immunesigGS_df$immune_sig, ITS_S),]
immunesigGS_R <- immunesigGS_df[is.element(immunesigGS_df$immune_sig, ITS_R),]

# ------------------------------------------------------------------------------

ITS_p_s_cancerstat <- list()

for(i in 1:length(ITS_S)){
  # i = 1
  ITSname <- ITS_S[i]
  listinput <- lapply(ITS_p_set, function(x){x[[ITSname]]})
  genelist = unique(unlist(listinput))
  listinput_df <- fromList(listinput)
  rownames(listinput_df) <- genelist
  df_stat <- listinput_df # * 22
  df_stat_cancer <- apply(df_stat, 1, function(x)names(x)[which(x == 1)])
  df_stat_cancer <- lapply(as.list(names(df_stat_cancer)), function(x){
    data.frame(gene_name = x, cancer = df_stat_cancer[[x]])
  })
  ITS_p_s_cancerstat[[i]] <- do.call(rbind,df_stat_cancer)
}

df_stat_Cancer_s <- unique(do.call(rbind,ITS_p_s_cancerstat))
df_stat_NumCancer_s <- data.frame(table(df_stat_Cancer_s$gene_name))
names(df_stat_NumCancer_s) <- c("gene_name", "num_Cancer")

ITS_p_s = list()
ITS_p_s_genestat <- list()

for(i in 1:length(ITS_S)){
  # i = 1
  ITSname <- ITS_S[i]
  listinput <- lapply(ITS_p_set, function(x){x[[ITSname]]})
  genelist = unique(unlist(listinput))
  listinput_df <- fromList(listinput)
  rownames(listinput_df) <- genelist
  # table(apply(listinput_df,1,sum))
  # apply(listinput_df,2,sum)
  # ITS_p_s[[i]] = unique(names(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)>=2)]))
  ITS_p_s_genestat[[i]] <- listinput_df
  
  
  df_stat <- listinput_df # * 22
  df_stat_row <- rev(sort(rowSums(sign(df_stat))))
  df_stat_row <- as.data.frame(df_stat_row)
  names(df_stat_row) = "num_Gene"
  df_stat_row$gene_name = rownames(df_stat_row)
  df_stat_row$immune_sig <- ITSname
  
  
  ITS_p_s[[i]] <- df_stat_row # [df_stat_row$num_Gene > 1,]
  
}
clu
ITSs_cancer_stat = lapply(ITS_p_s, function(x){
  data.frame(ratio_cancer = nrow(x[x$num_Gene >=2,])/nrow(x))
})
names(ITSs_cancer_stat)=ITS_S
ITSs_cancer_stat = do.call(rbind, ITSs_cancer_stat)
ITSs_cancer_stat$immune_sig = rownames(ITSs_cancer_stat)
summary(ITSs_cancer_stat[ITSs_cancer_stat$ratio_cancer<1,])
write.csv(ITSs_cancer_stat, "07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/ITSs_geneInAbove2cancer_stat.csv", quote = F, row.names = F)

ITSs_stat = lapply(ITS_p_set, function(x){
  # x=ITS_p_set[[1]]
  x = x[ITS_S]
  y = data.frame(table(unlist(x)))
  data.frame(ratio_ITS = nrow(y[y$Freq >=2,])/nrow(y))
})
names(ITSs_stat)=names(ITS_p_set)
ITSs_stat = do.call(rbind, ITSs_stat)
ITSs_stat$cancer = rownames(ITSs_stat)
summary(ITSs_stat[ITSs_stat$ratio_ITS<1,])
write.csv(ITSs_stat, "07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/ITSs_geneInAbove2ITS_stat.csv", quote = F, row.names = F)


ITS_p_s_df <- do.call(rbind, ITS_p_s)
ITS_p_s_df2 <- unique(data.frame(table(ITS_p_s_df$gene_name)))
names(ITS_p_s_df2) <- c("gene_name", "num_ITS")

# ------------------------------------------------------------------------------
ITS_p_s_df3 <- read.table("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/ITS_p_s_df.txt", header = T, sep = '\t')

dat_s <- merge(df_stat_NumCancer_s, ITS_p_s_df2,by="gene_name")
dat_s <- merge(unique(ITS_p_s_df3[c('gene_name', 'weighted_rate')]), dat_s, by="gene_name")
dat_s$density_Cancer <- dat_s$num_Cancer / 22
dat_s$density_ITS <- dat_s$num_ITS / length(ITS_S)
dat_s$weighted_rate_zscale = scale(dat_s$weighted_rate)

library(ggplot2)
library(dplyr)
library(viridis) # 使用viridis提供的翠绿色标度：scale_fill_viridis()
# library(ggpointdensity) # 绘制密度散点图

p_s <- ggplot(data = dat_s, mapping = aes(x = density_Cancer,
                                       y = density_ITS,
                                       color = weighted_rate)) + 
  geom_point() + #散点图
  scale_color_viridis(option="C", begin = 0, end = 1)  +
  scale_y_continuous(breaks = seq(0,1,0.2))+ 
  xlab("density_Cancer") + 
  theme(axis.title.x = element_text(size = 16,
                                    face = "bold", 
                                    vjust = 0.5, 
                                    hjust = 0.5))+
  ylab("density_ITS") + 
  theme(axis.title.y = element_text(size = 16,
                                    face = "bold", 
                                    vjust = 0.5, 
                                    hjust = 0.5))+
  theme_bw()+
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        text=element_text(size=12,  family="serif")) +
  ggtitle("Gene statistics in sensitive ITS")
        
p_s

# ------------------------------------------------------------------------------
ITS_p_r_cancerstat <- list()

for(i in 1:length(ITS_R)){
  # i = 1
  ITSname <- ITS_R[i]
  listinput <- lapply(ITS_p_set, function(x){x[[ITSname]]})
  genelist = unique(unlist(listinput))
  listinput_df <- fromList(listinput)
  rownames(listinput_df) <- genelist
  df_stat <- listinput_df # * 22
  df_stat_cancer <- apply(df_stat, 1, function(x)names(x)[which(x == 1)])
  df_stat_cancer <- lapply(as.list(names(df_stat_cancer)), function(x){
    data.frame(gene_name = x, cancer = df_stat_cancer[[x]])
  })
  ITS_p_r_cancerstat[[i]] <- do.call(rbind,df_stat_cancer)
}

df_stat_Cancer_r <- unique(do.call(rbind,ITS_p_r_cancerstat))
df_stat_NumCancer_r <- data.frame(table(df_stat_Cancer_r$gene_name))
names(df_stat_NumCancer_r) <- c("gene_name", "num_Cancer")



ITS_p_r = list()
ITS_p_r_genestat <- list()

for(i in 1:length(ITS_R)){
  # i = 1
  ITSname <- ITS_R[i]
  listinput <- lapply(ITS_p_set, function(x){x[[ITSname]]})
  genelist = unique(unlist(listinput))
  listinput_df <- fromList(listinput)
  rownames(listinput_df) <- genelist

  ITS_p_r_genestat[[i]] <- listinput_df
  
  
  df_stat <- listinput_df # * 22
  df_stat_row <- rev(sort(rowSums(sign(df_stat))))
  df_stat_row <- as.data.frame(df_stat_row)
  names(df_stat_row) = "num_Gene"
  df_stat_row$gene_name = rownames(df_stat_row)
  df_stat_row$immune_sig <- ITSname
  
  
  ITS_p_r[[i]] <- df_stat_row # [df_stat_row$num_Gene > 1,]
  
}

ITSr_cancer_stat = lapply(ITS_p_r, function(x){
  data.frame(ratio_cancer = nrow(x[x$num_Gene >=2,])/nrow(x))
})
names(ITSr_cancer_stat)=ITS_R
ITSr_cancer_stat = do.call(rbind, ITSr_cancer_stat)
ITSr_cancer_stat$immune_sig = rownames(ITSr_cancer_stat)
summary(ITSr_cancer_stat[ITSr_cancer_stat$ratio_cancer<1,])
write.csv(ITSr_cancer_stat, "07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/ITSr_geneInAbove2cancer_stat.csv", quote = F, row.names = F)


ITSr_stat = lapply(ITS_p_set, function(x){
  # x=ITS_p_set[[1]]
  x = x[ITS_R]
  y = data.frame(table(unlist(x)))
  data.frame(ratio_ITS = nrow(y[y$Freq >=2,])/nrow(y))
})
names(ITSr_stat)=names(ITS_p_set)
ITSr_stat = do.call(rbind, ITSr_stat)
ITSr_stat$cancer = rownames(ITSr_stat)
summary(ITSr_stat[ITSr_stat$ratio_ITS<1,])
write.csv(ITSr_stat, "07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/ITSr_geneInAbove2ITS_stat.csv", quote = F, row.names = F)

ITS_p_r_df <- do.call(rbind, ITS_p_r)
ITS_p_r_df2 <- unique(data.frame(table(ITS_p_r_df$gene_name)))
names(ITS_p_r_df2) <- c("gene_name", "num_ITS")

# ------------------------------------------------------------------------------
ITS_p_r_df3 <- read.table("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/ITS_p_r_df.txt", header = T, sep = '\t')

dat_r <- merge(df_stat_NumCancer_r, ITS_p_r_df2,by="gene_name")
dat_r <- merge(unique(ITS_p_r_df3[c('gene_name', 'weighted_rate')]), dat_r, by="gene_name")
dat_r$density_Cancer <- dat_r$num_Cancer / 22
dat_r$density_ITS <- dat_r$num_ITS / length(ITS_R)
dat_r$weighted_rate_zscale = scale(dat_r$weighted_rate)

library(ggplot2)
library(dplyr)
library(viridis) # 使用viridis提供的翠绿色标度：scale_fill_viridis()
# library(ggpointdensity) # 绘制密度散点图

p_r <- ggplot(data = dat_r, mapping = aes(x = density_Cancer,
                                       y = density_ITS,
                                       color = weighted_rate)) + 
  geom_point() + #散点图
  scale_color_viridis(option="C", begin = 0, end = 1)  +
  scale_y_continuous(breaks = seq(0,1,0.2))+ 
  xlab("density_Cancer") + 
  theme(axis.title.x = element_text(size = 16,
                                    face = "bold", 
                                    vjust = 0.5, 
                                    hjust = 0.5))+
  ylab("density_ITS") + 
  theme(axis.title.y = element_text(size = 16,
                                    face = "bold", 
                                    vjust = 0.5, 
                                    hjust = 0.5))+
  theme_bw()+
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        text=element_text(size=12,  family="serif")) +
  ggtitle("Gene statistics in resistant ITS")
        
p_r
# p_dot = p_s %>%
#   insert_bottom(p_r)
#   # insert_right(p_r)

# ------------------------------------------------------------------------------
# rank plot according to weighted_rate
library(ggrepel)
# dat_s
dat_s <- dat_s[order(dat_s$weighted_rate, decreasing = T), ]
dat_s$x = 1:nrow(dat_s)
dat_s$hits <- 'medium'
dat_s[dat_s$x < 30, ]$hits <- 'top30'
dat_s = dat_s %>%
  arrange(desc(weighted_rate)) %>%
  mutate("x" = row_number())

dat_s_shown <- dat_s[dat_s$x < 30, ]

p_s_rank <- ggplot(dat_s, aes(x = x, y = weighted_rate, color=hits))+
  geom_point(color="#B33333", size = 0.5)+
  theme_bw()+
  labs(x = "Ranked Genes", y = "Importance")+
  geom_text_repel(inherit.aes = F, 
                  data = dat_s_shown, 
                  aes(x = x, y = weighted_rate, label = gene_name, color = hits),
                  max.overlaps = getOption("ggrepel.max.overlaps", default = 30))+
scale_y_continuous(breaks = seq(0,1,0.2))



# dat_r
dat_r <- dat_r[order(dat_r$weighted_rate, decreasing = T), ]
dat_r$x = 1:nrow(dat_r)
dat_r$hits <- 'medium'
dat_r[dat_r$x < 30, ]$hits <- 'top30'
dat_r = dat_r %>%
  arrange(desc(weighted_rate)) %>%
  mutate("x" = row_number())

dat_r_shown <- dat_r[dat_r$x < 30, ]

p_r_rank <- ggplot(dat_r, aes(x = x, y = weighted_rate, color=hits))+
  geom_point(color="#338080", size = 0.5)+
  theme_bw()+
  labs(x = "Ranked Genes", y = "Importance")+
  geom_text_repel(inherit.aes = F, 
                  data = dat_r_shown, 
                  aes(x = x, y = weighted_rate, label = gene_name, color = hits),
                  max.overlaps = getOption("ggrepel.max.overlaps", default = 30))+
scale_y_continuous(breaks = seq(0,1,0.2))


p_s2 = p_s %>%
  insert_right(p_s_rank, width=0.5)
p_r2 =  p_r %>%
  insert_right(p_r_rank, width=0.5)

ggsave("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/generank_s.pdf", p_s2, height = 8, width = 12.5)
ggsave("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/generank_r.pdf", p_r2, height = 8, width = 12.5)

p_s2 =p_s / p_s_rank

p_r_rank <- ggplot(dat_r, aes(x = x, y = weighted_rate, color=hits))+
  geom_point(color="#338080", size = 0.5)+
  theme_bw()+
  labs(x = "Ranked Genes", y = "Importance")+
  geom_text_repel(inherit.aes = F, 
                  data = dat_r_shown, 
                  aes(x = x, y = weighted_rate, label = gene_name, color = hits),
                  max.overlaps = getOption("ggrepel.max.overlaps", default = 30))+
scale_y_continuous(breaks = seq(0,0.2,0.1))
p_r2 =p_r / p_r_rank
ggsave("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/generank_s_v2.pdf", p_s2, height = 8, width = 6)
ggsave("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/generank_r_v2.pdf", p_r2, height = 8, width = 6)


# ------------------------------------------------------------------------------
# GSEA enrichment result plot

# load("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/enrich/custom_GSEA_s_scaled.Rdata")
# res_s <- res
# load("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/enrich/custom_GSEA_r_scaled.Rdata")
# res_r <- res

# ridgeplot(res_s)
# gseaplot2(res_s, geneSetID = res_s$ID)

# ridgeplot(res_r)
# gseaplot2(res_r, geneSetID = res_r$ID)

# ------------------------------------------------------------------------------
res_s <- read.table("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/enrich/custom_GSEA_s_scaled.txt", sep = '\t', header = T)
res_r <- read.table("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/enrich/custom_GSEA_r_scaled.txt", sep = '\t', header = T)

GS_selected <- c("ITSp_vs_Sensitizer_PMID34980132",
                        "ITSp_vs_Resistor_PMID34980132",
                        "ITSp_vs_icb_enhancer_PMID34980132",
                        "ITSp_vs_icb_suppressor_PMID34980132",
                        "ITSp_vs_stimulatory_PMID29628290",
                        "ITSp_vs_inhibitor_PMID29628290",
                        "ITSp_vs_co_stim_inhib_PMID29628290",
                        "ITSp_vs_regulator_p_PMID34980132",
                        "ITSp_vs_regulator_n_PMID34980132")

res_s2 <- res_s[is.element(res_s$ID, GS_selected), ]
res_s2$type <- 's'
res_r2 <- res_r[is.element(res_r$ID, GS_selected), ]
res_r2$type <- 'r'

res = rbind(res_s2, res_r2)
res$signif <- ""
res[res$pvalue < 0.05, ]$signif =  "*"
# res[res$p.adjust < 0.2, ]$signif =  "*"
res[res$p.adjust < 0.1, ]$signif =  "**"
res[res$p.adjust < 0.01, ]$signif =  "***"

res$type <- factor(res$type, levels = c('s','r'))
res$ID <- factor(res$ID, levels = rev(GS_selected))


library(stringr)
res$ID <- gsub("pmid","PMID", gsub("Icb","ICB", str_to_title(gsub("ITSp_vs_", "", res$ID))))
p_heatmap = ggplot(res, aes(type, ID)) + 
  geom_tile(aes(fill = NES), colour = "grey", 
            size = 1,
            lwd = 2,
            linetype = 1)+
  scale_fill_gradient2(low = "#5C5DAF",mid = "white",high = "#EA2E2D", midpoint = 0) + 
  geom_text(aes(label=signif),col ="black",size = 8) +
  theme_bw() + 
  theme(axis.title.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.title.y=element_blank(), 
        axis.text.x = element_text(angle = 90, hjust = 1, size = 16), 
        axis.text.y = element_text(size = 12)) +
  labs(color = " * p < 0.05\n ** padj < 0.1\n *** padj < 0.01") 

ggsave("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/GSenrichment.pdf", p_heatmap, height = 5, width = 5)

p_bubble = ggplot(res, aes(type, ID)) + 
  geom_point(aes(color=NES, size=-log10(p.adjust)))+ # 
  scale_colour_gradient2(low="#5C5DAF",mid = 'white',high="#EA2E2D")+
  geom_text(aes(label=signif),col ="black",size = 8) +
  theme_bw() + 
  theme(axis.title.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.title.y=element_blank(), 
        axis.text.x = element_text(angle = 90, hjust = 1, size = 16), 
        axis.text.y = element_text(size = 12)) +
  labs(color =" * p < 0.05\n ** padj < 0.1\n *** padj < 0.01") 

ggsave("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/GSenrichment_bubble.pdf", p_bubble, height = 5, width = 5)


# rm(list=ls())
work_path <- "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/"
setwd(work_path)
options(stringsAsFactors = F)

# library(reshape2)
# library(ggplotify)
# library(grid)
# library(cowplot)
# library(parallel)
# library(patchwork)
# library(dplyr)
# library(stringr)
# library(data.table)
# require(GGally)
# library(plot3D)
# library(fmsb)
# library(readxl)
library(ggplot2)
library(dplyr)
library(UpSetR)
library(clusterProfiler)
library(enrichplot)
# source("06_results_summaryandplotting/scripts/plot_for_drug_function.R")
# source("06_results_summaryandplotting/scripts/ITS_analysis/ITS_gene_stat_function.R")
# source("06_results_summaryandplotting/scripts/immunesig_selection_function.R")
# source("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/scripts/drugITGSscore_lincs.R")
source("07_plotting_v2/scripts/functions/ITS_enrich.R")
source("07_plotting_v2/scripts/functions/ITS_function_detection.R")
source("06_results_summaryandplotting/scripts/immunesig_selection_function.R")
source("07_plotting_v2/scripts/functions/ITS_enrich.R")
source("07_plotting_v2/scripts/functions/ITS_function_detection.R")
source("07_plotting_v2/scripts/functions/ITS_g_cancer_plot_function.R")

dir.create("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/DT_enrich/")
# ------------------------------------------------------------------------------
# Drug target analysis
# ------------------------------------------------------------------------------
load( "07_plotting_v2/data/drugTarget_merged3DB.Rdata")

# ------------------------------------------------------------------------------
gs_dt = rbind(data.frame(gs_name = "drugTarget_all", SYMBOL = unique(drugTarget_all$genename)),
              data.frame(gs_name = "drugTarget_FDA", SYMBOL = unique(drugTarget_FDA$genename)),
              data.frame(gs_name = "drugTarget_FDAcancer", SYMBOL = unique(drugTarget_FDA[which(drugTarget_FDA$ifCancer_drugbank == 'yes'),]$genename)))

# gs_all = unique(rbind(gs_all, gs_dt))

ITS_p_s_df <- read.table("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/ITS_p_s_df.txt", sep = '\t', header = T)
tmp2 = as.vector(scale(ITS_p_s_df$weighted_rate))
names(tmp2) = ITS_p_s_df$gene_name

res = ITS_GS_enrich(tmp2,
                    gs=gs_dt,
                    pvalueCutoff = 1,
                    type = "custom",
                    method = 'GSEA')
data.frame(res)[c(1,6,7)]
res2=data.frame(res)
write.table(res2, "07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/DT_enrich/sITS_DT_enrich.txt", sep = '\t', quote = F, row.names = F)

p_s <- gseaplot2(res, geneSetID = c(3,1), pvalue_table = TRUE, subplots = 1:2)
p_s1 <- gseaplot2(res, geneSetID = c(2,3,1), pvalue_table = TRUE, subplots = 1:2)

p1 <- gseaplot(res, geneSetID = 2, by = "runningScore", title = res$Description[2])
p2 <- gseaplot(res, geneSetID = 3, by = "runningScore", title = res$Description[3])
p3 <- gseaplot(res, geneSetID = 1, by = "runningScore", title = res$Description[1])
p_s2 <- cowplot::plot_grid(p1, p2, p3, ncol=1, labels=LETTERS[1:3])

pp <- lapply(c(2,3,1), function(i) {
    anno <- res[i, c("NES", "pvalue", "p.adjust")]
    lab <- paste0(names(anno), "=",  round(anno, 3), collapse="\n")

    gsearank(res, i, res[i, 2]) + xlab(NULL) +ylab(NULL) +
        annotate("text", 10000, res[i, "enrichmentScore"] * .75, label = lab, hjust=0, vjust=0)
})
p_s3 <- plot_grid(plotlist=pp, ncol=1)

pdf("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/DT_enrich/sITS_DT_enrich.pdf", onefile = TRUE, width = 8, height = 8)
print(p_s)
dev.off()
pdf("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/DT_enrich/sITS_DT_enrich1.pdf", onefile = TRUE, width = 8, height = 8)
print(p_s1)
dev.off()
pdf("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/DT_enrich/sITS_DT_enrich2.pdf", onefile = TRUE, width = 10, height = 10)
print(p_s2)
dev.off()
pdf("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/DT_enrich/sITS_DT_enrich3.pdf", onefile = TRUE, width = 15, height = 10)
print(p_s3)
dev.off()



ITS_p_r_df <- read.table("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/ITS_p_r_df.txt", sep = '\t', header = T)
tmp2 = as.vector(scale(ITS_p_r_df$weighted_rate))
names(tmp2) = ITS_p_r_df$gene_name

res = ITS_GS_enrich(tmp2,
                    gs=gs_dt,
                    pvalueCutoff = 1,
                    type = "custom",
                    method = 'GSEA')
data.frame(res)[c(1,6,7)]
res2=data.frame(res)
write.table(res2, "07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/DT_enrich/rITS_DT_enrich.txt", sep = '\t', quote = F, row.names = F)

p_r <- gseaplot2(res, geneSetID = c(3,2), pvalue_table = TRUE, subplots = 1:2)
p_r1 <- gseaplot2(res, geneSetID = c(2,3,1), pvalue_table = TRUE, subplots = 1:2)

p1 <- gseaplot(res, geneSetID = 2, by = "runningScore", title = res$Description[2])
p2 <- gseaplot(res, geneSetID = 3, by = "runningScore", title = res$Description[3])
p3 <- gseaplot(res, geneSetID = 1, by = "runningScore", title = res$Description[1])
p_r2 <- cowplot::plot_grid(p1, p2, p3, ncol=1, labels=LETTERS[1:3])

pp <- lapply(c(2,3,1), function(i) {
    anno <- res[i, c("NES", "pvalue", "p.adjust")]
    lab <- paste0(names(anno), "=",  round(anno, 3), collapse="\n")

    gsearank(res, i, res[i, 2]) + xlab(NULL) +ylab(NULL) +
        annotate("text", 10000, res[i, "enrichmentScore"] * .75, label = lab, hjust=0, vjust=0)
})
p_r3 <- plot_grid(plotlist=pp, ncol=1)


pdf("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/DT_enrich/rITS_DT_enrich.pdf", onefile = TRUE, width = 8, height = 8)
print(p_r)
dev.off()
pdf("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/DT_enrich/rITS_DT_enrich1.pdf", onefile = TRUE, width = 8, height = 8)
print(p_r1)
dev.off()
pdf("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/DT_enrich/rITS_DT_enrich2.pdf", onefile = TRUE, width = 10, height = 10)
print(p_r2)
dev.off()
pdf("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/DT_enrich/rITS_DT_enrich3.pdf", onefile = TRUE, width = 15, height = 10)
print(p_r3)
dev.off()


# ------------------------------------------------------------------------------
# VENN plot
ITS_p_s_df <- read.table("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/ITS_p_s_df.txt", sep = '\t', header = T)
dt_s <- list(#drugTarget_all=unique(drugTarget_all$genename),
             drugTarget_FDA=unique(drugTarget_FDA$genename),
             drugTarget_FDAcancer=unique(drugTarget_FDA[which(drugTarget_FDA$ifCancer_drugbank == 'yes'),]$genename),
             ITS_s_gene=unique(ITS_p_s_df$gene_name))

dt_s <- lapply(dt_s, function(x)if(length(which(is.na(x)))>0){x[-which(is.na(x))]}else{x})

library(VennDiagram)
grid.newpage()
venn.plot <- venn.diagram(
  dt_s,
  filename = NULL,"07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/drug_merge3DB/DT_s_venn.tiff",
  col = "black",
  lty = "dotted",
  lwd = 4,
  fill = c("cornflowerblue", "green", "darkorchid1"), #"yellow", 
  alpha = 0.50,
  # label.col = c("orange", "white", "darkorchid4", "white", "white", "white",
  #               "white", "white", "darkblue", "white",
  #               "white", "white", "white", "darkgreen", "white"),
  cex = 1.5,
  fontfamily = "serif",
  fontface = "bold",
  cat.col =c("cornflowerblue", "green", "darkorchid1"), # c("darkblue", "darkgreen", "orange", "darkorchid4"),
  cat.cex = 1.5,
  cat.fontfamily = "serif")
dev.off()

pdf("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/drug_merge3DB/DT_s_venn.pdf", width = 5, height = 4)
grid.draw(venn.plot)
dev.off()

num_Gene = 20000
num_sITS = 3736+487+151
num_DT = 487+151+239+1816 

num_sITS_DT = 487+151
num_sITS_DTcancer = 151 
num_DTcancer = 151+239 


phyper(num_sITS_DT-1, num_DT, num_Gene-num_DT, num_sITS, lower.tail=F)
phyper(num_sITS_DTcancer-1, num_DTcancer, num_Gene-num_DTcancer, num_sITS, lower.tail=F)


# library(ggVennDiagram)
# # 可视化绘制
# library(ggplot2)
# ggVennDiagram(dt_s, category.names = c("drugTarget_all","drugTarget_FDA","drugTarget_FDAcancer", "ITS_s_gene"),
#               size=1,lty="longdash",color="gray60") + 
#   scale_fill_gradient(name="Count",low="#EC7D85",high = "#182F6F") +
#   hrbrthemes::theme_ipsum(base_family = "sans") +
#   labs(title = "Drug Targets in sensitive ITS") +
#   theme(plot.title = element_text(hjust = 0.5,vjust = .5,color = "black",face = 'bold',
#                                   size = 20, margin = margin(t = 1, b = 12)),
#         plot.subtitle = element_text(hjust = 0,vjust = .5,size=15),
#         plot.caption = element_text(face = 'bold',size = 12),
#         axis.text.x = element_blank(),
#         axis.text.y = element_blank(),
#         axis.title.x = element_blank(),
#         axis.title.y = element_blank())
# ------------------------------------------------------------------------------
ITS_p_r_df <- read.table("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/ITS_p_r_df.txt", sep = '\t', header = T)
dt_r <- list( # drugTarget_all=unique(drugTarget_all$genename),
             drugTarget_FDA=unique(drugTarget_FDA$genename),
             drugTarget_FDAcancer=unique(drugTarget_FDA[which(drugTarget_FDA$ifCancer_drugbank == 'yes'),]$genename),
             ITS_r_gene=unique(ITS_p_r_df$gene_name))

dt_r <- lapply(dt_r, function(x)if(length(which(is.na(x)))>0){x[-which(is.na(x))]}else{x})

library(VennDiagram)
grid.newpage()
venn.plot <- venn.diagram(
  dt_r,
  filename = NULL, # "07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/drug_merge3DB/DT_r_venn.tiff",
  col = "black",
  lty = "dotted",
  lwd = 4,
  fill = c("cornflowerblue", "green", "darkorchid1"),  #c("cornflowerblue", "green", "yellow", "darkorchid1"),
  alpha = 0.50,
  # label.col = c("orange", "white", "darkorchid4", "white", "white", "white",
  #               "white", "white", "darkblue", "white",
  #               "white", "white", "white", "darkgreen", "white"),
  # cex = 1.5,
  fontfamily = "serif",
  fontface = "bold",
  cat.col = c("cornflowerblue", "green", "darkorchid1"),  #c("darkblue", "darkgreen", "orange", "darkorchid4"),
  cat.cex = 1.5,
  cat.fontfamily = "serif")
dev.off()

pdf("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/drug_merge3DB/DT_r_venn.pdf", width = 5, height = 4)
grid.draw(venn.plot)
dev.off()

# ------------------------------------------------------------------------------
# drug stat
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# ITS_p_s_df
# d0 = length(unique(drugTarget_all$drugname_lower))
# dt1 = intersect(dt_s[['drugTarget_all']], dt_s[['ITS_s_gene']])
# d1 = length(unique(drugTarget_all[which(is.element(drugTarget_all$genename,dt1)),]$drugname_lower))

dt2 = intersect(dt_s[['drugTarget_FDA']], dt_s[['ITS_s_gene']])
d2 = length(unique(drugTarget_FDA[which(is.element(drugTarget_FDA$genename,dt2)),]$drugname_lower))

drugTarget_FDAcancer = drugTarget_FDA[drugTarget_FDA$ifCancer_drugbank=='yes', ]
dt3 = intersect(dt_s[['drugTarget_FDAcancer']], dt_s[['ITS_s_gene']])
d3 = length(unique(drugTarget_FDAcancer[which(is.element(drugTarget_FDAcancer$genename,dt2)),]$drugname_lower))


df_s <- data.frame(group = c("num_drugFDA","num_drugFDAcancer"), #"num_compound_all","num_compound",
                   value = c(d2,d3), # d0, d1,
                   type = 's')
df_s$percent <- df_s$value / df_s[df_s$group == 'num_drugFDA', ]$value
df_s$plotvalue = df_s$percent
# df_s[df_s$group == 'num_compound_all',]$plotvalue <- 1-df_s[df_s$group == 'num_compound', ]$percent
# df_s[df_s$group == 'num_compound',]$plotvalue <- df_s[df_s$group == 'num_compound', ]$percent - df_s[df_s$group == 'num_drugFDA', ]$percent
df_s[df_s$group == 'num_drugFDA',]$plotvalue <- df_s[df_s$group == 'num_drugFDA', ]$percent - df_s[df_s$group == 'num_drugFDAcancer', ]$percent


# ITS_p_r_df
# d0 = length(unique(drugTarget_all$drugname_lower))
# dt1 = intersect(dt_r[['drugTarget_all']], dt_r[['ITS_r_gene']])
# d1 = length(unique(drugTarget_all[which(is.element(drugTarget_all$genename,dt1)),]$drugname_lower))

dt2 = intersect(dt_r[['drugTarget_FDA']], dt_r[['ITS_r_gene']])
d2 = length(unique(drugTarget_FDA[which(is.element(drugTarget_FDA$genename,dt2)),]$drugname_lower))

drugTarget_FDAcancer = drugTarget_FDA[drugTarget_FDA$ifCancer_drugbank=='yes', ]
dt3 = intersect(dt_r[['drugTarget_FDAcancer']], dt_r[['ITS_r_gene']])
d3 = length(unique(drugTarget_FDAcancer[which(is.element(drugTarget_FDAcancer$genename,dt2)),]$drugname_lower))

df_r <- data.frame(group = c("num_drugFDA","num_drugFDAcancer"), # "num_compound_all","num_compound",
                   value = c(d2,d3), # d0, d1,
                   type = 'r')
df_r$percent <- df_r$value / df_r[df_r$group == 'num_drugFDA', ]$value
df_r$plotvalue = df_r$percent
# df_r[df_r$group == 'num_compound_all',]$plotvalue <- 1-df_r[df_r$group == 'num_compound', ]$percent
# df_r[df_r$group == 'num_compound',]$plotvalue <- df_r[df_r$group == 'num_compound', ]$percent - df_r[df_r$group == 'num_drugFDA', ]$percent
df_r[df_r$group == 'num_drugFDA',]$plotvalue <- df_r[df_r$group == 'num_drugFDA', ]$percent - df_r[df_r$group == 'num_drugFDAcancer', ]$percent

df_stat <- rbind(df_s, df_r)
df_stat$type <- factor(df_stat$type, levels = c('s','r'))
df_stat2 <- df_stat

p = ggplot(df_stat, aes( x = type, y=100 * plotvalue,fill = group))+
  #geom_col和geom_bar这两条命令都可以绘制堆叠柱形图
  geom_bar(position = "stack", stat = "identity", width = 0.6)+
  geom_text(aes(label=paste0(value,'\n(', 100 * round(plotvalue,4), '%)')),vjust=-2,size=4,color="black")+ 
  scale_colour_manual(values = c('green','yellow','red')) +
  theme_bw() + # 不要背景
  theme(axis.title.x=element_blank(), # 去掉 title
        axis.ticks.x=element_blank(), # 去掉x 轴
        axis.title.y=element_blank(), # 去掉 y 轴
        axis.text.x = element_text(angle = 90, hjust = 1, size = 18, face = "bold"), # 调整x轴文字，字体加粗
        axis.text.y = element_text(size = 14, face = "bold"))
p
ggsave("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/drug_merge3DB/drug_stat2.pdf", p, height = 5, width = 5)



# ------------------------------------------------------------------------------
# GSEA  ------------------------------------------------------------------------
# ------------------------------------------------------------------------------
filepath <- ("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/drug_merge3DB/GSEAenrich/ITSp_s_vs_drugtarget")
DT_drug_s <- read.table(paste0(filepath, "_stat.txt"), sep = '\t',header = T)
DT_drug_s$type = 's'
DT_drug_s$percent <- DT_drug_s$value / DT_drug_s[DT_drug_s$group == 'num_drug', ]$value
DT_drug_s$plotvalue = DT_drug_s$percent
DT_drug_s[DT_drug_s$group == 'num_drug',]$plotvalue <- 1-DT_drug_s[DT_drug_s$group == 'num_drugFDA', ]$percent
DT_drug_s[DT_drug_s$group == 'num_drugFDA',]$plotvalue <- DT_drug_s[DT_drug_s$group == 'num_drugFDA', ]$percent - DT_drug_s[DT_drug_s$group == 'num_drugFDAcancer', ]$percent

 
filepath <- ("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/drug_merge3DB/GSEAenrich/ITSp_r_vs_drugtarget")
DT_drug_r <- read.table(paste0(filepath, "_stat.txt"), sep = '\t',header = T)
DT_drug_r$type = 'r'
DT_drug_r$percent <- DT_drug_r$value / DT_drug_r[DT_drug_r$group == 'num_drug', ]$value
DT_drug_r$plotvalue = DT_drug_r$percent
DT_drug_r[DT_drug_r$group == 'num_drug',]$plotvalue <- 1-DT_drug_r[DT_drug_r$group == 'num_drugFDA', ]$percent
DT_drug_r[DT_drug_r$group == 'num_drugFDA',]$plotvalue <- DT_drug_r[DT_drug_s$group == 'num_drugFDA', ]$percent - DT_drug_r[DT_drug_r$group == 'num_drugFDAcancer', ]$percent

DT_drug_stat <- rbind(DT_drug_s, DT_drug_r)
DT_drug_stat$type <- factor(DT_drug_stat$type, levels = c('s','r'))
DT_drug_stat2 <- DT_drug_stat

p = ggplot(DT_drug_stat, aes( x = type, y=100 * plotvalue,fill = group))+
  #geom_col和geom_bar这两条命令都可以绘制堆叠柱形图
  geom_bar(position = "stack", stat = "identity", width = 0.6)+
  geom_text(aes(label=paste0(value,'\n(', 100 * round(plotvalue,4), '%)')),vjust=1,size=4,color="black")+ 
  scale_colour_manual(values = c('green','yellow','red')) +
  theme_bw() + # 不要背景
  theme(axis.title.x=element_blank(), # 去掉 title
        axis.ticks.x=element_blank(), # 去掉x 轴
        axis.title.y=element_blank(), # 去掉 y 轴
        axis.text.x = element_text(angle = 90, hjust = 1, size = 18, face = "bold"), # 调整x轴文字，字体加粗
        axis.text.y = element_text(size = 14, face = "bold"))

ggsave("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/drug_merge3DB/GSEAenrich/DT_drug_stat.pdf", p, height = 5, width = 7)




# ------------------------------------------------------------------------------
drugTarget_FDAcancer = drugTarget_FDA[drugTarget_FDA$ifCancer_drugbank=='yes', ]

res_s <- read.table("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/DT_enrich/sITS_DT_enrich.txt", sep = '\t', head = T)
res_s_core = strsplit(res_s$core_enrichment, split = '/')
names(res_s_core) = res_s$ID
sITS_FDA <- drugTarget_FDA[is.element(drugTarget_FDA$genename, res_s_core[[3]]),]
length(unique(sITS_FDA$drugname))
sort(table(sITS_FDA$genename))

sITS_FDAcancer <- drugTarget_FDAcancer[is.element(drugTarget_FDAcancer$genename, res_s_core[[1]]),]
length(unique(sITS_FDAcancer$drugname))
sort(table(sITS_FDAcancer$genename))

sITS_DTcancer_kegg <- ITS_GS_enrich(unique(sITS_FDAcancer$genename),
                                    pvalueCutoff = 0.05,
                                    type = "KEGG",
                                    method = 'ORA',
                                    savepath ="07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/DT_enrich/sITS_DTcancer_kegg")
data.frame(sITS_DTcancer_kegg)


length(intersect(tolower(unique(sITS_FDAcancer$drugname)), dat_s_shown$drugname))

res_r <- read.table("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/DT_enrich/rITS_DT_enrich.txt", sep = '\t', head = T)
res_r_core = strsplit(res_r$core_enrichment, split = '/')
names(res_r_core) = res_r$ID
rITS_FDA <- drugTarget_FDA[is.element(drugTarget_FDA$genename, res_r_core[[3]]),]
length(unique(rITS_FDA$drugname))
sort(table(rITS_FDA$genename))


rITS_FDAcancer <- drugTarget_FDAcancer[is.element(drugTarget_FDAcancer$genename, res_r_core[[1]]),]
length(unique(rITS_FDAcancer$drugname))
sort(table(rITS_FDAcancer$genename))
rITS_DTcancer_kegg <- ITS_GS_enrich(unique(rITS_FDAcancer$genename),
                                    pvalueCutoff = 0.05,
                                    type = "KEGG",
                                    method = 'ORA',
                                    savepath ="07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/DT_enrich/rITS_DTcancer_kegg")
data.frame(rITS_DTcancer_kegg)


num_Gene = 20000 
num_rITS = 2802+397+96 
num_DT = 397+96+294+1906
num_DTcancer = 96+294
num_rITS_DT = 397+96
num_rITS_DTcancer = 96

phyper(num_rITS_DT-1, num_DT, num_Gene-num_DT, num_rITS, lower.tail=F)
phyper(num_rITS_DTcancer-1, num_DTcancer, num_Gene-num_DTcancer, num_rITS, lower.tail=F)


# ------------------------------------------------------------------------------
# OG_OncoKB

# 2) OG-Drug -------------------------------------------------------------------
filepath <- ("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/drug_merge3DB/GSEAenrich/ITSp_s_vs_OG_oncoKB")

OG_drug_s <- read.table(paste0(filepath, ".txt"), sep = '\t', header = T)
df1 <- data.frame(table(OG_drug_s$genename))
names(df1) <- c('genename', 'num_drug')

OG_drugFDA_s <- read.table(paste0(filepath, "_FDA.txt"), sep = '\t', header = T)
df2 <- data.frame(table(OG_drugFDA_s$genename))
names(df2) <- c('genename', 'num_FDAdrug')


OG_drugFDAcancer_s <- OG_drugFDA_s[which(OG_drugFDA_s$ifCancer_drugbank == 'yes'), ]
df3 <- data.frame(table(OG_drugFDAcancer_s$genename))
names(df3) <- c('genename', 'num_FDAcancerdrug')
df_s <- merge(df1, df2, by = 'genename', all = T)
df_s <- merge(df_s, df3, by = 'genename', all = T)




# 2) OG-Drug -------------------------------------------------------------------
filepath <- ("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/drug_merge3DB/GSEAenrich/ITSp_r_vs_OG_oncoKB")

OG_drug_r <- read.table(paste0(filepath, ".txt"), sep = '\t', header = T)
df1 <- data.frame(table(OG_drug_r$genename))
names(df1) <- c('genename', 'num_drug')

OG_drugFDA_r <- read.table(paste0(filepath, "_FDA.txt"), sep = '\t', header = T)
df2 <- data.frame(table(OG_drugFDA_r$genename))
names(df2) <- c('genename', 'num_FDAdrug')


OG_drugFDAcancer_r <- OG_drugFDA_r[which(OG_drugFDA_r$ifCancer_drugbank == 'yes'), ]
df3 <- data.frame(table(OG_drugFDAcancer_r$genename))
names(df3) <- c('genename', 'num_FDAcancerdrug')
df_r <- merge(df1, df2, by = 'genename', all = T)
df_r <- merge(df_r, df3, by = 'genename', all = T)

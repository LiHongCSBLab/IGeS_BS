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
library(ggrepel)
library(patchwork)
library(aplot)
library(reshape2)

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



# ----------------------------------------------------------------------------
savepath = "07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/enrich/"

load(paste0(savepath, "custom_GSEA_s_scaled.Rdata"))
res2=data.frame(res)
gsname = "ITSp_vs_OG_oncoKB"

glist <- unlist(strsplit(res2[res2$ID == "ITSp_vs_OG_oncoKB",]$core_enrichment, "/"))


# ------------------------------------------------------------------------------
savepath = paste0("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/function_enrich/", gsname, "/")
filename = "gene_cancer_ITSrate_stat_s"
ITS_dc_df <- read.csv(paste0(savepath, filename, ".csv"), row.names = 1)
genelevel <- rownames(ITS_dc_df)
genelevel_s <- rownames(ITS_dc_df)

# ------------------------------------------------------------------------------
geneKEGGannot <- read.csv("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/function_enrich/ITSp_vs_OG_oncoKB/gene_cancer_ITSrate_stat_s_geneKEGGannot.csv", row.names = 1)



# 2) OG-Drug -------------------------------------------------------------------
filepath <- ("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/drug_merge3DB/GSEAenrich/ITSp_s_vs_OG_oncoKB")

drugstat_overall <- read.table(paste0(filepath, "_stat.txt"), sep = '\t', header = T)

OG_drug <- read.table(paste0(filepath, ".txt"), sep = '\t',header = T)
OG_drugFDA <- read.table(paste0(filepath, "_FDA.txt"), sep = '\t', header = T)
OG_drugFDAcancer <- OG_drugFDA[which(OG_drugFDA$ifCancer_drugbank == "yes"), ]



OG_drug_stat <- lapply(as.list(glist), function(gene){
    # gene=glist[1]
    if(is.element(gene, OG_drug$genename)) {
    dt1 <- OG_drug[which(OG_drug$genename == gene), ]
    dt2 <- OG_drugFDA[which(OG_drugFDA$genename == gene), ]
    dt3 <- dt2[which(dt2$ifCancer_drugbank == "yes"), ]
    df <- data.frame(genename = gene, 
                     NumDrug = length(unique(dt1$drugname_lower)),
                     NumDrugFDA = length(unique(dt2$drugname_lower)),
                     NumDrugFDAcancer = length(unique(dt3$drugname_lower)))
    }else{
    df <- data.frame(genename = gene, 
                     NumDrug = NA,
                     NumDrugFDA = NA,
                     NumDrugFDAcancer = NA)

    }
    return(df)
})
OG_drug_stat <- do.call(rbind, OG_drug_stat)
OG_drug_stat_s <- do.call(rbind, OG_drug_stat)





# ------------------------------------------------------------------------------
# plot v1 
library(corrplot)

plot_savepath <- "07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/function_enrich/ITSp_vs_OG_oncoKB/plot/"
dir.create(plot_savepath)

# 1) overall GSEA enrichment plot
pdf(paste0(plot_savepath, "overall_gseaplot_s.pdf"), width = 8, height = 4)
# gseaplot2(res,geneSetID = "ITSp_vs_OG_oncoKB")
gseaplot(res, geneSetID = "ITSp_vs_OG_oncoKB", by = "runningScore", title = "ITSp_vs_OG_oncoKB")
dev.off()


# 2) core enriched gene dot heatmap
filename = "gene_cancer_ITSrate_stat_s.pdf"
"07_plotting_v2/scripts/functions/ITS_g_cancer_plot_function.R"
colors <- colorRampPalette(c("blue", "white", "#B33333"))(20)
ITS_dc_df[is.na(ITS_dc_df)] = 0
pdf(paste0(plot_savepath, filename, ".pdf"),  width = 10, height = 10)
corrplot(as.matrix(ITS_dc_df), method = "pie", 
        col = colors, na.label = " ", tl.col = 'black')
dev.off()

# 3) annotate KEGG pathway of core enriched genes
genelevel

df_pw <- data.frame(KEGG_name = unique(geneKEGGannot$Description), KEGG_Index = 1:length(unique(geneKEGGannot$Description)))

df_g <- data.frame(gene_name = genelevel, gene_Index = 1:length(genelevel))
df_g_na <- data.frame(gene_name = "666",
                      gene_Index = (nrow(df_g)+1):nrow(df_pw))
df_g <- rbind(df_g, 
              df_g_na)



df_g_pw <- geneKEGGannot[c(3,4)]
names(df_g_pw) <- c("KEGG_name", "gene_name")

df_g_pw <- merge(df_g_pw, df_g, by = "gene_name")
df_g_pw <- merge(df_g_pw, df_pw, by = "KEGG_name")
df_g_pw <- df_g_pw[order(df_g_pw$gene_Index),]
df_g_pw <- df_g_pw[order(df_g_pw$KEGG_Index),]
df_g_pw$xmin <- -1
df_g_pw$xmax <- 1
# df_g_pw$KEGG_name <- factor(df_g_pw$KEGG_name, levels = unique(df_g_pw$KEGG_name))
df_g_pw$gene_Index <- factor(df_g_pw$gene_Index, levels = rev(1:length(genelevel)))

df_g$x <- -1
df_g$gene_Index <- factor(df_g$gene_Index, levels = rev(df_g$gene_Index))

df_pw$x <- 1
df_pw$KEGG_Index <- factor(df_pw$KEGG_Index, levels = unique(df_pw$KEGG_Index))

immunePW <- c("Cytokine-cytokine receptor interaction",
              "Chemokine signaling pathway",
              "Fc epsilon RI signaling pathway",
              "Fc gamma R-mediated phagocytosis",
              "Leukocyte transendothelial migration",
              "Natural killer cell mediated cytotoxicity",
              # "PD-L1 expression and PD-1 checkpoint pathway in cancer",
              "Primary immunodeficiency",
              "T cell receptor signaling pathway",
              "B cell receptor signaling pathway")
              # "Th1 and Th2 cell differentiation",
              # "Th17 cell differentiation")
signalingPW <- c("JAK-STAT signaling pathway",
              #  "NF-kappa B signaling pathway",
               "PI3K-Akt signaling pathway",
               "VEGF signaling pathway",
               "mTOR signaling pathway")

df_g_pw$KEGGtype <- "others"
df_g_pw[which(is.element(df_g_pw$KEGG_name, immunePW)),]$KEGGtype <- "immune"
df_g_pw[which(is.element(df_g_pw$KEGG_name, signalingPW)),]$KEGGtype <- "signaling"
df_g_pw$KEGGtype <- factor(df_g_pw$KEGGtype, levels = c("immune", "signaling", 'others'))
df_g_pw_s = df_g_pw
df_pw$KEGGtype <- "others"
df_pw[which(is.element(df_pw$KEGG_name, immunePW)),]$KEGGtype <- "immune"
df_pw[which(is.element(df_pw$KEGG_name, signalingPW)),]$KEGGtype <- "signaling"


p_geneKEGGannot1 <- ggplot()+
  geom_point(data=df_g, aes(x=x,y=gene_Index))  + 
  geom_text_repel(data = df_g, aes(x=x,y=gene_Index, label = gene_name))+
  geom_point(data=df_pw, aes(x=x,y=KEGG_Index)) +
  geom_text_repel(data = df_pw, aes(x=x,y=KEGG_Index, fill = KEGGtype, label = KEGG_name))+
  theme_void()
p_geneKEGGannot2 <- ggplot() +
  geom_point(data=df_g, aes(x=x,y=gene_Index))  + 
#   geom_text_repel(data = df_g, aes(x=x,y=gene_Index, label = gene_name))+
  geom_point(data=df_pw, aes(x=x,y=KEGG_Index)) +
#   geom_text_repel(data = df_pw, aes(x=x,y=KEGG_Index, fill = KEGGtype, label = KEGG_name))+
  geom_curve(data=df_g_pw,aes(x=xmin,y=gene_Index,
                              xend=xmax,yend=KEGG_Index,
                              color=KEGGtype #,
                              #size=rvalue
                              ),
                              curvature = 0)+
  theme_void() + 
  scale_color_manual(values=c("#b85315","#189164","#bababa")) 


p_geneKEGGannot = p_geneKEGGannot1 %>%
  insert_right(p_geneKEGGannot2, width = 0.2)
p_geneKEGGannot

ggsave(paste0(plot_savepath, "p_geneKEGGannot_s.pdf"), p_geneKEGGannot,  width = 10, height = 10)

# 4) annotate drug through target gene mapping on core enriched genes

OG_drug_stat$genename <- factor(OG_drug_stat$genename, levels = rev(genelevel))
dat <- OG_drug_stat
dat$NumDrugFDA <- dat$NumDrugFDA - dat$NumDrugFDAcancer
dat$NumDrug <- dat$NumDrug - dat$NumDrugFDA

df <- melt(dat[-which(names(dat) == "NumDrug")])

p_OG_drug <- ggplot(data = df, mapping = aes(x = genename, y = value, fill = variable)) + 
  geom_bar(stat = 'identity', position = 'fill') +
  theme_bw()+
  # geom_text(aes(y=value, label = value, size = 4)) +  
  labs(x = 'Sample',y = 'frequnency') +
  theme(axis.title =element_text(size = 16),
  axis.text =element_text(size = 14, color = 'black'),
  axis.text.x = element_text(angle = 45, hjust = 1))+
  coord_flip()
ggsave(paste0(plot_savepath, "p_OG_drug_s_v1.pdf"), p_OG_drug,  width = 5, height = 10)

dat_shown <- melt(OG_drug_stat[-which(names(OG_drug_stat) == "NumDrug")])
# dat_shown <- dat_shown[grep("FDA", dat_shown$variable), ]
p_OG_drug <- ggplot(data = df, mapping = aes(x = genename, y = value, fill = variable)) + 
  geom_bar(stat='identity',position='stack') +
  theme_bw()+
  # geom_text(aes(y=value, label = value, size = 4)) +  
  geom_text_repel(inherit.aes = F, 
                  data = dat_shown, 
                  aes(x = genename, y = value, label = value),
                  max.overlaps = getOption("ggrepel.max.overlaps", default = 30))+
  labs(x = 'Sample',y = 'frequnency') +
  theme(axis.title =element_text(size = 16),
  axis.text =element_text(size = 14, color = 'black'),
  axis.text.x = element_text(angle = 45, hjust = 1))+
  coord_flip()
p_OG_drug
ggsave(paste0(plot_savepath, "p_OG_drug_s_v2.pdf"), p_OG_drug,  width = 5, height = 10)

# 5) drug---target OG---KEGG pathway network






# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ITS_r
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

savepath = "07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/enrich/"

load(paste0(savepath, "custom_GSEA_r_scaled.Rdata"))
res2=data.frame(res)
gsname = "ITSp_vs_OG_oncoKB"

glist <- unlist(strsplit(res2[res2$ID == "ITSp_vs_OG_oncoKB",]$core_enrichment, "/"))


# ------------------------------------------------------------------------------
savepath = paste0("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/function_enrich/", gsname, "/")
filename = "gene_cancer_ITSrate_stat_r"
ITS_dc_df <- read.csv(paste0(savepath, filename, ".csv"), row.names = 1)
genelevel <- rownames(ITS_dc_df)
genelevel_r <- rownames(ITS_dc_df)

# ------------------------------------------------------------------------------
geneKEGGannot <- read.csv("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/function_enrich/ITSp_vs_OG_oncoKB/gene_cancer_ITSrate_stat_r_geneKEGGannot.csv", row.names = 1)



# 2) OG-Drug -------------------------------------------------------------------
filepath <- ("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/drug_merge3DB/GSEAenrich/ITSp_r_vs_OG_oncoKB")

drugstat_overall <- read.table(paste0(filepath, "_stat.txt"), sep = '\t', header = T)

OG_drug <- read.table(paste0(filepath, ".txt"), sep = '\t',header = T)
OG_drugFDA <- read.table(paste0(filepath, "_FDA.txt"), sep = '\t', header = T)
OG_drugFDAcancer <- OG_drugFDA[OG_drugFDA$ifCancer_drugbank == "yes", ]



OG_drug_stat <- lapply(as.list(glist), function(gene){
    # gene=glist[1]
    if(is.element(gene, OG_drug$genename)) {
    dt1 <- OG_drug[which(OG_drug$genename == gene), ]
    dt2 <- OG_drugFDA[which(OG_drugFDA$genename == gene), ]
    dt3 <- dt2[which(dt2$ifCancer_drugbank == "yes"), ]
    df <- data.frame(genename = gene, 
                     NumDrug = length(unique(dt1$drugname_lower)),
                     NumDrugFDA = length(unique(dt2$drugname_lower)),
                     NumDrugFDAcancer = length(unique(dt3$drugname_lower)))
    }else{
    df <- data.frame(genename = gene, 
                     NumDrug = NA,
                     NumDrugFDA = NA,
                     NumDrugFDAcancer = NA)

    }
    return(df)
})

OG_drug_stat <- do.call(rbind, OG_drug_stat)
OG_drug_stat_r <- do.call(rbind, OG_drug_stat)


# ------------------------------------------------------------------------------
# plot v1 
library(corrplot)

plot_savepath <- "07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/function_enrich/ITSp_vs_OG_oncoKB/plot/"
dir.create(plot_savepath)

# 1) overall GSEA enrichment plot
pdf(paste0(plot_savepath, "overall_gseaplot_r.pdf"), width = 8, height = 4)
# gseaplot2(res,geneSetID = "ITSp_vs_OG_oncoKB")
gseaplot(res, geneSetID = "ITSp_vs_OG_oncoKB", by = "runningScore", title = "ITSp_vs_OG_oncoKB")
dev.off()


# 2) core enriched gene dot heatmap
filename = "gene_cancer_ITSrate_stat_r.pdf"
# "07_plotting_v2/scripts/functions/ITS_g_cancer_plot_function.R"
colors <- colorRampPalette(c("blue", "white", "#B33333"))(20)
ITS_dc_df[is.na(ITS_dc_df)] = 0
pdf(paste0(plot_savepath, filename, ".pdf"),  width = 10, height = 10)
corrplot(as.matrix(ITS_dc_df), method = "pie", col = colors, na.label = " ", tl.col = 'black')
dev.off()

# 3) annotate KEGG pathway of core enriched genes
# genelevel

df_pw <- data.frame(KEGG_name = unique(geneKEGGannot$Description), KEGG_Index = 1:length(unique(geneKEGGannot$Description)))

df_g <- data.frame(gene_name = genelevel, gene_Index = 1:length(genelevel))
df_g_na <- data.frame(gene_name = "666",
                      gene_Index = (nrow(df_g)+1):nrow(df_pw))
df_g <- rbind(df_g, 
              df_g_na)



df_g_pw <- geneKEGGannot[c(3,4)]
names(df_g_pw) <- c("KEGG_name", "gene_name")

df_g_pw <- merge(df_g_pw, df_g, by = "gene_name")
df_g_pw <- merge(df_g_pw, df_pw, by = "KEGG_name")
df_g_pw <- df_g_pw[order(df_g_pw$gene_Index),]
df_g_pw <- df_g_pw[order(df_g_pw$KEGG_Index),]
df_g_pw$xmin <- -1
df_g_pw$xmax <- 1
# df_g_pw$KEGG_name <- factor(df_g_pw$KEGG_name, levels = unique(df_g_pw$KEGG_name))
df_g_pw$gene_Index <- factor(df_g_pw$gene_Index, levels = rev(1:length(genelevel)))

df_g$x <- -1
df_g$gene_Index <- factor(df_g$gene_Index, levels = rev(df_g$gene_Index))

df_pw$x <- 1
df_pw$KEGG_Index <- factor(df_pw$KEGG_Index, levels = unique(df_pw$KEGG_Index))

immunePW <- c("Cytokine-cytokine receptor interaction")
signalingPW <- c("Calcium signaling pathway",
                # "EGFR tyrosine kinase inhibitor resistance",
                # "Gap junction",
                "JAK-STAT signaling pathway",
                "MAPK signaling pathway",
                # "PI3K-Akt signaling pathway",
                # "Phospholipase D signaling pathway",
                # "Rap1 signaling pathway",
                # "Ras signaling pathway",
                "TGF-beta signaling pathway")             
df_g_pw$KEGGtype <- "others"
df_g_pw[which(is.element(df_g_pw$KEGG_name, immunePW)),]$KEGGtype <- "immune"
df_g_pw[which(is.element(df_g_pw$KEGG_name, signalingPW)),]$KEGGtype <- "signaling"
# df_g_pw$KEGGtype <- factor(df_g_pw$KEGGtype, levels = c("immune", "signaling", 'others'))
df_g_pw_r = df_g_pw

df_pw$KEGGtype <- "others"
df_pw[which(is.element(df_pw$KEGG_name, immunePW)),]$KEGGtype <- "immune"
df_pw[which(is.element(df_pw$KEGG_name, signalingPW)),]$KEGGtype <- "signaling"


p_geneKEGGannot1 <- ggplot()+
  geom_point(data=df_g, aes(x=x,y=gene_Index))  + 
  geom_text_repel(data = df_g, aes(x=x,y=gene_Index, label = gene_name)) +
  geom_point(data=df_pw, aes(x=x,y=KEGG_Index)) +
  geom_text_repel(data = df_pw, aes(x=x,y=KEGG_Index, fill = KEGGtype, label = KEGG_name))+
  theme_void()
p_geneKEGGannot2 <- ggplot() +
  geom_point(data=df_g, aes(x=x,y=gene_Index))  + 
#   geom_text_repel(data = df_g, aes(x=x,y=gene_Index, label = gene_name))+
  geom_point(data=df_pw, aes(x=x,y=KEGG_Index)) +
#   geom_text_repel(data = df_pw, aes(x=x,y=KEGG_Index, fill = KEGGtype, label = KEGG_name))+
  geom_curve(data=df_g_pw,aes(x=xmin,y=gene_Index,
                              xend=xmax,yend=KEGG_Index,
                              color=KEGGtype #,
                              #size=rvalue
                              ),
                              curvature = 0)+
  theme_void() + 
  scale_color_manual(values=c("#b85315","#189164","#bababa")) 


p_geneKEGGannot = p_geneKEGGannot1 %>%
  insert_right(p_geneKEGGannot2, width = 0.2)
p_geneKEGGannot

ggsave(paste0(plot_savepath, "p_geneKEGGannot_r.pdf"), p_geneKEGGannot,  width = 10, height = 10)

# 4) annotate drug through target gene mapping on core enriched genes

OG_drug_stat$genename <- factor(OG_drug_stat$genename, levels = rev(genelevel))
dat <- OG_drug_stat
dat$NumDrugFDA <- dat$NumDrugFDA - dat$NumDrugFDAcancer
dat$NumDrug <- dat$NumDrug - dat$NumDrugFDA

df <- melt(dat[-which(names(dat) == 'NumDrug')])

p_OG_drug <- ggplot(data = df, mapping = aes(x = genename, y = value, fill = variable)) + 
  geom_bar(stat = 'identity', position = 'fill') +
  theme_bw()+
  # geom_text(aes(y=value, label = value, size = 4)) +  
  labs(x = 'Sample',y = 'frequnency') +
  theme(axis.title =element_text(size = 16),
  axis.text =element_text(size = 14, color = 'black'),
  axis.text.x = element_text(angle = 45, hjust = 1))+
  coord_flip()
ggsave(paste0(plot_savepath, "p_OG_drug_r_v1.pdf"), p_OG_drug,  width = 5, height = 8)

dat_shown <- melt(OG_drug_stat[-which(names(OG_drug_stat) == 'NumDrug')])
# dat_shown <- dat_shown[grep("FDA", dat_shown$variable), ]
p_OG_drug <- ggplot(data = df, mapping = aes(x = genename, y = value, fill = variable)) + 
  geom_bar(stat='identity',position='stack') +
  theme_bw()+
  # geom_text(aes(y=value, label = value, size = 4)) +  
  geom_text_repel(inherit.aes = F, 
                  data = dat_shown, 
                  aes(x = genename, y = value, label = value),
                  max.overlaps = getOption("ggrepel.max.overlaps", default = 30))+
  labs(x = 'Sample',y = 'frequnency') +
  theme(axis.title =element_text(size = 16),
  axis.text =element_text(size = 14, color = 'black'),
  axis.text.x = element_text(angle = 45, hjust = 1))+
  coord_flip()
p_OG_drug
ggsave(paste0(plot_savepath, "p_OG_drug_r_v2.pdf"), p_OG_drug,  width = 5, height = 8)







# # ------------------------------------------------------------------------------
# # ------------------------------------------------------------------------------
# #  3) annotate KEGG pathway of core enriched genes (both S and R) --------------------
# # ------------------------------------------------------------------------------
# # ------------------------------------------------------------------------------

# #
# genelevel

# df_pw <- data.frame(KEGG_name = unique(geneKEGGannot$Description), KEGG_Index = 1:length(unique(geneKEGGannot$Description)))

# df_g <- data.frame(gene_name = genelevel, gene_Index = 1:length(genelevel))
# df_g_na <- data.frame(gene_name = "666",
#                       gene_Index = (nrow(df_g)+1):nrow(df_pw))
# df_g <- rbind(df_g, 
#               df_g_na)



# df_g_pw <- geneKEGGannot[c(3,4)]
# names(df_g_pw) <- c("KEGG_name", "gene_name")

# df_g_pw <- merge(df_g_pw, df_g, by = "gene_name")
# df_g_pw <- merge(df_g_pw, df_pw, by = "KEGG_name")
# df_g_pw <- df_g_pw[order(df_g_pw$gene_Index),]
# df_g_pw <- df_g_pw[order(df_g_pw$KEGG_Index),]
# df_g_pw$xmin <- -1
# df_g_pw$xmax <- 1
# # df_g_pw$KEGG_name <- factor(df_g_pw$KEGG_name, levels = unique(df_g_pw$KEGG_name))
# df_g_pw$gene_Index <- factor(df_g_pw$gene_Index, levels = rev(1:length(genelevel)))

# df_g$x <- -1
# df_g$gene_Index <- factor(df_g$gene_Index, levels = rev(df_g$gene_Index))

# df_pw$x <- 1
# df_pw$KEGG_Index <- factor(df_pw$KEGG_Index, levels = unique(df_pw$KEGG_Index))

# immunePW <- c("Cytokine-cytokine receptor interaction"# ,
#             #  "B cell receptor signaling pathway",
#             #  "Chemokine signaling pathway",
#             #   "Fc epsilon RI signaling pathway",
#             #  "Fc gamma R-mediated phagocytosis",
#             #  "Leukocyte transendothelial migration",
#             #  "Natural killer cell mediated cytotoxicity",
#             #  "PD-L1 expression and PD-1 checkpoint pathway in cancer",
#             #  "Primary immunodeficiency",
#             #  "T cell receptor signaling pathway",
#             #  "Th1 and Th2 cell differentiation",
#             # "Th17 cell differentiation"
#             )
# signalingPW <- c("JAK-STAT signaling pathway",
#                "NF-kappa B signaling pathway",
#                "PI3K-Akt signaling pathway",
#                "MAPK signaling pathway",
#                "Ras signaling pathway",
#                "Rap1 signaling pathway",
#                "Calcium signaling pathway",
#                "Phospholipase D signaling pathway",
#                "TGF-beta signaling pathway")             
# df_g_pw$KEGGtype <- "others"
# df_g_pw[which(is.element(df_g_pw$KEGG_name, immunePW)),]$KEGGtype <- "immune"
# df_g_pw[which(is.element(df_g_pw$KEGG_name, signalingPW)),]$KEGGtype <- "signaling"
# df_g_pw$KEGGtype <- factor(df_g_pw$KEGGtype, levels = c("immune", "signaling", 'others'))

# df_pw$KEGGtype <- "others"
# df_pw[which(is.element(df_pw$KEGG_name, immunePW)),]$KEGGtype <- "immune"
# df_pw[which(is.element(df_pw$KEGG_name, signalingPW)),]$KEGGtype <- "signaling"

# df_pw_s <- df_pw
# df_g_s <- df_g
# df_g_pw_s <- df_g_pw




# # 3) annotate KEGG pathway of core enriched genes
# genelevel

# df_pw <- data.frame(KEGG_name = unique(geneKEGGannot$Description), KEGG_Index = 1:length(unique(geneKEGGannot$Description)))

# df_g <- data.frame(gene_name = genelevel, gene_Index = 1:length(genelevel))
# df_g_na <- data.frame(gene_name = "666",
#                       gene_Index = (nrow(df_g)+1):nrow(df_pw))
# df_g <- rbind(df_g, 
#               df_g_na)



# df_g_pw <- geneKEGGannot[c(3,4)]
# names(df_g_pw) <- c("KEGG_name", "gene_name")

# df_g_pw <- merge(df_g_pw, df_g, by = "gene_name")
# df_g_pw <- merge(df_g_pw, df_pw, by = "KEGG_name")
# df_g_pw <- df_g_pw[order(df_g_pw$gene_Index),]
# df_g_pw <- df_g_pw[order(df_g_pw$KEGG_Index),]
# df_g_pw$xmin <- -1
# df_g_pw$xmax <- 1
# # df_g_pw$KEGG_name <- factor(df_g_pw$KEGG_name, levels = unique(df_g_pw$KEGG_name))
# df_g_pw$gene_Index <- factor(df_g_pw$gene_Index, levels = rev(1:length(genelevel)))

# df_g$x <- -1
# df_g$gene_Index <- factor(df_g$gene_Index, levels = rev(df_g$gene_Index))

# df_pw$x <- 1
# df_pw$KEGG_Index <- factor(df_pw$KEGG_Index, levels = unique(df_pw$KEGG_Index))

# immunePW <- c("Cytokine-cytokine receptor interaction")
# signalingPW <- c("Calcium signaling pathway",
#                 "EGFR tyrosine kinase inhibitor resistance",
#                 # "Gap junction",
#                 "JAK-STAT signaling pathway",
#                 "MAPK signaling pathway",
#                 "PI3K-Akt signaling pathway",
#                 "Phospholipase D signaling pathway",
#                 "Rap1 signaling pathway",
#                 "Ras signaling pathway",
#                 "TGF-beta signaling pathway")             
# df_g_pw$KEGGtype <- "others"
# df_g_pw[which(is.element(df_g_pw$KEGG_name, immunePW)),]$KEGGtype <- "immune"
# df_g_pw[which(is.element(df_g_pw$KEGG_name, signalingPW)),]$KEGGtype <- "signaling"
# df_g_pw$KEGGtype <- factor(df_g_pw$KEGGtype, levels = c("immune", "signaling", 'others'))

# df_pw$KEGGtype <- "others"
# df_pw[which(is.element(df_pw$KEGG_name, immunePW)),]$KEGGtype <- "immune"
# df_pw[which(is.element(df_pw$KEGG_name, signalingPW)),]$KEGGtype <- "signaling"


# df_pw_r <- df_pw
# df_g_r <- df_g
# df_g_pw_r <- df_g_pw

# # ------------------------------------------------------------------------------

# df_g_s$flag  <- 's'
# df_g_s <- df_g_s[- which(df_g_s$gene_name == '666'),]
# df_g_r$flag <- 'r'
# df_g_r <- df_g_r[- which(df_g_r$gene_name == '666'),]


# df_g <- rbind(df_g_s, df_g_r)
# df_g$gene_Index <- nrow(df_g):1
# df_g$gene_name<- factor(df_g$gene_name, levels = (df_g$gene_name))
# # df_g$gene_Index<- factor(df_g$gene_Index, levels = rev(df_g$gene_name))



# df_pw <- unique(rbind(df_pw_s, df_pw_r))
# df_pw<- df_pw[-grep('KEGG_Index', names(df_pw))]
# df_pw <- unique(df_pw)
# df_pw <- df_pw[order(df_pw$KEGGtype), ]
# df_pw$KEGGtype <- factor(df_pw$KEGGtype, levels = c("immune", "signaling", 'others'))
# df_pw$KEGG_name <- factor(df_pw$KEGG_name, levels = df_pw$KEGG_name)
# df_pw$KEGG_Index <- nrow(df_pw):1

# df_g_pw <- rbind(df_g_pw_s, df_g_pw_r)
# df_g_pw <- unique(df_g_pw[-c(3,4,8)])
# df_g_pw <- merge(df_g_pw, df_pw[c(1,4)], by = "KEGG_name")
# df_g_pw <- merge(df_g_pw, df_g[c(1,2)], by = "gene_name")


# df_g_pw$gene_name <- factor(df_g_pw$gene_name, levels = df_g$gene_name)

# df_g_pw$KEGGtype <- factor(df_g_pw$KEGGtype, levels = c("immune", "signaling", 'others'))


# p_geneKEGGannot1 <- ggplot()+
#   geom_point(data=df_g, aes(x=x,y=gene_Index))  + 
#   geom_text(data = df_g, aes(x=x,y=gene_Index, label = gene_name))+
# #   geom_text_repel(data = df_g, aes(x=x,y=gene_Index, label = gene_name))+
#   geom_point(data=df_pw, aes(x=x,y=KEGG_Index)) +
#   geom_text(data = df_pw, aes(x=x,y=KEGG_Index, fill = KEGGtype, label = KEGG_name))+
# #   geom_text_repel(data = df_pw, aes(x=x,y=KEGG_Index, fill = KEGGtype, label = KEGG_name))+
#   theme_void()
# p_geneKEGGannot2 <- ggplot() +
#   geom_point(data=df_g, aes(x=x,y=gene_Index))  + 
# #   geom_text_repel(data = df_g, aes(x=x,y=gene_Index, label = gene_name))+
#   geom_point(data=df_pw, aes(x=x,y=KEGG_Index)) +
# #   geom_text_repel(data = df_pw, aes(x=x,y=KEGG_Index, fill = KEGGtype, label = KEGG_name))+
#   geom_curve(data=df_g_pw,aes(x=xmin,y=gene_Index,
#                               xend=xmax,yend=KEGG_Index,
#                               color=KEGGtype #,
#                               #size=rvalue
#                               ),
#                               curvature = 0)+
#   theme_void() + 
#   scale_color_manual(values=c("#b85315","#189164","#bababa")) 


# # p_geneKEGGannot = p_geneKEGGannot1 %>%
# #   insert_right(p_geneKEGGannot2, width = 0.2)
# p_geneKEGGannot = p_geneKEGGannot1|p_geneKEGGannot2

# plot_savepath <- "07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/function_enrich/ITSp_vs_OG_oncoKB/plot/"

# ggsave(paste0(plot_savepath, "p_geneKEGGannot.pdf"), p_geneKEGGannot,  width = 15, height = 20)




# # ------------------------------------------------------------------------------
# # ------------------------------------------------------------------------------
# #  4) annotate drug through target gene mapping on core enriched genes (both S and R) --------------------
# # ------------------------------------------------------------------------------
# # ------------------------------------------------------------------------------

# # 4) annotate drug through target gene mapping on core enriched genes
# OG_drug_stat_s$genename <- factor(OG_drug_stat_s$genename, levels = rev(genelevel))
# dat <- OG_drug_stat_s
# dat$NumDrugFDA <- dat$NumDrugFDA - dat$NumDrugFDAcancer
# dat$NumDrug <- dat$NumDrug - dat$NumDrugFDA

# df_s <- melt(dat)


# OG_drug_stat_r$genename <- factor(OG_drug_stat_r$genename, levels = rev(genelevel))
# dat <- OG_drug_stat
# dat$NumDrugFDA <- dat$NumDrugFDA - dat$NumDrugFDAcancer
# dat$NumDrug <- dat$NumDrug - dat$NumDrugFDA

# df_r <- melt(dat)

# df <- unique(rbind(df_s, df_r))
# df$genename <- factor(df$genename, levels = rev(c(genelevel_s, genelevel_r)))

# p_OG_drug <- ggplot(data = df, mapping = aes(x = genename, y = value, fill = variable)) + 
#   geom_bar(stat = 'identity', position = 'fill') +
#   theme_bw()+
#   # geom_text(aes(y=value, label = value, size = 4)) +  
#   labs(x = 'Sample',y = 'frequnency') +
#   theme(axis.title =element_text(size = 16),
#   axis.text =element_text(size = 14, color = 'black'),
#   axis.text.x = element_text(angle = 45, hjust = 1))+
#   coord_flip()
# ggsave(paste0(plot_savepath, "p_OG_drug_v1.pdf"), p_OG_drug,  width = 5, height = 16)

# OG_drug_stat <- rbind(OG_drug_stat_s, OG_drug_stat_r)
# dat_shown <- melt(OG_drug_stat)
# # dat_shown <- dat_shown[grep("FDA", dat_shown$variable), ]
# p_OG_drug <- ggplot(data = df, mapping = aes(x = genename, y = value, fill = variable)) + 
#   geom_bar(stat='identity',position='stack') +
#   theme_bw()+
# #   geom_text(aes(y=value, label = value, size = 4)) +  
#   geom_text_repel(inherit.aes = F, 
#                   data = dat_shown, 
#                   aes(x = genename, y = value, label = value),
#                   max.overlaps = getOption("ggrepel.max.overlaps", default = 30))+
#   labs(x = 'Sample',y = 'frequnency') +
#   theme(axis.title =element_text(size = 16),
#   axis.text =element_text(size = 14, color = 'black'),
#   axis.text.x = element_text(angle = 45, hjust = 1))+
#   coord_flip()
# p_OG_drug
# ggsave(paste0(plot_savepath, "p_OG_drug_v2.pdf"), p_OG_drug,  width = 5, height = 16)

# -----------------------------------------------------------------------------=
# -----------------------------------------------------------------------------=
# 5) drug---target OG---KEGG pathway network
# -----------------------------------------------------------------------------=
# -----------------------------------------------------------------------------=
library(igraph)
library(devtools)
library(visNetwork)
library(dplyr)

# OG_drug <- read.table(paste0(filepath, ".txt"), sep = '\t',header = T)
# OG_drugFDA <- read.table(paste0(filepath, "_FDA.txt"), sep = '\t', header = T)
# OG_drugFDAcancer <- OG_drugFDA[OG_drugFDA$ifCancer_drugbank == "yes", ]
df_g_pw_s$flag = 's'
df_g_pw_r$flag = 'r'
df_g_pw = rbind(df_g_pw_s, df_g_pw_r)

names(OG_drug) <- c("drugname", "gene_name","MoA",  "database","drugname_lower")
df_OG_drug <- merge(unique(df_g_pw[c("KEGG_name","gene_name","KEGGtype")]),
                OG_drug, by = "gene_name")
write.table(df_OG_drug, paste0(plot_savepath, "OG_drug_net.txt"), sep = '\t', row.names = F, quote = F)

names(OG_drugFDA) <- c("drugname", "gene_name","MoA",  "database","drugname_lower","ifCancer_drugbank")
df_net <- merge(unique(df_g_pw[c("KEGG_name","gene_name","KEGGtype")]),
                OG_drugFDA, by = "gene_name")
write.table(df_net, paste0(plot_savepath, "OG_drugFDA_net.txt"), sep = '\t', row.names = F, quote = F)
df_net2 <- unique(df_net[c("drugname_lower","gene_name","KEGG_name","KEGGtype",
                     "ifCancer_drugbank")])
df_net2 <- df_net2[df_net2$KEGGtype != 'others',]


df1 <- unique(df_net2[c("drugname_lower", "gene_name")])
df2 <- unique(df_net2[c("KEGG_name", "gene_name")])
names(df1) <- c("N1", "N2")
names(df2) <- c("N1", "N2")
df <- rbind(df1, df2)

nodes = union(df$N1,df$N2)
net = igraph::graph_from_data_frame(d=df, vertices = unique(nodes), directed = F)


edge.color <- colorRampPalette(c("#D6D6D6", "#383838"), alpha=TRUE)
igraph::E(net)$color <- edge.color(igraph::ecount(net))


igraph::V(net)$color = "#E7F3FD"
igraph::V(net)[which(is.element(igraph::V(net)$name, df_net2$gene_name))]$color = "#a0c7dc"
igraph::V(net)[which(is.element(igraph::V(net)$name, df_net2[df_net2$ifCancer_drugbank == 'yes', ]$drugname_lower))]$color = "#f4d66cd0"
igraph::V(net)[which(is.element(igraph::V(net)$name, df_net2[is.na(df_net2$ifCancer_drugbank), ]$drugname_lower))]$color = "#b3ffd8"
igraph::V(net)[which(is.element(igraph::V(net)$name, df_net2[df_net2$KEGGtype == 'immune', ]$KEGG_name))]$color = "#be8866"
igraph::V(net)[which(is.element(igraph::V(net)$name, df_net2[df_net2$KEGGtype == 'signaling', ]$KEGG_name))]$color = "#f8e3d7"


igraph::V(net)$shape = "box"
igraph::V(net)[which(igraph::V(net)$name %in% df_net2$gene_name)]$shape = "round"
igraph::V(net)[which(igraph::V(net)$name %in% df_net2$KEGG_name)]$shape = "ellipse"

igraph::V(net)[which(is.element(igraph::V(net)$name, df_net2[df_net2$ifCancer_drugbank == 'yes', ]$drugname_lower))]$group = "cancer"
igraph::V(net)[which(is.element(igraph::V(net)$name, df_net2[is.na(df_net2$KEGGtype), ]$drugname_lower))]$group = "noncancer"

igraph::V(net)$font.size = 40

data <- toVisNetworkData(net)

visNetwork(nodes = data$nodes, edges = data$edges, width = "90%", height = "95vh")%>%
  visNodes(size = 40)%>%
  # visHierarchicalLayout(direction = "LR", levelSeparation = 500)%>% 
  visIgraphLayout(layout = "layout_on_sphere",type = "full",randomSeed = 123) %>%
  visOptions(highlightNearest = list(enabled = TRUE,  hideColor = "lightgrey", hover = T),
             nodesIdSelection =list(enabled = TRUE), selectedBy = "group") %>%
  # visConfigure(enabled = TRUE) %>%
  addFontAwesome() %>%
  visGroups(groupname = "noncancer", color = "#ff1717")%>%
  visGroups(groupname = "cancer", color = "#0631b0")%>%
  visLegend() %>%
  visInteraction(navigationButtons = TRUE) %>%
  visOptions(manipulation = TRUE) %>% 
  visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
  visSave(file = paste0(plot_savepath, "networkAnalysis.html"))





df <- df_net2[c("drugname_lower", "KEGG_name")]
names(df) <- c("N1", "N2")

nodes = union(df$N1,df$N2)
net = igraph::graph_from_data_frame(d=df, vertices = unique(nodes), directed = F)


edge.color <- colorRampPalette(c("#D6D6D6", "#383838"), alpha=TRUE)
igraph::E(net)$color <- edge.color(igraph::ecount(net))


igraph::V(net)$color = "#E7F3FD"
igraph::V(net)[which(is.element(igraph::V(net)$name, df_net2$gene_name))]$color = "#a0c7dc"
igraph::V(net)[which(is.element(igraph::V(net)$name, df_net2[df_net2$ifCancer_drugbank == 'yes', ]$drugname_lower))]$color = "#f4d66cd0"
igraph::V(net)[which(is.element(igraph::V(net)$name, df_net2[is.na(df_net2$ifCancer_drugbank), ]$drugname_lower))]$color = "#b3ffd8"
igraph::V(net)[which(is.element(igraph::V(net)$name, df_net2[df_net2$KEGGtype == 'immune', ]$KEGG_name))]$color = "#be8866"
igraph::V(net)[which(is.element(igraph::V(net)$name, df_net2[df_net2$KEGGtype == 'signaling', ]$KEGG_name))]$color = "#f8e3d7"


igraph::V(net)$shape = "box"
igraph::V(net)[which(igraph::V(net)$name %in% df_net2$gene_name)]$shape = "round"
igraph::V(net)[which(igraph::V(net)$name %in% df_net2$KEGG_name)]$shape = "ellipse"

igraph::V(net)[which(is.element(igraph::V(net)$name, df_net2[df_net2$ifCancer_drugbank == 'yes', ]$drugname_lower))]$group = "cancer"
igraph::V(net)[which(is.element(igraph::V(net)$name, df_net2[is.na(df_net2$KEGGtype), ]$drugname_lower))]$group = "noncancer"

igraph::V(net)$font.size = 40

data <- toVisNetworkData(net)

visNetwork(nodes = data$nodes, edges = data$edges, width = "90%", height = "95vh")%>%
  visNodes(size = 40)%>%
  # visHierarchicalLayout(direction = "LR", levelSeparation = 500)%>% 
  visIgraphLayout(layout = "layout_on_sphere",type = "full",randomSeed = 123) %>%
  visOptions(highlightNearest = list(enabled = TRUE,  hideColor = "lightgrey", hover = T),
             nodesIdSelection =list(enabled = TRUE), selectedBy = "group") %>%
  # visConfigure(enabled = TRUE) %>%
  addFontAwesome() %>%
  visGroups(groupname = "noncancer", color = "#ff1717")%>%
  visGroups(groupname = "cancer", color = "#0631b0")%>%
  visLegend() %>%
  visInteraction(navigationButtons = TRUE) %>%
  visOptions(manipulation = TRUE) %>% 
  visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
  visSave(file = paste0(plot_savepath, "networkAnalysis2.html"))

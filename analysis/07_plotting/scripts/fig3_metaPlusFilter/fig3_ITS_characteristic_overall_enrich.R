# rm(list=ls())
work_path <- "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/"

setwd(work_path)
options(stringsAsFactors = F)

library(reshape2)
library(ggplot2)
library(ggplotify)
library(grid)
library(cowplot)
library(parallel)
library(patchwork)
library(dplyr)
library(stringr)
library(data.table)
require(GGally)
library(plot3D)
library(fmsb)
library(readxl)
library(dplyr)

source("06_results_summaryandplotting/scripts/plot_for_drug_function.R")
source("06_results_summaryandplotting/scripts/ITS_analysis/ITS_gene_stat_function.R")
source("06_results_summaryandplotting/scripts/immunesig_selection_function.R")
source("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/scripts/drugITGSscore_lincs.R")


# ----------------------------------------------------------------------------

dir.create("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/")
# ----------------------------------------------------------------------------
# load meta analysis result

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

# LOAD GENES IN CCLE -------------------------------------------------------
geneExprInCCLE <- read.csv("03_drug_immuneSig_enrichment/data/CCLE_geneexpr/geneExprInCCLE.csv")

protein_coding_genes <- read.table("03_drug_immuneSig_enrichment/data/protein_coding_genes_hg38.txt", sep = '\t', header = T)
geneExprInCCLE_filtered <- inner_join(geneExprInCCLE, protein_coding_genes)
geneExprInCCLE_filtered <- geneExprInCCLE_filtered[geneExprInCCLE_filtered$rate < 0.2,]

# LOAD GENES IN LINCS ------------------------------------------------------
lincs_gene = read.csv("/picb/bigdata/project/FengFYM/drug_MoA_reanalysis/LINCS_data/GSE70138/file/GSE92742_Broad_LINCS_gene_info.txt", sep = "\t")


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
load( "07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITS_set.Rdata")



# ------------------------------------------------------------------------------
# 1. ITS varies in cancer types ---------------------------------------------------
# ------------------------------------------------------------------------------

library(UpSetR)
library(scales)
library(gridExtra)
blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )

#fromList
# listinput <- list(one = c(1, 2, 3, 5, 7, 8, 11, 12, 13), 
#                   two = c(1, 2, 4, 5, 10), 
#                   three = c(1, 5, 6, 7, 8, 9, 10, 12, 13))
dir.create("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITS_overall/")
pdf(paste0("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITS_overall/ITS_p_cancer_overlap.pdf"),
    width = 20,
    height = 10, 
    onefile = TRUE)

ITS_p_set = lapply(ITS_p_set,function(x){x[!is.na(names(x))]})
ITS_p_genestat = list()
dfstat = list()

for(i in 1:length(ITS_p_set[[1]])){
  # i = 2
  # i = 21
  ITSname <- names(ITS_p_set[[1]])[i]
  listinput <- lapply(ITS_p_set, function(x){x[[ITSname]]})
  genelist = unique(unlist(listinput))
  listinput_df <- fromList(listinput)
  rownames(listinput_df) <- genelist
  table(apply(listinput_df,1,sum))
  apply(listinput_df,2,sum)


  p = upset(fromList(listinput), 
            order.by = "freq",
            nsets = 22,
            nintersects = 30,
            scale.intersections = "identity")
  # View(upset)
  # dfpie = data.frame(group = c('gene_union', '50percent_cancer_shared'),
  #                    value = c(length(unique(unlist(listinput))),
  #                              length(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)>0.5*length(listinput))])))
  
  dfpie = data.frame(group = c('less50percent_cancer_shared', '50percent_cancer_shared'),
                     value = c(length(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)<=0.5*length(listinput))]),
                               length(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)>0.5*length(listinput))])))
  
  pie <- ggplot(dfpie, aes(x="", y=value, fill=group))+
    geom_bar( width = 1, stat = "identity", alpha = 0.5) + 
    coord_polar("y", start=0) + 
    guides(fill=FALSE) +
    scale_fill_brewer(palette="Dark2") +
    blank_theme +
    theme(axis.text.x=element_blank())+
    geom_text(aes(y = value/3 + c(0, cumsum(value)[-length(value)]), 
                  label = paste0(group, "\n", percent(value/sum(value)))), size=5)+
    ggtitle(names(ITS_p_set[[1]])[i])
  # pie
  
  library("ggplotify")
  # library(ggimage)
  p1 <- as.grob(as.ggplot(p))
  p2 <- as.grob(pie)
  grid.arrange(p2, p1, nrow = 1, widths = c(0.5, 2))

  
  ITS_p_genestat[[i]] <- listinput_df
  dfstat[[i]] <- data.frame(immune_sig = ITSname,
                            group = c('gene_cancerspecific', 'gene_2cancer','gene_3cancer','50percent_cancer_shared'),
                            num_ITSgene = nrow(listinput_df),
                            num_gene = c(length(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)==1)]),
                                      length(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)==2)]),
                                      length(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)==3)]),
                                      length(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)>=0.5*22)])),
                            value = c(length(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)==1)])/ nrow(listinput_df),
                                      length(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)==2)])/ nrow(listinput_df),
                                      length(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)==3)])/ nrow(listinput_df),
                                      length(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)>=0.5*22)])/nrow(listinput_df)))
  

}

dev.off()


names(ITS_p_genestat) <- names(ITS_p_set[[1]])
# names(dfstat) <- names(ITS_p_set[[1]])
save(ITS_p_genestat, dfstat, file = "07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITS_overall/ITS_p_genestat.Rdata")

# ------------------------------------------------------------------------------
load( "07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITS_overall/ITS_p_genestat.Rdata")
dfstat_all <- do.call(rbind, dfstat)
write.csv(dfstat_all, "07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITS_overall/ITS_p_genestat.csv", quote = F, row.names = F)


dfstat_all1 <- dfstat_all[dfstat_all$group == 'gene_cancerspecific',]
dfstat_all1 <- dfstat_all1[order(dfstat_all1$value, decreasing = T),]
dfstat_all$immune_sig <- factor(dfstat_all$immune_sig, levels = dfstat_all1$immune_sig)

dfstat_p <- ggplot(data = dfstat_all, aes(x = immune_sig, y = value, group=group)) + 
  geom_line(aes(color=group))+ 
  geom_point(aes(color=group)) +
  theme_bw() + 
  theme(axis.text.x=element_text(angle=90,size=8, vjust = 1, hjust = 0.5))
dfstat_p
ggsave("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITS_overall/ITS_p_overallstat_4lines.pdf",dfstat_p, width = 13, height = 10)

dfstat_all2 <- dfstat_all[dfstat_all$group != "gene_3cancer", ]
dfstat_p <- ggplot(data = dfstat_all2, aes(x = immune_sig, y = value, group=group)) + 
  geom_line(aes(color=group))+ 
  geom_point(aes(color=group)) +
  theme_bw() + 
  theme(axis.text.x=element_text(angle=90,size=8, vjust = 1, hjust = 0.5))
dfstat_p
ggsave("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITS_overall/ITS_p_overallstat_3lines.pdf", dfstat_p, width = 13, height = 10)


dfstat = list()
for(i in 1:length(ITS_p_set[[1]])){
  # i = 1
  # i = 21
  ITSname <- names(ITS_p_set[[1]])[i]
  listinput <- lapply(ITS_p_set, function(x){x[[ITSname]]})
  genelist = unique(unlist(listinput))
  listinput_df <- fromList(listinput)
  rownames(listinput_df) <- genelist
  (apply(listinput_df,1,sum))
  apply(listinput_df,2,sum)

  ITS_p_genestat[[i]] <- listinput_df
  dfstat[[i]] <- data.frame(immune_sig = ITSname,
                            group = c('gene_cancerspecific', 'gene_2cancer','gene_3cancer','50percent_cancer_shared'),
                            num_gene = c(length(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)==1)]),
                                      length(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)>=2)]),
                                      length(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)>=3)]),
                                      length(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)>= 0.5* length(listinput))])),
                            value = c(length(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)==1)])/ nrow(listinput_df),
                                      length(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)>=2)])/ nrow(listinput_df),
                                      length(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)>=3)])/ nrow(listinput_df),
                                      length(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)>=0.5*length(listinput))])/ nrow(listinput_df)))
  
}

save(ITS_p_genestat, dfstat, file = "07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITS_overall/ITS_p_genestat_binary.Rdata")

dfstat_all <- do.call(rbind, dfstat)
write.csv(dfstat_all, "07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITS_overall/ITS_p_genestat_binary.csv", quote = F, row.names = F)


dfstat_all1 <- dfstat_all[dfstat_all$group == 'gene_cancerspecific',]
dfstat_all1 <- dfstat_all1[order(dfstat_all1$value, decreasing = T),]
dfstat_all$immune_sig <- factor(dfstat_all$immune_sig, levels = dfstat_all1$immune_sig)
dfstat_all2 <- dfstat_all[is.element(dfstat_all$group, c('gene_cancerspecific', 'gene_2cancer')), ]
dfstat_p <- ggplot(data = dfstat_all2, aes(x = immune_sig, y = value, group=group)) + 
  geom_line(aes(color=group))+ 
  geom_point(aes(color=group)) +
  theme_bw() + 
  theme(axis.text.x=element_text(angle=90,size=8, vjust = 1, hjust = 0.5))
dfstat_p

ggsave("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITS_overall/ITS_p_overallstat_2lines.pdf",dfstat_p, width = 10, height = 5)




# ITS_n ------------------------------------------------------------------------
pdf(paste0("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITS_overall/ITS_n_cancer_overlap.pdf"),
    width = 20,
    height = 10, 
    onefile = TRUE)

ITS_n_set = lapply(ITS_n_set,function(x){x[!is.na(names(x))]})
ITS_n_genestat = list()
dfstat_n = list()

for(i in 1:length(ITS_n_set[[1]])){
  # i = 37
  
  ITSname <- names(ITS_n_set[[1]])[i]
  listinput <- lapply(ITS_n_set, function(x){x[[ITSname]]})
  genelist = unique(unlist(listinput))
  listinput_df <- fromList(listinput)
  rownames(listinput_df) <- genelist
  table(apply(listinput_df,1,sum))
  apply(listinput_df,2,sum)
  if(length(apply(listinput_df,1,sum)) > 0){
    p = upset(fromList(listinput), 
              order.by = "freq",
              nsets = 22,
              nintersects = 30,
              scale.intersections = "identity")
    # View(upset)
    dfpie = data.frame(group = c('gene_union', '50percent_cancer_shared'),
                       value = c(length(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)<=0.5*length(listinput))]),
                                 length(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)>0.5*length(listinput))])))
    
    
    pie <- ggplot(dfpie, aes(x="", y=value, fill=group))+
      geom_bar( width = 1, stat = "identity", alpha = 0.5) + 
      coord_polar("y", start=0) + 
      guides(fill=FALSE) +
      scale_fill_brewer(palette="Dark2") +
      blank_theme +
      theme(axis.text.x=element_blank())+
      geom_text(aes(y = value/3 + c(0, cumsum(value)[-length(value)]), 
                    label = percent(value/sum(value))), size=5)+
      ggtitle(names(ITS_n_set[[1]])[i])
    # pie
    
    library("ggplotify")
    # library(ggimage)
    p1 <- as.grob(as.ggplot(p))
    p2 <- as.grob(pie)
    grid.arrange(p2, p1, nrow = 1, widths = c(0.5, 2))

  ITS_n_genestat[[i]] <- listinput_df
  dfstat_n[[i]] <-  data.frame(immune_sig = ITSname,
                            group = c('gene_cancerspecific', 'gene_2cancer','gene_3cancer','50percent_cancer_shared'),
                            num_gene = c(length(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)==1)]),
                                      length(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)>=2)]),
                                      length(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)>=3)]),
                                      length(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)>= 0.5* length(listinput))])),
                            value = c(length(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)==1)])/ nrow(listinput_df),
                                      length(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)>=2)])/ nrow(listinput_df),
                                      length(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)>=3)])/ nrow(listinput_df),
                                      length(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)>=0.5*length(listinput))])/ nrow(listinput_df)))
  
  }
}


dev.off()

names(ITS_n_genestat) <- names(ITS_n_set[[1]][-41])
save(ITS_n_genestat, dfstat_n, file = "07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITS_overall/ITS_n_genestat.Rdata")
# names(dfstat_n) <- names(ITS_n_set[[1]])
dfstat_n_all <- do.call(rbind, dfstat_n)
write.csv(dfstat_n_all, "07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITS_overall/ITS_n_genestat.csv", quote = F, row.names = F)

dfstat_n_all1 <- dfstat_n_all[dfstat_n_all$group == 'gene_cancerspecific',]
dfstat_n_all1 <- dfstat_n_all1[order(dfstat_n_all1$value, decreasing = T),]
dfstat_n_all$immune_sig <- factor(dfstat_n_all$immune_sig, levels = dfstat_n_all1$immune_sig)

dfstat_n_plot <- ggplot(data = dfstat_n_all, aes(x = immune_sig, y = value, group=group)) + 
  geom_line(aes(color=group))+ 
  geom_point(aes(color=group)) +
  theme_bw() + 
  theme(axis.text.x=element_text(angle=90,size=8, vjust = 0.5, hjust = 0.5))
dfstat_n_plot
ggsave("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITS_overall/ITS_n_overallstat_4lines.pdf", dfstat_n_plot, width = 13, height = 10)

dfstat_n_all2 <- dfstat_n_all[dfstat_n_all$group != "gene_3cancer", ]
dfstat_n_plot <- ggplot(data = dfstat_n_all2, aes(x = immune_sig, y = value, group=group)) + 
  geom_line(aes(color=group))+ 
  geom_point(aes(color=group)) +
  theme_bw() + 
  theme(axis.text.x=element_text(angle=90,size=8, vjust = 1, hjust = 0.5))
dfstat_n_plot
ggsave("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITS_overall/ITS_n_overallstat_3lines.pdf", dfstat_n_plot, width = 13, height = 10)




dfstat = list()
for(i in 1:length(ITS_n_set[[1]])){
  # i = 1
  # i = 21
  ITSname <- names(ITS_n_set[[1]])[i]
  listinput <- lapply(ITS_n_set, function(x){x[[ITSname]]})
  genelist = unique(unlist(listinput))
  listinput_df <- fromList(listinput)
  rownames(listinput_df) <- genelist
  table(apply(listinput_df,1,sum))
  apply(listinput_df,2,sum)

  ITS_p_genestat[[i]] <- listinput_df
  dfstat[[i]] <- data.frame(immune_sig = ITSname,
                            group = c('gene_cancerspecific', 'gene_2cancer','gene_3cancer','50percent_cancer_shared'),
                            num_gene = c(length(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)==1)]),
                                      length(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)>=2)]),
                                      length(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)>=3)]),
                                      length(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)>= 0.5* length(listinput))])),
                            value = c(length(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)==1)])/ nrow(listinput_df),
                                      length(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)>=2)])/ nrow(listinput_df),
                                      length(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)>=3)])/ nrow(listinput_df),
                                      length(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)>=0.5*length(listinput))])/ nrow(listinput_df)))
  
}

dfstat_all <- do.call(rbind, dfstat)
write.csv(dfstat_all, "07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITS_overall/ITS_n_genestat_binary.csv", quote = F, row.names = F)


dfstat_all1 <- dfstat_all[dfstat_all$group == 'gene_cancerspecific',]
dfstat_all1 <- dfstat_all1[order(dfstat_all1$value, decreasing = T),]
dfstat_all$immune_sig <- factor(dfstat_all$immune_sig, levels = dfstat_all1$immune_sig)
dfstat_all2 <- dfstat_all[is.element(dfstat_all$group, c('gene_cancerspecific', 'gene_2cancer')), ]
dfstat_n <- ggplot(data = dfstat_all2, aes(x = immune_sig, y = value, group=group)) + 
  geom_line(aes(color=group))+ 
  geom_point(aes(color=group)) +
  theme_bw() + 
  theme(axis.text.x=element_text(angle=90,size=8, vjust = 1, hjust = 0.5))
dfstat_n

ggsave("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITS_overall/ITS_n_overallstat_2lines.pdf",dfstat_n, width = 13, height = 10)



# ------------------------------------------------------------------------------
# stat: sensitive-related and resistant-related ITS ----------------------------
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


# sensitive 

ITS_p_s_genestat = list()
dfstat_p_s = list()
for(i in 1:length(ITS_S)){
  # i = 1
  # i = 21
  ITSname <- ITS_S[i]
  listinput <- lapply(ITS_p_set, function(x){x[[ITSname]]})
  genelist = unique(unlist(listinput))
  listinput_df <- fromList(listinput)
  rownames(listinput_df) <- genelist
  table(apply(listinput_df,1,sum))
  apply(listinput_df,2,sum)

  ITS_p_s_genestat[[i]] <- listinput_df
  dfstat_p_s[[i]] <- data.frame(immune_sig = ITSname,
                            group = c('gene_cancerspecific', 'gene_2cancer','gene_3cancer','50percent_cancer_shared'),
                            num_ITSgene = nrow(listinput_df),
                            num_gene = c(length(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)==1)]),
                                      length(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)==2)]),
                                      length(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)==3)]),
                                      length(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)>=0.5*22)])),
                            value = c(length(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)==1)])/ nrow(listinput_df),
                                      length(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)==2)])/ nrow(listinput_df),
                                      length(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)==3)])/ nrow(listinput_df),
                                      length(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)>=0.5*22)])/nrow(listinput_df)))
  

}



names(ITS_p_s_genestat) <- ITS_p_s_genestat # names(ITS_p_set[[1]])
# names(dfstat) <- names(ITS_p_set[[1]])
save(ITS_p_s_genestat, dfstat, file = "07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITS_overall/ITS_p_s_genestat.Rdata")

dfstat_p_s_all <- do.call(rbind, dfstat_p_s)
write.csv(dfstat_p_s_all, "07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITS_overall/ITS_p_s_genestat.csv", quote = F, row.names = F)


dfstat_p_s_all1 <- dfstat_p_s_all[dfstat_p_s_all$group == 'gene_cancerspecific',]
dfstat_p_s_all1 <- dfstat_p_s_all1[order(dfstat_p_s_all1$value, decreasing = T),]
dfstat_p_s_all$immune_sig <- factor(dfstat_p_s_all$immune_sig, levels = dfstat_p_s_all1$immune_sig)

dfstat_p <- ggplot(data = dfstat_p_s_all, aes(x = immune_sig, y = value, group=group)) + 
  geom_line(aes(color=group))+ 
  geom_point(aes(color=group)) +
  theme_bw() + 
  theme(axis.text.x=element_text(angle=90,size=8, vjust = 1, hjust = 0.5))
dfstat_p
ggsave("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITS_overall/ITS_p_s_overallstat_4lines.pdf",dfstat_p, width = 13, height = 10)

dfstat_p_s_all2 <- dfstat_p_s_all[dfstat_p_s_all$group != "gene_3cancer", ]
dfstat_p <- ggplot(data = dfstat_p_s_all2, aes(x = immune_sig, y = value, group=group)) + 
  geom_line(aes(color=group))+ 
  geom_point(aes(color=group)) +
  theme_bw() + 
  theme(axis.text.x=element_text(angle=90,size=8, vjust = 1, hjust = 0.5))
dfstat_p
ggsave("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITS_overall/ITS_p_s_overallstat_3lines.pdf", dfstat_p, width = 13, height = 10)


dfstat_p_s = list()
for(i in 1:length(ITS_S)){
  # i = 1
  # i = 21
  ITSname <- ITS_S[i]
  listinput <- lapply(ITS_p_set, function(x){x[[ITSname]]})
  genelist = unique(unlist(listinput))
  listinput_df <- fromList(listinput)
  rownames(listinput_df) <- genelist
  table(apply(listinput_df,1,sum))
  apply(listinput_df,2,sum)

  ITS_p_genestat[[i]] <- listinput_df
  dfstat_p_s[[i]] <- data.frame(immune_sig = ITSname,
                            group = c('gene_cancerspecific', 'gene_2cancer','gene_3cancer','50percent_cancer_shared'),
                            num_gene = c(length(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)==1)]),
                                      length(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)>=2)]),
                                      length(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)>=3)]),
                                      length(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)>= 0.5* length(listinput))])),
                            value = c(length(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)==1)])/ nrow(listinput_df),
                                      length(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)>=2)])/ nrow(listinput_df),
                                      length(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)>=3)])/ nrow(listinput_df),
                                      length(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)>=0.5*length(listinput))])/ nrow(listinput_df)))
  
}

dfstat_p_s_all <- do.call(rbind, dfstat_p_s)
write.csv(dfstat_p_s_all, "07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITS_overall/ITS_p_s_genestat_binary.csv", quote = F, row.names = F)


dfstat_p_s_all1 <- dfstat_p_s_all[dfstat_p_s_all$group == 'gene_cancerspecific',]
dfstat_p_s_all1 <- dfstat_p_s_all1[order(dfstat_p_s_all1$value, decreasing = T),]
dfstat_p_s_all$immune_sig <- factor(dfstat_p_s_all$immune_sig, levels = dfstat_p_s_all1$immune_sig)
dfstat_p_s_all2 <- dfstat_p_s_all[is.element(dfstat_p_s_all$group, c('gene_cancerspecific', 'gene_2cancer')), ]
dfstat_p <- ggplot(data = dfstat_p_s_all2, aes(x = immune_sig, y = value, group=group)) + 
  geom_line(aes(color=group))+ 
  geom_point(aes(color=group)) +
  theme_bw() + 
  theme(axis.text.x=element_text(angle=90,size=8, vjust = 1, hjust = 0.5))
dfstat_p

ggsave("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITS_overall/ITS_p_s_overallstat_2lines.pdf",dfstat_p, width = 12, height = 5)



# resistant 
ITS_p_r_genestat = list()
dfstat_p_r = list()
for(i in 1:length(ITS_R)){
  # i = 1
  # i = 21
  ITSname <- ITS_R[i]
  listinput <- lapply(ITS_p_set, function(x){x[[ITSname]]})
  genelist = unique(unlist(listinput))
  listinput_df <- fromList(listinput)
  rownames(listinput_df) <- genelist
  table(apply(listinput_df,1,sum))
  apply(listinput_df,2,sum)

  ITS_p_r_genestat[[i]] <- listinput_df
  dfstat_p_r[[i]] <- data.frame(immune_sig = ITSname,
                            group = c('gene_cancerspecific', 'gene_2cancer','gene_3cancer','50percent_cancer_shared'),
                            num_ITSgene = nrow(listinput_df),
                            num_gene = c(length(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)==1)]),
                                      length(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)==2)]),
                                      length(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)==3)]),
                                      length(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)>=0.5*22)])),
                            value = c(length(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)==1)])/ nrow(listinput_df),
                                      length(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)==2)])/ nrow(listinput_df),
                                      length(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)==3)])/ nrow(listinput_df),
                                      length(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)>=0.5*22)])/nrow(listinput_df)))
  

}



names(ITS_p_r_genestat) <- ITS_R # names(ITS_p_set[[1]])
# names(dfstat) <- names(ITS_p_set[[1]])
save(ITS_p_r_genestat, dfstat, file = "07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITS_overall/ITS_p_r_genestat.Rdata")

dfstat_p_r_all <- do.call(rbind, dfstat_p_r)
write.csv(dfstat_p_r_all, "07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITS_overall/ITS_p_r_genestat.csv", quote = F, row.names = F)


dfstat_p_r_all1 <- dfstat_p_r_all[dfstat_p_r_all$group == 'gene_cancerspecific',]
dfstat_p_r_all1 <- dfstat_p_r_all1[order(dfstat_p_r_all1$value, decreasing = T),]
dfstat_p_r_all$immune_sig <- factor(dfstat_p_r_all$immune_sig, levels = dfstat_p_r_all1$immune_sig)

dfstat_p <- ggplot(data = dfstat_p_r_all, aes(x = immune_sig, y = value, group=group)) + 
  geom_line(aes(color=group))+ 
  geom_point(aes(color=group)) +
  theme_bw() + 
  theme(axis.text.x=element_text(angle=90,size=8, vjust = 1, hjust = 0.5))
dfstat_p
ggsave("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITS_overall/ITS_p_r_overallstat_4lines.pdf",dfstat_p, width = 13, height = 10)

dfstat_p_r_all2 <- dfstat_p_r_all[dfstat_p_r_all$group != "gene_3cancer", ]
dfstat_p <- ggplot(data = dfstat_p_r_all2, aes(x = immune_sig, y = value, group=group)) + 
  geom_line(aes(color=group))+ 
  geom_point(aes(color=group)) +
  theme_bw() + 
  theme(axis.text.x=element_text(angle=90,size=8, vjust = 1, hjust = 0.5))
dfstat_p
ggsave("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITS_overall/ITS_p_r_overallstat_3lines.pdf", dfstat_p, width = 13, height = 10)


dfstat_p_r = list()
for(i in 1:length(ITS_R)){
  # i = 1
  # i = 21
  ITSname <- ITS_R[i]
  listinput <- lapply(ITS_p_set, function(x){x[[ITSname]]})
  genelist = unique(unlist(listinput))
  listinput_df <- fromList(listinput)
  rownames(listinput_df) <- genelist
  table(apply(listinput_df,1,sum))
  apply(listinput_df,2,sum)

  ITS_p_genestat[[i]] <- listinput_df
  dfstat_p_r[[i]] <- data.frame(immune_sig = ITSname,
                            group = c('gene_cancerspecific', 'gene_2cancer','gene_3cancer','50percent_cancer_shared'),
                            num_gene = c(length(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)==1)]),
                                      length(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)>=2)]),
                                      length(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)>=3)]),
                                      length(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)>= 0.5* length(listinput))])),
                            value = c(length(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)==1)])/ nrow(listinput_df),
                                      length(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)>=2)])/ nrow(listinput_df),
                                      length(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)>=3)])/ nrow(listinput_df),
                                      length(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)>=0.5*length(listinput))])/ nrow(listinput_df)))
  
}

dfstat_p_r_all <- do.call(rbind, dfstat_p_r)
write.csv(dfstat_p_r_all, "07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITS_overall/ITS_p_r_genestat_binary.csv", quote = F, row.names = F)


dfstat_p_r_all1 <- dfstat_p_r_all[dfstat_p_r_all$group == 'gene_cancerspecific',]
dfstat_p_r_all1 <- dfstat_p_r_all1[order(dfstat_p_r_all1$value, decreasing = T),]
dfstat_p_r_all$immune_sig <- factor(dfstat_p_r_all$immune_sig, levels = dfstat_p_r_all1$immune_sig)
dfstat_p_r_all2 <- dfstat_p_r_all[is.element(dfstat_p_r_all$group, c('gene_cancerspecific', 'gene_2cancer')), ]
dfstat_p <- ggplot(data = dfstat_p_r_all2, aes(x = immune_sig, y = value, group=group)) + 
  geom_line(aes(color=group))+ 
  geom_point(aes(color=group)) +
  theme_bw() + 
  theme(axis.text.x=element_text(angle=90,size=8, vjust = 1, hjust = 0.5))
dfstat_p

ggsave("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITS_overall/ITS_p_r_overallstat_2lines.pdf",dfstat_p, width = 10, height = 5)
ggsave("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITS_overall/ITS_p_r_overallstat_2lines_v2.pdf",dfstat_p, width = 3.5, height = 5)



# ITS_N ------------------------------------------------------------------------


# sensitive 

ITS_n_s_genestat = list()
dfstat_n_s = list()
for(i in 1:length(ITS_S)){
  # i = 1
  # i = 21
  ITSname <- ITS_S[i]
  listinput <- lapply(ITS_n_set, function(x){x[[ITSname]]})
  genelist = unique(unlist(listinput))
  listinput_df <- fromList(listinput)
  rownames(listinput_df) <- genelist
  table(apply(listinput_df,1,sum))
  apply(listinput_df,2,sum)

  ITS_n_s_genestat[[i]] <- listinput_df
  dfstat_n_s[[i]] <- data.frame(immune_sig = ITSname,
                            group = c('gene_cancerspecific', 'gene_2cancer','gene_3cancer','50percent_cancer_shared'),
                            num_ITSgene = nrow(listinput_df),
                            num_gene = c(length(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)==1)]),
                                      length(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)==2)]),
                                      length(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)==3)]),
                                      length(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)>=0.5*22)])),
                            value = c(length(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)==1)])/ nrow(listinput_df),
                                      length(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)==2)])/ nrow(listinput_df),
                                      length(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)==3)])/ nrow(listinput_df),
                                      length(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)>=0.5*22)])/nrow(listinput_df)))
  

}



names(ITS_n_s_genestat) <- ITS_S # names(ITS_n_set[[1]])
# names(dfstat) <- names(ITS_n_set[[1]])
save(ITS_n_s_genestat, dfstat, file = "07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITS_overall/ITS_n_s_genestat.Rdata")

dfstat_n_s_all <- do.call(rbind, dfstat_n_s)
write.csv(dfstat_n_s_all, "07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITS_overall/ITS_n_s_genestat.csv", quote = F, row.names = F)


dfstat_n_s_all1 <- dfstat_n_s_all[dfstat_n_s_all$group == 'gene_cancerspecific',]
dfstat_n_s_all1 <- dfstat_n_s_all1[order(dfstat_n_s_all1$value, decreasing = T),]
dfstat_n_s_all$immune_sig <- factor(dfstat_n_s_all$immune_sig, levels = dfstat_n_s_all1$immune_sig)

dfstat_p <- ggplot(data = dfstat_n_s_all, aes(x = immune_sig, y = value, group=group)) + 
  geom_line(aes(color=group))+ 
  geom_point(aes(color=group)) +
  theme_bw() + 
  theme(axis.text.x=element_text(angle=90,size=8, vjust = 1, hjust = 0.5))
dfstat_p
ggsave("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITS_overall/ITS_n_s_overallstat_4lines.pdf",dfstat_p, width = 13, height = 10)

dfstat_n_s_all2 <- dfstat_n_s_all[dfstat_n_s_all$group != "gene_3cancer", ]
dfstat_p <- ggplot(data = dfstat_n_s_all2, aes(x = immune_sig, y = value, group=group)) + 
  geom_line(aes(color=group))+ 
  geom_point(aes(color=group)) +
  theme_bw() + 
  theme(axis.text.x=element_text(angle=90,size=8, vjust = 1, hjust = 0.5))
dfstat_p
ggsave("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITS_overall/ITS_n_s_overallstat_3lines.pdf", dfstat_p, width = 13, height = 10)


dfstat_n_s = list()
for(i in 1:length(ITS_S)){
  # i = 1
  # i = 21
  ITSname <- ITS_S[i]
  listinput <- lapply(ITS_n_set, function(x){x[[ITSname]]})
  genelist = unique(unlist(listinput))
  listinput_df <- fromList(listinput)
  rownames(listinput_df) <- genelist
  table(apply(listinput_df,1,sum))
  apply(listinput_df,2,sum)

  ITS_n_genestat[[i]] <- listinput_df
  dfstat_n_s[[i]] <- data.frame(immune_sig = ITSname,
                            group = c('gene_cancerspecific', 'gene_2cancer','gene_3cancer','50percent_cancer_shared'),
                            num_gene = c(length(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)==1)]),
                                      length(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)>=2)]),
                                      length(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)>=3)]),
                                      length(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)>= 0.5* length(listinput))])),
                            value = c(length(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)==1)])/ nrow(listinput_df),
                                      length(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)>=2)])/ nrow(listinput_df),
                                      length(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)>=3)])/ nrow(listinput_df),
                                      length(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)>=0.5*length(listinput))])/ nrow(listinput_df)))
  
}

dfstat_n_s_all <- do.call(rbind, dfstat_n_s)
write.csv(dfstat_n_s_all, "07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITS_overall/ITS_n_s_genestat_binary.csv", quote = F, row.names = F)


dfstat_n_s_all1 <- dfstat_n_s_all[dfstat_n_s_all$group == 'gene_cancerspecific',]
dfstat_n_s_all1 <- dfstat_n_s_all1[order(dfstat_n_s_all1$value, decreasing = T),]
dfstat_n_s_all$immune_sig <- factor(dfstat_n_s_all$immune_sig, levels = dfstat_n_s_all1$immune_sig)
dfstat_n_s_all2 <- dfstat_n_s_all[is.element(dfstat_n_s_all$group, c('gene_cancerspecific', 'gene_2cancer')), ]
dfstat_p <- ggplot(data = dfstat_n_s_all2, aes(x = immune_sig, y = value, group=group)) + 
  geom_line(aes(color=group))+ 
  geom_point(aes(color=group)) +
  theme_bw() + 
  theme(axis.text.x=element_text(angle=90,size=8, vjust = 1, hjust = 0.5))
dfstat_p

ggsave("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITS_overall/ITS_n_s_overallstat_2lines.pdf",dfstat_p, width = 13, height = 10)



# resistant 
ITS_n_r_genestat = list()
dfstat_n_r = list()
for(i in 1:length(ITS_R)){
  # i = 1
  # i = 21
  ITSname <- ITS_R[i]
  listinput <- lapply(ITS_n_set, function(x){x[[ITSname]]})
  genelist = unique(unlist(listinput))
  listinput_df <- fromList(listinput)
  rownames(listinput_df) <- genelist
  table(apply(listinput_df,1,sum))
  apply(listinput_df,2,sum)

  ITS_n_r_genestat[[i]] <- listinput_df
  dfstat_n_r[[i]] <- data.frame(immune_sig = ITSname,
                            group = c('gene_cancerspecific', 'gene_2cancer','gene_3cancer','50percent_cancer_shared'),
                            num_ITSgene = nrow(listinput_df),
                            num_gene = c(length(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)==1)]),
                                      length(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)==2)]),
                                      length(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)==3)]),
                                      length(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)>=0.5*22)])),
                            value = c(length(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)==1)])/ nrow(listinput_df),
                                      length(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)==2)])/ nrow(listinput_df),
                                      length(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)==3)])/ nrow(listinput_df),
                                      length(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)>=0.5*22)])/nrow(listinput_df)))
  

}



names(ITS_n_r_genestat) <- ITS_R# names(ITS_n_set[[1]])
# names(dfstat) <- names(ITS_n_set[[1]])
save(ITS_n_r_genestat, dfstat, file = "07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITS_overall/ITS_n_r_genestat.Rdata")

dfstat_n_r_all <- do.call(rbind, dfstat_n_r)
write.csv(dfstat_n_r_all, "07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITS_overall/ITS_n_r_genestat.csv", quote = F, row.names = F)


dfstat_n_r_all1 <- dfstat_n_r_all[dfstat_n_r_all$group == 'gene_cancerspecific',]
dfstat_n_r_all1 <- dfstat_n_r_all1[order(dfstat_n_r_all1$value, decreasing = T),]
dfstat_n_r_all$immune_sig <- factor(dfstat_n_r_all$immune_sig, levels = dfstat_n_r_all1$immune_sig)

dfstat_p <- ggplot(data = dfstat_n_r_all, aes(x = immune_sig, y = value, group=group)) + 
  geom_line(aes(color=group))+ 
  geom_point(aes(color=group)) +
  theme_bw() + 
  theme(axis.text.x=element_text(angle=90,size=8, vjust = 1, hjust = 0.5))
dfstat_p
ggsave("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITS_overall/ITS_n_r_overallstat_4lines.pdf",dfstat_p, width = 13, height = 10)

dfstat_n_r_all2 <- dfstat_n_r_all[dfstat_n_r_all$group != "gene_3cancer", ]
dfstat_p <- ggplot(data = dfstat_n_r_all2, aes(x = immune_sig, y = value, group=group)) + 
  geom_line(aes(color=group))+ 
  geom_point(aes(color=group)) +
  theme_bw() + 
  theme(axis.text.x=element_text(angle=90,size=8, vjust = 1, hjust = 0.5))
dfstat_p
ggsave("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITS_overall/ITS_n_r_overallstat_3lines.pdf", dfstat_p, width = 13, height = 10)


dfstat_n_r = list()
for(i in 1:length(ITS_R)){
  # i = 1
  # i = 21
  ITSname <- ITS_R[i]
  listinput <- lapply(ITS_n_set, function(x){x[[ITSname]]})
  genelist = unique(unlist(listinput))
  listinput_df <- fromList(listinput)
  rownames(listinput_df) <- genelist
  table(apply(listinput_df,1,sum))
  apply(listinput_df,2,sum)

  ITS_n_genestat[[i]] <- listinput_df
  dfstat_n_r[[i]] <- data.frame(immune_sig = ITSname,
                            group = c('gene_cancerspecific', 'gene_2cancer','gene_3cancer','50percent_cancer_shared'),
                            num_gene = c(length(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)==1)]),
                                      length(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)>=2)]),
                                      length(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)>=3)]),
                                      length(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)>= 0.5* length(listinput))])),
                            value = c(length(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)==1)])/ nrow(listinput_df),
                                      length(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)>=2)])/ nrow(listinput_df),
                                      length(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)>=3)])/ nrow(listinput_df),
                                      length(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)>=0.5*length(listinput))])/ nrow(listinput_df)))
  
}

dfstat_n_r_all <- do.call(rbind, dfstat_n_r)
write.csv(dfstat_n_r_all, "07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITS_overall/ITS_n_r_genestat_binary.csv", quote = F, row.names = F)


dfstat_n_r_all1 <- dfstat_n_r_all[dfstat_n_r_all$group == 'gene_cancerspecific',]
dfstat_n_r_all1 <- dfstat_n_r_all1[order(dfstat_n_r_all1$value, decreasing = T),]
dfstat_n_r_all$immune_sig <- factor(dfstat_n_r_all$immune_sig, levels = dfstat_n_r_all1$immune_sig)
dfstat_n_r_all2 <- dfstat_n_r_all[is.element(dfstat_n_r_all$group, c('gene_cancerspecific', 'gene_2cancer')), ]
dfstat_p <- ggplot(data = dfstat_n_r_all2, aes(x = immune_sig, y = value, group=group)) + 
  geom_line(aes(color=group))+ 
  geom_point(aes(color=group)) +
  theme_bw() + 
  theme(axis.text.x=element_text(angle=90,size=8, vjust = 1, hjust = 0.5))
dfstat_p

ggsave("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITS_overall/ITS_n_r_overallstat_2lines.pdf",dfstat_p, width = 13, height = 10)

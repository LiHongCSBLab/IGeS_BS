# rm(list=ls())
work_path <- "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/"
setwd(work_path)
options(stringsAsFactors = F)


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

# ----------------------------------------------------------------------------
dir.create("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/")
dir.create("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/")
dir.create("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITS_overall/")

# LOAD GENES IN LINCS ------------------------------------------------------
lincs_gene = read.csv("/picb/bigdata/project/FengFYM/drug_MoA_reanalysis/LINCS_data/GSE70138/file/GSE92742_Broad_LINCS_gene_info.txt", sep = "\t")

# ------------------------------------------------------------------------------
load( "07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITS_set.Rdata")


# ------------------------------------------------------------------------------
# ITS positive
# ------------------------------------------------------------------------------
# stat: sensitive-related and resistant-related ITS ----------------------------
# ------------------------------------------------------------------------------
# meta_result <- read.csv("05_immuneSig_ICBresponse/results/meta_results/res_metaRes.csv")
# meta_selected <- meta_result[meta_result$Meta_Pval < 0.05, ]
# meta_selected$flag = 'sensitive'
# meta_selected[meta_selected$OR < 1, ]$flag = 'resistant'
immunesig <- read.csv("05_immuneSig_ICBresponse/results/meta_results/res_metaRes.csv")
immunesig <- immunesig[immunesig$Meta_Pval < 0.05, ]
names(immunesig) <- c("immune_sig", names(immunesig)[-1])
immunesig$flag <- "sensitive"
immunesig[immunesig$OR < 1,]$flag <- "resistant"
immunesig$setindex = "01_sensitive"
immunesig[immunesig$flag == "resistant", ]$setindex = "02_resistant"
load("05_immuneSig_ICBresponse/results_r0.4_ITSssgsea_CombatRef/finalmodel_17datasets_ssgseanorm/metaAnalysis_filter/en/elasticnet_model_pred.Rdata")
en_model_para = m_en$learner.model$glmnet.fit
m_en_imp <- data.frame(en_model_para$beta[,which(en_model_para$lambda == m_en$learner.model$lambda.min)])
names(m_en_imp) <- "weight"
m_en_imp$immune_sig = row.names(m_en_imp)
immunesig = merge(immunesig, m_en_imp, by = "immune_sig")
immunesig =immunesig[immunesig$weight != 0 ,]
immunesig = immunesig[order(immunesig$flag),]

ITS_S = immunesig[immunesig$flag == 'sensitive', ]$immune_sig
ITS_R = immunesig[immunesig$flag == 'resistant', ]$immune_sig



# sensitive 
ITS_p_s_genestat <- list()
ITS_p_s_cancerspecific = list()
ITS_p_s_notspecific = list()

for(i in 1:length(ITS_S)){
  # i = 1
  ITSname <- ITS_S[i]
  listinput <- lapply(ITS_p_set, function(x){x[[ITSname]]})
  genelist = unique(unlist(listinput))
  listinput_df <- fromList(listinput)
  rownames(listinput_df) <- genelist
  # table(apply(listinput_df,1,sum))
  # apply(listinput_df,2,sum)
  # ITS_p_s_cancerspecific[[i]] =  unique(names(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)==1)]))
  # ITS_p_s_notspecific[[i]] = unique(names(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)>=2)]))
  ITS_p_s_genestat[[i]] <- listinput_df
  
  
  df_stat <- listinput_df # * 22
  df_stat_row <- rev(sort(rowSums(sign(df_stat))))
  df_stat_row <- as.data.frame(df_stat_row)
  names(df_stat_row) = "num_Gene"
  df_stat_row$gene_name = rownames(df_stat_row)
  df_stat_row$immune_sig <- ITSname
  
  
  ITS_p_s_notspecific[[i]] <- df_stat_row[df_stat_row$num_Gene > 1,]
  
  tmp <- df_stat_row[df_stat_row$num_Gene == 1,]
  if(nrow(tmp) > 0){
    tmp$cancer_type <- NA
    for(g in tmp$gene_name){
      tmp[tmp$gene_name == g,]$cancer_type <- names(listinput_df)[listinput_df[rownames(listinput_df)==g,]==1]
    }
    ITS_p_s_cancerspecific[[i]] <- tmp
  }else{
    ITS_p_s_cancerspecific[[i]] <- NULL
  }
  
}


df <- do.call(rbind, ITS_p_s_cancerspecific)

gs_s_specific <- list()
gs_s_nonspecific <- list()

for(i in seq(length(unique(df$gene_name)))) {
    g=unique(df$gene_name)[i]
    x = df[df$gene_name == g, ]
    if(length(unique(x$cancer_type))>1){
        gs_s_nonspecific[[i]] = x
    }else{
        gs_s_specific[[i]] = x
    }
}

gs_s_specific <- unique(do.call(rbind, gs_s_specific))
gs_s_nonspecific <- unique(do.call(rbind, gs_s_nonspecific))

# ITS_p_s_notspecific: gene in any ITS in >= 2 cancer type 
# ITS_p_s_cancerspecific: gene in any ITS in == 1 cancer type 
# # gs_s_nonspecific: ITS (which contain gene) in == 1 cancer type, but gene can show in > 2 ITS in > 2 cancer types
# # gs_s_specific: ITS (which contain gene) in == 1 cancer type, but gene can show in >2 ITS in ==2 cancer types

# Venn plot --------------------------------------------------------------------
# type1: ITS_p_s_notspecific: gene in any ITS in >= 2 cancer type 
ITS_p_s_notspecific_df <- do.call(rbind, ITS_p_s_notspecific)

# # type2: gs_nonspecific: ITS (which contain gene) in == 1 cancer type, but gene can show in > 2 ITS in > 2 cancer types
# gs_s_nonspecific
# # type3: gs_specificITS (which contain gene) in == 1 cancer type, but gene can show in >2 ITS in ==2 cancer types
# gs_s_specific
# type2 + type3 =ITS_p_s_cancerspecific

library(VennDiagram)
venn.diagram(x=list(type1_notspecific=unique(ITS_p_s_notspecific_df$gene_name), 
                     type2_ITSspecific_g_not=unique(gs_s_nonspecific$gene_name), 
                     type3_specific=unique(gs_s_specific$gene_name)), 
 filename = "07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITS_overall/ITS_p_s_gene_groups_Venn.png", 
 height = 450, width = 450, resolution =300, imagetype="png", 
 col="white", fill=c(colors()[616], colors()[38], colors()[468]), 
 alpha=c(0.6, 0.6, 0.6), lwd=c(1, 1, 1), cex=0.3, 
 cat.dist=c(-0.07, -0.07, -0.05), cat.pos=c(300, 60, 180), cat.cex=0.3) 


gene_p_s_specific = setdiff(unique(gs_s_specific$gene_name), unique(ITS_p_s_notspecific_df$gene_name))
gene_p_s_nonspecific = union(unique(gs_s_nonspecific$gene_name), unique(ITS_p_s_notspecific_df$gene_name))

ITS_p_s_cancerspecific_df = do.call(rbind, ITS_p_s_cancerspecific)
ITS_p_s_notspecific_df$cancer_type = 'pancancer'
ITS_p_s_genestat = rbind(ITS_p_s_notspecific_df, ITS_p_s_cancerspecific_df[])

ITS_p_s_genestat_specific = ITS_p_s_genestat[is.element(ITS_p_s_genestat$gene_name, gene_p_s_specific),]
ITS_p_s_genestat_nonspecific = ITS_p_s_genestat[is.element(ITS_p_s_genestat$gene_name, gene_p_s_nonspecific),]

length(intersect(unique(ITS_p_s_genestat_specific$gene_name), unique(ITS_p_s_genestat_nonspecific$gene_name)))

save(ITS_p_s_genestat, ITS_p_s_notspecific, 
     ITS_p_s_cancerspecific, gs_s_specific, gs_s_nonspecific, 
     ITS_p_s_genestat_specific, ITS_p_s_genestat_nonspecific,
     file = "07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITS_p_s_ITSstat.Rdata")

# overall-pie plot -------------------------------------------------------------
library(scales)
blank_theme <- theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )



library(RColorBrewer)
df_pie = data.frame(group = c("type1_notspecific",  "type3_specific"),
                    value = c(length(unique(gene_p_s_nonspecific)), 
                              length(unique(gene_p_s_specific))),
                    col = c(brewer.pal(3, "Set3")[1:2]))
pie <- ggplot(df_pie, aes(x="", y=value, fill=group))+
  geom_bar( width = 1, stat = "identity", alpha = 0.5) + 
  coord_polar("y", direction=1) +
  scale_fill_manual(values=brewer.pal(3, "Set3"))+
  geom_text(aes(y = sum(df_pie$value)-cumsum(df_pie$value)+df_pie$value/2, 
                label = paste(value,"\n", percent(value/sum(value)))), size=5)+
  blank_theme +
  theme(axis.text.x=element_blank())
pie
ggsave("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITS_overall/ITS_p_s_gene_groups.pdf", pie, width = 8, height = 8)



# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# resistant 
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

ITS_p_r_genestat <- list()
ITS_p_r_cancerspecific = list()
ITS_p_r_notspecific = list()

for(i in 1:length(ITS_R)){
  # i = 1
  ITSname <- ITS_R[i]
  listinput <- lapply(ITS_p_set, function(x){x[[ITSname]]})
  genelist = unique(unlist(listinput))
  listinput_df <- fromList(listinput)
  rownames(listinput_df) <- genelist
  # table(apply(listinput_df,1,sum))
  # apply(listinput_df,2,sum)
  # ITS_p_r_cancerspecific[[i]] =  unique(names(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)==1)]))
  # ITS_p_r_notspecific[[i]] = unique(names(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)>=2)]))
  ITS_p_r_genestat[[i]] <- listinput_df
  
  
  df_stat <- listinput_df # * 22
  df_stat_row <- rev(sort(rowSums(sign(df_stat))))
  df_stat_row <- as.data.frame(df_stat_row)
  names(df_stat_row) = "num_Gene"
  df_stat_row$gene_name = rownames(df_stat_row)
  df_stat_row$immune_sig <- ITSname
  
  
  ITS_p_r_notspecific[[i]] <- df_stat_row[df_stat_row$num_Gene > 1,]
  
  tmp <- df_stat_row[df_stat_row$num_Gene == 1,]
  if(nrow(tmp) > 0){
    tmp$cancer_type <- NA
    for(g in tmp$gene_name){
      tmp[tmp$gene_name == g,]$cancer_type <- names(listinput_df)[listinput_df[rownames(listinput_df)==g,]==1]
    }
    ITS_p_r_cancerspecific[[i]] <- tmp
  }else{
    ITS_p_r_cancerspecific[[i]] <- NULL
  }
  
}



df <- do.call(rbind, ITS_p_r_cancerspecific)

gs_r_specific <- list()
gs_r_nonspecific <- list()

for(i in seq(length(unique(df$gene_name)))) {
    g=unique(df$gene_name)[i]
    tmp = df[df$gene_name == g, ]
    if(length(unique(tmp$cancer_type))>1){
        gs_r_nonspecific[[i]] = tmp
    }else{
        gs_r_specific[[i]] = tmp
    }
}

gs_r_specific <- unique(do.call(rbind, gs_r_specific))
gs_r_nonspecific <- unique(do.call(rbind, gs_r_nonspecific))



# Venn plot --------------------------------------------------------------------
# type1: ITS_p_r_notspecific: gene in any ITS in >= 2 cancer type 
ITS_p_r_notspecific_df <- do.call(rbind, ITS_p_r_notspecific)

# # type2: gs_nonspecific: ITS (which contain gene) in == 1 cancer type, but gene can show in > 2 ITS in > 2 cancer types
# gs_r_nonspecific
# # type3: gs_specificITS (which contain gene) in == 1 cancer type, but gene can show in >2 ITS in ==2 cancer types
# gs_r_specific
# type2 + type3 =ITS_p_r_cancerspecific

library(VennDiagram)
venn.diagram(x=list(type1_notspecific=unique(ITS_p_r_notspecific_df$gene_name), 
                     type2_ITSspecific_g_not=unique(gs_r_nonspecific$gene_name), 
                     type3_specific=unique(gs_r_specific$gene_name)), 
 filename = "07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITS_overall/ITS_p_r_gene_groups_Venn.png", 
 height = 450, width = 450, resolution =300, imagetype="png", 
 col="white", fill=c(colors()[616], colors()[38], colors()[468]), 
 alpha=c(0.6, 0.6, 0.6), lwd=c(1, 1, 1), cex=0.3, 
 cat.dist=c(-0.07, -0.07, -0.05), cat.pos=c(300, 60, 180), cat.cex=0.3) 


gene_p_r_specific = setdiff(unique(gs_r_specific$gene_name), unique(ITS_p_r_notspecific_df$gene_name))
gene_p_r_nonspecific = union(unique(gs_r_nonspecific$gene_name), unique(ITS_p_r_notspecific_df$gene_name))

ITS_p_r_cancerspecific_df = do.call(rbind, ITS_p_r_cancerspecific)
ITS_p_r_notspecific_df$cancer_type = 'pancancer'
ITS_p_r_genestat = rbind(ITS_p_r_notspecific_df, ITS_p_r_cancerspecific_df[])

ITS_p_r_genestat_specific = ITS_p_r_genestat[is.element(ITS_p_r_genestat$gene_name, gene_p_r_specific),]
ITS_p_r_genestat_nonspecific = ITS_p_r_genestat[is.element(ITS_p_r_genestat$gene_name, gene_p_r_nonspecific),]


length(unique(ITS_p_r_genestat_specific$gene_name))
length(unique(ITS_p_r_genestat_nonspecific$gene_name))
length(intersect(unique(ITS_p_r_genestat_specific$gene_name), unique(ITS_p_r_genestat_nonspecific$gene_name)))

save(ITS_p_r_genestat, ITS_p_r_notspecific, 
     ITS_p_r_cancerspecific, gs_r_specific, gs_r_nonspecific, 
     ITS_p_r_genestat_specific, ITS_p_r_genestat_nonspecific,
     file = "07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITS_p_r_ITSstat.Rdata")

# overall-pie plot -------------------------------------------------------------
library(scales)
blank_theme <- theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )



library(RColorBrewer)
df_pie = data.frame(group = c("type1_notspecific",  "type3_specific"),
                    value = c(length(unique(gene_p_r_nonspecific)), 
                              length(unique(gene_p_r_specific))),
                    col = c(brewer.pal(3, "Set3")[1:2]))
pie <- ggplot(df_pie, aes(x="", y=value, fill=group))+
  geom_bar( width = 1, stat = "identity", alpha = 0.5) + 
  coord_polar("y", direction=1) +
  scale_fill_manual(values=brewer.pal(3, "Set3"))+
  geom_text(aes(y = sum(df_pie$value)-cumsum(df_pie$value)+df_pie$value/2, 
                label = paste(value,"\n", percent(value/sum(value)))), size=5)+
  blank_theme +
  theme(axis.text.x=element_blank())
pie
ggsave("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITS_overall/ITS_p_r_gene_groups.pdf", pie, width = 8, height = 8)




# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ITS negative
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# sensitive 
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

ITS_n_s_genestat <- list()
ITS_n_s_cancerspecific = list()
ITS_n_s_notspecific = list()

for(i in 1:length(ITS_S)){
  # i = 1
  ITSname <- ITS_S[i]
  listinput <- lapply(ITS_n_set, function(x){x[[ITSname]]})
  genelist = unique(unlist(listinput))
  if(length(genelist)>0){

  listinput_df <- fromList(listinput)
  rownames(listinput_df) <- genelist
  # table(apply(listinput_df,1,sum))
  # apply(listinput_df,2,sum)
  # ITS_n_s_cancerspecific[[i]] =  unique(names(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)==1)]))
  # ITS_n_s_notspecific[[i]] = unique(names(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)>=2)]))
  ITS_n_s_genestat[[i]] <- listinput_df
  
  
  df_stat <- listinput_df # * 22
  df_stat_row <- rev(sort(rowSums(sign(df_stat))))
  df_stat_row <- as.data.frame(df_stat_row)
  names(df_stat_row) = "num_Gene"
  df_stat_row$gene_name = rownames(df_stat_row)
  df_stat_row$immune_sig <- ITSname
  
  
  ITS_n_s_notspecific[[i]] <- df_stat_row[df_stat_row$num_Gene > 1,]
  
  tmp <- df_stat_row[df_stat_row$num_Gene == 1,]
  if(nrow(tmp) > 0){
    tmp$cancer_type <- NA
    for(g in tmp$gene_name){
      tmp[tmp$gene_name == g,]$cancer_type <- names(listinput_df)[listinput_df[rownames(listinput_df)==g,]==1]
    }
    ITS_n_s_cancerspecific[[i]] <- tmp
  }else{
    ITS_n_s_cancerspecific[[i]] <- NULL
  }
 }
}

df <- do.call(rbind, ITS_n_s_cancerspecific)

gs_s_specific <- list()
gs_s_nonspecific <- list()

for(i in seq(length(unique(df$gene_name)))) {
    g=unique(df$gene_name)[i]
    tmp = df[df$gene_name == g, ]
    if(length(unique(tmp$cancer_type))>1){
        gs_s_nonspecific[[i]] = tmp
    }else{
        gs_s_specific[[i]] = tmp
    }
}

gs_s_specific <- unique(do.call(rbind, gs_s_specific))
gs_s_nonspecific <- unique(do.call(rbind, gs_s_nonspecific))

# ITS_n_s_notspecific: gene in any ITS in >= 2 cancer type 
# ITS_n_s_cancerspecific: gene in any ITS in == 1 cancer type 
# # gs_s_nonspecific: ITS (which contain gene) in == 1 cancer type, but gene can show in > 2 ITS in > 2 cancer types
# # gs_s_specific: ITS (which contain gene) in == 1 cancer type, but gene can show in >2 ITS in ==2 cancer types

# Venn plot --------------------------------------------------------------------
# type1: ITS_n_s_notspecific: gene in any ITS in >= 2 cancer type 
ITS_n_s_notspecific_df <- do.call(rbind, ITS_n_s_notspecific)

# # type2: gs_nonspecific: ITS (which contain gene) in == 1 cancer type, but gene can show in > 2 ITS in > 2 cancer types
# gs_s_nonspecific
# # type3: gs_specificITS (which contain gene) in == 1 cancer type, but gene can show in >2 ITS in ==2 cancer types
# gs_s_specific
# type2 + type3 =ITS_n_s_cancerspecific

library(VennDiagram)
venn.diagram(x=list(type1_notspecific=unique(ITS_n_s_notspecific_df$gene_name), 
                     type2_ITSspecific_g_not=unique(gs_s_nonspecific$gene_name), 
                     type3_specific=unique(gs_s_specific$gene_name)), 
 filename = "07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITS_overall/ITS_n_s_gene_groups_Venn.png", 
 height = 450, width = 450, resolution =300, imagetype="png", 
 col="white", fill=c(colors()[616], colors()[38], colors()[468]), 
 alpha=c(0.6, 0.6, 0.6), lwd=c(1, 1, 1), cex=0.3, 
 cat.dist=c(-0.07, -0.07, -0.05), cat.pos=c(300, 60, 180), cat.cex=0.3) 


gene_n_s_specific = setdiff(unique(gs_s_specific$gene_name), unique(ITS_n_s_notspecific_df$gene_name))
gene_n_s_nonspecific = union(unique(gs_s_nonspecific$gene_name), unique(ITS_n_s_notspecific_df$gene_name))

ITS_n_s_cancerspecific_df = do.call(rbind, ITS_n_s_cancerspecific)
ITS_n_s_notspecific_df$cancer_type = 'pancancer'
ITS_n_s_genestat = rbind(ITS_n_s_notspecific_df, ITS_n_s_cancerspecific_df[])

ITS_n_s_genestat_specific = ITS_n_s_genestat[is.element(ITS_n_s_genestat$gene_name, gene_n_s_specific),]
ITS_n_s_genestat_nonspecific = ITS_n_s_genestat[is.element(ITS_n_s_genestat$gene_name, gene_n_s_nonspecific),]

length(unique(ITS_n_s_genestat_specific$gene_name))
length(unique(ITS_n_s_genestat_nonspecific$gene_name))
length(intersect(unique(ITS_n_s_genestat_specific$gene_name), unique(ITS_n_s_genestat_nonspecific$gene_name)))

save(ITS_n_s_genestat, ITS_n_s_notspecific, 
     ITS_n_s_cancerspecific, gs_s_specific, gs_s_nonspecific, 
     ITS_n_s_genestat_specific, ITS_n_s_genestat_nonspecific,
     file = "07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITS_n_s_ITSstat.Rdata")

# overall-pie plot -------------------------------------------------------------
library(scales)
blank_theme <- theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )



library(RColorBrewer)
df_pie = data.frame(group = c("type1_notspecific",  "type3_specific"),
                    value = c(length(unique(gene_n_s_nonspecific)), 
                              length(unique(gene_n_s_specific))),
                    col = c(brewer.pal(3, "Set3")[1:2]))
pie <- ggplot(df_pie, aes(x="", y=value, fill=group))+
  geom_bar( width = 1, stat = "identity", alpha = 0.5) + 
  coord_polar("y", direction=1) +
  scale_fill_manual(values=brewer.pal(3, "Set3"))+
  geom_text(aes(y = sum(df_pie$value)-cumsum(df_pie$value)+df_pie$value/2, 
                label = paste(value,"\n", percent(value/sum(value)))), size=5)+
  blank_theme +
  theme(axis.text.x=element_blank())
pie
ggsave("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITS_overall/ITS_n_s_gene_groups.pdf", pie, width = 8, height = 8)


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# resistant 
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

ITS_n_r_genestat <- list()
ITS_n_r_cancerspecific = list()
ITS_n_r_notspecific = list()

for(i in 1:length(ITS_R)){
  # i = 1
  ITSname <- ITS_R[i]
  listinput <- lapply(ITS_n_set, function(x){x[[ITSname]]})
  genelist = unique(unlist(listinput))
  if(length(genelist)>0){
  listinput_df <- fromList(listinput)
  rownames(listinput_df) <- genelist
  # table(apply(listinput_df,1,sum))
  # apply(listinput_df,2,sum)
  # ITS_n_r_cancerspecific[[i]] =  unique(names(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)==1)]))
  # ITS_n_r_notspecific[[i]] = unique(names(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)>=2)]))
  ITS_n_r_genestat[[i]] <- listinput_df
  
  
  df_stat <- listinput_df # * 22
  df_stat_row <- rev(sort(rowSums(sign(df_stat))))
  df_stat_row <- as.data.frame(df_stat_row)
  names(df_stat_row) = "num_Gene"
  df_stat_row$gene_name = rownames(df_stat_row)
  df_stat_row$immune_sig <- ITSname
  
  
  ITS_n_r_notspecific[[i]] <- df_stat_row[df_stat_row$num_Gene > 1,]
  
  tmp <- df_stat_row[df_stat_row$num_Gene == 1,]
  if(nrow(tmp) > 0){
    tmp$cancer_type <- NA
    for(g in tmp$gene_name){
      tmp[tmp$gene_name == g,]$cancer_type <- names(listinput_df)[listinput_df[rownames(listinput_df)==g,]==1]
    }
    ITS_n_r_cancerspecific[[i]] <- tmp
  }else{
    ITS_n_r_cancerspecific[[i]] <- NULL
   }
  }
}

# names(ITS_n_r_genestat) <- ITS_R
# names(ITS_n_r_cancerspecific) <- ITS_R
# names(ITS_n_r_notspecific) <- ITS_R
# save(ITS_n_r_genestat, ITS_n_r_cancerspecific, ITS_n_r_notspecific, 
#      file = "07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITS_n_r_ITSstat.Rdata")


df <- do.call(rbind, ITS_n_r_cancerspecific)

gs_r_specific <- list()
gs_r_nonspecific <- list()

for(i in seq(length(unique(df$gene_name)))) {
    g=unique(df$gene_name)[i]
    tmp = df[df$gene_name == g, ]
    if(length(unique(tmp$cancer_type))>1){
        gs_r_nonspecific[[i]] = tmp
    }else{
        gs_r_specific[[i]] = tmp
    }
}

gs_r_specific <- unique(do.call(rbind, gs_r_specific))
gs_r_nonspecific <- unique(do.call(rbind, gs_r_nonspecific))

# ITS_n_r_notspecific: gene in any ITS in >= 2 cancer type 
# ITS_n_r_cancerspecific: gene in any ITS in == 1 cancer type 
# # gs_nonspecific: ITS (which contain gene) in == 1 cancer type, but gene can show in > 2 ITS in > 2 cancer types
# # gs_specificITS (which contain gene) in == 1 cancer type, but gene can show in >2 ITS in ==2 cancer types

# Venn plot --------------------------------------------------------------------
# type1: ITS_n_r_notspecific: gene in any ITS in >= 2 cancer type 
ITS_n_r_notspecific_df <- do.call(rbind, ITS_n_r_notspecific)

# # type2: gs_nonspecific: ITS (which contain gene) in == 1 cancer type, but gene can show in > 2 ITS in > 2 cancer types
# gs_r_nonspecific
# # type3: gs_specificITS (which contain gene) in == 1 cancer type, but gene can show in >2 ITS in ==2 cancer types
# gs_r_specific
# type2 + type3 =ITS_n_r_cancerspecific

library(VennDiagram)
venn.diagram(x=list(type1_notspecific=unique(ITS_n_r_notspecific_df$gene_name), 
                     type2_ITSspecific_g_not=unique(gs_r_nonspecific$gene_name), 
                     type3_specific=unique(gs_r_specific$gene_name)), 
 filename = "07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITS_overall/ITS_n_r_gene_groups_Venn.png", 
 height = 450, width = 450, resolution =300, imagetype="png", 
 col="white", fill=c(colors()[616], colors()[38], colors()[468]), 
 alpha=c(0.6, 0.6, 0.6), lwd=c(1, 1, 1), cex=0.3, 
 cat.dist=c(-0.07, -0.07, -0.05), cat.pos=c(300, 60, 180), cat.cex=0.3) 


gene_n_r_specific = setdiff(unique(gs_r_specific$gene_name), unique(ITS_n_r_notspecific_df$gene_name))
gene_n_r_nonspecific = union(unique(gs_r_nonspecific$gene_name), unique(ITS_n_r_notspecific_df$gene_name))

ITS_n_r_cancerspecific_df = do.call(rbind, ITS_n_r_cancerspecific)
ITS_n_r_notspecific_df$cancer_type = 'pancancer'
ITS_n_r_genestat = rbind(ITS_n_r_notspecific_df, ITS_n_r_cancerspecific_df[])

ITS_n_r_genestat_specific = ITS_n_r_genestat[is.element(ITS_n_r_genestat$gene_name, gene_n_r_specific),]
ITS_n_r_genestat_nonspecific = ITS_n_r_genestat[is.element(ITS_n_r_genestat$gene_name, gene_n_r_nonspecific),]

length(unique(ITS_n_r_genestat_specific$gene_name))
length(unique(ITS_n_r_genestat_nonspecific$gene_name))
length(intersect(unique(ITS_n_r_genestat_specific$gene_name), unique(ITS_n_r_genestat_nonspecific$gene_name)))

save(ITS_n_r_genestat, ITS_n_r_notspecific, 
     ITS_n_r_cancerspecific, gs_r_specific, gs_r_nonspecific, 
     ITS_n_r_genestat_specific, ITS_n_r_genestat_nonspecific,
     file = "07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITS_n_r_ITSstat.Rdata")

# overall-pie plot -------------------------------------------------------------
library(scales)
blank_theme <- theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )



library(RColorBrewer)
df_pie = data.frame(group = c("type1_notspecific",  "type3_specific"),
                    value = c(length(unique(gene_n_r_nonspecific)), 
                              length(unique(gene_n_r_specific))),
                    col = c(brewer.pal(3, "Set3")[1:2]))
pie <- ggplot(df_pie, aes(x="", y=value, fill=group))+
  geom_bar( width = 1, stat = "identity", alpha = 0.5) + 
  coord_polar("y", direction=1) +
  scale_fill_manual(values=brewer.pal(3, "Set3"))+
  geom_text(aes(y = sum(df_pie$value)-cumsum(df_pie$value)+df_pie$value/2, 
                label = paste(value,"\n", percent(value/sum(value)))), size=5)+
  blank_theme +
  theme(axis.text.x=element_blank())
pie
ggsave("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITS_overall/ITS_n_r_gene_groups.pdf", pie, width = 8, height = 8)


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
library(UpSetR)
library(clusterProfiler)

# source("06_results_summaryandplotting/scripts/plot_for_drug_function.R")
# source("06_results_summaryandplotting/scripts/ITS_analysis/ITS_gene_stat_function.R")
# source("06_results_summaryandplotting/scripts/immunesig_selection_function.R")
# source("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/scripts/drugITGSscore_lincs.R")
source("07_plotting_v2/scripts/functions/ITS_enrich.R")
source("07_plotting_v2/scripts/functions/ITS_function_detection.R")



# ----------------------------------------------------------------------------
dir.create("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/")
dir.create("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITS_sharedgene/")
dir.create("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITS_sharedgene/enrich/")


# LOAD GENES IN LINCS ------------------------------------------------------
lincs_gene = read.csv("/picb/bigdata/project/FengFYM/drug_MoA_reanalysis/LINCS_data/GSE70138/file/GSE92742_Broad_LINCS_gene_info.txt", sep = "\t")

# ------------------------------------------------------------------------------
load("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITS_set.Rdata")

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

# ITS_S = ITS_S[ITS_S != "merck18_tide"]

# ------------------------------------------------------------------------------
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

names(ITS_p_s_genestat) <- ITS_S
names(ITS_p_s_cancerspecific) <- ITS_S
names(ITS_p_s_notspecific) <- ITS_S

# ------------------------------------------------------------------------------

ITS_p_s_notspecific_df <- do.call(rbind, ITS_p_s_notspecific)
ITS_p_s_notspecific_df2 <- unique(data.frame(table(ITS_p_s_notspecific_df$gene_name)))
names(ITS_p_s_notspecific_df2) <- c("gene_name", "num_ITS")
ITS_p_s_notspecific_df2 <- ITS_p_s_notspecific_df2[order(ITS_p_s_notspecific_df2$num_ITS, decreasing = T),]
ITS_p_s_notspecific_df2$rate <- ITS_p_s_notspecific_df2$num_ITS / length(ITS_S)
ITS_p_s_notspecific_df2$weighted_numITS <- NA
for(g in ITS_p_s_notspecific_df2$gene_name){
  ITS_p_s_notspecific_df2[ITS_p_s_notspecific_df2$gene_name==g,]$weighted_numITS <- sum(ITS_p_s_notspecific_df[ITS_p_s_notspecific_df$gene_name==g,]$num_Gene / length(ITS_p_set))
}
ITS_p_s_notspecific_df2$weighted_rate <- ITS_p_s_notspecific_df2$weighted_numITS / length(ITS_S)

write.table(ITS_p_s_notspecific_df2, "07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITS_sharedgene/enrich/ITS_p_s_notspecific_df.txt", quote = F, sep = '\t', row.names = F)


# enrich ----------------------------------------------------------------------
load("07_plotting_v2/data/custom_gs.Rdata")

library(org.Hs.eg.db)
library(topGO)
library(clusterProfiler)
library(msigdbr)
m_df <- msigdbr(species = "Homo sapiens")
#print(paste0("number of genes = ", length(x)))

# enrich ----------------------------------------------------------------------                        
genename <- bitr( ITS_p_s_notspecific_df2$gene_name, 
                  fromType = "SYMBOL", 
                  toType = c("ENTREZID", "SYMBOL"), 
                  OrgDb = org.Hs.eg.db)
names(genename) <- c("gene_name", "ENTREZID")
ITS_p_s_notspecific_df3 <- merge(ITS_p_s_notspecific_df2, genename, all.x=T, by = "gene_name")
ITS_p_s_notspecific_df3 <- ITS_p_s_notspecific_df3[order(ITS_p_s_notspecific_df3$weighted_rate, decreasing = T), ]

tmp2 = as.vector(ITS_p_s_notspecific_df3$weighted_rate)
names(tmp2) = ITS_p_s_notspecific_df3$gene_name

# enrich custom-ORA ---------------------------------------------------------------
savepath = "07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITS_sharedgene/enrich/"
dir.create(savepath)

res = ITS_GS_enrich(ITS_p_s_notspecific_df2$gene_name,
                    gs = gs_all,
                    pvalueCutoff = 1,
                    type = "custom",
                    method = 'ORA',
                    savepath = paste0(savepath, "custom_ORA_s"))
dotplot(res, showCategory = 10) 

# custom-Fisher etc ------------------------------------------------------------

res2 = ITS_GS_enrich(ITS_p_s_notspecific_df2$gene_name, 
                     gs=gs_all,
                     pvalueCutoff = 1,
                     type = "custom",
                     method = 'custom',
                     savepath = paste0(savepath, "custom_custom_s"))
dotplot(res, showCategory = 10) 

# enrich custom-GSEA -------------------------------------------------------------------
tmp2 = as.vector(ITS_p_s_notspecific_df3$weighted_rate)
names(tmp2) = ITS_p_s_notspecific_df3$gene_name
# x = tmp2
res3 = ITS_GS_enrich(tmp2,
                     gs=gs_all,
                     pvalueCutoff = 1,
                     type = "custom",
                     method = 'GSEA',
                     savepath = paste0(savepath, "custom_GSEA_s"))
dotplot(res, showCategory = 10) 

# ----------------------------------------------------------------------------
# functional enrich -----------------------------------------------------------------------
# ----------------------------------------------------------------------------
# enrich hallmark-ORA ---------------------------------------------------------------
res = ITS_GS_enrich(ITS_p_s_notspecific_df2$gene_name,
                    pvalueCutoff = 1,
                    type = "H",
                    method = 'ORA',
                    savepath = paste0(savepath, "hallmark_ORA_s"))
dotplot(res, showCategory = 10) 
cnetplot(res, showCategory = 5) 
cnetplot(res, categorySize="pvalue", foldChange=ITS_p_s_notspecific_df2$gene_name)

# enrich KEGG-ORA ---------------------------------------------------------------
res = ITS_GS_enrich(ITS_p_s_notspecific_df2$gene_name,
                    pvalueCutoff = 1,
                    type = "KEGG",
                    method = 'ORA',
                    savepath = paste0(savepath, "KEGG_ORA_s"))
dotplot(res, showCategory = 10) 

# enrich GOBP-ORA ---------------------------------------------------------------
res = ITS_GS_enrich(ITS_p_s_notspecific_df2$gene_name,
                    pvalueCutoff = 1,
                    type = "GOBP",
                    method = 'ORA',
                    savepath = paste0(savepath, "GOBP_ORA_s"))
dotplot(res, showCategory = 10) 
plotGOgraph(res)
cnetplot(res, showCategory = 5) 
cnetplot(res, categorySize="pvalue", foldChange=ITS_p_s_notspecific_df2$gene_name)


# KEGG GSEA --------------------------------------------------------------------
tmp2 = as.vector(ITS_p_s_notspecific_df3$weighted_rate)
names(tmp2) = ITS_p_s_notspecific_df3$ENTREZID
# x = tmp2
res3 = ITS_GS_enrich(tmp2,
                     pvalueCutoff = 1,
                     type = "KEGG",
                     method = 'GSEA',
                     savepath = paste0(savepath, "KEGG_GSEA_s"))
dotplot(res3, showCategory = 10) 
# hallmark GSEA --------------------------------------------------------------------
res3 = ITS_GS_enrich(tmp2,
                     pvalueCutoff = 1,
                     type = "H",
                     method = 'GSEA',
                     savepath = paste0(savepath, "hallmark_GSEA_s"))
dotplot(res3, showCategory = 10) 

# GOBP GSEA --------------------------------------------------------------------
res3 = ITS_GS_enrich(tmp2,
                     pvalueCutoff = 1,
                     type = "GOBP",
                     method = 'GSEA',
                     savepath = paste0(savepath, "GOBP_GSEA_s"))
dotplot(res3, showCategory = 10) 


# ------------------------------------------------------------------------------
# immune_sig group by prior knowledge
# ------------------------------------------------------------------------------
load("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSsimilarity/ITS_gs_similarity.Rdata")

# ITS_p_set_jcindex
# ITS_p_set_Ochiai

immunesigCL <- read.csv("03_drug_immuneSig_enrichment/data/ITS_hclust_Otsuka_Ochiai_TUMERIC_spearman_r_0.4_pan_metaP_0.05_genevote_0.5/SKCM_immunesig_ITS_index.txt", sep = '\t')
# immunesigCL <- immunesigCL[immunesigCL$Type == "TMEcellSig",]
immunesigCL <- immunesigCL[c("immune_sig", "Index")]
names(immunesigCL) <- c("immune_sig", "Index_CL")

# immunesigCL[grep("CD8+ T cell", immunesigCL$Index_CL),]$immune_sig
res = list()
res_func = list()
for(i in seq(length(unique(immunesigCL$Index_CL)[-33]))){
  # kword = "T cell"
  kword=unique(immunesigCL$Index_CL)[i]
  tmp <- ITS_p_s_notspecific_df[intersect(colnames(ITS_p_s_notspecific_df), immunesigCL[grep(kword, immunesigCL$Index_CL),]$immune_sig)]
  if(ncol(tmp)>0){
    
    tmpgene = names(rowSums(tmp)[rowSums(tmp) >= 0.8*ncol(tmp)])
    length(tmpgene)
    res_func[[i]] =ITSenrichment(ITS = tmpgene,
                                 method = "H")
    # res_func[[i]] = tmpgene_enrich[tmpgene_enrich$p.adjust < 0.05,]$ID
    
    res[[i]] <- ITS_GSenrichTest_res(x = tmpgene)
    res[[i]]$Index_CL = kword
  }
}

names(res_func) <- unique(immunesigCL$Index_CL)[-33]
save(res_func,resall,file="07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSsimilarity/temp_ITS_IndexCL_res.Rdata" )

resall <- do.call(rbind, res)
resall2 <- resall[resall$FishTestP < 0.05, ]
resall2[resall2$Index_CL == "Macrophage",]
res_func['Macrophage']

# # ------------------------------------------------------------------------------
# # cancer specific
# # ------------------------------------------------------------------------------

# names(ITS_p_s_cancerspecific) <- ITS_S
# ITS_p_s_cancerspecific_df <- fromList(ITS_p_s_cancerspecific)
# rownames(ITS_p_s_cancerspecific_df) <- unique(unlist(ITS_p_s_cancerspecific))
# head(rowSums(ITS_p_s_cancerspecific_df))

# ITS_p_s_cancerspecific_gs <- rownames(ITS_p_s_cancerspecific_df)
# ITS_p_s_cancerspecific_gs <- rownames(ITS_p_s_cancerspecific_df[which(rowSums(ITS_p_s_cancerspecific_df)>0.1*ncol(ITS_p_s_cancerspecific_df)),])

# ITS_GSenrichTest_res(x = ITS_p_s_cancerspecific_gs)

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

names(ITS_p_r_genestat) <- ITS_R
names(ITS_p_r_cancerspecific) <- ITS_R
names(ITS_p_r_notspecific) <- ITS_R

# ------------------------------------------------------------------------------

ITS_p_r_notspecific_df <- do.call(rbind, ITS_p_r_notspecific)
ITS_p_r_notspecific_df2 <- unique(data.frame(table(ITS_p_r_notspecific_df$gene_name)))
names(ITS_p_r_notspecific_df2) <- c("gene_name", "num_ITS")
ITS_p_r_notspecific_df2 <- ITS_p_r_notspecific_df2[order(ITS_p_r_notspecific_df2$num_ITS, decreasing = T),]
ITS_p_r_notspecific_df2$rate <- ITS_p_r_notspecific_df2$num_ITS / length(ITS_R)
ITS_p_r_notspecific_df2$weighted_numITS <- NA
for(g in ITS_p_r_notspecific_df2$gene_name){
  ITS_p_r_notspecific_df2[ITS_p_r_notspecific_df2$gene_name==g,]$weighted_numITS <- sum(ITS_p_r_notspecific_df[ITS_p_r_notspecific_df$gene_name==g,]$num_Gene / length(ITS_p_set))
}
ITS_p_r_notspecific_df2$weighted_rate <- ITS_p_r_notspecific_df2$weighted_numITS / length(ITS_R)

write.table(ITS_p_r_notspecific_df2, "07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITS_sharedgene/enrich/ITS_p_r_notspecific_df.txt", quote = F, sep = '\t', row.names = F)


# enrich ----------------------------------------------------------------------
load("07_plotting_v2/data/custom_gs.Rdata")

library(org.Hs.eg.db)
library(topGO)
library(clusterProfiler)
library(msigdbr)
m_df <- msigdbr(species = "Homo sapiens")
#print(paste0("number of genes = ", length(x)))

# enrich ----------------------------------------------------------------------                        
genename <- bitr( ITS_p_r_notspecific_df2$gene_name, 
                  fromType = "SYMBOL", 
                  toType = c("ENTREZID", "SYMBOL"), 
                  OrgDb = org.Hs.eg.db)
names(genename) <- c("gene_name", "ENTREZID")
ITS_p_r_notspecific_df3 <- merge(ITS_p_r_notspecific_df2, genename, all.x=T, by = "gene_name")
ITS_p_r_notspecific_df3 <- ITS_p_r_notspecific_df3[order(ITS_p_r_notspecific_df3$weighted_rate, decreasing = T), ]

tmp2 = as.vector(ITS_p_r_notspecific_df3$weighted_rate)
names(tmp2) = ITS_p_r_notspecific_df3$gene_name

# enrich custom-ORA ---------------------------------------------------------------
savepath = "07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITS_sharedgene/enrich/"
dir.create(savepath)

# custom-Fisher etc ------------------------------------------------------------

res2 = ITS_GS_enrich(ITS_p_r_notspecific_df2$gene_name, 
                     gs=gs_all,
                     pvalueCutoff = 1,
                     type = "custom",
                     method = 'custom',
                     savepath = paste0(savepath, "custom_custom_r"))


# custom-Fisher etc ------------------------------------------------------------

res2 = ITS_GS_enrich(ITS_p_r_notspecific_df2$gene_name, 
                     gs=gs_all,
                     pvalueCutoff = 1,
                     type = "custom",
                     method = 'custom',
                     savepath = paste0(savepath, "custom_custom_r"))



res = ITS_GS_enrich(ITS_p_r_notspecific_df2$gene_name,
                    gs = gs_all,
                    pvalueCutoff = 1,
                    type = "custom",
                    method = 'ORA',
                    savepath = paste0(savepath, "custom_ORA_r"))
dotplot(res, showCategory = 10) 

# enrich custom-GSEA -------------------------------------------------------------------
tmp2 = as.vector(ITS_p_r_notspecific_df3$weighted_rate)
names(tmp2) = ITS_p_r_notspecific_df3$gene_name
# x = tmp2
res3 = ITS_GS_enrich(tmp2,
                     gs=gs_all,
                     pvalueCutoff = 1,
                     type = "custom",
                     method = 'GSEA',
                     savepath = paste0(savepath, "custom_GSEA_r"))
dotplot(res3, showCategory = 10) 

# ----------------------------------------------------------------------------
# functional enrich -----------------------------------------------------------------------
# ----------------------------------------------------------------------------
# enrich hallmark-ORA ---------------------------------------------------------------
res = ITS_GS_enrich(ITS_p_r_notspecific_df2$gene_name,
                    pvalueCutoff = 1,
                    type = "H",
                    method = 'ORA',
                    savepath = paste0(savepath, "hallmark_ORA_r"))
dotplot(res, showCategory = 10) 
cnetplot(res, showCategory = 5) 
cnetplot(res, categorySize="pvalue", foldChange=ITS_p_r_notspecific_df2$gene_name)

# enrich KEGG-ORA ---------------------------------------------------------------
res = ITS_GS_enrich(ITS_p_r_notspecific_df2$gene_name,
                    pvalueCutoff = 1,
                    type = "KEGG",
                    method = 'ORA',
                    savepath = paste0(savepath, "KEGG_ORA_r"))
dotplot(res, showCategory = 10) 

# enrich GOBP-ORA ---------------------------------------------------------------
res = ITS_GS_enrich(ITS_p_r_notspecific_df2$gene_name,
                    pvalueCutoff = 1,
                    type = "GOBP",
                    method = 'ORA',
                    savepath = paste0(savepath, "GOBP_ORA_r"))
dotplot(res, showCategory = 10) 
plotGOgraph(res)
cnetplot(res, showCategory = 5) 
cnetplot(res, categorySize="pvalue", foldChange=ITS_p_r_notspecific_df2$gene_name)


# KEGG GSEA --------------------------------------------------------------------
tmp2 = as.vector(ITS_p_r_notspecific_df3$weighted_rate)
names(tmp2) = ITS_p_r_notspecific_df3$ENTREZID
# x = tmp2
res3 = ITS_GS_enrich(tmp2,
                     pvalueCutoff = 1,
                     type = "KEGG",
                     method = 'GSEA',
                     savepath = paste0(savepath, "KEGG_GSEA_r"))
dotplot(res3, showCategory = 10) 
# hallmark GSEA --------------------------------------------------------------------
res3 = ITS_GS_enrich(tmp2,
                     pvalueCutoff = 1,
                     type = "H",
                     method = 'GSEA',
                     savepath = paste0(savepath, "hallmark_GSEA_r"))
dotplot(res3, showCategory = 10) 

# GOBP GSEA --------------------------------------------------------------------
res3 = ITS_GS_enrich(tmp2,
                     pvalueCutoff = 1,
                     type = "GOBP",
                     method = 'GSEA',
                     savepath = paste0(savepath, "GOBP_GSEA_r"))
dotplot(res3, showCategory = 10) 







# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ITS_n



# sensitive 
ITS_n_s_genestat <- list()
ITS_n_s_cancerspecific = list()
ITS_n_s_notspecific = list()

for(i in 1:length(ITS_S)){
  # i = 1
  ITSname <- ITS_S[i]
  listinput <- lapply(ITS_n_set, function(x){x[[ITSname]]})
  genelist = unique(unlist(listinput))
  listinput_df <- fromList(listinput)
  rownames(listinput_df) <- genelist
  table(apply(listinput_df,1,sum))
  apply(listinput_df,2,sum)
  
  ITS_n_s_genestat[[i]] <- listinput_df
  
  ITS_n_s_cancerspecific[[i]] =  unique(names(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)==1)]))
  ITS_n_s_notspecific[[i]] = unique(names(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)>=2)]))
  
}


# names(ITS_n_s_cancerspecific) <- ITS_S
# ITS_n_s_cancerspecific_df <- fromList(ITS_n_s_cancerspecific)
# rownames(ITS_n_s_cancerspecific_df) <- unique(unlist(ITS_n_s_cancerspecific))
# head(rowSums(ITS_n_s_cancerspecific_df))
# ITS_n_s_cancerspecific_gs <- rownames(ITS_n_s_cancerspecific_df)
# ITS_n_s_cancerspecific_gs <- rownames(ITS_n_s_cancerspecific_df[which(rowSums(ITS_n_s_cancerspecific_df)>0.3*ncol(ITS_n_s_cancerspecific_df)),])
# ITS_GSenrichTest_res(x = ITS_n_s_cancerspecific_gs)

# ------------------------------------------------------------------------------

names(ITS_n_s_notspecific) <- ITS_S
ITS_n_s_notspecific_df <- fromList(ITS_n_s_notspecific)
rownames(ITS_n_s_notspecific_df) <- unique(unlist(ITS_n_s_notspecific))
head(rowSums(ITS_n_s_notspecific_df))

ITS_n_s_notspecific_gs <- rownames(ITS_n_s_notspecific_df)
ITS_n_s_notspecific_gs <- rownames(ITS_n_s_notspecific_df[which(rowSums(ITS_n_s_notspecific_df)>0.1*ncol(ITS_n_s_notspecific_df)),])


ITS_GSenrichTest_res(x = ITS_n_s_notspecific_gs)

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# resistant 
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

ITS_n_r_genestat = list()
ITS_n_r_cancerspecific = list()
ITS_n_r_notspecific = list()
for(i in 1:length(ITS_R)){
  # i = 1
  ITSname <- ITS_R[i]
  listinput <- lapply(ITS_n_set, function(x){x[[ITSname]]})
  genelist = unique(unlist(listinput))
  listinput_df <- fromList(listinput)
  rownames(listinput_df) <- genelist
  table(apply(listinput_df,1,sum))
  apply(listinput_df,2,sum)
  
  ITS_n_r_genestat[[i]] <- listinput_df
  
  ITS_n_r_cancerspecific[[i]] =  unique(names(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)==1)]))
  ITS_n_r_notspecific[[i]] = unique(names(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)>=2)]))
  
}

# names(ITS_n_r_cancerspecific) <- ITS_R
# ITS_n_r_cancerspecific_df <- fromList(ITS_n_r_cancerspecific)
# rownames(ITS_n_r_cancerspecific_df) <- unique(unlist(ITS_n_r_cancerspecific))
# head(rowSums(ITS_n_r_cancerspecific_df))
# ITS_n_r_cancerspecific_gs <- rownames(ITS_n_r_cancerspecific_df)
# ITS_n_r_cancerspecific_gs <- rownames(ITS_n_r_cancerspecific_df[which(rowSums(ITS_n_r_cancerspecific_df)>3),])
# ITS_GSenrichTest_res(x = ITS_n_r_cancerspecific_gs)

# ------------------------------------------------------------------------------

names(ITS_n_r_notspecific) <- ITS_R
ITS_n_r_notspecific_df <- fromList(ITS_n_r_notspecific)
rownames(ITS_n_r_notspecific_df) <- unique(unlist(ITS_n_r_notspecific))
head(rowSums(ITS_n_r_notspecific_df))

ITS_n_r_notspecific_gs <- rownames(ITS_n_r_notspecific_df)
ITS_n_r_notspecific_gs <- rownames(ITS_n_r_notspecific_df[which(rowSums(ITS_n_r_notspecific_df)>3),])


ITS_GSenrichTest_res(x = ITS_n_r_notspecific_gs)

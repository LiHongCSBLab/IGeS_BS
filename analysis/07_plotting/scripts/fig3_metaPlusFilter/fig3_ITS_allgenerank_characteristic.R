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
source("07_plotting_v2/scripts/functions/ITS_g_cancer_plot_function.R")



# ----------------------------------------------------------------------------
dir.create("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/")
dir.create("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/")

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
# ITS positive
# ------------------------------------------------------------------------------
# stat: sensitive-related ITS ----------------------------
# ------------------------------------------------------------------------------
# sensitive 
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


ITS_p_s_df <- do.call(rbind, ITS_p_s)
ITS_p_s_df2 <- unique(data.frame(table(ITS_p_s_df$gene_name)))
names(ITS_p_s_df2) <- c("gene_name", "num_ITS")
ITS_p_s_df2 <- ITS_p_s_df2[order(ITS_p_s_df2$num_ITS, decreasing = T),]
ITS_p_s_df2$rate <- ITS_p_s_df2$num_ITS / length(ITS_S)
ITS_p_s_df2$weighted_numITS <- NA
for(g in ITS_p_s_df2$gene_name){
  ITS_p_s_df2[ITS_p_s_df2$gene_name==g,]$weighted_numITS <- sum(ITS_p_s_df[ITS_p_s_df$gene_name==g,]$num_Gene / length(ITS_p_set))
}

ITS_p_s_df2$weighted_rate <- ITS_p_s_df2$weighted_numITS / length(ITS_S)
ITS_p_s_df2 <- ITS_p_s_df2[order(ITS_p_s_df2$weighted_rate, decreasing = T), ]
# ITS_p_s_df2


# enrich ----------------------------------------------------------------------
load("07_plotting_v2/data/custom_gs.Rdata")
load( "07_plotting_v2/data/drugTarget_merged3DB.Rdata")


library(org.Hs.eg.db)
library(topGO)
library(clusterProfiler)
library(msigdbr)
m_df <- msigdbr(species = "Homo sapiens")
#print(paste0("number of genes = ", length(x)))

# enrich ----------------------------------------------------------------------                        
genename <- bitr( ITS_p_s_df2$gene_name, 
                  fromType = "SYMBOL", 
                  toType = c("ENTREZID", "SYMBOL"), 
                  OrgDb = org.Hs.eg.db)
names(genename) <- c("gene_name", "ENTREZID")
ITS_p_s_df3 <- merge(ITS_p_s_df2, genename, all.x=T, by = "gene_name")
ITS_p_s_df3 <- ITS_p_s_df3[order(ITS_p_s_df3$weighted_rate, decreasing = T), ]

write.table(ITS_p_s_df3, "07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/ITS_p_s_df.txt", quote = F, sep = '\t', row.names = F)


# enrich custom-GSEA -------------------------------------------------------------------
# tmp2 = as.vector(ITS_p_s_df3$weighted_rate)
# ITS_p_s_df3 <- read.table("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/ITS_p_s_df.txt", sep = '\t', header = T)

tmp2 = as.vector(scale(ITS_p_s_df3$weighted_rate))
names(tmp2) = ITS_p_s_df3$gene_name

savepath = "07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/enrich/"
dir.create(savepath)

res = ITS_GS_enrich(tmp2,
                    gs=gs_all,
                    pvalueCutoff = 1,
                    type = "custom",
                    method = 'GSEA',
                    savepath = paste0(savepath, "custom_GSEA_s_scaled"))

gseaplot2(res, geneSetID =c(6,8))

data.frame(res)[c(1, 6,7)]
res2=data.frame(res)
res2[res2$ID == "ITSp_vs_OG_oncoKB",]

# plot - heatmap



dir.create("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/function_enrich/")

for(gsname in unique(res2$ID)){ #[11:15]
  # gsname="ITSp_vs_OG_oncoKB"
  # gsname="ITSp_vs_drugtarget"
  # ------------------------------------------------------------------------------
  # 1) functional enrichment
  # ITSp_vs_OG_oncoKB
  dir.create(paste0("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/function_enrich/", gsname, "/"))
  savepath = paste0("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/function_enrich/", gsname, "/")
  
  # ------------------------------------------------------------------------------
  # 0) extract gene list
  glist <- unlist(strsplit(res2[res2$ID == gsname,]$core_enrichment, "/"))
  
  # ------------------------------------------------------------------------------
  # 1) if these genes in original immune signatures
  try({
    length(glist)
    length(unique(immunesigGS_S[is.element(immunesigGS_S$gene_name, glist),]$gene_name))
    setdiff(glist, immunesigGS_S$gene_name)
    library(VennDiagram)
    venn.diagram(x=list(ITSgene=glist, 
                        originalSignaturesGenes=unique(immunesigGS_S$gene_name)), 
                 filename = paste0(savepath, "gene_S_Venn.png"), 
                 height = 450, width = 450, resolution =300, imagetype="png", 
                 col="white", fill=c(colors()[616], colors()[38]), 
                 alpha=c(0.6, 0.6), lwd=c(1, 1), cex=0.3, 
                 cat.dist=c(-0.07, -0.03), cat.pos=c(300, 60), cat.cex=0.3) 
  }, silent = T)
  
  # ------------------------------------------------------------------------------
  # 2) correlation between gene expression and original immune signatures 
  # paste0(07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITS_correlation_sgdotplot/", gsname)
  
  
  # ------------------------------------------------------------------------------
  # 3) functional enrichment
  
  # enrich hallmark-ORA ---------------------------------------------------------------
  try({
    
    enrichHallmark_res = ITS_GS_enrich(glist,
                                       pvalueCutoff = 0.1,
                                       type = "H",
                                       method = 'ORA',
                                       savepath = paste0(savepath,"hallmark_ORA_s"))
    enrichHallmark_res2 = data.frame(enrichHallmark_res)
    # unlist(strsplit(enrichHallmark_res2[enrichHallmark_res2$ID == "HALLMARK_KRAS_SIGNALING_UP",]$geneID, "/"))
    
    
    # if(nrow(enrichKEGG_res2)>10){
    #   pdf(paste0(savepath,"hallmark_ORA_s.pdf"), onefile = TRUE)
    #   print(dotplot(enrichHallmark_res, showCategory = 10))
    #   print(cnetplot(enrichHallmark_res, showCategory = 10))
    #   dev.off()
    # }else{
      pdf(paste0(savepath,"hallmark_ORA_s.pdf"), onefile = TRUE)
      print(dotplot(enrichHallmark_res, showCategory = nrow(enrichHallmark_res2)))
      print(cnetplot(enrichHallmark_res, showCategory = nrow(enrichHallmark_res2)))
      dev.off()
    # }
  }, silent = T)
  
  # enrich KEGG-ORA ---------------------------------------------------------------
  try({
    enrichKEGG_res = ITS_GS_enrich(glist,
                                   pvalueCutoff = 0.05,
                                   type = "KEGG",
                                   method = 'ORA',
                                   savepath = paste0(savepath,"KEGG_ORA_s"))
    
    enrichKEGG_res2 = data.frame(enrichKEGG_res)
    # unlist(strsplit(enrichKEGG_res2[enrichKEGG_res2$Description == "NF-kappa B signaling pathway",]$geneID, "/"))
    # unlist(strsplit(enrichKEGG_res2[enrichKEGG_res2$Description == "VEGF signaling pathway",]$geneID, "/"))
    
    # if(nrow(enrichKEGG_res2)>10){
    #   pdf(paste0(savepath,"KEGG_ORA_s.pdf"), onefile = TRUE)
    #   print(dotplot(enrichKEGG_res, showCategory = 10))
    #   print(cnetplot(enrichKEGG_res, showCategory = 10))
    #   dev.off()
    # }else{
      pdf(paste0(savepath,"KEGG_ORA_s.pdf"), onefile = TRUE)
      print(dotplot(enrichKEGG_res, showCategory = nrow(enrichKEGG_res2)) )
      print(cnetplot(enrichKEGG_res, showCategory = nrow(enrichKEGG_res2)) )
      dev.off()
    # }
  }, silent = T)
  
  # enrich GOBP-ORA ---------------------------------------------------------------
  try({
    enrichGOBP_res = ITS_GS_enrich(glist,
                                   pvalueCutoff = 0.05,
                                   type = "GOBP",
                                   method = 'ORA',
                                   savepath = paste0(savepath,"GOBP_ORA_s"))
    enrichGOBP_res2 = data.frame(enrichGOBP_res)
    
    # if(nrow(enrichGOBP_res2)>10){
    #   pdf(paste0(savepath,"GOBP_ORA_s.pdf"), onefile = TRUE)
    #   print(dotplot(enrichGOBP_res, showCategory = 10))
    #   print(cnetplot(enrichGOBP_res, showCategory = 10)) 
    #   dev.off()
      
    # }else{
      pdf(paste0(savepath,"GOBP_ORA_s.pdf"), onefile = TRUE)
      print(dotplot(enrichGOBP_res, showCategory = nrow(enrichGOBP_res2)))
      print(cnetplot(enrichGOBP_res, showCategory = nrow(enrichGOBP_res2))) 
      dev.off()
    # }
  }, silent = T)
  
  # enrich Reactome-ORA ---------------------------------------------------------------
  try({
    enrichReactome_res = ITS_GS_enrich(glist,
                                   pvalueCutoff = 0.05,
                                   type = "Reactome",
                                   method = 'ORA',
                                   savepath = paste0(savepath,"Reactome_ORA_s"))
    enrichReactome_res2 = data.frame(enrichReactome_res)
    
    # if(nrow(enrichReactome_res2)>10){
    #   pdf(paste0(savepath,"GOBP_ORA_s.pdf"), onefile = TRUE)
    #   print(dotplot(enrichGOBP_res, showCategory = 10))
    #   print(cnetplot(enrichGOBP_res, showCategory = 10)) 
    #   dev.off()
      
    # }else{
      pdf(paste0(savepath,"Reactome_ORA_s.pdf"), onefile = TRUE)
      print(dotplot(enrichReactome_res, showCategory = nrow(enrichReactome_res2)))
      print(cnetplot(enrichReactome_res, showCategory = nrow(enrichReactome_res2))) 
      dev.off()
    # }
  }, silent = T)
  
  # ------------------------------------------------------------------------------
  # 4) 这些基因在不同癌症中的ITS数目统计（饼图热图）.
  #    注释基因功能
  # source("07_plotting_v2/scripts/functions/ITS_g_cancer_plot_function.R")
  try({

  ITS_g_cancer_plot(ITS_set = ITS_p_set, 
                    genelist = glist,
                    type = 's',
                    enrichRes = enrichKEGG_res2,
                    savepath = savepath,
                    filename = "gene_cancer_ITSrate_stat_s")
  }, silent = T)

}


# ------------------------------------------------------------------------------
# Drug target analysis
# ------------------------------------------------------------------------------
load( "07_plotting_v2/data/drugTarget_merged3DB.Rdata")
# ------------------------------------------------------------------------------

dir.create("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/drug_merge3DB")
dir.create("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/drug_merge3DB/GSEAenrich/")
dir.create("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSn_sharedgene/")
dir.create("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSn_sharedgene/drug_merge3DB")
dir.create("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSn_sharedgene/drug_merge3DB/GSEAenrich/")

# 1) drugTarget_drugbank-Drug -------------------------------------------------------------------
filepath <- ("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/drug_merge3DB/GSEAenrich/ITSp_s_vs_drugtarget")

glist <- unlist(strsplit(res2[res2$ID == "ITSp_vs_drugtarget",]$core_enrichment, "/"))

# DT_druglist <- unique(drugTarget_all[is.element(drugTarget_all$genename, glist),]$drugname_lower)
# DT_drug <- drugTarget_all[is.element(drugTarget_all$drugname_lower, DT_druglist),]
DT_drug <- unique(drugTarget_all[is.element(drugTarget_all$genename, glist),])

write.table(DT_drug, paste0(filepath, ".txt"), sep = '\t', quote = F, row.names = F)

# DT_druglistFDA <- unique(drugTarget_FDA[is.element(drugTarget_FDA$genename, glist),]$drugname_lower)
# DT_drugFDA<- drugTarget_FDA[is.element(drugTarget_FDA$drugname_lower, DT_druglist),]
DT_drugFDA <- unique(drugTarget_FDA[is.element(drugTarget_FDA$genename, glist),])
write.table(DT_drugFDA, paste0(filepath, "_FDA.txt"), sep = '\t', quote = F, row.names = F)

DT_drug2 = unique(DT_drug[c('drugname_lower', 'genename')])
DT_drugFDA2 = unique(DT_drugFDA[c('drugname_lower', 'genename', 'ifCancer_drugbank')])

drugstat <- data.frame(group = c("num_drug",  "num_drugFDA", "num_drugFDAcancer"),
                       value = c(length(unique(DT_drug2$drugname_lower)),
                                 length(unique(DT_drugFDA2$drugname_lower)), 
                                 length(unique(DT_drugFDA2[which(DT_drugFDA2$ifCancer_drugbank == 'yes'), ]$drugname_lower))))

write.table(drugstat, paste0(filepath, "_stat.txt"), sep = '\t', quote = F, row.names = F)

# unique(DT_drugFDAcancer$name)
# sort(table(DT_drugFDAcancer$genename))


# 2) OG-Drug -------------------------------------------------------------------
filepath <- ("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/drug_merge3DB/GSEAenrich/ITSp_s_vs_OG_oncoKB")

glist <- unlist(strsplit(res2[res2$ID == "ITSp_vs_OG_oncoKB",]$core_enrichment, "/"))

# OG_druglist <- unique(drugTarget_all[is.element(drugTarget_all$genename, glist),]$drugname_lower)
# OG_drug <- drugTarget_all[is.element(drugTarget_all$drugname_lower, OG_druglist),]
OG_drug <- unique(drugTarget_all[is.element(drugTarget_all$genename, glist),])

write.table(OG_drug, paste0(filepath, ".txt"), sep = '\t', quote = F, row.names = F)

OG_druglistFDA <- unique(drugTarget_FDA[is.element(drugTarget_FDA$genename, glist),]$drugname_lower)
# OG_drugFDA<- drugTarget_FDA[is.element(drugTarget_FDA$drugname_lower, OG_druglist),]
OG_drugFDA<- unique(drugTarget_FDA[is.element(drugTarget_FDA$genename, glist),])


write.table(OG_drugFDA, paste0(filepath, "_FDA.txt"), sep = '\t', quote = F, row.names = F)

OG_drug2 = unique(OG_drug[c('drugname_lower', 'genename')])
OG_drugFDA2 = unique(OG_drugFDA[c('drugname_lower', 'genename', 'ifCancer_drugbank')])

drugstat <- data.frame(group = c("num_drug",  "num_drugFDA", "num_drugFDAcancer"),
                       value = c(length(unique(OG_drug2$drugname_lower)),
                                 length(unique(OG_drugFDA2$drugname_lower)), 
                                 length(unique(OG_drugFDA2[which(OG_drugFDA2$ifCancer_drugbank == 'yes'), ]$drugname_lower))))

write.table(drugstat, paste0(filepath, "_stat.txt"), sep = '\t', quote = F, row.names = F)




# 3) TSG-Drug -------------------------------------------------------------------
filepath <- ("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/drug_merge3DB/GSEAenrich/ITSp_s_vs_TSG_oncoKB")

glist <- unlist(strsplit(res2[res2$ID == "ITSp_vs_TSG_oncoKB",]$core_enrichment, "/"))


# TSG_druglist <- unique(drugTarget_all[is.element(drugTarget_all$genename, glist),]$drugname_lower)
# TSG_drug <- drugTarget_all[is.element(drugTarget_all$drugname_lower, TSG_druglist),]
TSG_drug <- unique(drugTarget_all[is.element(drugTarget_all$genename, glist),])


write.table(TSG_drug, paste0(filepath, ".txt"), sep = '\t', quote = F, row.names = F)

# TSG_druglistFDA <- unique(drugTarget_FDA[is.element(drugTarget_FDA$genename, glist),]$drugname_lower)
# TSG_drugFDA<- drugTarget_FDA[is.element(drugTarget_FDA$drugname_lower, TSG_druglist),]
TSG_drugFDA <- unique(drugTarget_FDA[is.element(drugTarget_FDA$genename, glist),])

write.table(TSG_drugFDA, paste0(filepath, "_FDA.txt"), sep = '\t', quote = F, row.names = F)

TSG_drug2 = unique(TSG_drug[c('drugname_lower', 'genename')])
TSG_drugFDA2 = unique(TSG_drugFDA[c('drugname_lower', 'genename', 'ifCancer_drugbank')])

drugstat <- data.frame(group = c("num_drug",  "num_drugFDA", "num_drugFDAcancer"),
                       value = c(length(unique(TSG_drug2$drugname_lower)),
                                 length(unique(TSG_drugFDA2$drugname_lower)), 
                                 length(unique(TSG_drugFDA2[which(TSG_drugFDA2$ifCancer_drugbank == 'yes'), ]$drugname_lower))))

write.table(drugstat, paste0(filepath, "_stat.txt"), sep = '\t', quote = F, row.names = F)





# # 3) TSG-Drug -------------------------------------------------------------------
# filepath <- ("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/drug_merge3DB/ORAenrich/ITSp_s_vs_TSG_oncoKB")

# glist <- unlist(strsplit(res2[res2$ID == "ITSp_vs_TSG_oncoKB",]$geneID, "/"))
# glist <- ITS_p_s_notspecific_df[is.element(ITS_p_s_notspecific_df$ENTREZID, glist), ]$gene_name


# TSG_druglist <- unique(drugTarget_all[is.element(drugTarget_all$genename, glist$gene_name),]$drugname_lower)
# TSG_drug <- drugTarget_all[is.element(drugTarget_all$drugname_lower, TSG_druglist),]

# write.table(TSG_drug, paste0(filepath, ".txt"), sep = '\t', quote = F, row.names = F)

# # TSG_druglistFDA <- unique(drugTarget_FDA[is.element(drugTarget_FDA$genename, glist),]$drugname_lower)
# # TSG_drugFDA<- drugTarget_FDA[is.element(drugTarget_FDA$drugname_lower, TSG_druglist),]
# TSG_drugFDA <- unique(drugTarget_FDA[is.element(drugTarget_FDA$genename, glist),])

# write.table(TSG_drugFDA, paste0(filepath, "_FDA.txt"), sep = '\t', quote = F, row.names = F)

# TSG_drug2 = unique(TSG_drug[c('drugname_lower', 'genename')])
# TSG_drugFDA2 = unique(TSG_drugFDA[c('drugname_lower', 'genename', 'ifCancer_drugbank')])

# drugstat <- data.frame(group = c("num_drug",  "num_drugFDA", "num_drugFDAcancer"),
#                     value = c(length(unique(TSG_drug2$drugname_lower)),
#                               length(unique(TSG_drugFDA2$drugname_lower)), 
#                               length(unique(TSG_drugFDA2[which(TSG_drugFDA2$ifCancer_drugbank == 'yes'), ]$drugname_lower))))

# write.table(drugstat, paste0(filepath, "_stat.txt"), sep = '\t', quote = F, row.names = F)




# ------------------------------------------------------------------------------
# stat: resistant-related ITS --------------------------------------------------
# ------------------------------------------------------------------------------

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


ITS_p_r_df <- do.call(rbind, ITS_p_r)
ITS_p_r_df2 <- unique(data.frame(table(ITS_p_r_df$gene_name)))
names(ITS_p_r_df2) <- c("gene_name", "num_ITS")
ITS_p_r_df2 <- ITS_p_r_df2[order(ITS_p_r_df2$num_ITS, decreasing = T),]
ITS_p_r_df2$rate <- ITS_p_r_df2$num_ITS / length(ITS_R)
ITS_p_r_df2$weighted_numITS <- NA
for(g in ITS_p_r_df2$gene_name){
  ITS_p_r_df2[ITS_p_r_df2$gene_name==g,]$weighted_numITS <- sum(ITS_p_r_df[ITS_p_r_df$gene_name==g,]$num_Gene / length(ITS_p_set))
}

ITS_p_r_df2$weighted_rate <- ITS_p_r_df2$weighted_numITS / length(ITS_R)
ITS_p_r_df2 <- ITS_p_r_df2[order(ITS_p_r_df2$weighted_rate, decreasing = T), ]
# head(ITS_p_r_df2)


# enrich ----------------------------------------------------------------------
load("07_plotting_v2/data/custom_gs.Rdata")


library(org.Hs.eg.db)
library(topGO)
library(clusterProfiler)
library(msigdbr)
m_df <- msigdbr(species = "Homo sapiens")
#print(paste0("number of genes = ", length(x)))

# enrich ----------------------------------------------------------------------                        
genename <- bitr( ITS_p_r_df2$gene_name, 
                  fromType = "SYMBOL", 
                  toType = c("ENTREZID", "SYMBOL"), 
                  OrgDb = org.Hs.eg.db)
names(genename) <- c("gene_name", "ENTREZID")
ITS_p_r_df3 <- merge(ITS_p_r_df2, genename, all.x=T, by = "gene_name")
ITS_p_r_df3 <- ITS_p_r_df3[order(ITS_p_r_df3$weighted_rate, decreasing = T), ]
write.table(ITS_p_r_df3, "07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/ITS_p_r_df.txt", quote = F, sep = '\t', row.names = F)


# enrich custom-GSEA -------------------------------------------------------------------
# tmp2 = as.vector(ITS_p_r_df3$weighted_rate)
tmp2 = as.vector(scale(ITS_p_r_df3$weighted_rate))
names(tmp2) = ITS_p_r_df3$gene_name

savepath = "07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/enrich/"
dir.create(savepath)

res = ITS_GS_enrich(tmp2,
                    gs=gs_all,
                    pvalueCutoff = 1,
                    type = "custom",
                    method = 'GSEA',
                    savepath = paste0(savepath, "custom_GSEA_r_scaled"))
data.frame(res)[c(1,6,7)]
res2=data.frame(res)
res2[res2$ID == "ITSp_vs_OG_oncoKB",]

# plot - heatmap



dir.create("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/function_enrich/")

for(gsname in unique(res2$ID)){
  #   gsname="ITSp_vs_OG_oncoKB"
  # gsname="ITSp_vs_drugtarget"
  # ------------------------------------------------------------------------------
  # 1) functional enrichment
  # ITSp_vs_OG_oncoKB
  dir.create(paste0("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/function_enrich/", gsname, "/"))
  savepath = paste0("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/function_enrich/", gsname, "/")
  
  # ------------------------------------------------------------------------------
  # 0) extract gene list
  glist <- unlist(strsplit(res2[res2$ID == gsname,]$core_enrichment, "/"))
  
  # ------------------------------------------------------------------------------
  # 1) if these genes in original immune signatures
  try({
    length(glist)
    length(unique(immunesigGS_R[is.element(immunesigGS_R$gene_name, glist),]$gene_name))
    setdiff(glist, immunesigGS_R$gene_name)
    library(VennDiagram)
    venn.diagram(x=list(ITSgene=glist, 
                        originalSignaturesGenes=unique(immunesigGS_R$gene_name)), 
                 filename = paste0(savepath, "gene_r_Venn.png"), 
                 height = 450, width = 450, resolution =300, imagetype="png", 
                 col="white", fill=c(colors()[616], colors()[38]), 
                 alpha=c(0.6, 0.6), lwd=c(1, 1), cex=0.3, 
                 cat.dist=c(-0.07, -0.03), cat.pos=c(300, 60), cat.cex=0.3) 
  }, silent = T)
  
  # ------------------------------------------------------------------------------
  # 2) correlation between gene expression and original immune signatures 
  # paste0(07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITS_correlation_sgdotplot/", gsname)
  
  
  # ------------------------------------------------------------------------------
  # 3) functional enrichment
  
  # enrich hallmark-ORA ---------------------------------------------------------------
  try({
    
    enrichHallmark_res = ITS_GS_enrich(glist,
                                       pvalueCutoff = 0.1,
                                       type = "H",
                                       method = 'ORA',
                                       savepath = paste0(savepath,"hallmark_ORA_r"))
    enrichHallmark_res2 = data.frame(enrichHallmark_res)
    enrichKEGG_res2 = data.frame(enrichKEGG_res)
    
    # if(nrow(enrichKEGG_res2)>10){
    #   pdf(paste0(savepath,"hallmark_ORA_r.pdf"), onefile = TRUE)
    #   print(dotplot(enrichHallmark_res, showCategory = 10))
    #   print(cnetplot(enrichHallmark_res, showCategory = 10))
    #   dev.off()
    # }else{
      pdf(paste0(savepath,"hallmark_ORA_r.pdf"), onefile = TRUE)
      print(dotplot(enrichHallmark_res, showCategory = nrow(enrichHallmark_res2)))
      print(cnetplot(enrichHallmark_res, showCategory = nrow(enrichHallmark_res2)))
      dev.off()
    # }
  }, silent = T)
  
  # enrich KEGG-ORA ---------------------------------------------------------------
  try({
    enrichKEGG_res = ITS_GS_enrich(glist,
                                   pvalueCutoff = 0.05,
                                   type = "KEGG",
                                   method = 'ORA',
                                   savepath = paste0(savepath,"KEGG_ORA_r"))
    
    enrichKEGG_res2 = data.frame(enrichKEGG_res)
    # unlist(strsplit(enrichKEGG_res2[enrichKEGG_res2$Description == "NF-kappa B signaling pathway",]$geneID, "/"))
    # unlist(strsplit(enrichKEGG_res2[enrichKEGG_res2$Description == "VEGF signaling pathway",]$geneID, "/"))
    
    # if(nrow(enrichKEGG_res2)>10){
    #   pdf(paste0(savepath,"KEGG_ORA_r.pdf"), onefile = TRUE)
    #   print(dotplot(enrichKEGG_res, showCategory = 10))
    #   print(cnetplot(enrichKEGG_res, showCategory = 10))
    #   dev.off()
    # }else{
      pdf(paste0(savepath,"KEGG_ORA_r.pdf"), onefile = TRUE)
      print(dotplot(enrichKEGG_res, showCategory = nrow(enrichKEGG_res2)) )
      print(cnetplot(enrichKEGG_res, showCategory = nrow(enrichKEGG_res2)) )
      dev.off()
    # }
  }, silent = T)
  
  # enrich GOBP-ORA ---------------------------------------------------------------
  try({
    enrichGOBP_res = ITS_GS_enrich(glist,
                                   pvalueCutoff = 0.05,
                                   type = "GOBP",
                                   method = 'ORA',
                                   savepath = paste0(savepath,"GOBP_ORA_r"))
    enrichGOBP_res2 = data.frame(enrichGOBP_res)
    
    # if(nrow(enrichGOBP_res2)>10){
    #   pdf(paste0(savepath,"GOBP_ORA_r.pdf"), onefile = TRUE)
    #   print(dotplot(enrichGOBP_res, showCategory = 10))
    #   print(cnetplot(enrichGOBP_res, showCategory = 10)) 
    #   dev.off()
      
    # }else{
      pdf(paste0(savepath,"GOBP_ORA_r.pdf"), onefile = TRUE)
      print(dotplot(enrichGOBP_res, showCategory = nrow(enrichGOBP_res2)))
      print(cnetplot(enrichGOBP_res, showCategory = nrow(enrichGOBP_res2))) 
      dev.off()
    # }
  }, silent = T)
  
  
  # enrich Reactome-ORA ---------------------------------------------------------------
  try({
    enrichReactome_res = ITS_GS_enrich(glist,
                                   pvalueCutoff = 0.05,
                                   type = "Reactome",
                                   method = 'ORA',
                                   savepath = paste0(savepath,"Reactome_ORA_r"))
    enrichReactome_res2 = data.frame(enrichReactome_res)
    
    # if(nrow(enrichReactome_res2)>10){
    #   pdf(paste0(savepath,"GOBP_ORA_r.pdf"), onefile = TRUE)
    #   print(dotplot(enrichGOBP_res, showCategory = 10))
    #   print(cnetplot(enrichGOBP_res, showCategory = 10)) 
    #   dev.off()
      
    # }else{
      pdf(paste0(savepath,"Reactome_ORA_r.pdf"), onefile = TRUE)
      print(dotplot(enrichReactome_res, showCategory = nrow(enrichReactome_res2)))
      print(cnetplot(enrichReactome_res, showCategory = nrow(enrichReactome_res2))) 
      dev.off()
    # }
  }, silent = T)
  
    # ------------------------------------------------------------------------------
  # 4) 这些基因在不同癌症中的ITS数目统计（饼图热图）.
  #    注释基因功能
  # source("07_plotting_v2/scripts/functions/ITS_g_cancer_plot_function.R")
  
  try({
  ITS_g_cancer_plot(ITS_set = ITS_p_set, 
                    genelist = glist,
                    type = 'r',
                    enrichRes = enrichKEGG_res2,
                    savepath = savepath,
                    filename = "gene_cancer_ITSrate_stat_r")
  }, silent = T)

}


# ------------------------------------------------------------------------------
# Drug target analysis
# ------------------------------------------------------------------------------
load( "07_plotting_v2/data/drugTarget_merged3DB.Rdata")
# ------------------------------------------------------------------------------

dir.create("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/drug_merge3DB")
dir.create("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/drug_merge3DB/GSEAenrich/")
dir.create("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSn_sharedgene/drug_merge3DB")
dir.create("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSn_sharedgene/drug_merge3DB/GSEAenrich/")

# 1) drugTarget_drugbank-Drug -------------------------------------------------------------------
filepath <- ("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/drug_merge3DB/GSEAenrich/ITSp_r_vs_drugtarget")

glist <- unlist(strsplit(res2[res2$ID == "ITSp_vs_drugtarget",]$core_enrichment, "/"))

# DT_druglist <- unique(drugTarget_all[is.element(drugTarget_all$genename, glist),]$drugname_lower)
# DT_drug <- drugTarget_all[is.element(drugTarget_all$drugname_lower, DT_druglist),]
DT_drug <- unique(drugTarget_all[is.element(drugTarget_all$genename, glist),])

write.table(DT_drug, paste0(filepath, ".txt"), sep = '\t', quote = F, row.names = F)

# DT_druglistFDA <- unique(drugTarget_FDA[is.element(drugTarget_FDA$genename, glist),]$drugname_lower)
# DT_drugFDA<- drugTarget_FDA[is.element(drugTarget_FDA$drugname_lower, DT_druglist),]
DT_drugFDA <- unique(drugTarget_FDA[is.element(drugTarget_FDA$genename, glist),])
write.table(DT_drugFDA, paste0(filepath, "_FDA.txt"), sep = '\t', quote = F, row.names = F)

DT_drug2 = unique(DT_drug[c('drugname_lower', 'genename')])
DT_drugFDA2 = unique(DT_drugFDA[c('drugname_lower', 'genename', 'ifCancer_drugbank')])

drugstat <- data.frame(group = c("num_drug",  "num_drugFDA", "num_drugFDAcancer"),
                       value = c(length(unique(DT_drug2$drugname_lower)),
                                 length(unique(DT_drugFDA2$drugname_lower)), 
                                 length(unique(DT_drugFDA2[which(DT_drugFDA2$ifCancer_drugbank == 'yes'), ]$drugname_lower))))

write.table(drugstat, paste0(filepath, "_stat.txt"), sep = '\t', quote = F, row.names = F)

# unique(DT_drugFDAcancer$name)
# sort(table(DT_drugFDAcancer$genename))


# 2) OG-Drug -------------------------------------------------------------------
filepath <- ("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/drug_merge3DB/GSEAenrich/ITSp_r_vs_OG_oncoKB")

glist <- unlist(strsplit(res2[res2$ID == "ITSp_vs_OG_oncoKB",]$core_enrichment, "/"))

# OG_druglist <- unique(drugTarget_all[is.element(drugTarget_all$genename, glist),]$drugname_lower)
# OG_drug <- drugTarget_all[is.element(drugTarget_all$drugname_lower, OG_druglist),]
OG_drug <- unique(drugTarget_all[is.element(drugTarget_all$genename, glist),])

write.table(OG_drug, paste0(filepath, ".txt"), sep = '\t', quote = F, row.names = F)

OG_druglistFDA <- unique(drugTarget_FDA[is.element(drugTarget_FDA$genename, glist),]$drugname_lower)
# OG_drugFDA<- drugTarget_FDA[is.element(drugTarget_FDA$drugname_lower, OG_druglist),]
OG_drugFDA<- unique(drugTarget_FDA[is.element(drugTarget_FDA$genename, glist),])


write.table(OG_drugFDA, paste0(filepath, "_FDA.txt"), sep = '\t', quote = F, row.names = F)

OG_drug2 = unique(OG_drug[c('drugname_lower', 'genename')])
OG_drugFDA2 = unique(OG_drugFDA[c('drugname_lower', 'genename', 'ifCancer_drugbank')])

drugstat <- data.frame(group = c("num_drug",  "num_drugFDA", "num_drugFDAcancer"),
                       value = c(length(unique(OG_drug2$drugname_lower)),
                                 length(unique(OG_drugFDA2$drugname_lower)), 
                                 length(unique(OG_drugFDA2[which(OG_drugFDA2$ifCancer_drugbank == 'yes'), ]$drugname_lower))))

write.table(drugstat, paste0(filepath, "_stat.txt"), sep = '\t', quote = F, row.names = F)




# 3) TSG-Drug -------------------------------------------------------------------
filepath <- ("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/drug_merge3DB/GSEAenrich/ITSp_r_vs_TSG_oncoKB")

glist <- unlist(strsplit(res2[res2$ID == "ITSp_vs_TSG_oncoKB",]$core_enrichment, "/"))


# TSG_druglist <- unique(drugTarget_all[is.element(drugTarget_all$genename, glist),]$drugname_lower)
# TSG_drug <- drugTarget_all[is.element(drugTarget_all$drugname_lower, TSG_druglist),]
TSG_drug <- unique(drugTarget_all[is.element(drugTarget_all$genename, glist),])


write.table(TSG_drug, paste0(filepath, ".txt"), sep = '\t', quote = F, row.names = F)

# TSG_druglistFDA <- unique(drugTarget_FDA[is.element(drugTarget_FDA$genename, glist),]$drugname_lower)
# TSG_drugFDA<- drugTarget_FDA[is.element(drugTarget_FDA$drugname_lower, TSG_druglist),]
TSG_drugFDA <- unique(drugTarget_FDA[is.element(drugTarget_FDA$genename, glist),])

write.table(TSG_drugFDA, paste0(filepath, "_FDA.txt"), sep = '\t', quote = F, row.names = F)

TSG_drug2 = unique(TSG_drug[c('drugname_lower', 'genename')])
TSG_drugFDA2 = unique(TSG_drugFDA[c('drugname_lower', 'genename', 'ifCancer_drugbank')])

drugstat <- data.frame(group = c("num_drug",  "num_drugFDA", "num_drugFDAcancer"),
                       value = c(length(unique(TSG_drug2$drugname_lower)),
                                 length(unique(TSG_drugFDA2$drugname_lower)), 
                                 length(unique(TSG_drugFDA2[which(TSG_drugFDA2$ifCancer_drugbank == 'yes'), ]$drugname_lower))))

write.table(drugstat, paste0(filepath, "_stat.txt"), sep = '\t', quote = F, row.names = F)



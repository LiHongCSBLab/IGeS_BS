# rm(list=ls())
work_path <- "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/"
setwd(work_path)
options(stringsAsFactors = F)
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
source("06_results_summaryandplotting/scripts/immunesig_selection_function.R")
source("07_plotting_v2/scripts/functions/ITS_g_cancer_plot_function.R")
source("07_plotting_v2/scripts/functions/ITS_enrich.R")
source("07_plotting_v2/scripts/functions/ITS_function_detection.R")

dir.create("07_plotting_v2/03_ITScharacteristic_metaPlusFilter")

# # ITS characteristics ----------------------------------------------------------
# # number of genes in the corresponding ITS -------------------------------------
# # ITS varies across cancer types -----------------------------------------------------

# files <- list.files("06_results_summaryandplotting/immuneSigGene_ITS_found_lincs_20220414/TUMERIC_spearman_r0.4_ITSproof_2/")

# ITSgenestat_selected <- lapply(as.list(files), function(f){
#   df =  as.data.frame(readxl::read_excel(paste0("06_results_summaryandplotting/immuneSigGene_ITS_found_lincs_20220414/TUMERIC_spearman_r0.4_ITSproof_2/", f), sheet = 1))
#   df = df[c("immune_sig", 
#             # "gene_in_OriginalSig",
#             # "gene_in_OriginalSig_expressedinCCLE",
#             "gene_in_ITS_p",
#             "gene_in_OriginalSig_and_ITS_p")]
#   #  "gene_in_OriginalSig_lincs_overlap", 
#   #  "gene_in_OriginalSig_lincs_all_overlap"
#   df_selected = merge(df, immunesigMeta[c('immune_sig','flag')], by = 'immune_sig')
#   df_selected$gene_in_ITS_not_OriginalSig = df_selected$gene_in_ITS_p - df_selected$gene_in_OriginalSig_and_ITS_p

#   df_selected = melt(df_selected[c("immune_sig", 
#                                   "gene_in_ITS_not_OriginalSig",
#                                   "gene_in_OriginalSig_and_ITS_p")])
#   df_selected$cancer = gsub(".xlsx","", f)
#   return(df_selected)
#   })

# ITSgenestat_selected2 = do.call(rbind, ITSgenestat_selected)

# ITSgenestat_selected2 %>% na.omit() %>% 
#   group_by(variable,immune_sig) %>% 
#   summarise(mean_value=mean(value),
#             sd_value=sd(value)) -> df1

# df1

# df1 %>% 
#   group_by(immune_sig) %>% 
#   mutate(new_col=cumsum(mean_value)) -> df2

# df2$variable<-factor(df2$variable,
#                 levels = unique(df2$variable))

# p_ITSgenestat <- ggplot(data=df2,aes(x=immune_sig,
#                     y=mean_value,
#                     fill=variable))+
#   geom_bar(position = "stack",stat="identity")+
#   geom_errorbar(aes(ymin=mean_value-sd_value,
#                     ymax=mean_value+sd_value),
#                 width=0.1)+
#   # scale_y_continuous(expand = c(0,0),
#   #                    limits = c(0,100))+
#   # scale_fill_material_d()+
#   theme_bw()+coord_flip()+
#   labs(x=NULL,y=NULL)

# p_all4 = p_all3 %>%
# insert_right(p_ITSgenestat, width  = .4)

# ggsave("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/meta_pancancer_genestat.pdf", p_all4, height = 45, width = 30)


# ITS function -----------------------------------------------------------------

ITS_p_s_df <- read.csv("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/ITS_p_s_df.txt", sep='\t')

  # ------------------------------------------------------------------------------
  # 1) functional enrichment
  # ITSp_vs_OG_oncoKB
  dir.create("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/ITS_s_function_enrich")
  savepath = "07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/ITS_s_function_enrich/"
  
  # run ORA analysis -----------------------------------------------------------
  # 0) extract gene list
  glist <- ITS_p_s_df$gene_name
  
  
  # enrich hallmark-ORA ---------------------------------------------------------------
  try({
    
    enrichHallmark_res = ITS_GS_enrich(glist,
                                       pvalueCutoff = 0.05,
                                       type = "H",
                                       method = 'ORA',
                                       savepath = paste0(savepath,"hallmark_ORA_s"))
    enrichHallmark_res2 = data.frame(enrichHallmark_res)
    # unlist(strsplit(enrichHallmark_res2[enrichHallmark_res2$ID == "HALLMARK_KRAS_SIGNALING_UP",]$geneID, "/"))
    
    
      pdf(paste0(savepath,"hallmark_ORA_s.pdf"), onefile = TRUE)
      print(dotplot(enrichHallmark_res, showCategory = nrow(enrichHallmark_res2)))
      print(cnetplot(enrichHallmark_res, showCategory = nrow(enrichHallmark_res2)))
      dev.off()
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
      pdf(paste0(savepath,"KEGG_ORA_s.pdf"), onefile = TRUE)
      print(dotplot(enrichKEGG_res, showCategory = nrow(enrichKEGG_res2)) )
      print(cnetplot(enrichKEGG_res, showCategory = nrow(enrichKEGG_res2)) )
      dev.off()
  }, silent = T)
  
  # enrich GOBP-ORA ---------------------------------------------------------------
  try({
    enrichGOBP_res = ITS_GS_enrich(glist,
                                   pvalueCutoff = 0.05,
                                   type = "GOBP",
                                   method = 'ORA',
                                   savepath = paste0(savepath,"GOBP_ORA_s"))
    enrichGOBP_res2 = data.frame(enrichGOBP_res)
      pdf(paste0(savepath,"GOBP_ORA_s.pdf"), onefile = TRUE)
      print(dotplot(enrichGOBP_res, showCategory = nrow(enrichGOBP_res2)))
      print(cnetplot(enrichGOBP_res, showCategory = nrow(enrichGOBP_res2))) 
      dev.off()
  }, silent = T)
  
  # enrich Reactome-ORA ---------------------------------------------------------------
  try({
    enrichReactome_res = ITS_GS_enrich(glist,
                                   pvalueCutoff = 0.05,
                                   type = "Reactome",
                                   method = 'ORA',
                                   savepath = paste0(savepath,"Reactome_ORA_s"))
    enrichReactome_res2 = data.frame(enrichReactome_res)
    
      pdf(paste0(savepath,"Reactome_ORA_s.pdf"), onefile = TRUE)
      print(dotplot(enrichReactome_res, showCategory = nrow(enrichReactome_res2)))
      print(cnetplot(enrichReactome_res, showCategory = nrow(enrichReactome_res2))) 
      dev.off()
  }, silent = T)
  
  # ------------------------------------------------------------------------------
  # 4) 这些基因在不同癌症中的ITS数目统计（饼图热图）.
  #    注释基因功能
  # source("07_plotting_v2/scripts/functions/ITS_g_cancer_plot_function.R")
  # try({

  # ITS_g_cancer_plot(ITS_set = ITS_p_set, 
  #                   genelist = glist,
  #                   type = 's',
  #                   enrichRes = enrichKEGG_res2,
  #                   savepath = savepath,
  #                   filename = "gene_cancer_ITSrate_stat_s")
  # # }, silent = T)

  # run GSEA analysis -----------------------------------------------------------
  # 0) extract gene list
  tmp2 = as.vector(scale(ITS_p_s_df$weighted_rate))
  names(tmp2) = ITS_p_s_df$gene_name
  
  # enrich hallmark-GSEA ---------------------------------------------------------------
  try({    

    enrichHallmark_res = ITS_GS_enrich(tmp2,
                                       pvalueCutoff = 0.05,
                                       type = "H",
                                       method = 'GSEA',
                                       savepath = paste0(savepath,"hallmark_GSEA_s"))
    enrichHallmark_res2 = data.frame(enrichHallmark_res)
    # unlist(strsplit(enrichHallmark_res2[enrichHallmark_res2$ID == "HALLMARK_KRAS_SIGNALING_UP",]$geneID, "/"))
      pdf(paste0(savepath,"hallmark_GSEA_s.pdf"), onefile = TRUE)
      print(dotplot(enrichHallmark_res, showCategory = nrow(enrichHallmark_res2)))
      print(cnetplot(enrichHallmark_res, showCategory = nrow(enrichHallmark_res2)))
      dev.off()
    
  }, silent = T)
  
  # enrich KEGG-GSEA---------------------------------------------------------------
  try({
    enrichKEGG_res = ITS_GS_enrich(tmp2,
                                   pvalueCutoff = 0.05,
                                   type = "KEGG",
                                   method = 'GSEA',
                                   savepath = paste0(savepath,"KEGG_GSEA_s"))
    
    enrichKEGG_res2 = data.frame(enrichKEGG_res)
    # unlist(strsplit(enrichKEGG_res2[enrichKEGG_res2$Description == "NF-kappa B signaling pathway",]$geneID, "/"))
    # unlist(strsplit(enrichKEGG_res2[enrichKEGG_res2$Description == "VEGF signaling pathway",]$geneID, "/"))
    
      pdf(paste0(savepath,"KEGG_GSEA_s.pdf"), onefile = TRUE)
      print(dotplot(enrichKEGG_res, showCategory = nrow(enrichKEGG_res2)) )
      print(cnetplot(enrichKEGG_res, showCategory = nrow(enrichKEGG_res2)) )
      dev.off()
  }, silent = T)
  
  # enrich GOBP-GSEA ---------------------------------------------------------------
  try({
    enrichGOBP_res = ITS_GS_enrich(tmp2,
                                   pvalueCutoff = 0.05,
                                   type = "GOBP",
                                   method = 'GSEA',
                                   savepath = paste0(savepath,"GOBP_GSEA_s"))
    enrichGOBP_res2 = data.frame(enrichGOBP_res)
    
      pdf(paste0(savepath,"GOBP_GSEA_s.pdf"), onefile = TRUE)
      print(dotplot(enrichGOBP_res, showCategory = nrow(enrichGOBP_res2)))
      print(cnetplot(enrichGOBP_res, showCategory = nrow(enrichGOBP_res2))) 
      dev.off()
  }, silent = T)
  
  # enrich Reactome-GSEA ---------------------------------------------------------------
  try({
    enrichReactome_res = ITS_GS_enrich(tmp2,
                                   pvalueCutoff = 0.05,
                                   type = "Reactome",
                                   method = 'GSEA',
                                   savepath = paste0(savepath,"Reactome_GSEA_s"))
    enrichReactome_res2 = data.frame(enrichReactome_res)
    
      pdf(paste0(savepath,"Reactome_GSEA_s.pdf"), onefile = TRUE)
      print(dotplot(enrichReactome_res, showCategory = nrow(enrichReactome_res2)))
      print(cnetplot(enrichReactome_res, showCategory = nrow(enrichReactome_res2))) 
      dev.off()
  }, silent = T)



# ------------------------------------------------------------------------------
ITS_p_r_df <- read.csv("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/ITS_p_r_df.txt", sep='\t')

  # ------------------------------------------------------------------------------
  # 1) functional enrichment
  # ITSp_vs_OG_oncoKB
  dir.create("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/ITS_r_function_enrich/")
  savepath = "07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/ITS_r_function_enrich/"
  
  # ------------------------------------------------------------------------------
  # 0) extract gene list
  glist <- ITS_p_r_df$gene_name
  
  # # ------------------------------------------------------------------------------
  # # 1) if these genes in original immune signatures
  # # try({
  #   length(glist)
  #   length(unique(immunesigGS_S[is.element(immunesigGS_S$gene_name, glist),]$gene_name))
  #   setdiff(glist, immunesigGS_S$gene_name)
  #   library(VennDiagram)
  #   venn.diagram(x=list(ITSgene=glist, 
  #                       originalSignaturesGenes=unique(immunesigGS_S$gene_name)), 
  #                filename = paste0(savepath, "gene_S_Venn.png"), 
  #                height = 450, width = 450, resolution =300, imagetype="png", 
  #                col="white", fill=c(colors()[616], colors()[38]), 
  #                alpha=c(0.6, 0.6), lwd=c(1, 1), cex=0.3, 
  #                cat.dist=c(-0.07, -0.03), cat.pos=c(300, 60), cat.cex=0.3) 
  # # }, silent = T)
  #     dev.off()

  # ------------------------------------------------------------------------------
  # 2) correlation between gene expression and original immune signatures 
  # paste0(07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITS_correlation_sgdotplot/", gsname)
  
  
  # ------------------------------------------------------------------------------
  # 3) functional enrichment
  
  # enrich hallmark-ORA ---------------------------------------------------------------
  try({
    
    enrichHallmark_res = ITS_GS_enrich(glist,
                                       pvalueCutoff = 0.05,
                                       type = "H",
                                       method = 'ORA',
                                       savepath = paste0(savepath,"hallmark_ORA_r"))
    enrichHallmark_res2 = data.frame(enrichHallmark_res)
    # unlist(strsplit(enrichHallmark_res2[enrichHallmark_res2$ID == "HALLMARK_KRAS_SIGNALING_UP",]$geneID, "/"))
      pdf(paste0(savepath,"hallmark_ORA_r.pdf"), onefile = TRUE)
      print(dotplot(enrichHallmark_res, showCategory = nrow(enrichHallmark_res2)))
      print(cnetplot(enrichHallmark_res, showCategory = nrow(enrichHallmark_res2)))
      dev.off()
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
      pdf(paste0(savepath,"KEGG_ORA_r.pdf"), onefile = TRUE)
      print(dotplot(enrichKEGG_res, showCategory = nrow(enrichKEGG_res2)) )
      print(cnetplot(enrichKEGG_res, showCategory = nrow(enrichKEGG_res2)) )
      dev.off()
  }, silent = T)
  
  # enrich GOBP-ORA ---------------------------------------------------------------
  try({
    enrichGOBP_res = ITS_GS_enrich(glist,
                                   pvalueCutoff = 0.05,
                                   type = "GOBP",
                                   method = 'ORA',
                                   savepath = paste0(savepath,"GOBP_ORA_r"))
    enrichGOBP_res2 = data.frame(enrichGOBP_res)
      pdf(paste0(savepath,"GOBP_ORA_r.pdf"), onefile = TRUE)
      print(dotplot(enrichGOBP_res, showCategory = nrow(enrichGOBP_res2)))
      print(cnetplot(enrichGOBP_res, showCategory = nrow(enrichGOBP_res2))) 
      dev.off()
  }, silent = T)
  
  # enrich Reactome-ORA ---------------------------------------------------------------
  try({
    enrichReactome_res = ITS_GS_enrich(glist,
                                   pvalueCutoff = 0.05,
                                   type = "Reactome",
                                   method = 'ORA',
                                   savepath = paste0(savepath,"Reactome_ORA_r"))
    enrichReactome_res2 = data.frame(enrichReactome_res)
      pdf(paste0(savepath,"Reactome_ORA_r.pdf"), onefile = TRUE)
      print(dotplot(enrichReactome_res, showCategory = nrow(enrichReactome_res2)))
      print(cnetplot(enrichReactome_res, showCategory = nrow(enrichReactome_res2))) 
      dev.off()
  }, silent = T)
  
  # # ------------------------------------------------------------------------------
  # # 4) 这些基因在不同癌症中的ITS数目统计（饼图热图）.
  # #    注释基因功能
  # # source("07_plotting_v2/scripts/functions/ITS_g_cancer_plot_function.R")
  # # try({

  # ITS_g_cancer_plot(ITS_set = ITS_p_set, 
  #                   genelist = glist,
  #                   type = 's',
  #                   enrichRes = enrichKEGG_res2,
  #                   savepath = savepath,
  #                   filename = "gene_cancer_ITSrate_stat_s")
  # # }, silent = T)


  # run GSEA analysis -----------------------------------------------------------
  # 0) extract gene list
  tmp2 = as.vector(scale(ITS_p_r_df$weighted_rate))
  names(tmp2) = ITS_p_r_df$gene_name
  
  # enrich hallmark-GSEA ---------------------------------------------------------------
  try({    

    enrichHallmark_res = ITS_GS_enrich(tmp2,
                                       pvalueCutoff = 0.05,
                                       type = "H",
                                       method = 'GSEA',
                                       savepath = paste0(savepath,"hallmark_GSEA_r"))
    enrichHallmark_res2 = data.frame(enrichHallmark_res)
    # unlist(strsplit(enrichHallmark_res2[enrichHallmark_res2$ID == "HALLMARK_KRAS_SIGNALING_UP",]$geneID, "/"))
      pdf(paste0(savepath,"hallmark_GSEA_r.pdf"), onefile = TRUE)
      print(dotplot(enrichHallmark_res, showCategory = nrow(enrichHallmark_res2)))
      print(cnetplot(enrichHallmark_res, showCategory = nrow(enrichHallmark_res2)))
      dev.off()
    
  }, silent = T)
  
  # enrich KEGG-GSEA---------------------------------------------------------------
  try({
    enrichKEGG_res = ITS_GS_enrich(glist,
                                   pvalueCutoff = 0.05,
                                   type = "KEGG",
                                   method = 'GSEA',
                                   savepath = paste0(savepath,"KEGG_GSEA_r"))
    
    enrichKEGG_res2 = data.frame(enrichKEGG_res)
    # unlist(strsplit(enrichKEGG_res2[enrichKEGG_res2$Description == "NF-kappa B signaling pathway",]$geneID, "/"))
    # unlist(strsplit(enrichKEGG_res2[enrichKEGG_res2$Description == "VEGF signaling pathway",]$geneID, "/"))
    
      pdf(paste0(savepath,"KEGG_GSEA_r.pdf"), onefile = TRUE)
      print(dotplot(enrichKEGG_res, showCategory = nrow(enrichKEGG_res2)) )
      print(cnetplot(enrichKEGG_res, showCategory = nrow(enrichKEGG_res2)) )
      dev.off()
  }, silent = T)
  
  # enrich GOBP-GSEA ---------------------------------------------------------------
  try({
    enrichGOBP_res = ITS_GS_enrich(glist,
                                   pvalueCutoff = 0.05,
                                   type = "GOBP",
                                   method = 'GSEA',
                                   savepath = paste0(savepath,"GOBP_GSEA_r"))
    enrichGOBP_res2 = data.frame(enrichGOBP_res)
    
      pdf(paste0(savepath,"GOBP_GSEA_r.pdf"), onefile = TRUE)
      print(dotplot(enrichGOBP_res, showCategory = nrow(enrichGOBP_res2)))
      print(cnetplot(enrichGOBP_res, showCategory = nrow(enrichGOBP_res2))) 
      dev.off()
  }, silent = T)
  
  # enrich Reactome-GSEA ---------------------------------------------------------------
  try({
    enrichReactome_res = ITS_GS_enrich(glist,
                                   pvalueCutoff = 0.05,
                                   type = "Reactome",
                                   method = 'GSEA',
                                   savepath = paste0(savepath,"Reactome_GSEA_r"))
    enrichReactome_res2 = data.frame(enrichReactome_res)
    
      pdf(paste0(savepath,"Reactome_GSEA_r.pdf"), onefile = TRUE)
      print(dotplot(enrichReactome_res, showCategory = nrow(enrichReactome_res2)))
      print(cnetplot(enrichReactome_res, showCategory = nrow(enrichReactome_res2))) 
      dev.off()
  }, silent = T)

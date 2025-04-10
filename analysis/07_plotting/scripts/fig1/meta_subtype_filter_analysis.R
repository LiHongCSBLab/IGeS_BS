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
library(tidyverse)



source("03_drug_immuneSig_enrichment/scripts/02_GSEA_geneset_preparation_function.R")
source("03_drug_immuneSig_enrichment/scripts/02_GSEA_geneset_merge_metaAnalysis_function.R")

source("06_results_summaryandplotting/scripts/immunesig_selection_function.R")
source("06_results_summaryandplotting/scripts/ITS_analysis/ITS_function_detection.R")
source("06_results_summaryandplotting/scripts/ITS_analysis/ITSexpr_purity_detection.R")
source("06_results_summaryandplotting/scripts/ITS_analysis/ITSexpr_purity_detection.R")
source("06_results_summaryandplotting/scripts/ITS_analysis/immuneModulatorDetection.R")



immunesig <- read.csv("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/05_immuneSig_ICBresponse/results/meta_results/res_metaRes.csv")
immunesig <- immunesig[immunesig$Meta_Pval < 0.05, ]
names(immunesig) <- c("immune_sig", names(immunesig)[-1])
immunesig$flag <- "sensitive"
immunesig[immunesig$OR < 1,]$flag <- "resistant"
immunesig$setindex = "01_sensitive"
immunesig[immunesig$flag == "resistant", ]$setindex = "02_resistant"

immuneSigInfoPath = "06_results_summaryandplotting/data/immene_cell_hieracy/"

immunesigInfo <- as.data.frame(readxl::read_excel(paste0(immuneSigInfoPath, "immune_sig_hieracy_new.xlsx"), sheet = 1))
immunesigInfo <- immunesigInfo[c(1,2,4,5:9)]
immunesigInfo$immune_sig = tolower(immunesigInfo$immune_sig)
immunesigInfo = immunesig_name_unify(immunesigInfo)

immunesigAnnot <- merge(unique(immunesig), immunesigInfo, by='immune_sig')
rownames(immunesigAnnot) = immunesigAnnot$immune_sig


# ------------------------------------------------------------------------------
load("05_immuneSig_ICBresponse/results/oriImmuneSig_res_all.Rdata")
merged_sig_scaled <- lapply(as.list(names(merged_sig)), function(x){
     df = merged_sig[[x]]
     print(x)
     df_scaled <- apply(df[-c(1:5)], 2, scale)
     df_scaled[is.nan(df_scaled)] <- 0
     df_scaled <- cbind(df[1:5], df_scaled)
     df_scaled <- df_scaled[!is.na(df_scaled$label), ]
     df_scaled$tag <- 1
     df_scaled[df_scaled$label == "NR",]$tag = 0
     df_scaled$tag <- factor(df_scaled$tag,levels = c(0,1), order=TRUE)
     return(df_scaled)
})

names(merged_sig_scaled) <- names(merged_sig)
merged_sig_scaled <- do.call(rbind,merged_sig_scaled)

rownames(merged_sig_scaled) <- paste0(merged_sig_scaled$dataset, merged_sig_scaled$patient)
res_df <-merged_sig_scaled[!is.na(merged_sig_scaled$label), ]
res_df = res_df[order(res_df$label), ]


xcell_remove <- xcell_sig_remover()
res_df <- res_df[-which(is.element(names(res_df), tolower(xcell_remove)))]

res_df[is.na(res_df)] = 0

res_df_annot <- res_df[c(1:5)]
res_df_metaselect <- res_df[is.element(names(res_df), immunesigAnnot$immune_sig)]
# res_df_metaselect <- cbind(res_df_annot, res_df_metaselect)

sig_cor = cor(res_df_metaselect, method='spearman')
immunesigAnnot = immunesigAnnot[order(immunesigAnnot$SubType), ]
sig_cor = sig_cor[immunesigAnnot$immune_sig, immunesigAnnot$immune_sig]
# pheatmap(sig_cor, cluster_row = F, annotation_row = immunesigAnnot[c( 'SubType','Type','flag')])


immunesigAnnot_s = immunesigAnnot[is.element(immunesigAnnot$flag, "sensitive"), ]
immunesigAnnot_s = immunesigAnnot_s[order(immunesigAnnot_s$SubType), ]
sig_cor_s = sig_cor[immunesigAnnot_s$immune_sig, immunesigAnnot_s$immune_sig]

sigs_s = lapply(as.list(seq(length(unique(immunesigAnnot_s$SubType)))), function(i){
     # i =2
     subtype = unique(immunesigAnnot_s$SubType)[i]
     print(subtype)
     # print(tmp)
     if(length(immunesigAnnot_s[immunesigAnnot_s$SubType == subtype,]$immune_sig) == 1){
     sigs = immunesigAnnot_s[is.element(immunesigAnnot_s$immune_sig, immunesigAnnot_s[immunesigAnnot_s$SubType == subtype,]$immune_sig), ]
     return(sigs)

     }else if(subtype == 'Miscellaneous'){
     sigs = immunesigAnnot_s[is.element(immunesigAnnot_s$immune_sig, immunesigAnnot_s[immunesigAnnot_s$SubType == subtype,]$immune_sig), ]
     return(sigs)

     }else {
     tmp = sig_cor_s[immunesigAnnot_s[immunesigAnnot_s$SubType == subtype,]$immune_sig, immunesigAnnot_s[immunesigAnnot_s$SubType == subtype,]$immune_sig]
     pheatmap(tmp, display_numbers=T, number_format = "%.2f",fontsize_number=12)
     # tmp
     sigs = immunesigAnnot_s[is.element(immunesigAnnot_s$immune_sig, colnames(tmp)), ]
     sigs = sigs[order(sigs$OR, decreasing = T), ]
     # print(sigs)
     return(head(sigs,1))
     }
})

sigs_s = do.call(rbind, sigs_s)
pheatmap(sig_cor[sigs_s$immune_sig, sigs_s$immune_sig], annotation_row = immunesigAnnot[c( 'SubType','Type','flag')],
         display_numbers=T, number_format = "%.2f",fontsize_number=5)

sig_cor_s = sig_cor[sigs_s$immune_sig, sigs_s$immune_sig]
sig_cor_s = melt(sig_cor_s[lower.tri(sig_cor_s)])
sig_cor_s = sig_cor_s[sig_cor_s$value !=0, ]


immunesigAnnot_r = immunesigAnnot[is.element(immunesigAnnot$flag, "resistant"), ]
immunesigAnnot_r = immunesigAnnot_r[order(immunesigAnnot_r$SubType), ]
sig_cor_r = sig_cor[immunesigAnnot_r$immune_sig, immunesigAnnot_r$immune_sig]

sigs_r = lapply(as.list(seq(length(unique(immunesigAnnot_r$SubType)))), function(i){
     # i = 5
     subtype = unique(immunesigAnnot_r$SubType)[i]
     print(subtype)
     # print(tmp)
     if(length(immunesigAnnot_r[immunesigAnnot_r$SubType == subtype,]$immune_sig) == 1){
     sigs = immunesigAnnot_r[is.element(immunesigAnnot_r$immune_sig, immunesigAnnot_r[immunesigAnnot_r$SubType == subtype,]$immune_sig), ]
     return(sigs)

     }else if(subtype == 'Miscellaneous'){
     sigs = immunesigAnnot_r[is.element(immunesigAnnot_r$immune_sig, immunesigAnnot_r[immunesigAnnot_r$SubType == subtype,]$immune_sig), ]
     return(sigs)

     }else{
     tmp = sig_cor_r[immunesigAnnot_r[immunesigAnnot_r$SubType == subtype,]$immune_sig, immunesigAnnot_r[immunesigAnnot_r$SubType == subtype,]$immune_sig]
     pheatmap(tmp)

     sigs = immunesigAnnot_r[is.element(immunesigAnnot_r$immune_sig, colnames(tmp)), ]
     sigs = sigs[order(sigs$OR, decreasing = T), ]
     # print(sigs)
     return(head(sigs,1))
     }
})

sigs_r = do.call(rbind, sigs_r)
pheatmap(sig_cor[sigs_r$immune_sig, sigs_r$immune_sig], annotation_row = immunesigAnnot[c( 'SubType','Type','flag')])


sig_select <- rbind(sigs_s,sigs_r)
pheatmap(sig_cor[sig_select$immune_sig, sig_select$immune_sig], annotation_row = immunesigAnnot[c( 'SubType','Type','flag')])

# Subtype
write.table(sig_select, "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/07_plotting_v2/01_metaAnalysis/table/metaSelected_immunesigAnnot_filter.txt", sep = '\t', quote = F, row.names = F)
# # SubType_detail
# write.table(sig_select, "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/07_plotting_v2/01_metaAnalysis/table/metaSelected_immunesigAnnot_subtypedetail_filter.txt", sep = '\t', quote = F, row.names = F)



dir.create("07_plotting_v2/01_metaAnalysis/")
meta_result <- read.csv("05_immuneSig_ICBresponse/results/meta_results/res_metaRes.csv")
meta_selected <- meta_result[meta_result$Meta_Pval < 0.05, ]
meta_selected$immune_sig = meta_selected$X
# meta_selected <- meta_selected[-c(grep("caf_tide", meta_selected$immune_sig),
#                                   grep("m2_tide", meta_selected$immune_sig)), ]

immunesigMeta = unique(meta_selected[c('immune_sig', 'OR', "lower_OR","upper_OR", "Meta_Pval")]) 
immunesigMeta$flag = 'sensitive'
immunesigMeta[immunesigMeta$OR < 1, ]$flag = 'resistant'

# Fig1c: mean number of genes in the signatures  --------------------------------------------
# number of genes in the signatures expressed in CCLE data --------------------------------------------
files <- list.files("06_results_summaryandplotting/immuneSigGene_ITS_found_lincs_20220414/TUMERIC_spearman_r0.4_ITSproof_2/")
files <- c(files, "SKCM_Primary.xlsx")

genestat_selected <- lapply(as.list(files), function(f){
  if(f == "SKCM_Primary.xlsx"){
    df =  as.data.frame(readxl::read_excel(paste0("06_results_summaryandplotting/immuneSigGene_ITS_found_lincs_20220414/CPE_spearman_r0.4_ITSproof_2/", f), sheet = 1))
    
  }else{
    df =  as.data.frame(readxl::read_excel(paste0("06_results_summaryandplotting/immuneSigGene_ITS_found_lincs_20220414/TUMERIC_spearman_r0.4_ITSproof_2/", f), sheet = 1))
  }
  df = df[c("immune_sig", 
            "gene_in_OriginalSig",
            "gene_in_OriginalSig_expressedinCCLE")]
  #  "gene_in_OriginalSig_lincs_overlap", 
  #  "gene_in_OriginalSig_lincs_all_overlap"
  df_selected = merge(df, immunesigMeta[c('immune_sig','flag')], by = 'immune_sig')
  df_selected$gene_in_OriginalSig_not_expressedinCCLE = df_selected$gene_in_OriginalSig - df_selected$gene_in_OriginalSig_expressedinCCLE
  
  df_selected$gene_in_OriginalSig_expressedinCCLE_rate = df_selected$gene_in_OriginalSig_expressedinCCLE / df_selected$gene_in_OriginalSig
  df_selected$gene_in_OriginalSig_not_expressedinCCLE_rate = df_selected$gene_in_OriginalSig_not_expressedinCCLE / df_selected$gene_in_OriginalSig

  df_selected = melt(df_selected[c("immune_sig", 
                                   "gene_in_OriginalSig",
                                   "gene_in_OriginalSig_not_expressedinCCLE",
                                   "gene_in_OriginalSig_expressedinCCLE",
                                   "gene_in_OriginalSig_expressedinCCLE_rate",
                                   "gene_in_OriginalSig_not_expressedinCCLE_rate")])
  df_selected$cancer = gsub(".xlsx","", f)
  return(df_selected)
})

genestat_selected2 = do.call(rbind, genestat_selected)


genestat_selected3 <- genestat_selected2[is.element(genestat_selected2$immune_sig, sig_select$immune_sig), ]

genestat_selected3 %>% na.omit() %>% 
  group_by(variable,immune_sig) %>% 
  summarise(mean_value=mean(value),
            sd_value=sd(value),
            min_value=min(value),
            max_value=max(value))  %>% 
  group_by(immune_sig) %>% 
  mutate(new_col=cumsum(mean_value)) -> df_ORI

df_ORI$variable<-factor(df_ORI$variable, levels = unique(df_ORI$variable))
df_ORI = as.data.frame(df_ORI)
df_ORI[df_ORI$variable == 'gene_in_OriginalSig_expressedinCCLE_rate',]
df_ORI_not_expressedinCCLE_rate = df_ORI[df_ORI$variable == 'gene_in_OriginalSig_not_expressedinCCLE_rate',]
nrow(df_ORI_not_expressedinCCLE_rate)
nrow(df_ORI_not_expressedinCCLE_rate[df_ORI_not_expressedinCCLE_rate$mean_value > 0.5,])
nrow(df_ORI_not_expressedinCCLE_rate[df_ORI_not_expressedinCCLE_rate$mean_value > 0.3,])


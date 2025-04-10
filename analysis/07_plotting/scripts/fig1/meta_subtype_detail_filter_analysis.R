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

sigs_s = lapply(as.list(seq(length(unique(immunesigAnnot_s$SubType_detail)))), function(i){
     # i = 44
     subtype = unique(immunesigAnnot_s$SubType_detail)[i]
     print(subtype)
     # print(tmp)
     if(length(immunesigAnnot_s[immunesigAnnot_s$SubType_detail == subtype,]$immune_sig) == 1){
     sigs = immunesigAnnot_s[is.element(immunesigAnnot_s$immune_sig, immunesigAnnot_s[immunesigAnnot_s$SubType_detail == subtype,]$immune_sig), ]
     # sigs
     return(sigs)

     }else if(subtype == 'Miscellaneous'){
     sigs = immunesigAnnot_s[is.element(immunesigAnnot_s$immune_sig, immunesigAnnot_s[immunesigAnnot_s$SubType_detail == subtype,]$immune_sig), ]
     return(sigs)

     }else {
     tmp = sig_cor_s[immunesigAnnot_s[immunesigAnnot_s$SubType_detail == subtype,]$immune_sig, immunesigAnnot_s[immunesigAnnot_s$SubType_detail == subtype,]$immune_sig]
     pheatmap(tmp, display_numbers=T, number_format = "%.2f",fontsize_number=12)
     # tmp
     sigs = immunesigAnnot_s[is.element(immunesigAnnot_s$immune_sig, colnames(tmp)), ]
     sigs = sigs[order(sigs$OR, decreasing = T), ]
     # print(sigs)
     return(head(sigs,1))
     }
})


sigs_s = do.call(rbind, sigs_s)
sigs_s = rbind(sigs_s, immunesigAnnot_s[immunesigAnnot_s$immune_sig == 't.cell.cd8._epic',])

pheatmap(sig_cor[sigs_s$immune_sig, sigs_s$immune_sig], annotation_row = immunesigAnnot[c( 'SubType','Type','flag')])





immunesigAnnot_r = immunesigAnnot[is.element(immunesigAnnot$flag, "resistant"), ]
immunesigAnnot_r = immunesigAnnot_r[order(immunesigAnnot_r$SubType), ]
sig_cor_r = sig_cor[immunesigAnnot_r$immune_sig, immunesigAnnot_r$immune_sig]

sigs_r = lapply(as.list(seq(length(unique(immunesigAnnot_r$SubType_detail)))), function(i){
     # i = 5
     subtype = unique(immunesigAnnot_r$SubType_detail)[i]
     print(subtype)
     # print(tmp)
     if(length(immunesigAnnot_r[immunesigAnnot_r$SubType_detail == subtype,]$immune_sig) == 1){
     sigs = immunesigAnnot_r[is.element(immunesigAnnot_r$immune_sig, immunesigAnnot_r[immunesigAnnot_r$SubType_detail == subtype,]$immune_sig), ]
     return(sigs)

     }else if(subtype == 'Miscellaneous'){
     sigs = immunesigAnnot_r[is.element(immunesigAnnot_r$immune_sig, immunesigAnnot_r[immunesigAnnot_r$SubType_detail == subtype,]$immune_sig), ]
     return(sigs)

     }else{
     tmp = sig_cor_r[immunesigAnnot_r[immunesigAnnot_r$SubType_detail == subtype,]$immune_sig, immunesigAnnot_r[immunesigAnnot_r$SubType_detail == subtype,]$immune_sig]
     pheatmap(tmp)

     sigs = immunesigAnnot_r[is.element(immunesigAnnot_r$immune_sig, colnames(tmp)), ]
     sigs = sigs[order(sigs$OR, decreasing = T), ]
     # print(sigs)
     return(head(sigs,1))
     }
})

sigs_r = do.call(rbind, sigs_r)
# pheatmap(sig_cor[sigs_r$immune_sig, sigs_r$immune_sig], annotation_row = immunesigAnnot[c( 'SubType','Type','flag')])


sig_select <- rbind(sigs_s,sigs_r)
pheatmap(sig_cor[sig_select$immune_sig, sig_select$immune_sig], annotation_row = immunesigAnnot[c( 'SubType','Type','flag')])

# Subtype
# write.table(sig_select, "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/07_plotting_v2/01_metaAnalysis/table/metaSelected_immunesigAnnot_filter.txt", sep = '\t', quote = F, row.names = F)
# Subtype_detail
write.table(sig_select, "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/07_plotting_v2/01_metaAnalysis/table/metaSelected_immunesigAnnot_subtypedetail_filter.txt", sep = '\t', quote = F, row.names = F)

sig_cor_select = sig_cor[sig_select$immune_sig, sig_select$immune_sig]
sig_cor_select[upper.tri(sig_cor_select)]=NA
sig_cor_select[diag(sig_cor_select)]=NA
sig_cor_select = reshape2::melt(sig_cor_select)
# View(sig_cor_select[!is.na(sig_cor_select$value),])
# View(sig_cor_select[which(sig_cor_select$value > 0.7 & sig_cor_select$value != 1),])

# detect if genes are highly overlapped ----------------------------------------

immuneSigGSpath = "02_tumorGene_immuneSig_correlation/data/immune_signatures/originalSignaturesGenes/"
savepath = "r1_drug_immuneSig_correlation_analysis/"

# LOAD original immune signature gene sets ---------------------------------
immunesigGS1 <- as.data.frame(read_excel(paste0(immuneSigGSpath, "immunesig_original_used.xlsx"), sheet = 1))
immunesigGS1 <- immunesigGS1[-grep("various", immunesigGS1$Genes), ]
immunesigGS1 <- immunesigGS1[-grep("SignatureMatrix", immunesigGS1$Genes), ]
immunesigGS1 <- immunesigGS1[-grep("ReferenceCellExpression", immunesigGS1$Genes), ]
immunesigGS1 <- immunesigGS1[-grep("mergedSignatures", immunesigGS1$Genes), ]
immunesigGS1 <- immunesigGS1[-grep("L22Matrix", immunesigGS1$Genes), ]

immunesigGS1 <- immunesigGS1[-c(2:5)]
immunesigGS1Name <- as.list(immunesigGS1$immune_sig)
immunesigGS1 <- lapply(immunesigGS1Name, function(x){
  y=immunesigGS1[immunesigGS1$immune_sig == x,][-1]
  unlist(y[-which(is.na(y))])
})
names(immunesigGS1) = unlist(immunesigGS1Name)

immunesigGS2 <- read_excel(paste0(immuneSigGSpath, "immunesig_original_used.xlsx"), sheet = 2)
cancer='SKCM'
immunesigGS2 <- as.data.frame(immunesigGS2[immunesigGS2$CancerType == cancer, ])
immunesigGS2Name <- as.list(immunesigGS2$Signatures)
immunesigGS2 <- lapply(immunesigGS2Name, function(x){
  y=immunesigGS2[immunesigGS2$Signatures == x,][-1]
  unlist(y[-which(is.na(y))])
})
names(immunesigGS2) = paste0(unlist(immunesigGS2Name),"_ConsensusTMEgsva")
immunesigGS = c(immunesigGS1, immunesigGS2)
immunesigGSname = data.frame(signame = names(immunesigGS), immune_sig = tolower(names(immunesigGS)))
immunesigGSname = immunesig_name_unify(immunesigGSname)
names(immunesigGS) = immunesigGSname$immune_sig



immunesigGS_selected = immunesigGS[intersect(names(immunesigGS), sig_select$immune_sig)]

jaccard_index <- function(x, y){
  length(intersect(x,y))/length(union(x,y))
}

Otsuka_Ochiai_coefficient <- function(x, y){
    length(intersect(x, y))/sqrt(length(x) *length(y))
}

maxoverlap <- function(x, y){
    max(length(intersect(x, y))/length(x), length(intersect(x, y))/length(y))
}
TMESig_p_jcindex <- list()
for(i in 1:length(immunesigGS_selected)){
  x = immunesigGS_selected[[i]]
  TMESig_p_jcindex[[i]] = sapply(immunesigGS_selected, function(y){
    jaccard_index(x,y)
    # length(intersect(x, y))
  })
}
names(TMESig_p_jcindex) = names(immunesigGS_selected)

# plot JC index to see the pattern -----------------------------------------
TMESig_p_jcindex = do.call(cbind, TMESig_p_jcindex)
TMESig_p_jcindex = TMESig_p_jcindex[sort(colnames(TMESig_p_jcindex)), ]
TMESig_p_jcindex = TMESig_p_jcindex[, sort(colnames(TMESig_p_jcindex))]


# pdf(paste0("06_results_summaryandplotting/data/immunesig_jaccard_index_spearman_r",s , "_ITSproof_",ITSproof,"/oriTMESig_jci.pdf"), 
#     width=18, height = 16)
# rownames(ITS_TMEcellSig_p_jcindex) = paste0(rownames(ITS_TMEcellSig_p_jcindex), " (", unlist(ng), ")")
    num_gene_immunesigGS_selected = lapply(immunesigGS_selected,length)
    ng = num_gene_immunesigGS_selected[colnames(TMESig_p_jcindex)]
#     colnames(TMESig_p_jcindex) = paste0(colnames(TMESig_p_jcindex), " (", unlist(ng), ")")

pheatmap::pheatmap(TMESig_p_jcindex,
                   cluster_cols = T,
                   cluster_rows = T,
                   display_numbers = T,
                   number_format = "%.2f",
                   fontsize = 5, 
                   annotation_row = immunesigAnnot[c( 'SubType','Type','flag')],
                   color = colorRampPalette(c("steelblue", "white", "firebrick3"))(20),
                   border_color = 'grey50' )
# dev.off()
TMESig_p_jcindex[diag(TMESig_p_jcindex)]=NA
TMESig_p_jcindex[upper.tri(TMESig_p_jcindex)]=NA
TMESig_p_jcindex = reshape2::melt(TMESig_p_jcindex)
View(TMESig_p_jcindex[which(TMESig_p_jcindex$value > 0.3),])

TMESig_p_Ochiai <- list()
for(i in 1:length(immunesigGS_selected)){
  x = immunesigGS_selected[[i]]
  TMESig_p_Ochiai[[i]] = sapply(immunesigGS_selected, function(y){
    Otsuka_Ochiai_coefficient(x,y)
    # length(intersect(x, y))
  })
}
names(TMESig_p_Ochiai) = names(immunesigGS_selected)

# plot JC index to see the pattern -----------------------------------------
TMESig_p_Ochiai = do.call(cbind, TMESig_p_Ochiai)
TMESig_p_Ochiai = TMESig_p_Ochiai[sort(colnames(TMESig_p_Ochiai)), ]
TMESig_p_Ochiai = TMESig_p_Ochiai[, sort(colnames(TMESig_p_Ochiai))]
TMESig_p_Ochiai[is.na(TMESig_p_Ochiai)] = 0


# pdf(paste0("06_results_summaryandplotting/data/immunesig_jaccard_index_spearman_r",s , "_ITSproof_",ITSproof,"/oriTMESig_Ochiai.pdf"), 
#     width=18, height = 16)
# rownames(ITS_TMEcellSig_p_Ochiai) = paste0(rownames(ITS_TMEcellSig_p_Ochiai), " (", unlist(ng), ")")
    num_gene_immunesigGS_selected = lapply(immunesigGS_selected,length)
    ng = num_gene_immunesigGS_selected[colnames(TMESig_p_Ochiai)]
#     colnames(TMESig_p_Ochiai) = paste0(colnames(TMESig_p_Ochiai), " (", unlist(ng), ")")

tmp = pheatmap::pheatmap(TMESig_p_Ochiai,
                   cluster_cols = T,
                   cluster_rows = T,
                   display_numbers = T,
                   number_format = "%.2f",
                   fontsize = 5, 
                   annotation_row = immunesigAnnot[c( 'SubType','Type','flag')],
                   color = colorRampPalette(c("steelblue", "white", "firebrick3"))(20),
                   border_color = 'grey50' )
# dev.off()
TMESig_p_Ochiai[upper.tri(TMESig_p_Ochiai)]=0
TMESig_p_Ochiai[diag(TMESig_p_Ochiai)]=NA

TMESig_p_Ochiai_df = reshape2::melt(TMESig_p_Ochiai)
View(TMESig_p_Ochiai_df[TMESig_p_Ochiai_df$value > 0.3,])

names(TMESig_p_jcindex) = c('V1','V2', 'JCindex')
names(TMESig_p_Ochiai_df) = c('V1','V2', 'Ochiaiindex')
names(sig_cor_select) = c('V1','V2', 'cor')
sig_cor_select2 = sig_cor_select[!is.na(sig_cor_select$cor),]
sig_cor_select2 = sig_cor_select2[which(sig_cor_select2$cor < 1),]
test = inner_join(inner_join(TMESig_p_jcindex, TMESig_p_Ochiai_df), sig_cor_select2)
View(test)

View(test[test$cor>0.9,])
View(test[test$cor>0.7 & test$Ochiaiindex > 0.2,])




# detect if gene can be found in CCLE ------------------------------------------
sig_select <- read.table("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/07_plotting_v2/01_metaAnalysis/table/metaSelected_immunesigAnnot_subtypedetail_filter.txt", sep = '\t', header=T)

genestat_selected3 <- genestat_selected2[is.element(genestat_selected2$immune_sig, sig_select$immune_sig), ]

genestat_selected3 %>% na.omit() %>% 
  group_by(variable,immune_sig) %>% 
  summarise(mean_value=mean(value),
            sd_value=sd(value),
            min_value=min(value),
            max_value=max(value)) -> df1


df1 %>% 
  group_by(immune_sig) %>% 
  mutate(new_col=cumsum(mean_value)) -> df_ORI

df_ORI$variable<-factor(df_ORI$variable, levels = unique(df_ORI$variable))
df_ORI = as.data.frame(df_ORI)
df_ORI[df_ORI$variable == 'gene_in_OriginalSig_expressedinCCLE_rate',]
df_ORI_not_expressedinCCLE_rate = df_ORI[df_ORI$variable == 'gene_in_OriginalSig_not_expressedinCCLE_rate',]
nrow(df_ORI_not_expressedinCCLE_rate)
nrow(df_ORI_not_expressedinCCLE_rate[df_ORI_not_expressedinCCLE_rate$mean_value > 0.5,])
nrow(df_ORI_not_expressedinCCLE_rate[df_ORI_not_expressedinCCLE_rate$mean_value > 0.3,])


# # single cohort

# sigs_ds = lapply(as.list(unique(res_df$dataset)), function(ds){
#      # ds = 'GSE145996'
#      print(ds)
#      res_df_ds = res_df[res_df$dataset == ds, ]

#      res_df_annot <- res_df_ds[c(1:5)]
#      res_df_metaselect <- res_df_ds[is.element(names(res_df_ds), immunesigAnnot$immune_sig)]


#      sig_cor = cor(res_df_metaselect, method='spearman')
#      immunesigAnnot = immunesigAnnot[order(immunesigAnnot$SubType), ]
#      sig_cor = sig_cor[immunesigAnnot$immune_sig, immunesigAnnot$immune_sig]
#      # pheatmap(sig_cor, cluster_row = F, annotation_row = immunesigAnnot[c( 'SubType','Type','flag')])


#      immunesigAnnot_s = immunesigAnnot[is.element(immunesigAnnot$flag, "sensitive"), ]
#      immunesigAnnot_s = immunesigAnnot_s[order(immunesigAnnot_s$SubType), ]
#      sig_cor_s = sig_cor[immunesigAnnot_s$immune_sig, immunesigAnnot_s$immune_sig]

#      sigs_s = lapply(as.list(seq(length(unique(immunesigAnnot_s$SubType)))), function(i){
#           # i = 4
#           subtype = unique(immunesigAnnot_s$SubType)[i]
#           print(subtype)
#           # print(tmp)
#           if(length(immunesigAnnot_s[immunesigAnnot_s$SubType == subtype,]$immune_sig) == 1){
#           sigs = immunesigAnnot_s[is.element(immunesigAnnot_s$immune_sig, immunesigAnnot_s[immunesigAnnot_s$SubType == subtype,]$immune_sig), ]
#           return(sigs)

#           }else if(subtype == 'Miscellaneous'){
#           sigs = immunesigAnnot_s[is.element(immunesigAnnot_s$immune_sig, immunesigAnnot_s[immunesigAnnot_s$SubType == subtype,]$immune_sig), ]
#           return(sigs)

#           }else {
#           tmp = sig_cor_s[immunesigAnnot_s[immunesigAnnot_s$SubType == subtype,]$immune_sig, immunesigAnnot_s[immunesigAnnot_s$SubType == subtype,]$immune_sig]
#           # pheatmap(tmp)

#           sigs = immunesigAnnot_s[is.element(immunesigAnnot_s$immune_sig, colnames(tmp)), ]
#           sigs = sigs[order(sigs$OR, decreasing = T), ]
#           # print(sigs)
#           return(head(sigs,1))
#           }
#      })

#      sigs_s = do.call(rbind, sigs_s)
#      # pheatmap(sig_cor[sigs_s$immune_sig, sigs_s$immune_sig], annotation_row = immunesigAnnot[c( 'SubType','Type','flag')])





#      immunesigAnnot_r = immunesigAnnot[is.element(immunesigAnnot$flag, "resistant"), ]
#      immunesigAnnot_r = immunesigAnnot_r[order(immunesigAnnot_r$SubType), ]
#      sig_cor_r = sig_cor[immunesigAnnot_r$immune_sig, immunesigAnnot_r$immune_sig]

#      sigs_r = lapply(as.list(seq(length(unique(immunesigAnnot_r$SubType)))), function(i){
#           # i = 5
#           subtype = unique(immunesigAnnot_r$SubType)[i]
#           print(subtype)
#           # print(tmp)
#           if(length(immunesigAnnot_r[immunesigAnnot_r$SubType == subtype,]$immune_sig) == 1){
#           sigs = immunesigAnnot_r[is.element(immunesigAnnot_r$immune_sig, immunesigAnnot_r[immunesigAnnot_r$SubType == subtype,]$immune_sig), ]
#           return(sigs)

#           }else if(subtype == 'Miscellaneous'){
#           sigs = immunesigAnnot_r[is.element(immunesigAnnot_r$immune_sig, immunesigAnnot_r[immunesigAnnot_r$SubType == subtype,]$immune_sig), ]
#           return(sigs)

#           }else{
#           tmp = sig_cor_r[immunesigAnnot_r[immunesigAnnot_r$SubType == subtype,]$immune_sig, immunesigAnnot_r[immunesigAnnot_r$SubType == subtype,]$immune_sig]
#           # pheatmap(tmp)

#           sigs = immunesigAnnot_r[is.element(immunesigAnnot_r$immune_sig, colnames(tmp)), ]
#           sigs = sigs[order(sigs$OR, decreasing = T), ]
#           # print(sigs)
#           return(head(sigs,1))
#           }
#      })

#      sigs_r = do.call(rbind, sigs_r)
#      # pheatmap(sig_cor[sigs_r$immune_sig, sigs_r$immune_sig], annotation_row = immunesigAnnot[c( 'SubType','Type','flag')])


#      sig_select <- rbind(sigs_s,sigs_r)
#      return(sig_select)
# })

# sigs_dsall = do.call(rbind, sigs_ds)
# table(sigs_dsall$immune_sig)

# unique( sigs_dsall$immune_sig)
# #     load("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/05_immuneSig_ICBresponse/results_r0.4_ITSssgsea_CombatRef/finalmodel_17datasets_ssgseanorm/metaAnalysis/en/elasticnet_model_pred.Rdata")

# #     en_model_para = m_en$learner.model$glmnet.fit
# #     m_en_imp <- data.frame(en_model_para$beta[,which(en_model_para$lambda == m_en$learner.model$lambda.min)])

# #     names(m_en_imp) <- "weight"
# #     m_en_imp$immune_sig = row.names(m_en_imp)
# #     immunesig = merge(immunesig, m_en_imp, by = "immune_sig")
# #     immunesig =immunesig[immunesig$weight != 0 ,]


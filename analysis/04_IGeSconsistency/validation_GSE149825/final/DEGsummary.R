rm(list=ls())
# Mapping drug DEGs to correlation results
# work_path <- "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/"
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
 
source("06_results_summaryandplotting/scripts/plot_for_drug_function.R")
source("06_results_summaryandplotting/scripts/immunesig_selection_function.R")

# ----------------------------------------------------------------------------
# drug-induced immune signatures NES
cancer = "SKCM"
tumor_subtype = "Metastatic"
purity_method = "TUMERIC"
dataset=70138
datatype="allgenes"
num_gene = 200


savepath = paste0("06_results_summaryandplotting/results_", purity_method,"/")
dir.create(savepath)



# ITSselected <- read.csv("06_results_summaryandplotting/data/immune_sig_selected_new/selected_ITS_0.4.txt", sep = '\t')
meta_res <- read.csv("05_immuneSig_ICBresponse/results/meta_results/res_metaRes.csv")
meta_res_selected <- meta_res[meta_res$Meta_Pval < 0.05, ]
sig_s = meta_res_selected[meta_res_selected$OR > 1,]
sig_r = meta_res_selected[meta_res_selected$OR < 1,]
sig_s = data.frame( immune_sig = sig_s$X, flag = "sensitive", setindex = '01_sensitive')
sig_r = data.frame( immune_sig = sig_r$X, flag = "resistant", setindex = '02_resistant')

load("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/05_immuneSig_ICBresponse/results_r0.4_ITSssgsea_CombatRef/finalmodel_17datasets_ssgseanorm/metaAnalysis_filter/en/elasticnet_model_pred.Rdata")

en_model_para = m_en$learner.model$glmnet.fit
m_en_imp <- data.frame(en_model_para$beta[,which(en_model_para$lambda == m_en$learner.model$lambda.min)])

names(m_en_imp) <- "weight"
m_en_imp$immune_sig = row.names(m_en_imp)


immunesig = rbind(sig_s, sig_r)
immunesig = merge(immunesig, m_en_imp, by = "immune_sig")
immunesig =immunesig[immunesig$weight != 0 ,]

load(paste0("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/02_tumorGene_immuneSig_correlation/results_selected/ITS_oe_original_spearman_TUMERIC_r0.4/for_GSVA/pcor/pcor_spearman_",cancer,"_",tumor_subtype,"_positive_200.Rdata"))
resultpath = paste0("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/result_GSE149825/")
DEGlist <- read.csv(paste0(resultpath, 'GSE149825_DEGanalysis.csv'))
DEGlist <- DEGlist[DEGlist$padj <0.1,]
DEGup <- DEGlist[DEGlist$log2FoldChange >0,]$ID
DEGdown <- DEGlist[DEGlist$log2FoldChange <0,]$ID
genelist_p_s_selected = genelist_p[immunesig[immunesig$flag == 'sensitive',]$immune_sig]
genelist_p_r_selected = genelist_p[immunesig[immunesig$flag == 'resistant',]$immune_sig]

num_sIGeS = t(as.data.frame(lapply(genelist_p_s_selected, length)))
num_intersect_sIGeS = t(as.data.frame(lapply(genelist_p_s_selected, function(x){length(intersect(x, toupper(DEGup)))})))
res_sIGeS = as.data.frame(cbind(num_intersect_sIGeS, num_sIGeS))
names(res_sIGeS) = c('intersect_DEGup_sIGeS', 'num_sIGeS')
res_sIGeS$numDEG = length(DEGup)
res_sIGeS$flag = 'sensitive'

num_rIGeS = t(as.data.frame(lapply(genelist_p_r_selected, length)))
num_intersect_rIGeS = t(as.data.frame(lapply(genelist_p_r_selected, function(x){length(intersect(x, toupper(DEGdown)))})))
res_rIGeS = as.data.frame(cbind(num_intersect_rIGeS, num_rIGeS))
names(res_rIGeS) = c('intersect_DEGdown_rIGeS', 'num_rIGeS')
res_sIGeS$numDEG = length(DEGdown)
res_rIGeS$flag = 'resistant'

write.csv(res_sIGeS, paste0(resultpath,'res_sIGeS.csv'), row.names=F, quote=F)
write.csv(res_rIGeS, paste0(resultpath,'res_rIGeS.csv'), row.names=F, quote=F)
write.csv(DEGlist, paste0(resultpath,'DEGlist_adjP0.1.csv'), row.names=F, quote=F)

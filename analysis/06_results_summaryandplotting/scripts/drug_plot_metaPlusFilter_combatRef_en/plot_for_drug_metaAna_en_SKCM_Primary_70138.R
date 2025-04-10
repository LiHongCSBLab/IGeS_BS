# rm(list=ls())
# Mapping drug DEGs to correlation results
# work_path <- "D:/veronica_files/lab/Projects/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/"
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
# library(readxl)
library(stringr)
library(data.table)


library(getopt)

command=matrix(c( 
  "weight_type",     "w", 1, "character", "Define weight type: model_weight or meta_weight",
  "type",            "t", 1, "character", "Define if weighted: unweighted, weighted or binary_weighted",
  "drugrank_method", "m", 1, "character", "Define integrated weight method: glmnet_weight, simple or sum_mean",
  "mergeITS",        "i", 1, "logical", "Define if merge ITS",
  "mergeITS_type",   "s", 1, "character", "Define the type for merge ITS: Type or Index",
  "ACAT_Pval",       "p", 1, "double",    "Define integrated weight method: glmnet_weight, simple or sum_mean",
  "help",            "h", 0, "logical",   "help file"),
  byrow=T,ncol=5)

args=getopt(command)
# -w model_weight -t unweighted -m simple

weight_type = args$weight_type
type = args$type
drugrank_method = args$drugrank_method
ACAT_Pval=args$ACAT_Pval
mergeITS = args$mergeITS
mergeITS_type = args$mergeITS_type


source("06_results_summaryandplotting/scripts/plot_for_drug_function_new.R")
source("06_results_summaryandplotting/scripts/immunesig_selection_function.R")
# source("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/scripts/drugITGSscore_lincs.R")

cancer = "SKCM"
tumor_subtype = "Primary"
purity_method = "CPE"
dataset=70138
datatype="allgenes"
num_gene = 200
cor_method = "spearman"
cluster_method = "mclust"
similarity_method = "Otsuka_Ochiai" # "jaccard_index",
num_gene = 200
s1 = 0.4
ITSproof = 2
genevote = 0.5
model = "m2"
# ACAT_Pval=0.05




immunesig <- read.csv("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/05_immuneSig_ICBresponse/results/meta_results/res_metaRes.csv")
immunesig <- immunesig[immunesig$Meta_Pval < 0.05, ]
names(immunesig) <- c("immune_sig", names(immunesig)[-1])
immunesig$flag <- "sensitive"
immunesig[immunesig$OR < 1,]$flag <- "resistant"
immunesig$setindex = "01_sensitive"
immunesig[immunesig$flag == "resistant", ]$setindex = "02_resistant"

load("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/05_immuneSig_ICBresponse/results_r0.4_ITSssgsea_CombatRef/finalmodel_17datasets_ssgseanorm/metaAnalysis_filter/en/elasticnet_model_pred.Rdata")

en_model_para = m_en$learner.model$glmnet.fit
m_en_imp <- data.frame(en_model_para$beta[,which(en_model_para$lambda == m_en$learner.model$lambda.min)])

names(m_en_imp) <- "weight"
m_en_imp$immune_sig = row.names(m_en_imp)
immunesig = merge(immunesig, m_en_imp, by = "immune_sig")
immunesig =immunesig[immunesig$weight != 0 ,]
# write.table(immunesig, '/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/05_immuneSig_ICBresponse/results_r0.4_ITSssgsea_CombatRef/finalmodel_17datasets_ssgseanorm/metaAnalysis_filter/en/elasticnet_model_weight.txt', quote=F, sep='\t')

dir.create("06_results_summaryandplotting/results_meta_filter_ITSssgsea_CombatRef_17datasets/")
savepath = paste0("06_results_summaryandplotting/results_meta_filter_ITSssgsea_CombatRef_17datasets/", purity_method,"_", model,"_meta_en_ACATP",ACAT_Pval,"/") #
dir.create(savepath)

drugITSscore_plot(cancer = cancer,
                   tumor_subtype = tumor_subtype,
                   purity_method = purity_method,
                   dataset = dataset,
                   datatype = datatype,
                   num_gene = num_gene,
                   immunesig = immunesig,
                   model = model,
                   weight_type = weight_type,
                   type = type,
                   drugrank_method = drugrank_method,
                   mergeITS = mergeITS,
                   mergeITS_type = mergeITS_type,
                   ACAT_Pval = ACAT_Pval,
                   ifProfile = TRUE, 
                   work_path = work_path,
                   savepath = savepath)


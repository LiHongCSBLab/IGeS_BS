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



# ----------------------------------------------------------------------------
# LOAD ITS gene sets
purity_method = "TUMERIC"

cancer = "SKCM"
tumor_subtype = "Metastatic"
datatype="allgenes"
num_gene = 200
purity_method = "TUMERIC"

cancerlist <- unique(gsub("spearman_","", 
                          gsub("_Metastatic_positive_200.Rdata","",
                               gsub("_Metastatic_negative_200.Rdata","",
                                    gsub("_Primary_positive_200.Rdata","",
                                         gsub("_Primary_negative_200.Rdata","",
                                              list.files("03_drug_immuneSig_enrichment/data/gene_immune_sig_TUMERIC_spearman/for_GSVA/pcor")))))))
cancerlist <- paste0(cancerlist, "_Primary")
cancerlist <- c("SKCM_Metastatic", cancerlist)

ITS_p_set <- lapply(as.list(cancerlist), function(cancer){
  num_gene = 200
  r = 0.4
  
  if(cancer == "SKCM_Primary"){
    purity_method = "CPE"
    ITSpath = paste0("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/02_tumorGene_immuneSig_correlation/results_selected/ITS_oe_original_spearman_",purity_method,"_r",r,"/for_GSVA/pcor/")
    load(paste0(ITSpath, "pcor_spearman_",cancer, "_positive_200.Rdata"))
  }else {
    purity_method = "TUMERIC"
    ITSpath = paste0("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/02_tumorGene_immuneSig_correlation/results_selected/ITS_oe_original_spearman_",purity_method,"_r",r,"/for_GSVA/pcor/")
    load(paste0(ITSpath, "pcor_spearman_",cancer, "_positive_200.Rdata"))
  }
  
  genelist_p_ori = genelist_p[meta_selected$immune_sig]
  
  genelist_p <- lapply(genelist_p_ori, function(x){
    intersect(x, geneExprInCCLE_filtered$gene_name)
  })
  names(genelist_p) = names(genelist_p_ori)
  return(genelist_p)
})
names(ITS_p_set) <- cancerlist


ITS_n_set <- lapply(as.list(cancerlist), function(cancer){
  num_gene = 200
  r = 0.4
  
  if(cancer == "SKCM_Primary"){
    purity_method = "CPE"
    ITSpath = paste0("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/02_tumorGene_immuneSig_correlation/results_selected/ITS_oe_original_spearman_",purity_method,"_r",r,"/for_GSVA/pcor/")
    load(paste0(ITSpath, "pcor_spearman_",cancer, "_negative_200.Rdata"))
  }else {
    purity_method = "TUMERIC"
    ITSpath = paste0("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/02_tumorGene_immuneSig_correlation/results_selected/ITS_oe_original_spearman_",purity_method,"_r",r,"/for_GSVA/pcor/")
    load(paste0(ITSpath, "pcor_spearman_",cancer, "_negative_200.Rdata"))
  }
  
  genelist_n_ori = genelist_n[meta_selected$immune_sig]
  genelist_n <- lapply(genelist_n_ori, function(x){
    intersect(x, geneExprInCCLE_filtered$gene_name)
  })
  names(genelist_n) = names(genelist_n_ori)
  return(genelist_n)
})

names(ITS_n_set) <- cancerlist

save(ITS_p_set, ITS_n_set, file = "07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITS_set.Rdata")
# save(ITS_p_set, ITS_n_set, file = "07_plotting_v2/data/ITS_set.Rdata")

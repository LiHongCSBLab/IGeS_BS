# Mapping drug DEGs to correlation results -------------------------------------
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
source("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/scripts/drugITGSscore_lincs.R")
source("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/scripts/drugITGSscore_GEOdataset.R")

# ----------------------------------------------------------------------------
# drug-induced immune signatures NES
# cancer = "LIHC"
tumor_subtype = "Primary"
purity_method = "TUMERIC"
dataset=70138
datatype="allgenes"
num_gene = 200

savepath  = "06_results_summaryandplotting/data/immune_sig_selected_new/"
dir.create(savepath)


immunesig <- immune_sig_filter_v2(auc_threshold = 0.6,
                                  immuSigproof = 2,
                                  ITSproof = 2, # how many datasets support the immune sigture to have prediction power
                                  savepath = savepath)





# ITS_removeredundancy_savepath = paste0("06_results_summaryandplotting/results_", purity_method,"/")
# dir.create(ITS_removeredundancy_savepath)
# cancerlist <- c( "UCEC","LIHC", "PRAD", "PAAD", "BRCA", "READ","LUAD","LUSC",
#                  "OV",  "SKCM","COAD", "DLBC", 
#                  "UCS", "HNSC", "BLCA","CESC","KIRC", "CHOL", "ESCA", "PCPG",
#                  "ACC","LGG", "GBM", "KICH",
#                  "SARC","STAD","TGCT","THYM", "MESO", "KIRP","THCA")

# for(cancer in cancerlist){
  
#   immunesig3 <- immune_sig_finefilter(
#     immunesig = immunesig2,
#     cancer = cancer,
#     tumor_subtype = tumor_subtype,
#     purity_method = purity_method,
#     datatype="allgenes",
#     num_gene = 200,
#     workpath = "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/", 
#     savepath = "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/" )
# }

# cancer="SKCM"
# tumor_subtype = "Metastatic"
# purity_method = "TUMERIC"
# dataset=70138
# datatype="allgenes"
# num_gene = 200


# immunesig3 <- immune_sig_finefilter(
#   immunesig = immunesig2,
#   cancer = cancer,
#   tumor_subtype = tumor_subtype,
#   purity_method = purity_method,
#   datatype="allgenes",
#   num_gene = 200,
#   workpath = "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/", 
#   savepath = "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/" )

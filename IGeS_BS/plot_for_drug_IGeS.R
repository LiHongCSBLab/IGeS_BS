# rm(list=ls())
# Mapping drug DEGs to correlation results
# work_path <- "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/IGeS_BS/IGeS_BS_tool/IGeS_BS/"
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
library(fs)

library(getopt)

command=matrix(c( 
  "dataset",          "a", 1, "character", "Define which dataset to use, 70138 or 92742",
  "cancer",           "c", 1, "character", "Define cancer type, e.g. LIHC",
  "tumor_subtype",    "t", 1, "character", "Primary or Metastatic. Current for Metastatics only support for SKCM",
  "purity_method",    "u", 1, "character", "Define method used for computing tumor purity. TUMERIC or CPE. Current forfor primary SKCM, only CPE method is supported",
  "ifProfile",        "f", 1, "logical",   "if IGeS profiling radar plot should be generated, TRUE or FALSE",
  "filepath",         "d", 1, "character",   "Define the path to IGeS_Profiling",
  "work_path",        "w", 1, "character",   "Define the path to work directory",
  "save_path",        "s", 1, "character",   "Define the path to save results",
  "help",              "h", 0, "logical",   "help file"),
  byrow=T,ncol=5)

args=getopt(command)
# runner:
# Rscript --vanilla plot_for_drug_IGeS.R -a 70138 -c LIHC -t Primary -u TUMERIC -f FALSE \
# -d "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/03_drug_immuneSig_enrichment/" \
# -w "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/IGeS_BS/IGeS_BS_tool/" \
# -s "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/IGeS_BS/IGeS_BS_tool/"> IGeS-BS.log 2>&1&


dataset=args$dataset
cancer = args$cancer
tumor_subtype = args$tumor_subtype
purity_method = args$purity_method
ifProfile = args$purity_method
filepath=args$filepath
work_path=args$work_path
save_path=args$save_path

# Check if optional parameter `ifProfile` was provided
if (!is.null(args$weight_type)) {
weight_type=args$weight_type
} 
if (!is.null(args$type)) {
type=args$type
} 
if (!is.null(args$ACAT_Pval)) {
ACAT_Pval=args$ACAT_Pval
} 

setwd(work_path)

# load functions
func1=path("IGeS_BS", "plot_for_drug_function_new.R")
func2=path("IGeS_BS", "immunesig_selection_function.R")
source(func1)
source(func2)
# source("./IGeS_BS/plot_for_drug_function_new.R")
# source("./IGeS_BS/immunesig_selection_function.R")


# load dependency files

immunesig <- read.csv("./Data/res_metaRes.csv")
immunesig <- immunesig[immunesig$Meta_Pval < 0.05, ]
names(immunesig) <- c("immune_sig", names(immunesig)[-1])
immunesig$flag <- "sensitive"
immunesig[immunesig$OR < 1,]$flag <- "resistant"
immunesig$setindex = "01_sensitive"
immunesig[immunesig$flag == "resistant", ]$setindex = "02_resistant"

load("./Data/IGeS_EN/elasticnet_model_pred.Rdata")
en_model_para = m_en$learner.model$glmnet.fit
m_en_imp <- data.frame(en_model_para$beta[,which(en_model_para$lambda == m_en$learner.model$lambda.min)])

names(m_en_imp) <- "weight"
m_en_imp$immune_sig = row.names(m_en_imp)
immunesig = merge(immunesig, m_en_imp, by = "immune_sig")
immunesig =immunesig[immunesig$weight != 0 ,]



# set up output folder
savepath=paste0(save_path, "./results/")
dir.create(savepath)

# example of  additional parameters
# dataset=70138
# cancer='LIHC'
# tumor_subtype = "Primary"
# purity_method = "TUMERIC"
# weight_type = "model_weight"
# type = "weighted"
# model = "m2"
# drugrank_method = "glmnet_weight"
# ACAT_Pval = 0.05
# datatype = "allgenes"
# num_gene = 200
# work_path2="/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/"03_drug_immuneSig_enrichment/"
filepath =  paste0(filepath, "results_", purity_method,"_meta_oe_original_final/")


drugITSscore_plot(cancer = cancer,
                  tumor_subtype = tumor_subtype,
                  purity_method = purity_method,
                  dataset = dataset,
                  immunesig = immunesig,
                  ifProfile = TRUE, 
                  filepath = filepath, 
                  savepath = savepath)


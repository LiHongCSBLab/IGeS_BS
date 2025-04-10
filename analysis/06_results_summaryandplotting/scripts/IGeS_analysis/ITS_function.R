# Mapping drug DEGs to correlation results
print("start computing")
work_path <- "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/"

setwd(work_path)
options(stringsAsFactors = F)
library(parallel)
library(reshape2)
library(stringr)
library(dplyr)


source("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/03_drug_immuneSig_enrichment/scripts/02_GSEA_geneset_preparation_function.R")
source("06_results_summaryandplotting/scripts/immunesig_selection_function.R")

# tumor_purity_methods = c("IHC", "CPE", "TUMERIC", "ABSOLUTE")
# for(tumor_purity_method in tumor_purity_methods){
tumor_purity_method= "TUMERIC"
cor_method = "spearman"
immunesig_path = "06_results_summaryandplotting/data/immune_sig_selected/"
savepath_allresult = paste0("02_tumorGene_immuneSig_correlation/results_selected/immunecell_TCGAgenes_result_", cor_method,"_", tumor_purity_method, "/")
dir.create(savepath_allresult)
savepath = paste0("03_drug_immuneSig_enrichment/data/gene_immune_sig_", tumor_purity_method, "_",cor_method,"/")


cancer = "SKCM"
sample_type = "Metastatic"
immunesig_confident_filter = 1.5
immunesig <- read.csv(paste0(immunesig_path,"/", cancer,"_",sample_type,".csv"))
immunesig <- immunesig[which(immunesig$confident >= immunesig_confident_filter ), ]

load(paste0("06_results_summaryandplotting/ITS_function/gobp_", cancer, "_", sample_type,"_p.Rdata"))
ITS_p_s <- ITSenrichRes_p[immunesig[immunesig$setindex == 1,]$immune_sig]
ITS_p_s <- do.call(rbind, ITS_p_s)

ITS_p_r <- ITSenrichRes_p[immunesig[immunesig$setindex == 2,]$immune_sig]

ITS_p_r <- do.call(rbind, ITS_p_r)
load(file= paste0("06_results_summaryandplotting/ITS_function/gobp_", 
                    cancer, "_", sample_type,"_n.Rdata"))
ITS_n_s <- ITSenrichRes_n[immunesig[immunesig$setindex == 1,]$immune_sig]
ITS_n_s <- do.call(rbind, ITS_n_s)
ITS_n_r <- ITSenrichRes_n[immunesig[immunesig$setindex == 2,]$immune_sig]
ITS_n_r <- do.call(rbind, ITS_n_r)

View(ITS_p_s[ITS_p_s$p.adjust < 0.05, ])
View(ITS_n_s[ITS_n_s$p.adjust < 0.05, ])

View(ITS_p_r[ITS_p_r$p.adjust < 0.05, ])
View(ITS_n_r[ITS_n_r$p.adjust < 0.05, ])


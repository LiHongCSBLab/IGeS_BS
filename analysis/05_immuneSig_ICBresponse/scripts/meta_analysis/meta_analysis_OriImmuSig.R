# rm(list=ls())
work_path <- "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/"
setwd(work_path)
options(stringsAsFactors = F)


library(reshape2)
library(dplyr)
library(stringr)
library(data.table)
library(ggplot2)
library(ggplotify)
library(grid)
library(cowplot)
library(patchwork)
require(GGally)
# library(plot3D)
# library(fmsb)

source("06_results_summaryandplotting/scripts/immunesig_selection_function.R")
source("03_drug_immuneSig_enrichment/scripts/merge_p.R")
source("05_immuneSig_ICBresponse/scripts/meta_analysis/meta_analysis_function.R")

dataset_annot <- read.csv("05_immuneSig_ICBresponse/dataset_annot.csv")
dataset_annot <- unique(dataset_annot[c("dataset", "cancer")])
dataset_annot$treatment_type <- c("aPD1", "aPD1", "aPD1", "aPD1", "aPDL1", 
                                  "aPD1", "aPD1", "aPD1_aPDL1", "aPD1","aCTLA4",
                                  "aPD1", "aPD1")
dataset_annot$dataset <- gsub("_outcome", "", dataset_annot$dataset)
dataset_annot$dataset <- gsub("_bladder_immunotherapy", "_dataset", dataset_annot$dataset)
dataset_annot$dataset <- gsub("Braun_dataset", "Braun_dataset_PD1", dataset_annot$dataset)

# ---------------------------------------------------------------------
icb_res_path = "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/05_immuneSig_ICBresponse/results/"
files = list.files(paste0(icb_res_path, "meta_analysis_files/"))
files = files[-grep("Lauss", files)]

resall = lapply(as.list(files), function(x){
  df = read.csv(paste0(icb_res_path, "meta_analysis_files/",x))
  df$dataset = gsub("_test_result.csv","",x)
  return(df)
})

resall = unique(do.call(rbind,resall))


xcell_cellfamily <- read.csv(paste0(work_path, "02_tumorGene_immuneSig_correlation/data/immune_signatures/originalSignaturesGenes/immune_sig_original_files/XCELL/cells_families.txt"), sep = '\t')
xcell_type1 <- xcell_cellfamily[is.element(xcell_cellfamily$Group, 
                                           c("B-cells",  "CD4+ memory T-cells", "CD4+ T-cells", "CD8+ T-cells",
                                             "DC",  "Endothelial cells", "Epithelial cells", "Fibroblasts", "Macrophages")), ]
xcell_type2 <- xcell_cellfamily[is.element(xcell_cellfamily$Group, "Parent"), ]
xcell_type3 <- xcell_type2[is.element(xcell_type2$Type, c("Lymphoid", "Myeloid")), ]
xcell_type4 <- xcell_type2[is.element(xcell_type2$Cells, c("Endothelial cells",  "Epithelial cells", "Fibroblasts")), ]
xcell_select <- rbind(xcell_type1, xcell_type3, xcell_type4)
xcell_remove <- xcell_cellfamily[!is.element(xcell_cellfamily$Cells, xcell_select$Cells),]
xcell_remove <- paste0(xcell_remove$Cells, "_xcell")
xcell_remove <- gsub(" ", ".", xcell_remove)

resall_immunesigname <- data.frame(immune_sig = tolower(resall$immunesig), immunesig = resall$immunesig)
resall_immunesigname = immunesig_name_unify(resall_immunesigname)
resall <- inner_join(resall_immunesigname, resall)
resall <- inner_join(dataset_annot, resall)
resall <- resall[-which(is.element(resall$immune_sig, tolower(xcell_remove))),]
resall <- unique(resall)
resall <- resall[-which(names(resall) == "X")]
write.csv(resall, paste0(icb_res_path, "meta_results/immunesig_icbResponse_all_results.csv"), quote = F)


source("05_immuneSig_ICBresponse/scripts/meta_analysis/meta_analysis_function.R")

res_all_combined_p <- ICB_combined_Pvalue(resall, 
                                        cancer = NULL,
                                        treatment = NULL,
                                        save_path = paste0(icb_res_path, "meta_results/"))
for(cancer in unique(dataset_annot$cancer)){
  res_skcm_combined_p <- ICB_combined_Pvalue(resall, 
                                          cancer =  cancer,
                                          treatment = NULL,
                                          save_path = paste0(icb_res_path, "meta_results/"))

}

# resall2 <- resall[resall$dataset != '07_IMvigor210_dataset' & resall$dataset != 'Braun_dataset_PD1',]
res_all_metagen <- ICB_metaAnalysis(resall = resall, 
                                    cancer = NULL,
                                    treatment = NULL,
                                    save_path = paste0(icb_res_path, "meta_results/"))
res_all_metagen[res_all_metagen$Meta_Pval < 0.05,]

for(cancer in unique(dataset_annot$cancer)){
res_skcm_metagen <- ICB_metaAnalysis(resall, 
                                    cancer = cancer,
                                    treatment = NULL,
                                    save_path =  paste0(icb_res_path, "meta_results/"))

}
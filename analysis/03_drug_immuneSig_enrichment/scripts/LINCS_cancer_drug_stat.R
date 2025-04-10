work_path <- "/picb/bigdata/project/FengFYM/drug_MoA_reanalysis/"
setwd(work_path)
options(stringsAsFactors=F)
library(ggplot2)
library(getopt)
library(prada)
library(rhdf5)
library(limma)
library(dplyr)

source("l1ktools-master/R/cmap/io.R")
options(stringsAsFactors = F)
options(quote = NULL)

cancer = "LIHC" 
dataset = "92742"
result_path <- "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/03_drug_immuneSig_enrichment/data/"

drug_index <- read.table(paste0(result_path, 'drug_index_LINCS.txt'), sep = '\t', header = T)
inst_GSE70138 <- read.csv("LINCS_data/GSE70138/file/GSE70138_Broad_LINCS_inst_info.txt", header = T, sep = "\t")
inst_GSE70138 <- inst_GSE70138[which(is.element(inst_GSE70138$pert_type, c("ctl_vehicle", "trt_cp")))]
# assign index to drugs
inst_GSE70138 <- inner_join(drug_index, inst_GSE70138)

cell_annot_selected <- read.table(paste0(result_path, 'cell_annot_selected_70138.txt'), sep = '\t',header = T)
inst_GSE70138 <- inner_join(cell_annot_selected, inst_GSE70138)

inst_GSE92742 <- read.csv("LINCS_data/GSE92742/file/GSE92742_Broad_LINCS_inst_info.txt", header=T, sep="\t")
inst_GSE92742[inst_GSE92742$pert_time == "24|4",]$pert_time <- "24"
inst_GSE92742[inst_GSE92742$pert_time == "6|4",]$pert_time <- "6"
inst_GSE92742 = inst_GSE92742[which(is.element(inst_GSE92742$pert_type, c("ctl_vehicle", "trt_cp"))), ]
inst_GSE92742 <- inner_join(drug_index, inst_GSE92742)
inst_GSE92742 = inst_GSE92742[-which(is.element(inst_GSE92742$pert_iname, c("PBS","H2O"))), ]

cell_annot_selected <- read.table(paste0(result_path, 'cell_annot_selected_92742.txt'), sep = '\t',header = T)
inst_GSE92742 <- inner_join(cell_annot_selected, inst_GSE92742)

drug_stat_70138 = data.frame(cancer = unique(c(
  "BRCA", "COAD", "LIHC", "LUAD", "PAAD", "PRAD", "READ", "SKCM", "BRCA", "COAD", "DLBC",
  "LAML", "LIHC", "LUAD", "LUSC", "OV", "PRAD", "READ", "SKCM", "STAD", "UCEC"
)), dataset = "70138", num_drug = 0)

drug_stat_92742 = data.frame(cancer = unique(c(
  "BRCA", "COAD", "LIHC", "LUAD", "PAAD", "PRAD", "READ", "SKCM", "BRCA", "COAD", "DLBC",
  "LAML", "LIHC", "LUAD", "LUSC", "OV", "PRAD", "READ", "SKCM", "STAD", "UCEC"
)), dataset = '92742', num_drug = 0)

for(dataset in c('70138', '92742')){
  if(dataset == '70138'){
    for(cancer in c("BRCA","COAD","LIHC","LUAD","PAAD","PRAD","READ","SKCM")){
      
      # dataset <- "70138"
      inst_GSE70138_cancer <- inst_GSE70138[inst_GSE70138$TCGA_cancer_type == cancer, ]
      
      geneinfo <- read.csv("LINCS_data/GSE70138/file/GSE92742_Broad_LINCS_gene_info.txt", sep = "\t", header = T)
      geneinfo_sig <- geneinfo[which(geneinfo$pr_is_lm == 1), ]
      
      drug_stat_70138[drug_stat_70138$cancer == cancer, ]$num_drug =
        length(unique(inst_GSE70138_cancer[-which(inst_GSE70138_cancer$pert_iname == "DMSO"), ]$drug_index))
      
    }
  } else if (dataset == '92742') {
    
    for(cancer in c("BRCA","COAD","DLBC","LAML","LIHC","LUAD","LUSC","OV","PRAD","READ","SKCM","STAD","UCEC")){
      
      inst_GSE92742_cancer <- inst_GSE92742[inst_GSE92742$TCGA_cancer_type == cancer, ]
      drug_stat_92742[drug_stat_92742$cancer == cancer,]$num_drug =
        length(unique(inst_GSE92742_cancer[-which(inst_GSE92742_cancer$pert_iname=='DMSO'),]$drug_index))
    }
  }
}

write.csv(drug_stat_70138, "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/03_drug_immuneSig_enrichment/data/drug_stat_70138.csv", quote = F, row.names = F)
write.csv(drug_stat_92742, '/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/03_drug_immuneSig_enrichment/data/drug_stat_92742.csv',quote = F, row.names=F)


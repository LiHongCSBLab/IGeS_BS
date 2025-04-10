# Tumor genes - infiltrated immune cell fraction
workpath <- "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/"
setwd(workpath)
library(parallel)
library(psych)
library(ppcor)
library(reshape2)
options(stringsAsFactors = F)

source("02_tumorGene_immuneSig_correlation/scripts/functions/02_immuneSig_generation_function.R")
source("02_tumorGene_immuneSig_correlation/scripts/functions/03_cor_gene_immuneSig_function.R")



# load immunesig files
immunesig_path <- "02_tumorGene_immuneSig_correlation/data/immune_signatures/"
immunesig_file <- "02_tumorGene_immuneSig_correlation/data/immune_signatures/immune_sigatures.txt"

patient_ID = read.csv("02_tumorGene_immuneSig_correlation/data/TCGA_patient_ID_tumor_type_selected.csv")

# load in integrated immune signature file
if(! file.exists(immunesig_file)){
  immune_sigatures <- immunesig_generator(immunesig_path)
}else{
  immune_sigatures <-  read.table(immunesig_file, sep = "\t", header = T)
}
dim(immune_sigatures)
immune_sigatures <- immune_sigatures[-which(is.element(names(immune_sigatures), "type"))]


# limit to protein_coding_genes ------------------------------------------------
protein_coding_genes <- read.table("03_drug_immuneSig_enrichment/data/protein_coding_genes_hg38.txt", sep = '\t', header = T)


# run correlation calculation in parallel mode ---------------------------------
project_ids <- c( "UCEC","LIHC", "PRAD", "PAAD", "BRCA", "READ","LUAD","LUSC",
                  "OV",  "SKCM","COAD", "DLBC", 
                  "UCS", "HNSC", "BLCA","CESC","KIRC", "CHOL", "ESCA", "PCPG",
                  "ACC","LGG", "GBM", "KICH",
                  "SARC","STAD","TGCT","THYM", "MESO", "KIRP","THCA") # "UVM", "LAML"

# spearman correlation with purity from TUMERIC
cancer_index <- seq(1, length(project_ids))

print("start parallel running for spearman correlation. ")
cores <- 20 # length(project_ids)
clust <- makeCluster(cores) # initialized clusters
clusterExport(clust, ls())  #load global parameters into each cluster
cancername <- parLapply(clust,
                        cancer_index,
                        immunecell_TCGAgene_cor_parallel_spearman_TUMERIC)
stopCluster(clust) 

cancer_index <- seq(1, length(project_ids))
print("start parallel running for TUMERIC purity pearson correlation. ")
cores <- 20 # length(project_ids)
clust <- makeCluster(cores) # initialized clusters
clusterExport(clust, ls())  #load global parameters into each cluster
cancername <- parLapply(clust,
                       cancer_index,
                       immunecell_TCGAgene_cor_parallel_TUMERIC)
stopCluster(clust) # closed clusters











# ------------------------------------------------------------------------------
project_ids <- c("UVM", "LAML")

# spearman correlation with purity from TUMERIC
cancer_index <- seq(1, length(project_ids))

print("start parallel running for spearman correlation. ")
cores <- 2 # length(project_ids)
clust <- makeCluster(cores) # initialized clusters
clusterExport(clust, ls())  #load global parameters into each cluster
cancername <- parLapply(clust,
                        cancer_index,
                        immunecell_TCGAgene_cor_parallel_spearman_TUMERIC)
stopCluster(clust) 

cancer_index <- seq(1, length(project_ids))
print("start parallel running for TUMERIC purity pearson correlation. ")
cores <- 2 # length(project_ids)
clust <- makeCluster(cores) # initialized clusters
clusterExport(clust, ls())  #load global parameters into each cluster
cancername <- parLapply(clust,
                       cancer_index,
                       immunecell_TCGAgene_cor_parallel_TUMERIC)
stopCluster(clust) # closed clusters


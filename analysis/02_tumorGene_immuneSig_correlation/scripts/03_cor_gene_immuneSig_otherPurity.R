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


  
# tcgaRNApath <- "/picb/bigdata/project/FengFYM/immune_targeted_combo/TCGAbiolinks_download_normalization/FPKM/"
# tcgaRNApath <- "/picb/bigdata/project/FengFYM/UCEC_Xena/EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena"
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


project_ids <- c( "UCEC","LIHC", "PRAD", "PAAD", "BRCA", "READ","LUAD","LUSC",
                  "OV",  "SKCM","COAD", "DLBC", "STAD",
                  "UCS", "HNSC", "BLCA","CESC","KIRC", "CHOL", "ESCA", "PCPG",
                  "ACC","LGG", "GBM", "KICH", "SARC","TGCT","THYM", "MESO", "KIRP","THCA") # "UVM", 
cancer_index <- seq(1, length(project_ids))
 


print("start parallel running for CPE purity spearman correlation. ")
cores <- 10 # length(project_ids)
clust <- makeCluster(cores) # initialized clusters
clusterExport(clust, ls())  #load global parameters into each cluster
cancername <- parLapply(clust,
                        cancer_index,
                        immunecell_TCGAgene_cor_parallel_spearman_CPE)
stopCluster(clust)

# cancer_index <- seq(1, length(project_ids))
# print("start parallel running for CPE purity pearson correlation. ")
# cores <- 18 # length(project_ids)
# clust <- makeCluster(cores) # initialized clusters
# clusterExport(clust, ls())  #load global parameters into each cluster
# cancername <- parLapply(clust,
#                         cancer_index,
#                         immunecell_TCGAgene_cor_parallel_CPE)
# stopCluster(clust)


# # timestart <- Sys.time()
print("start parallel running for IHC purity Spearman correlation. ")
# detectCores(logical = F)
cores <- 18 # length(project_ids)
clust <- makeCluster(cores) # initialized clusters
# clusterEvalQ(clust, library(SPEI)) 
clusterExport(clust, ls())  #load global parameters into each cluster
# # run function in a parallel manner

cancername <- parLapply(clust,
                        cancer_index,
                        immunecell_TCGAgene_cor_parallel_spearman_IHC)
stopCluster(clust) 


# print("start parallel running for IHC purity pearson correlation. ")
# # detectCores(logical = F)
# cores <- 18 # length(project_ids)
# clust <- makeCluster(cores) # initialized clusters
# # clusterEvalQ(clust, library(SPEI)) 
# clusterExport(clust, ls())  #load global parameters into each cluster
# # # run function in a parallel manner
# cancername <- parLapply(clust,
#                        cancer_index,
#                        immunecell_TCGAgene_cor_parallel_IHC)

# stopCluster(clust) # closed clusters
              
# # timeend <- Sys.time()
# # runningtime <- timeend - timestart
# # print(runningtime)


# ------------------------------------------------------------------------------
project_ids <- c("UVM", "LAML")
cancer_index <- seq(1, length(project_ids))
 
# # timestart <- Sys.time()
print("start parallel running for IHC purity pearson correlation. ")
# detectCores(logical = F)
cores <- 2 # length(project_ids)
clust <- makeCluster(cores) # initialized clusters
# clusterEvalQ(clust, library(SPEI)) 
clusterExport(clust, ls())  #load global parameters into each cluster
# # run function in a parallel manner

cancername <- parLapply(clust,
                        cancer_index,
                        immunecell_TCGAgene_cor_parallel_spearman_IHC)
stopCluster(clust) 


# print("start parallel running for IHC purity pearson correlation. ")
# # detectCores(logical = F)
# cores <- 2 # length(project_ids)
# clust <- makeCluster(cores) # initialized clusters
# # clusterEvalQ(clust, library(SPEI)) 
# clusterExport(clust, ls())  #load global parameters into each cluster
# # # run function in a parallel manner
# cancername <- parLapply(clust,
#                        cancer_index,
#                        immunecell_TCGAgene_cor_parallel_IHC)

# stopCluster(clust) # closed clusters
              
# # timeend <- Sys.time()
# # runningtime <- timeend - timestart
# # print(runningtime)


print("start parallel running for CPE purity pearson correlation. ")
cores <- 2 # length(project_ids)
clust <- makeCluster(cores) # initialized clusters
clusterExport(clust, ls())  #load global parameters into each cluster
cancername <- parLapply(clust,
                        cancer_index,
                        immunecell_TCGAgene_cor_parallel_spearman_CPE)
stopCluster(clust)

# cancer_index <- seq(1, length(project_ids))
# print("start parallel running for CPE purity pearson correlation. ")
# cores <- 2 # length(project_ids)
# clust <- makeCluster(cores) # initialized clusters
# clusterExport(clust, ls())  #load global parameters into each cluster
# cancername <- parLapply(clust,
#                         cancer_index,
#                         immunecell_TCGAgene_cor_parallel_CPE)
# stopCluster(clust)

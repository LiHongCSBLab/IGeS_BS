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

  setwd(workpath)
  library(psych)
  library(ppcor)
  library(reshape2)
  options(stringsAsFactors = F)
  # exprpath <- "/picb/bigdata/project/FengFYM/UCEC_Xena/EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena"
  # tcgaRNAmatrix <- as.data.frame(data.table::fread(exprpath))
  # colnames(tcgaRNAmatrix) <- substr(colnames(tcgaRNAmatrix),1 , 15)
  # limit to protein_coding_genes ------------------------------------------------
  protein_coding_genes <- read.table("03_drug_immuneSig_enrichment/data/protein_coding_genes_hg38.txt", sep = '\t', header = T)
  
  # sample_types <- read.table( paste0(workpath, "data/TCGA_cancertypes_subtypes.txt"), header=T, sep="\t")
  sample_types <- read.csv("02_tumorGene_immuneSig_correlation/data/TCGA_patient_ID_tumor_type_selected.csv")
  cancertypes <- unique(sample_types$TCGA_cancer_type)
  load("01_tumor_purity/data/processedData/TCGA_tumor_purity_merged.Rdata")
  tumor_purity <-purity_merged
  
  resultpath_TUMERIC2 <- "02_tumorGene_immuneSig_correlation/results/immunecell_TCGAgenes_result_spearman_TUMERIC/"
  dir.create(resultpath_TUMERIC2)

# limit to protein_coding_genes ------------------------------------------------
protein_coding_genes <- read.table("03_drug_immuneSig_enrichment/data/protein_coding_genes_hg38.txt", sep = '\t', header = T)

 immune_sig = immune_sigatures
 tumor_purity_method = "TUMERIC"
# run correlation calculation in parallel mode ---------------------------------
project_ids <- c( "BLCA",
                  "LIHC",
                  "LUAD",
                  "BRCA",
                  "LUSC",
                  "CESC",
                  "OV",
                  "COAD",
                  "PAAD",
                  "ESCA",
                  "PRAD",
                  "GBM",
                  "READ",
                  "HNSC",
                  "KIRC",
                  "STAD",
                  "KIRP",
                  "LGG",
                  "THCA",
                  "UCEC") 
  
sample_type = 'Primary'

res=list()
for(i in 1:length(project_ids) ){
  cancer = project_ids[i]
  processedDataPath <- '02_tumorGene_immuneSig_correlation/data/TCGA_processedData/'
  processedDataPath <- paste0(processedDataPath, cancer,'/')
    
    TPMMatPath = paste0(processedDataPath, "mixture_file_TCGA_", cancer, "_", sample_type, "_TPM.txt")
    mymatrix_filter <- read.csv(TPMMatPath, sep = '\t', row.names = 1)
    colnames(mymatrix_filter) <- gsub('\\.', '-',   colnames(mymatrix_filter))
    colnames(mymatrix_filter) <- substr(colnames(mymatrix_filter),1,15)
    mymatrix_filter <- mymatrix_filter[is.element(row.names(mymatrix_filter), protein_coding_genes$gene_name),]
    
    # matching files for correlation
    immune_sig$SampleID2 <- substr(immune_sig$ID, 1, 15)
    c.mat <- immune_sig[immune_sig$SampleID2 %in% unique(colnames(mymatrix_filter)),]
    c.mat <- c.mat[order(c.mat$SampleID2),]
    row.names(c.mat) <- c.mat$SampleID2
    
    t.mat <- mymatrix_filter[colnames(mymatrix_filter) %in% c.mat$SampleID2,]
    t.mat <- mymatrix_filter[, c.mat$SampleID2]
    
    tumor_purity_filter <- tumor_purity[is.element(tumor_purity$Sample_ID,  c.mat$SampleID2),]
    tumor_purity_filter2 <- unique(tumor_purity_filter[c('Sample_ID',tumor_purity_method)])
    names(tumor_purity_filter2) <- c('Sample_ID','purity')
    
    c.mat <- c.mat[which(is.element(unique(row.names(c.mat)), tumor_purity_filter2$Sample_ID)), ]
    c.mat <- c.mat[order(row.names(c.mat)),]
    t.mat <- t.mat[, which(is.element(unique(colnames(t.mat)), tumor_purity_filter2$Sample_ID)) ]
    t.mat <- t.mat[, order(colnames(t.mat))]
    t.mat <- log2(t.mat + 1)
    
    c.mat <- c.mat[- which(names( c.mat) %in% c("SampleID","SampleID2", 
                                                "PatientID",
                                                "TCGA.Study_immunelandscape", 
                                                "immunelandscape_TCGA.Study",
                                                "ID"))]
   res[[i]] <-data.frame(cancer=cancer, n=ncol(t.mat))
}

res$sample_type='Primary'

i=i+1
cancer='SKCM'
sample_type = 'Metastatic'
  processedDataPath <- '02_tumorGene_immuneSig_correlation/data/TCGA_processedData/'
  processedDataPath <- paste0(processedDataPath, cancer,'/')
    
    TPMMatPath = paste0(processedDataPath, "mixture_file_TCGA_", cancer, "_", sample_type, "_TPM.txt")
    mymatrix_filter <- read.csv(TPMMatPath, sep = '\t', row.names = 1)
    colnames(mymatrix_filter) <- gsub('\\.', '-',   colnames(mymatrix_filter))
    colnames(mymatrix_filter) <- substr(colnames(mymatrix_filter),1,15)
    mymatrix_filter <- mymatrix_filter[is.element(row.names(mymatrix_filter), protein_coding_genes$gene_name),]
    
    # matching files for correlation
    immune_sig$SampleID2 <- substr(immune_sig$ID, 1, 15)
    c.mat <- immune_sig[immune_sig$SampleID2 %in% unique(colnames(mymatrix_filter)),]
    c.mat <- c.mat[order(c.mat$SampleID2),]
    row.names(c.mat) <- c.mat$SampleID2
    
    t.mat <- mymatrix_filter[colnames(mymatrix_filter) %in% c.mat$SampleID2,]
    t.mat <- mymatrix_filter[, c.mat$SampleID2]
    
    tumor_purity_filter <- tumor_purity[is.element(tumor_purity$Sample_ID,  c.mat$SampleID2),]
    tumor_purity_filter2 <- unique(tumor_purity_filter[c('Sample_ID',tumor_purity_method)])
    names(tumor_purity_filter2) <- c('Sample_ID','purity')
    
    c.mat <- c.mat[which(is.element(unique(row.names(c.mat)), tumor_purity_filter2$Sample_ID)), ]
    c.mat <- c.mat[order(row.names(c.mat)),]
    t.mat <- t.mat[, which(is.element(unique(colnames(t.mat)), tumor_purity_filter2$Sample_ID)) ]
    t.mat <- t.mat[, order(colnames(t.mat))]
    t.mat <- log2(t.mat + 1)
    
    c.mat <- c.mat[- which(names( c.mat) %in% c("SampleID","SampleID2", 
                                                "PatientID",
                                                "TCGA.Study_immunelandscape", 
                                                "immunelandscape_TCGA.Study",
                                                "ID"))]
   res[[i]] <-data.frame(cancer=cancer, n=ncol(t.mat))
res=do.call(rbind,res)

res$sample_type='Primary'
res[21,3] = 'Metastatic'
sum(res$n)



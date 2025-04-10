# try normalized count with DESeq2 ---------------------------------------------
# rm(list=ls())
# Data set preprocessing -------------------------------------------------------

setwd("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/")
options(stringsAsFactors = F)
library(biomaRt)
library(GEOquery)
library(DESeq2)
dir.create("data/drug_pertubation_animal/GSE120500/")

GSE120500 <- getGEO('GSE120500', 
                    destdir = "data/drug_pertubation_animal/GSE120500/",  
                    GSEMatrix=TRUE)
show(GSE120500)
show(pData(phenoData(GSE120500[[1]]))[,c(1, 2, 43)])
metadata = pData(phenoData(GSE120500[[1]]))[,c(1, 2, 43)]
metadata$labels = c(paste0("L0",1:9),paste0("L", 10:12))
# load data
library(readxl)
counts <- as.data.frame(read_excel("data/drug_pertubation_animal/GSE120500/GSE120500_Mouse_Transcriptome_control_and_parpi.xls", 
                     sheet = 1))
row.names(counts) <- counts$Gene
counts <- counts[-1]

sf <- estimateSizeFactorsForMatrix(as.matrix(counts))
nCount_selected <- t(apply(counts, 1, function(x){x / sf}))

## filter out genes with low expression level
nCount_selected <- nCount_selected[apply(nCount_selected, 1, function(x){sum(x>1)>length(x)*0.9}),]
nCount_selected <- as.data.frame(nCount_selected)

# Gene name conversion ---------------------------------------------------------
genename = read.csv(
  paste0(
    '/picb/bigdata/project/FengFYM/p2_2_liver_cancer_specific_immuno_combo_therapy/',
    "03_mouse_data_analysis/data/mousedata_for_immuClustering_icbResponse/mouse2human_count_manually_genename_merge.txt"
  ),
  sep = '\t',
  header = T
)


nCount_selected$MGI_symbol  = rownames(nCount_selected)
nCount_selected = inner_join(genename, nCount_selected, by = "MGI_symbol")
nCount_selected = nCount_selected[!duplicated(nCount_selected$hgnc_symbol),]
rownames(nCount_selected) = nCount_selected$hgnc_symbol
nCount_selected = nCount_selected[-c(1,2,3)]

# rownames(nCount_selected) = toupper(rownames(nCount_selected))

# generating expression files for GSEA -----------------------------------------
# pre-, on- and post-AZD5363 biopsies (labeled as 1, 2, and 3, respectively) 
mat_gsea <- cbind(data.frame(NAME = rownames(nCount_selected)), 
                 data.frame(DESCRIPTION = rep(NA, nrow(nCount_selected))),
                 log2(1+nCount_selected[c(7:12,1:6)]))
write.table(mat_gsea, "data/drug_pertubation_animal/GSE120500/GSE120500RNASeq_treatedVsvehicle_nCounts_gsea.txt",
            row.names = F, quote = F, sep = "\t")

# generating class file for GSEA
grp <- list(paste(12, 2, 1, collapse = " "),
            paste("#", "treated", "vehicle", collapse = " "),
            paste(c(rep("treated", 6),
                    rep("vehicle", 6)),
                  collapse = " "))
grp <- do.call(rbind, grp)
write.table(grp, "data/drug_pertubation_animal/GSE120500/GSE120500RNASeq_treatedVsvehicle_nCounts_gsea.cls", 
            col.names = F, row.names = F, quote = F, sep = "\t")


# detect the number of DEG ---------------------------------------------------------
# counts
source("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/scripts/functions/DEGanalysis.R")

metadata = data.frame(sampleID = metadata$labels,
                      treatment = metadata$`treatment:ch1`)

rownames(metadata) <- metadata$sampleID


data <- counts[apply(counts, 1, function(x){sum(x>1)>length(x)*0.9}),]
dds <- DESeq(DESeqDataSetFromMatrix(countData = data,
                                    colData = metadata['treatment'],
                                    design= ~ treatment))
dataset = 'GSE120500'
cntrl = "Vehicle"
treat = "Olaparib"
padj_threshold = 0.1
logFC_abs_threshold = 1
resultpath = paste0("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/result_", dataset, '/')

DEGdetection(dds = dds, 
             dataset = dataset,
             cntrl = cntrl,
             treat = treat,
             padj_threshold = padj_threshold,
             logFC_abs_threshold = logFC_abs_threshold,
             resultpath = resultpath)


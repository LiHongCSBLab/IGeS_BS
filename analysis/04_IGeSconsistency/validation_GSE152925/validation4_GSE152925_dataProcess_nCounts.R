# ============================================================================= #
# DESeq2 side factor normalization counts
# ============================================================================= #
# rm(list=ls())
# Data preparation -------------------------------------------------------------

setwd("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/")
options(stringsAsFactors = F)
library(biomaRt)
library(GEOquery)
library(DESeq2)

dir.create("data/immune_combo_dataset/GSE152925/")

gse152925 <- getGEO('GSE152925', 
                    destdir = "data/immune_combo_dataset/GSE152925/",  
                    GSEMatrix=TRUE)
show(gse152925)
show(pData(phenoData(gse152925[[1]]))[,c(1, 8, 19, 45)])
show(pData(phenoData(gse152925[[2]]))[,c(1, 8, 19, 45)])
metadata1 = pData(phenoData(gse152925[[1]]))[,c(1, 8, 19, 45)]
metadata2 = pData(phenoData(gse152925[[2]]))[,c(1, 8, 19, 45)]
metadata_all = rbind(metadata1, metadata2)

# load data
counts <- read.csv(gzfile("data/immune_combo_dataset/GSE152925/GSE152925_B16_C57BL6_invivo_RawCount.txt.gz"), 
                        sep = '\t', row.names = 1)
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

metadata = metadata_all[is.element(metadata_all$`treatment:ch1`, 
                        c("Swainsonine +IgG", "IgG +PBS")), ]
nCount_selected <- nCount_selected[is.element(names(nCount_selected),  metadata$description)]


# generating expression files for GSEA -----------------------------------------
# pre-, on- and post-AZD5363 biopsies (labeled as 1, 2, and 3, respectively) 
mat_gsea <- cbind(data.frame(NAME = rownames(nCount_selected)), 
                 data.frame(DESCRIPTION = rep(NA, nrow(nCount_selected))),
                 log2(1+nCount_selected[c(4:6,1:3)]))
write.table(mat_gsea, "data/immune_combo_dataset/GSE152925/GSE152925RNASeq_treatedVsvehicle_nCounts_gsea.txt",
            row.names = F, quote = F, sep = "\t")

# generating class file for GSEA
grp <- list(paste(6, 2, 1, collapse = " "),
            paste("#", "treated", "vehicle", collapse = " "),
            paste(c(rep("treated", 3),
                    rep("vehicle", 3)),
                  collapse = " "))
grp <- do.call(rbind, grp)
write.table(grp, "data/immune_combo_dataset/GSE152925/GSE152925RNASeq_treatedVsvehicle_nCounts_gsea.cls", 
            col.names = F, row.names = F, quote = F, sep = "\t")



# detect the number of DEG ---------------------------------------------------------
# counts
source("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/scripts/functions/DEGanalysis.R")

metadata = data.frame(sampleID = metadata_all$description,
                      treatment = metadata_all$`treatment:ch1`)

rownames(metadata) <- metadata$sampleID
metadata <- metadata[is.element(metadata$sampleID, colnames(counts)), ]
counts <- counts[is.element(colnames(counts),  metadata$sampleID), ]

data <- counts[apply(counts, 1, function(x){sum(x>1)>length(x)*0.9}),]
dds <- DESeq(DESeqDataSetFromMatrix(countData = data,
                                    colData = metadata['treatment'],
                                    design= ~ treatment))
dataset = 'GSE152925'
cntrl = "IgG +PBS"
treat = "Swainsonine +IgG"
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


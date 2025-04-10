# ============================================================================= #
# DESeq2 side factor normalization counts
# ============================================================================= #
# rm(list=ls())
# Data preparation -------------------------------------------------------------
setwd("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/")
options(stringsAsFactors = F)

library(GEOquery)
library(dplyr)
library(DESeq2)

genename = read.csv(
  paste0(
    "/picb/bigdata/project/FengFYM/p2_2_liver_cancer_specific_immuno_combo_therapy/",
    "03_mouse_data_analysis/data/mousedata_for_immuClustering_icbResponse/mouse2human_count_manually_genename_merge.txt"
  ),
  sep = "\t",
  header = T
)

# Gene name conversion ---------------------------------------------------------
counts <- read.csv("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/data/immune_combo_dataset/GSE149825/GSE149825_SMAC_ICB_RNASeq_rawcount.csv", row.names = 1)

sf <- estimateSizeFactorsForMatrix(as.matrix(counts))
nCount_selected <- t(apply(counts, 1, function(x){x / sf}))

## filter out genes with low expression level
nCount_selected <- nCount_selected[apply(nCount_selected, 1, function(x){sum(x>1)>length(x)*0.9}),]
nCount_selected <- as.data.frame(nCount_selected)

nCount_selected$MGI_symbol  = rownames(nCount_selected)
nCount_selected = inner_join(genename, nCount_selected, by = "MGI_symbol")
nCount_selected = nCount_selected[!duplicated(nCount_selected$hgnc_symbol),]
rownames(nCount_selected) = nCount_selected$hgnc_symbol
nCount_selected = nCount_selected[-c(1,2,3)]


# generating expression files for GSEA
mat_gsea <- cbind(data.frame(NAME = rownames(nCount_selected)), 
                 data.frame(DESCRIPTION = rep(NA, nrow(nCount_selected))),
                 log2(1+nCount_selected[c(4:6, 1:3)]))
write.table(mat_gsea, "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/data/immune_combo_dataset/GSE149825/GSE149825RNASeq_ncounts_gsea.txt",
            row.names = F, quote = F, sep = "\t")

# generating class file for GSEA
grp <- list(paste(6, 2, 1, collapse = " "),
            paste("#", "treated", "vehicle", collapse = " "),
            paste(c(rep("treated", 3),
                    rep("vehicle", 3)),
                  collapse = " "))
grp <- do.call(rbind, grp)
write.table(grp, "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/data/immune_combo_dataset/GSE149825/GSE149825RNASeq_ncounts_gsea.cls", 
            col.names = F, row.names = F, quote = F, sep = "\t")


# detect the number of DEG ---------------------------------------------------------
# counts
source("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/scripts/functions/DEGanalysis.R")

metadata = data.frame(sampleID = colnames(counts),
                      treatment = c(rep('Vehicle',3) , rep('SMAC',3), rep('ICB',3), rep('SMAC.ICB',3)))
rownames(metadata) <- metadata$sampleID

data <- counts[apply(counts, 1, function(x){sum(x>1)>length(x)*0.9}),]
dds <- DESeq(DESeqDataSetFromMatrix(countData = data,
                                    colData = metadata['treatment'],
                                    design= ~ treatment))
dataset = 'GSE149825'
treat = "Vehicle"
cntrl = "SMAC"
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

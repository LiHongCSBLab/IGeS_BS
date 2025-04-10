# =============================================================================#
# DESeq2 size factor normalization counts
# =============================================================================#

# ------------------------------------------------------------------------------
setwd("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/")
options(stringsAsFactors = F)
library(DESeq2)

genename = read.csv(
  paste0(
    '/picb/bigdata/project/FengFYM/p2_2_liver_cancer_specific_immuno_combo_therapy/',
    "03_mouse_data_analysis/data/mousedata_for_immuClustering_icbResponse/mouse2human_count_manually_genename_merge.txt"
  ),
  sep = '\t',
  header = T
)

# Gene name conversion
counts <- read.csv("data/immune_combo_dataset/GSE114601/counts.raw.csv", row.names = 1)
metadata <- data.frame(sampleID = names(counts),
                       # s1795 s2467 s2503 s2521 s2596 s2617 s2625 s2669
                       treatment = c("antiPD1","combo","JQ1","antiPD1","vehicle","vehicle","JQ1","combo"))



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
                 log2(1+nCount_selected[which(is.element(names(nCount_selected), metadata[metadata$treatment == 'JQ1',]$sampleID))]),
                 log2(1+nCount_selected[which(is.element(names(nCount_selected), metadata[metadata$treatment == 'vehicle',]$sampleID))]))

write.table(mat_gsea, "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/data/immune_combo_dataset/GSE114601/GSE114601RNASeq_nCounts_gsea.txt",
            row.names = F, quote = F, sep = "\t")

# generating class file for GSEA
grp <- list(paste(4, 2, 1, collapse = " "),
            paste("#", "treated", "vehicle", collapse = " "),
            paste(c(rep("treated", 2),
                    rep("vehicle", 2)),
                  collapse = " "))
grp <- do.call(rbind, grp)
write.table(grp, "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/data/immune_combo_dataset/GSE114601/GSE114601RNASeq_nCounts_gsea.cls", 
            col.names = F, row.names = F, quote = F, sep = "\t")


# detect the number of DEG ---------------------------------------------------------
# counts
source("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/scripts/functions/DEGanalysis.R")
rownames(metadata) <- metadata$sampleID

data <- counts[apply(counts, 1, function(x){sum(x>1)>length(x)*0.9}),]
dds <- DESeq(DESeqDataSetFromMatrix(countData = data,
                                    colData = metadata['treatment'],
                                    design= ~ treatment))
dataset = 'GSE114601'
cntrl = "vehicle"
treat = "JQ1"
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

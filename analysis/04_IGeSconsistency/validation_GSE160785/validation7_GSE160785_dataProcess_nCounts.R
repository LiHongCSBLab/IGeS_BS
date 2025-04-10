# try normalized count with DESeq2 ---------------------------------------------
# rm(list=ls())
# Data set preprocessing -------------------------------------------------------
setwd("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/")
options(stringsAsFactors = F)
library(biomaRt)
library(GEOquery)
library(limma)
library(DESeq2)
dir.create("data/drug_pertubation_animal/GSE160785/")

# load series and platform data from GEO
gset <- getGEO("GSE160785",
               destdir = "data/drug_pertubation_animal/GSE160785/", 
               GSEMatrix =TRUE, getGPL=F)

metadata_all = pData(phenoData(gset[[1]]))[,c('title', 'geo_accession', 'treatment:ch1')]

counts <- read.csv(gzfile("data/drug_pertubation_animal/GSE160785/GSE160785_SZ9_raw_count_matrix.txt.gz"), 
                        sep = '\t', row.names = 1)

metadata = metadata_all[is.element(metadata_all$`treatment:ch1`, 
                        c("vehicle", "celecoxib")), ]
nCount_selected <- counts[is.element(names(counts),  metadata$title)]

sf <- estimateSizeFactorsForMatrix(as.matrix(nCount_selected))
nCount_selected <- t(apply(nCount_selected, 1, function(x){x / sf}))

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

nCount_selected$MusEnsembleID  = rownames(nCount_selected)
nCount_selected = inner_join(genename, nCount_selected, by = 'MusEnsembleID')
nCount_selected = nCount_selected[!duplicated(nCount_selected$hgnc_symbol),]
rownames(nCount_selected) = nCount_selected$hgnc_symbol
nCount_selected = nCount_selected[-c(1,2,3)]


# generating expression files for GSEA -----------------------------------------
metadata = metadata[order(metadata$`treatment:ch1`), ]
mat_gsea <- cbind(data.frame(NAME = rownames(nCount_selected)), 
                 data.frame(DESCRIPTION = rep(NA, nrow(nCount_selected))),
                 log2(1 + nCount_selected[metadata$title]))
write.table(mat_gsea, "data/drug_pertubation_animal/GSE160785/celecoxib_treatedVsvehicle_gsea_nCount.txt",
            row.names = F, quote = F, sep = "\t")

# generating class file for GSEA
grp <- list(paste(ncol(nCount_selected), 2, 1, collapse = " "),
            paste("#", "treated", "vehicle", collapse = " "),
            paste(c(rep("treated", nrow(metadata[metadata$`treatment:ch1` == 'celecoxib', ])),
                    rep("vehicle", nrow(metadata[metadata$`treatment:ch1` == 'vehicle', ]))),
                  collapse = " "))
grp <- do.call(rbind, grp)
write.table(grp, "data/drug_pertubation_animal/GSE160785/celecoxib_treatedVsvehicle_gsea_nCount.cls", 
            col.names = F, row.names = F, quote = F, sep = "\t")


# detect the number of DEG ---------------------------------------------------------
# counts
source("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/scripts/functions/DEGanalysis.R")

metadata = data.frame(sampleID = metadata_all$title,
                      treatment = metadata_all$`treatment:ch1`)

rownames(metadata) <- metadata$sampleID
metadata <- metadata[is.element(metadata$sampleID, colnames(counts)), ]
counts <- counts[is.element(colnames(counts),  metadata$sampleID), ]

data <- counts[apply(counts, 1, function(x){sum(x>1)>length(x)*0.9}),]
dds <- DESeq(DESeqDataSetFromMatrix(countData = data,
                                    colData = metadata['treatment'],
                                    design= ~ treatment))
dataset = 'GSE160785'
cntrl = "vehicle"
treat = "celecoxib"
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


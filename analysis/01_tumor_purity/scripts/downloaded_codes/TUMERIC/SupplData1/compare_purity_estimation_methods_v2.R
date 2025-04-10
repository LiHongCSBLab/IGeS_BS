# Example link for LUAD type data. To download data for another type, substitute LUAD with [cancer type] 
# The link below is log2(counts+1) (ensebml ids)  
# https://gdc.xenahubs.net/download/TCGA-LUAD/Xena_Matrices/TCGA-LUAD.htseq_counts.tsv.gz
# The link below in PANCANCER data in log2(TPM+0.001), ENSEMBL IDs
# https://toil.xenahubs.net/download/tcga_RSEM_gene_tpm.gz

library(tidyverse)
library(plyr)

### Read purity estimates of other methods
methods.purity = read.csv('data/purity_separate.csv',row.names = 1) # previous analysis results  (tumor purities from other methods)
methods.purity <- methods.purity[complete.cases(methods.purity), ] # only leave the samples that have all the purities available
# table with protein coding genes - their ENSEMBL IDs and HGNC symbols
# please download before running this script ...
protein.coding.genes.names <- read.csv("../data/proteinc_gene_name_conversion_table.csv", row.names = 1) 

### Read gene expression data
ctype = 'LUAD'
# ctype = 'HNSC'
# ctype = 'LUSC'
# ctype = 'BRCA'
# ctype = 'GBM'
# ctype = 'LGG' # can also combine normal samples from "GBM"
# ctype = 'SKCM' 
# ctype = 'THCA'

#  counts data
# please download before running this script ...
log_counts.path = paste0("../data/TCGA_RNA-seq/TCGA-", ctype, ".htseq_counts.tsv")
log2.counts <- read.csv(log_counts.path, sep = '\t', row.names = 1, check.names = FALSE)

# TPM data (use only for LUAD, HNSC and LUSC)
# please download before running this script ...
log2.tpm_pancan <- read.csv("../data/TCGA_RNA-seq/tcga_tpm_pancan_60k_10535.csv", sep='\t', row.names = 1, check.names = FALSE)
rownames(log2.tpm_pancan) <- sapply(strsplit(rownames(log2.tpm_pancan),"\\."), `[`, 1) # cutting version postfixes of ensemble ids

### data frame to matrix conversion + formatting
log2.counts = data.matrix(log2.counts) 
colnames(log2.counts) <- gsub('.{1}$', '', colnames(log2.counts)) # remove the last character in counts matrix (e.g. letter "A" representing vial)
rownames(log2.counts) <- sapply(strsplit(rownames(log2.counts),"\\."), `[`, 1) # cutting version postfixes of ensemble ids 
# leave only protein-coding genes 
log2.counts <- log2.counts[rownames(log2.counts) %in% protein.coding.genes.names$ensembl_id, ]

# get the names of samples with complete data from other methods and types of data
complete_samples_tumor <- intersect(rownames(methods.purity), colnames(log2.counts))

# move to linear scale, log2(x+1)
counts <- 2**log2.counts - 1

### choose normal and tumor samples
nsamples = 100 # limit the amount of mixed (tumor) samples

normal.counts <- (counts[, grep('.11$',colnames(counts)), drop=FALSE]) # get all possible normal samples (we're only using them for DeMixT)
tumor.counts <- (counts[, complete_samples_tumor, drop=FALSE])[, 1:nsamples] 
# names of genes in HGNC nomenclature
hgnc_counts_gene_names <- mapvalues(rownames(tumor.counts), 
                                 from=as.character(protein.coding.genes.names[protein.coding.genes.names$ensembl_id %in% rownames(tumor.counts), ]$ensembl_id), 
                                 to=as.character(protein.coding.genes.names[protein.coding.genes.names$ensembl_id %in% rownames(tumor.counts), ]$hgnc_sym))


# for GBM
# normal.counts.GBM <- normal.counts
# normal.counts <- normal.counts.GBM # for LGG, since GBM and LGG normal samples are ~same

## only for TPM values
tumor.log2.tpm <- log2.tpm_pancan[, complete_samples_tumor][, colnames(tumor.counts)]
tumor.log2.tpm <- tumor.log2.tpm[rownames(tumor.log2.tpm) %in% protein.coding.genes.names$ensembl_id, ]
tumor.tpm <- 2**tumor.log2.tpm - 0.001
hgnc_tpm_gene_names <- mapvalues(rownames(tumor.tpm),
          from=as.character(protein.coding.genes.names[protein.coding.genes.names$ensembl_id %in% rownames(tumor.tpm), ]$ensembl_id),
          to=as.character(protein.coding.genes.names[protein.coding.genes.names$ensembl_id %in% rownames(tumor.tpm), ]$hgnc_sym))

###########
### DeMixT-normalization
# quartile normalization provided by the authors of DeMixT in their email
library(DSS)

quart_normalize_for_demixt <- function(normal.counts, tumor.counts) {
  
  Count.matrix <- cbind(normal.counts, tumor.counts)
  newt <- Count.matrix
  colnames(newt)=NULL
  rownames(newt)=NULL
  
  normal.id <- colnames(normal.counts)
  tumor.id <- colnames(tumor.counts)
  
  designs=c(rep("1", length(normal.id)), rep("0", length(tumor.id)))
  seqData=newSeqCountSet(as.matrix(newt), designs)
  
  # Quartile normalization/total/median (median is provided by the authors, seems to work fine)
  seqData=estNormFactors(seqData, "quantile")
  k3=seqData@normalizationFactor
  mk3=median(k3)
  k3=k3/mk3
  
  temp <- newt
  
  for(i in 1:ncol(newt)){
    temp[,i] = temp[,i]/k3[i]
  }
  Count.matrix.normalized <- temp
  colnames(Count.matrix.normalized) <- colnames(Count.matrix)
  rownames(Count.matrix.normalized) <- rownames(Count.matrix)
  
  normal.counts.quartnorm <- Count.matrix.normalized[, colnames(normal.counts), drop=FALSE]
  tumor.counts.quartnorm <- Count.matrix.normalized[, colnames(tumor.counts), drop=FALSE]
  
  return(list("normal" = normal.counts.quartnorm, 
              "tumor" = tumor.counts.quartnorm))
}

quartnorm_counts_list <- quart_normalize_for_demixt(normal.counts, tumor.counts)
normal.counts.quartnorm <- quartnorm_counts_list$normal
tumor.counts.quartnorm <- quartnorm_counts_list$tumor

### DeMixT-method
# Note: author's used quartile-normalized counts data, although they say 
#   that theoretically DeMixT should accept any input expression as long as the distribution follows
#   log2-normal distribution    
#devtools::install_github("wwylab/DeMixTallmaterials/DeMixT_0.2")
library(DeMixT)

res <- DeMixT(data.Y = tumor.counts.quartnorm, 
      data.comp1 = normal.counts.quartnorm) 
demix.stromal <- res$pi[1,]
# saving results
write.csv(demix.stromal, paste0("../results/methods_benchmark/re-run/DeMixT_", ctype, "_stromal_counts_quartnorm.csv"))
###########
              
###########
### LinSeed
# devtools::install_github("ctlab/linseed")
library(linseed)

# normalization recommended by the package authors: linear space + normalized column-wise so that the columns (samples) have the same sum
expData <- tumor.counts
rownames(expData) <- hgnc_counts_gene_names
expData <- expData[!grepl("RPL|RPS", rownames(expData)), ] # removing RPL|RPS genes
expData[is.nan(expData)] <- 0
totalSum <- mean(colSums(expData))
expData <- apply(expData, 2, function(x) x * totalSum / sum(x))

lo <- LinseedObject$new(expData, topGenes=10000)

# Colinearity networks
# the dataset has to be in lin-scale and normalized
lo$calculatePairwiseLinearity()
lo$calculateSpearmanCorrelation()
lo$calculateSignificanceLevel(100)
# lo$significancePlot(0.01)

lo$filterDatasetByPval(0.01)

lo$setCellTypeNumber(2)
lo$project("filtered")
# lo$projectionPlot(color="filtered")
lo$smartSearchCorners(dataset="filtered", error="norm")
# deconvolution
lo$deconvolveByEndpoints()
# lo$selectGenes(100)
# lo$tsnePlot()

linseed.purity <- lo$proportions[2, ] # 2nd comp is purity judging by the values / good correlation
# saving results
write.csv(linseed.purity, paste0("../results/methods_benchmark/re-run/LinSeed_", ctype, "_purity_counts.csv"))
###########

###########
### CibersortX 
# uses its own online interface, requires non-log data and gene names in HGNC symbol nomenclature
# signature HNSCC matrix is found in Cibersort's paper supplementary 2e
norm_counts.cibersort <- tumor.tpm
norm_counts.cibersort.df <- rbind(setNames(data.frame(t(c("gene", colnames(norm_counts.cibersort))), stringsAsFactors = FALSE), c("gene", colnames(norm_counts.cibersort))), 
                                  setNames(data.frame(cbind(hgnc_tpm_gene_names, norm_counts.cibersort), check.names = FALSE, stringsAsFactors = FALSE), c("gene", colnames(norm_counts.cibersort))))
# saving results
write.table(norm_counts.cibersort.df, file = paste0("../data/TCGA_RNA-seq/", ctype, ".TPM.for_Cibersort.tsv"), sep='\t',
            quote=FALSE, col.names=FALSE, row.names=FALSE)

###########
### Loading analysis results and plotting

### Read purity estimates of other methods
methods.purity <- read.csv('../data/purity_separate.csv', row.names = 1) # previous analysis results
methods.purity <- methods.purity[complete.cases(methods.purity), ] # only leave the samples that have all the purities available

genomic.average.purity <- rowSums(methods.purity[, c('AbsCN.seq', 'ASCAT', 'PurBayes')]) / 3

methods = c('DeMixT', 'CiberSortX', 'LinSeed', 'ESTIMATE', 'TUMERIC')

# function to combine the purities from differene methods into one dataframe    
create_purity_df <- function(demix.purity, linseed.purity, cibersort.purity, ctype) {
  
  samples_names <- rownames(linseed.purity)
  
  all.purity = t(t(demix.purity))
  all.purity = cbind(all.purity, cibersort.purity)
  all.purity = cbind(all.purity, linseed.purity)
  all.purity = cbind(all.purity, methods.purity[samples_names, 'ESTIMATE'])
  all.purity = cbind(all.purity, methods.purity[samples_names, 'Purity'])
  all.purity <- cbind(all.purity, genomic.average.purity[samples_names])
  
  colnames(all.purity) = c(methods, 'genomic_avg')
  
  all.purity$type = ctype
  
  return(all.purity)
}


# LUAD results
ctype = 'LUAD'
demix.purity <- 1 - read.csv(paste0("../results/methods_benchmark/re-run/DeMixT_", ctype, "_stromal_counts_quartnorm.csv"), row.names=1)
linseed.purity <- read.csv(paste0("../results/methods_benchmark/re-run/LinSeed_", ctype, "_purity_counts.csv"), row.names=1)
cibersort.purity <- read.csv(paste0('../results/methods_benchmark/re-run/CIBERSORTx_', ctype, '_tpm_results.csv'), row.names=1)
rownames(cibersort.purity) <- cibersort.purity$Mixture; 
cibersort.purity <- select(cibersort.purity, Malignant.cells)

all.purity <- create_purity_df(demix.purity, linseed.purity, cibersort.purity, ctype)
assign(paste0("all.purity.", ctype), all.purity)

# HNSC results
ctype = 'HNSC'
demix.purity <- 1 - read.csv(paste0("../results/methods_benchmark/re-run/DeMixT_", ctype, "_stromal_counts_quartnorm.csv"), row.names=1)
linseed.purity <- read.csv(paste0("../results/methods_benchmark/re-run/LinSeed_", ctype, "_purity_counts.csv"), row.names=1)
cibersort.purity <- read.csv(paste0('../results/methods_benchmark/re-run/CIBERSORTx_', ctype, '_tpm_results.csv'), row.names=1)
rownames(cibersort.purity) <- cibersort.purity$Mixture; 
cibersort.purity <- select(cibersort.purity, Malignant.cells)

all.purity <- create_purity_df(demix.purity, linseed.purity, cibersort.purity, ctype)
assign(paste0("all.purity.", ctype), all.purity)

# LUSC results
ctype = 'LUSC'
demix.purity <- 1 - read.csv(paste0("../results/methods_benchmark/re-run/DeMixT_", ctype, "_stromal_counts_quartnorm.csv"), row.names=1)
linseed.purity <- read.csv(paste0("../results/methods_benchmark/re-run/LinSeed_", ctype, "_purity_counts.csv"), row.names=1)
cibersort.purity <- read.csv(paste0('../results/methods_benchmark/re-run/CIBERSORTx_', ctype, '_tpm_results.csv'), row.names=1)
rownames(cibersort.purity) <- cibersort.purity$Mixture; 
cibersort.purity <- select(cibersort.purity, Malignant.cells)

all.purity <- create_purity_df(demix.purity, linseed.purity, cibersort.purity, ctype)
assign(paste0("all.purity.", ctype), all.purity)

# BRCA results
ctype = 'BRCA'
demix.purity <- 1 - read.csv(paste0("../results/methods_benchmark/re-run/DeMixT_", ctype, "_stromal_counts_quartnorm.csv"), row.names=1)
linseed.purity <- read.csv(paste0("../results/methods_benchmark/re-run/LinSeed_", ctype, "_purity_counts.csv"), row.names=1)
cibersort.purity <- rep(NA, dim(linseed.purity)[1]) # this type doesn't have Cibersort results

all.purity <- create_purity_df(demix.purity, linseed.purity, cibersort.purity, ctype)
assign(paste0("all.purity.", ctype), all.purity)

# GBM results
ctype = 'GBM'
demix.purity <- 1 - read.csv(paste0("../results/methods_benchmark/re-run/DeMixT_", ctype, "_stromal_counts_quartnorm.csv"), row.names=1)
linseed.purity <- read.csv(paste0("../results/methods_benchmark/re-run/LinSeed_", ctype, "_purity_counts.csv"), row.names=1)
cibersort.purity <- rep(NA, dim(linseed.purity)[1]) # this type doesn't have Cibersort results

all.purity <- create_purity_df(demix.purity, linseed.purity, cibersort.purity, ctype)
assign(paste0("all.purity.", ctype), all.purity)

# LGG results
ctype = 'LGG'
linseed.purity <- read.csv(paste0("../results/methods_benchmark/re-run/LinSeed_", ctype, "_purity_counts.csv"), row.names=1)
demix.purity <- 1 - read.csv(paste0("../results/methods_benchmark/re-run/DeMixT_", ctype, "_stromal_counts_quartnorm.csv"), row.names=1)
cibersort.purity <- rep(NA, dim(linseed.purity)[1]) # this type doesn't have Cibersort results

all.purity <- create_purity_df(demix.purity, linseed.purity, cibersort.purity, ctype)
assign(paste0("all.purity.", ctype), all.purity)

# HNSC results
ctype = 'SKCM'
linseed.purity <- read.csv(paste0("../results/methods_benchmark/re-run/LinSeed_", ctype, "_purity_counts.csv"), row.names=1)
demix.purity <- rep(NA, dim(linseed.purity)[1]) # this type doesn't have Demixt results, cause there's only 1 normal sample
cibersort.purity <- rep(NA, dim(linseed.purity)[1]) # this type doesn't have Cibersort results

all.purity <- create_purity_df(demix.purity, linseed.purity, cibersort.purity, ctype)
assign(paste0("all.purity.", ctype), all.purity)

# THCA results
ctype = 'THCA'
linseed.purity <- read.csv(paste0("../results/methods_benchmark/re-run/LinSeed_", ctype, "_purity_counts.csv"), row.names=1)
demix.purity <- 1 - read.csv(paste0("../results/methods_benchmark/re-run/DeMixT_", ctype, "_stromal_counts_quartnorm.csv"), row.names=1)
cibersort.purity <- rep(NA, dim(linseed.purity)[1]) # this type doesn't have Cibersort results

all.purity <- create_purity_df(demix.purity, linseed.purity, cibersort.purity, ctype)
assign(paste0("all.purity.", ctype), all.purity)

# # saving the results
# for (ctype in ctypes) {
#   all.purity <- get(paste0('all.purity.', ctype))
#   write.table(x = cbind(all.purity, methods.purity[rownames(all.purity), c('AbsCN.seq', 'ASCAT', 'PurBayes', 'CPE')]),
#               file = paste0("../results/methods_benchmark/", ctype, "_purities_all_methods.tsv"), sep = '\t') # for saving all methods in 1 table
# }
  

### Plotting 
library(ggplot2)
#library(cowplot)
#theme_set(theme_cowplot())

# Download TCGA pan-cancer consortium ABSOLUTE purity values:
# https://gdc.cancer.gov/about-data/publications/pancanatlas
absolute = read.table('TCGA_mastercalls.abs_tables_JSedit.fixed.txt',header=T,sep="\t")[,c(2,4)]
absolute = absolute[substr(absolute[,1],1,4) == 'TCGA',]
absolute$sample = substr(absolute[,1],1,15)
colnames(absolute) = c('sample','ABSOLUTE')

cor.type <- function(f) {
  xl = read.table(f)[,-7]
  xl$sample = rownames(xl)
  xl = merge(xl,absolute)
  rownames(xl) = xl$sample
  xl = xl[,-1] # remove sample column
  cs = c()  
  ms = c('DeMixT','CiberSortX','LinSeed','TUMERIC')
  for(m in ms) {
    cs = rbind(cs,c('Method'=m,'Type'='ABS','Corr'=cor(xl,use='pairwise.complete.obs')[m,'ABSOLUTE']))
  }
  cs
}

cors = data.frame(cbind(cor.type('data/compare_purity/LUAD_purities_all_methods.tsv'),tumor='LUAD'))
cors = rbind(cors,cbind(cor.type('data/compare_purity/LUSC_purities_all_methods.tsv'),tumor='LUSC'))
cors = rbind(cors,cbind(cor.type('data/compare_purity/HNSC_purities_all_methods.tsv'),tumor='HNSC'))
cors = rbind(cors,cbind(cor.type('data/compare_purity/BRCA_purities_all_methods.tsv'),tumor='BRCA'))
cors = rbind(cors,cbind(cor.type('data/compare_purity/GBM_purities_all_methods.tsv'),tumor='GBM'))
cors = rbind(cors,cbind(cor.type('data/compare_purity/LGG_purities_all_methods.tsv'),tumor='LGG'))
cors = rbind(cors,cbind(cor.type('data/compare_purity/THCA_purities_all_methods.tsv'),tumor='THCA'))
cors = rbind(cors,cbind(cor.type('data/compare_purity/SKCM_purities_all_methods.tsv'),tumor='SKCM'))

ggplot(data=cors, aes(x=Method, y=as.numeric(as.character(Corr)), fill=tumor)) + geom_bar(stat="identity", position=position_dodge()) + ylab('Correlation with ABSOLUTE tumor purity (Pearson r)') + geom_hline(yintercept=c(0,0.25,0.5,0.75),linetype="dashed")


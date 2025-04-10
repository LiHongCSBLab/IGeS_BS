# Run the ipython and R markdown notebooks to reproduce the figures

# File and folder descriptions
purity_tumeric.csv : tumeric purity estimates

data/deconv_separate folder contains the separate deconvolution results from 20 cancers types, 5 BRCA subtypes and CPTAC OV and BRCA

data/output folder contains the outputs from the Figure1 ipython notebook which are the ligand-receptor product scores for pancancer and BRCA subtypes.

1370_pairs.txt : 1370 ligand-receptor pairs
1380_pairs.txt : 1380 ligand-receptor pairs

BRCA_pct_samp_fpkm1.csv: List of genes in BRCA by their expression in proportion of samples with fpkm>1
pct_samp_fpkm1.csv: List of genes in all cancer types by their expression in proportion of samples with fpkm>1

median_log2.csv: Median expression in BRCA subtypes [log2(fpkm+1)]

cancer_genes and stroma_genes: list of cancer and stroma specific genes from Tirosh et al. paper.
ESTIMATE_stroma_immune.txt: Stroma and immune specific genes from the ESTIMATE paper

cnv_exp_corr1.csv: correlation of expression vs cnv in cancer and stroma
gsea_cs.csv: log2(Cancer/Stroma) GSEA NES 

normal_median_logx+1.csv : Median expression in normal tissues [log2(fpkm+1)]

purity_separate.csv: purity estimates from different methods (AbsCN-seq, ASCAT, ESTIMATE, PurBayes and CPE)

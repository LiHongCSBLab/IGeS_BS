#!/bin/bash
# run under miniconda environment r_env
source /home/fengfangyoumin/miniconda3/bin/activate r_env

# cancer=$1
# sampleType=$2
# tumor_purity_method=$3
# method=$4
# sign=$5
# num_gene=$6
# dataset=$7



# for LIHC
# 01: run Rscript to generate expr file and cls file
    # parameter 1: "cancer",   "p", 1, "character", "the path of gene and fold-change",
    # parameter 2: "dataset",  "D", 1, "character", "dataset ID",

# 02 run enrichment analysis pipeline
# sample code:
# Rscript --vanilla 01_drug_expr_preparation_update.R -c SKCM -d 70138 > drug_expr_prep.log 2>&1 &

Rscript --vanilla 01_drug_expr_preparation_update.R -c LIHC -d 70138 > drug_expr_prep_LIHC_70138.log 2>&1 &
wait
bash 04_drug_immusig_GSEA_pipeline.sh LIHC Primary IHC pearson 70138  > drug_immusig_GSEA_LIHC_70138.log 2>&1 &
wait

Rscript --vanilla 01_drug_expr_preparation_update.R -c LIHC -d 92742 > drug_expr_prep_LIHC_92742.log 2>&1 &
wait
bash 04_drug_immusig_GSEA_pipeline.sh LIHC Primary IHC pearson 92742 > drug_immusig_GSEA_LIHC_92742.log 2>&1 &
wait
# clean files to release space
# rm -r /picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/03_drug_immuneSig_enrichment/results/drug_immunesig_IHC_GSEA_result/70138/LIHC
# rm -r /picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/03_drug_immuneSig_enrichment/data/drug_treated_expression_matrix/70138/LIHC

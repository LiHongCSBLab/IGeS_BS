#!/bin/bash
cancer=$1
sampleType=$2
tumor_purity_method=$3
method=$4
dataset=$5

# 02: run GSEA for each of cell line, then extract GSEA result
    # cancer=$1
    # sampleType=$2
    # tumor_purity_method=$3
    # method=$4
    # sign=$5
    # num_gene=$6
    # dataset=$7
# 04: merge enrichment score with different treatment time and dose
# 01: run Rscript to generate expr file and cls file
    # parameter 1: "cancer",   "p", 1, "character", "the path of gene and fold-change",
    # parameter 2: "dataset",  "D", 1, "character", "dataset ID"
    # parameter 3: "sign",     "s", 1, "character", "the work path for reading and saving data",
    # parameter 4: "num_gene", "n", 1, "character", "the path for saving result",

# Rscript --vanilla 03_summarize_drug_immuneSig_GSEA.R -c SKCM -d 70138 -s all -n 400 > drug_summarize.log 2>&1 &


# bash ./GSEA/02_drug_GSEA_update.sh SKCM Metastatic IHC spearman all 400 70138 > drug_GSEA.log 2>&1 &

# sign1=all
# num_gene1=400
# bash ./GSEA/02_drug_GSEA_update_978gene.sh ${cancer} ${sampleType} ${tumor_purity_method} ${method} ${sign1} ${num_gene1} ${dataset} > drug_GSEA_${cancer}_${dataset}_978gene.log 2>&1 &
# wait
# Rscript --vanilla 03_summarize_drug_immuneSig_GSEA.R -c ${cancer} -d ${dataset} -s ${sign1} -n ${num_gene1} -l 978genes > drug_summarize_${cancer}_${dataset}.log 2>&1 &
# wait

# sign2=positive
# num_gene2=200
# bash ./GSEA/02_drug_GSEA_update_978gene.sh ${cancer} ${sampleType} ${tumor_purity_method} ${method} ${sign2} ${num_gene2} ${dataset} > drug_GSEA_${cancer}_${dataset}_978gene.log 2>&1 &
# wait
# Rscript --vanilla 03_summarize_drug_immuneSig_GSEA.R -c ${cancer} -d ${dataset} -s ${sign3} -n ${num_gene3} -l 978genes > drug_summarize_${cancer}_${dataset}.log 2>&1 &
# wait


# sign3=negative
# num_gene3=200
# bash ./GSEA/02_drug_GSEA_update_978gene.sh ${cancer} ${sampleType} ${tumor_purity_method} ${method} ${sign3} ${num_gene3} ${dataset} > drug_GSEA_${cancer}_${dataset}_978gene.log 2>&1 &
# wait
# Rscript --vanilla 03_summarize_drug_immuneSig_GSEA.R -c ${cancer} -d ${dataset} -s ${sign3} -n ${num_gene3} -l 978genes > drug_summarize_${cancer}_${dataset}.log 2>&1 &
# wait


# sign1=all
# num_gene1=400
# bash ./GSEA/02_drug_GSEA_update.sh ${cancer} ${sampleType} ${tumor_purity_method} ${method} ${sign1} ${num_gene1} ${dataset} > drug_GSEA_${cancer}_${dataset}.log 2>&1 &
# wait
# Rscript --vanilla 03_summarize_drug_immuneSig_GSEA.R -c ${cancer} -t ${sampleType} -p ${tumor_purity_method} -d ${dataset} -s ${sign1} -n ${num_gene1} -l allgenes > drug_summarize_${cancer}_${dataset}.log 2>&1 &
# wait


sign2=positive
num_gene2=200
bash ./GSEA/02_drug_GSEA_update.sh ${cancer} ${sampleType} ${tumor_purity_method} ${method} ${sign2} ${num_gene2} ${dataset} > drug_GSEA_${cancer}_${dataset}.log 2>&1 #&
wait
Rscript --vanilla 03_summarize_drug_immuneSig_GSEA.R -c ${cancer} -t ${sampleType} -p ${tumor_purity_method} -d ${dataset} -s ${sign2} -n ${num_gene2} -l allgenes > drug_summarize_${cancer}_${dataset}.log 2>&1 # &
wait

# sign3=negative
# num_gene3=200
# # bash ./GSEA/02_drug_GSEA_update.sh ${cancer} ${sampleType} ${tumor_purity_method} ${method} ${sign3} ${num_gene3} ${dataset} > drug_GSEA_${cancer}_${dataset}.log 2>&1 &
# # wait
# Rscript --vanilla 03_summarize_drug_immuneSig_GSEA.R -c ${cancer} -t ${sampleType} -p ${tumor_purity_method} -d ${dataset} -s ${sign3} -n ${num_gene3} -l allgenes > drug_summarize_${cancer}_${dataset}.log 2>&1 &
# wait


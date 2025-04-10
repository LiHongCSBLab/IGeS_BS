#!/bin/bash

source /home/fengfangyoumin/miniconda3/bin/activate r_env

cancer=$1
sampleType=$2
tumor_purity_method=$3
dataset=$4
# sign=$5



sign3=negative
num_gene3=200
Rscript --vanilla 03_summarize_drug_immuneSig_GSEA.R -c ${cancer} -t ${sampleType} -p ${tumor_purity_method} -d ${dataset} -s ${sign3} -n ${num_gene3} -l allgenes > drug_summarize_${cancer}_${dataset}_${sign3}.log 2>&1 &
wait
 
# sign2=positive
# num_gene2=200
# Rscript --vanilla 03_summarize_drug_immuneSig_GSEA.R -c ${cancer} -t ${sampleType} -p ${tumor_purity_method} -d ${dataset} -s ${sign2} -n ${num_gene2} -l allgenes > drug_summarize_${cancer}_${dataset}.log 2>&1 &
# wait

# sign3=negative
# num_gene3=200
# Rscript --vanilla 03_summarize_drug_immuneSig_GSEA.R -c ${cancer} -t ${sampleType} -p ${tumor_purity_method} -d ${dataset} -s ${sign3} -n ${num_gene3} -l 978genes > drug_summarize_${cancer}_${dataset}.log 2>&1 &
# wait

# sign2=positive
# num_gene2=200
# Rscript --vanilla 03_summarize_drug_immuneSig_GSEA.R -c ${cancer} -t ${sampleType} -p ${tumor_purity_method} -d ${dataset} -s ${sign2} -n ${num_gene2} -l 978genes > drug_summarize_${cancer}_${dataset}.log 2>&1 &
# wait

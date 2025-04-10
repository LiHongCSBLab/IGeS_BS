#!/bin/bash
mkdir /picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/03_drug_immuneSig_enrichment/results/drug_immunesig_IHC_GSEA_result_files/

cancer=$1
sign=$2
num_gene=$3

# sign=all
# cancer=LIHC
# num_gene=400

resDir=/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/03_drug_immuneSig_enrichment/results/drug_immunesig_IHC_GSEA_result/${cancer}/${sign}_${num_gene}
desDir=/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/03_drug_immuneSig_enrichment/results/drug_immunesig_IHC_GSEA_result_files/${cancer}

if [ ! -d ${desDir} ]; then
    mkdir ${desDir}
    mkdir ${desDir}/${sign}_${num_gene}
fi

desDir=${desDir}/${sign}_${num_gene}

if [ ! -d ${desDir} ]; then
    mkdir ${desDir}
fi

for file in `ls ${resDir}`
do
  cp ${resDir}/${file}/treated_versus_DMSO.Gsea.*/gsea_report_for_treated_*.xls ${desDir}/gsea_report_for_treated_${file}.xls.
  cp ${resDir}/${file}/treated_versus_DMSO.Gsea.*/gsea_report_for_DMSO_*.xls ${desDir}/gsea_report_for_DMSO_${file}.xls
done



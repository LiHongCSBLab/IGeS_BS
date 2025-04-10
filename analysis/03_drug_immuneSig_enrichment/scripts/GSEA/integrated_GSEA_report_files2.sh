#!/bin/bash
mkdir /picb/bigdata/project/FengFYM/immune_targeted_combo/immunecell_DEG_cor/GSEA_GSVA/drug_cancer_intrinic_genes_result_files/

cancer=$1
sign=$2
num_gene=$3

# sign=all
# cancer=LIHC
# num_gene=400

resDir=/picb/bigdata/project/FengFYM/immune_targeted_combo/immunecell_DEG_cor/GSEA_GSVA/drug_cancer_intrinic_genes_result/${cancer}/${sign}_${num_gene}
desDir=/picb/bigdata/project/FengFYM/immune_targeted_combo/immunecell_DEG_cor/GSEA_GSVA/drug_cancer_intrinic_genes_result_files/${cancer}
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



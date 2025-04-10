#!/bin/bash
# run under miniconda environment r_env
source /home/fengfangyoumin/miniconda3/bin/activate r_env
cd /picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/03_drug_immuneSig_enrichment/scripts/

#####======================================================================#####
#####                      cancer with sampleTypes                         #####
#####======================================================================#####

# SKCM Primary Metastatic ------------------------------------------------------
cancer=SKCM

# Rscript --vanilla 01_drug_expr_preparation_update.R -c ${cancer} -d 70138 > drug_expr_prep_${cancer}_70138.log 2>&1 &
# wait

sampleType=Metastatic
tumor_purity_method=TUMERIC
# echo "start running for ${cancer} ${sampleType} ${tumor_purity_method}"
# wait
bash 04_drug_immusig_GSEA_pipeline.sh ${cancer} ${sampleType} ${tumor_purity_method} spearman 70138
bash 04_drug_immusig_GSEA_pipeline_978genes.sh ${cancer} ${sampleType} ${tumor_purity_method} spearman 70138
wait


# rm -r /picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/03_drug_immuneSig_enrichment/results_IHC/drug_immunesig_IHC_GSEA_result/70138/${cancer}
# rm -r /picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/03_drug_immuneSig_enrichment/results_TUMERIC/drug_immunesig_TUMERI_GSEA_result/70138/${cancer}
# rm -r /picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/03_drug_immuneSig_enrichment/results_CPE/drug_immunesig_CPE_GSEA_result/70138/${cancer}

# rm -r /picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/03_drug_immuneSig_enrichment/data/drug_treated_expression_matrix/70138/${cancer}
# rm -r /picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/03_drug_immuneSig_enrichment/data/drug_treated_expression_matrix/92742/${cancer}
# wait

#!/bin/bash
# run under miniconda environment r_env
source /home/fengfangyoumin/miniconda3/bin/activate r_env
cd /picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/03_drug_immuneSig_enrichment/scripts/


#####======================================================================#####
#####                      cancer with sampleTypes                         #####
#####======================================================================#####


# BRCA Primary -----------------------------------------------------------------

cancer=BRCA
# Rscript --vanilla 01_drug_expr_preparation_update.R -c ${cancer} -d 70138 > drug_expr_prep_${cancer}_70138.log 2>&1 &
# wait

sampleType=Primary
tumor_purity_method=TUMERIC
bash 04_drug_immusig_GSEA_pipeline.sh ${cancer} ${sampleType} ${tumor_purity_method} spearman 70138
wait
bash 04_drug_immusig_GSEA_pipeline_978genes.sh ${cancer} ${sampleType} ${tumor_purity_method} spearman 70138

# BRCA.Basal BRCA.Her2 BRCA.LumA BRCA.LumB BRCA.Normal -------------------------
# BRCA_sampleTypes=(BRCA.Basal BRCA.Her2 BRCA.LumA BRCA.LumB BRCA.Normal)

# for sampleType in ${BRCA_sampleTypes[@]}
# do

#     tumor_purity_method=TUMERIC
#     bash 04_drug_immusig_GSEA_pipeline.sh ${cancer} ${sampleType} ${tumor_purity_method} spearman 70138
#     wait
# done


# # rm -r /picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/03_drug_immuneSig_enrichment/results/drug_immunesig_IHC_GSEA_result/70138/${cancer}
# # rm -r /picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/03_drug_immuneSig_enrichment/data/drug_treated_expression_matrix/70138/${cancer}
# # rm -r /picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/03_drug_immuneSig_enrichment/results/drug_immunesig_IHC_GSEA_result/92742/${cancer}
# # rm -r /picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/03_drug_immuneSig_enrichment/data/drug_treated_expression_matrix/92742/${cancer}
# # wait


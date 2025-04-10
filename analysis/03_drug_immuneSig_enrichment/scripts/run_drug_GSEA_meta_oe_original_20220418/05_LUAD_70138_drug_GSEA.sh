#!/bin/bash
# run under miniconda environment r_env
source /home/fengfangyoumin/miniconda3/bin/activate r_env
cd /picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/03_drug_immuneSig_enrichment/scripts/



#####======================================================================#####
#####                      pan-cancer: Primary pcor                        #####
#####======================================================================#####


# cancer_names=(LUAD LUSC)
cancer=LUAD
sampleType=Primary
# Rscript --vanilla 01_drug_expr_preparation_update.R -c ${cancer} -d 70138 > drug_expr_prep_${cancer}_70138.log 2>&1 &
# wait
tumor_purity_method=TUMERIC
# echo "start running for ${cancer} ${sampleType} ${tumor_purity_method}"
# wait
bash 04_drug_immusig_GSEA_pipeline.sh ${cancer} ${sampleType} ${tumor_purity_method} spearman 70138
wait
bash 04_drug_immusig_GSEA_pipeline_978genes.sh ${cancer} ${sampleType} ${tumor_purity_method} spearman 70138


# cancer=LUSC
# sampleType=Primary
# # Rscript --vanilla 01_drug_expr_preparation_update.R -c ${cancer} -d 92742 > drug_expr_prep_${cancer}_92742.log 2>&1 &
# wait
# tumor_purity_method=TUMERIC
# # echo "start running for ${cancer} ${sampleType} ${tumor_purity_method}"
# # wait
# bash 04_drug_immusig_GSEA_pipeline.sh ${cancer} ${sampleType} ${tumor_purity_method} spearman 92742
# wait
# bash 04_drug_immusig_GSEA_pipeline_978genes.sh ${cancer} ${sampleType} ${tumor_purity_method} spearman 92742
# wait





# tumor_purity_method=CPE
# # echo "start running for ${cancer} ${sampleType} ${tumor_purity_method}"
    
# bash 04_drug_immusig_GSEA_pipeline_978genes.sh ${cancer} ${sampleType} ${tumor_purity_method} spearman 70138
# # wait
# bash 04_drug_immusig_GSEA_pipeline.sh ${cancer} ${sampleType} ${tumor_purity_method} spearman 70138
# wait

# tumor_purity_method=CPE
# # echo "start running for ${cancer} ${sampleType} ${tumor_purity_method}"
# bash 04_drug_immusig_GSEA_pipeline_978genes.sh ${cancer} ${sampleType} ${tumor_purity_method} spearman 92742
# # wait
# bash 04_drug_immusig_GSEA_pipeline.sh ${cancer} ${sampleType} ${tumor_purity_method} spearman 92742
# wait



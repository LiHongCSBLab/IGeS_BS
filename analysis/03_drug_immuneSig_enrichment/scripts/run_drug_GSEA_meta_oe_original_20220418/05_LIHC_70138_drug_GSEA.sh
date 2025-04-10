#!/bin/bash
# run under miniconda environment r_env
source /home/fengfangyoumin/miniconda3/bin/activate r_env
cd /picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/03_drug_immuneSig_enrichment/scripts/

wDir=/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/03_drug_immuneSig_enrichment/scripts/

# Use LIHC to test pipeline and the instructions
# this part has been run with testing_pipeline.sh.
# 01: run Rscript to generate expr file and cls file
    # parameter 1: "cancer",   "p", 1, "character", "the path of gene and fold-change",
    # parameter 2: "dataset",  "D", 1, "character", "dataset ID",

# 02 run enrichment analysis pipeline
# sample code:
# Rscript --vanilla 01_drug_expr_preparation_update.R -c SKCM -d 70138 > drug_expr_prep.log 2>&1 &

cancer=LIHC
sampleType=Primary

# Rscript --vanilla 01_drug_expr_preparation_update.R -c ${cancer} -d 70138 > drug_expr_prep_${cancer}_70138.log 2>&1 &
# wait

tumor_purity_method=TUMERIC
# # echo "start running for ${cancer} ${sampleType} ${tumor_purity_method}"
# wait
bash ${wDir}04_drug_immusig_GSEA_pipeline.sh ${cancer} ${sampleType} ${tumor_purity_method} spearman 70138
# bash ${wDir}04_drug_immusig_GSEA_pipeline_978genes.sh ${cancer} ${sampleType} ${tumor_purity_method} spearman 70138
wait


# Rscript --vanilla 01_drug_expr_preparation_update.R -c ${cancer} -d 92742 > drug_expr_prep_${cancer}_92742.log 2>&1 &
# wait

# tumor_purity_method=TUMERIC
# # # echo "start running for ${cancer} ${sampleType} ${tumor_purity_method}"
# # wait
# bash ${wDir}04_drug_immusig_GSEA_pipeline.sh ${cancer} ${sampleType} ${tumor_purity_method} spearman 92742
# bash ${wDir}04_drug_immusig_GSEA_pipeline_978genes.sh ${cancer} ${sampleType} ${tumor_purity_method} spearman 92742
# wait



# tumor_purity_method=CPE
# # # echo "start running for ${cancer} ${sampleType} ${tumor_purity_method}"
    
# bash 04_drug_immusig_GSEA_pipeline_978genes.sh ${cancer} ${sampleType} ${tumor_purity_method} spearman 70138
# # wait
# bash 04_drug_immusig_GSEA_pipeline.sh ${cancer} ${sampleType} ${tumor_purity_method} spearman 70138
# wait

# tumor_purity_method=IHC
# # # echo "start running for ${cancer} ${sampleType} ${tumor_purity_method}"
# bash 04_drug_immusig_GSEA_pipeline_978genes.sh ${cancer} ${sampleType} ${tumor_purity_method} spearman 70138
# # wait
# bash 04_drug_immusig_GSEA_pipeline.sh ${cancer} ${sampleType} ${tumor_purity_method} spearman 70138
# wait

# tumor_purity_method=CPE
# # # echo "start running for ${cancer} ${sampleType} ${tumor_purity_method}"
# bash 04_drug_immusig_GSEA_pipeline_978genes.sh ${cancer} ${sampleType} ${tumor_purity_method} spearman 92742
# # wait
# bash 04_drug_immusig_GSEA_pipeline.sh ${cancer} ${sampleType} ${tumor_purity_method} spearman 92742
# wait

# tumor_purity_method=IHC
# # # echo "start running for ${cancer} ${sampleType} ${tumor_purity_method}"
# bash 04_drug_immusig_GSEA_pipeline_978genes.sh ${cancer} ${sampleType} ${tumor_purity_method} spearman 92742
# # wait
# bash 04_drug_immusig_GSEA_pipeline.sh ${cancer} ${sampleType} ${tumor_purity_method} spearman 92742
# wait

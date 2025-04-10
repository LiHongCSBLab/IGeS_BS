#!/bin/bash
# run under miniconda environment r_env
source /home/fengfangyoumin/miniconda3/bin/activate r_env

nohup Rscript --vanilla plot_for_drug_metaAna_en_LUAD_70138.R -w model_weight -t weighted -m glmnet_weight -p 0.05 -i FALSE -s Type > meta_en_LUAD_70138.log 2>&1&
nohup Rscript --vanilla plot_for_drug_metaAna_en_othercancer_70138.R -w model_weight -t weighted -m glmnet_weight -p 0.05 -i FALSE -s Type > meta_en_othercancer_70138.log 2>&1&
nohup Rscript --vanilla plot_for_drug_metaAna_en_SKCM_Metastatic_70138.R -w model_weight -t weighted -m glmnet_weight -p 0.05 -i FALSE -s Type  > meta_en_MSKCM_70138.log 2>&1&
nohup Rscript --vanilla plot_for_drug_metaAna_en_SKCM_Primary_70138.R -w model_weight -t weighted -m glmnet_weight -p 0.05 -i FALSE -s Type  > meta_en_PSKCM_70138.log 2>&1&
# wait

# Rscript --vanilla plot_for_drug_metaAna_en_LUAD_70138.R -w meta_weight -t weighted -m sum_mean -p 0.05 -i FALSE -s Type > meta_en_LUAD_70138_metaweight_sum_mean.log 2>&1&
# Rscript --vanilla plot_for_drug_metaAna_en_othercancer_70138.R -w meta_weight -t weighted -m sum_mean -p 0.05 -i FALSE -s Type > meta_en_othercancer_70138_metaweight_sum_mean.log 2>&1&
# Rscript --vanilla plot_for_drug_metaAna_en_SKCM_Metastatic_70138.R -w meta_weight -t weighted -m sum_mean -p 0.05 -i FALSE -s Type  > meta_en_MSKCM_70138_metaweight_sum_mean.log 2>&1&
# Rscript --vanilla plot_for_drug_metaAna_en_SKCM_Primary_70138.R -w meta_weight -t weighted -m sum_mean -p 0.05 -i FALSE -s Type  > meta_en_PSKCM_70138_metaweight_sum_mean.log 2>&1&
# wait

# Rscript --vanilla plot_for_drug_metaAna_en_LUAD_70138.R -w meta_weight -t weighted -m simple -p 0.05 -i FALSE -s Type > meta_en_LUAD_70138_metaweight_simple.log 2>&1&
# Rscript --vanilla plot_for_drug_metaAna_en_othercancer_70138.R -w meta_weight -t weighted -m simple -p 0.05 -i FALSE -s Type > meta_en_othercancer_70138_metaweight_simple.log 2>&1&
# Rscript --vanilla plot_for_drug_metaAna_en_SKCM_Metastatic_70138.R -w meta_weight -t weighted -m simple -p 0.05 -i FALSE -s Type  > meta_en_MSKCM_70138_metaweight_simple.log 2>&1&
# Rscript --vanilla plot_for_drug_metaAna_en_SKCM_Primary_70138.R -w meta_weight -t weighted -m simple -p 0.05 -i FALSE -s Type  > meta_en_PSKCM_70138_metaweight_simple.log 2>&1&
# # wait

# # Rscript --vanilla plot_for_drug_metaAna_en_LUAD_70138.R -w model_weight -t weighted -m sum_mean -p 0.05 -i FALSE -s Type > meta_en_LUAD_70138_enweight_sum_mean.log 2>&1&
# # Rscript --vanilla plot_for_drug_metaAna_en_othercancer_70138.R -w model_weight -t weighted -m sum_mean -p 0.05 -i FALSE -s Type > meta_en_othercancer_70138_enweight_sum_mean.log 2>&1&
# # Rscript --vanilla plot_for_drug_metaAna_en_SKCM_Metastatic_70138.R -w model_weight -t weighted -m sum_mean -p 0.05 -i FALSE -s Type  > meta_en_MSKCM_70138_enweight_sum_mean.log 2>&1&
# # Rscript --vanilla plot_for_drug_metaAna_en_SKCM_Primary_70138.R -w model_weight -t weighted -m sum_mean -p 0.05 -i FALSE -s Type  > meta_en_PSKCM_70138_enweight_sum_mean.log 2>&1&
# # wait

# # Rscript --vanilla plot_for_drug_metaAna_en_LUAD_70138.R -w model_weight -t weighted -m simple -p 0.05 -i FALSE -s Type > meta_en_LUAD_70138_enweight_simple.log 2>&1&
# # Rscript --vanilla plot_for_drug_metaAna_en_othercancer_70138.R -w model_weight -t weighted -m simple -p 0.05 -i FALSE -s Type > meta_en_othercancer_70138_enweight_simple.log 2>&1&
# # Rscript --vanilla plot_for_drug_metaAna_en_SKCM_Metastatic_70138.R -w model_weight -t weighted -m simple -p 0.05 -i FALSE -s Type  > meta_en_MSKCM_70138_enweight_simple.log 2>&1&
# # Rscript --vanilla plot_for_drug_metaAna_en_SKCM_Primary_70138.R -w model_weight -t weighted -m simple -p 0.05 -i FALSE -s Type  > meta_en_PSKCM_70138_enweight_simple.log 2>&1&
# # wait



# # Rscript --vanilla plot_for_drug_metaAna_en_LUAD_70138.R -w model_weight -t unweighted -m simple -p 0.05 -i FALSE -s Type > meta_en_LUAD_70138_unweight_simple.log 2>&1&
# # Rscript --vanilla plot_for_drug_metaAna_en_othercancer_70138.R -w model_weight -t unweighted -m simple -p 0.05 -i FALSE -s Type > meta_en_othercancer_70138_unweight_simple.log 2>&1&
# # Rscript --vanilla plot_for_drug_metaAna_en_SKCM_Metastatic_70138.R -w model_weight -t unweighted -m simple -p 0.05 -i FALSE -s Type  > meta_en_MSKCM_70138_unweight_simple.log 2>&1&
# # Rscript --vanilla plot_for_drug_metaAna_en_SKCM_Primary_70138.R -w model_weight -t unweighted -m simple -p 0.05 -i FALSE -s Type  > meta_en_PSKCM_70138_unweight_simple.log 2>&1&
# # wait

# # Rscript --vanilla plot_for_drug_metaAna_en_LUAD_70138.R -w model_weight -t unweighted -m sum_mean -p 0.05 -i FALSE -s Type > meta_en_LUAD_70138_unweight_sum_mean.log 2>&1&
# # Rscript --vanilla plot_for_drug_metaAna_en_othercancer_70138.R -w model_weight -t unweighted -m sum_mean -p 0.05 -i FALSE -s Type > meta_en_othercancer_70138_unweight_sum_mean.log 2>&1&
# # Rscript --vanilla plot_for_drug_metaAna_en_SKCM_Metastatic_70138.R -w model_weight -t unweighted -m sum_mean -p 0.05 -i FALSE -s Type  > meta_en_MSKCM_70138_unweight_sum_mean.log 2>&1&
# # Rscript --vanilla plot_for_drug_metaAna_en_SKCM_Primary_70138.R -w model_weight -t unweighted -m sum_mean -p 0.05 -i FALSE -s Type  > meta_en_PSKCM_70138_unweight_sum_mean.log 2>&1&
# # wait


# Rscript --vanilla plot_for_drug_metaAna_en_LUAD_70138.R -w model_weight -t weighted -m glmnet_weight -p 0.05 -i TRUE -s Type> meta_en_LUAD_70138.log 2>&1&
# Rscript --vanilla plot_for_drug_metaAna_en_othercancer_70138.R -w model_weight -t weighted -m glmnet_weight -p 0.05  -i TRUE -s Type> meta_en_othercancer_70138.log 2>&1&
# Rscript --vanilla plot_for_drug_metaAna_en_SKCM_Metastatic_70138.R -w model_weight -t weighted -m glmnet_weight -p 0.05  -i TRUE -s Type > meta_en_MSKCM_70138.log 2>&1&
# Rscript --vanilla plot_for_drug_metaAna_en_SKCM_Primary_70138.R -w model_weight -t weighted -m glmnet_weight -p 0.05  -i TRUE -s Type > meta_en_PSKCM_70138.log 2>&1&
# wait

# Rscript --vanilla plot_for_drug_metaAna_en_LUAD_70138.R -w model_weight -t weighted -m sum_mean -p 0.05  -i TRUE -s Type > meta_en_LUAD_70138_enweight_sum_mean.log 2>&1&
# Rscript --vanilla plot_for_drug_metaAna_en_othercancer_70138.R -w model_weight -t weighted -m sum_mean -p 0.05 -i TRUE -s Type > meta_en_othercancer_70138_enweight_sum_mean.log 2>&1&
# Rscript --vanilla plot_for_drug_metaAna_en_SKCM_Metastatic_70138.R -w model_weight -t weighted -m sum_mean -p 0.05  -i TRUE -s Type > meta_en_MSKCM_70138_enweight_sum_mean.log 2>&1&
# Rscript --vanilla plot_for_drug_metaAna_en_SKCM_Primary_70138.R -w model_weight -t weighted -m sum_mean -p 0.05  -i TRUE -s Type > meta_en_PSKCM_70138_enweight_sum_mean.log 2>&1&
# wait

# Rscript --vanilla plot_for_drug_metaAna_en_LUAD_70138.R -w model_weight -t weighted -m simple -p 0.05 -i TRUE -s Type > meta_en_LUAD_70138_enweight_simple.log 2>&1&
# Rscript --vanilla plot_for_drug_metaAna_en_othercancer_70138.R -w model_weight -t weighted -m simple -p 0.05 -i TRUE -s Type > meta_en_othercancer_70138_enweight_simple.log 2>&1&
# Rscript --vanilla plot_for_drug_metaAna_en_SKCM_Metastatic_70138.R -w model_weight -t weighted -m simple -p 0.05 -i TRUE -s Type  > meta_en_MSKCM_70138_enweight_simple.log 2>&1&
# Rscript --vanilla plot_for_drug_metaAna_en_SKCM_Primary_70138.R -w model_weight -t weighted -m simple -p 0.05 -i TRUE -s Type  > meta_en_PSKCM_70138_enweight_simple.log 2>&1&
# wait


# Rscript --vanilla plot_for_drug_metaAna_en_LUAD_70138.R -w meta_weight -t weighted -m sum_mean -p 0.05 -i TRUE -s Type > meta_en_LUAD_70138_metaweight_sum_mean.log 2>&1&
# Rscript --vanilla plot_for_drug_metaAna_en_othercancer_70138.R -w meta_weight -t weighted -m sum_mean -p 0.05 -i TRUE -s Type > meta_en_othercancer_70138_metaweight_sum_mean.log 2>&1& 
# Rscript --vanilla plot_for_drug_metaAna_en_SKCM_Metastatic_70138.R -w meta_weight -t weighted -m sum_mean -p 0.05 -i TRUE -s Type  > meta_en_MSKCM_70138_metaweight_sum_mean.log 2>&1&
# Rscript --vanilla plot_for_drug_metaAna_en_SKCM_Primary_70138.R -w meta_weight -t weighted -m sum_mean -p 0.05 -i TRUE -s Type  > meta_en_PSKCM_70138_metaweight_sum_mean.log 2>&1&
# wait

# Rscript --vanilla plot_for_drug_metaAna_en_LUAD_70138.R -w meta_weight -t weighted -m simple -p 0.05 -i TRUE -s Type > meta_en_LUAD_70138_metaweight_simple.log 2>&1&
# Rscript --vanilla plot_for_drug_metaAna_en_othercancer_70138.R -w meta_weight -t weighted -m simple -p 0.05 -i TRUE -s Type > meta_en_othercancer_70138_metaweight_simple.log 2>&1&
# Rscript --vanilla plot_for_drug_metaAna_en_SKCM_Metastatic_70138.R -w meta_weight -t weighted -m simple -p 0.05 -i TRUE -s Type  > meta_en_MSKCM_70138_metaweight_simple.log 2>&1&
# Rscript --vanilla plot_for_drug_metaAna_en_SKCM_Primary_70138.R -w meta_weight -t weighted -m simple -p 0.05 -i TRUE -s Type  > meta_en_PSKCM_70138_metaweight_simple.log 2>&1&
# wait

# Rscript --vanilla plot_for_drug_metaAna_en_LUAD_70138.R -w model_weight -t unweighted -m simple -p 0.05 -i TRUE -s Type > meta_en_LUAD_70138_unweight_simple.log 2>&1&
# Rscript --vanilla plot_for_drug_metaAna_en_othercancer_70138.R -w model_weight -t unweighted -m simple -p 0.05 -i TRUE -s Type > meta_en_othercancer_70138_unweight_simple.log 2>&1&
# Rscript --vanilla plot_for_drug_metaAna_en_SKCM_Metastatic_70138.R -w model_weight -t unweighted -m simple -p 0.05 -i TRUE -s Type  > meta_en_MSKCM_70138_unweight_simple.log 2>&1&
# Rscript --vanilla plot_for_drug_metaAna_en_SKCM_Primary_70138.R -w model_weight -t unweighted -m simple -p 0.05 -i TRUE -s Type  > meta_en_PSKCM_70138_unweight_simple.log 2>&1&
# wait

# Rscript --vanilla plot_for_drug_metaAna_en_LUAD_70138.R -w model_weight -t unweighted -m sum_mean -p 0.05 -i TRUE -s Type > meta_en_LUAD_70138_unweight_sum_mean.log 2>&1&
# Rscript --vanilla plot_for_drug_metaAna_en_othercancer_70138.R -w model_weight -t unweighted -m sum_mean -p 0.05 -i TRUE -s Type > meta_en_othercancer_70138_unweight_sum_mean.log 2>&1&
# Rscript --vanilla plot_for_drug_metaAna_en_SKCM_Metastatic_70138.R -w model_weight -t unweighted -m sum_mean -p 0.05 -i TRUE -s Type  > meta_en_MSKCM_70138_unweight_sum_mean.log 2>&1&
# Rscript --vanilla plot_for_drug_metaAna_en_SKCM_Primary_70138.R -w model_weight -t unweighted -m sum_mean -p 0.05 -i TRUE -s Type  > meta_en_PSKCM_70138_unweight_sum_mean.log 2>&1&
# wait


# Rscript --vanilla plot_for_drug_metaAna_en_LUAD_70138.R -w model_weight -t weighted -m glmnet_weight -p 0.05 -i TRUE -s Index> meta_en_LUAD_70138.log 2>&1&
# Rscript --vanilla plot_for_drug_metaAna_en_othercancer_70138.R -w model_weight -t weighted -m glmnet_weight -p 0.05  -i TRUE -s Index> meta_en_othercancer_70138.log 2>&1&
# Rscript --vanilla plot_for_drug_metaAna_en_SKCM_Metastatic_70138.R -w model_weight -t weighted -m glmnet_weight -p 0.05  -i TRUE -s Index > meta_en_MSKCM_70138.log 2>&1&
# Rscript --vanilla plot_for_drug_metaAna_en_SKCM_Primary_70138.R -w model_weight -t weighted -m glmnet_weight -p 0.05  -i TRUE -s Index > meta_en_PSKCM_70138.log 2>&1&
# wait

# Rscript --vanilla plot_for_drug_metaAna_en_LUAD_70138.R -w model_weight -t weighted -m sum_mean -p 0.05  -i TRUE -s Index > meta_en_LUAD_70138_enweight_sum_mean.log 2>&1&
# Rscript --vanilla plot_for_drug_metaAna_en_othercancer_70138.R -w model_weight -t weighted -m sum_mean -p 0.05 -i TRUE -s Index > meta_en_othercancer_70138_enweight_sum_mean.log 2>&1&
# Rscript --vanilla plot_for_drug_metaAna_en_SKCM_Metastatic_70138.R -w model_weight -t weighted -m sum_mean -p 0.05  -i TRUE -s Index > meta_en_MSKCM_70138_enweight_sum_mean.log 2>&1&
# Rscript --vanilla plot_for_drug_metaAna_en_SKCM_Primary_70138.R -w model_weight -t weighted -m sum_mean -p 0.05  -i TRUE -s Index > meta_en_PSKCM_70138_enweight_sum_mean.log 2>&1&
# wait

# Rscript --vanilla plot_for_drug_metaAna_en_LUAD_70138.R -w model_weight -t weighted -m simple -p 0.05 -i TRUE -s Index > meta_en_LUAD_70138_enweight_simple.log 2>&1&
# Rscript --vanilla plot_for_drug_metaAna_en_othercancer_70138.R -w model_weight -t weighted -m simple -p 0.05 -i TRUE -s Index > meta_en_othercancer_70138_enweight_simple.log 2>&1&
# Rscript --vanilla plot_for_drug_metaAna_en_SKCM_Metastatic_70138.R -w model_weight -t weighted -m simple -p 0.05 -i TRUE -s Index  > meta_en_MSKCM_70138_enweight_simple.log 2>&1&
# Rscript --vanilla plot_for_drug_metaAna_en_SKCM_Primary_70138.R -w model_weight -t weighted -m simple -p 0.05 -i TRUE -s Index  > meta_en_PSKCM_70138_enweight_simple.log 2>&1&
# wait


# Rscript --vanilla plot_for_drug_metaAna_en_LUAD_70138.R -w meta_weight -t weighted -m sum_mean -p 0.05 -i TRUE -s Index > meta_en_LUAD_70138_metaweight_sum_mean.log 2>&1&
# Rscript --vanilla plot_for_drug_metaAna_en_othercancer_70138.R -w meta_weight -t weighted -m sum_mean -p 0.05 -i TRUE -s Index > meta_en_othercancer_70138_metaweight_sum_mean.log 2>&1&
# Rscript --vanilla plot_for_drug_metaAna_en_SKCM_Metastatic_70138.R -w meta_weight -t weighted -m sum_mean -p 0.05 -i TRUE -s Index  > meta_en_MSKCM_70138_metaweight_sum_mean.log 2>&1&
# Rscript --vanilla plot_for_drug_metaAna_en_SKCM_Primary_70138.R -w meta_weight -t weighted -m sum_mean -p 0.05 -i TRUE -s Index  > meta_en_PSKCM_70138_metaweight_sum_mean.log 2>&1&
# wait

# Rscript --vanilla plot_for_drug_metaAna_en_LUAD_70138.R -w meta_weight -t weighted -m simple -p 0.05 -i TRUE -s Index > meta_en_LUAD_70138_metaweight_simple.log 2>&1&
# Rscript --vanilla plot_for_drug_metaAna_en_othercancer_70138.R -w meta_weight -t weighted -m simple -p 0.05 -i TRUE -s Index > meta_en_othercancer_70138_metaweight_simple.log 2>&1&
# Rscript --vanilla plot_for_drug_metaAna_en_SKCM_Metastatic_70138.R -w meta_weight -t weighted -m simple -p 0.05 -i TRUE -s Index  > meta_en_MSKCM_70138_metaweight_simple.log 2>&1&
# Rscript --vanilla plot_for_drug_metaAna_en_SKCM_Primary_70138.R -w meta_weight -t weighted -m simple -p 0.05 -i TRUE -s Index  > meta_en_PSKCM_70138_metaweight_simple.log 2>&1&
# wait

# Rscript --vanilla plot_for_drug_metaAna_en_LUAD_70138.R -w model_weight -t unweighted -m simple -p 0.05 -i TRUE -s Index > meta_en_LUAD_70138_unweight_simple.log 2>&1&
# Rscript --vanilla plot_for_drug_metaAna_en_othercancer_70138.R -w model_weight -t unweighted -m simple -p 0.05 -i TRUE -s Index > meta_en_othercancer_70138_unweight_simple.log 2>&1&
# Rscript --vanilla plot_for_drug_metaAna_en_SKCM_Metastatic_70138.R -w model_weight -t unweighted -m simple -p 0.05 -i TRUE -s Index  > meta_en_MSKCM_70138_unweight_simple.log 2>&1&
# Rscript --vanilla plot_for_drug_metaAna_en_SKCM_Primary_70138.R -w model_weight -t unweighted -m simple -p 0.05 -i TRUE -s Index  > meta_en_PSKCM_70138_unweight_simple.log 2>&1&
# wait

# Rscript --vanilla plot_for_drug_metaAna_en_LUAD_70138.R -w model_weight -t unweighted -m sum_mean -p 0.05 -i TRUE -s Index > meta_en_LUAD_70138_unweight_sum_mean.log 2>&1&
# Rscript --vanilla plot_for_drug_metaAna_en_othercancer_70138.R -w model_weight -t unweighted -m sum_mean -p 0.05 -i TRUE -s Index > meta_en_othercancer_70138_unweight_sum_mean.log 2>&1&
# Rscript --vanilla plot_for_drug_metaAna_en_SKCM_Metastatic_70138.R -w model_weight -t unweighted -m sum_mean -p 0.05 -i TRUE -s Index  > meta_en_MSKCM_70138_unweight_sum_mean.log 2>&1&
# Rscript --vanilla plot_for_drug_metaAna_en_SKCM_Primary_70138.R -w model_weight -t unweighted -m sum_mean -p 0.05 -i TRUE -s Index  > meta_en_PSKCM_70138_unweight_sum_mean.log 2>&1&

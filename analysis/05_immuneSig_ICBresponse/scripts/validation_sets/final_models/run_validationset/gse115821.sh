Rscript --vanilla dataset_oriSig_validationset.R -s meta_mixed -f ranger_impurity -k 0.1 -d GSE115821_pre_outcome > logfiles/mixed_ranger_impurity_0.1_GSE115821.log 2>&1&
Rscript --vanilla dataset_oriSig_validationset.R -s meta_mixed -f ranger_impurity -k 0.2 -d GSE115821_pre_outcome > logfiles/mixed_ranger_impurity_0.2_GSE115821.log 2>&1&
Rscript --vanilla dataset_oriSig_validationset.R -s meta_mixed -f ranger_impurity -k 0.3 -d GSE115821_pre_outcome > logfiles/mixed_ranger_impurity_0.3_GSE115821.log 2>&1&
Rscript --vanilla dataset_oriSig_validationset.R -s meta_mixed -f ranger_impurity -k 0.5 -d GSE115821_pre_outcome > logfiles/mixed_ranger_impurity_0.5_GSE115821.log 2>&1&


Rscript --vanilla dataset_zscaled_oriSig_validationset.R -s meta_mixed -f ranger_impurity -k 0.1 -d GSE115821_pre_outcome > logfiles/mixed_ranger_impurity_0.1_GSE96619.log 2>&1&
Rscript --vanilla dataset_zscaled_oriSig_validationset.R -s meta_mixed -f ranger_impurity -k 0.2 -d GSE115821_pre_outcome > logfiles/mixed_ranger_impurity_0.2_GSE96619.log 2>&1&
Rscript --vanilla dataset_zscaled_oriSig_validationset.R -s meta_mixed -f ranger_impurity -k 0.3 -d GSE115821_pre_outcome > logfiles/mixed_ranger_impurity_0.3_GSE96619.log 2>&1&
Rscript --vanilla dataset_zscaled_oriSig_validationset.R -s meta_mixed -f ranger_impurity -k 0.5 -d GSE115821_pre_outcome > logfiles/mixed_ranger_impurity_0.5_GSE96619.log 2>&1&

Rscript --vanilla cancer_zscale_oriSig_validationset.R -s meta_mixed -f ranger_impurity -k 0.1 -d GSE115821_pre_outcome -c SKCM > logfiles/mixed_ranger_impurity_0.1_GSE96619_zscaled_SKCM.log >&1&
Rscript --vanilla cancer_zscale_oriSig_validationset.R -s meta_mixed -f ranger_impurity -k 0.2 -d GSE115821_pre_outcome -c SKCM > logfiles/mixed_ranger_impurity_0.2_GSE96619_zscaled_SKCM.log >&1&
Rscript --vanilla cancer_zscale_oriSig_validationset.R -s meta_mixed -f ranger_impurity -k 0.3 -d GSE115821_pre_outcome -c SKCM > logfiles/mixed_ranger_impurity_0.3_GSE96619_zscaled_SKCM.log >&1&
Rscript --vanilla cancer_zscale_oriSig_validationset.R -s meta_mixed -f ranger_impurity -k 0.5 -d GSE115821_pre_outcome -c SKCM > logfiles/mixed_ranger_impurity_0.5_GSE96619_zscaled_SKCM.log >&1&

Rscript --vanilla cancer_oriSig_validationset.R -s meta_mixed -f ranger_impurity -k 0.1 -d GSE115821_pre_outcome -c SKCM > logfiles/mixed_ranger_impurity_0.1_GSE96619_SKCM.log >&1&
Rscript --vanilla cancer_oriSig_validationset.R -s meta_mixed -f ranger_impurity -k 0.2 -d GSE115821_pre_outcome -c SKCM > logfiles/mixed_ranger_impurity_0.2_GSE96619_SKCM.log >&1&
Rscript --vanilla cancer_oriSig_validationset.R -s meta_mixed -f ranger_impurity -k 0.3 -d GSE115821_pre_outcome -c SKCM > logfiles/mixed_ranger_impurity_0.3_GSE96619_SKCM.log >&1&
Rscript --vanilla cancer_oriSig_validationset.R -s meta_mixed -f ranger_impurity -k 0.5 -d GSE115821_pre_outcome -c SKCM > logfiles/mixed_ranger_impurity_0.5_GSE96619_SKCM.log >&1&

# Rscript --vanilla dataset_oriSig_validationset.R -s mlrFilter -f ranger_impurity -k 0.1 -d GSE115821_pre_outcome > logfiles/mlrFilter_ranger_impurity_0.1_GSE115821.log 2>&1&
# Rscript --vanilla dataset_oriSig_validationset.R -s mlrFilter -f ranger_impurity -k 0.2 -d GSE115821_pre_outcome > logfiles/mlrFilter_ranger_impurity_0.2_GSE115821.log 2>&1&
# # Rscript --vanilla dataset_oriSig_validationset.R -s mlrFilter -f ranger_impurity -k 0.3 -d GSE115821_pre_outcome > logfiles/mlrFilter_ranger_impurity_0.3_GSE115821.log 2>&1&

# Rscript --vanilla dataset_oriSig_validationset.R -s mixed -f ranger_impurity -k 0.1 -d GSE115821_pre_outcome > logfiles/mixed_ranger_impurity_0.1_GSE115821.log 2>&1&
# Rscript --vanilla dataset_oriSig_validationset.R -s mixed -f ranger_impurity -k 0.2 -d GSE115821_pre_outcome > logfiles/mixed_ranger_impurity_0.2_GSE115821.log 2>&1&
# # Rscript --vanilla dataset_oriSig_validationset.R -s mixed -f ranger_impurity -k 0.3 -d GSE115821_pre_outcome > logfiles/mixed_ranger_impurity_0.3_GSE115821.log 2>&1&

# Rscript --vanilla dataset_zscaled_oriSig_validationset.R -s mlrFilter -f ranger_impurity -k 0.1 -d GSE115821_pre_outcome > logfiles/mlrFilter_ranger_impurity_0.1_GSE115821_zscaled.log 2>&1&
# Rscript --vanilla dataset_zscaled_oriSig_validationset.R -s mlrFilter -f ranger_impurity -k 0.2 -d GSE115821_pre_outcome > logfiles/mlrFilter_ranger_impurity_0.2_GSE115821_zscaled.log 2>&1&
# # Rscript --vanilla dataset_zscaled_oriSig_validationset.R -s mlrFilter -f ranger_impurity -k 0.3 -d GSE115821_pre_outcome > logfiles/mlrFilter_ranger_impurity_0.3_GSE115821_zscaled.log 2>&1&

# Rscript --vanilla dataset_zscaled_oriSig_validationset.R -s mixed -f ranger_impurity -k 0.1 -d GSE115821_pre_outcome > logfiles/mixed_ranger_impurity_0.1_GSE115821_zscaled.log 2>&1&
# Rscript --vanilla dataset_zscaled_oriSig_validationset.R -s mixed -f ranger_impurity -k 0.2 -d GSE115821_pre_outcome > logfiles/mixed_ranger_impurity_0.2_GSE115821_zscaled.log 2>&1&
# # Rscript --vanilla dataset_zscaled_oriSig_validationset.R -s mixed -f ranger_impurity -k 0.3 -d GSE115821_pre_outcome > logfiles/mixed_ranger_impurity_0.3_GSE115821_zscaled.log 2>&1&
# wait
# Rscript --vanilla cancer_zscale_oriSig_validationset.R -s mlrFilter -f ranger_impurity -k 0.1 -d GSE115821_pre_outcome -c SKCM > logfiles/mlrFilter_ranger_impurity_0.1_GSE115821_zscaled_SKCM.log >&1&
# Rscript --vanilla cancer_zscale_oriSig_validationset.R -s mlrFilter -f ranger_impurity -k 0.2 -d GSE115821_pre_outcome -c SKCM > logfiles/mlrFilter_ranger_impurity_0.2_GSE115821_zscaled_SKCM.log >&1&
# # Rscript --vanilla cancer_zscale_oriSig_validationset.R -s mlrFilter -f ranger_impurity -k 0.3 -d GSE115821_pre_outcome -c SKCM > logfiles/mlrFilter_ranger_impurity_0.3_GSE115821_zscaled_SKCM.log >&1&

# Rscript --vanilla cancer_zscale_oriSig_validationset.R -s mixed -f ranger_impurity -k 0.1 -d GSE115821_pre_outcome -c SKCM > logfiles/mixed_ranger_impurity_0.1_GSE115821_zscaled_SKCM.log >&1&
# Rscript --vanilla cancer_zscale_oriSig_validationset.R -s mixed -f ranger_impurity -k 0.2 -d GSE115821_pre_outcome -c SKCM > logfiles/mixed_ranger_impurity_0.2_GSE115821_zscaled_SKCM.log >&1&
# # Rscript --vanilla cancer_zscale_oriSig_validationset.R -s mixed -f ranger_impurity -k 0.3 -d GSE115821_pre_outcome -c SKCM > logfiles/mixed_ranger_impurity_0.3_GSE115821_zscaled_SKCM.log >&1&

# Rscript --vanilla cancer_oriSig_validationset.R -s mlrFilter -f ranger_impurity -k 0.1 -d GSE115821_pre_outcome -c SKCM > logfiles/mlrFilter_ranger_impurity_0.1_GSE115821_SKCM.log >&1&
# Rscript --vanilla cancer_oriSig_validationset.R -s mlrFilter -f ranger_impurity -k 0.2 -d GSE115821_pre_outcome -c SKCM > logfiles/mlrFilter_ranger_impurity_0.2_GSE115821_SKCM.log >&1&
# # Rscript --vanilla cancer_oriSig_validationset.R -s mlrFilter -f ranger_impurity -k 0.3 -d GSE115821_pre_outcome -c SKCM > logfiles/mlrFilter_ranger_impurity_0.3_GSE115821_SKCM.log >&1&

# Rscript --vanilla cancer_oriSig_validationset.R -s mixed -f ranger_impurity -k 0.1 -d GSE115821_pre_outcome -c SKCM > logfiles/mixed_ranger_impurity_0.1_GSE115821_SKCM.log >&1&
# Rscript --vanilla cancer_oriSig_validationset.R -s mixed -f ranger_impurity -k 0.2 -d GSE115821_pre_outcome -c SKCM > logfiles/mixed_ranger_impurity_0.2_GSE115821_SKCM.log >&1&
# # Rscript --vanilla cancer_oriSig_validationset.R -s mixed -f ranger_impurity -k 0.3 -d GSE115821_pre_outcome -c SKCM > logfiles/mixed_ranger_impurity_0.3_GSE115821_SKCM.log >&1&


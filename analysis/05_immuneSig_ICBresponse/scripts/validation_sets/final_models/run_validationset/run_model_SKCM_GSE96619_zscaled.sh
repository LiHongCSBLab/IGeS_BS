#! /bash/bin

Rscript --vanilla cancer_zscale_oriSig_validationset.R -s mlrFilter -f ranger_impurity -k 0.1 -d GSE96619_outcome -c SKCM > logfiles/mlrFilter_ranger_impurity_0.1_GSE96619_zscaled_SKCM.log >&1&
Rscript --vanilla cancer_zscale_oriSig_validationset.R -s mlrFilter -f ranger_impurity -k 0.3 -d GSE96619_outcome -c SKCM > logfiles/mlrFilter_ranger_impurity_0.3_GSE96619_zscaled_SKCM.log >&1&
Rscript --vanilla cancer_zscale_oriSig_validationset.R -s mlrFilter -f ranger_impurity -k 0.5 -d GSE96619_outcome -c SKCM > logfiles/mlrFilter_ranger_impurity_0.5_GSE96619_zscaled_SKCM.log >&1&
Rscript --vanilla cancer_zscale_oriSig_validationset.R -s mlrFilter -f ranger_impurity -k 0.7 -d GSE96619_outcome -c SKCM > logfiles/mlrFilter_ranger_impurity_0.7_GSE96619_zscaled_SKCM.log >&1&
Rscript --vanilla cancer_zscale_oriSig_validationset.R -s mlrFilter -f ranger_impurity -k 0.9 -d GSE96619_outcome -c SKCM > logfiles/mlrFilter_ranger_impurity_0.9_GSE96619_zscaled_SKCM.log >&1&

Rscript --vanilla cancer_zscale_oriSig_validationset.R -s mixed -f ranger_impurity -k 0.1 -d GSE96619_outcome -c SKCM > logfiles/mixed_ranger_impurity_0.1_GSE96619_zscaled_SKCM.log >&1&
Rscript --vanilla cancer_zscale_oriSig_validationset.R -s mixed -f ranger_impurity -k 0.3 -d GSE96619_outcome -c SKCM > logfiles/mixed_ranger_impurity_0.3_GSE96619_zscaled_SKCM.log >&1&
Rscript --vanilla cancer_zscale_oriSig_validationset.R -s mixed -f ranger_impurity -k 0.5 -d GSE96619_outcome -c SKCM > logfiles/mixed_ranger_impurity_0.5_GSE96619_zscaled_SKCM.log >&1&
Rscript --vanilla cancer_zscale_oriSig_validationset.R -s mixed -f ranger_impurity -k 0.7 -d GSE96619_outcome -c SKCM > logfiles/mixed_ranger_impurity_0.7_GSE96619_zscaled_SKCM.log >&1&
Rscript --vanilla cancer_zscale_oriSig_validationset.R -s mixed -f ranger_impurity -k 0.9 -d GSE96619_outcome -c SKCM > logfiles/mixed_ranger_impurity_0.9_GSE96619_zscaled_SKCM.log >&1&

wait

Rscript --vanilla cancer_zscale_oriSig_validationset.R -s metaAnalysis -d GSE96619_outcome -c SKCM > logfiles/metaAnalysis_GSE96619_zscaled_SKCM.log >&1&
Rscript --vanilla cancer_zscale_oriSig_validationset.R -s combinedP -d GSE96619_outcome -c SKCM > logfiles/combinedP_GSE96619_zscaled_SKCM.log >&1&
Rscript --vanilla cancer_zscale_oriSig_validationset.R -s voting -d GSE96619_outcome -c SKCM > logfiles/voting_GSE96619_zscaled_SKCM.log >&1&

Rscript --vanilla cancer_zscale_oriSig_validationset.R -s mlrFilter -f FSelector_symmetrical.uncertainty -d GSE96619_outcome -c SKCM > logfiles/mlrFilter_FSelector_symmetrical.uncertainty_GSE96619_zscaled_SKCM.log >&1&
wait

Rscript --vanilla cancer_zscale_oriSig_validationset.R -s mlrFilter -f auc -k 0.1 -d GSE96619_outcome -c SKCM > logfiles/mlrFilter_auc_0.1_GSE96619_zscaled_SKCM.log >&1&
Rscript --vanilla cancer_zscale_oriSig_validationset.R -s mlrFilter -f auc -k 0.3 -d GSE96619_outcome -c SKCM > logfiles/mlrFilter_auc_0.3_GSE96619_zscaled_SKCM.log >&1&
Rscript --vanilla cancer_zscale_oriSig_validationset.R -s mlrFilter -f auc -k 0.5 -d GSE96619_outcome -c SKCM > logfiles/mlrFilter_auc_0.5_GSE96619_zscaled_SKCM.log >&1&
Rscript --vanilla cancer_zscale_oriSig_validationset.R -s mlrFilter -f auc -k 0.7 -d GSE96619_outcome -c SKCM > logfiles/mlrFilter_auc_0.7_GSE96619_zscaled_SKCM.log >&1&
Rscript --vanilla cancer_zscale_oriSig_validationset.R -s mlrFilter -f auc -k 0.9 -d GSE96619_outcome -c SKCM > logfiles/mlrFilter_auc_0.9_GSE96619_zscaled_SKCM.log >&1&
wait

Rscript --vanilla cancer_zscale_oriSig_validationset.R -s mlrFilter -f praznik_JMIM -k 0.1 -d GSE96619_outcome -c SKCM > logfiles/mlrFilter_praznik_JMIM_0.1_GSE96619_zscaled_SKCM.log >&1&
Rscript --vanilla cancer_zscale_oriSig_validationset.R -s mlrFilter -f praznik_JMIM -k 0.3 -d GSE96619_outcome -c SKCM > logfiles/mlrFilter_praznik_JMIM_0.3_GSE96619_zscaled_SKCM.log >&1&
Rscript --vanilla cancer_zscale_oriSig_validationset.R -s mlrFilter -f praznik_JMIM -k 0.5 -d GSE96619_outcome -c SKCM > logfiles/mlrFilter_praznik_JMIM_0.5_GSE96619_zscaled_SKCM.log >&1&
Rscript --vanilla cancer_zscale_oriSig_validationset.R -s mlrFilter -f praznik_JMIM -k 0.7 -d GSE96619_outcome -c SKCM > logfiles/mlrFilter_praznik_JMIM_0.7_GSE96619_zscaled_SKCM.log >&1&
Rscript --vanilla cancer_zscale_oriSig_validationset.R -s mlrFilter -f praznik_JMIM -k 0.9 -d GSE96619_outcome -c SKCM > logfiles/mlrFilter_praznik_JMIM_0.9_GSE96619_zscaled_SKCM.log >&1&
wait

Rscript --vanilla cancer_zscale_oriSig_validationset.R -s mlrFilter -f praznik_MRMR -k 0.1 -d GSE96619_outcome -c SKCM > logfiles/mlrFilter_praznik_MRMR_0.1_GSE96619_zscaled_SKCM.log >&1&
Rscript --vanilla cancer_zscale_oriSig_validationset.R -s mlrFilter -f praznik_MRMR -k 0.3 -d GSE96619_outcome -c SKCM > logfiles/mlrFilter_praznik_MRMR_0.3_GSE96619_zscaled_SKCM.log >&1&
Rscript --vanilla cancer_zscale_oriSig_validationset.R -s mlrFilter -f praznik_MRMR -k 0.5 -d GSE96619_outcome -c SKCM > logfiles/mlrFilter_praznik_MRMR_0.5_GSE96619_zscaled_SKCM.log >&1&
Rscript --vanilla cancer_zscale_oriSig_validationset.R -s mlrFilter -f praznik_MRMR -k 0.7 -d GSE96619_outcome -c SKCM > logfiles/mlrFilter_praznik_MRMR_0.7_GSE96619_zscaled_SKCM.log >&1&
Rscript --vanilla cancer_zscale_oriSig_validationset.R -s mlrFilter -f praznik_MRMR -k 0.9 -d GSE96619_outcome -c SKCM > logfiles/mlrFilter_praznik_MRMR_0.9_GSE96619_zscaled_SKCM.log >&1&

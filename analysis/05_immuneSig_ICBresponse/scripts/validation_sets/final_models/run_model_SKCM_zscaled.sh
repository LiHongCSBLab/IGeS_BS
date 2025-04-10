#! /bash/bin


Rscript --vanilla cancer_zscale_oriSig_merge_prediction_mlr_finalmodel.R -s mlrFilter -f ranger_impurity -k 0.1 -c SKCM > logfiles/mlrFilter_ranger_impurity_0.1_zscaled_SKCM.log >&1&
Rscript --vanilla cancer_zscale_oriSig_merge_prediction_mlr_finalmodel.R -s mlrFilter -f ranger_impurity -k 0.2 -c SKCM > logfiles/mlrFilter_ranger_impurity_0.2_zscaled_SKCM.log >&1&
Rscript --vanilla cancer_zscale_oriSig_merge_prediction_mlr_finalmodel.R -s mlrFilter -f ranger_impurity -k 0.3 -c SKCM > logfiles/mlrFilter_ranger_impurity_0.3_zscaled_SKCM.log >&1&
Rscript --vanilla cancer_zscale_oriSig_merge_prediction_mlr_finalmodel.R -s mlrFilter -f ranger_impurity -k 0.5 -c SKCM > logfiles/mlrFilter_ranger_impurity_0.5_zscaled_SKCM.log >&1&
# Rscript --vanilla cancer_zscale_oriSig_merge_prediction_mlr_finalmodel.R -s mlrFilter -f ranger_impurity -k 0.7 -c SKCM > logfiles/mlrFilter_ranger_impurity_0.7_zscaled_SKCM.log >&1&
# Rscript --vanilla cancer_zscale_oriSig_merge_prediction_mlr_finalmodel.R -s mlrFilter -f ranger_impurity -k 0.9 -c SKCM > logfiles/mlrFilter_ranger_impurity_0.9_zscaled_SKCM.log >&1&

Rscript --vanilla cancer_zscale_oriSig_merge_prediction_mlr_finalmodel.R -s mixed -f ranger_impurity -k 0.1 -c SKCM > logfiles/mixed_ranger_impurity_0.1_zscaled_SKCM.log >&1&
Rscript --vanilla cancer_zscale_oriSig_merge_prediction_mlr_finalmodel.R -s mixed -f ranger_impurity -k 0.2 -c SKCM > logfiles/mixed_ranger_impurity_0.2_zscaled_SKCM.log >&1&
Rscript --vanilla cancer_zscale_oriSig_merge_prediction_mlr_finalmodel.R -s mixed -f ranger_impurity -k 0.3 -c SKCM > logfiles/mixed_ranger_impurity_0.3_zscaled_SKCM.log >&1&
Rscript --vanilla cancer_zscale_oriSig_merge_prediction_mlr_finalmodel.R -s mixed -f ranger_impurity -k 0.5 -c SKCM > logfiles/mixed_ranger_impurity_0.5_zscaled_SKCM.log >&1&
# Rscript --vanilla cancer_zscale_oriSig_merge_prediction_mlr_finalmodel.R -s mixed -f ranger_impurity -k 0.7 -c SKCM > logfiles/mixed_ranger_impurity_0.7_zscaled_SKCM.log >&1&
# Rscript --vanilla cancer_zscale_oriSig_merge_prediction_mlr_finalmodel.R -s mixed -f ranger_impurity -k 0.9 -c SKCM > logfiles/mixed_ranger_impurity_0.9_zscaled_SKCM.log >&1&
wait

Rscript --vanilla cancer_zscale_oriSig_merge_prediction_mlr_finalmodel.R -s metaAnalysis -c SKCM > logfiles/metaAnalysis_zscaled_SKCM.log >&1&
Rscript --vanilla cancer_zscale_oriSig_merge_prediction_mlr_finalmodel.R -s combinedP -c SKCM > logfiles/combinedP_zscaled_SKCM.log >&1&
Rscript --vanilla cancer_zscale_oriSig_merge_prediction_mlr_finalmodel.R -s voting -c SKCM > logfiles/voting_zscaled_SKCM.log >&1&

Rscript --vanilla cancer_zscale_oriSig_merge_prediction_mlr_finalmodel.R -s mlrFilter -f FSelector_symmetrical.uncertainty -c SKCM > logfiles/mlrFilter_FSelector_symmetrical.uncertainty_zscaled_SKCM.log >&1&
wait

Rscript --vanilla cancer_zscale_oriSig_merge_prediction_mlr_finalmodel.R -s mlrFilter -f auc -k 0.1 -c SKCM > logfiles/mlrFilter_auc_0.1_zscaled_SKCM.log >&1&
Rscript --vanilla cancer_zscale_oriSig_merge_prediction_mlr_finalmodel.R -s mlrFilter -f auc -k 0.3 -c SKCM > logfiles/mlrFilter_auc_0.3_zscaled_SKCM.log >&1&
Rscript --vanilla cancer_zscale_oriSig_merge_prediction_mlr_finalmodel.R -s mlrFilter -f auc -k 0.5 -c SKCM > logfiles/mlrFilter_auc_0.5_zscaled_SKCM.log >&1&
# Rscript --vanilla cancer_zscale_oriSig_merge_prediction_mlr_finalmodel.R -s mlrFilter -f auc -k 0.7 -c SKCM > logfiles/mlrFilter_auc_0.7_zscaled_SKCM.log >&1&
# Rscript --vanilla cancer_zscale_oriSig_merge_prediction_mlr_finalmodel.R -s mlrFilter -f auc -k 0.9 -c SKCM > logfiles/mlrFilter_auc_0.9_zscaled_SKCM.log >&1&
wait

Rscript --vanilla cancer_zscale_oriSig_merge_prediction_mlr_finalmodel.R -s mlrFilter -f praznik_JMIM -k 0.1 -c SKCM > logfiles/mlrFilter_praznik_JMIM_0.1_zscaled_SKCM.log >&1&
Rscript --vanilla cancer_zscale_oriSig_merge_prediction_mlr_finalmodel.R -s mlrFilter -f praznik_JMIM -k 0.3 -c SKCM > logfiles/mlrFilter_praznik_JMIM_0.3_zscaled_SKCM.log >&1&
Rscript --vanilla cancer_zscale_oriSig_merge_prediction_mlr_finalmodel.R -s mlrFilter -f praznik_JMIM -k 0.5 -c SKCM > logfiles/mlrFilter_praznik_JMIM_0.5_zscaled_SKCM.log >&1&
# Rscript --vanilla cancer_zscale_oriSig_merge_prediction_mlr_finalmodel.R -s mlrFilter -f praznik_JMIM -k 0.7 -c SKCM > logfiles/mlrFilter_praznik_JMIM_0.7_zscaled_SKCM.log >&1&
# Rscript --vanilla cancer_zscale_oriSig_merge_prediction_mlr_finalmodel.R -s mlrFilter -f praznik_JMIM -k 0.9 -c SKCM > logfiles/mlrFilter_praznik_JMIM_0.9_zscaled_SKCM.log >&1&
wait

Rscript --vanilla cancer_zscale_oriSig_merge_prediction_mlr_finalmodel.R -s mlrFilter -f praznik_MRMR -k 0.1 -c SKCM > logfiles/mlrFilter_praznik_MRMR_0.1_zscaled_SKCM.log >&1&
Rscript --vanilla cancer_zscale_oriSig_merge_prediction_mlr_finalmodel.R -s mlrFilter -f praznik_MRMR -k 0.3 -c SKCM > logfiles/mlrFilter_praznik_MRMR_0.3_zscaled_SKCM.log >&1&
Rscript --vanilla cancer_zscale_oriSig_merge_prediction_mlr_finalmodel.R -s mlrFilter -f praznik_MRMR -k 0.5 -c SKCM > logfiles/mlrFilter_praznik_MRMR_0.5_zscaled_SKCM.log >&1&
# Rscript --vanilla cancer_zscale_oriSig_merge_prediction_mlr_finalmodel.R -s mlrFilter -f praznik_MRMR -k 0.7 -c SKCM > logfiles/mlrFilter_praznik_MRMR_0.7_zscaled_SKCM.log >&1&
# Rscript --vanilla cancer_zscale_oriSig_merge_prediction_mlr_finalmodel.R -s mlrFilter -f praznik_MRMR -k 0.9 -c SKCM > logfiles/mlrFilter_praznik_MRMR_0.9_zscaled_SKCM.log >&1&

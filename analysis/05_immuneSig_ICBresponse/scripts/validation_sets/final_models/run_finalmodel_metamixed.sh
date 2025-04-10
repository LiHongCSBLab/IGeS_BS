#! /bin/bash
Rscript --vanilla dataset_oriSig_merge_prediction_mlr_finalmodel.R -s meta_mixed -f ranger_impurity -k 0.1 > logfiles/meta_mixed_ranger_impurity_0.1.log 2>&1&
Rscript --vanilla dataset_oriSig_merge_prediction_mlr_finalmodel.R -s meta_mixed -f ranger_impurity -k 0.2 > logfiles/meta_mixed_ranger_impurity_0.2.log 2>&1&
Rscript --vanilla dataset_oriSig_merge_prediction_mlr_finalmodel.R -s meta_mixed -f ranger_impurity -k 0.3 > logfiles/meta_mixed_ranger_impurity_0.3.log 2>&1&
Rscript --vanilla dataset_oriSig_merge_prediction_mlr_finalmodel.R -s meta_mixed -f ranger_impurity -k 0.5 > logfiles/meta_mixed_ranger_impurity_0.5.log 2>&1&
# Rscript --vanilla dataset_oriSig_merge_prediction_mlr_finalmodel.R -s meta_mixed -f ranger_impurity -k 0.7 > logfiles/meta_mixed_ranger_impurity_0.7.log 2>&1&
# Rscript --vanilla dataset_oriSig_merge_prediction_mlr_finalmodel.R -s meta_mixed -f ranger_impurity -k 0.9 > logfiles/meta_mixed_ranger_impurity_0.9.log 2>&1&
wait

Rscript --vanilla dataset_zscale_oriSig_merge_prediction_mlr_finalmodel.R -s meta_mixed -f ranger_impurity -k 0.1 > logfiles/meta_mixed_ranger_impurity_0.1_zscaled.log 2>&1&
Rscript --vanilla dataset_zscale_oriSig_merge_prediction_mlr_finalmodel.R -s meta_mixed -f ranger_impurity -k 0.2 > logfiles/meta_mixed_ranger_impurity_0.2_zscaled.log 2>&1&
Rscript --vanilla dataset_zscale_oriSig_merge_prediction_mlr_finalmodel.R -s meta_mixed -f ranger_impurity -k 0.3 > logfiles/meta_mixed_ranger_impurity_0.3_zscaled.log 2>&1&
Rscript --vanilla dataset_zscale_oriSig_merge_prediction_mlr_finalmodel.R -s meta_mixed -f ranger_impurity -k 0.5 > logfiles/meta_mixed_ranger_impurity_0.5_zscaled.log 2>&1&
# Rscript --vanilla dataset_zscale_oriSig_merge_prediction_mlr_finalmodel.R -s meta_mixed -f ranger_impurity -k 0.7 > logfiles/meta_mixed_ranger_impurity_0.7_zscaled.log 2>&1&
# Rscript --vanilla dataset_zscale_oriSig_merge_prediction_mlr_finalmodel.R -s meta_mixed -f ranger_impurity -k 0.9 > logfiles/meta_mixed_ranger_impurity_0.9_zscaled.log 2>&1&
wait
Rscript --vanilla cancer_oriSig_merge_prediction_mlr_finalmodel.R -s meta_mixed -f ranger_impurity -k 0.1 -c SKCM > logfiles/meta_mixed_ranger_impurity_0.1_SKCM.log >&1&
Rscript --vanilla cancer_oriSig_merge_prediction_mlr_finalmodel.R -s meta_mixed -f ranger_impurity -k 0.2 -c SKCM > logfiles/meta_mixed_ranger_impurity_0.2_SKCM.log >&1&
Rscript --vanilla cancer_oriSig_merge_prediction_mlr_finalmodel.R -s meta_mixed -f ranger_impurity -k 0.3 -c SKCM > logfiles/meta_mixed_ranger_impurity_0.3_SKCM.log >&1&
Rscript --vanilla cancer_oriSig_merge_prediction_mlr_finalmodel.R -s meta_mixed -f ranger_impurity -k 0.5 -c SKCM > logfiles/meta_mixed_ranger_impurity_0.5_SKCM.log >&1&
# Rscript --vanilla cancer_oriSig_merge_prediction_mlr_finalmodel.R -s meta_mixed -f ranger_impurity -k 0.7 -c SKCM > logfiles/meta_mixed_ranger_impurity_0.7_SKCM.log >&1&
# Rscript --vanilla cancer_oriSig_merge_prediction_mlr_finalmodel.R -s meta_mixed -f ranger_impurity -k 0.9 -c SKCM > logfiles/meta_mixed_ranger_impurity_0.9_SKCM.log >&1&
wait


Rscript --vanilla cancer_zscale_oriSig_merge_prediction_mlr_finalmodel.R -s meta_mixed -f ranger_impurity -k 0.1 -c SKCM > logfiles/meta_mixed_ranger_impurity_0.1_zscaled_SKCM.log >&1&
Rscript --vanilla cancer_zscale_oriSig_merge_prediction_mlr_finalmodel.R -s meta_mixed -f ranger_impurity -k 0.2 -c SKCM > logfiles/meta_mixed_ranger_impurity_0.2_zscaled_SKCM.log >&1&
Rscript --vanilla cancer_zscale_oriSig_merge_prediction_mlr_finalmodel.R -s meta_mixed -f ranger_impurity -k 0.3 -c SKCM > logfiles/meta_mixed_ranger_impurity_0.3_zscaled_SKCM.log >&1&
Rscript --vanilla cancer_zscale_oriSig_merge_prediction_mlr_finalmodel.R -s meta_mixed -f ranger_impurity -k 0.5 -c SKCM > logfiles/meta_mixed_ranger_impurity_0.5_zscaled_SKCM.log >&1&
# Rscript --vanilla cancer_zscale_oriSig_merge_prediction_mlr_finalmodel.R -s meta_mixed -f ranger_impurity -k 0.7 -c SKCM > logfiles/meta_mixed_ranger_impurity_0.7_zscaled_SKCM.log >&1&
# Rscript --vanilla cancer_zscale_oriSig_merge_prediction_mlr_finalmodel.R -s meta_mixed -f ranger_impurity -k 0.9 -c SKCM > logfiles/meta_mixed_ranger_impurity_0.9_zscaled_SKCM.log >&1&
wait

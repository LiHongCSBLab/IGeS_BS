#! /bash/bin


Rscript --vanilla dataset_oriSig_merge_prediction_mlr_finalmodel.R -s mlrFilter -f ranger_impurity -k 0.1 > logfiles/mlrFilter_ranger_impurity_0.1.log 2>&1&
Rscript --vanilla dataset_oriSig_merge_prediction_mlr_finalmodel.R -s mlrFilter -f ranger_impurity -k 0.2 > logfiles/mlrFilter_ranger_impurity_0.2.log 2>&1&
Rscript --vanilla dataset_oriSig_merge_prediction_mlr_finalmodel.R -s mlrFilter -f ranger_impurity -k 0.3 > logfiles/mlrFilter_ranger_impurity_0.3.log 2>&1&
Rscript --vanilla dataset_oriSig_merge_prediction_mlr_finalmodel.R -s mlrFilter -f ranger_impurity -k 0.5 > logfiles/mlrFilter_ranger_impurity_0.5.log 2>&1&
# Rscript --vanilla dataset_oriSig_merge_prediction_mlr_finalmodel.R -s mlrFilter -f ranger_impurity -k 0.7 > logfiles/mlrFilter_ranger_impurity_0.7.log 2>&1&
# Rscript --vanilla dataset_oriSig_merge_prediction_mlr_finalmodel.R -s mlrFilter -f ranger_impurity -k 0.9 > logfiles/mlrFilter_ranger_impurity_0.9.log 2>&1&

Rscript --vanilla dataset_oriSig_merge_prediction_mlr_finalmodel.R -s mixed -f ranger_impurity -k 0.1 > logfiles/mixed_ranger_impurity_0.1.log 2>&1&
Rscript --vanilla dataset_oriSig_merge_prediction_mlr_finalmodel.R -s mixed -f ranger_impurity -k 0.2 > logfiles/mixed_ranger_impurity_0.2.log 2>&1&
Rscript --vanilla dataset_oriSig_merge_prediction_mlr_finalmodel.R -s mixed -f ranger_impurity -k 0.3 > logfiles/mixed_ranger_impurity_0.3.log 2>&1&
Rscript --vanilla dataset_oriSig_merge_prediction_mlr_finalmodel.R -s mixed -f ranger_impurity -k 0.5 > logfiles/mixed_ranger_impurity_0.5.log 2>&1&
# Rscript --vanilla dataset_oriSig_merge_prediction_mlr_finalmodel.R -s mixed -f ranger_impurity -k 0.7 > logfiles/mixed_ranger_impurity_0.7.log 2>&1&
# Rscript --vanilla dataset_oriSig_merge_prediction_mlr_finalmodel.R -s mixed -f ranger_impurity -k 0.9 > logfiles/mixed_ranger_impurity_0.9.log 2>&1&
wait

Rscript --vanilla dataset_oriSig_merge_prediction_mlr_finalmodel.R -s metaAnalysis > logfiles/metaAnalysis.log 2>&1&
Rscript --vanilla dataset_oriSig_merge_prediction_mlr_finalmodel.R -s combinedP > logfiles/combinedP.log 2>&1&
Rscript --vanilla dataset_oriSig_merge_prediction_mlr_finalmodel.R -s voting > logfiles/voting.log 2>&1&

Rscript --vanilla dataset_oriSig_merge_prediction_mlr_finalmodel.R -s mlrFilter -f FSelector_symmetrical.uncertainty > logfiles/mlrFilter_FSelector_symmetrical.uncertainty.log 2>&1&
wait

Rscript --vanilla dataset_oriSig_merge_prediction_mlr_finalmodel.R -s mlrFilter -f auc -k 0.1 > logfiles/mlrFilter_auc_0.1.log 2>&1&
Rscript --vanilla dataset_oriSig_merge_prediction_mlr_finalmodel.R -s mlrFilter -f auc -k 0.3 > logfiles/mlrFilter_auc_0.3.log 2>&1&
Rscript --vanilla dataset_oriSig_merge_prediction_mlr_finalmodel.R -s mlrFilter -f auc -k 0.5 > logfiles/mlrFilter_auc_0.5.log 2>&1&
# Rscript --vanilla dataset_oriSig_merge_prediction_mlr_finalmodel.R -s mlrFilter -f auc -k 0.7 > logfiles/mlrFilter_auc_0.7.log 2>&1&
# Rscript --vanilla dataset_oriSig_merge_prediction_mlr_finalmodel.R -s mlrFilter -f auc -k 0.9 > logfiles/mlrFilter_auc_0.9.log 2>&1&
wait

Rscript --vanilla dataset_oriSig_merge_prediction_mlr_finalmodel.R -s mlrFilter -f praznik_JMIM -k 0.1 > logfiles/mlrFilter_praznik_JMIM_0.1.log 2>&1&
Rscript --vanilla dataset_oriSig_merge_prediction_mlr_finalmodel.R -s mlrFilter -f praznik_JMIM -k 0.3 > logfiles/mlrFilter_praznik_JMIM_0.3.log 2>&1&
Rscript --vanilla dataset_oriSig_merge_prediction_mlr_finalmodel.R -s mlrFilter -f praznik_JMIM -k 0.5 > logfiles/mlrFilter_praznik_JMIM_0.5.log 2>&1&
# Rscript --vanilla dataset_oriSig_merge_prediction_mlr_finalmodel.R -s mlrFilter -f praznik_JMIM -k 0.7 > logfiles/mlrFilter_praznik_JMIM_0.7.log 2>&1&
# Rscript --vanilla dataset_oriSig_merge_prediction_mlr_finalmodel.R -s mlrFilter -f praznik_JMIM -k 0.9 > logfiles/mlrFilter_praznik_JMIM_0.9.log 2>&1&
wait

Rscript --vanilla dataset_oriSig_merge_prediction_mlr_finalmodel.R -s mlrFilter -f praznik_MRMR -k 0.1 > logfiles/mlrFilter_praznik_MRMR_0.1.log 2>&1&
Rscript --vanilla dataset_oriSig_merge_prediction_mlr_finalmodel.R -s mlrFilter -f praznik_MRMR -k 0.3 > logfiles/mlrFilter_praznik_MRMR_0.3.log 2>&1&
Rscript --vanilla dataset_oriSig_merge_prediction_mlr_finalmodel.R -s mlrFilter -f praznik_MRMR -k 0.5 > logfiles/mlrFilter_praznik_MRMR_0.5.log 2>&1&
# Rscript --vanilla dataset_oriSig_merge_prediction_mlr_finalmodel.R -s mlrFilter -f praznik_MRMR -k 0.7 > logfiles/mlrFilter_praznik_MRMR_0.7.log 2>&1&
# Rscript --vanilla dataset_oriSig_merge_prediction_mlr_finalmodel.R -s mlrFilter -f praznik_MRMR -k 0.9 > logfiles/mlrFilter_praznik_MRMR_0.9.log 2>&1&

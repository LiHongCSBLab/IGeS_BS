#! /bash/bin


Rscript --vanilla dataset_zscaled_oriSig_validationset.R -s mlrFilter -f ranger_impurity -k 0.1 -d Liu_dataset_CTLA4Naive_outcome > logfiles/mlrFilter_ranger_impurity_0.1_Liu_zscaled.log 2>&1&
Rscript --vanilla dataset_zscaled_oriSig_validationset.R -s mlrFilter -f ranger_impurity -k 0.2 -d Liu_dataset_CTLA4Naive_outcome > logfiles/mlrFilter_ranger_impurity_0.2_Liu_zscaled.log 2>&1&
Rscript --vanilla dataset_zscaled_oriSig_validationset.R -s mlrFilter -f ranger_impurity -k 0.3 -d Liu_dataset_CTLA4Naive_outcome > logfiles/mlrFilter_ranger_impurity_0.3_Liu_zscaled.log 2>&1&
Rscript --vanilla dataset_zscaled_oriSig_validationset.R -s mlrFilter -f ranger_impurity -k 0.5 -d Liu_dataset_CTLA4Naive_outcome > logfiles/mlrFilter_ranger_impurity_0.5_Liu_zscaled.log 2>&1&
Rscript --vanilla dataset_zscaled_oriSig_validationset.R -s mlrFilter -f ranger_impurity -k 0.7 -d Liu_dataset_CTLA4Naive_outcome > logfiles/mlrFilter_ranger_impurity_0.7_Liu_zscaled.log 2>&1&
Rscript --vanilla dataset_zscaled_oriSig_validationset.R -s mlrFilter -f ranger_impurity -k 0.9 -d Liu_dataset_CTLA4Naive_outcome > logfiles/mlrFilter_ranger_impurity_0.9_Liu_zscaled.log 2>&1&

Rscript --vanilla dataset_zscaled_oriSig_validationset.R -s mixed -f ranger_impurity -k 0.1 -d Liu_dataset_CTLA4Naive_outcome > logfiles/mixed_ranger_impurity_0.1_Liu_zscaled.log 2>&1&
Rscript --vanilla dataset_zscaled_oriSig_validationset.R -s mixed -f ranger_impurity -k 0.2 -d Liu_dataset_CTLA4Naive_outcome > logfiles/mixed_ranger_impurity_0.2_Liu_zscaled.log 2>&1&
Rscript --vanilla dataset_zscaled_oriSig_validationset.R -s mixed -f ranger_impurity -k 0.3 -d Liu_dataset_CTLA4Naive_outcome > logfiles/mixed_ranger_impurity_0.3_Liu_zscaled.log 2>&1&
Rscript --vanilla dataset_zscaled_oriSig_validationset.R -s mixed -f ranger_impurity -k 0.5 -d Liu_dataset_CTLA4Naive_outcome > logfiles/mixed_ranger_impurity_0.5_Liu_zscaled.log 2>&1&
Rscript --vanilla dataset_zscaled_oriSig_validationset.R -s mixed -f ranger_impurity -k 0.7 -d Liu_dataset_CTLA4Naive_outcome > logfiles/mixed_ranger_impurity_0.7_Liu_zscaled.log 2>&1&
Rscript --vanilla dataset_zscaled_oriSig_validationset.R -s mixed -f ranger_impurity -k 0.9 -d Liu_dataset_CTLA4Naive_outcome > logfiles/mixed_ranger_impurity_0.9_Liu_zscaled.log 2>&1&
wait

Rscript --vanilla dataset_zscaled_oriSig_validationset.R -s metaAnalysis -d Liu_dataset_CTLA4Naive_outcome > logfiles/metaAnalysis_Liu_zscaled.log 2>&1&
Rscript --vanilla dataset_zscaled_oriSig_validationset.R -s combinedP -d Liu_dataset_CTLA4Naive_outcome > logfiles/combinedP_Liu_zscaled.log 2>&1&
Rscript --vanilla dataset_zscaled_oriSig_validationset.R -s voting -d Liu_dataset_CTLA4Naive_outcome > logfiles/voting_Liu_zscaled.log 2>&1&

Rscript --vanilla dataset_zscaled_oriSig_validationset.R -s mlrFilter -f FSelector_symmetrical.uncertainty -d Liu_dataset_CTLA4Naive_outcome > logfiles/mlrFilter_FSelector_symmetrical.uncertainty_Liu_zscaled.log 2>&1&
wait

Rscript --vanilla dataset_zscaled_oriSig_validationset.R -s mlrFilter -f auc -k 0.1 -d Liu_dataset_CTLA4Naive_outcome > logfiles/mlrFilter_auc_0.1_Liu_zscaled.log 2>&1&
Rscript --vanilla dataset_zscaled_oriSig_validationset.R -s mlrFilter -f auc -k 0.3 -d Liu_dataset_CTLA4Naive_outcome > logfiles/mlrFilter_auc_0.3_Liu_zscaled.log 2>&1&
Rscript --vanilla dataset_zscaled_oriSig_validationset.R -s mlrFilter -f auc -k 0.5 -d Liu_dataset_CTLA4Naive_outcome > logfiles/mlrFilter_auc_0.5_Liu_zscaled.log 2>&1&
Rscript --vanilla dataset_zscaled_oriSig_validationset.R -s mlrFilter -f auc -k 0.7 -d Liu_dataset_CTLA4Naive_outcome > logfiles/mlrFilter_auc_0.7_Liu_zscaled.log 2>&1&
Rscript --vanilla dataset_zscaled_oriSig_validationset.R -s mlrFilter -f auc -k 0.9 -d Liu_dataset_CTLA4Naive_outcome > logfiles/mlrFilter_auc_0.9_Liu_zscaled.log 2>&1&
wait

Rscript --vanilla dataset_zscaled_oriSig_validationset.R -s mlrFilter -f praznik_JMIM -k 0.1 -d Liu_dataset_CTLA4Naive_outcome > logfiles/mlrFilter_praznik_JMIM_0.1_Liu_zscaled.log 2>&1&
Rscript --vanilla dataset_zscaled_oriSig_validationset.R -s mlrFilter -f praznik_JMIM -k 0.3 -d Liu_dataset_CTLA4Naive_outcome > logfiles/mlrFilter_praznik_JMIM_0.3_Liu_zscaled.log 2>&1&
Rscript --vanilla dataset_zscaled_oriSig_validationset.R -s mlrFilter -f praznik_JMIM -k 0.5 -d Liu_dataset_CTLA4Naive_outcome > logfiles/mlrFilter_praznik_JMIM_0.5_Liu_zscaled.log 2>&1&
Rscript --vanilla dataset_zscaled_oriSig_validationset.R -s mlrFilter -f praznik_JMIM -k 0.7 -d Liu_dataset_CTLA4Naive_outcome > logfiles/mlrFilter_praznik_JMIM_0.7_Liu_zscaled.log 2>&1&
Rscript --vanilla dataset_zscaled_oriSig_validationset.R -s mlrFilter -f praznik_JMIM -k 0.9 -d Liu_dataset_CTLA4Naive_outcome > logfiles/mlrFilter_praznik_JMIM_0.9_Liu_zscaled.log 2>&1&
wait

Rscript --vanilla dataset_zscaled_oriSig_validationset.R -s mlrFilter -f praznik_MRMR -k 0.1 -d Liu_dataset_CTLA4Naive_outcome > logfiles/mlrFilter_praznik_MRMR_0.1_Liu_zscaled.log 2>&1&
Rscript --vanilla dataset_zscaled_oriSig_validationset.R -s mlrFilter -f praznik_MRMR -k 0.3 -d Liu_dataset_CTLA4Naive_outcome > logfiles/mlrFilter_praznik_MRMR_0.3_Liu_zscaled.log 2>&1&
Rscript --vanilla dataset_zscaled_oriSig_validationset.R -s mlrFilter -f praznik_MRMR -k 0.5 -d Liu_dataset_CTLA4Naive_outcome > logfiles/mlrFilter_praznik_MRMR_0.5_Liu_zscaled.log 2>&1&
Rscript --vanilla dataset_zscaled_oriSig_validationset.R -s mlrFilter -f praznik_MRMR -k 0.7 -d Liu_dataset_CTLA4Naive_outcome > logfiles/mlrFilter_praznik_MRMR_0.7_Liu_zscaled.log 2>&1&
Rscript --vanilla dataset_zscaled_oriSig_validationset.R -s mlrFilter -f praznik_MRMR -k 0.9 -d Liu_dataset_CTLA4Naive_outcome > logfiles/mlrFilter_praznik_MRMR_0.9_Liu_zscaled.log 2>&1&

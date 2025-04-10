#! /bash/bin

Rscript --vanilla cancer_zscale_oriSig_validationset.R -s mlrFilter -f ranger_impurity -k 0.1 -d 08_Lauss_dataset_outcome -c SKCM > logfiles/mlrFilter_ranger_impurity_0.1_Lauss_zscaled_SKCM.log >&1&
Rscript --vanilla cancer_zscale_oriSig_validationset.R -s mlrFilter -f ranger_impurity -k 0.3 -d 08_Lauss_dataset_outcome -c SKCM > logfiles/mlrFilter_ranger_impurity_0.3_Lauss_zscaled_SKCM.log >&1&
Rscript --vanilla cancer_zscale_oriSig_validationset.R -s mlrFilter -f ranger_impurity -k 0.5 -d 08_Lauss_dataset_outcome -c SKCM > logfiles/mlrFilter_ranger_impurity_0.5_Lauss_zscaled_SKCM.log >&1&
Rscript --vanilla cancer_zscale_oriSig_validationset.R -s mlrFilter -f ranger_impurity -k 0.7 -d 08_Lauss_dataset_outcome -c SKCM > logfiles/mlrFilter_ranger_impurity_0.7_Lauss_zscaled_SKCM.log >&1&
Rscript --vanilla cancer_zscale_oriSig_validationset.R -s mlrFilter -f ranger_impurity -k 0.9 -d 08_Lauss_dataset_outcome -c SKCM > logfiles/mlrFilter_ranger_impurity_0.9_Lauss_zscaled_SKCM.log >&1&


Rscript --vanilla cancer_zscale_oriSig_validationset.R -s mixed -f ranger_impurity -k 0.1 -d 08_Lauss_dataset_outcome -c SKCM > logfiles/mixed_ranger_impurity_0.1_Lauss_zscaled_SKCM.log >&1&
Rscript --vanilla cancer_zscale_oriSig_validationset.R -s mixed -f ranger_impurity -k 0.3 -d 08_Lauss_dataset_outcome -c SKCM > logfiles/mixed_ranger_impurity_0.3_Lauss_zscaled_SKCM.log >&1&
Rscript --vanilla cancer_zscale_oriSig_validationset.R -s mixed -f ranger_impurity -k 0.5 -d 08_Lauss_dataset_outcome -c SKCM > logfiles/mixed_ranger_impurity_0.5_Lauss_zscaled_SKCM.log >&1&
Rscript --vanilla cancer_zscale_oriSig_validationset.R -s mixed -f ranger_impurity -k 0.7 -d 08_Lauss_dataset_outcome -c SKCM > logfiles/mixed_ranger_impurity_0.7_Lauss_zscaled_SKCM.log >&1&
Rscript --vanilla cancer_zscale_oriSig_validationset.R -s mixed -f ranger_impurity -k 0.9 -d 08_Lauss_dataset_outcome -c SKCM > logfiles/mixed_ranger_impurity_0.9_Lauss_zscaled_SKCM.log >&1&

wait


Rscript --vanilla cancer_zscale_oriSig_validationset.R -s metaAnalysis -d 08_Lauss_dataset_outcome -c SKCM > logfiles/metaAnalysis_Lauss_zscaled_SKCM.log >&1&
Rscript --vanilla cancer_zscale_oriSig_validationset.R -s combinedP -d 08_Lauss_dataset_outcome -c SKCM > logfiles/combinedP_Lauss_zscaled_SKCM.log >&1&
Rscript --vanilla cancer_zscale_oriSig_validationset.R -s voting -d 08_Lauss_dataset_outcome -c SKCM > logfiles/voting_Lauss_zscaled_SKCM.log >&1&

Rscript --vanilla cancer_zscale_oriSig_validationset.R -s mlrFilter -f FSelector_symmetrical.uncertainty -d 08_Lauss_dataset_outcome -c SKCM > logfiles/mlrFilter_FSelector_symmetrical.uncertainty_Lauss_zscaled_SKCM.log >&1&
wait

Rscript --vanilla cancer_zscale_oriSig_validationset.R -s mlrFilter -f auc -k 0.1 -d 08_Lauss_dataset_outcome -c SKCM > logfiles/mlrFilter_auc_0.1_Lauss_zscaled_SKCM.log >&1&
Rscript --vanilla cancer_zscale_oriSig_validationset.R -s mlrFilter -f auc -k 0.3 -d 08_Lauss_dataset_outcome -c SKCM > logfiles/mlrFilter_auc_0.3_Lauss_zscaled_SKCM.log >&1&
Rscript --vanilla cancer_zscale_oriSig_validationset.R -s mlrFilter -f auc -k 0.5 -d 08_Lauss_dataset_outcome -c SKCM > logfiles/mlrFilter_auc_0.5_Lauss_zscaled_SKCM.log >&1&
Rscript --vanilla cancer_zscale_oriSig_validationset.R -s mlrFilter -f auc -k 0.7 -d 08_Lauss_dataset_outcome -c SKCM > logfiles/mlrFilter_auc_0.7_Lauss_zscaled_SKCM.log >&1&
Rscript --vanilla cancer_zscale_oriSig_validationset.R -s mlrFilter -f auc -k 0.9 -d 08_Lauss_dataset_outcome -c SKCM > logfiles/mlrFilter_auc_0.9_Lauss_zscaled_SKCM.log >&1&
wait

Rscript --vanilla cancer_zscale_oriSig_validationset.R -s mlrFilter -f praznik_JMIM -k 0.1 -d 08_Lauss_dataset_outcome -c SKCM > logfiles/mlrFilter_praznik_JMIM_0.1_Lauss_zscaled_SKCM.log >&1&
Rscript --vanilla cancer_zscale_oriSig_validationset.R -s mlrFilter -f praznik_JMIM -k 0.3 -d 08_Lauss_dataset_outcome -c SKCM > logfiles/mlrFilter_praznik_JMIM_0.3_Lauss_zscaled_SKCM.log >&1&
Rscript --vanilla cancer_zscale_oriSig_validationset.R -s mlrFilter -f praznik_JMIM -k 0.5 -d 08_Lauss_dataset_outcome -c SKCM > logfiles/mlrFilter_praznik_JMIM_0.5_Lauss_zscaled_SKCM.log >&1&
Rscript --vanilla cancer_zscale_oriSig_validationset.R -s mlrFilter -f praznik_JMIM -k 0.7 -d 08_Lauss_dataset_outcome -c SKCM > logfiles/mlrFilter_praznik_JMIM_0.7_Lauss_zscaled_SKCM.log >&1&
Rscript --vanilla cancer_zscale_oriSig_validationset.R -s mlrFilter -f praznik_JMIM -k 0.9 -d 08_Lauss_dataset_outcome -c SKCM > logfiles/mlrFilter_praznik_JMIM_0.9_Lauss_zscaled_SKCM.log >&1&
wait

Rscript --vanilla cancer_zscale_oriSig_validationset.R -s mlrFilter -f praznik_MRMR -k 0.1 -d 08_Lauss_dataset_outcome -c SKCM > logfiles/mlrFilter_praznik_MRMR_0.1_Lauss_zscaled_SKCM.log >&1&
Rscript --vanilla cancer_zscale_oriSig_validationset.R -s mlrFilter -f praznik_MRMR -k 0.3 -d 08_Lauss_dataset_outcome -c SKCM > logfiles/mlrFilter_praznik_MRMR_0.3_Lauss_zscaled_SKCM.log >&1&
Rscript --vanilla cancer_zscale_oriSig_validationset.R -s mlrFilter -f praznik_MRMR -k 0.5 -d 08_Lauss_dataset_outcome -c SKCM > logfiles/mlrFilter_praznik_MRMR_0.5_Lauss_zscaled_SKCM.log >&1&
Rscript --vanilla cancer_zscale_oriSig_validationset.R -s mlrFilter -f praznik_MRMR -k 0.7 -d 08_Lauss_dataset_outcome -c SKCM > logfiles/mlrFilter_praznik_MRMR_0.7_Lauss_zscaled_SKCM.log >&1&
Rscript --vanilla cancer_zscale_oriSig_validationset.R -s mlrFilter -f praznik_MRMR -k 0.9 -d 08_Lauss_dataset_outcome -c SKCM > logfiles/mlrFilter_praznik_MRMR_0.9_Lauss_zscaled_SKCM.log >&1&

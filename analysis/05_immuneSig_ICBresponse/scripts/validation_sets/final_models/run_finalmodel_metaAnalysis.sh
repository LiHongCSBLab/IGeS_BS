#! /bin/bash
Rscript --vanilla dataset_zscale_oriSig_merge_prediction_mlr_finalmodel.R -s metaAnalysis > logfiles/metaAnalysis_zscaled.log 2>&1&
Rscript --vanilla dataset_oriSig_merge_prediction_mlr_finalmodel.R -s metaAnalysis > logfiles/metaAnalysis_unscaled.log 2>&1&

Rscript --vanilla cancer_oriSig_merge_prediction_mlr_finalmodel.R -s metaAnalysis -c SKCM > logfiles/metaAnalysis_unscaled_SKCM.log >&1&
Rscript --vanilla cancer_zscale_oriSig_merge_prediction_mlr_finalmodel.R -s metaAnalysis -c SKCM > logfiles/metaAnalysis_zscaled_SKCM.log >&1&


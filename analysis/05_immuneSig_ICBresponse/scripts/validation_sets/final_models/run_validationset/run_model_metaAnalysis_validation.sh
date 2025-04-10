#! /bin/bash
Rscript --vanilla dataset_oriSig_validationset.R -s metaAnalysis -d 08_Lauss_dataset_outcome > logfiles/metaAnalysis_Lauss.log 2>&1&
Rscript --vanilla dataset_zscaled_oriSig_validationset.R -s metaAnalysis -d 08_Lauss_dataset_outcome > logfiles/metaAnalysis_Lauss_zscaled.log 2>&1&
Rscript --vanilla cancer_oriSig_validationset.R -s metaAnalysis -d 08_Lauss_dataset_outcome -c SKCM > logfiles/metaAnalysis_Lauss_SKCM.log >&1&
Rscript --vanilla cancer_zscale_oriSig_validationset.R -s metaAnalysis -d 08_Lauss_dataset_outcome -c SKCM > logfiles/metaAnalysis_Lauss_zscaled_SKCM.log >&1&

Rscript --vanilla dataset_oriSig_validationset.R -s metaAnalysis -d Liu_dataset_CTLA4Naive_outcome > logfiles/metaAnalysis_Liu.log 2>&1&
Rscript --vanilla dataset_zscaled_oriSig_validationset.R -s metaAnalysis -d Liu_dataset_CTLA4Naive_outcome > logfiles/metaAnalysis_Liu_zscaled.log 2>&1&
Rscript --vanilla cancer_oriSig_validationset.R -s metaAnalysis -d Liu_dataset_CTLA4Naive_outcome -c SKCM > logfiles/metaAnalysis_Liu_SKCM.log >&1&
Rscript --vanilla cancer_zscale_oriSig_validationset.R -s metaAnalysis -d Liu_dataset_CTLA4Naive_outcome -c SKCM > logfiles/metaAnalysis_Liu_zscaled_SKCM.log >&1&

Rscript --vanilla dataset_oriSig_validationset.R -s metaAnalysis -d GSE96619_outcome > logfiles/metaAnalysis_GSE96619.log 2>&1&
Rscript --vanilla dataset_zscaled_oriSig_validationset.R -s metaAnalysis -d GSE96619_outcome > logfiles/metaAnalysis_GSE96619_zscaled.log 2>&1&
Rscript --vanilla cancer_oriSig_validationset.R -s metaAnalysis -d GSE96619_outcome -c SKCM > logfiles/metaAnalysis_GSE96619_SKCM.log >&1&
Rscript --vanilla cancer_zscale_oriSig_validationset.R -s metaAnalysis -d GSE96619_outcome -c SKCM > logfiles/metaAnalysis_GSE96619_zscaled_SKCM.log >&1&

Rscript --vanilla dataset_oriSig_validationset.R -s metaAnalysis -d GSE115821_pre_outcome > logfiles/metaAnalysis_GSE115821.log 2>&1&
Rscript --vanilla dataset_zscaled_oriSig_validationset.R -s metaAnalysis -d GSE115821_pre_outcome > logfiles/metaAnalysis_GSE115821_zscaled.log 2>&1&
Rscript --vanilla cancer_oriSig_validationset.R -s metaAnalysis -d GSE115821_pre_outcome -c SKCM > logfiles/metaAnalysis_GSE115821_SKCM.log >&1&
Rscript --vanilla cancer_zscale_oriSig_validationset.R -s metaAnalysis -d GSE115821_pre_outcome -c SKCM > logfiles/metaAnalysis_GSE115821_zscaled_SKCM.log >&1&



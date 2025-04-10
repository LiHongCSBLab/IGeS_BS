#!/bin/bash
Rscript --vanilla ICB_predictor_evaluation_metaAnalysis.R > meta_analysis_preparation.log 2>&1&
wait
Rscript --vanilla /picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/05_immuneSig_ICBresponse/scripts/meta_analysis/meta_analysis_OriImmuSig.R > meta_for_oriImmuSig.log 2>&1&



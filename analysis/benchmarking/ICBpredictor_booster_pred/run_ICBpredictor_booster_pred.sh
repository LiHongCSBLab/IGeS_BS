#!/bin/bash

source /home/fengfangyoumin/miniconda3/bin/activate r_env

cd /picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/07_plotting_v2/scripts/revision/ICBpredictor_booster_pred
mkdir logfiles

Rscript --vanilla ICBpredictor_booster_pred.R > logfiles/ICBpredictor_booster_pred2.log 2>&1

# Rscript --vanilla ICBpredictor_booster_pred_summary.R > logfiles/ICBpredictor_booster_pred_summary.log 2>&1&


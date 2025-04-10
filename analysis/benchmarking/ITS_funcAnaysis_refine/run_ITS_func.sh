#!/bin/bash

source /home/fengfangyoumin/miniconda3/bin/activate r_env

cd /picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/07_plotting_v2/scripts/revision/ITS_funcAnaysis_refine
mkdir logfiles

Rscript --vanilla fig3_ITS_funcAnaysis_refine.R > logfiles/fig3_ITS_funcAnaysis_refine.log 2>&1


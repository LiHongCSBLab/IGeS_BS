#!/bin/bash

source /home/fengfangyoumin/miniconda3/bin/activate r_env

cd /picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/07_plotting_v2/scripts/revision/CM-Drug-main/
mkdir logfiles


# Rscript --vanilla CM-Drug_part1.R > logfiles/CM-Drug_part1.log 2>&1&
# Rscript --vanilla CM-Drug_part2.R > logfiles/CM-Drug_part2.log 2>&1&
# Rscript --vanilla CM-Drug_LIHC.R > logfiles/CM-Drug_LIHC.log 2>&1
# Rscript --vanilla CM-Drug_LIHC_70138.R > logfiles/CM-Drug_LIHCC_70138.log 2>&1
Rscript --vanilla CM-Drug_LIHC_92742.R > logfiles/CM-Drug_LIHC_92742.log 2>&1


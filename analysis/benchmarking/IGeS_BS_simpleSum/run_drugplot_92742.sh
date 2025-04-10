#!/bin/bash
# run under miniconda environment r_env
source /home/fengfangyoumin/miniconda3/bin/activate r_env


Rscript --vanilla plot_for_drug_metaAna_en_othercancer_92742.R -w model_weight -t unweighted -m simple -p 0.05 -i FALSE -s Type > meta_en_othercancer_92742_unweight_simple.log 2>&1&

#!/bin/bash

source /home/fengfangyoumin/miniconda3/bin/activate r_env

cd /picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/07_plotting_v2/scripts/revision/IGeS_EN_shuffle
mkdir logfiles


echo "Running script1.R... Output will be saved to $LOG1"
Rscript --vanilla ITSssgsea_CombatRef_17datasets.R > logfiles/ITSssgsea_CombatRef_17datasets.log 2>&1

# Ensure script1.R completed successfully
if [ $? -ne 0 ]; then
    echo "Error: ITSssgsea_CombatRef_17datasets.R failed to execute."
    exit 1
fi

echo "ITSssgsea_CombatRef_17datasets.R completed. "


Rscript --vanilla ITSssgsea_CombatRef_17datasets_finalmodel_metaPlusFilter.R > logfiles/ITSssgsea_CombatRef_17datasets_finalmodel_metaPlusFilter.log 2>&1

# Ensure script2.R completed successfully
if [ $? -ne 0 ]; then
    echo "Error: ITSssgsea_CombatRef_17datasets_finalmodel_metaPlusFilter.R failed to execute. "
    exit 1
fi

echo "Both R scripts executed successfully."

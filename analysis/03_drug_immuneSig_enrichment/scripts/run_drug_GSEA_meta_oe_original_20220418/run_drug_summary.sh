

#!/bin/bash
# run under miniconda environment r_env
source /home/fengfangyoumin/miniconda3/bin/activate r_env
# cd /picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/03_drug_immuneSig_enrichment/scripts/

bash 05_OV_92742_drug_GSEA.sh > logfiles/OV_92742.log 2>&1&
bash 05_PAAD_STAD_DLBC_drug_GSEA.sh > logfiles/PAAD_STAD_DLBC.log 2>&1&
bash 05_LUSC_92742_drug_GSEA.sh > logfiles/LUSC_92742.log 2>&1&

bash 05_BRCA_70138_drug_GSEA.sh > logfiles/BRCA_70138.log 2>&1&
bash 05_COAD_70138_drug_GSEA.sh > logfiles/COAD_70138.log 2>&1&
bash 05_LIHC_70138_drug_GSEA.sh > logfiles/LIHC_70138.log 2>&1&
bash 05_LUAD_70138_drug_GSEA.sh > logfiles/LUAD_70138.log 2>&1&
bash 05_PRAD_70138_drug_GSEA.sh > logfiles/PRAD_70138.log 2>&1&
bash 05_READ_70138_drug_GSEA.sh > logfiles/READ_70138.log 2>&1&
bash 05_SKCMMetastatic_70138_drug_GSEA.sh > logfiles/SKCMMetastatic_70138.log 2>&1&
bash 05_SKCMPrimary_70138_drug_GSEA.sh > logfiles/SKCMPrimary_70138.log 2>&1&
wait

bash 05_LIHC_92742_drug_GSEA.sh > logfiles/LIHC_92742.log 2>&1&
bash 05_LUAD_92742_drug_GSEA.sh > logfiles/LUAD_92742.log 2>&1&
bash 05_BRCA_92742_drug_GSEA.sh > logfiles/BRCA_92742.log 2>&1&
bash 05_COAD_92742_drug_GSEA.sh > logfiles/COAD_92742.log 2>&1&
bash 05_PRAD_92742_drug_GSEA.sh > logfiles/PRAD_92742.log 2>&1&
bash 05_READ_92742_drug_GSEA.sh > logfiles/READ_92742.log 2>&1&
bash 05_SKCMMetastatic_92742_drug_GSEA.sh > logfiles/SKCMMetastatic_92742.log 2>&1&
bash 05_SKCMPrimary_92742_drug_GSEA.sh > logfiles/SKCMPrimary_92742.log 2>&1&
bash 05_UCEC_92742_drug_GSEA.sh > logfiles/UCEC_92742.log 2>&1&

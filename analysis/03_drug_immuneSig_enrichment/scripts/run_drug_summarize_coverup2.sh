#!/bin/bash
# run under miniconda environment r_env
source /home/fengfangyoumin/miniconda3/bin/activate r_env


#####======================================================================#####
#####                      pan-cancer: Primary cor                         #####
#####======================================================================#####

# "BRCA" "DLBC" "LAML" "LIHC" "LUAD" "LUSC" "OV"   "PRAD" "SKCM" "STAD" "UCEC" "COAD" "READ" "PAAD"

# cancer=LIHC
# sampleType=Primary
# tumor_purity_method=TUMERIC
# dataset=70138

# run for both positive and negative directory, using drug_summarize_coverup.sh
# nohup bash drug_summarize_coverup.sh PRAD Primary TUMERIC 70138 &

# run for only negative diectory, using: drug_summarize_coverup2.sh
nohup bash drug_summarize_coverup2.sh OV Primary TUMERIC 92742 &

nohup bash drug_summarize_coverup2.sh UCEC Primary TUMERIC 92742 &



# wait

# nohup bash drug_summarize_coverup2.sh BRCA Primary TUMERIC 70138 &
# nohup bash drug_summarize_coverup2.sh BRCA BRCA.Basal TUMERIC 70138 &
# nohup bash drug_summarize_coverup2.sh BRCA BRCA.Normal TUMERIC 70138 &
# nohup bash drug_summarize_coverup2.sh BRCA BRCA.LumB TUMERIC 70138 &
# wait
# nohup bash drug_summarize_coverup2.sh BRCA BRCA.LumA TUMERIC 70138 &
# nohup bash drug_summarize_coverup2.sh BRCA BRCA.Her2 TUMERIC 70138 &
# wait

# nohup bash drug_summarize_coverup2.sh BRCA Primary TUMERIC 92742 &
# nohup bash drug_summarize_coverup2.sh COAD Primary TUMERIC 92742 &
# nohup bash drug_summarize_coverup2.sh LUAD Primary TUMERIC 92742 &
# nohup bash drug_summarize_coverup2.sh SKCM Metastatic TUMERIC 92742 &


# # nohup bash drug_summarize_coverup2.sh SKCM Primary TUMERIC 70138 &


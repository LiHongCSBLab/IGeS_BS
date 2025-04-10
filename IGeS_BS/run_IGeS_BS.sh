#!/bin/bash


path=/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/

# RUN FOR GSE70138 dataset
dataset=70138

# Define the array of cancers (this was missing)
cancers=("LIHC" "BRCA" "COAD" "READ" "LUAD" "PAAD" "PRAD")  # You can add more: ("LIHC" "BRCA" "LUAD")
cancer_type="Primary"

# Loop through the cancers
for cancer in "${cancers[@]}"
do 
    nohup Rscript --vanilla plot_for_drug_IGeS.R \
        -a "$dataset" \
        -c "$cancer" \
        -t "$cancer_type" \
        -u TUMERIC \
        -f TRUE \
        -d "${path}/03_drug_immuneSig_enrichment/" \
        -w "${path}/IGeS_BS/IGeS_BS_tool/" \
        -s "${path}/IGeS_BS/" > "logfiles/IGeS-BS_${cancer}_${cancer_type}_${dataset}.log" 2>&1 &
done


cancer="SKCM"
cancer_type="Metastatic"
nohup Rscript --vanilla plot_for_drug_IGeS.R \
        -a "$dataset" \
        -c "$cancer" \
        -t "$cancer_type" \
        -u TUMERIC \
        -f TRUE \
        -d "${path}/03_drug_immuneSig_enrichment/" \
        -w "${path}/IGeS_BS/IGeS_BS_tool/" \
        -s "${path}/IGeS_BS/" > "logfiles/IGeS-BS_${cancer}_${cancer_type}_${dataset}.log" 2>&1 &

cancer="SKCM"
cancer_type="Primary"
nohup Rscript --vanilla plot_for_drug_IGeS.R \
        -a "$dataset" \
        -c "$cancer" \
        -t "$cancer_type" \
        -u CPE \
        -f TRUE \
        -d "${path}/03_drug_immuneSig_enrichment/" \
        -w "${path}/IGeS_BS/IGeS_BS_tool/" \
        -s "${path}/IGeS_BS/" > "logfiles/IGeS-BS_${cancer}_${cancer_type}_${dataset}.log" 2>&1 &


# RUN FOR GSE92742 dataset
dataset=92742

# Define the array of cancers (this was missing)
cancers=("LIHC" "BRCA" "COAD" "READ" "LUAD" "LUSC" "OV" "PRAD" "STAD" "UCEC")  # You can add more: ("LIHC" "BRCA" "LUAD")
cancer_type="Primary"
# Loop through the cancers
for cancer in "${cancers[@]}"
do 
    nohup Rscript --vanilla plot_for_drug_IGeS.R \
        -a "$dataset" \
        -c "$cancer" \
        -t ${cancer_type} \
        -u TUMERIC \
        -f TRUE \
        -d "${path}/03_drug_immuneSig_enrichment/" \
        -w "${path}/IGeS_BS/IGeS_BS_tool/" \
        -s "${path}/IGeS_BS/" > "logfiles/IGeS-BS_${cancer}_${cancer_type}_${dataset}.log" 2>&1 &

done

cancer="SKCM"
cancer_type="Metastatic"
nohup Rscript --vanilla plot_for_drug_IGeS.R \
        -a "$dataset" \
        -c "$cancer" \
        -t "$cancer_type" \
        -u TUMERIC \
        -f TRUE \
        -d "${path}/03_drug_immuneSig_enrichment/" \
        -w "${path}/IGeS_BS/IGeS_BS_tool/" \
        -s "${path}/IGeS_BS/" > "logfiles/IGeS-BS_${cancer}_${cancer_type}_${dataset}.log" 2>&1 &

cancer="SKCM"
cancer_type="Primary"
nohup Rscript --vanilla plot_for_drug_IGeS.R \
        -a "$dataset" \
        -c "$cancer" \
        -t "$cancer_type" \
        -u CPE \
        -f TRUE \
        -d "${path}/03_drug_immuneSig_enrichment/" \
        -w "${path}/IGeS_BS/IGeS_BS_tool/" \
        -s "${path}/IGeS_BS/"  > "logfiles/IGeS-BS_${cancer}_${cancer_type}_${dataset}.log" 2>&1 &


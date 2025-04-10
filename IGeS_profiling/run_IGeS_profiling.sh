#!/bin/bash
# run under miniconda environment r_env
source conda activate IGeS_BS_env

#####======================================================================#####
#####                      IGeS Profiling config                           #####
#####======================================================================#####


wDir=./IGeS_BS_tool/ # <------------Please change it accordingly
wDir=/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/IGeS_BS/IGeS_BS_tool
cd $wDir


java_path=/path/to/java # <------------Please change it accordingly
java_path=/picb/extprog/bin/linux-x86/java


# output path
output=${wDir}/IGeS_profiling_result # <------------Please change it accordingly
mkdir $output

cancer=LIHC # <------------Please change it accordingly
sampleType=Primary # <------------Please change it accordingly
tumor_purity_method=TUMERIC # <------------Please change it accordingly

dataset=70138 # <------------Please change it accordingly

#####======================================================================#####
#####                    time tracking funciton                            #####
#####======================================================================#####
# Function to print the elapsed time in a user-friendly format
print_time() {
  START_TIME=$1
  END_TIME=$2
  ELAPSED_TIME=$((END_TIME - START_TIME))
  printf "Elapsed time: %02d:%02d:%02d\n" \
         $(($ELAPSED_TIME / 3600)) $(($ELAPSED_TIME % 3600 / 60)) $(($ELAPSED_TIME % 60))
}

# Function to show a simple progress bar
show_progress_bar() {
  local duration=$1
  local interval=1
  local n=$((duration / interval))
  
  for i in $(seq 1 $n); do
    printf "#"
    sleep $interval
  done
  echo ""
}

#####======================================================================#####
#####======================================================================#####
# set script location
software=${wDir}/IGeS_profiling/GSEA/
geneset_path=${wDir}/Data/IGeS/
# make dir for log files
mkdir ${wDir}/logfiles


#####======================================================================#####
echo "Step 01"
echo "Extracting compound expression profiles for cancer type: ${cancer}_${sampleType} in dataset: ${dataset}"

START_TIME=$(date +%s)

Rscript --vanilla ${wDir}/IGeS_profiling/01_drug_expr_preparation_update.R -c ${cancer} -d 70138 -w ${wDir} -r ${output} > ${wDir}/logfiles/drug_expr_prep_${cancer}_${dataset}.log 2>&1 

END_TIME=$(date +%s)
print_time $START_TIME $END_TIME
show_progress_bar $((END_TIME - START_TIME))

wait
echo "Extraction is done, proceeding to next step."

#####======================================================================#####

echo "Step 02"
echo "Starting computing IGeS-GSEA  for each compound in each pertubation condition for ${cancer}_${sampleType} in  dataset GSE${dataset}"


START_TIME=$(date +%s)

bash ${wDir}/IGeS_profiling/GSEA/02_drug_GSEA_update.sh ${cancer} ${sampleType} ${tumor_purity_method} ${dataset} ${wDir} ${output} ${software} ${geneset_path} ${java_path} > ${wDir}/logfiles/drug_GSEA_${cancer}_${dataset}.log 2>&1

END_TIME=$(date +%s)
print_time $START_TIME $END_TIME
show_progress_bar $((END_TIME - START_TIME))
wait

echo "IGeS-GSEA update completed for ${cancer}_${sampleType}. Proceeding to summarize results."

#####======================================================================#####

echo "Step 03"
echo "Starting to summarize compound-IGeS GSEA results from different pertubation condition to cancer level"
echo "- ${cancer}_${sampleType}in  dataset GSE${dataset}"

START_TIME=$(date +%s)

Rscript --vanilla ${wDir}/IGeS_profiling/03_summarize_drug_immuneSig_GSEA.R -a ${dataset} -c ${cancer} -t ${sampleType} -u ${tumor_purity_method} -w ${wDir} -s ${output} > ${wDir}/logfiles/drug_summarize_${cancer}_${dataset}.log 2>&1

END_TIME=$(date +%s)
print_time $START_TIME $END_TIME
show_progress_bar $((END_TIME - START_TIME))

wait
echo "IGeS-Profiling is done. Please check output in:"
echo "${output}"

#####======================================================================#####
#####======================================================================#####
#####======================================================================#####
# clean folder after summarize compound-IGeS GSEA results 
# rm -r ${output}/drug_immunesig_IHC_GSEA_result/70138/${cancer}
# rm -r ${output}/data/drug_treated_expression_matrix/70138/${cancer}
# wait



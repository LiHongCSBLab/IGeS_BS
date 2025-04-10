#!/bin/bash

#$ -S /bin/bash
#$ -N gsea
#$ -pe make 1
cancer=$1
sampleType=$2
tumor_purity_method=$3
dataset=$4
wDir=$5
expPath=$6
software=$7
geneset_path=$8
java_path=$9

sign=positive
num_gene=200
method=spearman

# test example:
# cancer=SKCM
# sampleType=Metastatic
# method=spearman
# tumor_purity_method=TUMERIC
# sign=positive
# num_gene=200
# dataset=70138

# wDir=/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/03_drug_immuneSig_enrichment
# software=/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/03_drug_immuneSig_enrichment/scripts/GSEA
# geneset_path=/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/IGeS_BS/IGeS_BS_tool/Data/IGeS/

mkdir ${wDir}/results_${tumor_purity_method}_meta_oe_original/
mkdir ${wDir}/results_${tumor_purity_method}_meta_oe_original/drug_immunesig_${tumor_purity_method}_GSEA_result
mkdir ${wDir}/results_${tumor_purity_method}_meta_oe_original/drug_immunesig_${tumor_purity_method}_GSEA_result/${dataset}
outPath=results_${tumor_purity_method}_meta_oe_original/drug_immunesig_${tumor_purity_method}_GSEA_result/${dataset}/${cancer}_${sampleType}


mkdir ${wDir}/results_${tumor_purity_method}_meta_oe_original/drug_immunesig_${tumor_purity_method}_GSEA_result_files
mkdir ${wDir}/results_${tumor_purity_method}_meta_oe_original/drug_immunesig_${tumor_purity_method}_GSEA_result_files/${dataset}
mkdir ${wDir}/results_${tumor_purity_method}_meta_oe_original/drug_immunesig_${tumor_purity_method}_GSEA_result_files/${dataset}/${cancer}_${sampleType}
desDir=${wDir}/results_${tumor_purity_method}_meta_oe_original/drug_immunesig_${tumor_purity_method}_GSEA_result_files/${dataset}/${cancer}_${sampleType}




if [ ${tumor_purity_method} == CPE ]; then
  geneset=${geneset_path}/${method}_${cancer}_${sampleType}_${sign}_${num_gene}.gmt
else
  geneset=${geneset_path}/${method}_${cancer}_${sampleType}_${sign}_${num_gene}.gmt
fi 

exp=${expPath}/data/drug_treated_expression_matrix/${dataset}/${cancer}/
grp=${expPath}/data/drug_treated_design/for_GSEA/${dataset}/${cancer}/

pairs=treated_versus_DMSO


if [ ! -d ${wDir}/${exp} ]; then
    mkdir ${wDir}/${exp}
fi

if [ ! -d ${wDir}/${grp} ]; then
    mkdir ${wDir}/${grp}
fi

if [ ! -d ${wDir}/${outPath} ]; then
    mkdir ${wDir}/${outPath}
fi
if [ ! -d ${wDir}/${outPath}/${sign}_${num_gene} ]; then
    mkdir ${wDir}/${outPath}/${sign}_${num_gene}
fi
if [ ! -d ${wDir}/${outPath}/${sign}_${num_gene} ]; then
    mkdir ${wDir}/${outPath}/${sign}_${num_gene}
fi

if [ ! -d ${desDir} ]; then
  mkdir ${desDir}
fi

if [ ! -d ${desDir}/${sign}_${num_gene} ]; then
  mkdir ${desDir}/${sign}_${num_gene}
fi
# desDir=${desDir}/${sign}_${num_gene}


i=0
for drug in `ls ${wDir}/${exp}`
do
  # filelist[$i]="$file" ${tmp%.*}
  outDir=${wDir}/${outPath}/${sign}_${num_gene}/${drug}
  if [ ! -d ${outDir} ]; then
        mkdir ${outDir}
    fi
  
  exp_path=${wDir}/${exp}/${drug}
  grp_path=${wDir}/${grp}/${drug}
  
  j=0
  for file in `ls ${exp_path}`
  # for ((j=0;j<=${#exp_path[@]};j++))
  do
    res=${exp_path}/${file}
    cls="${grp_path}/${file%.*}.cls#${pairs}"
    # gmt=${wDir}/${geneset}
    gmt=${geneset}
    
    out=${outDir}/${file%.*}
    if [ ! -d ${out} ]; then
        mkdir ${out}
    fi

    # nohup /picb/extprog/bin/linux-x86/java \
    nohup ${java_path}/java \
    -cp ${software}/gsea-3.0.jar \
    -Xmx2g xtools.gsea.Gsea  \
    -gmx "${gmt}" \
    -res "${res}" \
    -cls "${cls}" \
    -collapse false \
    -mode Max_probe \
    -norm meandiv \
    -nperm 1000  \
    -scoring_scheme weighted \
    -permute gene_set \
    -rpt_label ${pairs} \
    -out "${out}" \
    -include_only_symbols true \
    -make_sets true \
    -plot_top_x 0 \
    -rnd_seed timestamp \
    -set_max 600 \
    -set_min 5 \
    -zip_report false \
    -gui false & # > drug${i}.log 2>& 1 &
    # wait
    ((j++))
  done
  wait


  # copy files to new directory
  desSaveDir=${desDir}/${sign}_${num_gene}/${drug}
  if [ ! -d ${desSaveDir} ]; then
    mkdir ${desSaveDir}
  fi

  for file in `ls ${outDir}`
  do
    cp ${outDir}/${file}/treated_versus_DMSO.Gsea.*/gsea_report_for_treated_*.xls ${desSaveDir}/gsea_report_for_treated_${file}.xls
    cp ${outDir}/${file}/treated_versus_DMSO.Gsea.*/gsea_report_for_DMSO_*.xls ${desSaveDir}/gsea_report_for_DMSO_${file}.xls
  done

  wait
  ((i++))
done


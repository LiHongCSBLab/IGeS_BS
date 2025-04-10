#!/bin/bash

#$ -S /bin/bash
#$ -N gsea
#$ -pe make 1
cancer=$1
sampleType=$2
tumor_purity_method=$3
method=$4
sign=$5
num_gene=$6

dataset=$7
cell=$8

# test example:
cancer=SKCM
sampleType=Metastatic
method=pearson
tumor_purity_method=IHC
sign=all
num_gene=400
dataset=70138
cell=A375

wDir=/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/03_drug_immuneSig_enrichment

software=/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/03_drug_immuneSig_enrichment/scripts/GSEA


pairs=treated_versus_DMSO
geneset=data/gene_immune_sig_${tumor_purity_method}/for_GSEA/pcor/${method}_${cancer}_${sampleType}_${sign}_${num_gene}.gmt

exp=data/drug_treated_expression_matrix/${cancer}/${dataset}_${cell}
grp=data/drug_treated_design/for_GSEA/${cancer}/${dataset}_${cell}

mkdir ${wDir}/results/drug_immunesig_${tumor_purity_method}_GSEA_result/
mkdir ${wDir}/results/drug_immunesig_${tumor_purity_method}_GSEA_result/${cancer}
outPath=results/drug_immunesig_${tumor_purity_method}_GSEA_result/${cancer}/${dataset}_${cell}

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


# echo ${wDir}${outPath}
i=0
for file in `ls ${wDir}/${exp}`
do
  # filelist[$i]="$file" ${tmp%.*}
  outDir[$i]=${wDir}/${outPath}/${sign}_${num_gene}/${file%.*}
  exp_path[$i]=${wDir}/${exp}/$file
  grp_path[$i]=${wDir}/${grp}/${file%.*}.cls
  ((i++))
done


for ((i=0;i<=${#exp_path[@]};i++))
# for ((i=0;i<=100;i++))
do 
    # while ((i > 100))
    # do
    res=${exp_path[$i]}
    cls="${grp_path[$i]}#${pairs}"
    gmt=${wDir}/${geneset}
    
    out=${outDir[$i]}
    if [ ! -d ${out} ]; then
        mkdir ${out}
    fi

    nohup /picb/extprog/bin/linux-x86/java \
    -cp ${software}/gsea-3.0.jar \
    -Xmx512m xtools.gsea.Gsea  \
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
    wait
done

wait


# echo "out/n"
# echo ${outDir[1]}
# echo "exp/n"
# echo ${exp_path[1]}
# echo "grp/n"
# echo ${grp_path[1]}

resDir=${wDir}/${outPath}/${sign}_${num_gene} #/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/03_drug_immuneSig_enrichment/results/drug_immunesig_IHC_GSEA_result/${cancer}/${sign}_${num_gene}

mkdir /picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/03_drug_immuneSig_enrichment/results/drug_immunesig_IHC_GSEA_result_files/
desDir=/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/03_drug_immuneSig_enrichment/results/drug_immunesig_IHC_GSEA_result_files/${cancer}

if [ ! -d ${desDir} ]; then
    mkdir ${desDir}
    mkdir ${desDir}/${dataset}_${cell}
    mkdir ${desDir}/${dataset}_${cell}/${sign}_${num_gene}

fi

desDir=${desDir}/${dataset}_${cell}/${sign}_${num_gene}

if [ ! -d ${desDir} ]; then
    mkdir ${desDir}
fi

for file in `ls ${resDir}`
do
  cp ${resDir}/${file}/treated_versus_DMSO.Gsea.*/gsea_report_for_treated_*.xls ${desDir}/gsea_report_for_treated_${file}.xls
  cp ${resDir}/${file}/treated_versus_DMSO.Gsea.*/gsea_report_for_DMSO_*.xls ${desDir}/gsea_report_for_DMSO_${file}.xls
done

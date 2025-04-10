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


# test example:
# cancer=SKCM
# sampleType=Metastatic
# method=spearman
# tumor_purity_method=TUMERIC
# sign=positive
# num_gene=200
# dataset=70138

wDir=/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/03_drug_immuneSig_enrichment

software=/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/03_drug_immuneSig_enrichment/scripts/GSEA


pairs=treated_versus_DMSO

s1=0.4
ITSproof=2
genevote=0.5
cluster_method=hclust
similarity_method=Otsuka_Ochiai

if [ ${tumor_purity_method} == CPE ]; then
  # # geneset=data/gene_immune_sig_${tumor_purity_method}/for_GSEA/cor/${method}_${cancer}_${sampleType}_${sign}_${num_gene}.gmt
  # # geneset=data/gene_immune_sig_cor/for_GSEA/${method}_${cancer}_${sampleType}_${sign}_${num_gene}.gmt
  # geneset=data/gene_immune_sig_TUMERIC_new/for_GSEA/cor/${method}_${cancer}_${sampleType}_${sign}_${num_gene}.gmt
  # geneset=data/ITS_${cluster_method}_${similarity_method}_${tumor_purity_method}_${method}_r_${s1}_ITSproof_${ITSproof}_genevote_${genevote}/for_GSEA/cor/${method}_${cancer}_${sampleType}_${sign}_${num_gene}.gmt
  # geneset=/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/03_drug_immuneSig_enrichment/data/ITS_mclust_Otsuka_Ochiai_TUMERIC_spearman_r_0.4_metaP_0.05_genevote_0.5/for_GSEA/cor/${method}_${cancer}_${sampleType}_${sign}_${num_gene}.gmt
  geneset=/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/02_tumorGene_immuneSig_correlation/results_selected/ITS_oe_original_spearman_CPE_r0.4/for_GSEA/pcor/${method}_${cancer}_${sampleType}_${sign}_${num_gene}.gmt
  # geneset=/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/03_drug_immuneSig_enrichment/data/ITS_hclust_Otsuka_Ochiai_TUMERIC_spearman_r_0.4_metaP_0.05_genevote_0.5/for_GSEA/cor/${method}_${cancer}_${sampleType}_${sign}_${num_gene}.gmt
else
  # # geneset=data/gene_immune_sig_${tumor_purity_method}/for_GSEA/pcor/${method}_${cancer}_${sampleType}_${sign}_${num_gene}.gmt
  # geneset=data/gene_immune_sig_${tumor_purity_method}_new/for_GSEA/pcor/${method}_${cancer}_${sampleType}_${sign}_${num_gene}.gmt
  # merged ITS
  # geneset=data/ITS_${cluster_method}_${similarity_method}_${tumor_purity_method}_${method}_r_${s1}_ITSproof_${ITSproof}_genevote_${genevote}/for_GSEA/pcor/${method}_${cancer}_${sampleType}_${sign}_${num_gene}.gmt
  # geneset=/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/02_tumorGene_immuneSig_correlation/results_selected/immunecell_TCGAgenes_result_${method}_${tumor_purity_method}_r${s1}/for_GSEA/pcor/${method}_${cancer}_${sampleType}_${sign}_${num_gene}.gmt
  # geneset=/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/03_drug_immuneSig_enrichment/data/ITS_mclust_Otsuka_Ochiai_TUMERIC_spearman_r_0.4_metaP_0.05_genevote_0.5/for_GSEA/pcor/${method}_${cancer}_${sampleType}_${sign}_${num_gene}.gmt
  geneset=/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/02_tumorGene_immuneSig_correlation/results_selected/ITS_oe_original_spearman_TUMERIC_r0.4/for_GSEA/pcor/${method}_${cancer}_${sampleType}_${sign}_${num_gene}.gmt
  # geneset=/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/03_drug_immuneSig_enrichment/data/ITS_hclust_Otsuka_Ochiai_TUMERIC_spearman_r_0.4_metaP_0.05_genevote_0.5/for_GSEA/pcor/${method}_${cancer}_${sampleType}_${sign}_${num_gene}.gmt
fi 

exp=data/drug_treated_expression_matrix/${dataset}/${cancer}/
grp=data/drug_treated_design/for_GSEA/${dataset}/${cancer}/

mkdir ${wDir}/results_${tumor_purity_method}_meta_oe_original/
mkdir ${wDir}/results_${tumor_purity_method}_meta_oe_original/drug_immunesig_${tumor_purity_method}_GSEA_result
mkdir ${wDir}/results_${tumor_purity_method}_meta_oe_original/drug_immunesig_${tumor_purity_method}_GSEA_result/${dataset}
outPath=results_${tumor_purity_method}_meta_oe_original/drug_immunesig_${tumor_purity_method}_GSEA_result/${dataset}/${cancer}_${sampleType}



mkdir /picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/03_drug_immuneSig_enrichment/results_${tumor_purity_method}_meta_oe_original/drug_immunesig_${tumor_purity_method}_GSEA_result_files
mkdir /picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/03_drug_immuneSig_enrichment/results_${tumor_purity_method}_meta_oe_original/drug_immunesig_${tumor_purity_method}_GSEA_result_files/${dataset}
mkdir /picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/03_drug_immuneSig_enrichment/results_${tumor_purity_method}_meta_oe_original/drug_immunesig_${tumor_purity_method}_GSEA_result_files/${dataset}/${cancer}_${sampleType}
desDir=/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/03_drug_immuneSig_enrichment/results_${tumor_purity_method}_meta_oe_original/drug_immunesig_${tumor_purity_method}_GSEA_result_files/${dataset}/${cancer}_${sampleType}


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

    nohup /picb/extprog/bin/linux-x86/java \
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


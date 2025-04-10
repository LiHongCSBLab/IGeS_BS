
# ------------------------------------------------------------------------------
# DESeq2 size factor normalized counts
# ------------------------------------------------------------------------------

#!/bin/bash

#$ -S /bin/bash
#$ -N gsea
#$ -pe make 1

cancer=SKCM
sampleType=Metastatic
method=spearman
tumor_purity_method=TUMERIC
num_gene=200
dataset=70138
pairs=treated_versus_vehicle

s1=0.4
ITSproof=2
 

wDir=/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/03_drug_immuneSig_enrichment
software=/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/03_drug_immuneSig_enrichment/scripts/GSEA


exp=/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/data/immune_combo_dataset/GSE149825/GSE149825RNASeq_ncounts_gsea.txt
grp=/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/data/immune_combo_dataset/GSE149825/GSE149825RNASeq_ncounts_gsea.cls

# mkdir ${wDir}/results_${tumor_purity_method}/
# mkdir ${wDir}/results_${tumor_purity_method}/drug_immunesig_${tumor_purity_method}_GSEA_result
# mkdir ${wDir}/results_${tumor_purity_method}/drug_immunesig_${tumor_purity_method}_GSEA_result/${dataset}

mkdir /picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/result_GSE149825
mkdir /picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/result_GSE149825/drug_immunesig_${tumor_purity_method}_GSEA_result_r${s1}_meta_oe_original
outPath=/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/result_GSE149825/drug_immunesig_${tumor_purity_method}_GSEA_result_r${s1}_meta_oe_original/



sign=positive


if [ ${tumor_purity_method} == none ]; then
  # geneset=data/gene_immune_sig_${tumor_purity_method}/for_GSEA/cor/${method}_${cancer}_${sampleType}_${sign}_${num_gene}.gmt
  # geneset=data/gene_immune_sig_${tumor_purity_method}_new/for_GSEA/cor/${method}_${cancer}_${sampleType}_${sign}_${num_gene}.gmt
  geneset=/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/02_tumorGene_immuneSig_correlation/results_selected/ITS_oe_original_spearman_CPE_r0.4/for_GSEA/pcor/${method}_${cancer}_${sampleType}_${sign}_${num_gene}.gmt
else
#   geneset=data/gene_immune_sig_${tumor_purity_method}_${method}_r${s1}_meta_oe_original/for_GSEA/pcor/${method}_${cancer}_${sampleType}_${sign}_${num_gene}.gmt
  # geneset=data/gene_immune_sig_${tumor_purity_method}_new/for_GSEA/pcor/${method}_${cancer}_${sampleType}_${sign}_${num_gene}.gmt
#   geneset=/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/02_tumorGene_immuneSig_correlation/results_selected/ITS_oe_original_${method}_${tumor_purity_method}_r${s1}/for_GSEA/pcor/${method}_${cancer}_${sampleType}_${sign}_${num_gene}.gmt
  geneset=/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/02_tumorGene_immuneSig_correlation/results_selected/ITS_oe_original_spearman_TUMERIC_r0.4/for_GSEA/pcor/${method}_${cancer}_${sampleType}_${sign}_${num_gene}.gmt
fi

gmt=${geneset}


if [ ! -d ${outPath} ]; then
    mkdir ${outPath}
fi
if [ ! -d ${outPath}/${sign}_${num_gene} ]; then
    mkdir ${outPath}/${sign}_${num_gene}
fi
if [ ! -d ${outPath}/${sign}_${num_gene} ]; then
    mkdir ${outPath}/${sign}_${num_gene}
fi

outDir=${outPath}/${sign}_${num_gene}/

nohup /picb/extprog/bin/linux-x86/java \
    -cp ${software}/gsea-3.0.jar \
    -Xmx512m xtools.gsea.Gsea  \
    -gmx "${gmt}" \
    -res "${exp}" \
    -cls "${grp}" \
    -collapse false \
    -mode Max_probe \
    -norm meandiv \
    -nperm 1000  \
    -scoring_scheme weighted \
    -permute gene_set \
    -rpt_label ${pairs} \
    -out "${outDir}" \
    -include_only_symbols true \
    -make_sets true \
    -plot_top_x 0 \
    -rnd_seed timestamp \
    -set_max 600 \
    -set_min 5 \
    -zip_report false \
    -gui false > GSE149825.log 2>& 1 &


sign=negative


if [ ${tumor_purity_method} == none ]; then
  # geneset=data/gene_immune_sig_${tumor_purity_method}/for_GSEA/cor/${method}_${cancer}_${sampleType}_${sign}_${num_gene}.gmt
  # geneset=data/gene_immune_sig_${tumor_purity_method}_new/for_GSEA/cor/${method}_${cancer}_${sampleType}_${sign}_${num_gene}.gmt
  geneset=/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/02_tumorGene_immuneSig_correlation/results_selected/ITS_oe_original_spearman_CPE_r0.4/for_GSEA/pcor/${method}_${cancer}_${sampleType}_${sign}_${num_gene}.gmt
else
#   geneset=data/gene_immune_sig_${tumor_purity_method}_${method}_r${s1}_meta_oe_original/for_GSEA/pcor/${method}_${cancer}_${sampleType}_${sign}_${num_gene}.gmt
  # geneset=data/gene_immune_sig_${tumor_purity_method}_new/for_GSEA/pcor/${method}_${cancer}_${sampleType}_${sign}_${num_gene}.gmt
#   geneset=/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/02_tumorGene_immuneSig_correlation/results_selected/ITS_oe_original_${method}_${tumor_purity_method}_r${s1}/for_GSEA/pcor/${method}_${cancer}_${sampleType}_${sign}_${num_gene}.gmt
  geneset=/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/02_tumorGene_immuneSig_correlation/results_selected/ITS_oe_original_spearman_TUMERIC_r0.4/for_GSEA/pcor/${method}_${cancer}_${sampleType}_${sign}_${num_gene}.gmt
fi

gmt=${geneset}

# gmt=${wDir}/${geneset}

if [ ! -d ${outPath} ]; then
    mkdir ${outPath}
fi
if [ ! -d ${outPath}/${sign}_${num_gene} ]; then
    mkdir ${outPath}/${sign}_${num_gene}
fi
if [ ! -d ${outPath}/${sign}_${num_gene} ]; then
    mkdir ${outPath}/${sign}_${num_gene}
fi

outDir=${outPath}/${sign}_${num_gene}/

nohup /picb/extprog/bin/linux-x86/java \
    -cp ${software}/gsea-3.0.jar \
    -Xmx512m xtools.gsea.Gsea  \
    -gmx "${gmt}" \
    -res "${exp}" \
    -cls "${grp}" \
    -collapse false \
    -mode Max_probe \
    -norm meandiv \
    -nperm 1000  \
    -scoring_scheme weighted \
    -permute gene_set \
    -rpt_label ${pairs} \
    -out "${outDir}" \
    -include_only_symbols true \
    -make_sets true \
    -plot_top_x 0 \
    -rnd_seed timestamp \
    -set_max 600 \
    -set_min 5 \
    -zip_report false \
    -gui false > GSE149825.log 2>& 1 &



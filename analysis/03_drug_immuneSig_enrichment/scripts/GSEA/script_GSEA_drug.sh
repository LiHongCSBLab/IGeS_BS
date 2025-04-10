#!/bin/bash

#$ -S /bin/bash
#$ -N gsea
#$ -pe make 1

wDir=/picb/bigdata/project/FengFYM/immune_targeted_combo/immunecell_DEG_cor/GSEA_GSVA
cancer=$1	# TCGA cancer names
cell=$2		# cancer cell line names
drug=$3  	# drug names
dataset=$4	# GSE datasets
treatment_time=$5 # treatment time
# cancer=LIHC 
num_gene=100

exp=drug_treated_expression_matrix/${drug}_${dataset}_${cell}_${treatment_time}.txt
grp=drug_treated_design/for_GSEA/${drug}_${dataset}_${cell}_${treatment_time}.cls
pairs=treated_versus_DMSO
software=/picb/bigdata/project/FengFYM/immune_targeted_combo/immunecell_DEG_cor/scripts/immunesig_drug_GSEA_GSVA/GSEA

sign=negative 
geneset=gene_immune_sig/for_GSEA/${cancer}_${sign}_${num_gene}.gmt
outDir=/picb/bigdata/project/FengFYM/immune_targeted_combo/immunecell_DEG_cor/GSEA_GSVA/drug_immunesig_GSEA_result/${cancer}_${sign}_${num_gene}

res=${wDir}/${exp}
cls="${wDir}/${grp}#${pairs}"
gmt=${wDir}/${geneset}

if [ ! -d ${outDir} ]; then
	mkdir ${outDir}
fi

/picb/extprog/bin/linux-x86/java \
-cp ${software}/gsea-3.0.jar \
-Xmx8g xtools.gsea.Gsea  \
-gmx ${gmt} \
-res ${res} \
-cls ${cls} \
-collapse false \
-mode Max_probe \
-norm meandiv \
-nperm 1000  \
-scoring_scheme weighted \
-permute gene_set \
-rpt_label ${pairs} \
-out ${outDir} \
-include_only_symbols true \
-make_sets true \
-plot_top_x 0 \
-rnd_seed timestamp \
-set_max 600 \
-set_min 5 \
-zip_report false \
-gui false


sign=positive 
geneset=gene_immune_sig/for_GSEA/${cancer}_${sign}_${num_gene}.gmt
outDir=/picb/bigdata/project/FengFYM/immune_targeted_combo/immunecell_DEG_cor/GSEA_GSVA/drug_immunesig_GSEA_result/${cancer}_${sign}_${num_gene}

res=${wDir}/${exp}
cls="${wDir}/${grp}#${pairs}"
gmt=${wDir}/${geneset}

if [ ! -d ${outDir} ]; then
	mkdir ${outDir}
fi

/picb/extprog/bin/linux-x86/java \
-cp ${software}/gsea-3.0.jar \
-Xmx8g xtools.gsea.Gsea  \
-gmx ${gmt} \
-res ${res} \
-cls ${cls} \
-collapse false \
-mode Max_probe \
-norm meandiv \
-nperm 1000  \
-scoring_scheme weighted \
-permute gene_set \
-rpt_label ${pairs} \
-out ${outDir} \
-include_only_symbols true \
-make_sets true \
-plot_top_x 0 \
-rnd_seed timestamp \
-set_max 600 \
-set_min 5 \
-zip_report false \
-gui false

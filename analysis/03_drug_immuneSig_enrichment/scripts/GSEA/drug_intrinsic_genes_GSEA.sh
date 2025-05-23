#!/bin/bash

#$ -S /bin/bash
#$ -N gsea
#$ -pe make 1
cancer=$1
sign=$2
num_gene=$3
#sign=all
#cancer=COAD_MSI_H
#num_gene=400

wDir=/picb/bigdata/project/FengFYM/immune_targeted_combo/immunecell_DEG_cor/GSEA_GSVA

software=/picb/bigdata/project/FengFYM/immune_targeted_combo/immunecell_DEG_cor/scripts/immunesig_drug_GSEA_GSVA/GSEA


pairs=treated_versus_DMSO
geneset=gene_immune_sig/cancer_intrinic_genes.gmt

exp=drug_treated_expression_matrix/${cancer}
grp=drug_treated_design/for_GSEA/${cancer}
#mkdir drug_cancer_intrinic_genes_result
outPath=drug_cancer_intrinic_genes_result/${cancer}

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
do 
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
    -out ${out} \
    -include_only_symbols true \
    -make_sets true \
    -plot_top_x 0 \
    -rnd_seed timestamp \
    -set_max 600 \
    -set_min 5 \
    -zip_report false \
    -gui false > drug${i}.log 2>& 1 &
    wait

done

wait


# echo "out/n"
#echo ${outDir[1]}
# echo "exp/n"
# echo ${exp_path[1]}
# echo "grp/n"
# echo ${grp_path[1]}

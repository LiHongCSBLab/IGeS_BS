#!/bin/bash

#$ -S /bin/bash
#$ -N gsea
#$ -pe make 1

#wDir=/picb/bigdata/project/linping/yang/RNAseq/GSEA
wDir=$1
exp=$2		# expression matrix
grp=$3		# group info
pairs=$4  	# pairs to compare
geneset=$5	# which database to be used
outDir=$6

software=/home/linping/software

if [ $# -ne 6 ]; then
	echo "Usage: bash script.GSEA.sh [wDir] [exp] [grp] [pairsToCompare] [geneset] [outDir]"
	echo "Example: bash script.GSEA.sh wDir expr.txt grp.cls c1_versus_c2 hallmarks.gmt outDir"
	echo "Example: bash script.GSEA.sh wDir expr.txt grp.cls c0_versus_c1 hallmarks.gmt outDir"
	exit 1 
fi

#go=${wDir}/c5.bp.v6.0.symbols.gmt
#kegg=${wDir}/c2.cp.kegg.v6.0.symbols.gmt
res=${wDir}/${exp}
cls="${wDir}/${grp}#${pairs}"
gmt=${wDir}/${geneset}

#################
# set input file background DataBase, expMatrix, groupCls
#	-gmx ${gmt} -res ${res} -cls ${cls} \
################
# algorithm part 1
#	-collapse false -mode Max_probe -norm meandiv -nperm 1000  -scoring_scheme weighted -permute phenotype\
################
# set output
#	-rpt_label ${subDir}_${database} -out ${wDir}/All \
################ 
# algorithm part 2
#	-include_only_symbols true -make_sets true -plot_top_x 50 -rnd_seed timestamp -set_max 500 -set_min 5 -zip_report false -gui false
###############

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
-plot_top_x 200 \
-rnd_seed timestamp \
-set_max 600 \
-set_min 5 \
-zip_report false \
-gui false



#!/bin/bash
# run under miniconda environment r_env
source /home/fengfangyoumin/miniconda3/bin/activate r_env
cd /picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/03_drug_immuneSig_enrichment/scripts/


# cancer=SKCM
# sampleType=Metastatic
# tumor_purity_method=TUMERIC

# method=spearman
# dataset=92742

# sign2=positive
# num_gene2=200
# Rscript --vanilla 03_summarize_drug_immuneSig_GSEA.R -c ${cancer} -t ${sampleType} -p ${tumor_purity_method} -d ${dataset} -s ${sign2} -n ${num_gene2} -l allgenes > drug_summarize_${cancer}_${dataset}.log 2>&1 &
# wait

# sign3=negative
# num_gene3=200
# Rscript --vanilla 03_summarize_drug_immuneSig_GSEA.R -c ${cancer} -t ${sampleType} -p ${tumor_purity_method} -d ${dataset} -s ${sign3} -n ${num_gene3} -l allgenes > drug_summarize_${cancer}_${dataset}.log 2>&1 &
# wait

# cancer=BRCA
# sampleType=Primary
# tumor_purity_method=TUMERIC

# method=spearman
# dataset=92742

# sign2=positive
# num_gene2=200
# Rscript --vanilla 03_summarize_drug_immuneSig_GSEA.R -c ${cancer} -t ${sampleType} -p ${tumor_purity_method} -d ${dataset} -s ${sign2} -n ${num_gene2} -l allgenes > drug_summarize_${cancer}_${dataset}.log 2>&1 &
# wait

# sign3=negative
# num_gene3=200
# Rscript --vanilla 03_summarize_drug_immuneSig_GSEA.R -c ${cancer} -t ${sampleType} -p ${tumor_purity_method} -d ${dataset} -s ${sign3} -n ${num_gene3} -l allgenes > drug_summarize_${cancer}_${dataset}.log 2>&1 &
# wait


# cancer=COAD
# sign2=positive
# num_gene2=200
Rscript --vanilla 03_summarize_drug_immuneSig_GSEA.R -c COAD -t Primary -p TUMERIC -d 92742 -s positive -n 200 -l allgenes > drug_summarize_COAD_92742.log 2>&1 &
Rscript --vanilla 03_summarize_drug_immuneSig_GSEA.R -c COAD -t Primary -p TUMERIC -d 92742 -s negative -n 200 -l allgenes > drug_summarize_COAD_92742.log 2>&1 &
# wait

# sign3=negative
# num_gene3=200
# Rscript --vanilla 03_summarize_drug_immuneSig_GSEA.R -c ${cancer} -t ${sampleType} -p ${tumor_purity_method} -d ${dataset} -s ${sign3} -n ${num_gene3} -l allgenes > drug_summarize_${cancer}_${dataset}.log 2>&1 &
# wait

# cancer=LUAD
# sign2=positive
# num_gene2=200
Rscript --vanilla 03_summarize_drug_immuneSig_GSEA.R -c LUAD -t Primary -p TUMERIC -d 92742 -s positive -n 200 -l allgenes > drug_summarize_LUAD_92742.log 2>&1 &
Rscript --vanilla 03_summarize_drug_immuneSig_GSEA.R -c LUAD -t Primary -p TUMERIC -d 92742 -s negative -n 200 -l allgenes > drug_summarize_LUAD_92742.log 2>&1 &
# wait

# sign3=negative
# num_gene3=200
# Rscript --vanilla 03_summarize_drug_immuneSig_GSEA.R -c ${cancer} -t ${sampleType} -p ${tumor_purity_method} -d ${dataset} -s ${sign3} -n ${num_gene3} -l allgenes > drug_summarize_${cancer}_${dataset}.log 2>&1 &
# wait

cancer=PRAD
sign2=positive
num_gene2=200
Rscript --vanilla 03_summarize_drug_immuneSig_GSEA.R -c PRAD -t Primary -p TUMERIC -d 92742 -s positive -n 200 -l allgenes > drug_summarize_PRAD_92742.log 2>&1 &
Rscript --vanilla 03_summarize_drug_immuneSig_GSEA.R -c PRAD -t Primary -p TUMERIC -d 92742 -s negative -n 200 -l allgenes > drug_summarize_PRAD_92742.log 2>&1 &
# Rscript --vanilla 03_summarize_drug_immuneSig_GSEA.R -c ${cancer} -t ${sampleType} -p ${tumor_purity_method} -d ${dataset} -s ${sign2} -n ${num_gene2} -l allgenes > drug_summarize_${cancer}_${dataset}.log 2>&1 &
# # wait

# sign3=negative
# num_gene3=200
# Rscript --vanilla 03_summarize_drug_immuneSig_GSEA.R -c ${cancer} -t ${sampleType} -p ${tumor_purity_method} -d ${dataset} -s ${sign3} -n ${num_gene3} -l allgenes > drug_summarize_${cancer}_${dataset}.log 2>&1 &
# # wait

# cancer=READ
# sign2=positive
# num_gene2=200
Rscript --vanilla 03_summarize_drug_immuneSig_GSEA.R -c READ -t Primary -p TUMERIC -d 92742 -s positive -n 200 -l allgenes > drug_summarize_READ_92742.log 2>&1 &
Rscript --vanilla 03_summarize_drug_immuneSig_GSEA.R -c READ -t Primary -p TUMERIC -d 92742 -s negative -n 200 -l allgenes > drug_summarize_READ_92742.log 2>&1 &

# Rscript --vanilla 03_summarize_drug_immuneSig_GSEA.R -c ${cancer} -t ${sampleType} -p ${tumor_purity_method} -d ${dataset} -s ${sign2} -n ${num_gene2} -l allgenes > drug_summarize_${cancer}_${dataset}.log 2>&1 &
# # wait

# sign3=negative
# num_gene3=200
# Rscript --vanilla 03_summarize_drug_immuneSig_GSEA.R -c ${cancer} -t ${sampleType} -p ${tumor_purity_method} -d ${dataset} -s ${sign3} -n ${num_gene3} -l allgenes > drug_summarize_${cancer}_${dataset}.log 2>&1 &
# # wait

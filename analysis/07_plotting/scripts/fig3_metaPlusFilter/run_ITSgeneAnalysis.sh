#!/bin/bash

source /home/fengfangyoumin/miniconda3/bin/activate r_env
mkdir logfiles
# Rscript --vanilla fig3_ITS_preparation.R > logfiles/fig3_ITS_preparation.log 2>&1&

Rscript --vanilla fig3_ITS_characteristic_overall_enrich.R > logfiles/fig3_ITS_characteristic_overall_enrich.log 2>&1&
wait
Rscript --vanilla fig3_ITS_allgenerank_characteristic.R > logfiles/fig3_ITS_allgenerank_characteristic.log 2>&1&
wait
Rscript --vanilla fig3_ITS_allgenerank_characteristic_plot.R > logfiles/fig3_ITS_allgenerank_characteristic_plot.log 2>&1&
wait
Rscript --vanilla fig3_ITS_allgenerank_overallDT_plot.R > logfiles/fig3_ITS_allgenerank_overallDT_plot.log 2>&1&
wait
Rscript --vanilla fig3_ITS_allgenerank_drugEnrichonITS_plot.R > logfiles/fig3_ITS_allgenerank_drugEnrichonITS_plot.log 2>&1&
wait
Rscript --vanilla fig3_ITS_allgenerank_OGgsea_plot.R > logfiles/fig3_ITS_allgenerank_OGgsea_plot.log 2>&1&
wait
# Rscript --vanilla fig3_ITS_funcAnaysis.R > logfiles/fig3_ITS_funcAnaysis.log 2>&1&
Rscript --vanilla fig3_suppl_ITS_cancerspecific_characteristic.R > logfiles/fig3_suppl_ITS_cancerspecific_characteristic.log 2>&1&
wait
Rscript --vanilla fig3_ITS_similarity_analysis.R > logfiles/fig3_ITS_similarity_analysis.log 2>&1&
wait
Rscript --vanilla fig3_OriImmuSig_similarity_analysis.R > logfiles/fig3_OriImmuSig_similarity_analysis.log 2>&1&
wait
Rscript --vanilla fig3_suppl_ITSgene_cor_detection.R > logfiles/fig3_suppl_ITSgene_cor_detection.log 2>&1&
wait
Rscript --vanilla fig3_suppl_ITSgene_enrich.R > logfiles/fig3_suppl_ITSgene_enrich.log 2>&1&


Rscript --vanilla ITS_characteristic_2groups/fig3_ITS_gene_groups.R > logfiles/fig3_ITS_gene_groups.log 2>&1&
wait
Rscript --vanilla ITS_characteristic_2groups/fig3_ITS_sharedgene_characteristic.R > logfiles/fig3_ITS_sharedgene_characteristic.log 2>&1&
wait
Rscript --vanilla ITS_characteristic_2groups/fig3_ITSn_sharedgene_characteristic.R > logfiles/fig3_ITSn_sharedgene_characteristic.log 2>&1&
wait
Rscript --vanilla ITS_characteristic_2groups/fig3_ITS_sharedgene_drug_LINCS.R > logfiles/fig3_ITS_sharedgene_drug_LINCS.log 2>&1&
wait
Rscript --vanilla ITS_characteristic_2groups/fig3_ITS_sharedgene_drug_merged3DB.R > logfiles/fig3_ITS_sharedgene_drug_merged3DB.log 2>&1&
wait
Rscript --vanilla ITS_characteristic_2groups/fig3_ITS_sharedgene_drugAnalysis.R > logfiles/fig3_ITS_sharedgene_drugAnalysis.log 2>&1&



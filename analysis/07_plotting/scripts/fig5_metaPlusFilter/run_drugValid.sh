#!/bin/bash
source /home/fengfangyoumin/miniconda3/bin/activate r_env


Rscript --vanilla fig5_drug_summary_stat.R &
wait
Rscript --vanilla fig5_drugrank_meta_en_summary.R &
wait
Rscript --vanilla fig5_drugITS_summary_new_meta_en.R &
wait
Rscript --vanilla fig5_drugITS_meta_en_proof_enrich.R &
wait
Rscript --vanilla fig5_drugITS_meta_en_proof_enrich_plot.R &
wait
nohup Rscript --vanilla fig5_meta_en_rank_dotplot.R &
nohup Rscript --vanilla fig5_meta_en_rank_drugFDA_dotplot.R &
nohup Rscript --vanilla fig5_meta_en_rank_drugFDA_filterMOA_dotplot.R &
nohup Rscript --vanilla fig5_meta_en_rank_drugFDAcancer_dotplot.R &

nohup Rscript --vanilla fig5_drugrank_meta_en_diffTest_lincsPortal2.R &
nohup Rscript --vanilla fig5_drugrank_meta_en_alldrugFDA_diffTest_lincsPortal2.R &
nohup Rscript --vanilla fig5_drugrank_meta_en_alldrugFDAcancer_diffTest_lincsPortal2.R &
nohup Rscript --vanilla fig5_drugrank_meta_en_diffTest_isgold_lincsPortal2.R &
nohup Rscript --vanilla fig5_drugrank_meta_en_alldrugFDA_diffTest_isgold_lincsPortal2.R &


nohup Rscript --vanilla fig5_drugITS_shared_drug_meta_en.R &
nohup Rscript --vanilla fig5_drugITS_shared_drug_all_meta_en.R &
nohup Rscript --vanilla fig5_drugITS_shared_alldrugFDA_meta_en.R &
nohup Rscript --vanilla fig5_drugITS_shared_alldrugFDA_all_meta_en.R &

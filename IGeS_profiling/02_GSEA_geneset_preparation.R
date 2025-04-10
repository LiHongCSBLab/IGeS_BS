# Mapping drug DEGs to correlation results
print("start computing")
work_path <- "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/"

setwd(work_path)
options(stringsAsFactors = F)
library(parallel)
library(reshape2)
library(stringr)
library(dplyr)


source("03_drug_immuneSig_enrichment/scripts/02_GSEA_geneset_preparation_function.R")
source("03_drug_immuneSig_enrichment/scripts/02_GSEA_geneset_merge_function.R")
 
source("06_results_summaryandplotting/scripts/immunesig_selection_function.R")
source("06_results_summaryandplotting/scripts/ITS_analysis/ITS_function_detection.R")
# source("06_results_summaryandplotting/scripts/ITS_analysis/ITSexpr_purity_detection.R")
source("06_results_summaryandplotting/scripts/ITS_analysis/ITSexpr_purity_detection.R")
source("06_results_summaryandplotting/scripts/ITS_analysis/immuneModulatorDetection.R")
# tumor_purity_methods = c("IHC", "CPE", "TUMERIC", "ABSOLUTE")
# for(tumor_purity_method in tumor_purity_methods){
tumor_purity_method= "TUMERIC"
cor_method = "spearman"
immunesig_path = "06_results_summaryandplotting/data/immune_sig_selected/"
savepath_allresult = paste0("02_tumorGene_immuneSig_correlation/results_selected/")

# savepath_allresult = paste0("02_tumorGene_immuneSig_correlation/results_selected/immunecell_TCGAgenes_result_", cor_method,"_", tumor_purity_method, "_r0.4/")
dir.create(savepath_allresult)
savepath = paste0("03_drug_immuneSig_enrichment/data/")
# savepath = paste0("03_drug_immuneSig_enrichment/data/gene_immune_sig_", tumor_purity_method, "_",cor_method,"_r0.4/")
dir.create(savepath)
# dir.create(paste0(savepath, "for_GSVA/"))
# dir.create(paste0(savepath, "for_GSVA/cor"))
# dir.create(paste0(savepath, "for_GSVA/pcor"))
# dir.create(paste0(savepath, "for_GSEA/"))
# dir.create(paste0(savepath, "for_GSEA/cor"))
# dir.create(paste0(savepath, "for_GSEA/pcor"))

merged_ITS_savepath = paste0("03_drug_immuneSig_enrichment/data/")


resultpath = paste0("02_tumorGene_immuneSig_correlation/results/immunecell_TCGAgenes_result_",cor_method,"_", tumor_purity_method, "/")
filelist = list.files(resultpath)

filelist1 = filelist[grep('pcor', filelist)][-c(2:6, 9:11, 25, 26)]
filelist2 = filelist[-grep('pcor', filelist)][-8]


for (file1 in filelist1) { # [c(1, 14, 20, 21, 28, 29)]
  # file1 = "pcor_spearman_OV_Primary_TUMERIC.Rdata"
  # file1 = "pcor_spearman_SKCM_Metastatic_TUMERIC.Rdata"
  # file1 = "pcor_spearman_LUAD_Primary_TUMERIC.Rdata" 

  tmp = gsub("Rdata","", gsub("pcor_spearman_","",file1))
  tmp = gsub("\\.","",unlist(strsplit(tmp,split="_")))
  cancer = tmp[1]
  sample_type = tmp[2]
  gene_immune_pcor_sig(file_name = file1,
                       immunesig_path = immunesig_path,
                       resultpath = resultpath,
                       savepath_allresult = savepath_allresult,
                       savepath = savepath,
                       cancer = cancer,
                       sample_type = sample_type,
                       tumor_purity_method = tumor_purity_method,
                       method = cor_method,
                       s1 = 0.4,
                       s2 = -0.4,
                       padj = 0.05,
                      #  immunesig_confident_filter = 0.5,
                       auc_threshold = 0.6,
                       immuSigproof = 2, # how many datasets support the immune sigture to have prediction power
                       ITSproof = 2,
                       num_gene_p = 200,
                       num_gene_n = 200)

  # ITS_reduncancy_detection(cancer = cancer,
  #                          tumor_subtype = sample_type,
  #                          purity_method = "TUMERIC",
  #                          cor_method = "spearman",
  #                          num_gene = 200,
  #                          s = 0.4,
  #                          ITSproof = 2,
  #                          g_vote_threshold = 0.5,
  #                          ITSpath,
  #                          immuneSigInfoPath,
  #                          savepath = savepath)

  gene_immune_pcor_sig(file_name = file1,
                       immunesig_path = immunesig_path,
                       resultpath = resultpath,
                       savepath_allresult = savepath_allresult,
                       savepath = savepath,
                       cancer = cancer,
                       sample_type = sample_type,
                       tumor_purity_method = tumor_purity_method,
                       method = cor_method,
                       s1 = 0.3,
                       s2 = -0.3,
                       padj = 0.05,
                      #  immunesig_confident_filter = 0.5,
                       auc_threshold = 0.6,
                       immuSigproof = 2, # how many datasets support the immune sigture to have prediction power
                       ITSproof = 2,
                       num_gene_p = 200,
                       num_gene_n = 200)


} 


for (file1 in filelist1) {
    # file1 = "pcor_spearman_SKCM_Metastatic_TUMERIC.Rdata"

  tmp = gsub("Rdata","", gsub("pcor_spearman_","",file1))
  tmp = gsub("\\.","",unlist(strsplit(tmp,split="_")))
  cancer = tmp[1]
  sample_type = tmp[2]

  ITS_fileparation(savepath_allresult = savepath_allresult,
                                  cancer = cancer,
                                  sample_type = sample_type,
                                  method = cor_method,
                                  tumor_purity_method = tumor_purity_method,
                                  s1 = 0.4,
                                  num_gene_p = 200,
                                  num_gene_n = 200)

  ITS_fileparation(savepath_allresult = savepath_allresult,
                                  cancer = cancer,
                                  sample_type = sample_type,
                                  method = cor_method,
                                  tumor_purity_method = tumor_purity_method,
                                  s1 = 0.3,
                                  num_gene_p = 200,
                                  num_gene_n = 200)

}


for (file1 in filelist1) {
    # file1 = "pcor_spearman_SKCM_Metastatic_TUMERIC.Rdata"

  tmp = gsub("Rdata","", gsub("pcor_spearman_","",file1))
  tmp = gsub("\\.","",unlist(strsplit(tmp,split="_")))
  cancer = tmp[1]
  sample_type = tmp[2]

  ITS_fileparation_oe_original(savepath_allresult = savepath_allresult,
                                  cancer = cancer,
                                  sample_type = sample_type,
                                  method = cor_method,
                                  tumor_purity_method = tumor_purity_method,
                                  s1 = 0.4,
                                  num_gene_p = 200,
                                  num_gene_n = 200)

  ITS_fileparation_oe_original(savepath_allresult = savepath_allresult,
                                  cancer = cancer,
                                  sample_type = sample_type,
                                  method = cor_method,
                                  tumor_purity_method = tumor_purity_method,
                                  s1 = 0.3,
                                  num_gene_p = 200,
                                  num_gene_n = 200)

}



for (file1 in filelist1) {
    # file1 = "pcor_spearman_SKCM_Metastatic_TUMERIC.Rdata"

  tmp = gsub("Rdata","", gsub("pcor_spearman_","",file1))
  tmp = gsub("\\.","",unlist(strsplit(tmp,split="_")))
  cancer = tmp[1]
  sample_type = tmp[2]
  gene_immune_pcor_sig_finefilter(file_name = file1,
                                  immunesig_path = immunesig_path,
                                  resultpath = resultpath,
                                  savepath_allresult = savepath_allresult,
                                  savepath = savepath,
                                  cancer = cancer,
                                  sample_type = sample_type,
                                  tumor_purity_method = tumor_purity_method,
                                  method = cor_method,
                                  s1 = 0.4,
                                  s2 = -0.4,
                                  padj = 0.05,
                                  #  immunesig_confident_filter = 0.5,
                                  auc_threshold = 0.6,
                                  immuSigproof = 2, # how many datasets support the immune sigture to have prediction power
                                  ITSproof = 2,
                                  num_gene_p = 200,
                                  num_gene_n = 200)

  gene_immune_pcor_sig_finefilter(file_name = file1,
                                  immunesig_path = immunesig_path,
                                  resultpath = resultpath,
                                  savepath_allresult = savepath_allresult,
                                  savepath = savepath,
                                  cancer = cancer,
                                  sample_type = sample_type,
                                  tumor_purity_method = tumor_purity_method,
                                  method = cor_method,
                                  s1 = 0.3,
                                  s2 = -0.3,
                                  padj = 0.05,
                                  #  immunesig_confident_filter = 0.5,
                                  auc_threshold = 0.6,
                                  immuSigproof = 2, # how many datasets support the immune sigture to have prediction power
                                  ITSproof = 2,
                                  num_gene_p = 200,
                                  num_gene_n = 200)

}



tumor_purity_method= "CPE"
savepath_allresult = paste0("02_tumorGene_immuneSig_correlation/results_selected/")

# savepath_allresult = paste0("02_tumorGene_immuneSig_correlation/results_selected/immunecell_TCGAgenes_result_", cor_method,"_", tumor_purity_method, "_r0.4/")
dir.create(savepath_allresult)
savepath = paste0("03_drug_immuneSig_enrichment/data/")
# savepath = paste0("03_drug_immuneSig_enrichment/data/gene_immune_sig_", tumor_purity_method, "_",cor_method,"_r0.4/")
dir.create(savepath)
dir.create(paste0(savepath, "for_GSVA/"))
dir.create(paste0(savepath, "for_GSVA/cor"))
dir.create(paste0(savepath, "for_GSVA/pcor"))
dir.create(paste0(savepath, "for_GSEA/"))
dir.create(paste0(savepath, "for_GSEA/cor"))
dir.create(paste0(savepath, "for_GSEA/pcor"))


resultpath = paste0("02_tumorGene_immuneSig_correlation/results/immunecell_TCGAgenes_result_",cor_method,"_", tumor_purity_method, "/")
filelist = list.files(resultpath)

filelist1 = filelist[grep('pcor', filelist)][-c(3:7, 10:12, 25, 26)]
filelist2 = filelist[-grep('pcor', filelist)][-8]

for (file1 in filelist1) { # [c(1, 14, 20, 21, 28, 29)]
  # file1 = "pcor_spearman_SKCM_Primary_CPE.Rdata"   

  tmp = gsub("Rdata","", gsub("pcor_spearman_","",file1))
  tmp = gsub("\\.","",unlist(strsplit(tmp,split="_")))
  cancer = tmp[1]
  sample_type = tmp[2]
  gene_immune_pcor_sig(file_name = file1,
                       immunesig_path = immunesig_path,
                       resultpath = resultpath,
                       savepath_allresult = savepath_allresult,
                       savepath = savepath,
                       cancer = cancer,
                       sample_type = sample_type,
                       tumor_purity_method = tumor_purity_method,
                       method = cor_method,
                       s1 = 0.4,
                       s2 = -0.4,
                       padj = 0.05,
                      #  immunesig_confident_filter = 0.5,
                       auc_threshold = 0.6,
                       immuSigproof = 2, # how many datasets support the immune sigture to have prediction power
                       ITSproof = 2,
                       num_gene_p = 200,
                       num_gene_n = 200)

  # ITS_reduncancy_detection(cancer = cancer,
  #                          tumor_subtype = sample_type,
  #                          purity_method = "TUMERIC",
  #                          cor_method = "spearman",
  #                          num_gene = 200,
  #                          s = 0.4,
  #                          ITSproof = 2,
  #                          g_vote_threshold = 0.5,
  #                          ITSpath,
  #                          immuneSigInfoPath,
  #                          savepath = savepath)

  gene_immune_pcor_sig(file_name = file1,
                       immunesig_path = immunesig_path,
                       resultpath = resultpath,
                       savepath_allresult = savepath_allresult,
                       savepath = savepath,
                       cancer = cancer,
                       sample_type = sample_type,
                       tumor_purity_method = tumor_purity_method,
                       method = cor_method,
                       s1 = 0.3,
                       s2 = -0.3,
                       padj = 0.05,
                      #  immunesig_confident_filter = 0.5,
                       auc_threshold = 0.6,
                       immuSigproof = 2, # how many datasets support the immune sigture to have prediction power
                       ITSproof = 2,
                       num_gene_p = 200,
                       num_gene_n = 200)


} 


for (file1 in filelist1) {
  # file1 = "pcor_spearman_SKCM_Primary_CPE.Rdata"   

  tmp = gsub("Rdata","", gsub("pcor_spearman_","",file1))
  tmp = gsub("\\.","",unlist(strsplit(tmp,split="_")))
  cancer = tmp[1]
  sample_type = tmp[2]

  ITS_fileparation(savepath_allresult = savepath_allresult,
                                  cancer = cancer,
                                  sample_type = sample_type,
                                  method = cor_method,
                                  tumor_purity_method = tumor_purity_method,
                                  s1 = 0.4,
                                  num_gene_p = 200,
                                  num_gene_n = 200)

  ITS_fileparation(savepath_allresult = savepath_allresult,
                                  cancer = cancer,
                                  sample_type = sample_type,
                                  method = cor_method,
                                  tumor_purity_method = tumor_purity_method,
                                  s1 = 0.3,
                                  num_gene_p = 200,
                                  num_gene_n = 200)

}


for (file1 in filelist1) {
  # file1 = "pcor_spearman_SKCM_Primary_CPE.Rdata"   

  tmp = gsub("Rdata","", gsub("pcor_spearman_","",file1))
  tmp = gsub("\\.","",unlist(strsplit(tmp,split="_")))
  cancer = tmp[1]
  sample_type = tmp[2]

  ITS_fileparation_oe_original(savepath_allresult = savepath_allresult,
                                  cancer = cancer,
                                  sample_type = sample_type,
                                  method = cor_method,
                                  tumor_purity_method = tumor_purity_method,
                                  s1 = 0.4,
                                  num_gene_p = 200,
                                  num_gene_n = 200)

  ITS_fileparation_oe_original(savepath_allresult = savepath_allresult,
                                  cancer = cancer,
                                  sample_type = sample_type,
                                  method = cor_method,
                                  tumor_purity_method = tumor_purity_method,
                                  s1 = 0.3,
                                  num_gene_p = 200,
                                  num_gene_n = 200)

}



for (file1 in filelist1) {
  # file1 = "pcor_spearman_SKCM_Primary_CPE.Rdata"   

  tmp = gsub("Rdata","", gsub("pcor_spearman_","",file1))
  tmp = gsub("\\.","",unlist(strsplit(tmp,split="_")))
  cancer = tmp[1]
  sample_type = tmp[2]
  gene_immune_pcor_sig_finefilter(file_name = file1,
                                  immunesig_path = immunesig_path,
                                  resultpath = resultpath,
                                  savepath_allresult = savepath_allresult,
                                  savepath = savepath,
                                  cancer = cancer,
                                  sample_type = sample_type,
                                  tumor_purity_method = tumor_purity_method,
                                  method = cor_method,
                                  s1 = 0.4,
                                  s2 = -0.4,
                                  padj = 0.05,
                                  #  immunesig_confident_filter = 0.5,
                                  auc_threshold = 0.6,
                                  immuSigproof = 2, # how many datasets support the immune sigture to have prediction power
                                  ITSproof = 2,
                                  num_gene_p = 200,
                                  num_gene_n = 200)

  gene_immune_pcor_sig_finefilter(file_name = file1,
                                  immunesig_path = immunesig_path,
                                  resultpath = resultpath,
                                  savepath_allresult = savepath_allresult,
                                  savepath = savepath,
                                  cancer = cancer,
                                  sample_type = sample_type,
                                  tumor_purity_method = tumor_purity_method,
                                  method = cor_method,
                                  s1 = 0.3,
                                  s2 = -0.3,
                                  padj = 0.05,
                                  #  immunesig_confident_filter = 0.5,
                                  auc_threshold = 0.6,
                                  immuSigproof = 2, # how many datasets support the immune sigture to have prediction power
                                  ITSproof = 2,
                                  num_gene_p = 200,
                                  num_gene_n = 200)

}



# tumor_purity_method= "IHC"
# cor_method = "spearman"
# immunesig_path = "06_results_summaryandplotting/data/immune_sig_selected/"
# savepath_allresult = paste0("02_tumorGene_immuneSig_correlation/results_selected/immunecell_TCGAgenes_result_", cor_method,"_", tumor_purity_method, "/")
# dir.create(savepath_allresult)
# savepath = paste0("03_drug_immuneSig_enrichment/data/gene_immune_sig_", tumor_purity_method, "_",cor_method,"/")
# dir.create(savepath)
# dir.create(paste0(savepath, "for_GSVA/"))
# dir.create(paste0(savepath, "for_GSVA/cor"))
# dir.create(paste0(savepath, "for_GSVA/pcor"))
# dir.create(paste0(savepath, "for_GSEA/"))
# dir.create(paste0(savepath, "for_GSEA/cor"))
# dir.create(paste0(savepath, "for_GSEA/pcor"))


# resultpath = paste0("02_tumorGene_immuneSig_correlation/results/immunecell_TCGAgenes_pcor_",cor_method,"_", tumor_purity_method, "/")
# filelist = list.files(resultpath)

# filelist1 = filelist[grep('pcor', filelist)]
# filelist2 = filelist[-grep('pcor', filelist)]
# for (file1 in filelist1) {
#   tmp = gsub("Rdata","", gsub("pcor_spearman_","",file1))
#   tmp = unlist(strsplit(tmp,split="_"))
#   cancer = tmp[1]
#   sample_type = tmp[2]
#   gene_immune_pcor_sig(file_name = file1,
#                        immunesig_path = immunesig_path,
#                        resultpath = resultpath,
#                        savepath_allresult = savepath_allresult,
#                        savepath=savepath,
#                        cancer = cancer,
#                        sample_type = sample_type,
#                        tumor_purity_method = tumor_purity_method,
#                        method = cor_method,
#                        s1 = 0.4,
#                        s2 = -0.4,
#                        immunesig_confident_filter = 0.5,
#                        num_gene_p = 200,
#                        num_gene_n = 200 )
  
# }

# for (file2 in filelist2) {
#   tmp = gsub("Rdata","", gsub("pcor_spearman_","",file2))
#   tmp = unlist(strsplit(tmp,split="_"))
#   cancer = tmp[1]
#   sample_type = tmp[2]
#   gene_immune_pcor_sig(file_name = file2,
#                        immunesig_path = immunesig_path,
#                        resultpath = resultpath,
#                        savepath_allresult = savepath_allresult,
#                        savepath=savepath,
#                        cancer = cancer,
#                        sample_type = sample_type,
#                        tumor_purity_method = tumor_purity_method,
#                        method = cor_method,
#                        s1 = 0.4,
#                        s2 = -0.4,
#                        immunesig_confident_filter = 0.5,
#                        num_gene_p = 200,
#                        num_gene_n = 200 )
# }
# # }


# rm(list=ls())
work_path <- "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/"

setwd(work_path)
options(stringsAsFactors = F)
dir.create("07_plotting_v2/05_drugITSprofiler_CombatRef_metaPlusFilter/drugITSres_metaPlusFilter/")
# inst_GSE70138_is_gold <- read.csv("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/06_results_summaryandplotting/results_meta_filter_ITSssgsea_CombatRef_17datasets/data/LINCS_is_gold/inst_GSE70138_is_gold.csv", header = T, sep = '\t')

filepath <- "06_results_summaryandplotting/results_meta_filter_ITSssgsea_CombatRef_17datasets/TUMERIC_m2_meta_en_ACATP0.05/70138/allgenes/"
resfilepath <- list.files(filepath)

filenames <- list.files(paste0(filepath, "SKCM_Metastatic/drugITSRank/"))
filenames <- filenames[grep('.txt', filenames)]
# filenames <- filenames[-grep('unweight', filenames)]

for(f in filenames){
  #  f=filenames[3]
  filepath <- "06_results_summaryandplotting/results_meta_filter_ITSssgsea_CombatRef_17datasets/TUMERIC_m2_meta_en_ACATP0.05/70138/allgenes/"
  resfilepath <- list.files(filepath)

  res70138 <- list()
  for(i in 1:length(resfilepath)){
    # i=8
    # res_weighted_en <- read.csv(paste0(filepath, resfilepath[i], "/weighted/glmnet_weight/glmnet_weight_analysis_res.txt"), header = T, sep = '\t')
    # res_weighted_en$rank_tile <- res_weighted_en$rank / nrow(res_weighted_en)
    # res_weighted_en$cancer <- paste0(resfilepath[i],"70138")
    res_weighted_en <- read.csv(paste0(filepath, resfilepath[i], "/drugITSRank/", f), header = T, sep = '\t')
    res_weighted_en$rank_tile <- res_weighted_en$rank / nrow(res_weighted_en)
    res_weighted_en$cancer <- resfilepath[i]
    res_weighted_en$dataset <- "70138"
    
    res70138[[i]] <- res_weighted_en
  }
  res70138 <- do.call(rbind, res70138)
  res70138_SKCMprimary <-  read.csv(paste0("06_results_summaryandplotting/results_meta_filter_ITSssgsea_CombatRef_17datasets/CPE_m2_meta_en_ACATP0.05/70138/allgenes/SKCM_Primary/drugITSRank/", f), header = T, sep = '\t')
  res70138_SKCMprimary$rank_tile <- res70138_SKCMprimary$rank / nrow(res70138_SKCMprimary)
  res70138_SKCMprimary$cancer <- "SKCM_Primary"
  res70138_SKCMprimary$dataset = '70138'
  res70138 <- rbind(res70138, res70138_SKCMprimary)
  
  write.table(res70138, paste0("07_plotting_v2/05_drugITSprofiler_CombatRef_metaPlusFilter/drugITSres_metaPlusFilter/meta_en_ACATP0.05_70138_", f), sep = '\t', quote = F, row.names = F)
  # }
  
  
  filepath <- "06_results_summaryandplotting/results_meta_filter_ITSssgsea_CombatRef_17datasets/TUMERIC_m2_meta_en_ACATP0.05/92742/allgenes/"
  resfilepath <- list.files(filepath)
  
  # filenames <- list.files(paste0(filepath, "SKCM_Metastatic/drugITSRank/"))
  # filenames <- filenames[grep('.txt', filenames)]
  # filenames <- filenames[-grep('unweight', filenames)]
  
  # for(f in filenames){
  #  f=filenames[1]
  res92742 <- list()
  for(i in 1:length(resfilepath)){
    # i=9
    # res_weighted_en <- read.csv(paste0(filepath, resfilepath[i], "/weighted/glmnet_weight/glmnet_weight_analysis_res.txt"), header = T, sep = '\t')
    # res_weighted_en$rank_tile <- res_weighted_en$rank / nrow(res_weighted_en)
    # res_weighted_en$cancer <- paste0(resfilepath[i],"92742")
    res_weighted_en <- read.csv(paste0(filepath, resfilepath[i], "/drugITSRank/", f), header = T, sep = '\t')
    res_weighted_en$rank_tile <- res_weighted_en$rank / nrow(res_weighted_en)
    res_weighted_en$cancer <- resfilepath[i]
    res_weighted_en$dataset <- "92742"
    res92742[[i]] <- res_weighted_en
  }
  res92742 <- do.call(rbind, res92742)
  res92742_SKCMprimary <-  read.csv(paste0("06_results_summaryandplotting/results_meta_filter_ITSssgsea_CombatRef_17datasets/CPE_m2_meta_en_ACATP0.05/92742/allgenes/SKCM_Primary/drugITSRank/", f), header = T, sep = '\t')
  res92742_SKCMprimary$rank_tile <- res92742_SKCMprimary$rank / nrow(res92742_SKCMprimary)
  res92742_SKCMprimary$cancer <- "SKCM_Primary"
  res92742_SKCMprimary$dataset = '92742'
  
  res92742 <- rbind(res92742, res92742_SKCMprimary)
  
  write.table(res92742, paste0("07_plotting_v2/05_drugITSprofiler_CombatRef_metaPlusFilter/drugITSres_metaPlusFilter/meta_en_ACATP0.05_92742_", f), sep = '\t', quote = F, row.names = F)
  
  
  res <- rbind(res70138, res92742)
  write.table(res, paste0("07_plotting_v2/05_drugITSprofiler_CombatRef_metaPlusFilter/drugITSres_metaPlusFilter/meta_en_ACATP0.05_",f), sep = '\t', quote = F, row.names = F)
  # write.table(res[res$drug_index == "drug90",], "07_plotting_v2/05_drugITSprofiler_CombatRef_metaPlusFilter/crizotinib_rank/crizotinib_rank_ITSprofiler_meta_en_ACATP0.05.txt", sep = '\t', quote = F, row.names = F)
  
}

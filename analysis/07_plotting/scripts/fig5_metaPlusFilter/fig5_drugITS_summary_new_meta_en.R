# rm(list=ls())
work_path <- "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/"

setwd(work_path)
options(stringsAsFactors = F)



drug_proof <- read.csv("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/data/FDAapproved_p3Trial_IOtherapy/drug_list_proof_cancertype.txt", header = T, sep = '\t')

druglist <- unique(drug_proof[drug_proof$if_with_ICB == 'yes', ]$drug_index)
druglist <- druglist[!is.na(druglist)]


resfilepath <- list.files("06_results_summaryandplotting/results_meta_filter_ITSssgsea_CombatRef_17datasets/")
resfilepath <- resfilepath[grep("m2_meta_", resfilepath)]
resfilepath <- resfilepath[grep("ACATP0.05", resfilepath)]
resfilepath <- gsub("CPE_","",resfilepath)
resfilepath <- unique(gsub("TUMERIC_","",resfilepath))

for(respath in resfilepath){
  
  # respath=resfilepath[1]
    dir.create(paste0("07_plotting_v2/05_drugITSprofiler_CombatRef_metaPlusFilter/", respath, "_metaPlusFilter/"))
  cancerlist1 <- list.files(paste0("06_results_summaryandplotting/results_meta_filter_ITSssgsea_CombatRef_17datasets/TUMERIC_m2_meta_en_ACATP0.05/70138/allgenes/"))
  cancerlist1 <- c(cancerlist1, "SKCM_Primary")
  
  res_summary_70138 = list()
  
  for(n in seq(length(cancerlist1))){
    # n = 1 
    cancer = cancerlist1[n]
    
    # i = 1 
    
    if(cancer == "SKCM_Primary"){
      
      
      respath1 = paste0("CPE_", respath)
      filenames <- list.files(paste0("06_results_summaryandplotting/results_meta_filter_ITSssgsea_CombatRef_17datasets/",respath1,"/70138/allgenes/", cancer, "/drugITSRank/"))
      filenames <- filenames[grep('.txt',filenames)]
      if(length(grep("unweight",  filenames))>0)filenames <- filenames[-grep('unweight',filenames)]
      
      df1 <- lapply(as.list(filenames), function(f){
        res <- read.csv(paste0("06_results_summaryandplotting/results_meta_filter_ITSssgsea_CombatRef_17datasets/",respath1,"/70138/allgenes/", cancer, "/drugITSRank/", f), header = T, sep = '\t')
        res$rank_tile <- res$rank / nrow(res)
        res2 = res[is.element(res$drug_index, druglist), ]
        res2$strategy = gsub("meta_en_ACATP0.05_","",  gsub("_res.txt","", paste0( "",respath1,"_", f)))
        return(res2)
      })
      res_summary <- do.call(rbind, df1)
      res_summary$cancer_type <- strsplit(cancer,'_')[[1]][1]
      res_summary$cancer_subtype <- strsplit(cancer,'_')[[1]][2]
      
      res_summary_70138[[n]] <- res_summary[order(res_summary$drug_index), ]
      
    }else{
      
      respath2 = paste0("TUMERIC_", respath)
      filenames <- list.files(paste0("06_results_summaryandplotting/results_meta_filter_ITSssgsea_CombatRef_17datasets/",respath2,"/70138/allgenes/", cancer, "/drugITSRank/"))
      filenames <- filenames[grep('.txt',filenames)]
      if(length(grep("unweight",  filenames))>0)filenames <- filenames[-grep('unweight',filenames)]
      
      df1 <- lapply(as.list(filenames), function(f){
        res <- read.csv(paste0("06_results_summaryandplotting/results_meta_filter_ITSssgsea_CombatRef_17datasets/",respath2,"/70138/allgenes/", cancer, "/drugITSRank/", f), header = T, sep = '\t')
        res$rank_tile <- res$rank / nrow(res)
        res2 = res[is.element(res$drug_index, druglist), ]
        res2$strategy = gsub("meta_en_ACATP0.05_","",  gsub("_res.txt","", paste0( "",respath2,"_", f)))
        
        return(res2)
      })
      res_summary<- do.call(rbind, df1)
      res_summary$cancer_type <- strsplit(cancer,'_')[[1]][1]
      res_summary$cancer_subtype <- strsplit(cancer,'_')[[1]][2]
      res_summary_70138[[n]] <- res_summary[order(res_summary$drug_index), ]
    }
  }
  
  res_summary_70138 <- do.call(rbind, res_summary_70138)
  
  
  
  cancerlist2 <- list.files(paste0("06_results_summaryandplotting/results_meta_filter_ITSssgsea_CombatRef_17datasets/TUMERIC_m2_meta_en_ACATP0.05/92742/allgenes/"))
  cancerlist2 <- c(cancerlist2, "SKCM_Primary")
  
  res_summary_92742 = list()
  
  for(n in seq(length(cancerlist2))){
    # n = 1 
    cancer = cancerlist2[n]
    
    # i = 1 
    
    if(cancer == "SKCM_Primary"){
      
      
      respath1 = paste0("CPE_", respath)
      filenames <- list.files(paste0("06_results_summaryandplotting/results_meta_filter_ITSssgsea_CombatRef_17datasets/",respath1,"/92742/allgenes/", cancer, "/drugITSRank/"))
      filenames <- filenames[grep('.txt',filenames)]
      if(length(grep("unweight",  filenames))>0)filenames <- filenames[-grep('unweight',filenames)]
      
      df1 <- lapply(as.list(filenames), function(f){
        res <- read.csv(paste0("06_results_summaryandplotting/results_meta_filter_ITSssgsea_CombatRef_17datasets/",respath1,"/92742/allgenes/", cancer, "/drugITSRank/", f), header = T, sep = '\t')
        res$rank_tile <- res$rank / nrow(res)
        res2 = res[is.element(res$drug_index, druglist), ]
        res2$strategy = gsub("meta_en_ACATP0.05_","",  gsub("_res.txt","", paste0(  "",respath1,"_", f)))
        return(res2)
      })
      res_summary <- do.call(rbind, df1)
      res_summary$cancer_type <- strsplit(cancer,'_')[[1]][1]
      res_summary$cancer_subtype <- strsplit(cancer,'_')[[1]][2]
      
      res_summary_92742[[n]] <- res_summary[order(res_summary$drug_index), ]
      
    }else{
      
      respath2 = paste0("TUMERIC_", respath)
      
      filenames <- list.files(paste0("06_results_summaryandplotting/results_meta_filter_ITSssgsea_CombatRef_17datasets/",respath2,"/92742/allgenes/", cancer, "/drugITSRank/"))
      filenames <- filenames[grep('.txt',filenames)]
      if(length(grep("unweight",  filenames))>0)filenames <- filenames[-grep('unweight',filenames)]
      
      df1 <- lapply(as.list(filenames), function(f){
        res <- read.csv(paste0("06_results_summaryandplotting/results_meta_filter_ITSssgsea_CombatRef_17datasets/",respath2,"/92742/allgenes/", cancer, "/drugITSRank/", f), header = T, sep = '\t')
        res$rank_tile <- res$rank / nrow(res)
        res2 = res[is.element(res$drug_index, druglist), ]
        res2$strategy = gsub("meta_en_ACATP0.05_","",  gsub("_res.txt","", paste0( "",respath2,"_", f)))
        
        return(res2)
      })
      res_summary<- do.call(rbind, df1)
      res_summary$cancer_type <- strsplit(cancer,'_')[[1]][1]
      res_summary$cancer_subtype <- strsplit(cancer,'_')[[1]][2]
      res_summary_92742[[n]] <- res_summary[order(res_summary$drug_index), ]
    }
  }
  
  res_summary_92742 <- do.call(rbind, res_summary_92742)
  res_summary_70138$dataset = '70138'
  res_summary_92742$dataset = '92742'
  res_summary_all <- rbind(res_summary_70138, res_summary_92742)
  
  
  write.table(res_summary_all, paste0("07_plotting_v2/05_drugITSprofiler_CombatRef_metaPlusFilter/", respath, "_metaPlusFilter/res_summary_ACATP0.05_proof.txt"), sep = '\t', row.names = F, quote = F)
  
  drug_proof2 <- unique(drug_proof[-c(1,5,9)])
  names(drug_proof2)[7] <- 'cancer_type'
  res_summary_matchProof <- merge(res_summary_all, drug_proof2, by = c("drug_index","cancer_type"))
  
  write.table(res_summary_matchProof, paste0("07_plotting_v2/05_drugITSprofiler_CombatRef_metaPlusFilter/", respath, "_metaPlusFilter/res_summary_ACATP0.05_matchProof.txt"), sep = '\t', row.names = F, quote = F)
}

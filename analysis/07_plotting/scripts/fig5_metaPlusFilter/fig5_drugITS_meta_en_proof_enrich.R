# rm(list=ls())
work_path <- "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/"

setwd(work_path)
options(stringsAsFactors = F)
# dir.create("06_results_summaryandplotting/results_meta_filter_ITSssgsea_CombatRef_17datasets/res_summary/")
dir.create("07_plotting_v2/05_drugITSprofiler_CombatRef_metaPlusFilter/res_enrich_metaPlusFilter/")



library(org.Hs.eg.db)
library(topGO)
library(clusterProfiler)
library(msigdbr)
library(enrichplot)


drug_proof <- read.table("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/data/FDAapproved_p3Trial_IOtherapy/drug_list_proof_cancertype.txt", header = T, sep = '\t')

for(rankthreshold in rev(seq(0.1,1,0.1))) {
  # rankthreshold = 1

  resfilepath <- list.files("06_results_summaryandplotting/results_meta_filter_ITSssgsea_CombatRef_17datasets/")
  resfilepath <- resfilepath[grep("m2_meta_", resfilepath)]
  resfilepath <- resfilepath[grep("ACATP0.05", resfilepath)]
  resfilepath <- gsub("CPE_","",resfilepath)
  resfilepath <- unique(gsub("TUMERIC_","",resfilepath))

  for(respath in resfilepath){
    # respath = resfilepath[1]
    dir.create(paste0("07_plotting_v2/05_drugITSprofiler_CombatRef_metaPlusFilter/res_enrich_metaPlusFilter/",respath,"/"))
    cancerlist1 <- list.files(paste0("06_results_summaryandplotting/results_meta_filter_ITSssgsea_CombatRef_17datasets/TUMERIC_m2_meta_en_ACATP0.05/70138/allgenes/"))
    cancerlist1 <- c(cancerlist1, "SKCM_Primary")
    
    res_summary_70138 = list()
    
    for(n in seq(length(cancerlist1))){
      # n = 9

      cancer = cancerlist1[n]
      
      # i = 1 
      
      if(cancer == "SKCM_Primary"){
        
        respath1 = paste0("CPE_", respath)
        
        filenames <- list.files(paste0("06_results_summaryandplotting/results_meta_filter_ITSssgsea_CombatRef_17datasets/",respath1,"//70138/allgenes/", cancer, "/drugITSRank/"))
        filenames <- filenames[grep('.txt',filenames)]
       if(length(grep('unweighted',filenames)) != 0) filenames <- filenames[-grep('unweighted',filenames)]
        
        df1 <- lapply(as.list(filenames), function(f){
          # f=filenames[1]
          res <- read.csv(paste0("06_results_summaryandplotting/results_meta_filter_ITSssgsea_CombatRef_17datasets/",respath1,"//70138/allgenes/", cancer, "/drugITSRank/", f), header = T, sep = '\t')
          res$rank_tile <- res$rank / nrow(res)
          
          # res = res[res$u>0 & res$v < 0,]
          # res$rank2 <- seq(1, nrow(res), 1)
          # # res_score = as.vector(scale(res$rank_tile)) * -1
          res_score = as.vector((res$score) + 10e-10)
          names(res_score) = res$drug_index
          res_score = head(res_score, length(res_score)* rankthreshold)
          
          druglist <- unique(drug_proof[drug_proof$if_with_ICB == 'yes', ])
          druglist <- druglist[!is.na(druglist$drug_index), ]
          druglist <- druglist[druglist$proof != 'animal',]
          
          
          druglist1 =  data.frame(proof = 'proof', drug_index = unique(druglist$drug_index))
          
          res_enrich1 <- GSEA(geneList = res_score, 
                              pvalueCutoff = 1, 
                              minGSSize = 0,
                              maxGSSize = 10000,
                              TERM2GENE = druglist1)
          # data.frame(proof = 'proof', drug_index = druglist))
          res_enrich1_2 = as.data.frame(summary(res_enrich1))
          # res_enrich1_2
          
          druglist2 =  unique(druglist[c('Cancer_TCGA', 'drug_index')])
          res_enrich2 <- GSEA(geneList = res_score, 
                              pvalueCutoff = 1, 
                              minGSSize = 0,
                              maxGSSize = 10000,
                              TERM2GENE = druglist2)
          # data.frame(proof = 'proof', drug_index = druglist))
          
          res_enrich2_2 = as.data.frame(summary(res_enrich2))
          # res_enrich2_2
          enrich_ID = res_enrich2_2[res_enrich2_2$NES >0, ]
          enrich_ID = enrich_ID[order(enrich_ID$NES, decreasing = T), ]
          # gseaplot2(res_enrich2, geneSetID = enrich_ID$ID, pvalue_table = F, color = c("#E495A5", "#86B875", "#7DB0DD"))    
          # gseaplot2(res_enrich2, geneSetID = enrich_ID$ID, pvalue_table = F, color = c("#E495A5", "#86B875", "#7DB0DD"))    
          
          druglist <- unique(drug_proof[drug_proof$if_with_ICB == 'yes', ])
          druglist <- druglist[!is.na(druglist$drug_index), ]
          druglist3 =  unique(druglist[c('proof', 'drug_index')])
          res_enrich3 <- GSEA(geneList = res_score, 
                              pvalueCutoff = 1, 
                              minGSSize = 0,
                              maxGSSize = 10000,
                              TERM2GENE = druglist3)
          # data.frame(proof = 'proof', drug_index = druglist))
          res_enrich3_2 = as.data.frame(summary(res_enrich3))
          # res_enrich3_2
          druglist4 = druglist[is.element(druglist$proof, c("1","2","1-2","1/2a","1b/2","1/2","1b")), ]
          druglist4 =  data.frame(proof = 'ClinicalTrial_beforeIII', drug_index = unique(druglist4$drug_index))
          res_enrich4 <- GSEA(geneList = res_score, 
                              pvalueCutoff = 1, 
                              minGSSize = 0,
                              maxGSSize = 10000,
                              TERM2GENE = druglist4)
          # data.frame(proof = 'proof', drug_index = druglist))
          res_enrich4_2 = as.data.frame(summary(res_enrich4))
          res_enrich = rbind(res_enrich1_2, rbind(res_enrich2_2, res_enrich3_2, res_enrich4_2))
          res_enrich$strategy = gsub("meta_en_ACATP0.05_","",  gsub("_res.txt","", paste0( "",respath1,"/_", f)))
          return(res_enrich)
        })
        
        res_summary <- do.call(rbind, df1)
        res_summary$cancer_type <- strsplit(cancer,'_')[[1]][1]
        res_summary$cancer_subtype <- strsplit(cancer,'_')[[1]][2]
              
        res_summary_70138[[n]] <- res_summary
        
      }else{
        
        respath2 = paste0("TUMERIC_", respath)
        
        filenames <- list.files(paste0("06_results_summaryandplotting/results_meta_filter_ITSssgsea_CombatRef_17datasets/",respath2,"//70138/allgenes/", cancer, "/drugITSRank/"))
        filenames <- filenames[grep('.txt',filenames)]
       if(length(grep('unweighted',filenames)) != 0) filenames <- filenames[-grep('unweighted',filenames)]
        
        df1 <- lapply(as.list(filenames), function(f){
          #  f=filenames[4]
          # for(f in filenames){
          res <- read.csv(paste0("06_results_summaryandplotting/results_meta_filter_ITSssgsea_CombatRef_17datasets/",respath2,"//70138/allgenes/", cancer, "/drugITSRank/", f), header = T, sep = '\t')
          res$rank_tile <- res$rank / nrow(res)
          # res = res[res$u>0 & res$v < 0,]
          # res$rank2 <- seq(1, nrow(res), 1)

          # res_score = as.vector(scale(res$rank_tile)) * -1
          res_score = as.vector((res$score) + 10e-10)
          names(res_score) = res$drug_index
          res_score = head(res_score, length(res_score)* rankthreshold)
          
          druglist <- unique(drug_proof[drug_proof$if_with_ICB == 'yes', ])
          druglist <- druglist[!is.na(druglist$drug_index), ]
          druglist <- druglist[druglist$proof != 'animal',]
          
          
          druglist1 =  data.frame(proof = 'proof', drug_index = unique(druglist$drug_index))
          
          res_enrich1 <- GSEA(geneList = res_score, 
                              pvalueCutoff = 1, 
                              minGSSize = 0,
                              maxGSSize = 10000,
                              TERM2GENE = druglist1)
          # data.frame(proof = 'proof', drug_index = druglist))
          res_enrich1_2 = as.data.frame(summary(res_enrich1))
          # res_enrich1_2
          gseaplot2(res_enrich1, geneSetID = 1, pvalue_table = F, color = c("#E495A5", "#86B875", "#7DB0DD"))    

          druglist2 =  unique(druglist[c('Cancer_TCGA', 'drug_index')])
          res_enrich2 <- GSEA(geneList = res_score, 
                              pvalueCutoff = 1, 
                              minGSSize = 0,
                              maxGSSize = 10000,
                              TERM2GENE = druglist2)
          # data.frame(proof = 'proof', drug_index = druglist))
          
          res_enrich2_2 = as.data.frame(summary(res_enrich2))
          # res_enrich2_2
          enrich_ID = res_enrich2_2[res_enrich2_2$NES >0, ]
          enrich_ID = enrich_ID[order(enrich_ID$NES, decreasing = T), ]
          # gseaplot2(res_enrich2, geneSetID = enrich_ID$ID, pvalue_table = F, color = c("#E495A5", "#86B875", "#7DB0DD"))    
          gseaplot2(res_enrich2, geneSetID = 'SKCM', pvalue_table = F, color = c("#E495A5", "#86B875", "#7DB0DD"))    
          
          druglist <- unique(drug_proof[drug_proof$if_with_ICB == 'yes', ])
          druglist <- druglist[!is.na(druglist$drug_index), ]
          druglist3 =  unique(druglist[c('proof', 'drug_index')])
          res_enrich3 <- GSEA(geneList = res_score, 
                              pvalueCutoff = 1, 
                              minGSSize = 0,
                              maxGSSize = 10000,
                              TERM2GENE = druglist3)
          # data.frame(proof = 'proof', drug_index = druglist))
          res_enrich3_2 = as.data.frame(summary(res_enrich3))
          # res_enrich3_2
          
          druglist4 = druglist[is.element(druglist$proof, c("1","2","1-2","1/2a","1b/2","1/2","1b")), ]
          druglist4 =  data.frame(proof = 'ClinicalTrial_beforeIII', drug_index = unique(druglist4$drug_index))
          res_enrich4 <- GSEA(geneList = res_score, 
                              pvalueCutoff = 1, 
                              minGSSize = 0,
                              maxGSSize = 10000,
                              TERM2GENE = druglist4)
          # data.frame(proof = 'proof', drug_index = druglist))
          res_enrich4_2 = as.data.frame(summary(res_enrich4))
          res_enrich = rbind(res_enrich1_2, rbind(res_enrich2_2, res_enrich3_2, res_enrich4_2))
          res_enrich$strategy = gsub("meta_en_ACATP0.05_","",  gsub("_res.txt","", paste0("",respath2,"/_", f)))


          return(res_enrich)
        })
        
        res_summary<- do.call(rbind, df1)
        res_summary$cancer_type <- strsplit(cancer,'_')[[1]][1]
        res_summary$cancer_subtype <- strsplit(cancer,'_')[[1]][2]
        res_summary_70138[[n]] <- res_summary
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

        filenames <- list.files(paste0("06_results_summaryandplotting/results_meta_filter_ITSssgsea_CombatRef_17datasets/",respath1,"//92742/allgenes/", cancer, "/drugITSRank/"))
        filenames <- filenames[grep('.txt',filenames)]
       if(length(grep('unweighted',filenames)) != 0) filenames <- filenames[-grep('unweighted',filenames)]
        
        df1 <- lapply(as.list(filenames), function(f){
          res <- read.csv(paste0("06_results_summaryandplotting/results_meta_filter_ITSssgsea_CombatRef_17datasets/",respath1,"//92742/allgenes/", cancer, "/drugITSRank/", f), header = T, sep = '\t')
          res$rank_tile <- res$rank / nrow(res)

          res_score = as.vector((res$score) + 10e-10)
          names(res_score) = res$drug_index
          res_score = head(res_score, length(res_score) * rankthreshold)
          
          druglist <- unique(drug_proof[drug_proof$if_with_ICB == 'yes', ])
          druglist <- druglist[!is.na(druglist$drug_index), ]
          druglist <- druglist[druglist$proof != 'animal',]
          
          
          druglist1 =  data.frame(proof = 'proof', drug_index = unique(druglist$drug_index))
          
          res_enrich1 <- GSEA(geneList = res_score, 
                              pvalueCutoff = 1, 
                              minGSSize = 0,
                              maxGSSize = 10000,
                              TERM2GENE = druglist1)
          # data.frame(proof = 'proof', drug_index = druglist))
          res_enrich1_2 = as.data.frame(summary(res_enrich1))
          # res_enrich1_2
          
          druglist2 =  unique(druglist[c('Cancer_TCGA', 'drug_index')])
          res_enrich2 <- GSEA(geneList = res_score, 
                              pvalueCutoff = 1, 
                              minGSSize = 0,
                              maxGSSize = 10000,
                              TERM2GENE = druglist2)
          # data.frame(proof = 'proof', drug_index = druglist))
          
          res_enrich2_2 = as.data.frame(summary(res_enrich2))
          # res_enrich2_2
          enrich_ID = res_enrich2_2[res_enrich2_2$NES >0, ]
          enrich_ID = enrich_ID[order(enrich_ID$NES, decreasing = T), ]
          # gseaplot2(res_enrich2, geneSetID = enrich_ID$ID, pvalue_table = F, color = c("#E495A5", "#86B875", "#7DB0DD"))    
          
          druglist <- unique(drug_proof[drug_proof$if_with_ICB == 'yes', ])
          druglist <- druglist[!is.na(druglist$drug_index), ]
          druglist3 =  unique(druglist[c('proof', 'drug_index')])
          res_enrich3 <- GSEA(geneList = res_score, 
                              pvalueCutoff = 1, 
                              minGSSize = 0,
                              maxGSSize = 10000,
                              TERM2GENE = druglist3)
          # data.frame(proof = 'proof', drug_index = druglist))
          res_enrich3_2 = as.data.frame(summary(res_enrich3))
          # res_enrich3_2        
          druglist4 = druglist[is.element(druglist$proof, c("1","2","1-2","1/2a","1b/2","1/2","1b")), ]
          druglist4 =  data.frame(proof = 'ClinicalTrial_beforeIII', drug_index = unique(druglist4$drug_index))
          res_enrich4 <- GSEA(geneList = res_score, 
                              pvalueCutoff = 1, 
                              minGSSize = 0,
                              maxGSSize = 10000,
                              TERM2GENE = druglist4)
          # data.frame(proof = 'proof', drug_index = druglist))
          res_enrich4_2 = as.data.frame(summary(res_enrich4))
          res_enrich = rbind(res_enrich1_2, rbind(res_enrich2_2, res_enrich3_2, res_enrich4_2))
          res_enrich$strategy = gsub("meta_en_ACATP0.05_","",  gsub("_res.txt","", paste0(  "",respath1,"/_", f)))

          return(res_enrich)
        })
        res_summary <- do.call(rbind, df1)
        res_summary$cancer_type <- strsplit(cancer,'_')[[1]][1]
        res_summary$cancer_subtype <- strsplit(cancer,'_')[[1]][2]
        
        res_summary_92742[[n]] <- res_summary
        
      }else{
        
        respath1 = paste0("TUMERIC_", respath)

        filenames <- list.files(paste0("06_results_summaryandplotting/results_meta_filter_ITSssgsea_CombatRef_17datasets/",respath2,"//92742/allgenes/", cancer, "/drugITSRank/"))
        filenames <- filenames[grep('.txt',filenames)]
       if(length(grep('unweighted',filenames)) != 0) filenames <- filenames[-grep('unweighted',filenames)]
        
        df1 <- lapply(as.list(filenames), function(f){
          res <- read.csv(paste0("06_results_summaryandplotting/results_meta_filter_ITSssgsea_CombatRef_17datasets/",respath2,"//92742/allgenes/", cancer, "/drugITSRank/", f), header = T, sep = '\t')
          res$rank_tile <- res$rank / nrow(res)
          res_score = as.vector((res$score) + 10e-10)
          names(res_score) = res$drug_index
          res_score = head(res_score, length(res_score) * rankthreshold)
          
          druglist <- unique(drug_proof[drug_proof$if_with_ICB == 'yes', ])
          druglist <- druglist[!is.na(druglist$drug_index), ]
          druglist <- druglist[druglist$proof != 'animal',]
          
          
          druglist1 =  data.frame(proof = 'proof', drug_index = unique(druglist$drug_index))
          
          res_enrich1 <- GSEA(geneList = res_score, 
                              pvalueCutoff = 1, 
                              minGSSize = 0,
                              maxGSSize = 10000,
                              TERM2GENE = druglist1)
          # data.frame(proof = 'proof', drug_index = druglist))
          res_enrich1_2 = as.data.frame(summary(res_enrich1))
          # res_enrich1_2
          
          druglist2 =  unique(druglist[c('Cancer_TCGA', 'drug_index')])
          res_enrich2 <- GSEA(geneList = res_score, 
                              pvalueCutoff = 1, 
                              minGSSize = 0,
                              maxGSSize = 10000,
                              TERM2GENE = druglist2)
          # data.frame(proof = 'proof', drug_index = druglist))
          
          res_enrich2_2 = as.data.frame(summary(res_enrich2))
          # res_enrich2_2
          enrich_ID = res_enrich2_2[res_enrich2_2$NES >0, ]
          enrich_ID = enrich_ID[order(enrich_ID$NES, decreasing = T), ]
          # gseaplot2(res_enrich2, geneSetID = enrich_ID$ID, pvalue_table = F, color = c("#E495A5", "#86B875", "#7DB0DD"))    
          
          druglist <- unique(drug_proof[drug_proof$if_with_ICB == 'yes', ])
          druglist <- druglist[!is.na(druglist$drug_index), ]
          druglist3 =  unique(druglist[c('proof', 'drug_index')])
          res_enrich3 <- GSEA(geneList = res_score, 
                              pvalueCutoff = 1, 
                              minGSSize = 0,
                              maxGSSize = 10000,
                              TERM2GENE = druglist3)
          # data.frame(proof = 'proof', drug_index = druglist))
          res_enrich3_2 = as.data.frame(summary(res_enrich3))
          # res_enrich3_2
          druglist4 = druglist[is.element(druglist$proof, c("1","2","1-2","1/2a","1b/2","1/2","1b")), ]
          druglist4 =  data.frame(proof = 'ClinicalTrial_beforeIII', drug_index = unique(druglist4$drug_index))
          res_enrich4 <- GSEA(geneList = res_score, 
                              pvalueCutoff = 1, 
                              minGSSize = 0,
                              maxGSSize = 10000,
                              TERM2GENE = druglist4)
          # data.frame(proof = 'proof', drug_index = druglist))
          res_enrich4_2 = as.data.frame(summary(res_enrich4))
          res_enrich = rbind(res_enrich1_2, rbind(res_enrich2_2, res_enrich3_2, res_enrich4_2))
          res_enrich$strategy = gsub("meta_en_ACATP0.05_","",  gsub("_res.txt","", paste0("",respath2,"/_", f)))

          return(res_enrich)
          
        })
        res_summary<- do.call(rbind, df1)
        res_summary$cancer_type <- strsplit(cancer,'_')[[1]][1]
        res_summary$cancer_subtype <- strsplit(cancer,'_')[[1]][2]
        res_summary_92742[[n]] <- res_summary
      }
    }
    
    res_summary_92742 <- do.call(rbind, res_summary_92742)
    res_summary_70138$dataset = '70138'
    res_summary_92742$dataset = '92742'
    res_summary_all <- rbind(res_summary_70138, res_summary_92742)
    
    
    write.table(res_summary_all, paste0("07_plotting_v2/05_drugITSprofiler_CombatRef_metaPlusFilter/res_enrich_metaPlusFilter/",respath,"/",rankthreshold,".txt"), sep = '\t', row.names = F, quote = F)
    
  }
}


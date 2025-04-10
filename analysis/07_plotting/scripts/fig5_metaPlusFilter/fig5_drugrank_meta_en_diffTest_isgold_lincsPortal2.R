# rm(list=ls())
work_path <- "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/"

setwd(work_path)
options(stringsAsFactors = F)
set.seed(1234)
dir.create('07_plotting_v2/05_drugITSprofiler_CombatRef_metaPlusFilter/drugrank_diff_isgold_lincsPortal2_metaPlusFilter/')
dir.create('07_plotting_v2/05_drugITSprofiler_CombatRef_metaPlusFilter/drugrank_diff_isgold_lincsPortal2_metaPlusFilter/70138/')
dir.create('07_plotting_v2/05_drugITSprofiler_CombatRef_metaPlusFilter/drugrank_diff_isgold_lincsPortal2_metaPlusFilter/92742/')

inst_GSE70138_is_gold <- read.csv("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/06_results_summaryandplotting/data/LINCS_is_gold/inst_GSE70138_is_gold.csv", header = T, sep = '\t')
inst_GSE92742_is_gold <- read.csv("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/06_results_summaryandplotting/data//LINCS_is_gold/inst_GSE92742_is_gold.csv", header = T, sep = '\t')

cell_70138 <- read.csv("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/03_drug_immuneSig_enrichment/data/cell_annot_selected_70138.txt", sep = '\t', header = T)
cell_92742 <- read.csv("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/03_drug_immuneSig_enrichment/data/cell_annot_selected_92742.txt", sep = '\t', header = T)


drug_proof <- read.table("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/data/FDAapproved_p3Trial_IOtherapy/drug_list_proof_cancertype.txt", header = T, sep = '\t')
druglist <- unique(drug_proof[drug_proof$if_with_ICB == 'yes', ])
druglist <- druglist[!is.na(druglist$drug_index), ]
druglist <- druglist[druglist$proof != 'animal',]


# load( "07_plotting_v2/data/drugTarget_merged3DB.Rdata")
# drug_FDAcancer <- unique(drugTarget_FDA[which(drugTarget_FDA$ifCancer_drugbank == 'yes'),])
# # drug_FDAcancer <- data.frame(drugname = unique(DTFDA_drugbank$name), drugname_lower = tolower(unique(DTFDA_drugbank$name)))
# # drug_FDAcancer <- data.frame(drugname = unique(DTFDAcancer_drugbank$name), drugname_lower = tolower(unique(DTFDAcancer_drugbank$name)))



# drugMoA <- read.csv('/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/03_drug_immuneSig_enrichment/data/drug_moa_merged_cleaned.txt', sep = '\t', header = T)
# drugMoA <- read.csv('/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/s3_lincs_drug_MoA/lincs_portal2.0/drugMoA_lincsPortal2.csv',header = F)
# drugMoA <- do.call(rbind, 
#                   lapply(as.list(drugMoA$V2), function(x){
#                   tmp = drugMoA[drugMoA$V2 == x,]
#                   df = data.frame(drug_index = x, MoA = unique(unlist(tmp[-c(1,2)])))
#                   return(df) }))
# drugMoA <- unique(drugMoA[drugMoA$MoA != "", ])
# write.csv(drugMoA, '/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/s3_lincs_drug_MoA/lincs_portal2.0/drugMoA_lincsPortal2_2.csv', quote = F, row.names = F)
drugMoA <- read.csv('/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/s3_lincs_drug_MoA/lincs_portal2.0/drugMoA_lincsPortal2_2.csv',header = T)
drugTarget <- read.csv('/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/s3_lincs_drug_MoA/lincs_portal2.0/drugtarget_lincsPortal2.csv', header = T)
drug_maxfda <- read.csv('/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/s3_lincs_drug_MoA/lincs_portal2.0/drug_max_fda_phase_lincsPortal2.csv', header = T)
drug_maxfda34 <- unique(drug_maxfda[is.element(drug_maxfda$max_fda_phase,c(3,4)), ])

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# 70138  -----------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

filepath <- "06_results_summaryandplotting/results_meta_filter_ITSssgsea_CombatRef_17datasets/TUMERIC_m2_meta_en_ACATP0.05/70138/allgenes/"
resfilepath <- list.files(filepath)

filenames <- list.files(paste0(filepath, "SKCM_Metastatic/drugITSRank/"))
filenames <- filenames[grep('.txt', filenames)]
filenames <- filenames[grep('_weighted_', filenames)]

for(f in filenames){
  #  f=filenames[1]
  dir.create(paste0('07_plotting_v2/05_drugITSprofiler_CombatRef_metaPlusFilter/drugrank_diff_isgold_lincsPortal2_metaPlusFilter/70138/', gsub("_res.txt","", f)))
  savepath = paste0('07_plotting_v2/05_drugITSprofiler_CombatRef_metaPlusFilter/drugrank_diff_isgold_lincsPortal2_metaPlusFilter/70138/', gsub("_res.txt","", f),'/')
  
  for(i in 1:length(resfilepath)){
    # i = 8
    
    cancer = strsplit(resfilepath[i],'_')[[1]][1]
    print(cancer)
    
    res_weighted_en <- read.csv(paste0(filepath, resfilepath[i], "/drugITSRank/", f), header = T, sep = '\t')
    inst_select <- inst_GSE70138_is_gold[is.element(inst_GSE70138_is_gold$cell_id, cell_70138[cell_70138$TCGA_cancer_type == cancer,]$cell_id), ]
    res_weighted_en <- merge(res_weighted_en, inst_select, by = 'pert_iname')
    res_weighted_en <- unique(res_weighted_en[!is.element(names(res_weighted_en),names(inst_GSE70138_is_gold)[-which(names(inst_GSE70138_is_gold) == 'pert_iname')])])
    
    res_weighted_en$rank_tile <- res_weighted_en$rank / nrow(res_weighted_en)
    res_weighted_en$cancer <- resfilepath[i]
    res_weighted_en$dataset <- "70138"
    
    
    for(proofType in c(cancer, 'FDA', 'ClinicalTrial_III', 'ClinicalTrial_beforeIII','ClinicalTrial_afterIII', 'proof')){
      
      res_fda <- res_weighted_en # [is.element(res_weighted_en$pert_iname, drug_FDAcancer$drugname_lower), ]
      res_fda$label = 'others'
      
      if(is.element(proofType, c( 'FDA', 'ClinicalTrial_III'))){
        res_fda[is.element(res_fda$drug_index, unique(druglist[is.element(druglist$proof, proofType), ]$drug_index)), ]$label = 'fda_combo'
      }else if(proofType == 'ClinicalTrial_afterIII'){
        proof = c('FDA', 'ClinicalTrial_III')
        res_fda[is.element(res_fda$drug_index, unique(druglist[is.element(druglist$proof, proof), ]$drug_index)), ]$label = 'fda_combo'
      }else if(proofType == 'ClinicalTrial_beforeIII'){
        proof = c("1","1-2","1b","2","1/2a","1b/2","1/2")
        res_fda[is.element(res_fda$drug_index, unique(druglist[is.element(druglist$proof, proof), ]$drug_index)), ]$label = 'fda_combo'
      }else if(proofType == 'proof'){
        res_fda[is.element(res_fda$drug_index, unique(druglist$drug_index)), ]$label = 'fda_combo'
    }else if(is.element(proofType, unique(druglist$Cancer_TCGA))) {  
      if(length(intersect(res_fda$drug_index, unique(druglist[druglist$Cancer_TCGA == cancer,]$drug_index)))>0) {  
        res_fda[is.element(res_fda$drug_index, unique(druglist[druglist$Cancer_TCGA == cancer,]$drug_index)), ]$label = 'fda_combo'
      }
    }


      
      
      if(length(unique(res_fda$label))>1){
        
        res_fda <- res_fda[order(res_fda$score, decreasing = T), ]
        res_fda$rank2 = 1:nrow(res_fda)
        res_fda$rank2_tile = res_fda$rank2 / nrow(res_fda)
        
        pdf(paste0(savepath, cancer, "_", proofType,'.pdf'), height = 20, width = 10)
        par(mfrow = c(5,2))
        plot(density(res_fda[res_fda$label == 'others',]$score), main = 'all')
        if(nrow(res_fda[res_fda$label == 'fda_combo',])>1) lines(density(res_fda[res_fda$label == 'fda_combo',]$score), col = 'red')
        boxplot(score~label, data=res_fda, main = 'all')
        
        dt1 = data.frame(index = 'all',
                         num_drug = nrow(res_fda), 
                         num_proof = table(res_fda$label)[1],
                         num_others = table(res_fda$label)[2],
                         mean_proof = mean(res_fda[res_fda$label == 'fda_combo',]$score),
                         mean_others = mean(res_fda[res_fda$label == 'others',]$score),
                         effectsize = mean(res_fda[res_fda$label == 'fda_combo',]$score) - mean(res_fda[res_fda$label == 'others',]$score),
                         pval = wilcox.test(res_fda[res_fda$label == 'fda_combo',]$score, 
                                            res_fda[res_fda$label == 'others',]$score, 
                                            alternative = 'greater')$p.value)
        
        
        
        tmp1 = res_fda[res_fda$label == 'fda_combo',]
        tmp1 = unique(tmp1)
        dim(tmp1)
        res_fda1 = merge(res_fda, drugTarget[c('drug_index','target')], by = 'drug_index')
        
        sharetarget = intersect(unique(res_fda1[res_fda1$label == 'fda_combo',]$target), unique(res_fda1[res_fda1$label == 'others',]$target))
        dst = unique(res_fda1[is.element(res_fda1$target, sharetarget),]$drug_index)
        tmp2 = res_fda1[res_fda1$label == 'others',]
        tmp2 = tmp2[!is.element(tmp2$drug_index, dst), ]
        dim(tmp2)
        tmp2 = unique(tmp2[!is.element(names(tmp2),'target')])
        tmp2 = unique(tmp2)
        dim(tmp2)
        tmp2 = tmp2[!is.element(tmp2$drug_index, unique(druglist$drug_index)),]
        dim(tmp2)
        
        res_fda1 = merge(res_fda, drugMoA[c('drug_index','MoA')], by = 'drug_index')
        shareMoA = intersect(unique(res_fda1[res_fda1$label == 'fda_combo',]$MoA), unique(res_fda1[res_fda1$label == 'others',]$MoA))
        dst = unique(res_fda1[is.element(res_fda1$MoA, shareMoA),]$drug_index)
        # tmp2 = res_fda1[res_fda1$label == 'others',]
        tmp2 = tmp2[!is.element(tmp2$drug_index, dst), ]
        dim(tmp2)
        tmp2 = unique(tmp2[!is.element(names(tmp2),'MoA')])
        tmp2 = unique(tmp2)
        dim(tmp2)
        tmp2 = tmp2[!is.element(tmp2$drug_index, unique(druglist$drug_index)),]
        dim(tmp2)
        
        res1 = unique(rbind(tmp1, tmp2))
        res1 = res1[order(res1$score, decreasing = T),]
        res1$rank2 = 1:nrow(res1)
        res1$rank2_tile = res1$rank2 / nrow(res1)
        
        dt2 = data.frame(index = 'all_filterMoA',
                         num_drug = nrow(res1), 
                         num_proof = table(res1$label)[1],
                         num_others = table(res1$label)[2],
                         mean_proof = mean(res1[res1$label == 'fda_combo',]$score),
                         mean_others = mean(res1[res1$label == 'others',]$score),
                         effectsize = mean(res1[res1$label == 'fda_combo',]$score) - mean(res1[res1$label == 'others',]$score),
                         pval = wilcox.test(res1[res1$label == 'fda_combo',]$score, 
                                            res1[res1$label == 'others',]$score, 
                                            alternative = 'greater')$p.value)
        
        plot(density(res1[res1$label == 'others',]$score), main='all_filterMoA')
        if(nrow(res1[res1$label == 'fda_combo',])>1) lines(density(res1[res1$label == 'fda_combo',]$score), col = 'red')
        boxplot(score~label, data=res1, main='all_filterMoA')
        
        
        res_fda2 <- res_fda1[is.element(res_fda1$drug_index, unique(drug_maxfda34$drug_index)), ]
        # res_fda2$drugname_lower <- res_fda2$pert_iname
        res_fda2 <- res_fda2[order(res_fda2$score, decreasing = T), ]
        
        plot(density(res_fda2[res_fda2$label == 'others',]$score), main = 'FDAcancer')
        if(nrow(res_fda2[res_fda2$label == 'fda_combo',])>1) lines(density(res_fda2[res_fda2$label == 'fda_combo',]$score), col = 'red')
        boxplot(score~label, data=res_fda2,main = 'FDAcancer')
        
        dt3 = data.frame(index = 'proof',
                         num_drug = nrow(res_fda2), 
                         num_proof = table(res_fda2$label)[1],
                         num_others = table(res_fda2$label)[2],
                         mean_proof = mean(res_fda2[res_fda2$label == 'fda_combo',]$score),
                         mean_others = mean(res_fda2[res_fda2$label == 'others',]$score),
                         effectsize = mean(res_fda2[res_fda2$label == 'fda_combo',]$score) - mean(res_fda2[res_fda2$label == 'others',]$score),
                         pval = wilcox.test(res_fda2[res_fda2$label == 'fda_combo',]$score, 
                                            res_fda2[res_fda2$label == 'others',]$score, 
                                            alternative = 'greater')$p.value)
        
        
        # tmp1 = res_fda2[res_fda2$label == 'fda_combo',]
        # tmp1 = unique(tmp1)
        
        res_fda3 = merge(res_fda2, drugTarget[c('drug_index','target')], by = 'drug_index')
        
        sharetarget = intersect(unique(res_fda3[res_fda3$label == 'fda_combo',]$target), unique(res_fda3[res_fda3$label == 'others',]$target))
        dst = unique(res_fda3[is.element(res_fda3$target, sharetarget),]$drug_index)
        tmp2 = res_fda3[res_fda3$label == 'others',]
        tmp2 = tmp2[!is.element(tmp2$drug_index, dst), ]
        dim(tmp2)
        tmp2 = unique(tmp2[!is.element(names(tmp2),'target')])
        tmp2 = unique(tmp2)
        dim(tmp2)
        tmp2 = tmp2[!is.element(tmp2$drug_index, unique(druglist$drug_index)),]
        dim(tmp2)
        
        res_fda3 = merge(res_fda2, drugMoA[c('drug_index','MoA')], by = 'drug_index')
        shareMoA = intersect(unique(res_fda3[res_fda3$label == 'fda_combo',]$MoA), unique(res_fda3[res_fda3$label == 'others',]$MoA))
        dst = unique(res_fda3[is.element(res_fda3$MoA, shareMoA),]$drug_index)
        # tmp2 = res_fda3[res_fda3$label == 'others',]
        tmp2 = tmp2[!is.element(tmp2$drug_index, dst), ]
        dim(tmp2)
        tmp2 = unique(tmp2[!is.element(names(tmp2),'MoA')])
        tmp2 = unique(tmp2)
        dim(tmp2)
        tmp2 = tmp2[!is.element(tmp2$drug_index, unique(druglist$drug_index)),]
        dim(tmp2)
        
        res2 = unique(rbind(tmp1, tmp2))
        res2 = res2[order(res2$score, decreasing = T),]
        res2$rank2 = 1:nrow(res2)
        res2$rank2_tile = res2$rank2 / nrow(res2)
        
        dt4 = data.frame(index = 'proof_filterMoA',
                         num_drug = nrow(res2), 
                         num_proof = table(res2$label)[1],
                         num_others = table(res2$label)[2],
                         mean_proof = mean(res2[res2$label == 'fda_combo',]$score),
                         mean_others = mean(res2[res2$label == 'others',]$score),
                         effectsize = mean(res2[res2$label == 'fda_combo',]$score) - mean(res2[res2$label == 'others',]$score),
                         pval = wilcox.test(res2[res2$label == 'fda_combo',]$score, 
                                            res2[res2$label == 'others',]$score, 
                                            alternative = 'greater')$p.value)
        
        plot(density(res2[res2$label == 'others',]$score), main='FDAcancer_filterMoA')
        if(nrow(res2[res2$label == 'fda_combo',])>1) lines(density(res2[res2$label == 'fda_combo',]$score), col = 'red')
        boxplot(score~label, data=res2, main='FDAcancer_filterMoA')
        
        # boxplot(rank2_tile~label, data=res2)
        # plot(density(res2[res2$label == 'others',]$rank2_tile))
        # lines(density(res2[res2$label == 'fda_combo',]$rank2_tile), col = 'red')
        
        if(nrow(tmp2)>nrow(tmp1)){
          
          set.seed(1234)
          testp = lapply(as.list(1:1000), function(i){
            
            tmpindex = sample( 1:nrow(tmp2), nrow(tmp1), replace = FALSE, prob = NULL)
            tmp3 = tmp2[tmpindex,]
            res3 = rbind(tmp1, tmp3)
            res3 = res3[order(res3$score, decreasing = T),]
            res3$rank2 = 1:nrow(res3)
            res3$rank2_tile = res3$rank2 / nrow(res3)
            
            dt5 = data.frame(index = 'proof_filterMoA_sample',
                             num_drug = nrow(res3), 
                             num_proof = table(res3$label)[1],
                             num_others = table(res3$label)[2],
                             mean_proof = mean(res3[res3$label == 'fda_combo',]$score),
                             mean_others = mean(res3[res3$label == 'others',]$score),
                             effectsize = mean(res3[res3$label == 'fda_combo',]$score) - mean(res3[res3$label == 'others',]$score),
                             pval = wilcox.test(res3[res3$label == 'fda_combo',]$score, 
                                                res3[res3$label == 'others',]$score, 
                                                alternative = 'greater')$p.value)
            
            return(dt5)
            
          })
          testp = do.call(rbind, testp)
          write.csv(testp, paste0(savepath, cancer, "_", proofType,'_diffPermTest.csv'), row.names = F, quote = F)
          
          plot(density(testp$pval), main = paste0('p<0.05: ', nrow(testp[testp$pval < 0.05,])))
          abline(v=0.05)
          
        }
        
        dev.off()
        
        
        dt = rbind(dt1, dt2, dt3, dt4)
        write.csv(dt, paste0(savepath, cancer, "_", proofType,'_diffTest.csv'), row.names = F, quote = F)
        
      }
    }
  }
  
  res70138_SKCMprimary <-  read.csv(paste0("06_results_summaryandplotting/results_meta_filter_ITSssgsea_CombatRef_17datasets/CPE_m2_meta_en_ACATP0.05/70138/allgenes/SKCM_Primary/drugITSRank/", f), header = T, sep = '\t')
  inst_select <- inst_GSE70138_is_gold[is.element(inst_GSE70138_is_gold$cell_id, cell_70138[cell_70138$TCGA_cancer_type == cancer,]$cell_id), ]
  res70138_SKCMprimary <- merge(res70138_SKCMprimary, inst_select, by = 'pert_iname')
  res70138_SKCMprimary <- unique(res70138_SKCMprimary[!is.element(names(res70138_SKCMprimary),names(inst_GSE70138_is_gold)[-which(names(inst_GSE70138_is_gold) == 'pert_iname')])])

  res70138_SKCMprimary$rank_tile <- res70138_SKCMprimary$rank / nrow(res70138_SKCMprimary)
  res70138_SKCMprimary$cancer <- "SKCM_Primary"
  res70138_SKCMprimary$dataset = '70138'
  
  res_fda <- res70138_SKCMprimary
  cancer = "SKCM"

  for(proofType in c(cancer, 'FDA', 'ClinicalTrial_III', 'ClinicalTrial_beforeIII','ClinicalTrial_afterIII', 'proof')){
    
    res_fda$label = 'others'
    
    if(is.element(proofType, c( 'FDA', 'ClinicalTrial_III'))){
      res_fda[is.element(res_fda$drug_index, unique(druglist[is.element(druglist$proof, proofType), ]$drug_index)), ]$label = 'fda_combo'
    }else if(proofType == 'ClinicalTrial_afterIII'){
      proof = c('FDA', 'ClinicalTrial_III')
      res_fda[is.element(res_fda$drug_index, unique(druglist[is.element(druglist$proof, proof), ]$drug_index)), ]$label = 'fda_combo'
    }else if(proofType == 'ClinicalTrial_beforeIII'){
      proof = c("1","1-2","1b","2","1/2a","1b/2","1/2")
      res_fda[is.element(res_fda$drug_index, unique(druglist[is.element(druglist$proof, proof), ]$drug_index)), ]$label = 'fda_combo'
    }else if(proofType == 'proof'){
      res_fda[is.element(res_fda$drug_index, unique(druglist$drug_index)), ]$label = 'fda_combo'
    }else if(is.element(proofType, unique(druglist$Cancer_TCGA))) {  
      if(length(intersect(res_fda$drug_index, unique(druglist[druglist$Cancer_TCGA == cancer,]$drug_index)))>0) {  
        res_fda[is.element(res_fda$drug_index, unique(druglist[druglist$Cancer_TCGA == cancer,]$drug_index)), ]$label = 'fda_combo'
      }
    }


      
      if(length(unique(res_fda$label))>1){
        
        res_fda <- res_fda[order(res_fda$score, decreasing = T), ]
        res_fda$rank2 = 1:nrow(res_fda)
        res_fda$rank2_tile = res_fda$rank2 / nrow(res_fda)
        
        pdf(paste0(savepath, cancer, "_", proofType,'.pdf'), height = 20, width = 10)
        par(mfrow = c(5,2))
        plot(density(res_fda[res_fda$label == 'others',]$score), main = 'all')
        if(nrow(res_fda[res_fda$label == 'fda_combo',])>1) lines(density(res_fda[res_fda$label == 'fda_combo',]$score), col = 'red')
        boxplot(score~label, data=res_fda, main = 'all')
        
        dt1 = data.frame(index = 'all',
                         num_drug = nrow(res_fda), 
                         num_proof = table(res_fda$label)[1],
                         num_others = table(res_fda$label)[2],
                         mean_proof = mean(res_fda[res_fda$label == 'fda_combo',]$score),
                         mean_others = mean(res_fda[res_fda$label == 'others',]$score),
                         effectsize = mean(res_fda[res_fda$label == 'fda_combo',]$score) - mean(res_fda[res_fda$label == 'others',]$score),
                         pval = wilcox.test(res_fda[res_fda$label == 'fda_combo',]$score, 
                                            res_fda[res_fda$label == 'others',]$score, 
                                            alternative = 'greater')$p.value)
        
        
        
        tmp1 = res_fda[res_fda$label == 'fda_combo',]
        tmp1 = unique(tmp1)
        dim(tmp1)
        res_fda1 = merge(res_fda, drugTarget[c('drug_index','target')], by = 'drug_index')
        
        sharetarget = intersect(unique(res_fda1[res_fda1$label == 'fda_combo',]$target), unique(res_fda1[res_fda1$label == 'others',]$target))
        dst = unique(res_fda1[is.element(res_fda1$target, sharetarget),]$drug_index)
        tmp2 = res_fda1[res_fda1$label == 'others',]
        tmp2 = tmp2[!is.element(tmp2$drug_index, dst), ]
        dim(tmp2)
        tmp2 = unique(tmp2[!is.element(names(tmp2),'target')])
        tmp2 = unique(tmp2)
        dim(tmp2)
        tmp2 = tmp2[!is.element(tmp2$drug_index, unique(druglist$drug_index)),]
        dim(tmp2)
        
        res_fda1 = merge(res_fda, drugMoA[c('drug_index','MoA')], by = 'drug_index')
        shareMoA = intersect(unique(res_fda1[res_fda1$label == 'fda_combo',]$MoA), unique(res_fda1[res_fda1$label == 'others',]$MoA))
        dst = unique(res_fda1[is.element(res_fda1$MoA, shareMoA),]$drug_index)
        # tmp2 = res_fda1[res_fda1$label == 'others',]
        tmp2 = tmp2[!is.element(tmp2$drug_index, dst), ]
        dim(tmp2)
        tmp2 = unique(tmp2[!is.element(names(tmp2),'MoA')])
        tmp2 = unique(tmp2)
        dim(tmp2)
        tmp2 = tmp2[!is.element(tmp2$drug_index, unique(druglist$drug_index)),]
        dim(tmp2)
        
        res1 = unique(rbind(tmp1, tmp2))
        res1 = res1[order(res1$score, decreasing = T),]
        res1$rank2 = 1:nrow(res1)
        res1$rank2_tile = res1$rank2 / nrow(res1)
        
        dt2 = data.frame(index = 'all_filterMoA',
                         num_drug = nrow(res1), 
                         num_proof = table(res1$label)[1],
                         num_others = table(res1$label)[2],
                         mean_proof = mean(res1[res1$label == 'fda_combo',]$score),
                         mean_others = mean(res1[res1$label == 'others',]$score),
                         effectsize = mean(res1[res1$label == 'fda_combo',]$score) - mean(res1[res1$label == 'others',]$score),
                         pval = wilcox.test(res1[res1$label == 'fda_combo',]$score, 
                                            res1[res1$label == 'others',]$score, 
                                            alternative = 'greater')$p.value)
        
        plot(density(res1[res1$label == 'others',]$score), main='all_filterMoA')
        if(nrow(res1[res1$label == 'fda_combo',])>1) lines(density(res1[res1$label == 'fda_combo',]$score), col = 'red')
        boxplot(score~label, data=res1, main='all_filterMoA')
        
        
        res_fda2 <- res_fda1[is.element(res_fda1$drug_index, unique(drug_maxfda34$drug_index)), ]
        # res_fda2$drugname_lower <- res_fda2$pert_iname
        res_fda2 <- res_fda2[order(res_fda2$score, decreasing = T), ]
        
        plot(density(res_fda2[res_fda2$label == 'others',]$score), main = 'FDAcancer')
        if(nrow(res_fda2[res_fda2$label == 'fda_combo',])>1) lines(density(res_fda2[res_fda2$label == 'fda_combo',]$score), col = 'red')
        boxplot(score~label, data=res_fda2,main = 'FDAcancer')
        
        dt3 = data.frame(index = 'proof',
                         num_drug = nrow(res_fda2), 
                         num_proof = table(res_fda2$label)[1],
                         num_others = table(res_fda2$label)[2],
                         mean_proof = mean(res_fda2[res_fda2$label == 'fda_combo',]$score),
                         mean_others = mean(res_fda2[res_fda2$label == 'others',]$score),
                         effectsize = mean(res_fda2[res_fda2$label == 'fda_combo',]$score) - mean(res_fda2[res_fda2$label == 'others',]$score),
                         pval = wilcox.test(res_fda2[res_fda2$label == 'fda_combo',]$score, 
                                            res_fda2[res_fda2$label == 'others',]$score, 
                                            alternative = 'greater')$p.value)
        
        
        # tmp1 = res_fda2[res_fda2$label == 'fda_combo',]
        # tmp1 = unique(tmp1)
        
        res_fda3 = merge(res_fda2, drugTarget[c('drug_index','target')], by = 'drug_index')
        
        sharetarget = intersect(unique(res_fda3[res_fda3$label == 'fda_combo',]$target), unique(res_fda3[res_fda3$label == 'others',]$target))
        dst = unique(res_fda3[is.element(res_fda3$target, sharetarget),]$drug_index)
        tmp2 = res_fda3[res_fda3$label == 'others',]
        tmp2 = tmp2[!is.element(tmp2$drug_index, dst), ]
        dim(tmp2)
        tmp2 = unique(tmp2[!is.element(names(tmp2),'target')])
        tmp2 = unique(tmp2)
        dim(tmp2)
        tmp2 = tmp2[!is.element(tmp2$drug_index, unique(druglist$drug_index)),]
        dim(tmp2)
        
        res_fda3 = merge(res_fda2, drugMoA[c('drug_index','MoA')], by = 'drug_index')
        shareMoA = intersect(unique(res_fda3[res_fda3$label == 'fda_combo',]$MoA), unique(res_fda3[res_fda3$label == 'others',]$MoA))
        dst = unique(res_fda3[is.element(res_fda3$MoA, shareMoA),]$drug_index)
        # tmp2 = res_fda3[res_fda3$label == 'others',]
        tmp2 = tmp2[!is.element(tmp2$drug_index, dst), ]
        dim(tmp2)
        tmp2 = unique(tmp2[!is.element(names(tmp2),'MoA')])
        tmp2 = unique(tmp2)
        dim(tmp2)
        tmp2 = tmp2[!is.element(tmp2$drug_index, unique(druglist$drug_index)),]
        dim(tmp2)
        
        res2 = unique(rbind(tmp1, tmp2))
        res2 = res2[order(res2$score, decreasing = T),]
        res2$rank2 = 1:nrow(res2)
        res2$rank2_tile = res2$rank2 / nrow(res2)
        
        dt4 = data.frame(index = 'proof_filterMoA',
                         num_drug = nrow(res2), 
                         num_proof = table(res2$label)[1],
                         num_others = table(res2$label)[2],
                         mean_proof = mean(res2[res2$label == 'fda_combo',]$score),
                         mean_others = mean(res2[res2$label == 'others',]$score),
                         effectsize = mean(res2[res2$label == 'fda_combo',]$score) - mean(res2[res2$label == 'others',]$score),
                         pval = wilcox.test(res2[res2$label == 'fda_combo',]$score, 
                                            res2[res2$label == 'others',]$score, 
                                            alternative = 'greater')$p.value)
        
        plot(density(res2[res2$label == 'others',]$score), main='FDAcancer_filterMoA')
        if(nrow(res2[res2$label == 'fda_combo',])>1) lines(density(res2[res2$label == 'fda_combo',]$score), col = 'red')
        boxplot(score~label, data=res2, main='FDAcancer_filterMoA')
        
        # boxplot(rank2_tile~label, data=res2)
        # plot(density(res2[res2$label == 'others',]$rank2_tile))
        # lines(density(res2[res2$label == 'fda_combo',]$rank2_tile), col = 'red')
        
        if(nrow(tmp2)>nrow(tmp1)){
          
          set.seed(1234)
          testp = lapply(as.list(1:1000), function(i){
            
            tmpindex = sample( 1:nrow(tmp2), nrow(tmp1), replace = FALSE, prob = NULL)
            tmp3 = tmp2[tmpindex,]
            res3 = rbind(tmp1, tmp3)
            res3 = res3[order(res3$score, decreasing = T),]
            res3$rank2 = 1:nrow(res3)
            res3$rank2_tile = res3$rank2 / nrow(res3)
            
            dt5 = data.frame(index = 'proof_filterMoA_sample',
                             num_drug = nrow(res3), 
                             num_proof = table(res3$label)[1],
                             num_others = table(res3$label)[2],
                             mean_proof = mean(res3[res3$label == 'fda_combo',]$score),
                             mean_others = mean(res3[res3$label == 'others',]$score),
                             effectsize = mean(res3[res3$label == 'fda_combo',]$score) - mean(res3[res3$label == 'others',]$score),
                             pval = wilcox.test(res3[res3$label == 'fda_combo',]$score, 
                                                res3[res3$label == 'others',]$score, 
                                                alternative = 'greater')$p.value)
            
            return(dt5)
            
          })
          testp = do.call(rbind, testp)
          write.csv(testp, paste0(savepath, "SKCM_Primary_", proofType,'_diffPermTest.csv'), row.names = F, quote = F)
          
          plot(density(testp$pval), main = paste0('p<0.05: ', nrow(testp[testp$pval < 0.05,])))
          abline(v=0.05)
          
        }
        
        dev.off()
        
        
        dt = rbind(dt1, dt2, dt3, dt4)
        write.csv(dt, paste0(savepath, "SKCM_Primary_", proofType,'_diffTest.csv'), row.names = F, quote = F)
        
    }
  }
}

  
  # ------------------------------------------------------------------------------
  # ------------------------------------------------------------------------------
  # 92742  -----------------------------------------------------------------------
  # ------------------------------------------------------------------------------
  # ------------------------------------------------------------------------------
  
  filepath <- "06_results_summaryandplotting/results_meta_filter_ITSssgsea_CombatRef_17datasets/TUMERIC_m2_meta_en_ACATP0.05/92742/allgenes/"
  resfilepath <- list.files(filepath)
  
  filenames <- list.files(paste0(filepath, "SKCM_Metastatic/drugITSRank/"))
  filenames <- filenames[grep('.txt', filenames)]
filenames <- filenames[grep('_weighted_', filenames)]
  
 dir.create('07_plotting_v2/05_drugITSprofiler_CombatRef_metaPlusFilter/drugrank_diff_isgold_lincsPortal2_metaPlusFilter/92742/')
 
for(f in filenames){

  dir.create(paste0('07_plotting_v2/05_drugITSprofiler_CombatRef_metaPlusFilter/drugrank_diff_isgold_lincsPortal2_metaPlusFilter/92742/', gsub("_res.txt","", f)))
  savepath = paste0('07_plotting_v2/05_drugITSprofiler_CombatRef_metaPlusFilter/drugrank_diff_isgold_lincsPortal2_metaPlusFilter/92742/', gsub("_res.txt","", f),'/')
  
  for(i in 1:length(resfilepath)){
    # i = 11
    
    cancer = strsplit(resfilepath[i],'_')[[1]][1]
    print(cancer)
    res_weighted_en <- read.csv(paste0(filepath, resfilepath[i], "/drugITSRank/", f), header = T, sep = '\t')
    inst_select <- inst_GSE92742_is_gold[is.element(inst_GSE92742_is_gold$cell_id, cell_92742[cell_92742$TCGA_cancer_type == cancer,]$cell_id), ]
    res_weighted_en <- merge(res_weighted_en, inst_select, by = 'pert_iname')
    res_weighted_en <- unique(res_weighted_en[!is.element(names(res_weighted_en),names(inst_GSE92742_is_gold)[-which(names(inst_GSE92742_is_gold) == 'pert_iname')])])
    
    res_weighted_en$rank_tile <- res_weighted_en$rank / nrow(res_weighted_en)
    res_weighted_en$cancer <- resfilepath[i]
    res_weighted_en$dataset <- "92742"
    
    
    for(proofType in c(cancer, 'FDA', 'ClinicalTrial_III', 'ClinicalTrial_beforeIII','ClinicalTrial_afterIII', 'proof')){
      
      res_fda <- res_weighted_en # [is.element(res_weighted_en$pert_iname, drug_FDAcancer$drugname_lower), ]
      res_fda$label = 'others'
      
      if(is.element(proofType, c( 'FDA', 'ClinicalTrial_III'))){
        res_fda[is.element(res_fda$drug_index, unique(druglist[is.element(druglist$proof, proofType), ]$drug_index)), ]$label = 'fda_combo'
      }else if(proofType == 'ClinicalTrial_afterIII'){
        proof = c('FDA', 'ClinicalTrial_III')
        res_fda[is.element(res_fda$drug_index, unique(druglist[is.element(druglist$proof, proof), ]$drug_index)), ]$label = 'fda_combo'
      }else if(proofType == 'ClinicalTrial_beforeIII'){
        proof = c("1","1-2","1b","2","1/2a","1b/2","1/2")
        res_fda[is.element(res_fda$drug_index, unique(druglist[is.element(druglist$proof, proof), ]$drug_index)), ]$label = 'fda_combo'
      }else if(proofType == 'proof'){
        res_fda[is.element(res_fda$drug_index, unique(druglist$drug_index)), ]$label = 'fda_combo'
      }else if(is.element(proofType, unique(druglist$Cancer_TCGA))) {  
        if(length(intersect(res_fda$drug_index, unique(druglist[druglist$Cancer_TCGA == cancer,]$drug_index)))>0) {  
        res_fda[is.element(res_fda$drug_index, unique(druglist[druglist$Cancer_TCGA == cancer,]$drug_index)), ]$label = 'fda_combo'
        }
      }
      
      
      if(length(unique(res_fda$label))>1){
        
        res_fda <- res_fda[order(res_fda$score, decreasing = T), ]
        res_fda$rank2 = 1:nrow(res_fda)
        res_fda$rank2_tile = res_fda$rank2 / nrow(res_fda)
        
        pdf(paste0(savepath, cancer, "_", proofType,'.pdf'), height = 20, width = 10)
        par(mfrow = c(5,2))
        plot(density(res_fda[res_fda$label == 'others',]$score), main = 'all')
        if(nrow(res_fda[res_fda$label == 'fda_combo',])>1) lines(density(res_fda[res_fda$label == 'fda_combo',]$score), col = 'red')
        boxplot(score~label, data=res_fda, main = 'all')
        
        dt1 = data.frame(index = 'all',
                         num_drug = nrow(res_fda), 
                         num_proof = table(res_fda$label)[1],
                         num_others = table(res_fda$label)[2],
                         mean_proof = mean(res_fda[res_fda$label == 'fda_combo',]$score),
                         mean_others = mean(res_fda[res_fda$label == 'others',]$score),
                         effectsize = mean(res_fda[res_fda$label == 'fda_combo',]$score) - mean(res_fda[res_fda$label == 'others',]$score),
                         pval = wilcox.test(res_fda[res_fda$label == 'fda_combo',]$score, 
                                            res_fda[res_fda$label == 'others',]$score, 
                                            alternative = 'greater')$p.value)
        
        
        
        tmp1 = res_fda[res_fda$label == 'fda_combo',]
        tmp1 = unique(tmp1)
        dim(tmp1)
        res_fda1 = merge(res_fda, drugTarget[c('drug_index','target')], by = 'drug_index')
        
        sharetarget = intersect(unique(res_fda1[res_fda1$label == 'fda_combo',]$target), unique(res_fda1[res_fda1$label == 'others',]$target))
        dst = unique(res_fda1[is.element(res_fda1$target, sharetarget),]$drug_index)
        tmp2 = res_fda1[res_fda1$label == 'others',]
        tmp2 = tmp2[!is.element(tmp2$drug_index, dst), ]
        dim(tmp2)
        tmp2 = unique(tmp2[!is.element(names(tmp2),'target')])
        tmp2 = unique(tmp2)
        dim(tmp2)
        tmp2 = tmp2[!is.element(tmp2$drug_index, unique(druglist$drug_index)),]
        dim(tmp2)
        
        res_fda1 = merge(res_fda, drugMoA[c('drug_index','MoA')], by = 'drug_index')
        shareMoA = intersect(unique(res_fda1[res_fda1$label == 'fda_combo',]$MoA), unique(res_fda1[res_fda1$label == 'others',]$MoA))
        dst = unique(res_fda1[is.element(res_fda1$MoA, shareMoA),]$drug_index)
        # tmp2 = res_fda1[res_fda1$label == 'others',]
        tmp2 = tmp2[!is.element(tmp2$drug_index, dst), ]
        dim(tmp2)
        tmp2 = unique(tmp2[!is.element(names(tmp2),'MoA')])
        tmp2 = unique(tmp2)
        dim(tmp2)
        tmp2 = tmp2[!is.element(tmp2$drug_index, unique(druglist$drug_index)),]
        dim(tmp2)
        
        res1 = unique(rbind(tmp1, tmp2))
        res1 = res1[order(res1$score, decreasing = T),]
        res1$rank2 = 1:nrow(res1)
        res1$rank2_tile = res1$rank2 / nrow(res1)
        
        dt2 = data.frame(index = 'all_filterMoA',
                         num_drug = nrow(res1), 
                         num_proof = table(res1$label)[1],
                         num_others = table(res1$label)[2],
                         mean_proof = mean(res1[res1$label == 'fda_combo',]$score),
                         mean_others = mean(res1[res1$label == 'others',]$score),
                         effectsize = mean(res1[res1$label == 'fda_combo',]$score) - mean(res1[res1$label == 'others',]$score),
                         pval = wilcox.test(res1[res1$label == 'fda_combo',]$score, 
                                            res1[res1$label == 'others',]$score, 
                                            alternative = 'greater')$p.value)
        
        plot(density(res1[res1$label == 'others',]$score), main='all_filterMoA')
        if(nrow(res1[res1$label == 'fda_combo',])>1 ) lines(density(res1[res1$label == 'fda_combo',]$score), col = 'red')
        boxplot(score~label, data=res1, main='all_filterMoA')
        
        
        res_fda2 <- res_fda1[is.element(res_fda1$drug_index, unique(drug_maxfda34$drug_index)), ]
        # res_fda2$drugname_lower <- res_fda2$pert_iname
        res_fda2 <- res_fda2[order(res_fda2$score, decreasing = T), ]
        
        plot(density(res_fda2[res_fda2$label == 'others',]$score), main = 'FDAcancer')
        if(nrow(res_fda2[res_fda2$label == 'fda_combo',])>1) lines(density(res_fda2[res_fda2$label == 'fda_combo',]$score), col = 'red')
        boxplot(score~label, data=res_fda2,main = 'FDAcancer')
        
        dt3 = data.frame(index = 'proof',
                         num_drug = nrow(res_fda2), 
                         num_proof = table(res_fda2$label)[1],
                         num_others = table(res_fda2$label)[2],
                         mean_proof = mean(res_fda2[res_fda2$label == 'fda_combo',]$score),
                         mean_others = mean(res_fda2[res_fda2$label == 'others',]$score),
                         effectsize = mean(res_fda2[res_fda2$label == 'fda_combo',]$score) - mean(res_fda2[res_fda2$label == 'others',]$score),
                         pval = wilcox.test(res_fda2[res_fda2$label == 'fda_combo',]$score, 
                                            res_fda2[res_fda2$label == 'others',]$score, 
                                            alternative = 'greater')$p.value)
        
        
        # tmp1 = res_fda2[res_fda2$label == 'fda_combo',]
        # tmp1 = unique(tmp1)
        
        res_fda3 = merge(res_fda2, drugTarget[c('drug_index','target')], by = 'drug_index')
        
        sharetarget = intersect(unique(res_fda3[res_fda3$label == 'fda_combo',]$target), unique(res_fda3[res_fda3$label == 'others',]$target))
        dst = unique(res_fda3[is.element(res_fda3$target, sharetarget),]$drug_index)
        tmp2 = res_fda3[res_fda3$label == 'others',]
        tmp2 = tmp2[!is.element(tmp2$drug_index, dst), ]
        dim(tmp2)
        tmp2 = unique(tmp2[!is.element(names(tmp2),'target')])
        tmp2 = unique(tmp2)
        dim(tmp2)
        tmp2 = tmp2[!is.element(tmp2$drug_index, unique(druglist$drug_index)),]
        dim(tmp2)
        
        res_fda3 = merge(res_fda2, drugMoA[c('drug_index','MoA')], by = 'drug_index')
        shareMoA = intersect(unique(res_fda3[res_fda3$label == 'fda_combo',]$MoA), unique(res_fda3[res_fda3$label == 'others',]$MoA))
        dst = unique(res_fda3[is.element(res_fda3$MoA, shareMoA),]$drug_index)
        # tmp2 = res_fda3[res_fda3$label == 'others',]
        tmp2 = tmp2[!is.element(tmp2$drug_index, dst), ]
        dim(tmp2)
        tmp2 = unique(tmp2[!is.element(names(tmp2),'MoA')])
        tmp2 = unique(tmp2)
        dim(tmp2)
        tmp2 = tmp2[!is.element(tmp2$drug_index, unique(druglist$drug_index)),]
        dim(tmp2)
        
        res2 = unique(rbind(tmp1, tmp2))
        res2 = res2[order(res2$score, decreasing = T),]
        res2$rank2 = 1:nrow(res2)
        res2$rank2_tile = res2$rank2 / nrow(res2)
        
        dt4 = data.frame(index = 'proof_filterMoA',
                         num_drug = nrow(res2), 
                         num_proof = table(res2$label)[1],
                         num_others = table(res2$label)[2],
                         mean_proof = mean(res2[res2$label == 'fda_combo',]$score),
                         mean_others = mean(res2[res2$label == 'others',]$score),
                         effectsize = mean(res2[res2$label == 'fda_combo',]$score) - mean(res2[res2$label == 'others',]$score),
                         pval = wilcox.test(res2[res2$label == 'fda_combo',]$score, 
                                            res2[res2$label == 'others',]$score, 
                                            alternative = 'greater')$p.value)
        
        plot(density(res2[res2$label == 'others',]$score), main='FDAcancer_filterMoA')
        if(nrow(res2[res2$label == 'fda_combo',])>1) lines(density(res2[res2$label == 'fda_combo',]$score), col = 'red')
        boxplot(score~label, data=res2, main='FDAcancer_filterMoA')
        
        # boxplot(rank2_tile~label, data=res2)
        # plot(density(res2[res2$label == 'others',]$rank2_tile))
        # lines(density(res2[res2$label == 'fda_combo',]$rank2_tile), col = 'red')
        
        if(nrow(tmp2)>nrow(tmp1)){
          
          set.seed(1234)
          testp = lapply(as.list(1:1000), function(i){
            
            tmpindex = sample( 1:nrow(tmp2), nrow(tmp1), replace = FALSE, prob = NULL)
            tmp3 = tmp2[tmpindex,]
            res3 = rbind(tmp1, tmp3)
            res3 = res3[order(res3$score, decreasing = T),]
            res3$rank2 = 1:nrow(res3)
            res3$rank2_tile = res3$rank2 / nrow(res3)
            
            dt5 = data.frame(index = 'proof_filterMoA_sample',
                             num_drug = nrow(res3), 
                             num_proof = table(res3$label)[1],
                             num_others = table(res3$label)[2],
                             mean_proof = mean(res3[res3$label == 'fda_combo',]$score),
                             mean_others = mean(res3[res3$label == 'others',]$score),
                             effectsize = mean(res3[res3$label == 'fda_combo',]$score) - mean(res3[res3$label == 'others',]$score),
                             pval = wilcox.test(res3[res3$label == 'fda_combo',]$score, 
                                                res3[res3$label == 'others',]$score, 
                                                alternative = 'greater')$p.value)
            
            return(dt5)
            
          })
          testp = do.call(rbind, testp)
          write.csv(testp, paste0(savepath, cancer, "_", proofType,'_diffPermTest.csv'), row.names = F, quote = F)
          
          plot(density(testp$pval), main = paste0('p<0.05: ', nrow(testp[testp$pval < 0.05,])))
          abline(v=0.05)
          
        }
        
        dev.off()
        
        
        dt = rbind(dt1, dt2, dt3, dt4)
        write.csv(dt, paste0(savepath, cancer, "_", proofType,'_diffTest.csv'), row.names = F, quote = F)
        
      }
    }
  }
  
  
  
  
  res92742_SKCMprimary <-  read.csv(paste0("06_results_summaryandplotting/results_meta_filter_ITSssgsea_CombatRef_17datasets/CPE_m2_meta_en_ACATP0.05/92742/allgenes/SKCM_Primary/drugITSRank/", f), header = T, sep = '\t')
  inst_select <- inst_GSE92742_is_gold[is.element(inst_GSE92742_is_gold$cell_id, cell_92742[cell_92742$TCGA_cancer_type == cancer,]$cell_id), ]
  res92742_SKCMprimary <- merge(res92742_SKCMprimary, inst_select, by = 'pert_iname')
  res92742_SKCMprimary <- unique(res92742_SKCMprimary[!is.element(names(res92742_SKCMprimary),names(inst_GSE92742_is_gold)[-which(names(inst_GSE92742_is_gold) == 'pert_iname')])])
  
  res92742_SKCMprimary$rank_tile <- res92742_SKCMprimary$rank / nrow(res92742_SKCMprimary)
  res92742_SKCMprimary$cancer <- "SKCM_Primary"
  res92742_SKCMprimary$dataset = '92742'
  cancer = 'SKCM'
  
  for(proofType in c(cancer, 'FDA', 'ClinicalTrial_III', 'ClinicalTrial_beforeIII','ClinicalTrial_afterIII', 'proof')){
    
    res_fda <- res92742_SKCMprimary # [is.element(res_weighted_en$pert_iname, drug_FDAcancer$drugname_lower), ]
    res_fda$label = 'others'
    
    if(is.element(proofType, c( 'FDA', 'ClinicalTrial_III'))){
      res_fda[is.element(res_fda$drug_index, unique(druglist[is.element(druglist$proof, proofType), ]$drug_index)), ]$label = 'fda_combo'
    }else if(proofType == 'ClinicalTrial_afterIII'){
      proof = c('FDA', 'ClinicalTrial_III')
      res_fda[is.element(res_fda$drug_index, unique(druglist[is.element(druglist$proof, proof), ]$drug_index)), ]$label = 'fda_combo'
    }else if(proofType == 'ClinicalTrial_beforeIII'){
      proof = c("1","1-2","1b","2","1/2a","1b/2","1/2")
      res_fda[is.element(res_fda$drug_index, unique(druglist[is.element(druglist$proof, proof), ]$drug_index)), ]$label = 'fda_combo'
    }else if(proofType == 'proof'){
      res_fda[is.element(res_fda$drug_index, unique(druglist$drug_index)), ]$label = 'fda_combo'
    }else if(is.element(proofType, unique(druglist$Cancer_TCGA))) {  
        if(length(intersect(res_fda$drug_index, unique(druglist[druglist$Cancer_TCGA == cancer,]$drug_index)))>0) {  
        res_fda[is.element(res_fda$drug_index, unique(druglist[druglist$Cancer_TCGA == cancer,]$drug_index)), ]$label = 'fda_combo'
        }
      }
    
      
      if(length(unique(res_fda$label))>1){
        
        res_fda <- res_fda[order(res_fda$score, decreasing = T), ]
        res_fda$rank2 = 1:nrow(res_fda)
        res_fda$rank2_tile = res_fda$rank2 / nrow(res_fda)
        
        pdf(paste0(savepath, cancer, "_", proofType,'.pdf'), height = 20, width = 10)
        par(mfrow = c(5,2))
        plot(density(res_fda[res_fda$label == 'others',]$score), main = 'all')
        if(nrow(res_fda[res_fda$label == 'fda_combo',])>1) lines(density(res_fda[res_fda$label == 'fda_combo',]$score), col = 'red')
        boxplot(score~label, data=res_fda, main = 'all')
        
        dt1 = data.frame(index = 'all',
                         num_drug = nrow(res_fda), 
                         num_proof = table(res_fda$label)[1],
                         num_others = table(res_fda$label)[2],
                         mean_proof = mean(res_fda[res_fda$label == 'fda_combo',]$score),
                         mean_others = mean(res_fda[res_fda$label == 'others',]$score),
                         effectsize = mean(res_fda[res_fda$label == 'fda_combo',]$score) - mean(res_fda[res_fda$label == 'others',]$score),
                         pval = wilcox.test(res_fda[res_fda$label == 'fda_combo',]$score, 
                                            res_fda[res_fda$label == 'others',]$score, 
                                            alternative = 'greater')$p.value)
        
        
        
        tmp1 = res_fda[res_fda$label == 'fda_combo',]
        tmp1 = unique(tmp1)
        dim(tmp1)
        res_fda1 = merge(res_fda, drugTarget[c('drug_index','target')], by = 'drug_index')
        
        sharetarget = intersect(unique(res_fda1[res_fda1$label == 'fda_combo',]$target), unique(res_fda1[res_fda1$label == 'others',]$target))
        dst = unique(res_fda1[is.element(res_fda1$target, sharetarget),]$drug_index)
        tmp2 = res_fda1[res_fda1$label == 'others',]
        tmp2 = tmp2[!is.element(tmp2$drug_index, dst), ]
        dim(tmp2)
        tmp2 = unique(tmp2[!is.element(names(tmp2),'target')])
        tmp2 = unique(tmp2)
        dim(tmp2)
        tmp2 = tmp2[!is.element(tmp2$drug_index, unique(druglist$drug_index)),]
        dim(tmp2)
        
        res_fda1 = merge(res_fda, drugMoA[c('drug_index','MoA')], by = 'drug_index')
        shareMoA = intersect(unique(res_fda1[res_fda1$label == 'fda_combo',]$MoA), unique(res_fda1[res_fda1$label == 'others',]$MoA))
        dst = unique(res_fda1[is.element(res_fda1$MoA, shareMoA),]$drug_index)
        # tmp2 = res_fda1[res_fda1$label == 'others',]
        tmp2 = tmp2[!is.element(tmp2$drug_index, dst), ]
        dim(tmp2)
        tmp2 = unique(tmp2[!is.element(names(tmp2),'MoA')])
        tmp2 = unique(tmp2)
        dim(tmp2)
        tmp2 = tmp2[!is.element(tmp2$drug_index, unique(druglist$drug_index)),]
        dim(tmp2)
        
        res1 = unique(rbind(tmp1, tmp2))
        res1 = res1[order(res1$score, decreasing = T),]
        res1$rank2 = 1:nrow(res1)
        res1$rank2_tile = res1$rank2 / nrow(res1)
        
        dt2 = data.frame(index = 'all_filterMoA',
                         num_drug = nrow(res1), 
                         num_proof = table(res1$label)[1],
                         num_others = table(res1$label)[2],
                         mean_proof = mean(res1[res1$label == 'fda_combo',]$score),
                         mean_others = mean(res1[res1$label == 'others',]$score),
                         effectsize = mean(res1[res1$label == 'fda_combo',]$score) - mean(res1[res1$label == 'others',]$score),
                         pval = wilcox.test(res1[res1$label == 'fda_combo',]$score, 
                                            res1[res1$label == 'others',]$score, 
                                            alternative = 'greater')$p.value)
        
        plot(density(res1[res1$label == 'others',]$score), main='all_filterMoA')
        if(nrow(res1[res1$label == 'fda_combo',])>1) lines(density(res1[res1$label == 'fda_combo',]$score), col = 'red')
        boxplot(score~label, data=res1, main='all_filterMoA')
        
        
        res_fda2 <- res_fda1[is.element(res_fda1$drug_index, unique(drug_maxfda34$drug_index)), ]
        # res_fda2$drugname_lower <- res_fda2$pert_iname
        res_fda2 <- res_fda2[order(res_fda2$score, decreasing = T), ]
        
        plot(density(res_fda2[res_fda2$label == 'others',]$score), main = 'FDAcancer')
        if(nrow(res_fda2[res_fda2$label == 'fda_combo',])>1) lines(density(res_fda2[res_fda2$label == 'fda_combo',]$score), col = 'red')
        boxplot(score~label, data=res_fda2,main = 'FDAcancer')
        
        dt3 = data.frame(index = 'proof',
                         num_drug = nrow(res_fda2), 
                         num_proof = table(res_fda2$label)[1],
                         num_others = table(res_fda2$label)[2],
                         mean_proof = mean(res_fda2[res_fda2$label == 'fda_combo',]$score),
                         mean_others = mean(res_fda2[res_fda2$label == 'others',]$score),
                         effectsize = mean(res_fda2[res_fda2$label == 'fda_combo',]$score) - mean(res_fda2[res_fda2$label == 'others',]$score),
                         pval = wilcox.test(res_fda2[res_fda2$label == 'fda_combo',]$score, 
                                            res_fda2[res_fda2$label == 'others',]$score, 
                                            alternative = 'greater')$p.value)
        
        
        # tmp1 = res_fda2[res_fda2$label == 'fda_combo',]
        # tmp1 = unique(tmp1)
        
        res_fda3 = merge(res_fda2, drugTarget[c('drug_index','target')], by = 'drug_index')
        
        sharetarget = intersect(unique(res_fda3[res_fda3$label == 'fda_combo',]$target), unique(res_fda3[res_fda3$label == 'others',]$target))
        dst = unique(res_fda3[is.element(res_fda3$target, sharetarget),]$drug_index)
        tmp2 = res_fda3[res_fda3$label == 'others',]
        tmp2 = tmp2[!is.element(tmp2$drug_index, dst), ]
        dim(tmp2)
        tmp2 = unique(tmp2[!is.element(names(tmp2),'target')])
        tmp2 = unique(tmp2)
        dim(tmp2)
        tmp2 = tmp2[!is.element(tmp2$drug_index, unique(druglist$drug_index)),]
        dim(tmp2)
        
        res_fda3 = merge(res_fda2, drugMoA[c('drug_index','MoA')], by = 'drug_index')
        shareMoA = intersect(unique(res_fda3[res_fda3$label == 'fda_combo',]$MoA), unique(res_fda3[res_fda3$label == 'others',]$MoA))
        dst = unique(res_fda3[is.element(res_fda3$MoA, shareMoA),]$drug_index)
        # tmp2 = res_fda3[res_fda3$label == 'others',]
        tmp2 = tmp2[!is.element(tmp2$drug_index, dst), ]
        dim(tmp2)
        tmp2 = unique(tmp2[!is.element(names(tmp2),'MoA')])
        tmp2 = unique(tmp2)
        dim(tmp2)
        tmp2 = tmp2[!is.element(tmp2$drug_index, unique(druglist$drug_index)),]
        dim(tmp2)
        
        res2 = unique(rbind(tmp1, tmp2))
        res2 = res2[order(res2$score, decreasing = T),]
        res2$rank2 = 1:nrow(res2)
        res2$rank2_tile = res2$rank2 / nrow(res2)
        
        dt4 = data.frame(index = 'proof_filterMoA',
                         num_drug = nrow(res2), 
                         num_proof = table(res2$label)[1],
                         num_others = table(res2$label)[2],
                         mean_proof = mean(res2[res2$label == 'fda_combo',]$score),
                         mean_others = mean(res2[res2$label == 'others',]$score),
                         effectsize = mean(res2[res2$label == 'fda_combo',]$score) - mean(res2[res2$label == 'others',]$score),
                         pval = wilcox.test(res2[res2$label == 'fda_combo',]$score, 
                                            res2[res2$label == 'others',]$score, 
                                            alternative = 'greater')$p.value)
        
        plot(density(res2[res2$label == 'others',]$score), main='FDAcancer_filterMoA')
        if(nrow(res2[res2$label == 'fda_combo',])>1) lines(density(res2[res2$label == 'fda_combo',]$score), col = 'red')
        boxplot(score~label, data=res2, main='FDAcancer_filterMoA')
        
        # boxplot(rank2_tile~label, data=res2)
        # plot(density(res2[res2$label == 'others',]$rank2_tile))
        # lines(density(res2[res2$label == 'fda_combo',]$rank2_tile), col = 'red')
        
        if(nrow(tmp2)>nrow(tmp1)){
          
          set.seed(1234)
          testp = lapply(as.list(1:1000), function(i){
            
            tmpindex = sample( 1:nrow(tmp2), nrow(tmp1), replace = FALSE, prob = NULL)
            tmp3 = tmp2[tmpindex,]
            res3 = rbind(tmp1, tmp3)
            res3 = res3[order(res3$score, decreasing = T),]
            res3$rank2 = 1:nrow(res3)
            res3$rank2_tile = res3$rank2 / nrow(res3)
            
            dt5 = data.frame(index = 'proof_filterMoA_sample',
                             num_drug = nrow(res3), 
                             num_proof = table(res3$label)[1],
                             num_others = table(res3$label)[2],
                             mean_proof = mean(res3[res3$label == 'fda_combo',]$score),
                             mean_others = mean(res3[res3$label == 'others',]$score),
                             effectsize = mean(res3[res3$label == 'fda_combo',]$score) - mean(res3[res3$label == 'others',]$score),
                             pval = wilcox.test(res3[res3$label == 'fda_combo',]$score, 
                                                res3[res3$label == 'others',]$score, 
                                                alternative = 'greater')$p.value)
            
            return(dt5)
            
          })
          testp = do.call(rbind, testp)
          write.csv(testp, paste0(savepath, "SKCM_Primary_", proofType,'_diffPermTest.csv'), row.names = F, quote = F)
          
          plot(density(testp$pval), main = paste0('p<0.05: ', nrow(testp[testp$pval < 0.05,])))
          abline(v=0.05)
          
        }
        
        dev.off()
        
        
        dt = rbind(dt1, dt2, dt3, dt4)
        write.csv(dt, paste0(savepath, "SKCM_Primary_", proofType,'_diffTest.csv'), row.names = F, quote = F)
        
    }
  }
}

  
  
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# merge files and plot ---------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
library(ggplot2)
library(ggsci)
library(RColorBrewer)
library(patchwork)


# 70138
filepath <- '07_plotting_v2/05_drugITSprofiler_CombatRef_metaPlusFilter/drugrank_diff_isgold_lincsPortal2_metaPlusFilter/70138/'
filenames <- list.files(filepath)
# filenames=filenames[-c(1,2)]
for(i in seq(length(filenames))){
  # i = 2
  f=filenames[i]
  
  cancertype <- list.files(paste0(filepath, f))
  cancertype <- do.call(c, lapply(strsplit(cancertype,'_'),function(x)x[1]))
  cancertype = c(unique(cancertype),'SKCM_Primary')

  dir.create('07_plotting_v2/05_drugITSprofiler_CombatRef_metaPlusFilter/drugrank_diff_isgold_lincsPortal2_metaPlusFilter_summary/')
  dir.create('07_plotting_v2/05_drugITSprofiler_CombatRef_metaPlusFilter/drugrank_diff_isgold_lincsPortal2_metaPlusFilter_summary/70138/')

  ressavepath = paste0('07_plotting_v2/05_drugITSprofiler_CombatRef_metaPlusFilter/drugrank_diff_isgold_lincsPortal2_metaPlusFilter_summary/70138/', f,'/')
  dir.create(ressavepath)

  savepath = paste0(filepath, f,'/')
  resfilename = list.files(savepath)
  res70138 = list()

  for(j in seq(length(cancertype))){
    # j = 9
    cancer = cancertype[j]
      
    if(cancer == 'SKCM'){
      fname <- resfilename[grep('SKCM', resfilename)]
      fname <- fname[-grep('SKCM_Primary', fname)]
      fname <- fname[grep('_diffTest.csv', fname)]
    }else if(cancer == 'SKCM_Primary'){
      fname <- resfilename[grep('SKCM_Primary', resfilename)]
      fname <- fname[grep('_diffTest.csv', fname)]
    }else{
      fname <- resfilename[grep(cancer, resfilename)]
      fname <- fname[grep('_diffTest.csv', fname)]
    }
        
    res = lapply(as.list(fname), function(x){
      df = read.csv(paste0(savepath, x))
      df$cancer = cancer
      df$proofType = gsub("_diffTest.csv", "", x)
      df$shape = gsub(paste0(cancer, '_'), '', df$proofType)
      # df$shape = gsub(paste0(cancer, '_Primary_'), '', df$shape)
      
      return(df)
    })

    res =  do.call(rbind, res)   
    res[which(res$shape == gsub('_Primary','', cancer)), ]$shape = 'cancer'

    res70138[[j]] = res
  }

  res70138 = do.call(rbind, res70138)
  write.csv(res70138, paste0(ressavepath, "difftest_summary.csv"), quote = F, row.names = F)

  # plot
  lapply(as.list(unique(res70138$index)), function(x){

    pdata = res70138[res70138$index == x,]
    p = lapply(as.list(unique(res70138$shape)), function(y){

    pdata2 = pdata[pdata$shape == y, ]
    p = ggplot(data = pdata2, aes(x=effectsize, y= -log10(pval), colour=cancer,fill=shape)) + # , shape=shape
      geom_point( aes(size = num_proof)) +
      geom_hline(aes(yintercept = -log10(0.05)), colour = 'grey') +
      geom_vline(aes(xintercept = 0), colour = 'grey') +
      scale_color_brewer(palette = "Spectral") +
      scale_fill_brewer(palette = "Spectral") +
      theme_bw() + 
      scale_color_igv()+
        ggtitle(y)

      ggsave(paste0(ressavepath, x, "_", y, ".pdf"), p , width = 8, height = 6)      
      return(p)
    })
    
    p2 = (p[[1]]+p[[2]]) / (p[[3]] + p[[4]]) / (p[[5]]+p[[6]])
    ggsave(paste0(ressavepath, x, ".pdf"), p2 , width = 15, height = 18)
  })

}


# 92742
filepath <- '07_plotting_v2/05_drugITSprofiler_CombatRef_metaPlusFilter/drugrank_diff_isgold_lincsPortal2_metaPlusFilter/92742/'
filenames <- list.files(filepath)

for(i in seq(length(filenames))){
  # i = 2
  f=filenames[i]
  
  cancertype <- list.files(paste0(filepath, f))
  cancertype <- do.call(c, lapply(strsplit(cancertype,'_'),function(x)x[1]))
  cancertype = c(unique(cancertype),'SKCM_Primary')

  dir.create('07_plotting_v2/05_drugITSprofiler_CombatRef_metaPlusFilter/drugrank_diff_isgold_lincsPortal2_metaPlusFilter_summary/')
  dir.create('07_plotting_v2/05_drugITSprofiler_CombatRef_metaPlusFilter/drugrank_diff_isgold_lincsPortal2_metaPlusFilter_summary/92742/')

  ressavepath = paste0('07_plotting_v2/05_drugITSprofiler_CombatRef_metaPlusFilter/drugrank_diff_isgold_lincsPortal2_metaPlusFilter_summary/92742/', f,'/')
  dir.create(ressavepath)

  savepath = paste0(filepath, f,'/')
  resfilename = list.files(savepath)
  res92742 = list()

  for(j in seq(length(cancertype))){
    # j = 9
    cancer = cancertype[j]
      
    if(cancer == 'SKCM'){
      fname <- resfilename[grep('SKCM', resfilename)]
      fname <- fname[-grep('SKCM_Primary', fname)]
      fname <- fname[grep('_diffTest.csv', fname)]
    }else if(cancer == 'SKCM_Primary'){
      fname <- resfilename[grep('SKCM_Primary', resfilename)]
      fname <- fname[grep('_diffTest.csv', fname)]
    }else{
      fname <- resfilename[grep(cancer, resfilename)]
      fname <- fname[grep('_diffTest.csv', fname)]
    }
        
    res = lapply(as.list(fname), function(x){
      df = read.csv(paste0(savepath, x))
      df$cancer = cancer
      df$proofType = gsub("_diffTest.csv", "", x)
      df$shape = gsub(paste0(cancer, '_'), '', df$proofType)
      # df$shape = gsub(paste0(cancer, '_Primary_'), '', df$shape)
      
      return(df)
    })

    res =  do.call(rbind, res)   
    if(nrow(res[which(res$shape == gsub('_Primary','', cancer)), ])>0)res[which(res$shape == gsub('_Primary','', cancer)), ]$shape = 'cancer'

    res92742[[j]] = res
  }

  res92742 = do.call(rbind, res92742)
  write.csv(res92742, paste0(ressavepath, "difftest_summary.csv"), quote = F, row.names = F)

  # plot

  lapply(as.list(unique(res92742$index)), function(x){

    pdata = res92742[res92742$index == x,]
    p = lapply(as.list(unique(res92742$shape)), function(y){

    pdata2 = pdata[pdata$shape == y, ]
    p = ggplot(data = pdata2, aes(x=effectsize, y= -log10(pval), colour=cancer,fill=shape)) + # , shape=shape
      geom_point( aes(size = num_proof)) +
      geom_hline(aes(yintercept = -log10(0.05)), colour = 'grey') +
      geom_vline(aes(xintercept = 0), colour = 'grey') +
      scale_color_brewer(palette = "Spectral") +
      scale_fill_brewer(palette = "Spectral") +
      theme_bw() + 
      scale_color_igv()+
      ggtitle(y)
      ggsave(paste0(ressavepath, x, "_", y, ".pdf"), p , width = 8, height = 6)      
      return(p)
    })
    
    p2 = (p[[1]]+p[[2]]) / (p[[3]] + p[[4]]) / (p[[5]]+p[[6]])
    ggsave(paste0(ressavepath, x, ".pdf"), p2 , width = 15, height = 18)
  })

}


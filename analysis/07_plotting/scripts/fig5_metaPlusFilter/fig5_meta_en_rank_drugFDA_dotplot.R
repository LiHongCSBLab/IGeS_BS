# rm(list=ls())
work_path <- "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/"
setwd(work_path)
options(stringsAsFactors = F)
library(dplyr)
library(ggrepel)

dir.create("07_plotting_v2/05_drugITSprofiler_CombatRef_metaPlusFilter/drugITSres_metaPlusFilter_drugFDA_dotplot/")

drugMoA <- read.csv('/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/s3_lincs_drug_MoA/lincs_portal2.0/drugMoA_lincsPortal2_2.csv',header = T)
drugTarget <- read.csv('/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/s3_lincs_drug_MoA/lincs_portal2.0/drugtarget_lincsPortal2.csv', header = T)
drug_maxfda <- read.csv('/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/s3_lincs_drug_MoA/lincs_portal2.0/drug_max_fda_phase_lincsPortal2.csv', header = T)
drug_maxfda34 <- unique(drug_maxfda[is.element(drug_maxfda$max_fda_phase, c(3,4)), ])
drug_maxfda34$drugname_lower <- tolower(drug_maxfda34$pert_iname)

load( "07_plotting_v2/data/drugTarget_merged3DB.Rdata")
drug_FDA <- unique(drugTarget_FDA[c(1,5,6)])
drug_FDA$drug_index = NA
drug_FDA <- unique(drug_FDA[c(1,4,3,2)])
names(drug_FDA) <- names(drug_maxfda34)
drug_FDA$max_fda_phase = 'approved'
drug_maxfda34 <- rbind(drug_maxfda34,drug_FDA)

# ------------------------------------------------------------------------------
# run for dataset GSE70138

dataset = '70138'
filepath <- "06_results_summaryandplotting/results_meta_filter_ITSssgsea_CombatRef_17datasets/TUMERIC_m2_meta_en_ACATP0.05/70138/allgenes/"
cancerlist <- list.files(filepath)

filenames <- list.files(paste0(filepath, "SKCM_Metastatic/drugITSRank/"))
filenames <- filenames[grep('.txt', filenames)]
if(length(grep("unweight",  filenames))>0)filenames <- filenames[-grep('unweight',filenames)]

for(f in filenames){
  # f=filenames[3]
  for(cancer in cancerlist){
    
    # cancer = 'SKCM_Metastatic' 
    # filenames <- list.files("07_plotting_v2/05_drugITSprofiler_CombatRef_metaPlusFilter/drugITSres_metaPlusFilter_drugFDA_dotplot/")
    # drugrank <- read.csv("07_plotting_v2/05_drugITSprofiler_CombatRef_metaPlusFilter/drugITSres_metaPlusFilter_drugFDA_dotplot/meta_en_ACATP0.05_meta_weight_weighted_sum_mean_res.txt", sep = '\t', header = T)
    # dat <- drugrank[drugrank$cancer == cancer & drugrank$dataset == dataset, ]
    
    cancer_type = strsplit(cancer,'_')[[1]][1]
    cancer_subtype = strsplit(cancer,'_')[[1]][2]
    
    drugrankFile = paste0("06_results_summaryandplotting/results_meta_filter_ITSssgsea_CombatRef_17datasets/TUMERIC_m2_meta_en_ACATP0.05/",dataset,"/allgenes/",cancer,"/drugITSRank/",f)
    drugrank <- read.csv(drugrankFile, sep = '\t', header = T)
    drugrank$rank_tile = drugrank$rank / nrow(drugrank)
    drugrank$drugname_lower <- tolower(drugrank$pert_iname)
    dat <- drugrank[is.element(drugrank$drugname_lower, drug_maxfda34$drugname_lower), ]
    
    drug_matchproof <- read.csv("07_plotting_v2/05_drugITSprofiler_CombatRef_metaPlusFilter/m2_meta_en_ACATP0.05_metaPlusFilter/res_summary_ACATP0.05_matchProof.txt", sep = '\t', header = T)
    drug_matchproof <- drug_matchproof[drug_matchproof$proof!='animal',] 
    drug_matchproof <- drug_matchproof[drug_matchproof$if_with_ICB == 'yes', ]
    dat = dat %>%
      arrange(desc(score)) %>%
      mutate("x" = row_number())
    dat_shown <- drug_matchproof[grep(gsub("_res.txt","", f), drug_matchproof$strategy), ]
    drug_shown <- unique(dat_shown[dat_shown$cancer_type == cancer_type & 
                                     dat_shown$cancer_subtype == cancer_subtype, ]$drug_index)
    
    dat_shown <- dat[is.element(dat$drug_index,drug_shown), ]
    if(nrow(dat_shown) > 0){
      
      dat_shown$proof <- "withProof"
      if(nrow(dat_shown[dat_shown$rank_tile < 0.1,])>0)dat_shown[dat_shown$rank_tile < 0.1,]$proof <- "withProof_top10%"
      if(nrow(dat_shown[dat_shown$rank_tile < 0.05,])>0)dat_shown[dat_shown$rank_tile < 0.05,]$proof <- "withProof_top5%"
      # dat_shown2 <- dat[dat$rank_tile < 0.05,]
      # dat_shown2$proof = NA
      # dat_shown <- rbind(dat_shown,dat_shown2)
      
      # ------------------------------------------------------------------------------
      # dot rank plot, 
      # marking top 0.05%, FDA approved anti-cancer and other diseases treatment seperately dat.
      
      
      p_rank <- ggplot(dat, aes(x = x, y = score))+
        geom_point(color="#1762f7", size = 0.5)+
        theme_bw()+
        labs(x = "Ranked drugs", y = "Potential")+
        geom_text_repel(inherit.aes = F, 
                        data = dat_shown, 
                        aes(x = x, y = score, label = pert_iname, color = proof),
                        max.overlaps = getOption("ggrepel.max.overlaps", default = 200))+
        # scale_y_continuous(breaks = seq(0,1,0.2)) +
        ggtitle(paste0(cancer,"_",dataset))
      # p_rank
    }else{
      p_rank <- ggplot(dat, aes(x = x, y = score))+
        geom_point(color="#1762f7", size = 0.5)+
        theme_bw()+
        labs(x = "Ranked drugs", y = "Potential")+
        # scale_y_continuous(breaks = seq(0,1,0.2)) +
        ggtitle(paste0(cancer,"_",dataset))
      
    }
    
    ggsave(paste0("07_plotting_v2/05_drugITSprofiler_CombatRef_metaPlusFilter/drugITSres_metaPlusFilter_drugFDA_dotplot/", gsub("_res.txt", "", f), "_", cancer, "_",dataset, ".pdf"), 
           p_rank, height = 4, width = 5)
  }
  
  
  cancer = 'SKCM_Primary' 
  # filenames <- list.files("07_plotting_v2/05_drugITSprofiler_CombatRef_metaPlusFilter/drugITSres_metaPlusFilter_drugFDA_dotplot/")
  # drugrank <- read.csv("07_plotting_v2/05_drugITSprofiler_CombatRef_metaPlusFilter/drugITSres_metaPlusFilter_drugFDA_dotplot/meta_en_ACATP0.05_meta_weight_weighted_sum_mean_res.txt", sep = '\t', header = T)
  # dat <- drugrank[drugrank$cancer == cancer & drugrank$dataset == dataset, ]
  
  cancer_type = strsplit(cancer,'_')[[1]][1]
  cancer_subtype = strsplit(cancer,'_')[[1]][2]
  
  drugrankFile = paste0("06_results_summaryandplotting/results_meta_filter_ITSssgsea_CombatRef_17datasets/CPE_m2_meta_en_ACATP0.05/",dataset,"/allgenes/",cancer,"/drugITSRank/", f)
  drugrank <- read.csv(drugrankFile, sep = '\t', header = T)
  drugrank$rank_tile = drugrank$rank / nrow(drugrank)
    drugrank$drugname_lower <- tolower(drugrank$pert_iname)
    dat <- drugrank[is.element(drugrank$drugname_lower, drug_maxfda34$drugname_lower), ]
  
  drug_matchproof <- read.csv("07_plotting_v2/05_drugITSprofiler_CombatRef_metaPlusFilter/m2_meta_en_ACATP0.05_metaPlusFilter/res_summary_ACATP0.05_matchProof.txt", sep = '\t', header = T)
    drug_matchproof <- drug_matchproof[drug_matchproof$proof!='animal',] 
    drug_matchproof <- drug_matchproof[drug_matchproof$if_with_ICB == 'yes', ]
  dat = dat %>%
    arrange(desc(score)) %>%
    mutate("x" = row_number())
  dat_shown <- drug_matchproof[grep(gsub("_res.txt","", f), drug_matchproof$strategy), ]
  drug_shown <- unique(dat_shown[dat_shown$cancer_type == cancer_type & 
                                   dat_shown$cancer_subtype == cancer_subtype, ]$drug_index)
  
  dat_shown <- dat[is.element(dat$drug_index,drug_shown), ]
  dat_shown$proof <- "withProof"
  if(nrow(dat_shown[dat_shown$rank_tile < 0.1,])>0)dat_shown[dat_shown$rank_tile < 0.1,]$proof <- "withProof_top10%"
  if(nrow(dat_shown[dat_shown$rank_tile < 0.05,])>0)dat_shown[dat_shown$rank_tile < 0.05,]$proof <- "withProof_top5%"
  
  p_rank <- ggplot(dat, aes(x = x, y = score))+
    geom_point(color="#1762f7", size = 0.5)+
    theme_bw()+
    labs(x = "Ranked drugs", y = "Potential")+
    geom_text_repel(inherit.aes = F, 
                    data = dat_shown, 
                    aes(x = x, y = score, label = pert_iname, color = proof),
                    max.overlaps = getOption("ggrepel.max.overlaps", default = 200))+
    # scale_y_continuous(breaks = seq(0,1,0.2)) +
    ggtitle(paste0(cancer,"_",dataset))
  # p_rank
  ggsave(paste0("07_plotting_v2/05_drugITSprofiler_CombatRef_metaPlusFilter/drugITSres_metaPlusFilter_drugFDA_dotplot/", gsub("_res.txt", "", f), "_", cancer, "_",dataset, ".pdf"), 
         p_rank, height = 4, width = 5)
}


# ------------------------------------------------------------------------------
# run for dataset GSE92742

dataset = '92742'
filepath <- "06_results_summaryandplotting/results_meta_filter_ITSssgsea_CombatRef_17datasets/TUMERIC_m2_meta_en_ACATP0.05/92742/allgenes/"
cancerlist <- list.files(filepath)

filenames <- list.files(paste0(filepath, "SKCM_Metastatic/drugITSRank/"))
filenames <- filenames[grep('.txt', filenames)]
if(length(grep("unweight",  filenames))>0)filenames <- filenames[-grep('unweight',filenames)]

for(f in filenames){
  # f=filenames[2]
  for(cancer in cancerlist){
    
    # cancer = 'SKCM_Metastatic' 
    # filenames <- list.files("07_plotting_v2/05_drugITSprofiler_CombatRef_metaPlusFilter/drugITSres_metaPlusFilter_drugFDA_dotplot/")
    # drugrank <- read.csv("07_plotting_v2/05_drugITSprofiler_CombatRef_metaPlusFilter/drugITSres_metaPlusFilter_drugFDA_dotplot/meta_en_ACATP0.05_meta_weight_weighted_sum_mean_res.txt", sep = '\t', header = T)
    # dat <- drugrank[drugrank$cancer == cancer & drugrank$dataset == dataset, ]
    
    cancer_type = strsplit(cancer,'_')[[1]][1]
    cancer_subtype = strsplit(cancer,'_')[[1]][2]
    
    drugrankFile = paste0("06_results_summaryandplotting/results_meta_filter_ITSssgsea_CombatRef_17datasets/TUMERIC_m2_meta_en_ACATP0.05/",dataset,"/allgenes/",cancer,"/drugITSRank/",f)
    drugrank <- read.csv(drugrankFile, sep = '\t', header = T)
    drugrank$rank_tile = drugrank$rank / nrow(drugrank)
    drugrank$drugname_lower <- tolower(drugrank$pert_iname)
    dat <- drugrank[is.element(drugrank$drugname_lower, drug_maxfda34$drugname_lower), ]
    
    drug_matchproof <- read.csv("07_plotting_v2/05_drugITSprofiler_CombatRef_metaPlusFilter/m2_meta_en_ACATP0.05_metaPlusFilter/res_summary_ACATP0.05_matchProof.txt", sep = '\t', header = T)
    drug_matchproof <- drug_matchproof[drug_matchproof$proof!='animal',] 
    drug_matchproof <- drug_matchproof[drug_matchproof$if_with_ICB == 'yes', ]
    dat = dat %>%
      arrange(desc(score)) %>%
      mutate("x" = row_number())
    dat_shown <- drug_matchproof[grep(gsub("_res.txt","", f), drug_matchproof$strategy), ]
    drug_shown <- unique(dat_shown[dat_shown$cancer_type == cancer_type & 
                                     dat_shown$cancer_subtype == cancer_subtype, ]$drug_index)
    
    dat_shown <- dat[is.element(dat$drug_index,drug_shown), ]
    if(nrow(dat_shown) > 0){
      
      dat_shown$proof <- "withProof"
      if(nrow(dat_shown[dat_shown$rank_tile < 0.1,])>0)dat_shown[dat_shown$rank_tile < 0.1,]$proof <- "withProof_top10%"
      if(nrow(dat_shown[dat_shown$rank_tile < 0.05,])>0)dat_shown[dat_shown$rank_tile < 0.05,]$proof <- "withProof_top5%"
      # dat_shown2 <- dat[dat$rank_tile < 0.05,]
      # dat_shown2$proof = NA
      # dat_shown <- rbind(dat_shown,dat_shown2)
      
      # ------------------------------------------------------------------------------
      # dot rank plot, 
      # marking top 0.05%, FDA approved anti-cancer and other diseases treatment seperately dat.
      
      
      p_rank <- ggplot(dat, aes(x = x, y = score))+
        geom_point(color="#1762f7", size = 0.5)+
        theme_bw()+
        labs(x = "Ranked drugs", y = "Potential")+
        geom_text_repel(inherit.aes = F, 
                        data = dat_shown, 
                        aes(x = x, y = score, label = pert_iname, color = proof),
                        max.overlaps = getOption("ggrepel.max.overlaps", default = 200))+
        # scale_y_continuous(breaks = seq(0,1,0.2)) +
        ggtitle(paste0(cancer,"_",dataset))
      # p_rank
    }else{
      p_rank <- ggplot(dat, aes(x = x, y = score))+
        geom_point(color="#1762f7", size = 0.5)+
        theme_bw()+
        labs(x = "Ranked drugs", y = "Potential")+
        # scale_y_continuous(breaks = seq(0,1,0.2)) +
        ggtitle(paste0(cancer,"_",dataset))
      
    }
    ggsave(paste0("07_plotting_v2/05_drugITSprofiler_CombatRef_metaPlusFilter/drugITSres_metaPlusFilter_drugFDA_dotplot/", gsub("_res.txt", "", f), "_", cancer, "_",dataset, ".pdf"), 
           p_rank, height = 4, width = 5)
  }
  
  
  cancer = 'SKCM_Primary' 
  # filenames <- list.files("07_plotting_v2/05_drugITSprofiler_CombatRef_metaPlusFilter/drugITSres_metaPlusFilter_drugFDA_dotplot/")
  # drugrank <- read.csv("07_plotting_v2/05_drugITSprofiler_CombatRef_metaPlusFilter/drugITSres_metaPlusFilter_drugFDA_dotplot/meta_en_ACATP0.05_meta_weight_weighted_sum_mean_res.txt", sep = '\t', header = T)
  # dat <- drugrank[drugrank$cancer == cancer & drugrank$dataset == dataset, ]
  
  cancer_type = strsplit(cancer,'_')[[1]][1]
  cancer_subtype = strsplit(cancer,'_')[[1]][2]
  
  drugrankFile = paste0("06_results_summaryandplotting/results_meta_filter_ITSssgsea_CombatRef_17datasets/CPE_m2_meta_en_ACATP0.05/",dataset,"/allgenes/",cancer,"/drugITSRank/", f)
  drugrank <- read.csv(drugrankFile, sep = '\t', header = T)
  drugrank$rank_tile = drugrank$rank / nrow(drugrank)
    drugrank$drugname_lower <- tolower(drugrank$pert_iname)
    dat <- drugrank[is.element(drugrank$drugname_lower, drug_maxfda34$drugname_lower), ]
  
  drug_matchproof <- read.csv("07_plotting_v2/05_drugITSprofiler_CombatRef_metaPlusFilter/m2_meta_en_ACATP0.05_metaPlusFilter/res_summary_ACATP0.05_matchProof.txt", sep = '\t', header = T)
    drug_matchproof <- drug_matchproof[drug_matchproof$proof!='animal',] 
    drug_matchproof <- drug_matchproof[drug_matchproof$if_with_ICB == 'yes', ]
  dat = dat %>%
    arrange(desc(score)) %>%
    mutate("x" = row_number())
  dat_shown <- drug_matchproof[grep(gsub("_res.txt","", f), drug_matchproof$strategy), ]
  drug_shown <- unique(dat_shown[dat_shown$cancer_type == cancer_type & 
                                   dat_shown$cancer_subtype == cancer_subtype, ]$drug_index)
  
  dat_shown <- dat[is.element(dat$drug_index,drug_shown), ]
  dat_shown$proof <- "withProof"
  if(nrow(dat_shown[dat_shown$rank_tile < 0.1,])>0)dat_shown[dat_shown$rank_tile < 0.1,]$proof <- "withProof_top10%"
  if(nrow(dat_shown[dat_shown$rank_tile < 0.05,])>0)dat_shown[dat_shown$rank_tile < 0.05,]$proof <- "withProof_top5%"
  
  p_rank <- ggplot(dat, aes(x = x, y = score))+
    geom_point(color="#1762f7", size = 0.5)+
    theme_bw()+
    labs(x = "Ranked drugs", y = "Potential")+
    geom_text_repel(inherit.aes = F, 
                    data = dat_shown, 
                    aes(x = x, y = score, label = pert_iname, color = proof),
                    max.overlaps = getOption("ggrepel.max.overlaps", default = 200))+
    # scale_y_continuous(breaks = seq(0,1,0.2)) +
    ggtitle(paste0(cancer,"_",dataset))
  # p_rank
  ggsave(paste0("07_plotting_v2/05_drugITSprofiler_CombatRef_metaPlusFilter/drugITSres_metaPlusFilter_drugFDA_dotplot/", gsub("_res.txt", "", f), "_", cancer, "_",dataset, ".pdf"), 
         p_rank, height = 4, width = 5)
}


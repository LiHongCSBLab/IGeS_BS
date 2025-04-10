
work_path <- "/picb/bigdata/project/FengFYM/"
setwd(work_path)
options(stringsAsFactors = F)

library(getopt)

library(dplyr)
library(reshape2)
library(parallel)
library(patchwork)
# library(readxl)
library(stringr)

lincs_cell_line_weight <- read.csv("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/s1_cellline_selection/TCGA_110CL/correlation_celllines_sum.csv")
names(lincs_cell_line_weight) <- c("cancertype","cell_id","median_cor","mean_cor","sd_cor","max_cor","min_cor")

drug_index <- read.table("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/03_drug_immuneSig_enrichment/data/drug_index_LINCS.txt", sep = "\t", header = T)

dataset="70138" 
sign="positive"
num_gene=200
tumor_purity_method='TUMERIC'
sampleType='Primary'

res_70138=list()
cancerlist=c("BRCA","COAD","LIHC","LUAD","PAAD","PRAD","READ")
for( j in 1:length(cancerlist)){
  
  cancer=cancerlist[j]
  resPath <- paste0("p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/03_drug_immuneSig_enrichment/results_",tumor_purity_method, "_meta_oe_original/drug_immunesig_",tumor_purity_method, "_GSEA_result/")
  savePath <- paste0("p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/03_drug_immuneSig_enrichment/results_",tumor_purity_method, "_meta_oe_original/")
  savePath_allresult = paste0(savePath, "/drug_immunesig_",tumor_purity_method, "_GSEA_result_all/")
  saveAllPath_mergedresult = paste0(savePath, "/drug_immunesig_",tumor_purity_method, "_GSEA_result_summarized/")
  
  
  cancer_type = paste0(cancer,"_", sampleType)
  resPath <- paste0(resPath, "/", dataset, "/",cancer_type,"/", sign, "_", num_gene, "/")
  drugs <- list.files(resPath)
  
  length(drugs)
  res=list()
  for(i in 1:length(drugs)){
    drug=drugs[i]
    # drug = 'drug101'
    resPath2 <- paste0(resPath, "/", drug, "/")
    files <- list.files(resPath2)
    if(length(files) != 0){
      
      # files_treated <- grep("gsea_report_for_treated_", files, value = T)
      # files_DMSO <- grep("gsea_report_for_DMSO_", files, value = T)
      druginfo <- gsub(".xls", "",
                       gsub("gsea_report_for_DMSO_", "",
                            gsub("gsea_report_for_treated_", "", files)
                       ))
      
      druginfo <- strsplit(druginfo, "_", fixed = TRUE)
      druginfo <- do.call(rbind, druginfo)
      
      druginfo <- as.data.frame(druginfo)
      druginfo$drug_index <- drug
      druginfo$dataset <- dataset
      druginfo$files <- files
      druginfo <- druginfo[c("drug_index", "V1","V2","V3","dataset","files")]
      names(druginfo) <- c("drug_index", "cell_id", "treatment_time", "treatment_dose","dataset","files")
      res[[i]] = druginfo
    }
  }
  
  res_70138[[j]]=do.call(rbind, res)
  res_70138[[j]]$cancer=cancer
  res_70138[[j]]$sampleType=sampleType
}

cancer='SKCM'
sampleType='Primary'
tumor_purity_method='CPE'
j=j+1

resPath <- paste0("p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/03_drug_immuneSig_enrichment/results_",tumor_purity_method, "_meta_oe_original/drug_immunesig_",tumor_purity_method, "_GSEA_result/")
savePath <- paste0("p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/03_drug_immuneSig_enrichment/results_",tumor_purity_method, "_meta_oe_original/")
savePath_allresult = paste0(savePath, "/drug_immunesig_",tumor_purity_method, "_GSEA_result_all/")
saveAllPath_mergedresult = paste0(savePath, "/drug_immunesig_",tumor_purity_method, "_GSEA_result_summarized/")


cancer_type = paste0(cancer,"_", sampleType)

resPath <- paste0(resPath, "/", dataset, "/",cancer_type,"/", sign, "_", num_gene, "/")

drugs <- list.files(resPath)

length(drugs)
res=list()
for(i in 1:length(drugs)){
  drug=drugs[i]
  # drug = 'drug101'
  resPath2 <- paste0(resPath, "/", drug, "/")
  files <- list.files(resPath2)
  if(length(files) != 0){
    
    # files_treated <- grep("gsea_report_for_treated_", files, value = T)
    # files_DMSO <- grep("gsea_report_for_DMSO_", files, value = T)
    druginfo <- gsub(".xls", "",
                     gsub("gsea_report_for_DMSO_", "",
                          gsub("gsea_report_for_treated_", "", files)
                     ))
    
    druginfo <- strsplit(druginfo, "_", fixed = TRUE)
    druginfo <- do.call(rbind, druginfo)
    
    druginfo <- as.data.frame(druginfo)
    druginfo$drug_index <- drug
    druginfo$dataset <- dataset
    druginfo$files <- files
    druginfo <- druginfo[c("drug_index", "V1","V2","V3","dataset","files")]
    names(druginfo) <- c("drug_index", "cell_id", "treatment_time", "treatment_dose","dataset","files")
    res[[i]] = druginfo
  }
}
res_70138[[j]]=do.call(rbind, res)
res_70138[[j]]$cancer=cancer
res_70138[[j]]$sampleType=sampleType

cancer='SKCM'
sampleType='Metastatic'
tumor_purity_method='TUMERIC'
j=j+1

resPath <- paste0("p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/03_drug_immuneSig_enrichment/results_",tumor_purity_method, "_meta_oe_original/drug_immunesig_",tumor_purity_method, "_GSEA_result/")
savePath <- paste0("p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/03_drug_immuneSig_enrichment/results_",tumor_purity_method, "_meta_oe_original/")
savePath_allresult = paste0(savePath, "/drug_immunesig_",tumor_purity_method, "_GSEA_result_all/")
saveAllPath_mergedresult = paste0(savePath, "/drug_immunesig_",tumor_purity_method, "_GSEA_result_summarized/")


cancer_type = paste0(cancer,"_", sampleType)

resPath <- paste0(resPath, "/", dataset, "/",cancer_type,"/", sign, "_", num_gene, "/")

drugs <- list.files(resPath)

length(drugs)
res=list()
for(i in 1:length(drugs)){
  drug=drugs[i]
  # drug = 'drug101'
  resPath2 <- paste0(resPath, "/", drug, "/")
  files <- list.files(resPath2)
  if(length(files) != 0){
    
    # files_treated <- grep("gsea_report_for_treated_", files, value = T)
    # files_DMSO <- grep("gsea_report_for_DMSO_", files, value = T)
    druginfo <- gsub(".xls", "",
                     gsub("gsea_report_for_DMSO_", "",
                          gsub("gsea_report_for_treated_", "", files)
                     ))
    
    druginfo <- strsplit(druginfo, "_", fixed = TRUE)
    druginfo <- do.call(rbind, druginfo)
    
    druginfo <- as.data.frame(druginfo)
    druginfo$drug_index <- drug
    druginfo$dataset <- dataset
    druginfo$files <- files
    druginfo <- druginfo[c("drug_index", "V1","V2","V3","dataset","files")]
    names(druginfo) <- c("drug_index", "cell_id", "treatment_time", "treatment_dose","dataset","files")
    res[[i]] = druginfo
  }
}
res_70138[[j]]=do.call(rbind, res)
res_70138[[j]]$cancer=cancer
res_70138[[j]]$sampleType=sampleType
res_70138_all=do.call(rbind,res_70138)
write.csv(res_70138_all,'/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/03_drug_immuneSig_enrichment/usedLINCS_drugExpressionMatix_70138.csv', row.names=F, quote=F)
res_70138_all <- read.csv('/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/03_drug_immuneSig_enrichment/usedLINCS_drugExpressionMatix_70138.csv')
length(unique(res_70138_all$cell_id))
length(unique(res_70138_all$drug_index))
head(res_70138_all)
nrow(unique(res_70138_all[c("drug_index","files")]))




dataset="92742" 
sign="positive"
num_gene=200
tumor_purity_method='TUMERIC'
sampleType='Primary'

res_92742=list()
cancerlist=c("BRCA","COAD","LIHC","LUAD","LUSC","OV","PRAD","READ","STAD","UCEC")
for( j in 1:length(cancerlist)){
  
  cancer=cancerlist[j]
  resPath <- paste0("p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/03_drug_immuneSig_enrichment/results_",tumor_purity_method, "_meta_oe_original/drug_immunesig_",tumor_purity_method, "_GSEA_result/")
  savePath <- paste0("p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/03_drug_immuneSig_enrichment/results_",tumor_purity_method, "_meta_oe_original/")
  savePath_allresult = paste0(savePath, "/drug_immunesig_",tumor_purity_method, "_GSEA_result_all/")
  saveAllPath_mergedresult = paste0(savePath, "/drug_immunesig_",tumor_purity_method, "_GSEA_result_summarized/")
  
  
  cancer_type = paste0(cancer,"_", sampleType)
  resPath <- paste0(resPath, "/", dataset, "/",cancer_type,"/", sign, "_", num_gene, "/")
  drugs <- list.files(resPath)
  
  length(drugs)
  res=list()
  for(i in 1:length(drugs)){
    drug=drugs[i]
    # drug = 'drug101'
    resPath2 <- paste0(resPath, "/", drug, "/")
    files <- list.files(resPath2)
    if(length(files) != 0){
      
      # files_treated <- grep("gsea_report_for_treated_", files, value = T)
      # files_DMSO <- grep("gsea_report_for_DMSO_", files, value = T)
      druginfo <- gsub(".xls", "",
                       gsub("gsea_report_for_DMSO_", "",
                            gsub("gsea_report_for_treated_", "", files)
                       ))
      
      druginfo <- strsplit(druginfo, "_", fixed = TRUE)
      druginfo <- do.call(rbind, druginfo)
      
      druginfo <- as.data.frame(druginfo)
      druginfo$drug_index <- drug
      druginfo$dataset <- dataset
      druginfo$files <- files
      druginfo <- druginfo[c("drug_index", "V1","V2","V3","dataset","files")]
      names(druginfo) <- c("drug_index", "cell_id", "treatment_time", "treatment_dose","dataset","files")
      res[[i]] = druginfo
    }
  }
  
  res_92742[[j]]=do.call(rbind, res)
  res_92742[[j]]$cancer=cancer
  res_92742[[j]]$sampleType=sampleType
}

# save(res_92742,'/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/03_drug_immuneSig_enrichment/usedLINCS_drugExpressionMatix_92742.Rdata')

cancer='SKCM'
sampleType='Primary'
tumor_purity_method='CPE'
j=j+1

resPath <- paste0("p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/03_drug_immuneSig_enrichment/results_",tumor_purity_method, "_meta_oe_original/drug_immunesig_",tumor_purity_method, "_GSEA_result/")
savePath <- paste0("p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/03_drug_immuneSig_enrichment/results_",tumor_purity_method, "_meta_oe_original/")
savePath_allresult = paste0(savePath, "/drug_immunesig_",tumor_purity_method, "_GSEA_result_all/")
saveAllPath_mergedresult = paste0(savePath, "/drug_immunesig_",tumor_purity_method, "_GSEA_result_summarized/")


cancer_type = paste0(cancer,"_", sampleType)

resPath <- paste0(resPath, "/", dataset, "/",cancer_type,"/", sign, "_", num_gene, "/")

drugs <- list.files(resPath)

length(drugs)
res=list()

for(i in 1:length(drugs)){
  drug=drugs[i]
  # drug = 'drug101'
  resPath2 <- paste0(resPath, "/", drug, "/")
  files <- list.files(resPath2)
  if(length(files) != 0){
    
    # files_treated <- grep("gsea_report_for_treated_", files, value = T)
    # files_DMSO <- grep("gsea_report_for_DMSO_", files, value = T)
    druginfo <- gsub(".xls", "",
                     gsub("gsea_report_for_DMSO_", "",
                          gsub("gsea_report_for_treated_", "", files)
                     ))
    
    druginfo <- strsplit(druginfo, "_", fixed = TRUE)
    druginfo <- do.call(rbind, druginfo)
    
    druginfo <- as.data.frame(druginfo)
    druginfo$drug_index <- drug
    druginfo$dataset <- dataset
    druginfo$files <- files
    druginfo <- druginfo[c("drug_index", "V1","V2","V3","dataset","files")]
    names(druginfo) <- c("drug_index", "cell_id", "treatment_time", "treatment_dose","dataset","files")
    res[[i]] = druginfo
  }
}
res_92742[[j]]=do.call(rbind, res)
res_92742[[j]]$cancer=cancer
res_92742[[j]]$sampleType=sampleType




cancer='SKCM'
sampleType='Metastatic'
tumor_purity_method='TUMERIC'
j=j+1

resPath <- paste0("p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/03_drug_immuneSig_enrichment/results_",tumor_purity_method, "_meta_oe_original/drug_immunesig_",tumor_purity_method, "_GSEA_result/")
savePath <- paste0("p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/03_drug_immuneSig_enrichment/results_",tumor_purity_method, "_meta_oe_original/")
savePath_allresult = paste0(savePath, "/drug_immunesig_",tumor_purity_method, "_GSEA_result_all/")
saveAllPath_mergedresult = paste0(savePath, "/drug_immunesig_",tumor_purity_method, "_GSEA_result_summarized/")


cancer_type = paste0(cancer,"_", sampleType)

resPath <- paste0(resPath, "/", dataset, "/",cancer_type,"/", sign, "_", num_gene, "/")

drugs <- list.files(resPath)

length(drugs)
res=list()
for(i in 1:length(drugs)){
  drug=drugs[i]
  # drug = 'drug101'
  resPath2 <- paste0(resPath, "/", drug, "/")
  files <- list.files(resPath2)
  if(length(files) != 0){
    
    # files_treated <- grep("gsea_report_for_treated_", files, value = T)
    # files_DMSO <- grep("gsea_report_for_DMSO_", files, value = T)
    druginfo <- gsub(".xls", "",
                     gsub("gsea_report_for_DMSO_", "",
                          gsub("gsea_report_for_treated_", "", files)
                     ))
    
    druginfo <- strsplit(druginfo, "_", fixed = TRUE)
    druginfo <- do.call(rbind, druginfo)
    
    druginfo <- as.data.frame(druginfo)
    druginfo$drug_index <- drug
    druginfo$dataset <- dataset
    druginfo$files <- files
    druginfo <- druginfo[c("drug_index", "V1","V2","V3","dataset","files")]
    names(druginfo) <- c("drug_index", "cell_id", "treatment_time", "treatment_dose","dataset","files")
    res[[i]] = druginfo
  }
}
res_92742[[j]]=do.call(rbind, res)
res_92742[[j]]$cancer=cancer
res_92742[[j]]$sampleType=sampleType

res_927428_all=do.call(rbind,res_92742)
write.csv(res_927428_all,'/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/03_drug_immuneSig_enrichment/usedLINCS_drugExpressionMatix_92742.csv', row.names=F, quote=F)

print(paste0('num_cellline: ' ,length(unique(res_927428_all$cell_id))))
print(paste0('num_drug: ' ,length(unique(res_927428_all$drug_index))))
print(paste0('num_expression: ' ,nrow(unique(res_927428_all[c("drug_index","files")]))))

res_70138_all <- read.csv('/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/03_drug_immuneSig_enrichment/usedLINCS_drugExpressionMatix_70138.csv')
res_927428_all <- read.csv('/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/03_drug_immuneSig_enrichment/usedLINCS_drugExpressionMatix_92742.csv')
length(unique(c(res_70138_all$drug_index, res_927428_all$drug_index)))



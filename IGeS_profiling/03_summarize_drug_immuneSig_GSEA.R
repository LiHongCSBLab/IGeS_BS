# summarizing NES and P-value for multiple treatment time and dose
# 1) NES score merge method: 
#    OCTAD (https://github.com/Bin-Chen-Lab/RGES)
#
#
# 2) p.value merge method: 
#    Fisher"s method: https://github.com/reimandlab/ActivePathways/blob/master/R/merge_p.r
#    (DIDN"T USE)EmpiricalBrownsMethod: https://github.com/IlyaLab/CombiningDependentPvaluesUsingEBM/blob/master/R/EmpiricalBrownsMethod/R/ebm.R
#    ACAT: ACAT: A Fast and Powerful p Value Combination Method for Rare-Variant Analysis in Sequencing Studies
#    https://github.com/yaowuliu/ACAT/blob/master/R/ACAT.R
#    codes from these three are downloaded from: https://github.com/reimandlab/ActivePathways/blob/master/R/merge_p.r
#    additionally, we compute if num_sig(FDR<0.05) /num_sig > 0.5 as a naive method
# rm(list=ls())
options(stringsAsFactors = F)

library(getopt)

command=matrix(c( 
  "dataset",          "a", 1, "character", "Define which dataset to use, 70138 or 92742",
  "cancer",           "c", 1, "character", "Define cancer type, e.g. LIHC",
  "tumor_subtype",    "t", 1, "character", "Primary or Metastatic. Current for Metastatics only support for SKCM",
  "purity_method",    "u", 1, "character", "Define method used for computing tumor purity. TUMERIC or CPE. Current forfor primary SKCM, only CPE method is supported",
  "work_path",        "w", 1, "character",   "Define the path to work directory",
  "savePath",          "s", 1, "character",   "Define the path to save results",
  "help",              "h", 0, "logical",   "help file"),
  byrow=T,ncol=5)



args=getopt(command)

cancer = args$cancer
sampleType = args$tumor_subtype
tumor_purity_method = args$tumor_purity_method
dataset = args$dataset

work_path=args$work_path
savePath=args$savePath

# cancer='SKCM'
# sampleType = 'Primary'
# tumor_purity_method = 'CPE'
# cancer='LIHC'
# sampleType = 'Primary'
# tumor_purity_method = 'TUMERIC'
# dataset = 92742

# work_path <- "../IGeS_BS/IGeS_BS_tool"
setwd(work_path)

library(dplyr)
library(reshape2)
library(parallel)
library(patchwork)
library(stringr)
library(fs)


func1=path("./IGeS_profiling/", "merge_p.R")
source(func1)
func2=path("./IGeS_profiling/", "getsRGES.R")
source(func2)
# source("p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/03_drug_immuneSig_enrichment/scripts/merge_p.R")
# source("p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/03_drug_immuneSig_enrichment/scripts/getsRGES.R")


merge_result <- function(cancer,
                         sampleType,
                         dataset,
                         LINCS_type='allgenes' ,
                         tumor_purity_method = tumor_purity_method,
                         sign= 'positive',
                         num_gene=200,
                         resPath,
                         savePath_allresult,
                         saveAllPath_mergedresult) {
  
  #   cancer="SKCM"
  #   dataset="70138" 
  #   # cell="A375"
  #   sign="all"
  #   num_gene=400
  #   resPath = resPath
  #   savePath = savePath
  lincs_cell_line_weight_path=path("./Data/", "correlation_celllines_sum.csv")
  lincs_cell_line_weight <- read.csv(lincs_cell_line_weight_path)
  names(lincs_cell_line_weight) <- c("cancertype","cell_id","median_cor","mean_cor","sd_cor","max_cor","min_cor")
  
  drug_index_path=path("./Data/", "drug_index_LINCS.txt")
  drug_index <- read.table(drug_index_path, sep = "\t", header = T)
  
  cancer_type = paste0(cancer,"_", sampleType)
  
  if(LINCS_type == "allgenes"){
    resPath <- paste0(resPath, "/", dataset, "/",cancer_type,"/", sign, "_", num_gene, "/")
  }else if(LINCS_type == "978genes"){
    resPath <- paste0(resPath, "/", dataset, "/",cancer,"/", sign, "_", num_gene, "/")
    
  }
  
  drugs <- list.files(resPath)
  
  for(drug in drugs){
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
      
      # druginfo$indexs <- paste(druginfo$dataset, druginfo$cellline, druginfo$treatment_time, sep = "_")
      
      res = lapply(files, function(file) {
        if(file.info(paste0(resPath2, file))$size != 0){
          df <- read.csv(paste0(resPath2, file), sep = "\t", stringsAsFactors = F)
          if(nrow(df)>0){
            df <- df[c(1, 5:9)]
            df$immune_sig <- tolower(df$NAME)
            df$files <- file
            return(df)
          }
        }
      })
      res <- do.call(rbind, res)
      
      
      if(!is.null(res)) {
        
        res <- inner_join(druginfo, res)
        
        if(length(grep(".311", res$cell_id)) != 0 ){ 
          res[grep(".311", res$cell_id),]$cell_id = gsub(".311","", res[grep(".311", res$cell_id),]$cell_id)}
        if(length(grep(".101", res$cell_id)) != 0){ 
          res[grep(".101", res$cell_id),]$cell_id = gsub(".101","", res[grep(".101", res$cell_id),]$cell_id)}
        
        res_merged <- data.frame(drug_index = drug, 
                                 # cell_id = cell,
                                 sig = unique(res$immune_sig),
                                 mNES = 0, 
                                 NESmean = 0,
                                 NESn = 0,
                                 NESmedian = 0,
                                 NESsd = 0,
                                 NEScoef=0,
                                 mPvalueNaive = 0,
                                 mPvalue_fisher = 0,
                                 mPvalue_ACAT = 0)
        
        for (sig in unique(res$immune_sig)) {
          
          tmp <- res[res$immune_sig == sig, ]
          if(all(!is.na(unique(tmp$NES)))) {
            
            tmp <- tmp[order(tmp$treatment_dose),]
            tmp <- tmp[!is.na(tmp$NES),]
            # print(nrow(tmp[tmp$FDR.q.val < 0.1, ]))
            
            # merge multiple NES score
            # ref <- tmp[tmp$treatment_dose == 10 & tmp$treatment_time == 24, ]
            # target <- tmp[-which(tmp$treatment_dose == 10 & tmp$treatment_time == 24), ]
            mNES <- merge_NESmodified(
              resfile = tmp,
              cancer = cancer,
              lincs_cell_line_weight = lincs_cell_line_weight
            )
            
            res_merged[which(res_merged$sig == sig),]$mNES <- mNES$mNES
            res_merged[which(res_merged$sig == sig),]$NESn <- mNES$n
            res_merged[which(res_merged$sig == sig),]$NESmean <- mNES$mean
            res_merged[which(res_merged$sig == sig),]$NESmedian <- mNES$median
            res_merged[which(res_merged$sig == sig),]$NESmedian <- mNES$median
            res_merged[which(res_merged$sig == sig),]$NESsd <- mNES$sd
            res_merged[which(res_merged$sig == sig),]$NEScoef <- mNES$coef
            
            # merge multiple P values
            merged_FDRp_naive <- nrow(tmp[tmp$FDR.q.val<0.05,])/nrow(tmp)
            merged_FDRp_fisher <- merge_p_values(tmp$FDR.q.val, method = "Fisher")
            # merged_p_ebm <- brownsMethod(tmp$FDR.q.val)
            # tmp[tmp$NOM.p.val == 1, ]$NOM.p.val = 0.999
            # tmp[tmp$NOM.p.val == 0, ]$NOM.p.val = 10^-20
            
            # if (unique(tmp$FDR.q.val) == 0) {
            #   tmp$FDR.q.val = tmp$FDR.q.val + 10^-10
            #   merged_FDRp_ACAT <- ACAT(Pvals = tmp$FDR.q.val, weights = NULL)
            
            # } else {
            tmp[tmp$FDR.q.val == 0, ]$FDR.q.val = tmp[tmp$FDR.q.val == 0, ]$FDR.q.val + 10^-20
            tmp[tmp$FDR.q.val == 1, ]$FDR.q.val = tmp[tmp$FDR.q.val == 1, ]$FDR.q.val - 10^-20
            
            # weights = as.numeric(tmp$treatment_dose)
            # weights = weights / sum(weights)
            # # dose = as.numeric(tmp[tmp$FDR.q.val != 0 & tmp$FDR.q.val != 1, ]$treatment_dose)
            # # weights=rep(0.5/length(dose), length(dose))
            # # weights[which(dose==10)] = 0.5
            
            # # merged_FDRp_ACAT <- ACAT(Pvals = tmp[tmp$FDR.q.val != 0 & tmp$FDR.q.val != 1, ]$FDR.q.val, weights = weights)
            # # merged_FDRp_ACAT <- ACAT(Pvals = tmp$FDR.q.val, weights = weights)
            merged_FDRp_ACAT <- ACAT(Pvals = tmp$FDR.q.val, weights = NULL)
            # print(summary(tmp$FDR.q.val))
            # }
            
            res_merged[which(res_merged$sig == sig),]$mPvalueNaive <- merged_FDRp_naive
            res_merged[which(res_merged$sig == sig),]$mPvalue_fisher <- merged_FDRp_fisher
            res_merged[which(res_merged$sig == sig),]$mPvalue_ACAT <- merged_FDRp_ACAT
          }
        }
        # unique(res[res$NES > 0 & res$FDR.q.val < 0.05, ]$immune_sig) -> a1
        # unique(res[res$NES < 0 & res$FDR.q.val < 0.05, ]$immune_sig) -> a2
        # unique(sig1[sig1$mNES > 0,]$sig) -> b1
        # unique(sig1[sig1$mNES < 0, ]$sig) -> b2
        
        dir.create(paste0(savePath_allresult, "/", dataset))
        dir.create(paste0(saveAllPath_mergedresult, "/", dataset))
        dir.create(paste0(savePath_allresult, "/", dataset, "/", cancer_type, "/"))
        dir.create(paste0(saveAllPath_mergedresult, "/", dataset, "/", cancer_type, "/"))
        dir.create(paste0(savePath_allresult, "/", dataset, "/", cancer_type, "/",sign,"_",num_gene,"/"))
        dir.create(paste0(saveAllPath_mergedresult, "/", dataset, "/", cancer_type, "/",sign,"_",num_gene,"/"))
        
        write.csv(res, paste0(savePath_allresult, "/", dataset, "/", cancer_type, "/",sign,"_",num_gene,"/", drug, "_result_all.csv"), 
                  quote = F, row.names=F)
        write.csv(res_merged, paste0(saveAllPath_mergedresult,  "/",dataset, "/", cancer_type, "/",sign,"_",num_gene,"/", drug, "_result_merged.csv"), 
                  quote = F, row.names=F)
        
      }
    }
  }
}

# plot together locally


# resPath <- paste0("p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/03_drug_immuneSig_enrichment/results_",tumor_purity_method, "_meta_oe_original/drug_immunesig_",tumor_purity_method, "_GSEA_result_files/")
# savePath <- paste0("p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/03_drug_immuneSig_enrichment/results_",tumor_purity_method, "_meta_oe_original/")
resPath <- paste0(savePath, "drug_immunesig_",tumor_purity_method, "_GSEA_result_files/")
dir.create(savePath)
dir.create(paste0(savePath, "/drug_immunesig_",tumor_purity_method, "_GSEA_result_all/"))
dir.create(paste0(savePath, "/drug_immunesig_",tumor_purity_method, "_GSEA_result_summarized/"))
savePath_allresult = paste0(savePath, "/drug_immunesig_",tumor_purity_method, "_GSEA_result_all/")
saveAllPath_mergedresult = paste0(savePath, "/drug_immunesig_",tumor_purity_method, "_GSEA_result_summarized/")

merge_result(
  cancer = cancer,
  sampleType=sampleType,
  LINCS_type=LINCS_type,
  tumor_purity_method=tumor_purity_method,
  dataset = dataset,
  sign = sign,
  num_gene = num_gene,
  resPath = resPath,
  savePath_allresult = savePath_allresult,
  saveAllPath_mergedresult = saveAllPath_mergedresult
) 



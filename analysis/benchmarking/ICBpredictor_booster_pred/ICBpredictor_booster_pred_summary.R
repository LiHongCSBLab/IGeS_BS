# summarizing NES and P-value for multiple treatment time and dose
# 1) NES score merge method: 
#    OCTAD (https://github.com/Bin-Chen-Lab/RGES)
#
#
# rm(list=ls())
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


source("p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/03_drug_immuneSig_enrichment/scripts/merge_p.R")
source("p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/03_drug_immuneSig_enrichment/scripts/getsRGES.R")



merge_result <- function(cancer,
                         sampleType,
                         dataset,
                         LINCS_type,
                         tumor_purity_method = tumor_purity_method,
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
  
  lincs_cell_line_weight <- read.csv("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/s1_cellline_selection/TCGA_110CL/correlation_celllines_sum.csv")
  names(lincs_cell_line_weight) <- c("cancertype","cell_id","median_cor","mean_cor","sd_cor","max_cor","min_cor")
  
  drug_index <- read.table("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/03_drug_immuneSig_enrichment/data/drug_index_LINCS.txt", sep = "\t", header = T)

  cancer_type = paste0(cancer,"_", sampleType)
  # resPath <- paste0("p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/07_plotting_v2/ICBpredictor_booster_pred/")
  resPath <- paste0(resPath, "/", dataset, "/",cancer, "/")

  drugs <- list.files(resPath)
  
  for(drug in drugs){
    # drug = 'drug2'
    resPath2 <- paste0(resPath, "/", drug, "/")
    files <- list.files(resPath2)
    if(length(files) != 0){
    druginfo <- strsplit(files, "_", fixed = TRUE)
    druginfo <- do.call(rbind, druginfo)
    
    druginfo <- as.data.frame(druginfo)

    druginfo$drug_index <- drug
    druginfo$dataset <- dataset
    druginfo$files <- files
    druginfo <- druginfo[c("drug_index", "V1","V2","V3","dataset","files")]
    names(druginfo) <- c("drug_index", "cell_id", "treatment_time", "treatment_dose","dataset","files")
    
    # druginfo$indexs <- paste(druginfo$dataset, druginfo$cellline, druginfo$treatment_time, sep = "_")
    

    res = lapply(files, function(file) {
        f1='TIDE_TIP_result/TIDE_TIP_pred.csv'
        f2='IPS_result/IPS_pred.csv'
        f3='sig_OE_results/sig_OE_res.csv'

        if(file.info(paste0(resPath2, file, '/', f1))$size != 0){
            df1 <- read.csv(paste0(resPath2, file, '/', f1), stringsAsFactors = F)
        }
        if(file.info(paste0(resPath2, file, '/', f2))$size != 0){
            df2 <- read.csv(paste0(resPath2, file, '/', f2), stringsAsFactors = F)
        }
        # if(file.info(paste0(resPath2, file, '/', f3))$size != 0){
        #     df3 <- read.csv(paste0(resPath2, file, '/', f3), stringsAsFactors = F)
        # }
            
        df = rbind(df1, df2)#, df3)
        names(df) = c('immune_sig', 'NES')
        df$files <- file
        return(df)
    #   if(file.info(paste0(resPath2, file))$size != 0){
    #     df <- read.csv(paste0(resPath2, file), sep = "\t", stringsAsFactors = F)
    #     if(nrow(df)>0){
    #       df <- df[c(1, 5:9)]
    #       df$immune_sig <- tolower(df$NAME)
    #       df$files <- file
    #       return(df)
    #     }
    #   }
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
        
      }
    }

    dir.create(paste0(savePath_allresult, "/", dataset))
    dir.create(paste0(saveAllPath_mergedresult, "/", dataset))
    dir.create(paste0(savePath_allresult, "/", dataset, "/", cancer_type, "/"))
    dir.create(paste0(saveAllPath_mergedresult, "/", dataset, "/", cancer_type, "/"))
    
    write.csv(res, paste0(savePath_allresult, "/", dataset, "/", cancer_type, "/", drug, "_result_all.csv"), 
              quote = F, row.names=F)
    write.csv(res_merged, paste0(saveAllPath_mergedresult,  "/",dataset, "/", cancer_type, "/", drug, "_result_merged.csv"), 
              quote = F, row.names=F)
              
    }
    }
  }
}


cancer='LIHC'
sampleType = 'Primary'
tumor_purity_method = 'TUMERIC'
dataset = 92742
LINCS_type = 'allgenes' 


resPath <- paste0("p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/07_plotting_v2/ICBpredictor_booster_pred/")
savePath <- paste0("p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/07_plotting_v2/ICBpredictor_booster_pred/")
dir.create(savePath)
dir.create(paste0(savePath, "/result_all/"))
dir.create(paste0(savePath, "/result_summarized/"))
savePath_allresult = paste0(savePath, "/result_all/")
saveAllPath_mergedresult = paste0(savePath, "/result_summarized/")
    
merge_result(
      cancer = cancer,
      sampleType=sampleType,
      LINCS_type=LINCS_type,
      tumor_purity_method=tumor_purity_method,
      dataset = dataset,
      resPath = resPath,
      savePath_allresult = savePath_allresult,
      saveAllPath_mergedresult = saveAllPath_mergedresult
    ) 


# merge results in one file
  cancer_type = paste0(cancer,"_", sampleType)

files = list.files(paste0(saveAllPath_mergedresult,  "/",dataset, "/", cancer_type, "/"))
res_all = do.call(rbind, lapply(as.list(files), function(x){
                    read.csv(paste0(saveAllPath_mergedresult,  "/",dataset, "/", cancer_type, "/", x))
                }))

write.csv(res_all, paste0(saveAllPath_mergedresult,  "/",dataset, "/", cancer_type, "_all_result_merged.csv"), quote = F, row.names=F)
tmp=res_all[res_all$sig == 'TIDE', ]
head(tmp[order(tmp$mNES, decreasing=T), ],20)
write.csv(res_all, paste0(saveAllPath_mergedresult,  "/",dataset, "/", cancer_type, "_TIDE_merged.csv"), quote = F, row.names=F)

tmp=res_all[res_all$sig == 'IPS', ]
head(tmp[order(tmp$mNES, decreasing=T), ],20)
write.csv(res_all, paste0(saveAllPath_mergedresult,  "/",dataset, "/", cancer_type, "_IPS_merged.csv"), quote = F, row.names=F)

tmp=res_all[res_all$sig == 'TIP_signature', ]
head(tmp[order(tmp$mNES, decreasing=T), ],20)
write.csv(res_all, paste0(saveAllPath_mergedresult,  "/",dataset, "/", cancer_type, "_TIP_merged.csv"), quote = F, row.names=F)

tmp=res_all[res_all$sig == 'resu', ]
head(tmp[order(tmp$mNES, decreasing=T), ],20)
write.csv(res_all, paste0(saveAllPath_mergedresult,  "/",dataset, "/", cancer_type, "_OE_resu_merged.csv"), quote = F, row.names=F)

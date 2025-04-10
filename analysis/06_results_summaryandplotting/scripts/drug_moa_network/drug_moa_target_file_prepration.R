print("prepare files for ploting")
work_path <- "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/"

setwd(work_path)
options(stringsAsFactors = F)
library(parallel)
library(reshape2)
library(stringr)
library(dplyr)

# drug_index <- read.csv('06_results_summaryandplotting/data/drug_index_LINCS.txt',sep = "\t")
# # 1) select drug with both MoA and target genes
# drugtarget <- read.csv("06_results_summaryandplotting/data/drugTargetMerge.txt",sep = "\t")
# drugMoA <- read.csv("06_results_summaryandplotting/data/drugMoA_lincsPortal2.csv", sep = ',', header = F)
# drugMoA <- drugMoA[c(1,2,3)]
# drugMoA <- drugMoA[!is.na(drugMoA$V3),]
# names(drugMoA) <-c("pert_iname", "drug_index", "MoA")

# drugTargetMoA <- inner_join(drugMoA, drugtarget)
# drugTargetMoA <- unique(drugTargetMoA[-1])
# drugTargetMoA <- inner_join(drug_index, drugTargetMoA)
# drugMoA <- unique(drugTargetMoA[(c("pert_iname","drug_index","MoA"))])
# drugTarget <- unique(drugTargetMoA[(c("pert_iname","drug_index","target"))])

# write.table(drugTargetMoA, "06_results_summaryandplotting/data/drugTargetMoA_merged.txt", sep = '\t', row.names = F, quote = F)
# write.table(drugMoA, "06_results_summaryandplotting/data/drugTarget_final.txt", sep = '\t', row.names = F, quote = F)
# write.table(drugTarget, "06_results_summaryandplotting/data/drugMoA_final.txt", sep = '\t', row.names = F, quote = F)

# length(unique(drugTargetMoA$drug_index))
# [1] 1295
# length(unique(drugTargetMoA$MoA))
# [1] 301
# length(unique(drugTargetMoA$target))
# [1] 1584


# 2) for each cancer type, extract target genes correlation analysis results 
# with immune signatures
drugTargetMoA <- read.table("06_results_summaryandplotting/data/drugTargetMoA_merged.txt", sep = '\t', header = T)
dir.create("06_results_summaryandplotting/results/targets_MoA_immusig/")

for(tumor_purity_method in c("IHC")){ # , "TUMERIC", "CPE"
  
  # tumor_purity_method= "TUMERIC"
  method = "pearson"
  
  dir.create(paste0("06_results_summaryandplotting/results/targets_MoA_immusig/", tumor_purity_method))
purity_method = 'TUMERIC'
if(purity_method == 'IHC'){
  resultpath = paste0("02_tumorGene_immuneSig_correlation/results/immunecell_TCGAgenes_pcor_pearson_", tumor_purity_method, "/")

}else{
  resultpath = paste0("02_tumorGene_immuneSig_correlation/results/immunecell_TCGAgenes_result_pearson_", tumor_purity_method, "/")

}




  filelist = list.files(resultpath)
  
  cancerlist1 = filelist[grep('pcor', filelist)]
  cancerlist2 = filelist[-grep('pcor', filelist)]
   
  for (file_name in cancerlist1) {
    savefile = gsub(paste0("pcor_",method,"_"),"",gsub(".Rdata","", file_name))
    savepath = paste0("06_results_summaryandplotting/results/targets_MoA_immusig/", tumor_purity_method,"/", savefile)
    
    file = paste0(resultpath, file_name)
    load(file)
    
    result_all <- lapply(cor_result, function(x){
      # x=cor_result[[1]]
      x = data.frame(r = unlist(x$cor),
                     p.value = unlist(x$p.value),
                     adj.p = unlist(x$p.adj))
      x$target = rownames(x)
      x = inner_join(drugTargetMoA,x)
      return(x)
    })
    
    
    result_significant_loose <- lapply(cor_result, function(x){
      # x=cor_result[[1]]
      x = data.frame(r = unlist(x$cor),
                     p.value = unlist(x$p.value),
                     adj.p = unlist(x$p.adj))
      x$target = rownames(x)
      x = inner_join(drugTargetMoA,x)
      # x = inner_join(drugTarget['target'],x)
      x <- x[x$adj.p < 0.05,]
      
      return(x)
    })
    
    result_significant_strict <- lapply(cor_result, function(x){
      # x=cor_result[[1]]
      x = data.frame(r = unlist(x$cor),
                     p.value = unlist(x$p.value),
                     adj.p = unlist(x$p.adj))
      x$target = rownames(x)
      x = inner_join(drugTargetMoA,x)
      # x = inner_join(drugTarget['target'],x)
      x <- x[x$adj.p < 0.05,]
      s1 <- summary(x[x$r > 0,]$r)[5]
      s2 <- summary(x[x$r < 0,]$r)[2]
      y <- unique(x[(x$r > s1 | x$r < s2) & x$adj.p < 0.05,])
      # y <- y[!is.na(y$r),]
      return(y)
    })
    
    save(result_all, 
         result_significant_loose, 
         result_significant_strict,
         file = paste0(savepath,".Rdata"))
    
  }
  
  for(file_name in cancerlist2){
    
    savefile = gsub(paste0("cor_",method,"_"),"",gsub(".Rdata","", file_name))
    savepath = paste0("06_results_summaryandplotting/results/targets_MoA_immusig/", tumor_purity_method,"/", savefile)
    # dir.create(savepath)
    
    file = paste0(resultpath, file_name)
    load(file)
    
    
    result_all <- lapply(cor_result2, function(x){
      x$target = rownames(x)
      x = inner_join(drugTargetMoA, x)
      
    })
    
    result_significant_strict <- lapply(cor_result2, function(x){
      x$target = rownames(x)
      x = inner_join(drugTargetMoA,x)
      
      x <- x[x$adj.p < 0.05,]
    })
    result_significant_strict <- lapply(cor_result2, function(x){
      x$target = rownames(x)
      x = inner_join(drugTargetMoA,x)
      
      x <- x[x$adj.p < 0.05,]
      s1 <- summary(x[x$r > 0,]$r)[5]
      s2 <- summary(x[x$r < 0,]$r)[2]
      y <- unique(x[(x$r > s1 | x$r < s2) & x$adj.p < 0.05,])
      # y <- y[!is.na(y$r),]
      return(y)
    })
    
    save(result_all, 
         result_significant_loose, 
         result_significant_strict,
         file = paste0(savepath,".Rdata"))
    
  }
}



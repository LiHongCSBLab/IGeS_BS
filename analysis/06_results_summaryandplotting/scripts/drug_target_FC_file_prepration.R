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


# 3) for each cell line, extract differentially expressed target genes and FC values
source("03_drug_immuneSig_enrichment/scripts/merge_p.R")
source("03_drug_immuneSig_enrichment/scripts/getsRGES.R")


mergeFC_r1 <- function(drug,
                       cancer,
                       datatype,
                       dataset,
                       lincs_cell_line_weight,
                       drugDEG_path,
                       drugMergeFCpath,
                       drugtargetMergeFCpath){
  
  
  
  drugfilepath = paste0(drugDEG_path, drug,"/")
  drugfiles = list.files(drugfilepath)
  # drugfiles = drugfiles[grep(cell, drugfiles)]
  
  
  allgeneFC = lapply(as.list(drugfiles), function(file){
    filepath = paste0(drugDEG_path, drug, "/", file)
    df = read.csv(filepath)
    df = df[order(df$X),]
    df$FC = 2^df$logFC
    df$filename = file
    file_info = data.frame(t(unlist(strsplit(gsub(".csv","",file), "_", fixed=TRUE))))
    names(file_info) = c("cell_id","treatment_time","treatment_dose")
    file_info$filename = file
    df=merge(df, file_info, by="filename")
    return(df)
  })
  names(allgeneFC) = drugfiles
  res = do.call(rbind, allgeneFC)
  names(res) = c("files", "gene", "logFC", "AveExpr","t", "P.Value","adj.P.Val",
                 "B","NES","cell_id","treatment_time","treatment_dose")
  res$drug_index = drug
  
  merged_geneFC = lapply(as.list(unique(res$gene)), function(gene){
    tmp <- res[res$gene == gene, ]
    tmp <- tmp[order(tmp$treatment_dose),]
    tmp <- tmp[!is.na(tmp$NES),]
    # print(nrow(tmp[tmp$FDR.q.val < 0.1, ]))
    
    mNES <- merge_NESmodified(
      resfile = tmp,
      lincs_cell_line_weight = lincs_cell_line_weight
    )
    mNES$gene = gene
    # merge multiple P values
    # merged_FDRp_naive <- nrow(tmp[tmp$adj.P.Val<0.05,])/nrow(tmp)
    mNES$merged_p_fisher <- merge_p_values(tmp$P.Value, method = "Fisher")
    mNES$merged_p_ACAT <- ACAT(Pvals=tmp[tmp$P.Value != 0 & tmp$P.Value != 1, ]$P.Value)
    mNES$merged_FDRp_fisher <- merge_p_values(tmp$adj.P.Val, method = "Fisher")
    mNES$merged_FDRp_ACAT <- ACAT(Pvals=tmp[tmp$adj.P.Val != 0 & tmp$adj.P.Val != 1, ]$adj.P.Val)
    
    return(mNES)
  })
  
  allgeneFC = do.call(rbind,merged_geneFC)
  allgeneFC$logFC = log2(allgeneFC$mNES)
  allgeneFC = allgeneFC[c("gene","logFC", "merged_p_ACAT", "merged_FDRp_ACAT")]  
  names(allgeneFC) = c("gene","logFC", "P.Value", "adj.P.Val")
  allgeneFC$drug = drug
  
  drugMergeFCpath = paste0(drugMergeFCpath, "/", datatype, "/", dataset, "/", cancer, "/")
  write.csv(allgeneFC, paste0(drugMergeFCpath, drug, "_allgeneFC.csv"), row.names = F, quote = F)
  return(allgeneFC)
  
}

lincs_cell_line_weight <- read.csv("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/s1_cellline_selection/TCGA_110CL/correlation_celllines_sum.csv")
names(lincs_cell_line_weight) <- c("cancertype","cell_id","median_cor","mean_cor","sd_cor","max_cor","min_cor")


# cell = 'HEPG2'

for(dataset in c("70138", "92742")){
  for(datatype in c("allgenes", "978genes")){
    
    
    # datatype = "allgenes"
    # dataset = "70138"
    drugDEG_path = "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/s2_drug_inducted_transcriptome_change/"
    
    # cancer_list <- list.files(paste0(drugDEG_path, "/", datatype, "/", dataset, "/"))
    if(dataset == "70138"){
        cancer_list <- c("LIHC","SKCM","LUAD","PAAD","PRAD","COAD","BRCA","READ")
    }else if(dataset == "92742"){
        cancer_list <- c("LIHC","SKCM","LUAD","PAAD","DLBC","LUSC","OV","STAD","UCEC","PRAD","BRCA","COAD","READ")
    }
    
    for(cancer in cancer_list){
      
      # cancer = "LIHC"

      drugDEG_path = paste0(drugDEG_path, "/", datatype, "/", dataset, "/", cancer, "/")
      druglist = list.files(drugDEG_path)
      druglistpath = paste0(drugDEG_path, druglist)
      
      
      drugMergeFCpath = "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/s2_drug_inducted_transcriptome_change/mergeFC_r1/"
      dir.create(drugMergeFCpath)
      dir.create(paste0(drugMergeFCpath, "/", datatype))
      dir.create(paste0(drugMergeFCpath, "/", datatype, "/", dataset))
      dir.create(paste0(drugMergeFCpath, "/", datatype, "/", dataset, "/", cancer, "/"))
      
      drugTargetMergeFCpath = "06_results_summaryandplotting/results/drugtargetMergeFC/"
      dir.create(drugTargetMergeFCpath)
      dir.create(paste0(drugTargetMergeFCpath, "/", datatype))
      dir.create(paste0(drugTargetMergeFCpath, "/", datatype, "/", dataset))
      dir.create(paste0(drugTargetMergeFCpath, "/", datatype, "/", dataset, "/", cancer, "/"))
      drugtargetMergeFCpath = paste0(drugMergeFCpath, "/", datatype, "/", dataset, "/", cancer, "/")
      
      allFCmerged <- lapply(as.list(druglist), function(drug){    
        allgeneFC = mergeFC_r1(drug = drug,
                               cancer = cancer, 
                               datatype = datatype, 
                               dataset = dataset, 
                               lincs_cell_line_weight = lincs_cell_line_weight,
                               drugDEG_path = drugDEG_path,
                               drugMergeFCpath = drugMergeFCpath,
                               drugtargetMergeFCpath = drugtargetMergeFCpath)
        
      })
      names(allFCmerged) <- druglist
      
      drugtargetFC <- lapply(allFCmerged, function(allgeneFC){    
      names(allgeneFC) = c("target","logFC", "P.Value", "adj.P.Val","drug_index")
      dtFC = inner_join(drugTargetMoA, allgeneFC)
      return(dtFC)
      })
      
      drugtargetFC_loose <- lapply(drugtargetFC, function(dtFC){    
        dtFC = dtFC[dtFC$adj.P.Val < 0.1, ]
        return(dtFC)
      })
      
      drugtargetFC_strict <- lapply(drugtargetFC, function(dtFC){    
        dtFC = dtFC[dtFC$adj.P.Val < 0.1, ]
        dtFC = dtFC[abs(dtFC$logFC) > 0.5 & dtFC$adj.P.Val < 0.1, ]
        return(dtFC)
      })
      
      save(drugtargetFC, 
           drugtargetFC_loose, 
           drugtargetFC_strict,
           file = paste0(drugTargetMergeFCpath, "drugMoATarget_FC.Rdata"))

    }
  }
}



# select drug-target-FC-immunesig

# dt = do.call(rbind, drugtargetFC)

# res = lapply(result_all, function(ig){
#     dt_ig = inner_join(ig, dt)
#     #dt_ig = dt_ig[-which(is.na(dt_ig$MoA)),]
# })

# names(res) = names(result_all)
# res2 = do.call(rbind,res)
# res2 = res2[res2$MoA != "",]

# res3 = res2[res2$adj.p < 0.1 & res2$ adj.P.Val <0.1,]
# rownames(res3[res3$r > 0.4 & res3$logFC>0,])
# rownames(res3[res3$r < -0.3 & res3$logFC>0,])

# Mapping drug DEGs to correlation results
# work_path <- "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/"
work_path <- "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/"

setwd(work_path)
options(stringsAsFactors = F)

library(reshape2)
library(ggplot2)
library(ggplotify)
library(grid)
library(cowplot)
library(parallel)
library(patchwork)
library(dplyr)
library(stringr)
library(data.table)
require(GGally)
library(plot3D)
library(fmsb)

source("06_results_summaryandplotting/scripts/plot_for_drug_function.R")
source("06_results_summaryandplotting/scripts/immunesig_selection_function.R")
source("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/scripts/drugITGSscore_lincs.R")
source("03_drug_immuneSig_enrichment/scripts/02_GSEA_geneset_preparation_function.R")
source("06_results_summaryandplotting/scripts/immunesig_selection_function.R")
source("06_results_summaryandplotting/scripts/ITS_analysis/ITSexpr_purity_detection.R")
source("06_results_summaryandplotting/scripts/ITS_analysis/ITS_function_detection.R")

# correlation: genes in original genesets  with tumor purity ---------------------

gene_purity_correlation <- function(cancer = "SKCM",
                                    sample_type = "Metastatic",
                                    tumor_purity_method = "TUMERIC",
                                    datatype="allgenes",
                                    num_gene = 200,
                                    r_threshold = -0.2, 
                                    workpath = "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/", 
                                    immuneSigGSpath = "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/02_tumorGene_immuneSig_correlation/data/immune_signatures/originalSignaturesGenes/",
                                    ITSpath = "03_drug_immuneSig_enrichment/data/gene_immune_sig_TUMERIC_spearman/for_GSVA/pcor/",
                                    savepath = "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/"
){
  
  cancer = "SKCM"
  sample_type = "Metastatic"
  tumor_purity_method = "TUMERIC"
  datatype="allgenes"
  num_gene = 200
  workpath = work_path
  immuneSigGSpath = "02_tumorGene_immuneSig_correlation/data/immune_signatures/originalSignaturesGenes/"
  savepath = "r1_drug_immuneSig_correlation_analysis/"
  
  
  # if gene expression has high negative correlation witg purity--------------------------------
  processedDataPath <- '02_tumorGene_immuneSig_correlation/data/TCGA_processedData/'
  processedDataPath <- paste0(processedDataPath, cancer,'/')
  TPMMatPath = paste0(processedDataPath, "mixture_file_TCGA_", cancer,"_",
                      sample_type, "_TPM.txt")
  mymatrix_filter <- read.csv(TPMMatPath, sep = '\t', row.names = 1)
  colnames(mymatrix_filter) <- gsub('\\.', '-',   colnames(mymatrix_filter))
  colnames(mymatrix_filter) <- substr(colnames(mymatrix_filter),1,15)
  
  load("01_tumor_purity/data/processedData/TCGA_tumor_purity_merged.Rdata")
  tumor_purity <- purity_merged
  tumor_purity = tumor_purity[-which(duplicated(tumor_purity$Sample_ID)),]
  
  tumor_purity_filter <- tumor_purity[is.element(tumor_purity$Sample_ID,  colnames(mymatrix_filter)),]
  tumor_purity_filter2 <- unique(tumor_purity_filter[c('Sample_ID',tumor_purity_method)])
  names(tumor_purity_filter2) <- c('Sample_ID','purity')
  
  
  library(readxl)
  library(dplyr)
  
  dir.create("06_results_summaryandplotting/immuneSigGene_found_lincs")
  dir.create(paste0("06_results_summaryandplotting/immuneSigGene_found_lincs/", purity_method))
  
  # LOAD GENES IN CCLE -------------------------------------------------------
  geneExprInCCLE <- read.csv("03_drug_immuneSig_enrichment/data/CCLE_geneexpr/geneExprInCCLE.csv")
  
  protein_coding_genes <- read.table("03_drug_immuneSig_enrichment/data/protein_coding_genes_hg38.txt", sep = '\t', header = T)
  geneExprInCCLE_filtered <- inner_join(geneExprInCCLE, protein_coding_genes)
  geneExprInCCLE_filtered <- geneExprInCCLE_filtered[geneExprInCCLE_filtered$rate < 0.2,]
  
  # LOAD GENES IN LINCS ------------------------------------------------------
  lincs_gene = read.csv("/picb/bigdata/project/FengFYM/drug_MoA_reanalysis/LINCS_data/GSE70138/file/GSE92742_Broad_LINCS_gene_info.txt", sep = "\t")
  
  # LOAD drug target GENES ---------------------------------------------------
  drugtarget <- read.csv('/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/s3_lincs_drug_MoA/drugTargetMerge.txt', sep = '\t')
  drugtarget <- data.frame(gene_name = unique(drugtarget$target))
  drugtarget_OASIS <- read.csv('/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/s3_lincs_drug_MoA/drugtarget_OASIS_TIDE.csv')
  
  # LOAD original immune signature gene sets ---------------------------------
  immunesigGS1 <- as.data.frame(read_excel(paste0(immuneSigGSpath, "immunesig_original_used.xlsx"), sheet = 1))
  immunesigGS1 <- immunesigGS1[-grep("various", immunesigGS1$Genes), ]
  immunesigGS1 <- immunesigGS1[-grep("SignatureMatrix", immunesigGS1$Genes), ]
  immunesigGS1 <- immunesigGS1[-grep("ReferenceCellExpression", immunesigGS1$Genes), ]
  immunesigGS1 <- immunesigGS1[-grep("mergedSignatures", immunesigGS1$Genes), ]
  immunesigGS1 <- immunesigGS1[-grep("L22Matrix", immunesigGS1$Genes), ]
  
  immunesigGS1 <- immunesigGS1[-c(2:5)]
  immunesigGS1Name <- as.list(immunesigGS1$immune_sig)
  immunesigGS1 <- lapply(immunesigGS1Name, function(x){
    y=immunesigGS1[immunesigGS1$immune_sig == x,][-1]
    unlist(y[-which(is.na(y))])
  })
  names(immunesigGS1) = unlist(immunesigGS1Name)
  
  immunesigGS2 <- read_excel(paste0(immuneSigGSpath, "immunesig_original_used.xlsx"), sheet = 2)
  immunesigGS2 <- as.data.frame(immunesigGS2[immunesigGS2$CancerType == cancer, ])
  immunesigGS2Name <- as.list(immunesigGS2$Signatures)
  immunesigGS2 <- lapply(immunesigGS2Name, function(x){
    y=immunesigGS2[immunesigGS2$Signatures == x,][-1]
    unlist(y[-which(is.na(y))])
  })
  names(immunesigGS2) = paste0(unlist(immunesigGS2Name),"_ConsensusTMEgsva")
  immunesigGS = c(immunesigGS1, immunesigGS2)
  immunesigGSname = data.frame(signame = names(immunesigGS), immune_sig = tolower(names(immunesigGS)))
  immunesigGSname = immunesig_name_unify(immunesigGSname)
  names(immunesigGS) = immunesigGSname$immune_sig
  
  
  immunesigGS_purity_cor <- lapply(immunesigGS, function(gs){
    # for(i in seq(length(gs))){
    # gs = immunesigGS[["cd8..t.cells_xcell"]]
    gs_filtered <- lapply(as.list(gs), function(g){
      # g=gs[i]
      if(is.element(g, rownames(mymatrix_filter))){
        gexpr = as.data.frame(t(data.frame(mymatrix_filter[rownames(mymatrix_filter) == g ,])))
        gexpr$Sample_ID = gsub('\\.', '-',   rownames(gexpr))
        df <- inner_join(gexpr, tumor_purity_filter2)
        names(df) <- c("gene","sample","purity")
        df <- df[!is.na(df$purity), ]
        
        r = cor.test(log2(1+df$gene), df$purity, method = "pearson")
        
        if(r$p.value < 0.05 & r$estimate < r_threshold) {
          g=NULL
        }else {
          g=g
        }
      }
      return(g)
    })
    gs_filtered <- unique(unlist(gs_filtered))
    if(is.null(gs_filtered)){
      return(NULL)
    }else {
      return(gs_filtered)
    }
  })
  
  res = lapply(as.list(names(immunesigGS)), function(x){
    gs1 = immunesigGS[[x]]
    gs2 = immunesigGS_purity_cor[[x]]
    res = data.frame(gene_in_sig= length(gs1), 
                     gene_cor_purity = length(gs2),
                     ratio = length(gs2)/length(gs1))
  })
  names(res) = names(immunesigGS)
  res = do.call(rbind, res)
  
  
  # function dection -----------------------------------------------------------
  print("6. conduct function analysis for each ITS.")
  # 02_tumorGene_immuneSig_correlation/scripts/correlation_dection_plot.R
  print("KEGG analysis ############################")
  
  enrichRes_kegg=list()

  for(i in 320:length(immunesigGS_purity_cor)){
    if(is.null(immunesigGS_purity_cor[[i]])){
      enrichRes_kegg[[i]] <- NULL
    }else if(immunesigGS_purity_cor[[i]] == "merged"){
      enrichRes_kegg[[i]] <- NULL
    }else {
      enrichRes_kegg[[i]] <-ITSenrichment(ITS = immunesigGS_purity_cor[[i]],
                                                 method ="kegg")
    }
  }
  names(enrichRes_kegg) = names(immunesigGS_purity_cor)
  
  print("GO-BP analysis ############################")
  enrichRes_gobp=list()
  
  for(i in 320:length(immunesigGS_purity_cor)) {
    
    if(is.null(immunesigGS_purity_cor[[i]])) {
      enrichRes_gobp[[i]] <- NULL
    }else if( is.element(immunesigGS_purity_cor[[i]], "merged")) {
      enrichRes_gobp[[i]] <- NULL
    } else {
      enrichRes_gobp[[i]] <- ITSenrichment(ITS = immunesigGS_purity_cor[[i]],
                                                  method ="gobp")
    }
  }
  
  names(enrichRes_gobp) = names(immunesigGS_purity_cor)
   
  
  dir.create(paste0(immuneSigGSpath, "/immuneSigGene_purity_correlation/"))
  save(immunesigGS, 
       immunesigGS_purity_cor, 
       res, 
       enrichRes_gobp, 
       enrichRes_kegg, 
       file = paste0(immuneSigGSpath, "/immuneSigGene_purity_correlation/",
                     cancer, "_", sample_type,".Rdata"))
  
}







# correlation: gene in ITSs with tumor purity ----------------------------------
ITSexpr_purity_detection <- function(genelist_p,
                                     genelist_n,
                                     r_threshold = -0.2,
                                     cancer = 'LIHC',
                                     sample_type = 'Primary',
                                     tumor_purity_method = "TUMERIC") {
  
  cancer = "SKCM"
  sample_type = "Metastatic"
  tumor_purity_method = "TUMERIC"
  datatype="allgenes"
  num_gene = 200
  workpath = work_path
  immuneSigGSpath = "02_tumorGene_immuneSig_correlation/data/immune_signatures/originalSignaturesGenes/"
  savepath = "r1_drug_immuneSig_correlation_analysis/"
  ITSpath = "03_drug_immuneSig_enrichment/data/gene_immune_sig_TUMERIC_spearman/for_GSVA/pcor/"
  load(paste0(ITSpath, "spearman_",cancer,"_", sample_type, "_positive_200.Rdata"))
  load(paste0(ITSpath, "spearman_",cancer,"_", sample_type, "_negative_200.Rdata"))
  
  
  
  # if gene expression has high negative correlation witg purity--------------------------------
  processedDataPath <- '02_tumorGene_immuneSig_correlation/data/TCGA_processedData/'
  processedDataPath <- paste0(processedDataPath, cancer,'/')
  TPMMatPath = paste0(processedDataPath, "mixture_file_TCGA_", cancer,"_",
                      sample_type, "_TPM.txt")
  mymatrix_filter <- read.csv(TPMMatPath, sep = '\t', row.names = 1)
  colnames(mymatrix_filter) <- gsub('\\.', '-',   colnames(mymatrix_filter))
  colnames(mymatrix_filter) <- substr(colnames(mymatrix_filter),1,15)
  
  load("01_tumor_purity/data/processedData/TCGA_tumor_purity_merged.Rdata")
  tumor_purity <- purity_merged
  tumor_purity = tumor_purity[-which(duplicated(tumor_purity$Sample_ID)),]
  
  tumor_purity_filter <- tumor_purity[is.element(tumor_purity$Sample_ID,  colnames(mymatrix_filter)),]
  tumor_purity_filter2 <- unique(tumor_purity_filter[c('Sample_ID',tumor_purity_method)])
  names(tumor_purity_filter2) <- c('Sample_ID','purity')
  
  
  # for(j in seq(length(genelist_p))){
  #   gs = genelist_p[[j]]
  
  genelist_p_filtered <- lapply(genelist_p, function(gs){
    # for(i in seq(length(gs))){
    # gs = genelist_p[[22]]
    gs_filtered <- lapply(as.list(gs), function(g){
      # g=gs[i]
      if(is.element(g, rownames(mymatrix_filter))){
        gexpr = as.data.frame(t(data.frame(mymatrix_filter[rownames(mymatrix_filter) == g ,])))
        gexpr$Sample_ID = gsub('\\.', '-',   rownames(gexpr))
        df <- inner_join(gexpr, tumor_purity_filter2)
        names(df) <- c("gene","sample","purity")
        df <- df[!is.na(df$purity), ]
        
        r = cor.test(log2(1+df$gene), df$purity, method = "pearson")
        
        if(r$p.value < 0.05 & r$estimate < r_threshold) {
          g=NULL
        }else {
          g=g
        }
      }
      return(g)
    })
    gs_filtered <- unique(unlist(gs_filtered))
    return(gs_filtered)
  })
  
  
  res_p = lapply(as.list(names(genelist_p)), function(x){
    gs1 = genelist_p[[x]]
    gs2 = genelist_p_filtered[[x]]
    res = data.frame(gene_in_sig= length(gs1), 
                     gene_cor_purity = length(gs2),
                     ratio = length(gs2)/length(gs1))
  })
  
  names(res_p) = names(genelist_p)
  res_p = do.call(rbind, res_p)
  
  
  # function dection -----------------------------------------------------------
  print("6. conduct function analysis for each ITS.")
  # 02_tumorGene_immuneSig_correlation/scripts/correlation_dection_plot.R
  print("KEGG analysis ############################")
  
  enrichRes_p_kegg=list()
  for(i in seq(length(genelist_p_filtered))){
    enrichRes_p_kegg[[i]] <-ITSenrichment(ITS = genelist_p_filtered[[i]],
                                                 method ="kegg")
  }
  
  names(enrichRes_p_kegg) = names(genelist_p_filtered)
  
  print("GO-BP analysis ############################")
  enrichRes_p_gobp=list()
  for(i in seq(length(genelist_p_filtered))){
    enrichRes_p_gobp[[i]] <-ITSenrichment(ITS = genelist_p_filtered[[i]],
                                                 method ="gobp")
  }
  
  names(enrichRes_p_gobp) = names(genelist_p_filtered)
  #   save(enrichRes, 
  #        file= paste0("06_results_summaryandplotting/ITS_function/gobp_", 
  #                     cancer, "_", sample_type,"_p.Rdata"))
  
  
  
  dir.create(paste0(immuneSigGSpath, "/ITSGene_purity_correlation/"))
  save(genelist_p, genelist_p_filtered, res_p, enrichRes_p_gobp, enrichRes_p_kegg, 
       file = paste0(immuneSigGSpath, "/ITSGene_purity_correlation/",
                    cancer, "_", sample_type,"_p.Rdata"))
  
  
  
  genelist_n_filtered <- lapply(genelist_n, function(gs){
    # for(i in seq(length(gs))){
    
    gs_filtered <- lapply(as.list(gs), function(g){
      # g=gs[i]
      if(is.element(g, rownames(mymatrix_filter))){
        gexpr = as.data.frame(t(data.frame(mymatrix_filter[rownames(mymatrix_filter) == g ,])))
        gexpr$Sample_ID = gsub('\\.', '-',   rownames(gexpr))
        df <- inner_join(gexpr, tumor_purity_filter2)
        names(df) <- c("gene","sample","purity")
        df <- df[!is.na(df$purity), ]
        
        r = cor.test(log2(1+df$gene), df$purity, method = "pearson")
        
        if(r$p.value < 0.05 & r$estimate < r_threshold) {
          g=NULL
        }else {
          g=g
        }
      }
      return(g)
    })
    gs_filtered <- unique(unlist(gs_filtered))
    return(gs_filtered)
  })
    
    
    res_n = lapply(as.list(names(genelist_n)), function(x){
      gs1 = genelist_n[[x]]
      gs2 = genelist_n_filtered[[x]]
      res = data.frame(gene_in_sig= length(gs1), 
                       gene_cor_purity = length(gs2),
                       ratio = length(gs2)/length(gs1))
    })
    names(res_n) = names(genelist_n)
    res_n = do.call(rbind, res_n)
    
    
    # function dection -----------------------------------------------------------
    print("6. conduct function analysis for each ITS.")
    
    print("KEGG analysis ############################")
    
    enrichRes_n_kegg=list()
    for(i in seq(length(genelist_n_filtered))){
      enrichRes_n_kegg[[i]] <- ITSenrichment(ITS = genelist_n_filtered[[i]],
                                                    method ="kegg")
    }
    names(enrichRes_n_kegg) = names(genelist_n_filtered)
    
    
    print("GO-BP analysis ############################")
    
    # 02_tumorGene_immuneSig_correlation/scripts/correlation_dection_plot.R
    enrichRes_n_gobp=list()
    for(i in seq(length(genelist_n_filtered))){
      enrichRes_n_gobp[[i]] <-ITSenrichment(ITS = genelist_n_filtered[[i]],
                                                   method ="gobp")
    }
    names(enrichRes_n_gobp) = names(genelist_n_filtered)
    #   save(enrichRes, 
    #        file= paste0("06_results_summaryandplotting/ITS_function/gobp_", 
    #                     cancer, "_", sample_type,"_p.Rdata"))
    
    
    dir.create(immuneSigGSpath, "/ITSGene_purity_correlation/")
    save(genelist_n, genelist_n_filtered, res_n, enrichRes_n_gobp, enrichRes_n_kegg, 
         file= paste0(immuneSigGSpath, "/ITSGene_purity_correlation/",
                      cancer, "_", sample_type,"_n.Rdata"))
    
    
  
  genelist_filter <- list(genelist_p_filtered, genelist_n_filtered)
  return(genelist_filter)
  
}


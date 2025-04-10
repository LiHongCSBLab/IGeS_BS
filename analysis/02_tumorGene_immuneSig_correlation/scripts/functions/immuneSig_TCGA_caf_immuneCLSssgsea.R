# raw count normalization for pathway enrichment analysis
# code written by Ping Lin, modified by FYM Feng
## rawCount: raw count of your RNAseq expression data
## method: select analysis method, currently support "GSEA" and "GSVA",

nCountGenerator <- function(rawCount,
                            method = c("GSVA", "GSEA"),
                            output = NULL) {
  
  
  if (!require(DESeq2)) {
    BiocManager::install("DESeq2")
  }  
  
  
  sf <- estimateSizeFactorsForMatrix(rawCount)
  nCounts <- t(apply(rawCount, 1, function(x){x / sf}))
  
  ## filter out genes with low expression level
  data <- nCounts[apply(nCounts, 1, function(x){sum(x>1)>length(x)*0.9}),]
  if(method == "GSVA"){
    
    return(data)
    
  }else if(method == "GSEA"){
    
    if(is.null(output)){
      
      print(" ERROR: prefix of output must be provided! ")
      
    }else{
      
      ## data: your normalized expression data (no log-scaled data)
      ## Please make sure the rownames of your data are the same as the gene names in the gene set you provided to GSEA
      logExpr <- log2(data + 1)
      logExpr.gsea <- as.data.frame(logExpr)
      
      ## add NAME and DESCRIPTION
      logExpr.gsea <- cbind(rownames(logExpr.gsea),
                            rownames(logExpr.gsea),
                            logExpr.gsea)
      colnames(logExpr.gsea) <- c("NAME", "DESCRIPTION", colnames(logExpr))
      
      ## output gsea expr data
      write.table(logExpr.gsea,
                  file = paste0(output, "logExpr.gsea.txt"),
                  sep = "\t", row.names = F, col.names = T)
      
      ## clsFile
      clsFile <- paste0(output, "grpInfo.cls")
      numGrp <- 2
      yourGrpOr <- c("R", "NR")
      cat(ncol(logExpr), "2 1\n", file = clsFile, append = F)
      cat("#", yourGrpOr, file = clsFile, append = T)
      cat("\n", file = clsFile, append = T)
      cat(paste(colnames(logExpr), collapse = " "), "\n",
          file = clsFile, append = T)
      
    }
  }
}



# original code conduct a median scale in python 
#   def median_scale(data, clip=None):
#     c_data = (data - data.median()) / data.mad()
#     if clip is not None:
#         return c_data.clip(-clip, clip)
#     return c_data

median_scale <- function(x){
  return(t(apply(x, 1, function(x) (x - median(na.omit(x)))/mad(na.omit(x)))))
}

#  Gene Enrichment Analysis with `GSVA` package for immune signatures
## data: output nCount matrix with `nCountGenerator(..., method = "GSVA", ...),
##       RNAseq raw counts data -> "Possion"
## cancer: TCGA cancer name abbreivation
## method: select method for enrichement analysis, support "gsva" and "ssgsea"
## geneset: several pathways or gene sets have been provided:
##          kegg: pahtway from KEGG
##          bgobp: Biological process from Gene Ontology
##          hallmark: cancer hallmark pathway
##          reactome: pathway from Reactome database
## output: prefix of output

TCGA_caf_immuneCLS_consistency <- function(cancer,
                                           sample_type,
                                           method = c("gsva", "ssgsea"),
                                           kcdf = c("Poisson", "Gaussian"), 
                                           min.sz=1, 
                                           max.sz=1000,
                                           sig_dir,
                                           processedDataPath,
                                           originalSigPath, 
                                           savePath
) {
  
  ## load R package
  if (!require(GSVA)) {
    BiocManager::install("GSVA")
  }
  library(stringr)
  library(ggplot2)
  library(plyr)
  library(patchwork)
  library(gridExtra)
  
  ## please make sure that the rownames of data
  ## are the same with the gene names in the gene set
  ## that you want to run GSVA
  # rownames(data) <- .....
  
  ## load ImmuneSignature --------------------------------------------------------------
  #   library(readxl)
  #   sig_dir = "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/02_tumorGene_immuneSig_correlation/data/immune_signatures/"
  #   immuneSigGeneSet <- read_excel(paste0(sig_dir, "originalSignaturesGenes/immunesig_original.xlsx"), sheet = 1)
  #   info <- as.data.frame(immuneSigGeneSet[c(1:5)])
  #   immuneSigGeneSet <- as.data.frame(immuneSigGeneSet[-c(2:5)])
  #   immuneSigGeneSetList <- list()
  #   for(i in seq(nrow(immuneSigGeneSet))){
  #       g = immuneSigGeneSet[i,][-1]
  #       immuneSigGeneSetList[[i]] = g[!is.na(g)]
  #   }
  #   names(immuneSigGeneSetList) = immuneSigGeneSet[,1]
  #   save(immuneSigGeneSetList, file = paste0(sig_dir, "originalSignaturesGenes/immuneSigGeneSet_all.Rdata"))
  
  #   gs1 = immuneSigGeneSetList[grep("TCGA_caf_immuneCLS", names(immuneSigGeneSetList))]
  #   save(gs1, file = paste0(sig_dir, "originalSignaturesGenes/immuneSigGeneSetList_TCGA_caf_immuneCLS.Rdata"))
  #   otherImmuneSigGeneSet = info[c(grep("Senbabaoglu", info$OriginalPublications),
  #                                  grep("Attractors", info$OriginalPublications),
  #                                  grep("ICR", info$OriginalPublications),
  #                                  grep("Bindea", info$OriginalPublications)), ]$Signatures
  
  #   gs2 = immuneSigGeneSetList[which(is.element(names(immuneSigGeneSetList), otherImmuneSigGeneSet))]
  #   save(gs2, file = paste0(sig_dir, "originalSignaturesGenes/immuneSigGeneSetList_sig160_othersigs.Rdata"))
  
  
  #   sig_dir = "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/02_tumorGene_immuneSig_correlation/data/immune_signatures/"
  
  load(paste0(sig_dir, "originalSignaturesGenes/immuneSigGeneSetList_TCGA_caf_immuneCLS.Rdata"))
  
  
  # cancer = 'LIHC'
  cancer_type = paste0('TCGA-', cancer)
  
  processedDataPath <- paste0(processedDataPath, cancer,'/')
  dir.create(processedDataPath)
  
  dir.create(paste0(savePath, "/TCGA_caf_immuneCLS_result/"))
  dir.create(paste0(savePath, "/TCGA_caf_immuneCLS_result/raw_output/"))
  dir.create(paste0(savePath, "/TCGA_caf_immuneCLS_result/scaled_output/"))
  dir.create(paste0(savePath, "/TCGA_caf_immuneCLS_result/consistency_plot/"))
  dir.create(paste0(savePath, "/TCGA_caf_immuneCLS_result/consistency_plot_scaled/"))
  
  
  # Prepare TCGA expression matrix --------------------------------------------
  
  # load in full expression matrix
  # drop out low expression genes
  # save to certain path
  
  TPMMatPath = paste0(processedDataPath, "mixture_file_TCGA_", cancer,"_",sample_type, "_TPM.txt")
  TPM_filtered <- read.csv(TPMMatPath, sep = '\t', row.names = 1)
  colnames(TPM_filtered) <- gsub('\\.', '-',   colnames(TPM_filtered))
  data = apply(TPM_filtered, 2, function(a) log2(a+1))
  
  ## gsva
  if(method == "gsva"){
    
    esOut <- gsva(data,
                  gs = gs1, 
                  method="gsva",
                  mx.diff=FALSE, kcdf=kcdf,
                  verbose=FALSE, parallel.sz=1,
                  min.sz=min.sz, max.sz=max.sz,
                  ssgsea.norm=FALSE)
    
  }else if (method == 'ssgsea'){
    
    esOut <- gsva(data,
                  gs=gs1, 
                  method="ssgsea",
                  mx.diff=FALSE, kcdf=kcdf,
                  verbose=FALSE, parallel.sz=1,
                  min.sz=min.sz, max.sz=max.sz,
                  ssgsea.norm=TRUE)
    
  }
  save(esOut, 
       file = paste0(savePath,"/TCGA_caf_immuneCLS_result/raw_output/", method, "_", cancer,"_",sample_type,"_immueSig.Rdata"))
  
  # original code conduct a median scale in python 
  #   def median_scale(data, clip=None):
  #     c_data = (data - data.median()) / data.mad()
  #     if clip is not None:
  #         return c_data.clip(-clip, clip)
  #     return c_data
  
  #   median_scale <- function(x){
  #     return(t(apply(x, 1, function(x) (x - median(na.omit(x)))/mad(na.omit(x)))))
  #   }
  
  
  esOut_scaled <- median_scale(esOut)
  save(esOut_scaled, 
       file = paste0(savePath,"/TCGA_caf_immuneCLS_result/scaled_output/", method, "_", cancer,"_",sample_type,"_immueSig.Rdata"))
  
  
  #   originalSigPath <- "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/02_tumorGene_immuneSig_correlation/data/immune_signatures/"
  precompute_sig <- read.csv(paste0(originalSigPath, "cluster_signatures_cancercell/signatures-tcga.tsv"), header = T, sep = "\t", row.names=1)
  x <- str_replace_all(colnames(precompute_sig), pattern='\\.', replacement='-')
  colnames(precompute_sig) <- x
  
  sig = esOut
  colnames(sig) <- substr(colnames(sig),1,12)
  rownames(sig) <- gsub("_TCGA_caf_immuneCLS", "", rownames(sig))
  sum(colnames(sig) %in% colnames(precompute_sig))
  # 
  # Then we need to get these into a comparible form #
  #
  
  sharedIDs <- intersect(colnames(precompute_sig), colnames(sig))
  genesets  <- intersect(rownames(precompute_sig), rownames(sig))
  length(sharedIDs)
  length(genesets)
  
  if(length(sharedIDs) > 0 ){
    
    sig2 <- sig[genesets,sharedIDs]
    precompute_sig2 <- precompute_sig[genesets,sharedIDs]
    sig3 = as.data.frame(t(sig2))
    sig3$type = "recompute"
    precompute_sig3 = as.data.frame(t(precompute_sig2))
    precompute_sig3$type = "downloaded"
    df = rbind(sig3, precompute_sig3)
    
    b <- median_scale(sig2)
    b2 = as.data.frame(t(b))
    b2$type = "recompute"
    precompute_sig3 = as.data.frame(t(precompute_sig2))
    precompute_sig3$type = "downloaded"
    df2 = rbind(b2, precompute_sig3)
    
    
    # and replot #
    p1=list()
    p2=list()
    
    
    for(i in 1:nrow(sig2)){
      
      cor1s = cor.test(x=as.numeric(sig2[i,]), 
                       y=as.numeric(precompute_sig2[i,]), 
                       method = "spearman")
      
      p1_1 = qplot(x=as.numeric(sig2[i,]), y=as.numeric(precompute_sig2[i,])) + 
        geom_abline(intercept = 0, slope = 1) +
        theme_bw()+
        ggtitle(paste0(rownames(sig2)[i], 
                       "\n spearman's r = ", round(cor1s$estimate,3), 
                       "\n P value = ", signif(cor1s$p.value,3)))
      
      df1_2 = df[c(colnames(df)[i], "type")]
      names(df1_2) = c('immunesig', 'type')
      mu1 <- ddply(df1_2, "type", summarise, grp.mean=mean(immunesig))
      p1_2 <- ggplot(df1_2, aes(x=immunesig, color=type)) +
        geom_density()+
        geom_vline(data=mu1, aes(xintercept=grp.mean, color=type),
                   linetype="dashed")+
        theme_bw()
      
      
      
      cor1_2s = cor.test(x=as.numeric(b[i,]), 
                         y=as.numeric(precompute_sig2[i,]), 
                         method = "spearman")
      
      p2_1 = qplot(x=as.numeric(b[i,]), y=as.numeric(precompute_sig2[i,])) + 
        geom_abline(intercept = 0, slope = 1) +
        theme_bw()+
        ggtitle(paste0(rownames(sig2)[i], "\n spearman's r = ", round(cor1_2s$estimate,3), 
                       "\n P value = ", signif(cor1_2s$p.value,3)))      
      
      df2_2 = df2[c(colnames(df2)[i], "type")]
      names(df2_2) = c('immunesig', 'type')
      mu2 <- ddply(df2_2, "type", summarise, grp.mean=mean(immunesig))
      p2_2 <- ggplot(df2_2, aes(x=immunesig, color=type)) +
        geom_density()+
        geom_vline(data=mu2, aes(xintercept=grp.mean, color=type),
                   linetype="dashed")+
        theme_bw()
      
      p1[[i]] = p1_1 / p1_2
      p2[[i]] = p2_1 / p2_2
      
    }
    
    
    pdf(paste0(savePath, "/TCGA_caf_immuneCLS_result/consistency_plot/", 
               cancer, "_", sample_type,".pdf"),
        width = 5,
        height = 8, 
        onefile = TRUE)
    for(i in seq(length(p1))){
      print(p1[[i]])
    }
    dev.off()
    
    
    pdf(paste0(savePath, "/TCGA_caf_immuneCLS_result/consistency_plot_scaled/", 
               cancer, "_", sample_type,"_scaled.pdf"), 
        width = 5,
        height = 8, 
        onefile = TRUE)
    
    for(i in seq(length(p2))){
      print(p2[[i]])
    }
    dev.off()
    
  }  
}




TCGA_caf_immuneCLS <- function(data, 
                               cancer,
                               sample_type,
                               method = c("gsva", "ssgsea"),
                               kcdf = c("Poisson", "Gaussian"), 
                               min.sz=1, 
                               max.sz=1000,
                               sig_dir,
                               savePath) {
  
  ## load R package
  if (!require(GSVA)) {
    BiocManager::install("GSVA")
  }
  library(stringr)
  library(ggplot2)
  library(plyr)
  library(patchwork)
  library(gridExtra)
  
  ## please make sure that the rownames of data
  ## are the same with the gene names in the gene set
  ## that you want to run GSVA
  # rownames(data) <- .....
  
  
  #   sig_dir = "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/02_tumorGene_immuneSig_correlation/data/immune_signatures/"
  
  load(paste0(sig_dir, "originalSignaturesGenes/immuneSigGeneSetList_TCGA_caf_immuneCLS.Rdata"))
  
  
  # cancer = 'LIHC'
  
  dir.create(paste0(savePath, "/TCGA_caf_immuneCLS_result/"))
  dir.create(paste0(savePath, "/TCGA_caf_immuneCLS_result/raw_output/"))
  dir.create(paste0(savePath, "/TCGA_caf_immuneCLS_result/scaled_output/"))
  
  
  # Prepare TCGA expression matrix --------------------------------------------
  
  # load in full expression matrix
  # drop out low expression genes
  # save to certain path
  
  data = apply(data, 2, function(a) log2(a+1))
  
  ## gsva
  if(method == "gsva"){
    
    esOut <- gsva(data,
                  gs = gs1, 
                  method="gsva",
                  mx.diff=FALSE, kcdf=kcdf,
                  verbose=FALSE, parallel.sz=1,
                  min.sz=min.sz, max.sz=max.sz,
                  ssgsea.norm=FALSE)
    
  }else if (method == 'ssgsea'){
    
    esOut <- gsva(data,
                  gs=gs1, 
                  method="ssgsea",
                  mx.diff=FALSE, kcdf=kcdf,
                  verbose=FALSE, parallel.sz=1,
                  min.sz=min.sz, max.sz=max.sz,
                  ssgsea.norm=TRUE)
    
  }
  save(esOut, 
       file = paste0(savePath,"/TCGA_caf_immuneCLS_result/raw_output/", method, "_", cancer,"_",sample_type,"_immueSig.Rdata"))
  
  
  esOut_scaled <- median_scale(esOut)
  save(esOut_scaled, 
       file = paste0(savePath,"/TCGA_caf_immuneCLS_result/scaled_output/", method, "_", cancer,"_",sample_type,"_immueSig.Rdata"))
}

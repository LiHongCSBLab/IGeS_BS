zscore.rows2<-function(x){
  return(t(apply(x, 1, function(x) (x - median(na.omit(x)))/sd(na.omit(x)))))
}


immuneSig_sig160_consistency <- function(cancer,
                                         sample_type,
                                         processedDataPath,
                                         sourceCodePath,
                                         originalSigPath, 
                                         savePath){
  
  library(stringr)
  library(ggplot2)
  library(plyr)
  library(patchwork)
  
  # cancer = 'LIHC'
  cancer_type = paste0('TCGA-', cancer)
  
  processedDataPath <- paste0(processedDataPath, cancer,'/')
  dir.create(processedDataPath)
  
  dir.create(paste0(savePath, "/sig160_result/"))
  dir.create(paste0(savePath, "/sig160_result/raw_output/"))
  dir.create(paste0(savePath, "/sig160_result/zscore_output/"))
  dir.create(paste0(savePath, "/sig160_result/consistency_plot/"))
  dir.create(paste0(savePath, "/sig160_result/consistency_plot_zscore/"))
  # Prepare TCGA expression matrix --------------------------------------------
  
  # load in full expression matrix
  # drop out low expression genes
  # save to certain path
  
  TPMMatPath = paste0(processedDataPath, "mixture_file_TCGA_", cancer,"_",sample_type, "_TPM.txt")
  TPM_filtered <- read.csv(TPMMatPath, sep = '\t', row.names = 1)
  colnames(TPM_filtered) <- gsub('\\.', '-',   colnames(TPM_filtered))
  
  
  # apply the transforms #
  # first the log2
  datSubsetTransformed <- apply(TPM_filtered, 2, function(a) log2(a+1))
  
  dim(datSubsetTransformed)
  
  # then the median scale
  datSubsetTransformed <- t(apply(datSubsetTransformed, 1, function(a) a - median(a, na.rm=T)))
  
  dim(datSubsetTransformed)
  
  load(paste0(sourceCodePath, 'Immune-Subtype-Clustering-master/shiny-app/Immune-Subtype-Clustering/data/comparative_immuneSigs_geneLists4.rda'))
  
  sigs1_2_geneIDs2<-as.character(na.omit(sigs1_2_eg2[[1]]$probe))
  for(i in 2:length(sigs1_2_eg2)){
    sigs1_2_geneIDs2<-c(sigs1_2_geneIDs2,as.character(na.omit(sigs1_2_eg2[[i]]$probe)))
  }
  sigs1_2_geneIDs2<-unique(sigs1_2_geneIDs2)   ## 2652 unique
  length(sigs1_2_geneIDs2)
  
  
  source(paste0(sourceCodePath, 'Immune-Subtype-Clustering-master/Scripts/ImmuneSigs68_function.R'))
  scores <- ImmuneSigs_function(datSubsetTransformed, 
                                sigs1_2_eg2,
                                sigs12_weighted_means,
                                sigs12_module_weights,
                                sigs1_2_names2,
                                sigs1_2_type2)
  
  write.csv(scores, paste0(savePath, "/sig160_result/raw_output/", cancer, "_", sample_type, ".csv"),
            row.names = T, quote = F)
  scores_z <- zscore.rows2(scores)
  write.csv(scores_z, paste0(savePath, "/sig160_result/zscore_output/", cancer, "_", sample_type, ".csv"),
            row.names = T, quote = F)
  
  reportedScores <- read.table(paste0(originalSigPath, 'Scores_160_Signatures.tsv'), header=T, sep='\t', stringsAsFactors = F)
  
  x <- str_replace_all(colnames(reportedScores), pattern='\\.', replacement='-')
  colnames(reportedScores) <- x
  
  # did we get them matched?
  sum(colnames(scores) %in% colnames(reportedScores))
  # 
  # Then we need to get these into a comparible form #
  #
  
  sharedIDs <- intersect(colnames(reportedScores), colnames(scores))
  genesets  <- intersect(reportedScores$SetName, rownames(scores))
  length(sharedIDs)
  length(genesets)
  if(length(sharedIDs) > 0 ){
    
    scores2 <- scores[genesets,sharedIDs]
    reportedScores2 <- reportedScores[genesets,sharedIDs]
    
    #   dim(scores2)
    #   dim(reportedScores2)
    
    #   scores2[1:5,1:5]  
    #   reportedScores2[1:5,1:5]
    
    #   diag(cor(t(scores2), t(reportedScores2), method = "spearman"))
    #   diag(cor(t(scores2), t(reportedScores2), method = "pearson"))
    
    
    representativeSig = c("CSF1_response",
                          "LIexpression_score",
                          "TGFB_score_21050467",
                          "Module3_IFN_score",
                          "CHANG_CORE_SERUM_RESPONSE_UP")
    sample <- which(is.element(rownames(scores2), representativeSig))
    
    # comparing the raw computed scores with reported scores #
    cor1s = cor.test(x=as.numeric(scores2[sample[1],]), 
                     y=as.numeric(reportedScores2[sample[1],]), 
                     method = "spearman")
    
    p1_1 = qplot(x=as.numeric(scores2[sample[1],]), y=as.numeric(reportedScores2[sample[1],])) + 
      geom_abline(intercept = 0, slope = 1)+
      theme_bw() +
      ggtitle(paste0(rownames(scores2)[sample[1]], "\n spearman's r = ", round(cor1s$estimate,3), "\n P value = ", signif(cor1s$p.value,3)))
    
    cor2s = cor.test(x=as.numeric(scores2[sample[2],]), 
                     y=as.numeric(reportedScores2[sample[2],]), 
                     method = "spearman")
    p2_1 = qplot(x = as.numeric(scores2[sample[2], ]), y = as.numeric(reportedScores2[sample[2], ])) +
      geom_abline(intercept = 0, slope = 1) +
      theme_bw()+      
      ggtitle(label = paste0(rownames(scores2)[sample[2]], "\n spearman's r = ", round(cor2s$estimate,3), "\n P value = ", signif(cor2s$p.value,3)))
    
    cor3s = cor.test(x=as.numeric(scores2[sample[3],]), 
                     y=as.numeric(reportedScores2[sample[3],]), 
                     method = "spearman")
    p3_1 = qplot(x=as.numeric(scores2[sample[3],]), y=as.numeric(reportedScores2[sample[3],])) + 
      geom_abline(intercept = 0, slope = 1) +
      theme_bw()+
      ggtitle(label = paste0(rownames(scores2)[sample[3]], "\n spearman's r = ", round(cor3s$estimate,3), "\n P value = ", signif(cor3s$p.value,3)))
    
    cor4s = cor.test(x=as.numeric(scores2[sample[4],]), 
                     y=as.numeric(reportedScores2[sample[4],]), 
                     method = "spearman")
    p4_1 = qplot(x=as.numeric(scores2[sample[4],]), y=as.numeric(reportedScores2[sample[4],])) + 
      geom_abline(intercept = 0, slope = 1) +
      theme_bw()+
      ggtitle(label = paste0(rownames(scores2)[sample[4]], "\n spearman's r = ", round(cor4s$estimate,3), "\n P value = ", signif(cor4s$p.value,3)))
    
    cor5s = cor.test(x=as.numeric(scores2[sample[5],]), 
                     y=as.numeric(reportedScores2[sample[5],]), 
                     method = "spearman")
    p5_1 = qplot(x=as.numeric(scores2[sample[5],]), y=as.numeric(reportedScores2[sample[5],])) + 
      geom_abline(intercept = 0, slope = 1) +
      theme_bw()+
      ggtitle(label = paste0(rownames(scores2)[sample[5]], "\n spearman's r = ", round(cor5s$estimate,3), "\n P value = ", signif(cor5s$p.value,3)))
    
    
    scores3 = as.data.frame(t(scores2))
    scores3$type = "recompute"
    reportedScores3 = as.data.frame(t(reportedScores2))
    reportedScores3$type = "downloaded"
    df = rbind(scores3, reportedScores3)
    
    mu1 <- ddply(df, "type", summarise, grp.mean=mean(CSF1_response))
    p1_3 <- ggplot(df, aes(x=CSF1_response, color=type)) +
      geom_density()+
      geom_vline(data=mu1, aes(xintercept=grp.mean, color=type),
                 linetype="dashed")+
      theme_bw()
    
    mu2 <- ddply(df, "type", summarise, grp.mean=mean(LIexpression_score))
    p2_3  <- ggplot(df, aes(x=LIexpression_score, color=type)) +
      geom_density()+
      geom_vline(data=mu2, aes(xintercept=grp.mean, color=type),
                 linetype="dashed")+
      theme_bw()
    
    mu3 <- ddply(df, "type", summarise, grp.mean=mean(TGFB_score_21050467))
    p3_3  <- ggplot(df, aes(x=TGFB_score_21050467, color=type)) +
      geom_density()+
      geom_vline(data=mu3, aes(xintercept=grp.mean, color=type),
                 linetype="dashed")+
      theme_bw()
    
    mu4 <- ddply(df, "type", summarise, grp.mean=mean(Module3_IFN_score))
    p4_3  <- ggplot(df, aes(x=Module3_IFN_score, color=type)) +
      geom_density()+
      geom_vline(data=mu4, aes(xintercept=grp.mean, color=type),
                 linetype="dashed")+
      theme_bw()
    
    mu5 <- ddply(df, "type", summarise, grp.mean=mean(CHANG_CORE_SERUM_RESPONSE_UP))
    p5_3  <- ggplot(df, aes(x=CHANG_CORE_SERUM_RESPONSE_UP, color=type)) +
      geom_density()+
      geom_vline(data=mu5, aes(xintercept=grp.mean, color=type),
                 linetype="dashed")+
      theme_bw()
    
    p1 = p1_1 | p2_1 | p3_1  | p4_1 | p5_1
    p3 = p1_3 | p2_3 | p3_3  | p4_3 | p5_3
    
    p = p1/p3 + 
      plot_layout(widths = c(2, 1), heights = unit(c(18, 1), c('cm', 'null')))
    ggsave(paste0(savePath, "/sig160_result/consistency_plot/", cancer, "_", sample_type, ".pdf"),
           p, width = 40, height = 15)
    

      # now let's apply the sample-level normalization #
    b <- zscore.rows2(scores2)
    # and replot #
    # p1_2 = qplot(x=as.numeric(b[sample[1],]), y=as.numeric(reportedScores2[sample[1],])) + geom_abline()
    # p2_2 = qplot(x=as.numeric(b[sample[2],]), y=as.numeric(reportedScores2[sample[2],])) + geom_abline()
    # p3_2 = qplot(x=as.numeric(b[sample[3],]), y=as.numeric(reportedScores2[sample[3],])) + geom_abline()
    # p4_2 = qplot(x=as.numeric(b[sample[4],]), y=as.numeric(reportedScores2[sample[4],])) + geom_abline()
    # p5_2 = qplot(x=as.numeric(b[sample[5],]), y=as.numeric(reportedScores2[sample[5],])) + geom_abline()

        # comparing the raw computed scores with reported scores #
    cor1s = cor.test(x=as.numeric(b[sample[1],]), 
                     y=as.numeric(reportedScores2[sample[1],]), 
                     method = "spearman")
    
    p1_2 = qplot(x=as.numeric(b[sample[1],]), y=as.numeric(reportedScores2[sample[1],])) + 
      geom_abline(intercept = 0, slope = 1)+
      theme_bw() +
      ggtitle(paste0(rownames(b)[sample[1]], "\n spearman's r = ", round(cor1s$estimate,3), "\n P value = ", signif(cor1s$p.value,3)))
    
    cor2s = cor.test(x=as.numeric(b[sample[2],]), 
                     y=as.numeric(reportedScores2[sample[2],]), 
                     method = "spearman")
    p2_2 = qplot(x = as.numeric(b[sample[2], ]), y = as.numeric(reportedScores2[sample[2], ])) +
      geom_abline(intercept = 0, slope = 1) +
      theme_bw()+      
      ggtitle(label = paste0(rownames(b)[sample[2]], "\n spearman's r = ", round(cor2s$estimate,3), "\n P value = ", signif(cor2s$p.value,3)))
    
    cor3s = cor.test(x=as.numeric(b[sample[3],]), 
                     y=as.numeric(reportedScores2[sample[3],]), 
                     method = "spearman")
    p3_2 = qplot(x=as.numeric(b[sample[3],]), y=as.numeric(reportedScores2[sample[3],])) + 
      geom_abline(intercept = 0, slope = 1) +
      theme_bw()+
      ggtitle(label = paste0(rownames(b)[sample[3]], "\n spearman's r = ", round(cor3s$estimate,3), "\n P value = ", signif(cor3s$p.value,3)))
    
    cor4s = cor.test(x=as.numeric(b[sample[4],]), 
                     y=as.numeric(reportedScores2[sample[4],]), 
                     method = "spearman")
    p4_2 = qplot(x=as.numeric(b[sample[4],]), y=as.numeric(reportedScores2[sample[4],])) + 
      geom_abline(intercept = 0, slope = 1) +
      theme_bw()+
      ggtitle(label = paste0(rownames(b)[sample[4]], "\n spearman's r = ", round(cor4s$estimate,3), "\n P value = ", signif(cor4s$p.value,3)))
    
    cor5s = cor.test(x=as.numeric(b[sample[5],]), 
                     y=as.numeric(reportedScores2[sample[5],]), 
                     method = "spearman")
    p5_2 = qplot(x=as.numeric(b[sample[5],]), y=as.numeric(reportedScores2[sample[5],])) + 
      geom_abline(intercept = 0, slope = 1) +
      theme_bw()+
      ggtitle(label = paste0(rownames(b)[sample[5]], "\n spearman's r = ", round(cor5s$estimate,3), "\n P value = ", signif(cor5s$p.value,3)))

    b3 = as.data.frame(t(b))
    b3$type = "recompute"
    reportedScores3 = as.data.frame(t(reportedScores2))
    reportedScores3$type = "downloaded"
    df = rbind(b3, reportedScores3)
    
    mu1 <- ddply(df, "type", summarise, grp.mean=mean(CSF1_response))
    p1_4 <- ggplot(df, aes(x=CSF1_response, color=type)) +
      geom_density()+
      geom_vline(data=mu1, aes(xintercept=grp.mean, color=type),
                 linetype="dashed")+
      theme_bw()
    
    mu2 <- ddply(df, "type", summarise, grp.mean=mean(LIexpression_score))
    p2_4  <- ggplot(df, aes(x=LIexpression_score, color=type)) +
      geom_density()+
      geom_vline(data=mu2, aes(xintercept=grp.mean, color=type),
                 linetype="dashed")+
      theme_bw()
    
    mu3 <- ddply(df, "type", summarise, grp.mean=mean(TGFB_score_21050467))
    p3_4  <- ggplot(df, aes(x=TGFB_score_21050467, color=type)) +
      geom_density()+
      geom_vline(data=mu3, aes(xintercept=grp.mean, color=type),
                 linetype="dashed")+
      theme_bw()
    
    mu4 <- ddply(df, "type", summarise, grp.mean=mean(Module3_IFN_score))
    p4_4  <- ggplot(df, aes(x=Module3_IFN_score, color=type)) +
      geom_density()+
      geom_vline(data=mu4, aes(xintercept=grp.mean, color=type),
                 linetype="dashed")+
      theme_bw()
    
    mu5 <- ddply(df, "type", summarise, grp.mean=mean(CHANG_CORE_SERUM_RESPONSE_UP))
    p5_4  <- ggplot(df, aes(x=CHANG_CORE_SERUM_RESPONSE_UP, color=type)) +
      geom_density()+
      geom_vline(data=mu5, aes(xintercept=grp.mean, color=type),
                 linetype="dashed")+
      theme_bw()

    p2 = p1_2 | p2_2 | p3_2  | p4_2 | p5_2
    p4 = p1_4 | p2_4 | p3_4  | p4_4 | p5_4
    p = p2/p4 + 
      plot_layout(widths = c(2, 1), heights = unit(c(18, 1), c('cm', 'null')))
    ggsave(paste0(savePath, "/sig160_result/consistency_plot_zscore/", cancer, "_", sample_type, ".pdf"),
           p, width = 40, height = 15)
  }
}






immuneSig_sig160 <- function(data,
                             cancer,
                             sample_type,
                             sourceCodePath,
                             savePath){
  
  library(stringr)
  library(ggplot2)
  library(plyr)
  library(patchwork)
  
  
  dir.create(paste0(savePath, "/sig160_result/"))
  dir.create(paste0(savePath, "/sig160_result/raw_output/"))
  dir.create(paste0(savePath, "/sig160_result/zscore_output/"))
  
  # Prepare gene expression matrix --------------------------------------------
  
  # apply the transforms 
  ## first the log2
  datSubsetTransformed <- apply(data, 2, function(a) log2(a+1))
  
  dim(datSubsetTransformed)
  
  ## then the median scale
  datSubsetTransformed <- t(apply(datSubsetTransformed, 1, function(a) a - median(a, na.rm=T)))
  
  dim(datSubsetTransformed)
  
  load(paste0(sourceCodePath, 'Immune-Subtype-Clustering-master/shiny-app/Immune-Subtype-Clustering/data/comparative_immuneSigs_geneLists4.rda'))
  
  sigs1_2_geneIDs2<-as.character(na.omit(sigs1_2_eg2[[1]]$probe))
  for(i in 2:length(sigs1_2_eg2)){
    sigs1_2_geneIDs2<-c(sigs1_2_geneIDs2,as.character(na.omit(sigs1_2_eg2[[i]]$probe)))
  }
  sigs1_2_geneIDs2<-unique(sigs1_2_geneIDs2)   ## 2652 unique
  length(sigs1_2_geneIDs2)
  
  source(paste0(sourceCodePath, 'Immune-Subtype-Clustering-master/Scripts/ImmuneSigs68_function.R'))
  scores <- ImmuneSigs_function(datSubsetTransformed, 
                                sigs1_2_eg2,
                                sigs12_weighted_means,
                                sigs12_module_weights,
                                sigs1_2_names2,
                                sigs1_2_type2)
  
  write.csv(scores, paste0(savePath, "/sig160_result/raw_output/", cancer, "_", sample_type, ".csv"),
            row.names = T, quote = F)
  scores_z <- zscore.rows2(scores)
  write.csv(scores_z, paste0(savePath, "/sig160_result/zscore_output/", cancer, "_", sample_type, ".csv"),
            row.names = T, quote = F)
  
}





immuneSig_sig160_othersigs_consistency <- function(cancer,
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
  
  load(paste0(sig_dir, "originalSignaturesGenes/immuneSigGeneSetList_sig160_othersigs.Rdata"))
  
  
  # cancer = 'LIHC'
  cancer_type = paste0('TCGA-', cancer)
  
  processedDataPath <- paste0(processedDataPath, cancer,'/')
  dir.create(processedDataPath)
  
  dir.create(paste0(savePath, "/sig160ssgsea_result/"))
  dir.create(paste0(savePath, "/sig160ssgsea_result/raw_output/"))
  dir.create(paste0(savePath, "/sig160ssgsea_result/scaled_output/"))
  dir.create(paste0(savePath, "/sig160ssgsea_result/consistency_plot/"))
  dir.create(paste0(savePath, "/sig160ssgsea_result/consistency_plot_scaled/"))
  
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
                  gs = gs2, 
                  method="gsva",
                  mx.diff=FALSE, kcdf=kcdf,
                  verbose=FALSE, parallel.sz=1,
                  min.sz=min.sz, max.sz=max.sz,
                  ssgsea.norm=FALSE)
    
  }else if (method == 'ssgsea'){
    
    esOut <- gsva(data,
                  gs=gs2, 
                  method="ssgsea",
                  mx.diff=FALSE, kcdf=kcdf,
                  verbose=FALSE, parallel.sz=1,
                  min.sz=min.sz, max.sz=max.sz,
                  ssgsea.norm=TRUE)
    
  }
  save(esOut, 
       file = paste0(savePath,"/sig160ssgsea_result/raw_output/", method, "_", cancer,"_",sample_type,"_immueSig.Rdata"))
  
  
  #   median_scale <- function(x){
  #     return(t(apply(x, 1, function(x) (x - median(na.omit(x)))/mad(na.omit(x)))))
  #   }
  
  
  #   esOut_scaled <- median_scale(esOut)
  esOut_scaled <- zscore.rows2(esOut)
  
  save(esOut_scaled, 
       file = paste0(savePath,"/sig160ssgsea_result/scaled_output/", method, "_", cancer,"_",sample_type,"_immueSig.Rdata"))
  
  
  
  
  #   originalSigPath <- '/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/02_tumorGene_immuneSig_correlation/data/immune_signatures/'
  precompute_sig <- read.table(paste0(originalSigPath, 'Scores_160_Signatures.tsv'), header=T, sep='\t', stringsAsFactors = F)
  
  x <- str_replace_all(colnames(precompute_sig), pattern='\\.', replacement='-')
  colnames(precompute_sig) <- x
  
  
  sig = esOut[grep("_sig160",rownames( esOut)),]
  #   colnames(sig) <- substr(colnames(sig),1,12)
  rownames(sig) <- gsub("_sig160", "", rownames(sig))
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
    
    # b <- median_scale(sig2)
    b <- zscore.rows2(sig2)
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
        ggtitle(paste0(rownames(sig2)[i], 
                       "\n spearman's r = ", round(cor1_2s$estimate,3), 
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
    
    
    pdf(paste0(savePath, "/sig160ssgsea_result/consistency_plot/", 
               cancer, "_", sample_type,".pdf"),
        width = 5,
        height = 8, 
        onefile = TRUE)
    for(i in seq(length(p1))){
      print(p1[[i]])
    }
    dev.off()
    
    
    pdf(paste0(savePath, "/sig160ssgsea_result/consistency_plot_scaled/", 
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






immuneSig_sig160_othersigs <- function(data,
                                       cancer,
                                       sample_type, 
                                       method = c("gsva", "ssgsea"),
                                       kcdf = c("Poisson", "Gaussian"), 
                                       min.sz=1, 
                                       max.sz=1000,
                                       sig_dir,
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
  load(paste0(sig_dir, "originalSignaturesGenes/immuneSigGeneSetList_sig160_othersigs.Rdata"))
  
  
  dir.create(paste0(savePath, "/sig160ssgsea_result/"))
  dir.create(paste0(savePath, "/sig160ssgsea_result/raw_output/"))
  dir.create(paste0(savePath, "/sig160ssgsea_result/scaled_output/"))
  
  # Prepare TCGA expression matrix --------------------------------------------
  data = apply(data, 2, function(a) log2(a+1))
  
  ## gsva
  if(method == "gsva"){
    
    esOut <- gsva(data,
                  gs = gs2, 
                  method="gsva",
                  mx.diff=FALSE, kcdf=kcdf,
                  verbose=FALSE, parallel.sz=1,
                  min.sz=min.sz, max.sz=max.sz,
                  ssgsea.norm=FALSE)
    
  }else if (method == 'ssgsea'){
    
    esOut <- gsva(data,
                  gs=gs2, 
                  method="ssgsea",
                  mx.diff=FALSE, kcdf=kcdf,
                  verbose=FALSE, parallel.sz=1,
                  min.sz=min.sz, max.sz=max.sz,
                  ssgsea.norm=TRUE)
    
  }
  save(esOut, 
       file = paste0(savePath,"/sig160ssgsea_result/raw_output/", method, "_", cancer,"_",sample_type,"_immueSig.Rdata"))
  
  #   esOut_scaled <- median_scale(esOut)
  esOut_scaled <- zscore.rows2(esOut)
  
  save(esOut_scaled, 
       file = paste0(savePath,"/sig160ssgsea_result/scaled_output/", method, "_", cancer,"_",sample_type,"_immueSig.Rdata"))
}



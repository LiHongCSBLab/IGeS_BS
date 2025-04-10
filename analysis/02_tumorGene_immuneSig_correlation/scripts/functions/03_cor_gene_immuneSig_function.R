
immunecell_TCGAgene_cor_parallel_IHC <- function(cancer_index){
  
  setwd(workpath)
  library(psych)
  library(ppcor)
  library(reshape2)
  options(stringsAsFactors = F)
  # exprpath <- "/picb/bigdata/project/FengFYM/UCEC_Xena/EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena"
  # tcgaRNAmatrix <- as.data.frame(data.table::fread(exprpath))
  # colnames(tcgaRNAmatrix) <- substr(colnames(tcgaRNAmatrix),1 , 15)
  # limit to protein_coding_genes ------------------------------------------------
  protein_coding_genes <- read.table("03_drug_immuneSig_enrichment/data/protein_coding_genes_hg38.txt", sep = '\t', header = T)
  
  # sample_types <- read.table( paste0(workpath, "data/TCGA_cancertypes_subtypes.txt"), header=T, sep="\t")
  sample_types <- read.csv("02_tumorGene_immuneSig_correlation/data/TCGA_patient_ID_tumor_type_selected.csv")
  cancertypes <- unique(sample_types$TCGA_cancer_type)
  cancer = project_ids[cancer_index]
  
  load("01_tumor_purity/data/processedData/TCGA_tumor_purity_merged.Rdata")
  tumor_purity <- purity_merged
  
  resultpath_IHC1 <- "02_tumorGene_immuneSig_correlation/results/immunecell_TCGAgenes_pcor_pearson_IHC/"
  dir.create(resultpath_IHC1)
  outlog_IHC1 <- immunecell_TCGAgene_pcor(cancer = cancer,
                                          sample_types = sample_types,
                                          protein_coding_genes = protein_coding_genes,
                                          tumor_purity = tumor_purity,
                                          tumor_purity_method = "IHC",
                                          immune_sig = immune_sigatures,
                                          method = 'pearson',
                                          resultpath = resultpath_IHC1)
  print(paste0("finished calculating for ", outlog_IHC1))
  
  
  return(project_ids[cancer_index])
} 



immunecell_TCGAgene_cor_parallel_TUMERIC <- function(cancer_index){
  
  setwd(workpath)
  library(psych)
  library(ppcor)
  library(reshape2)
  options(stringsAsFactors = F)
  # exprpath <- "/picb/bigdata/project/FengFYM/UCEC_Xena/EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena"
  # tcgaRNAmatrix <- as.data.frame(data.table::fread(exprpath))
  # colnames(tcgaRNAmatrix) <- substr(colnames(tcgaRNAmatrix),1 , 15)
  # limit to protein_coding_genes ------------------------------------------------
  protein_coding_genes <- read.table("03_drug_immuneSig_enrichment/data/protein_coding_genes_hg38.txt", sep = '\t', header = T)
  
  # sample_types <- read.table( paste0(workpath, "data/TCGA_cancertypes_subtypes.txt"), header=T, sep="\t")
  sample_types <- read.csv("02_tumorGene_immuneSig_correlation/data/TCGA_patient_ID_tumor_type_selected.csv")
  cancertypes <- unique(sample_types$TCGA_cancer_type)
  load("01_tumor_purity/data/processedData/TCGA_tumor_purity_merged.Rdata")
  tumor_purity <-purity_merged
  cancer = project_ids[cancer_index]
  
  
  resultpath_TUMERIC1 <- "02_tumorGene_immuneSig_correlation/results/immunecell_TCGAgenes_result_pearson_TUMERIC/"
  dir.create(resultpath_TUMERIC1)
  outlog_TUMERIC1 <- immunecell_TCGAgene_pcor(cancer = cancer,
                                              sample_types = sample_types,
                                              protein_coding_genes = protein_coding_genes,
                                              tumor_purity = tumor_purity,
                                              tumor_purity_method = "TUMERIC",
                                              immune_sig = immune_sigatures,
                                              method = 'pearson',
                                              resultpath = resultpath_TUMERIC1)
  print(paste0("finished calculating for ", outlog_TUMERIC1))
  
  return(project_ids[cancer_index])
} 


immunecell_TCGAgene_cor_parallel_CPE <- function(cancer_index){
  
  setwd(workpath)
  library(psych)
  library(ppcor)
  library(reshape2)
  options(stringsAsFactors = F)
  # exprpath <- "/picb/bigdata/project/FengFYM/UCEC_Xena/EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena"
  # tcgaRNAmatrix <- as.data.frame(data.table::fread(exprpath))
  # colnames(tcgaRNAmatrix) <- substr(colnames(tcgaRNAmatrix),1 , 15)
  # limit to protein_coding_genes ------------------------------------------------
  protein_coding_genes <- read.table("03_drug_immuneSig_enrichment/data/protein_coding_genes_hg38.txt", sep = '\t', header = T)
  
  # sample_types <- read.table( paste0(workpath, "data/TCGA_cancertypes_subtypes.txt"), header=T, sep="\t")
  sample_types <- read.csv("02_tumorGene_immuneSig_correlation/data/TCGA_patient_ID_tumor_type_selected.csv")
  cancertypes <- unique(sample_types$TCGA_cancer_type)
  load("01_tumor_purity/data/processedData/TCGA_tumor_purity_merged.Rdata")
  tumor_purity <-purity_merged
  cancer = project_ids[cancer_index]
  
  resultpath_CPE1 <- "02_tumorGene_immuneSig_correlation/results/immunecell_TCGAgenes_result_pearson_CPE/"
  dir.create(resultpath_CPE1)
  outlog_CPE <- immunecell_TCGAgene_pcor(cancer = cancer,
                                         sample_types = sample_types,
                                              protein_coding_genes = protein_coding_genes,
                                         tumor_purity = tumor_purity,
                                         tumor_purity_method = "CPE",
                                         immune_sig = immune_sigatures,
                                         method = 'pearson',
                                         resultpath = resultpath_CPE1)
  print(paste0("finished calculating for ", outlog_CPE))
  
  
  return(project_ids[cancer_index])
} 



immunecell_TCGAgene_cor_parallel_spearman_TUMERIC <- function(cancer_index){
  
  setwd(workpath)
  library(psych)
  library(ppcor)
  library(reshape2)
  options(stringsAsFactors = F)
  # exprpath <- "/picb/bigdata/project/FengFYM/UCEC_Xena/EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena"
  # tcgaRNAmatrix <- as.data.frame(data.table::fread(exprpath))
  # colnames(tcgaRNAmatrix) <- substr(colnames(tcgaRNAmatrix),1 , 15)
  # limit to protein_coding_genes ------------------------------------------------
  protein_coding_genes <- read.table("03_drug_immuneSig_enrichment/data/protein_coding_genes_hg38.txt", sep = '\t', header = T)
  
  # sample_types <- read.table( paste0(workpath, "data/TCGA_cancertypes_subtypes.txt"), header=T, sep="\t")
  sample_types <- read.csv("02_tumorGene_immuneSig_correlation/data/TCGA_patient_ID_tumor_type_selected.csv")
  cancertypes <- unique(sample_types$TCGA_cancer_type)
  load("01_tumor_purity/data/processedData/TCGA_tumor_purity_merged.Rdata")
  tumor_purity <-purity_merged
  cancer = project_ids[cancer_index]
  
  resultpath_TUMERIC2 <- "02_tumorGene_immuneSig_correlation/results/immunecell_TCGAgenes_result_spearman_TUMERIC/"
  dir.create(resultpath_TUMERIC2)
  
  outlog_TUMERIC2 <- immunecell_TCGAgene_pcor(cancer = cancer,
                                              sample_types = sample_types,
                                              protein_coding_genes = protein_coding_genes,
                                              tumor_purity = tumor_purity,
                                              tumor_purity_method = "TUMERIC",
                                              immune_sig = immune_sigatures,
                                              method = 'spearman',
                                              resultpath = resultpath_TUMERIC2)
  print(paste0("finished calculating for ", outlog_TUMERIC2))
  
  # resultpath_ABSOLUTE1 <- "02_tumorGene_immuneSig_correlation/results/immunecell_TCGAgenes_result_pearson_ABSOLUTE/"
  # dir.create(resultpath_ABSOLUTE1)
  # outlog_ABSOLUTE1 <- immunecell_TCGAgene_pcor(cancer = cancer,
  #                                              sample_types = sample_types,
  #                                              tcgaRNAmatrix = tcgaRNAmatrix,
  #                                              tumor_purity = tumor_purity,
  #                                              tumor_purity_method = "ABSOLUTE",
  #                                              immune_sig = immune_sigatures,
  #                                              method = 'pearson',
  #                                              resultpath = resultpath_ABSOLUTE1)
  # print(paste0("finished calculating for ", outlog_ABSOLUTE1))
  # 
  
  # resultpath_IHC2 <- "02_tumorGene_immuneSig_correlation/results/immunecell_TCGAgenes_pcor_spearman_IHC/"
  # dir.create(resultpath_IHC2)
  
  # outlog_IHC2 <- immunecell_TCGAgene_pcor(cancer = cancer,
  #                                         sample_types = sample_types,
  #                                         # tcgaRNAmatrix = tcgaRNAmatrix,
  #                                         tumor_purity = tumor_purity,
  #                                         tumor_purity_method = "IHC",
  #                                         immune_sig = immune_sigatures,
  #                                         method = 'spearman',
  #                                         resultpath = resultpath_IHC2)
  # print(paste0("finished calculating for ", outlog_IHC2))
  
  
  return(project_ids[cancer_index])
} 


immunecell_TCGAgene_cor_parallel_spearman_IHC <- function(cancer_index){
  
  setwd(workpath)
  library(psych)
  library(ppcor)
  library(reshape2)
  options(stringsAsFactors = F)
  # exprpath <- "/picb/bigdata/project/FengFYM/UCEC_Xena/EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena"
  # tcgaRNAmatrix <- as.data.frame(data.table::fread(exprpath))
  # colnames(tcgaRNAmatrix) <- substr(colnames(tcgaRNAmatrix),1 , 15)
  
  # limit to protein_coding_genes ------------------------------------------------
  protein_coding_genes <- read.table("03_drug_immuneSig_enrichment/data/protein_coding_genes_hg38.txt", sep = '\t', header = T)
  # sample_types <- read.table( paste0(workpath, "data/TCGA_cancertypes_subtypes.txt"), header=T, sep="\t")
  sample_types <- read.csv("02_tumorGene_immuneSig_correlation/data/TCGA_patient_ID_tumor_type_selected.csv")
  cancertypes <- unique(sample_types$TCGA_cancer_type)
  load("01_tumor_purity/data/processedData/TCGA_tumor_purity_merged.Rdata")
  tumor_purity <-purity_merged
  cancer = project_ids[cancer_index]
  
  resultpath_IHC2 <- "02_tumorGene_immuneSig_correlation/results/immunecell_TCGAgenes_pcor_spearman_IHC/"
  dir.create(resultpath_IHC2)
  
  outlog_IHC2 <- immunecell_TCGAgene_pcor(cancer = cancer,
                                          sample_types = sample_types,
                                              protein_coding_genes = protein_coding_genes,
                                          tumor_purity = tumor_purity,
                                          tumor_purity_method = "IHC",
                                          immune_sig = immune_sigatures,
                                          method = 'spearman',
                                          resultpath = resultpath_IHC2)
  print(paste0("finished calculating for ", outlog_IHC2))
  
  return(project_ids[cancer_index])
} 

immunecell_TCGAgene_cor_parallel_spearman_CPE <- function(cancer_index){
  
  setwd(workpath)
  library(psych)
  library(ppcor)
  library(reshape2)
  options(stringsAsFactors = F)
  # exprpath <- "/picb/bigdata/project/FengFYM/UCEC_Xena/EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena"
  # tcgaRNAmatrix <- as.data.frame(data.table::fread(exprpath))
  # colnames(tcgaRNAmatrix) <- substr(colnames(tcgaRNAmatrix),1 , 15)
  # limit to protein_coding_genes ------------------------------------------------
  protein_coding_genes <- read.table("03_drug_immuneSig_enrichment/data/protein_coding_genes_hg38.txt", sep = '\t', header = T)
  
  # sample_types <- read.table( paste0(workpath, "data/TCGA_cancertypes_subtypes.txt"), header=T, sep="\t")
  sample_types <- read.csv("02_tumorGene_immuneSig_correlation/data/TCGA_patient_ID_tumor_type_selected.csv")
  cancertypes <- unique(sample_types$TCGA_cancer_type)
  load("01_tumor_purity/data/processedData/TCGA_tumor_purity_merged.Rdata")
  tumor_purity <-purity_merged
  cancer = project_ids[cancer_index]
  
  resultpath_CPE1 <- "02_tumorGene_immuneSig_correlation/results/immunecell_TCGAgenes_result_pearson_CPE/"
  dir.create(resultpath_CPE1)
  outlog_CPE1 <- immunecell_TCGAgene_pcor(cancer = cancer,
                                          sample_types = sample_types,
                                              protein_coding_genes = protein_coding_genes,
                                          tumor_purity = tumor_purity,
                                          tumor_purity_method = "CPE",
                                          immune_sig = immune_sigatures,
                                          method = 'spearman',
                                          resultpath = resultpath_CPE1)
  print(paste0("finished calculating for ", outlog_CPE1))
  
  return(project_ids[cancer_index])
} 




immunecell_TCGAgene_pcor <- function(cancer,
                                     sample_types,
                                     protein_coding_genes,
                                     tumor_purity,
                                     tumor_purity_method = c("IHC", "CPE", "TUMERIC", "ABSOLUTE"),
                                     immune_sig,
                                     method = c('pearson','spearman'),
                                     resultpath){
  library(ppcor)
  library(reshape2)
  #for(cancer in project_ids){
  print(paste0("start calculation for ", cancer))
  immune_sig$SampleID2 <- substr(immune_sig$ID, 1, 15)
  print("test Step 1 ")
  
  tumor_purity = tumor_purity[-which(duplicated(tumor_purity$Sample_ID)),]
  #  cancer = 'LIHC'
  sample_type = 'Primary'
  # prepare expression matrix
  # mymatrix_filter <- expr_matrix_preparation(cancer,
  #                                            sample_types,
  #                                            tcgaRNAmatrix)
  # 
  # Prepare TCGA expression matrix --------------------------------------------
  
  processedDataPath <- '02_tumorGene_immuneSig_correlation/data/TCGA_processedData/'
  processedDataPath <- paste0(processedDataPath, cancer,'/')
  
  # when running BRCA, also computing for subtypes ----------------------------
  if(cancer == 'BRCA'){
    
    sample_type = 'Primary'
    
    TPMMatPath = paste0(processedDataPath, "mixture_file_TCGA_", cancer, "_", sample_type, "_TPM.txt")
    mymatrix_filter <- read.csv(TPMMatPath, sep = '\t', row.names = 1)
    colnames(mymatrix_filter) <- gsub('\\.', '-',   colnames(mymatrix_filter))
    colnames(mymatrix_filter) <- substr(colnames(mymatrix_filter),1,15)
    mymatrix_filter <- mymatrix_filter[is.element(row.names(mymatrix_filter), protein_coding_genes$gene_name),]
    
    # matching files for correlation
    immune_sig$SampleID2 <- substr(immune_sig$ID, 1, 15)
    c.mat <- immune_sig[immune_sig$SampleID2 %in% unique(colnames(mymatrix_filter)),]
    c.mat <- c.mat[order(c.mat$SampleID2),]
    row.names(c.mat) <- c.mat$SampleID2
    
    t.mat <- mymatrix_filter[colnames(mymatrix_filter) %in% c.mat$SampleID2,]
    t.mat <- mymatrix_filter[, c.mat$SampleID2]
    
    tumor_purity_filter <- tumor_purity[is.element(tumor_purity$Sample_ID,  c.mat$SampleID2),]
    tumor_purity_filter2 <- unique(tumor_purity_filter[c('Sample_ID',tumor_purity_method)])
    names(tumor_purity_filter2) <- c('Sample_ID','purity')
    
    c.mat <- c.mat[which(is.element(unique(row.names(c.mat)), tumor_purity_filter2$Sample_ID)), ]
    c.mat <- c.mat[order(row.names(c.mat)),]
    t.mat <- t.mat[, which(is.element(unique(colnames(t.mat)), tumor_purity_filter2$Sample_ID)) ]
    t.mat <- t.mat[, order(colnames(t.mat))]
    t.mat <- log2(t.mat + 1)
    
    c.mat <- c.mat[- which(names( c.mat) %in% c("SampleID","SampleID2", 
                                                "PatientID",
                                                "TCGA.Study_immunelandscape", 
                                                "immunelandscape_TCGA.Study",
                                                "ID"))]
    print("test Step 2 ")
    
    if(nrow(tumor_purity_filter2[!is.na(tumor_purity_filter2$purity),]) > 10 ){
      row.names(tumor_purity_filter2) <- tumor_purity_filter2$Sample_ID
      tpf = tumor_purity_filter2[order(tumor_purity_filter2$Sample_ID),]
      
      # partial correlation
      print("calculating correlation after filtering low expressed genes")  
      # timestart <- Sys.time()
      
      print("test Step 3 ")
      
      
      cor_result <- lapply(as.list(c.mat), function(x){
        # [grep('CIBERSORT.ABS', colnames(c.mat))]
        if(length(which(is.na(x)))<0.5*length(x)){
          
          immuneSig = data.frame(immuneSig=x);
          # purity = data.frame(purity = tpf[,2]);
          tmp = cbind(immuneSig, t(as.matrix(t.mat)), tpf[2]);
          tmp = tmp[!is.na(tmp[,1]),];
          tmp = tmp[!is.na(tmp$purity),];
          
          cor_result <- sapply(tmp[-c(grep('immuneSig',colnames(tmp)),
                                      grep('purity',colnames(tmp)))],
                               function(y){
                                 tryCatch({pcor.test(x = tmp$immuneSig, 
                                                     y = y,
                                                     z = tmp$purity, 
                                                     method = method)},
                                          error = function(e){
                                            return(data.frame(estimate = NA,
                                                              p.value = NA,
                                                              statistic = NA,
                                                              n = NA, gp = NA, Method = NA))})
                               });
          cor_result <- as.data.frame(t(cor_result));
          cor_result <- cor_result[c('estimate','p.value')];
          names(cor_result) <- c("cor","p.value");
          cor_result$p.adj <- p.adjust(cor_result$p.value, method="fdr")
          return(cor_result)
        }
      })  
      
      save(cor_result, file = paste0(resultpath, "pcor_", method,"_", cancer,
                                     "_",sample_type, "_", tumor_purity_method,".Rdata"))
      
    }else{
      
      # timestart <- Sys.time()
      print("test Step 3 ")
      
      cor_result = lapply(c.mat, function(x){  corr.test(x, t(as.matrix(t.mat)), 
                                                         use = "pairwise", 
                                                         method = method, 
                                                         adjust = "none", 
                                                         ci = TRUE)})
      
      cor_result2 <- list()
      for(i in 1:length(cor_result)){
        cor_result[[i]]$ci$adj.p <- p.adjust(cor_result[[i]]$ci$p,method="fdr")
        cor_result2[[i]] <- cor_result[[i]]$ci
        cor_result2[[i]]$immune_sig <- names(cor_result)[i]
        cor_result2[[i]]$gene_name <- rownames(t.mat)
      }
      names(cor_result2) <- names(cor_result)
      
      save(cor_result2, file = paste0(resultpath, "cor_", method,"_", cancer,"_",sample_type, ".Rdata"))
      
    }
    
    
    sample_types <- sample_types[sample_types$TCGA_cancer_type == cancer, ]
    
    subtypes <- na.omit(unique(sample_types$TCGA.Subtype))
    
    for(type in subtypes){
      
      #type = subtypes[1]
      # prepare expression matrix
      
      cancer_samples <- sample_types[sample_types$TCGA.Subtype == type,]
      mymatrix <- mymatrix_filter[colnames(mymatrix_filter) %in% cancer_samples$TGCA_barcode]
      
      # matching files for correlation
      c.mat <- immune_sig[immune_sig$SampleID2 %in% unique(colnames(mymatrix)),]
      c.mat <- c.mat[order(c.mat$SampleID2),]
      row.names(c.mat) <- c.mat$SampleID2
      
      t.mat <- mymatrix[colnames(mymatrix) %in% c.mat$SampleID2,]
      t.mat <- mymatrix[, c.mat$SampleID2]
      
      tumor_purity_filter <- tumor_purity[is.element(tumor_purity$Sample_ID,  c.mat$SampleID2),]
      tumor_purity_filter2 <- unique(tumor_purity_filter[c('Sample_ID',tumor_purity_method)])
      names(tumor_purity_filter2) <- c('Sample_ID','purity')
      
      c.mat <- c.mat[which(is.element(unique(row.names(c.mat)), tumor_purity_filter2$Sample_ID)), ]
      c.mat <- c.mat[order(row.names(c.mat)),]
      t.mat <- t.mat[, which(is.element(unique(colnames(t.mat)), tumor_purity_filter2$Sample_ID)) ]
      t.mat <- t.mat[, order(colnames(t.mat))]
      t.mat <- log2(t.mat + 1)
      
      c.mat <- c.mat[- which(names( c.mat) %in% c("SampleID","SampleID2", "PatientID",
                                                  "TCGA.Study_immunelandscape", "immunelandscape_TCGA.Study",
                                                  "ID"))]
      print("test Step 2 ")
      
      if(nrow(tumor_purity_filter2[!is.na(tumor_purity_filter2$purity),]) > 10 ){
        row.names(tumor_purity_filter2) <- tumor_purity_filter2$Sample_ID
        tpf = tumor_purity_filter2[order(tumor_purity_filter2$Sample_ID),]
        # partial correlation
        print("test Step 3 ")
        
        
        cor_result <- lapply(as.list(c.mat), function(x){
          # [grep('CIBERSORT.ABS', colnames(c.mat))]
          if(length(which(is.na(x)))<0.5*length(x)){
            
            immuneSig = data.frame(immuneSig=x);
            # purity = data.frame(purity = tpf[,2]);
            tmp = cbind(immuneSig, t(as.matrix(t.mat)), tpf[2]);
            tmp = tmp[!is.na(tmp[,1]),];
            tmp = tmp[!is.na(tmp$purity),];
            
            cor_result <- sapply(tmp[-c(grep('immuneSig',colnames(tmp)),
                                        grep('purity',colnames(tmp)))],
                                 function(y){
                                   tryCatch({pcor.test(x = tmp$immuneSig, 
                                                       y = y,
                                                       z = tmp$purity, 
                                                       method = method)},
                                            error = function(e){
                                              return(data.frame(estimate = NA,
                                                                p.value = NA,
                                                                statistic = NA,
                                                                n = NA, gp = NA, Method = NA))})
                                 });
            cor_result <- as.data.frame(t(cor_result));
            cor_result <- cor_result[c('estimate','p.value')];
            names(cor_result) <- c("cor","p.value");
            cor_result$p.adj <- p.adjust(cor_result$p.value, method="fdr")
            return(cor_result)
          }
        })   
        
        save(cor_result, file = paste0(resultpath, "pcor_", method,"_", cancer,
                                       "_", type, "_", tumor_purity_method,".Rdata"))
      }else{
        
        cor_result = lapply(c.mat, function(x){corr.test(x, t(as.matrix(t.mat)), 
                                                         use = "pairwise", 
                                                         method = method, 
                                                         adjust = "none", 
                                                         ci = TRUE)})
        
        cor_result2 <- list()
        print("test Step 3 ")
        
        for(i in 1:length(cor_result)){
          cor_result[[i]]$ci$adj.p <- p.adjust(cor_result[[i]]$ci$p,method="fdr")
          cor_result2[[i]] <- cor_result[[i]]$ci
          cor_result2[[i]]$immune_sig <- names(cor_result)[i]
          cor_result2[[i]]$gene_name <- rownames(t.mat)
        }
        names(cor_result2) <- names(cor_result)
        
        save(cor_result2, file = paste0(resultpath, "cor_", method,"_", cancer,
                                        "_",type, ".Rdata"))
        
      }
    }
    
  }else if(is.element(cancer, c('COAD','READ'))){
    
    # when running COAD and READ, also computing for different MSI status ----------------------------
    
    sample_type = 'Primary'  
    TPMMatPath = paste0(processedDataPath, "mixture_file_TCGA_", cancer,"_",
                        sample_type, "_TPM.txt")
    mymatrix_filter <- read.csv(TPMMatPath, sep = '\t', row.names = 1)
    colnames(mymatrix_filter) <- gsub('\\.', '-',   colnames(mymatrix_filter))
    colnames(mymatrix_filter) <- substr(colnames(mymatrix_filter),1,15)
    mymatrix_filter <- mymatrix_filter[is.element(row.names(mymatrix_filter), protein_coding_genes$gene_name),]

    # matching files for correlation
    immune_sig$SampleID2 <- substr(immune_sig$ID, 1, 15)
    c.mat <- immune_sig[immune_sig$SampleID2 %in% unique(colnames(mymatrix_filter)),]
    c.mat <- c.mat[order(c.mat$SampleID2),]
    row.names(c.mat) <- c.mat$SampleID2
    
    t.mat <- mymatrix_filter[colnames(mymatrix_filter) %in% c.mat$SampleID2,]
    t.mat <- mymatrix_filter[, c.mat$SampleID2]
    
    tumor_purity_filter <- tumor_purity[is.element(tumor_purity$Sample_ID,  c.mat$SampleID2),]
    tumor_purity_filter2 <- unique(tumor_purity_filter[c('Sample_ID',tumor_purity_method)])
    names(tumor_purity_filter2) <- c('Sample_ID','purity')
    
    c.mat <- c.mat[which(is.element(unique(row.names(c.mat)), tumor_purity_filter2$Sample_ID)), ]
    c.mat <- c.mat[order(row.names(c.mat)),]
    t.mat <- t.mat[, which(is.element(unique(colnames(t.mat)), tumor_purity_filter2$Sample_ID)) ]
    t.mat <- t.mat[, order(colnames(t.mat))]
    t.mat <- log2(t.mat + 1)
    
    c.mat <- c.mat[- which(names( c.mat) %in% c("SampleID","SampleID2", 
                                                "PatientID",
                                                "TCGA.Study_immunelandscape", 
                                                "immunelandscape_TCGA.Study",
                                                "ID"))]
    print("test Step 2 ")
    
    if(nrow(tumor_purity_filter2[!is.na(tumor_purity_filter2$purity),]) > 10 ){
      row.names(tumor_purity_filter2) <- tumor_purity_filter2$Sample_ID
      tpf = tumor_purity_filter2[order(tumor_purity_filter2$Sample_ID),]
      
      # partial correlation
      print("test Step 3 ")
      
      
      cor_result <- lapply(as.list(c.mat), function(x){
        # [grep('CIBERSORT.ABS', colnames(c.mat))]
        if(length(which(is.na(x)))<0.5*length(x)){
          
          immuneSig = data.frame(immuneSig=x);
          # purity = data.frame(purity = tpf[,2]);
          tmp = cbind(immuneSig, t(as.matrix(t.mat)), tpf[2]);
          tmp = tmp[!is.na(tmp[,1]),];
          tmp = tmp[!is.na(tmp$purity),];
          
          cor_result <- sapply(tmp[-c(grep('immuneSig',colnames(tmp)),
                                      grep('purity',colnames(tmp)))],
                               function(y){
                                 tryCatch({pcor.test(x = tmp$immuneSig, 
                                                     y = y,
                                                     z = tmp$purity, 
                                                     method = method)},
                                          error = function(e){
                                            return(data.frame(estimate = NA,
                                                              p.value = NA,
                                                              statistic = NA,
                                                              n = NA, gp = NA, 
                                                              Method = NA))})
                               });
          cor_result <- as.data.frame(t(cor_result));
          cor_result <- cor_result[c('estimate','p.value')];
          names(cor_result) <- c("cor","p.value");
          cor_result$p.adj <- p.adjust(cor_result$p.value, method="fdr")
          return(cor_result)
        }
      })        
      save(cor_result, file = paste0(resultpath, "pcor_", method,"_", cancer,
                                     "_",sample_type, "_", tumor_purity_method,".Rdata"))
      
    }else{
      
      # timestart <- Sys.time()
      print("test Step 3 ")
      
      cor_result = lapply(c.mat, function(x){  corr.test(x, t(as.matrix(t.mat)), 
                                                         use = "pairwise", 
                                                         method = method, 
                                                         adjust = "none", 
                                                         ci = TRUE)})
      
      cor_result2 <- list()
      for(i in 1:length(cor_result)){
        cor_result[[i]]$ci$adj.p <- p.adjust(cor_result[[i]]$ci$p,method="fdr")
        cor_result2[[i]] <- cor_result[[i]]$ci
        cor_result2[[i]]$immune_sig <- names(cor_result)[i]
        cor_result2[[i]]$gene_name <- rownames(t.mat)
      }
      names(cor_result2) <- names(cor_result)
      
      save(cor_result2, file = paste0(resultpath, "cor_", method,"_", cancer,"_",sample_type, ".Rdata"))
      
    }
    
    
    sample_types <- sample_types[sample_types$TCGA_cancer_type == cancer, ]
    subtype <- na.omit(unique(sample_types$MSI_status))
    subtype <- subtype[-which(subtype == "Indeterminate")]
    
    
    for(type in subtype){
      # type = subtype[1]
      cancer_samples <- sample_types[which(sample_types$MSI_status == type),]
      mymatrix <- mymatrix_filter[colnames(mymatrix_filter) %in% cancer_samples$TGCA_barcode]
      
      # matching files for correlation
      c.mat <- immune_sig[immune_sig$SampleID2 %in% unique(colnames(mymatrix)),]
      c.mat <- c.mat[order(c.mat$SampleID2),]
      row.names(c.mat) <- c.mat$SampleID2
      
      t.mat <- mymatrix[colnames(mymatrix) %in% c.mat$SampleID2,]
      t.mat <- mymatrix[, c.mat$SampleID2]
      
      tumor_purity_filter <- tumor_purity[is.element(tumor_purity$Sample_ID,  c.mat$SampleID2),]
      tumor_purity_filter2 <- unique(tumor_purity_filter[c('Sample_ID',tumor_purity_method)])
      names(tumor_purity_filter2) <- c('Sample_ID','purity')
      
      
      c.mat <- c.mat[which(is.element(unique(row.names(c.mat)), tumor_purity_filter2$Sample_ID)), ]
      c.mat <- c.mat[order(row.names(c.mat)),]
      t.mat <- t.mat[, which(is.element(unique(colnames(t.mat)), tumor_purity_filter2$Sample_ID)) ]
      t.mat <- t.mat[, order(colnames(t.mat))]
      t.mat <- log2(t.mat + 1)
      
      c.mat <- c.mat[- which(names( c.mat) %in% c("SampleID","SampleID2", 
                                                  "PatientID",
                                                  "TCGA.Study_immunelandscape", 
                                                  "immunelandscape_TCGA.Study",
                                                  "ID"))]
      print("test Step 2 ")
      
      if(nrow(tumor_purity_filter2[!is.na(tumor_purity_filter2$purity),]) > 10 ){
        
        row.names(tumor_purity_filter2) <- tumor_purity_filter2$Sample_ID
        tpf = tumor_purity_filter2[order(tumor_purity_filter2$Sample_ID),]
        
        # partial correlation
        print("calculating correlation after filtering low expressed genes")  
        # timestart <- Sys.time()
        
        print("test Step 3 ")
        
        cor_result <- lapply(as.list(c.mat), function(x){
          # [grep('CIBERSORT.ABS', colnames(c.mat))]
          if(length(which(is.na(x)))<0.5*length(x)){
            
            immuneSig = data.frame(immuneSig=x);
            # purity = data.frame(purity = tpf[,2]);
            tmp = cbind(immuneSig, t(as.matrix(t.mat)), tpf[2]);
            tmp = tmp[!is.na(tmp[,1]),];
            tmp = tmp[!is.na(tmp$purity),];
            
            cor_result <- sapply(tmp[-c(grep('immuneSig',colnames(tmp)),
                                        grep('purity',colnames(tmp)))],
                                 function(y){
                                   tryCatch({pcor.test(x = tmp$immuneSig, 
                                                       y = y,
                                                       z = tmp$purity, 
                                                       method = method)},
                                            error = function(e){
                                              return(data.frame(estimate = NA,
                                                                p.value = NA,
                                                                statistic = NA,
                                                                n = NA, gp = NA, 
                                                                Method = NA))})
                                 });
            cor_result <- as.data.frame(t(cor_result));
            cor_result <- cor_result[c('estimate','p.value')];
            names(cor_result) <- c("cor","p.value");
            cor_result$p.adj <- p.adjust(cor_result$p.value, method="fdr")
            return(cor_result)
          }
        })  
        
        save(cor_result, file = paste0(resultpath, "pcor_", method,"_", cancer,
                                       "_",type, "_", tumor_purity_method,".Rdata"))
      }else{
        
        cor_result = lapply(c.mat, function(x){  corr.test(x, t(as.matrix(t.mat)), 
                                                           use = "pairwise", 
                                                           method = method, 
                                                           adjust = "none", 
                                                           ci = TRUE)})
        
        print("test Step 3 ")
        cor_result2 <- list()
        
        for(i in 1:length(cor_result)){
          cor_result[[i]]$ci$adj.p <- p.adjust(cor_result[[i]]$ci$p,method="fdr")
          cor_result2[[i]] <- cor_result[[i]]$ci
          cor_result2[[i]]$immune_sig <- names(cor_result)[i]
          cor_result2[[i]]$gene_name <- rownames(t.mat)
        }
        names(cor_result2) <- names(cor_result)
        
        save(cor_result2, file = paste0(resultpath, "cor_", method,"_", cancer,
                                        "_",type, ".Rdata"))
        
      }
    }
    
  }else if(cancer == 'SKCM'){
    # SKCM has both metastatic and primary, compute for both seperately --------------------
    for(sample_type in c('Primary','Metastatic')){
      
      TPMMatPath = paste0(processedDataPath, "mixture_file_TCGA_", cancer,"_",
                          sample_type, "_TPM.txt")
      
      mymatrix_filter <- read.csv(TPMMatPath, sep = '\t', row.names = 1)
      colnames(mymatrix_filter) <- gsub('\\.', '-',   colnames(mymatrix_filter))
      colnames(mymatrix_filter) <- substr(colnames(mymatrix_filter),1,15)
      mymatrix_filter <- mymatrix_filter[is.element(row.names(mymatrix_filter), protein_coding_genes$gene_name),]

      # matching files for correlation
      immune_sig$SampleID2 <- substr(immune_sig$ID, 1, 15)
      c.mat <- immune_sig[immune_sig$SampleID2 %in% unique(colnames(mymatrix_filter)),]
      c.mat <- c.mat[order(c.mat$SampleID2),]
      row.names(c.mat) <- c.mat$SampleID2
      
      t.mat <- mymatrix_filter[colnames(mymatrix_filter) %in% c.mat$SampleID2,]
      t.mat <- mymatrix_filter[, c.mat$SampleID2]
      
      tumor_purity_filter <- tumor_purity[is.element(tumor_purity$Sample_ID,  c.mat$SampleID2),]
      tumor_purity_filter2 <- unique(tumor_purity_filter[c('Sample_ID',tumor_purity_method)])
      names(tumor_purity_filter2) <- c('Sample_ID','purity')
      
      c.mat <- c.mat[which(is.element(unique(row.names(c.mat)), tumor_purity_filter2$Sample_ID)), ]
      c.mat <- c.mat[order(row.names(c.mat)),]
      t.mat <- t.mat[, which(is.element(unique(colnames(t.mat)), tumor_purity_filter2$Sample_ID)) ]
      t.mat <- t.mat[, order(colnames(t.mat))]
      t.mat <- log2(t.mat + 1)
      
      c.mat <- c.mat[- which(names( c.mat) %in% c("SampleID","SampleID2", 
                                                  "PatientID",
                                                  "TCGA.Study_immunelandscape", 
                                                  "immunelandscape_TCGA.Study",
                                                  "ID"))]
      print("test Step 2 ")
      
      if(nrow(tumor_purity_filter2[!is.na(tumor_purity_filter2$purity),]) > 10 ){
        row.names(tumor_purity_filter2) <- tumor_purity_filter2$Sample_ID
        tpf = tumor_purity_filter2[order(tumor_purity_filter2$Sample_ID),]
        
        # partial correlation
        print("calculating correlation after filtering low expressed genes")  
        # timestart <- Sys.time()
        
        print("test Step 3 ")
        
        cor_result <- lapply(as.list(c.mat), function(x){
          # [grep('timer', colnames(c.mat))]
          if(length(which(is.na(x)))<0.5*length(x)){
            
            immuneSig = data.frame(immuneSig=x);
            # purity = data.frame(purity = tpf[,2]);
            tmp = cbind(immuneSig, t(as.matrix(t.mat)), tpf[2]);
            tmp = tmp[!is.na(tmp[,1]),];
            tmp = tmp[!is.na(tmp$purity),];
            
            cor_result <- sapply(tmp[-c(grep('immuneSig',colnames(tmp)),
                                        grep('purity',colnames(tmp)))],
                                 function(y){
                                   tryCatch({pcor.test(x = tmp$immuneSig, 
                                                       y = y,
                                                       z = tmp$purity, 
                                                       method = method)},
                                            error = function(e){
                                              return(data.frame(estimate = NA,
                                                                p.value = NA,
                                                                statistic = NA,
                                                                n = NA, gp = NA, 
                                                                Method = NA))})
                                 });
            cor_result <- as.data.frame(t(cor_result));
            cor_result <- cor_result[c('estimate','p.value')];
            names(cor_result) <- c("cor","p.value");
            cor_result$p.adj <- p.adjust(cor_result$p.value, method="fdr")
            return(cor_result)
          }
        })  
        
        
        save(cor_result, file = paste0(resultpath, "pcor_", method,"_", cancer,
                                       "_",sample_type, "_", tumor_purity_method,".Rdata"))
        
      }else{
        
        cor_result = lapply(c.mat, function(x){  corr.test(x, t(as.matrix(t.mat)),
                                                           use = "pairwise", 
                                                           method = method, 
                                                           adjust = "none", 
                                                           ci = TRUE)})
        print("test Step 3 ")
        
        cor_result2 <- list()
        for(i in 1:length(cor_result)){
          cor_result[[i]]$ci$adj.p <- p.adjust(cor_result[[i]]$ci$p,method="fdr")
          cor_result2[[i]] <- cor_result[[i]]$ci
          cor_result2[[i]]$immune_sig <- names(cor_result)[i]
          cor_result2[[i]]$gene_name <- rownames(t.mat)
        }
        names(cor_result2) <- names(cor_result)
        
        save(cor_result2, file = paste0(resultpath, "cor_", method,"_", cancer,
                                        "_",sample_type, ".Rdata"))
        
      }
    }
    
  }else{
    
    sample_type = 'Primary'
    TPMMatPath = paste0(processedDataPath, "mixture_file_TCGA_", cancer,"_",
                        sample_type, "_TPM.txt")
    mymatrix_filter <- read.csv(TPMMatPath, sep = '\t', row.names = 1)
    colnames(mymatrix_filter) <- gsub('\\.', '-',   colnames(mymatrix_filter))
    colnames(mymatrix_filter) <- substr(colnames(mymatrix_filter),1,15)
    mymatrix_filter <- mymatrix_filter[is.element(row.names(mymatrix_filter), protein_coding_genes$gene_name),]

    # matching files for correlation
    immune_sig$SampleID2 <- substr(immune_sig$ID, 1, 15)
    c.mat <- immune_sig[immune_sig$SampleID2 %in% unique(colnames(mymatrix_filter)),]
    c.mat <- c.mat[order(c.mat$SampleID2),]
    row.names(c.mat) <- c.mat$SampleID2
    
    t.mat <- mymatrix_filter[colnames(mymatrix_filter) %in% c.mat$SampleID2,]
    t.mat <- mymatrix_filter[, c.mat$SampleID2]
    
    tumor_purity_filter <- tumor_purity[is.element(tumor_purity$Sample_ID,  c.mat$SampleID2),]
    tumor_purity_filter2 <- unique(tumor_purity_filter[c('Sample_ID',tumor_purity_method)])
    names(tumor_purity_filter2) <- c('Sample_ID','purity')
    
    c.mat <- c.mat[which(is.element(unique(row.names(c.mat)), tumor_purity_filter2$Sample_ID)), ]
    c.mat <- c.mat[order(row.names(c.mat)),]
    t.mat <- t.mat[, which(is.element(unique(colnames(t.mat)), tumor_purity_filter2$Sample_ID)) ]
    t.mat <- t.mat[, order(colnames(t.mat))]
    t.mat <- log2(t.mat + 1)
    
    c.mat <- c.mat[- which(names( c.mat) %in% c("SampleID","SampleID2", 
                                                "PatientID","type",
                                                "TCGA.Study_immunelandscape", 
                                                "immunelandscape_TCGA.Study",
                                                "ID"))]
    print("test Step 2 ")
    
    if(nrow(tumor_purity_filter2[!is.na(tumor_purity_filter2$purity),]) > 10 ){
      row.names(tumor_purity_filter2) <- tumor_purity_filter2$Sample_ID
      tpf = tumor_purity_filter2[order(tumor_purity_filter2$Sample_ID),]
      
      # partial correlation
      print("calculating correlation after filtering low expressed genes")  
      # timestart <- Sys.time()
      
      print("test Step 3 ")
      
      cor_result <- lapply(as.list(c.mat), function(x){
        # [grep('CIBERSORT.ABS', colnames(c.mat))]
        # x = as.list(c.mat)[[103]]
        if(length(which(is.na(x)))<0.5*length(x)){
          
          immuneSig = data.frame(immuneSig=x);
          # purity = data.frame(purity = tpf[,2]);
          tmp = cbind(immuneSig, t(as.matrix(t.mat)), tpf[2]);
          tmp = tmp[!is.na(tmp[,1]),];
          tmp = tmp[!is.na(tmp$purity),];
          # dim(tmp)
          cor_result <- sapply(tmp[-c(grep('immuneSig',colnames(tmp)),
                                      grep('purity',colnames(tmp)))],
                               function(y){
                                 tryCatch({pcor.test(x = tmp$immuneSig, 
                                                     y = y,
                                                     z = tmp$purity, 
                                                     method = method)},
                                          error = function(e){
                                            return(data.frame(estimate = NA,
                                                              p.value = NA,
                                                              statistic = NA,
                                                              n = NA, gp = NA, 
                                                              Method = NA))})
                               });
          cor_result <- as.data.frame(t(cor_result));
          cor_result <- cor_result[c('estimate','p.value')];
          names(cor_result) <- c("cor","p.value");
          cor_result$p.adj <- p.adjust(cor_result$p.value, method="fdr")
          return(cor_result)
        }
      })  
      
      # timeend <- Sys.time()
      # runningtime <- timeend-timestart
      # print(runningtime)
      # total runningtime estimation for one cancer type: 
      # runningtime * ncol(c.mat) / 2 / 60
      # in total about 5.3 hour for one cancer
      
      save(cor_result, file = paste0(resultpath, "pcor_", method,"_", cancer,
                                     "_",sample_type, "_", tumor_purity_method,".Rdata"))
      
    }else{
      
      # timestart <- Sys.time()
      print("test Step 3 ")
      
      cor_result = lapply(c.mat, function(x){  corr.test(x, t(as.matrix(t.mat)), 
                                                         use = "pairwise", 
                                                         method = method, 
                                                         adjust = "none", 
                                                         ci = TRUE)})
      # timeend <- Sys.time()
      # runningtime <- timeend-timestart
      # print(runningtime)
      # total runningtime estimation for one cancer type: 
      # runningtime * ncol(c.mat) / 2 / 60
      # in total about 5.3 hour for one cancer
      
      cor_result2 <- list()
      for(i in 1:length(cor_result)){
        cor_result[[i]]$ci$adj.p <- p.adjust(cor_result[[i]]$ci$p,method="fdr")
        cor_result2[[i]] <- cor_result[[i]]$ci
        cor_result2[[i]]$immune_sig <- names(cor_result)[i]
        cor_result2[[i]]$gene_name <- rownames(t.mat)
      }
      names(cor_result2) <- names(cor_result)
      
      save(cor_result2, file = paste0(resultpath, "cor_", method,"_", cancer,"_",sample_type, ".Rdata"))
      
    }
  }
  
  return(cancer)
}
 

immunecell_TCGAgene_cor <- function(cancer,
                                    sample_types,
                                    protein_coding_genes,
                                    sample_type = 'Primary',
                                    method = c('pearson','spearman'),
                                    immune_sig,
                                    resultpath){
  library(psych)
  library(reshape2)
  #for(cancer in project_ids){
  print(paste0("start calculation for ", cancer))
  
  
  
  
  TPMMatPath = paste0(processedDataPath, "mixture_file_TCGA_", cancer,"_",sample_type, "_TPM.txt")
  mymatrix_filter <- read.csv(TPMMatPath, sep = '\t', row.names = 1)
  colnames(mymatrix_filter) <- gsub('\\.', '-',   colnames(mymatrix_filter))
  colnames(mymatrix_filter) <- substr(colnames(mymatrix_filter),1,15)
  mymatrix_filter <- mymatrix_filter[is.element(row.names(mymatrix_filter), protein_coding_genes$gene_name),]

  # matching files for correlation
  immune_sig$SampleID2 <- substr(immune_sig$ID, 1, 15)
  c.mat <- immune_sig[immune_sig$SampleID2 %in% unique(colnames(mymatrix_filter)),]
  c.mat <- c.mat[order(c.mat$SampleID2),]
  row.names(c.mat) <- c.mat$SampleID2
  
  t.mat <- mymatrix_filter[colnames(mymatrix_filter) %in% c.mat$SampleID2,]
  t.mat <- mymatrix_filter[, c.mat$SampleID2]
  t.mat <- log2(t.mat + 1)
  tumor_purity_filter <- tumor_purity[is.element(tumor_purity$Sample_ID,  c.mat$SampleID2),]
  tumor_purity_filter2 <- unique(tumor_purity_filter[c('Sample_ID',tumor_purity_method)])
  names(tumor_purity_filter2) <- c('Sample_ID','purity')
  
  
  c.mat <- c.mat[- which(names( c.mat) %in% c("SampleID","SampleID2", "PatientID",
                                              "TCGA.Study_immunelandscape", 
                                              "immunelandscape_TCGA.Study",
                                              "ID"))]
  
  
  print("calculating correlation after filtering low expressed genes")  
  
  # tmp=cor(as.matrix(c.mat[-c(1,2)]), t(as.matrix(t.mat)),use="pairwise.complete.obs")
  
  cor_result = lapply(c.mat, function(x){  corr.test(x, t(as.matrix(t.mat)), use = "pairwise", 
                                                     method = method, adjust = "none", ci = TRUE)})
  
  cor_result2 <- list()
  for(i in 1:length(cor_result)){
    cor_result[[i]]$ci$adj.p <- p.adjust(cor_result[[i]]$ci$p,method="fdr")
    cor_result2[[i]] <- cor_result[[i]]$ci
    cor_result2[[i]]$immune_sig <- names(cor_result)[i]
    cor_result2[[i]]$gene_name <- rownames(t.mat)
  }
  names(cor_result2) <- names(cor_result)
  
  save(cor_result2, file = paste0(resultpath, "cor_", method,"_", cancer, ".Rdata"))
  
  return(cancer)
}


expr_matrix_preparation <- function(cancer,
                                    sample_types,
                                    tcgaRNAmatrix ){
  
  cancer_samples <- sample_types[sample_types$TCGA_cancer_type == cancer,]
  mymatrix <- tcgaRNAmatrix[colnames(tcgaRNAmatrix) %in% c("sample",cancer_samples$TGCA_barcode)]
  #  mymatrix <-  aggregate(mymatrix[-which(colnames(mymatrix)=="sample")],by=list(mymatrix$sample),median)
  #  rownames(mymatrix) <- mymatrix$Group.1
  mymatrix <- mymatrix[!duplicated(mymatrix$sample),]
  rownames(mymatrix) <- mymatrix$sample
  mymatrix <- mymatrix[-1]
  if(length(which(is.na(rowMeans(mymatrix))))>0){
    mymatrix <- mymatrix[-which(is.na(rowMeans(mymatrix))),]
  }
  
  return(mymatrix)
  
}


# Annotate gene names

gene_ensembl_mapping <- function(genelist, resultpath){  
  library(biomaRt)
  
  ensembl <- useMart(biomart="ensembl",dataset="hsapiens_gene_ensembl")
  result_gene_ensembl <- getBM(values = genelist,
                               filters = "ensembl_gene_id", 
                               mart = ensembl,
                               attributes = c("ensembl_gene_id","hgnc_symbol"))
  write.table(result_gene_ensembl, paste0(resultpath, "data/result_gene_ensembl.txt"), sep="\t",row.names=F, quote=F)
  
  return(result_gene_ensembl)
}







geneExprInCCLE <- function(x){
  # ccle_expr <- data.table::fread("03_drug_immuneSig_enrichment/data/CCLE_geneexpr/CCLE_RNAseq_reads.csv",sep=',')
  # rownames(ccle_expr) <- ccle_expr$V1
  # ccle_expr <- as.data.frame(t(ccle_expr))
  # colnames(ccle_expr) <- ccle_expr[1,]
  # ccle_expr <- ccle_expr[-1,]
  # write.csv(ccle_expr, "03_drug_immuneSig_enrichment/data/CCLE_geneexpr/CCLE_RNAseq_reads_processed.csv", quote = F)
  ccle_expr <- read.csv("03_drug_immuneSig_enrichment/data/CCLE_geneexpr/CCLE_RNAseq_reads_processed.csv", row.names = 1)
  ccle_expr_rowname <- row.names(ccle_expr)
  ccle_expr_rowname_gene <- unlist(lapply(as.list(ccle_expr_rowname), function(x)unlist(strsplit(x, fixed = T, split = " ("))[1]))
  ccle_expr_rowname_geneID <- unlist(lapply(as.list(ccle_expr_rowname), function(x)gsub(")","",unlist(strsplit(x, fixed = T, split = " ("))[2])))
  rownames(ccle_expr) <- ccle_expr_rowname
  ccle_expr$X = ccle_expr_rowname 
  ccle_expr$gene_name = ccle_expr_rowname_gene
  ccle_expr$gene_id = ccle_expr_rowname_geneID
  geneExprInCCLE <- ccle_expr[c("X","gene_name","gene_id")]
  geneExprInCCLEstat <- apply(ccle_expr[-which(is.element(colnames(ccle_expr), c("X","gene_name","gene_id")))],1,
                              function(x)length(which(x == 0)))
  # countsMat <- ccle_expr[-which(is.element(colnames(ccle_expr), c("X","gene_name","gene_id")))]
  # countsMatkeep <- rowSums(countsMat>0) >= floor(0.8*ncol(countsMat))# 
  # countsMatkeep <- countsMat[countsMatkeep,] 
  
  geneExprInCCLEstat <- data.frame(geneExprInCCLEstat)
  geneExprInCCLEstat$X <- rownames(geneExprInCCLEstat) 
  geneExprInCCLE <- inner_join(geneExprInCCLE, geneExprInCCLEstat)
  geneExprInCCLE$rate <- geneExprInCCLE$geneExprInCCLEstat / ncol(ccle_expr[-which(is.element(colnames(ccle_expr), c("X","gene_name","gene_id")))])              
  write.csv(geneExprInCCLE, "03_drug_immuneSig_enrichment/data/CCLE_geneexpr/geneExprInCCLE.csv", row.names = F, quote = F)
  
}




gene_immune_pcor_sig <- function(file_name,
                                 immunesig_path,
                                 resultpath,
                                 savepath_allresult,
                                 savepath, 
                                 cancer,
                                 sample_type,
                                 tumor_purity_method,
                                 method = "spearman",
                                 s1 = 0.4,
                                 s2 = -0.4,
                                 padj = 0.05, # 10e-10,
                                 auc_threshold = 0.6,
                                 immuSigproof = 2, # how many datasets support the immune sigture to have prediction power
                                 ITSproof = 2,
                                #  immunesig_confident_filter = 0.5,
                                 num_gene_p = 200,
                                 num_gene_n = 200) {
  
  library(reshape2)
  library(stringr)
  library(dplyr)

  savepath_allresult = paste0(savepath_allresult, "/immunecell_TCGAgenes_result_", cor_method,"_", tumor_purity_method, "_r",s1,"/") # ,"_ITSproof_",ITSproof
 # savepath_allresult = paste0(savepath_allresult, "/immunecell_TCGAgenes_result_", cor_method,"_", tumor_purity_method, "_r",s1,"_ITSproof_",ITSproof,"/")
  dir.create(savepath_allresult)
  savepath = paste0(savepath, "gene_immune_sig_", tumor_purity_method, "_",cor_method,"_r",s1,"_ITSproof_",ITSproof,"/")
  dir.create(savepath)

  dir.create(paste0(savepath, "for_GSVA/"))
  dir.create(paste0(savepath, "for_GSVA/cor"))
  dir.create(paste0(savepath, "for_GSVA/pcor"))
  dir.create(paste0(savepath, "for_GSEA/"))
  dir.create(paste0(savepath, "for_GSEA/cor"))
  dir.create(paste0(savepath, "for_GSEA/pcor"))

  # method = "spearman"
  # cancer = 'LIHC'
  # sample_type = 'Primary'
  # tumor_purity_method = 'TUMERIC'
  
  # Filter criteria: -----------------------------------------------------------
  # 1. significance: p.adj < 0.05
  # 2. protein-coding genes
  # 3. express in above 80% samples in original TCGA dataset
  # 3. expressed in 80% of cell lines in CCLE
  # 4. high correlation (positive or negative): 
  #     positive correlated: cor > *
  #     negative correlated: cor < -*
  # ----------------------------------------------------------------------------
  print("# --------------------------------------------------------------------")
  print("1. load in filter files")
  # protein_coding_genes 
  # hg38 <- readGFF("/picb/bigdata/project/FengFYM/s7_annotation_files/gencode.v22.annotation.gtf")
  # anno <- setDT(hg38)
  # anno <- as.data.frame(anno[gene_type == "protein_coding", ])
  # protein_coding_genes <- unique(anno[c("gene_name","gene_id","gene_type")])
  # protein_coding_genes$gene_id <- substr(protein_coding_genes$gene_id, 1, 15)
  # write.table(protein_coding_genes, "03_drug_immuneSig_enrichment/data/protein_coding_genes_hg38.txt", sep = '\t', row.names = F, quote = F)
  protein_coding_genes <- read.table("03_drug_immuneSig_enrichment/data/protein_coding_genes_hg38.txt", sep = '\t', header = T)
  
  # express in above 80% samples in original TCGA dataset
  processedDataPath <- '02_tumorGene_immuneSig_correlation/data/TCGA_processedData/'
  processedDataPath <- paste0(processedDataPath, cancer,'/')
  TPMMatPath = paste0(processedDataPath, "mixture_file_TCGA_", cancer, "_", sample_type, "_TPM.txt")
  TPMMat <- read.csv(TPMMatPath, sep = '\t', row.names = 1)
  TPMMatkeep <- rowSums(TPMMat>0) >= floor(0.8*ncol(TPMMat))
  TPMMatkeep <- TPMMat[TPMMatkeep,] 
  genetokeep <- data.frame(gene_name = rownames(TPMMatkeep))
  protein_coding_genes <- inner_join(genetokeep, protein_coding_genes)
  
  
  # gene expressed in > 80% of CCLE cell lines ---------------------------------
  # file generated with geneExprInCCLE()
  geneExprInCCLE <- read.csv("03_drug_immuneSig_enrichment/data/CCLE_geneexpr/geneExprInCCLE.csv")
  geneExprInCCLE_filted <- geneExprInCCLE[geneExprInCCLE$rate < 0.2,]
  
  # file = paste0(resultpath, "pcor_", method,"_", cancer,
  #               "_",sample_type, "_", tumor_purity_method,".Rdata")
  
  file = paste0(resultpath, file_name)
  cancer = unlist(strsplit(file_name, split = "_"))[3]
  sample_type = unlist(strsplit(file_name, split = "_"))[4]
  tumor_purity_method = gsub('.Rdata', '', unlist(strsplit(file_name, split = "_"))[5])
  
  print("2. load correlation files")
  
  load(file)
  
  genelist <- lapply(cor_result, function(x){
    
    #summary(x$r)
    # x = cor_result[[325]]
    x = data.frame(r = unlist(x$cor),
                   p.value = unlist(x$p.value),
                   adj.p = unlist(x$p.adj))
    # keep proteincoding genes
    x1 <- x[which(is.element(rownames(x), protein_coding_genes$gene_name)),]
    # keep genes expressed in >80% cell lines
    x1 <- x1[which(is.element(rownames(x1), geneExprInCCLE_filted$gene_name)),]
    
    # x1 <- x1[x1$adj.p < padj,]
     
    # s2 <- summary(x[x$r < 0,]$r)[2]
    y <- unique(x1[which(x1$adj.p < padj),])
    y <- y[!is.na(y$r),]
    y$gene_name = rownames(y)
    
    if(nrow(y) == 0){
      g_p = NULL
      g_n = NULL
      return(list(g_p, g_n))
    }else{
      y_p <- y[which(y$r > s1),]
      y_n <- y[which(y$r < s2),]
      
      y_p <- y_p[order(y_p$r, decreasing = T), ]
      y_n <- y_n[order(y_n$r, decreasing = F), ]
      
      #length(y)
      if(nrow(y_p) ==0){
        g_p = NULL
      }else{
        if(nrow(y_p) <= num_gene_p){
          
          g_p <- na.omit(unique(y_p$gene_name))
          
        }else{
          g_p <- na.omit(unique(y_p$gene_name[1:num_gene_p]))
        }
        
      }
      
      if(nrow(y_n) ==0){
        g_n = NULL
      }else{
        
        if(nrow(y) <= num_gene_n) {
          g_n <- y_n$gene_name
        }else{
          g_n <- na.omit(unique(y_n$gene_name[1:num_gene_n]))
        }
      }
      
      # g_p <- paste(g_p, collapse = "\t")
      return(list(p=g_p, n=g_n))
    }
  })
  
  names(genelist) = names(cor_result)
  
  genelist_p <- lapply(genelist, function(x){ x[[1]]  })
  genelist_n <- lapply(genelist, function(x){ x[[2]]  })
  #   plot(density(unlist(lapply(genelist_p,length))))
  #   lines(density(unlist(lapply(genelist_n,length))), col='red')
  lapply(genelist_p, length)
  lapply(genelist_n, length)
  

  # drop out xcell non-TME related cells
    print("4.drop out xcell non-TME related cells.")
  # prefilter out xcell unrelated cell type
  xcell_cellfamily <- read.csv("02_tumorGene_immuneSig_correlation/data/immune_signatures/originalSignaturesGenes/immune_sig_original_files/XCELL/cells_families.txt", sep = '\t')
  xcell_type1 <- xcell_cellfamily[is.element(xcell_cellfamily$Group, 
                                           c("B-cells",
                                             "CD4+ memory T-cells",
                                             "CD4+ T-cells",
                                             "CD8+ T-cells",
                                             "DC",
                                             "Endothelial cells",
                                             "Epithelial cells",
                                             "Fibroblasts",
                                             "Macrophages")), ]

  xcell_type2 <- xcell_cellfamily[is.element(xcell_cellfamily$Group, 
                                           "Parent"), ]
  xcell_type3 <- xcell_type2[is.element(xcell_type2$Type, 
                                      c("Lymphoid", 
                                        "Myeloid")), ]
  xcell_type4 <- xcell_type2[is.element(xcell_type2$Cells, 
                                      c("Endothelial cells", 
                                        "Epithelial cells",
                                        "Fibroblasts")), ]
  xcell_select <- rbind(xcell_type1, xcell_type3, xcell_type4)
  xcell_remove <- xcell_cellfamily[!is.element(xcell_cellfamily$Cells, xcell_select$Cells),]
  xcell_remove <- paste0(xcell_remove$Cells, "_xcell")
  xcell_remove <- gsub(" ", ".", xcell_remove)
  
  genelist_p <- genelist_p[-which(is.element(names(genelist_p), xcell_remove))]
  genelist_n <- genelist_n[-which(is.element(names(genelist_n), xcell_remove))]

  # function analysis
  print("4. save all gene sets before further filtering.")
  

  save(genelist, file = paste0(savepath_allresult, "/pcor_", 
                               method,"_", cancer, "_",sample_type, "_all_", num_gene_p, ".Rdata"))
  save(genelist_p, file = paste0(savepath_allresult, "/pcor_", 
                                 method,"_", cancer, "_",sample_type, "_positive_", num_gene_p, ".Rdata"))
  save(genelist_n, file = paste0(savepath_allresult, "/pcor_", 
                                 method,"_", cancer, "_",sample_type, "_negative_", num_gene_n, ".Rdata"))
}
                                 
                                   


ITS_fileparation <- function(savepath_allresult,
                                 cancer,
                                 sample_type,
                                 method = "spearman",
                                 tumor_purity_method,
                                 s1 = 0.4,
                                 num_gene_p = 200,
                                 num_gene_n = 200) {
  
  library(reshape2)
  library(stringr)
  library(dplyr)

  savepath_allresult = paste0(savepath_allresult, "/immunecell_TCGAgenes_result_", method,"_", tumor_purity_method, "_r",s1,"/") # ,"_ITSproof_",ITSproof
  dir.create(paste0(savepath_allresult, "for_GSEA/"))
  dir.create(paste0(savepath_allresult, "for_GSEA/pcor"))
  
  load(paste0(savepath_allresult, "/pcor_", 
                                 method,"_", cancer, "_",sample_type, "_positive_", num_gene_p, ".Rdata"))
  load(paste0(savepath_allresult, "/pcor_", 
                                 method,"_", cancer, "_",sample_type, "_negative_", num_gene_n, ".Rdata"))
  
  ITSname = data.frame(immunesig = names(genelist_p), immune_sig = tolower(names(genelist_p)))
  ITSname = immunesig_name_unify(ITSname)
  names(genelist_p) = ITSname$immune_sig

  ITSname = data.frame(immunesig = names(genelist_n), immune_sig = tolower(names(genelist_n)))
  ITSname = immunesig_name_unify(ITSname)
  names(genelist_n) = ITSname$immune_sig

  # generate gmt files for GSEA analysis ---------------------------------------- 
  print("9. save ITS for GSEA analysis.")

  genelist_p2 <- sapply(genelist_p, function(x)paste(x, collapse = "\t"))
  genelist_n2 <- sapply(genelist_n, function(x)paste(x, collapse = "\t"))
  
  genelist_p2 <- cbind(matrix(NA, length(genelist_p2)), genelist_p2)
  genelist_n2 <- cbind(matrix(NA, length(genelist_n2)), genelist_n2)
  
  if(nrow(genelist_p2) != 0){
    write.table(genelist_p2, paste0(savepath_allresult, "for_GSEA/pcor/", 
                                    method,"_", cancer, "_",sample_type, "_positive_", num_gene_p, ".gmt"), 
                sep='\t', row.names=T, col.names = F, quote=F)
  }
  if(nrow(genelist_n2) != 0){
    write.table(genelist_n2, paste0(savepath_allresult, "for_GSEA/pcor/", 
                                    method,"_", cancer, "_",sample_type, "_negative_", num_gene_n, ".gmt"), 
                sep='\t', row.names=T, col.names = F, quote=F)
  }
  
  print("End of analysis")
  print("# --------------------------------------------------------------------")
  
}


ITS_fileparation_oe_original <- function(savepath_allresult,
                                 cancer,
                                 sample_type,
                                 method = "spearman",
                                 tumor_purity_method,
                                 s1 = 0.4,
                                 num_gene_p = 200,
                                 num_gene_n = 200) {
  
  library(reshape2)
  library(stringr)
  library(dplyr)
  library(readxl)


  savepath_allresult2 = paste0( savepath_allresult, "/ITS_oe_original_", method,"_", tumor_purity_method, "_r",s1,"/") # ,"_ITSproof_",ITSproof
  dir.create(paste0(savepath_allresult2, "/"))
  dir.create(paste0(savepath_allresult2, "for_GSEA/"))
  dir.create(paste0(savepath_allresult2, "for_GSEA/pcor"))
  dir.create(paste0(savepath_allresult2, "for_GSVA/"))
  dir.create(paste0(savepath_allresult2, "for_GSVA/pcor"))

  savepath_allresult = paste0(savepath_allresult, "/immunecell_TCGAgenes_result_", method,"_", tumor_purity_method, "_r",s1,"/") # ,"_ITSproof_",ITSproof
  dir.create(paste0(savepath_allresult, "for_GSEA/"))
  dir.create(paste0(savepath_allresult, "for_GSEA/pcor"))
  
  load(paste0(savepath_allresult, "/pcor_", 
                                 method,"_", cancer, "_",sample_type, "_positive_", num_gene_p, ".Rdata"))
  load(paste0(savepath_allresult, "/pcor_", 
                                 method,"_", cancer, "_",sample_type, "_negative_", num_gene_n, ".Rdata"))
  

  ITSname = data.frame(immunesig = names(genelist_p), immune_sig = tolower(names(genelist_p)))
  ITSname = immunesig_name_unify(ITSname)
  names(genelist_p) = ITSname$immune_sig

  ITSname = data.frame(immunesig = names(genelist_n), immune_sig = tolower(names(genelist_n)))
  ITSname = immunesig_name_unify(ITSname)
  names(genelist_n) = ITSname$immune_sig

  

  # replace tumor sig with original genelists ---------------------------------------------------------------------------- 
  immuneSigInfoPath = "06_results_summaryandplotting/data/immene_cell_hieracy/"
  
  immunesigInfo <- as.data.frame(read_excel(paste0(immuneSigInfoPath, "immune_sig_hieracy_new.xlsx"), sheet = 1))
  immunesigInfo <- immunesigInfo[c(1,4,5:9)]
  immunesigInfo$immune_sig = tolower(immunesigInfo$immune_sig)
  immunesigInfo = immunesig_name_unify(immunesigInfo)
  
  immunesigInfoSelected = immunesigInfo # [is.element(immunesigInfo$immune_sig, ITSselected2$immune_sig), ]
  tumorsig = immunesigInfoSelected[is.element(immunesigInfoSelected$Type, c("Tumor cell signatures", "Melanoma cell state")), ]
  
  
  immunesigInfoSelected = immunesigInfoSelected[!is.element(immunesigInfoSelected$Type, c("Tumor cell signatures", "Melanoma cell state")), ]
  
  TMEcellSig = immunesigInfoSelected[immunesigInfoSelected$Classification == "TME",]
  immuneSig = immunesigInfoSelected[immunesigInfoSelected$Classification == "Immune signatures",]
  icb_predictor = immunesigInfoSelected[immunesigInfoSelected$Classification == "ICB predictor",]
  ICBresisSig = immunesigInfoSelected[immunesigInfoSelected$Classification == "ICB resistant signatures",]
  TMESig = immunesigInfoSelected[immunesigInfoSelected$Classification == "TME signatures",]
  
  immuneSig = rbind(immuneSig, rbind(icb_predictor, TMESig))
  
  # LOAD original immune signature gene sets ---------------------------------
  immuneSigGSpath = "02_tumorGene_immuneSig_correlation/data/immune_signatures/originalSignaturesGenes/"
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
  
  
  # ICB resistant signatures ---------------------------------------------------------------
  if(length(intersect(names(genelist_p), ICBresisSig$immune_sig)) > 1){
    
    icb_res_path = "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/05_immuneSig_ICBresponse/results/"
    # resall <- read.csv(paste0(icb_res_path, "immunesig_icbResponse_all_results.csv"))
    resall <- read.csv("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/05_immuneSig_ICBresponse/results/meta_results/res_metaRes.csv")
    names(resall) <- c("immune_sig", names(resall)[-1])
    resall_oe <- resall[grep("oe_", resall$immune_sig),]
    resall_oe <- resall_oe[-grep("oe_tme", resall_oe$immune_sig),]
    resall_oe <- resall_oe[-grep("ayers", resall_oe$immune_sig),]
    resall_oe <- resall_oe[-grep("oe_co.culture.screen.hits10", resall_oe$immune_sig),]
    
    oe_kept = union(resall_oe$immune_sig, union(tumorsig$immune_sig, TMESig$immune_sig))


    ITS_ICBresisSig_p_new = immunesigGS[intersect(names(immunesigGS), oe_kept)]
    ITS_ICBresisSig_n_new <- lapply(ITS_ICBresisSig_p_new, function(x){NULL })
    
  }
  
    genelist_p1 = genelist_p[!is.element(names(genelist_p), names(ITS_ICBresisSig_p_new))]
    genelist_p = c(genelist_p1, ITS_ICBresisSig_p_new)
 
    genelist_n1 = genelist_n[!is.element(names(genelist_n), names(ITS_ICBresisSig_n_new))]
    genelist_n = c(genelist_n1, ITS_ICBresisSig_n_new)
  
    save(genelist_p, file = paste0(savepath_allresult2, "/for_GSVA/pcor/pcor_", 
                                  method,"_", cancer, "_",sample_type, "_positive_", num_gene_p, ".Rdata"))
    save(genelist_n, file = paste0(savepath_allresult2, "/for_GSVA/pcor//pcor_", 
                                  method,"_", cancer, "_",sample_type, "_negative_", num_gene_n, ".Rdata"))

  # generate gmt files for GSEA analysis ---------------------------------------- 
  print("9. save ITS for GSEA analysis.")

  genelist_p2 <- sapply(genelist_p, function(x)paste(x, collapse = "\t"))
  genelist_n2 <- sapply(genelist_n, function(x)paste(x, collapse = "\t"))
  
  genelist_p2 <- cbind(matrix(NA, length(genelist_p2)), genelist_p2)
  genelist_n2 <- cbind(matrix(NA, length(genelist_n2)), genelist_n2)
  
  if(nrow(genelist_p2) != 0){
    write.table(genelist_p2, paste0(savepath_allresult2, "/for_GSEA/pcor/", 
                                    method,"_", cancer, "_",sample_type, "_positive_", num_gene_p, ".gmt"), 
                sep='\t', row.names=T, col.names = F, quote=F)
  }
  if(nrow(genelist_n2) != 0){
    write.table(genelist_n2, paste0(savepath_allresult2, "/for_GSEA/pcor/", 
                                    method,"_", cancer, "_",sample_type, "_negative_", num_gene_n, ".gmt"), 
                sep='\t', row.names=T, col.names = F, quote=F)
  }
  
  print("End of analysis")
  print("# --------------------------------------------------------------------")
  
}





gene_immune_pcor_sig_finefilter <- function(file_name,
                                            immunesig_path,
                                            resultpath,
                                            savepath_allresult,
                                            savepath, 
                                            cancer,
                                            sample_type,
                                            tumor_purity_method,
                                            method = "spearman",
                                            s1 = 0.4,
                                            s2 = -0.4,
                                            padj = 0.05, # 10e-10,
                                            auc_threshold = 0.6,
                                            immuSigproof = 2, # how many datasets support the immune sigture to have prediction power
                                            ITSproof = 2,
                                            #  immunesig_confident_filter = 0.5,
                                            num_gene_p = 200,
                                            num_gene_n = 200) {
  
  library(reshape2)
  library(stringr)
  library(dplyr)

  savepath_allresult = paste0(savepath_allresult, "/immunecell_TCGAgenes_result_", cor_method,"_", tumor_purity_method, "_r",s1,"_ITSproof_",ITSproof,"/")
  dir.create(savepath_allresult)
  savepath = paste0(savepath, "gene_immune_sig_", tumor_purity_method, "_",cor_method,"_r",s1,"_ITSproof_",ITSproof,"/")
  dir.create(savepath)

  dir.create(paste0(savepath, "for_GSVA/"))
  dir.create(paste0(savepath, "for_GSVA/cor"))
  dir.create(paste0(savepath, "for_GSVA/pcor"))
  dir.create(paste0(savepath, "for_GSEA/"))
  dir.create(paste0(savepath, "for_GSEA/cor"))
  dir.create(paste0(savepath, "for_GSEA/pcor"))

  file = paste0(resultpath, file_name)
  cancer = unlist(strsplit(file_name, split = "_"))[3]
  sample_type = unlist(strsplit(file_name, split = "_"))[4]
  tumor_purity_method = gsub('.Rdata', '', unlist(strsplit(file_name, split = "_"))[5])

  load(paste0(savepath_allresult, "/pcor_", method,"_", cancer, "_",sample_type, "_positive_", num_gene_p, ".Rdata"))
  load(paste0(savepath_allresult, "/pcor_", method,"_", cancer, "_",sample_type, "_negative_", num_gene_n, ".Rdata"))

  # drop ITS with less than 3 genes --------------------------------------------
  print("6. drop ITS with less than 3 genes.")
  drop_sig =  data.frame(table(c(rownames(data.frame(which(lapply(genelist_p,length) <= 3))),
                                 rownames(data.frame(which(lapply(genelist_n,length) <= 3))))))
  drop_sig = drop_sig[drop_sig$Freq == 2,]
  
  if(nrow(drop_sig)>0){
    
    # lapply(genelist_p[is.element(names(genelist_p), drop_sig$Var1)],length)
    # lapply(genelist_n[is.element(names(genelist_n), drop_sig$Var1)],length)
    
    genelist_p <- genelist_p[!is.element(names(genelist_p), drop_sig$Var1)]
    genelist_n <- genelist_n[!is.element(names(genelist_n), drop_sig$Var1)]
    
    #   plot(density(unlist(lapply(genelist_p,length))))
    #   lines(density(unlist(lapply(genelist_n,length))), col='red')
  }
  
  lapply(genelist_p,length)
  lapply(genelist_n,length)
  


  # function dection -----------------------------------------------------------
  print("7. conduct function analysis for each ITS.")
  print("GO-BP analysis ############################")
  
  

  # correlation between gene expression and purity------------------------------
  # genelist_filtered <- ITSexpr_purity_dection(genelist_p = genelist_p,
  #                                             genelist_n = genelist_n,
  #                                             r_threshold = -0.4,
  #                                             cancer = cancer,
  #                                             sample_type = sample_type,
  #                                             tumor_purity_method = tumor_purity_method)
  # genelist_p <- genelist_filtered[[1]]
  # genelist_n <- genelist_filtered[[2]]
  # function dection confirm ITS -----------------------------------------------
  # 02_tumorGene_immuneSig_correlation/scripts/correlation_dection_plot.R
  # ITSenrichment()
  #                                           
  #
  #   plot(density(unlist(lapply(genelist_p,length))))
  #   lines(density(unlist(lapply(genelist_n,length))), col='red')
  
  # merge ITSs: ---------------------------------------------------------------- 
  # 1. whose original immune signatures are describing same cells or biological
  # meaning.
  # 
  # 2. whose jaccard index are high (like > 0.7)
  # 
  # replace ITS: the original  signatures that is no need to regeneration ITS -----------------------
  # 
  # 

      # ----------------------------------------------------------------------------
  # select ITS whose original immunesig is either ICB predictive or strongly relate to survival  
  # ----------------------------------------------------------------------------
  print("5. filter out ITSs failed for either ICB predictive analysis or survival analysis.")
   
  # immunesig <- read.csv(paste0(immunesig_path,"/", cancer,"_",sample_type,".csv"))
  # immunesig <- immunesig[which(immunesig$confident >= immunesig_confident_filter ), ]

# ITSproof=4
  # immunesig <- immune_sig_filter_v2(auc_threshold = auc_threshold,
  #                                     immuSigproof = immuSigproof, # how many datasets support the immune sigture to have prediction power
  #                                     ITSproof = ITSproof, # how many datasets support the immune sigture to have prediction power
  #                                     workpath = "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/",
  #                                     savepath  = "06_results_summaryandplotting/data/immune_sig_selected_new/")

  immunesig <- read.table(paste0("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/06_results_summaryandplotting/data/immune_sig_selected_new/selected_ITS_",s1,".txt"),
                             sep = "\t", header = T)
  
  immunesig <- immunesig[immunesig$Freq >= ITSproof, ]

  names(genelist_p) <- tolower( names(genelist_p))
  names(genelist_n) <- tolower( names(genelist_n))
  # length(setdiff(immunesig$immune_sig, names(genelist_p)))
  signame = data.frame(x =  names(genelist_p), immune_sig = names(genelist_p))
  signame = immunesig_name_unify(signame)
  # length(setdiff(immunesig$immune_sig, signame$immune_sig))
  names(genelist_p) <- signame$immune_sig
  names(genelist_n) <- signame$immune_sig
  genelist_p <- genelist_p[is.element(names(genelist_p), immunesig$immune_sig)]
  genelist_n <- genelist_n[is.element(names(genelist_n), immunesig$immune_sig)]
  lapply(genelist_p,length)
  lapply(genelist_n,length)

  print("8. save ITS as Rdata.")
  
  save(genelist_p, file = paste0(savepath, "for_GSVA/pcor/", 
                                 method,"_", cancer, "_",sample_type, "_positive_", num_gene_p, ".Rdata"))
  save(genelist_n, file = paste0(savepath, "for_GSVA/pcor/", 
                                 method,"_", cancer, "_",sample_type, "_negative_", num_gene_n, ".Rdata"))
  
  
  # generate gmt files for GSEA analysis ---------------------------------------- 
  print("9. save ITS for GSEA analysis.")
  
  genelist_p2 <- sapply(genelist_p, function(x)paste(x, collapse = "\t"))
  genelist_n2 <- sapply(genelist_n, function(x)paste(x, collapse = "\t"))
  
  genelist_p2 <- cbind(matrix(NA, length(genelist_p2)), genelist_p2)
  genelist_n2 <- cbind(matrix(NA, length(genelist_n2)), genelist_n2)
  
  if(nrow(genelist_p2) != 0){
    write.table(genelist_p2, paste0(savepath, "for_GSEA/pcor/", 
                                    method,"_", cancer, "_",sample_type, "_positive_", num_gene_p, ".gmt"), 
                sep='\t', row.names=T, col.names = F, quote=F)
  }
  if(nrow(genelist_n2) != 0){
    write.table(genelist_n2, paste0(savepath, "for_GSEA/pcor/", 
                                    method,"_", cancer, "_",sample_type, "_negative_", num_gene_n, ".gmt"), 
                sep='\t', row.names=T, col.names = F, quote=F)
  }
  
  print("End of analysis")
  print("# --------------------------------------------------------------------")
  
}



gene_immune_cor_sig <- function(file_name,
                                immunesig_path,
                                resultpath,
                                savepath_allresult,
                                savepath, 
                                cancer,
                                sample_type,
                                tumor_purity_method,
                                method = "spearman",
                                s1 = 0.4,
                                s2 = -0.4, 
                                padj = 0.05,
                                auc_threshold = 0.6,
                                immuSigproof = 2, # how many datasets support the immune sigture to have prediction power
                                ITSproof = 2,
                                # immunesig_confident_filter = 0.5,
                                num_gene_p = 200,
                                num_gene_n = 200 ) {
  library(reshape2)
  library(stringr)
  library(dplyr)

  print("# --------------------------------------------------------------------")
  print("1. load in filter files")
  
  # 1. protein_coding_genes 
  protein_coding_genes <- read.table("03_drug_immuneSig_enrichment/data/protein_coding_genes_hg38.txt", sep = '\t', header = T)
  
  # express in above 80% samples in original TCGA dataset
  processedDataPath <- '02_tumorGene_immuneSig_correlation/data/TCGA_processedData/'
  processedDataPath <- paste0(processedDataPath, cancer,'/')
  TPMMatPath = paste0(processedDataPath, "mixture_file_TCGA_", cancer, "_", sample_type, "_TPM.txt")
  TPMMat <- read.csv(TPMMatPath, sep = '\t', row.names = 1)
  TPMMatkeep <- rowSums(TPMMat>0) >= floor(0.8*ncol(TPMMat))
  TPMMatkeep <- TPMMat[TPMMatkeep,] 
  genetokeep <- data.frame(gene_name = rownames(TPMMatkeep))
  protein_coding_genes <- inner_join(genetokeep, protein_coding_genes)
  
   
  # 2. gene expressed in CCLE-cell lines
  geneExprInCCLE <- read.csv("03_drug_immuneSig_enrichment/data/CCLE_geneexpr/geneExprInCCLE.csv")
  geneExprInCCLE_filted <- geneExprInCCLE[geneExprInCCLE$rate < 0.2,]
  
  
  file = paste0(resultpath, file_name)
  cancer = unlist(strsplit(file_name, split = "_"))[3]
  sample_type = gsub(".Rdata", "", unlist(strsplit(file_name, split = "_"))[4])
  
  print("2. load correlation files")
  load(file)
  
  genelist <- lapply(cor_result2, function(x){
    #summary(x$r)
    # x = cor_result2[[1]]
    x1 <- x[which(is.element(x$gene_name, protein_coding_genes$gene_name)),]
    # keep genes expressed in >80% cell lines
    x1 <- x1[which(is.element(x1$gene_name, geneExprInCCLE_filted$gene_name)),]
    
    # s2 <- summary(x[x$r < 0,]$r)[2]
    y <- unique(x1[which(x1$adj.p < padj),])
    y <- y[!is.na(y$r),]

    if(nrow(y) == 0){
      g_p = NULL
      g_n = NULL
      return(list(g_p, g_n))
    }else{
      
      y_p <- y[which(y$r > s1),]
      y_n <- y[which(y$r < s2),]
      
      y_p <- y_p[order(y_p$r, decreasing = T), ]
      y_n <- y_n[order(y_n$r, decreasing = F), ]
      
      #length(y)
      if(nrow(y_p) ==0){
        g_p = NULL
      }else{
        if(nrow(y_p) <= num_gene_p){
          g_p <- na.omit(unique(y_p$gene_name))
          
        }else{
          g_p <- na.omit(unique(y_p$gene_name[1:num_gene_p]))
        }
      }
      
      if(nrow(y_n) ==0){
        g_n = NULL
      }else{
        if(nrow(y) <= num_gene_n) {
          g_n <- y_n$gene_name
        }else{
          if(nrow(y_n) <= num_gene_n){
            
            g_n <- na.omit(unique(y_n$gene_name))
            
          }else{
            g_n <- na.omit(unique(y_n$gene_name[1:num_gene_n]))
          }
        }
      }
      
      # g_p <- paste(g_p, collapse = "\t")
      return(list(g_p, g_n))
    }
  })
  
  names(genelist) = names(cor_result2)
  
  genelist_p <- lapply(genelist, function(x){ x[[1]]  })
  genelist_n <- lapply(genelist, function(x){ x[[2]]  })
  
  names(genelist) = names(cor_result2)
  
  genelist_p <- lapply(genelist, function(x){ x[[1]]  })
  genelist_n <- lapply(genelist, function(x){ x[[2]]  })
  
  save(genelist, file = paste0(savepath_allresult, "/cor_", 
                               method,"_", cancer, "_",sample_type, "_all_", num_gene_p, ".Rdata"))
  save(genelist_p, file = paste0(savepath_allresult, "/cor_", 
                                 method,"_", cancer, "_",sample_type, "_positive_", num_gene_p, ".Rdata"))
  save(genelist_n, file = paste0(savepath_allresult, "/cor_", 
                                 method,"_", cancer, "_",sample_type, "_negative_", num_gene_n, ".Rdata"))
  
  # ----------------------------------------------------------------------------
  # select ITS whose original immunesig is either ICB predictive or strongly relate to survival  
  # ----------------------------------------------------------------------------
  print("4. filter out ITSs failed for either ICB predictive analysis or survival analysis.")
  
  # immunesig <- read.csv(paste0(immunesig_path,"/", cancer,"_",sample_type,".csv"))
  # immunesig <- immunesig[which(immunesig$confident >= immunesig_confident_filter ), ]

  immunesig <- immune_sig_filter_v2(auc_threshold = auc_threshold,
                                    immuSigproof = immuSigproof, # how many datasets support the immune sigture to have prediction power
                                    ITSproof = ITSproof, # how many datasets support the immune sigture to have prediction power
                                    workpath = "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/",
                                    savepath  = "06_results_summaryandplotting/data/immune_sig_selected_new/")

  names(genelist_p) <- tolower( names(genelist_p))
  names(genelist_n) <- tolower( names(genelist_n))
  # length(setdiff(immunesig$immune_sig, names(genelist_p)))
  signame = data.frame(x =  names(genelist_p), immune_sig = names(genelist_p))
  signame = immunesig_name_unify(signame)
  # length(setdiff(immunesig$immune_sig, signame$immune_sig))
  names(genelist_p) <- signame$immune_sig
  names(genelist_n) <- signame$immune_sig
  genelist_p <- genelist_p[is.element(names(genelist_p), immunesig$immune_sig)]
  genelist_n <- genelist_n[is.element(names(genelist_n), immunesig$immune_sig)]
  # lapply(genelist_p,length)
  # lapply(genelist_n,length)
  
  #   plot(density(unlist(lapply(genelist_p,length))))
  #   lines(density(unlist(lapply(genelist_n,length))), col='red')
  # drop ITS with less than 3 genes --------------------------------------------
  print("5. drop ITS with less than 3 genes.")
  drop_sig =  data.frame(table(c(rownames(data.frame(which(lapply(genelist_p,length) <= 3))),
                                 rownames(data.frame(which(lapply(genelist_n,length) <= 3))))))
  drop_sig = drop_sig[drop_sig$Freq == 2,]
  if(nrow(drop_sig)>0){
    
    # lapply(genelist_p[is.element(names(genelist_p), drop_sig$Var1)],length)
    # lapply(genelist_n[is.element(names(genelist_n), drop_sig$Var1)],length)
    
    genelist_p <- genelist_p[!is.element(names(genelist_p), drop_sig$Var1)]
    genelist_n <- genelist_n[!is.element(names(genelist_n), drop_sig$Var1)]
    # lapply(genelist_p,length)
    # lapply(genelist_n,length)
    
    # plot(density(unlist(lapply(genelist_p,length))))
    # lines(density(unlist(lapply(genelist_n,length))), col='red')
  }
  

  # correlation between gene expression and purity------------------------------
  # ITSexpr_purity_dection(genelist_p,
  #                        genelist_n,
  #                        r_threshold,
  #                        cancer = 'LIHC',
  #                        sample_type = 'Primary',
  #                        tumor_purity_method = "TUMERIC")

  save(genelist_p, file = paste0(savepath, "for_GSVA/cor/", method, "_", cancer, "_",sample_type, "_positive_", num_gene_p, ".Rdata"))
  save(genelist_n, file = paste0(savepath, "for_GSVA/cor/", method, "_", cancer, "_",sample_type, "_negative_", num_gene_n, ".Rdata"))

  # generate gmt files for GSEA analysis ---------------------------------------- 
  print("7. save ITS for GSEA analysis.")
  
  genelist2 <- sapply(genelist, function(x)paste(x, collapse = "\t"))
  genelist_p2 <- sapply(genelist_p, function(x)paste(x, collapse = "\t"))
  genelist_n2 <- sapply(genelist_n, function(x)paste(x, collapse = "\t"))
  
  genelist2 <- cbind(matrix(NA, length(genelist2)), genelist2)
  genelist_p2 <- cbind(matrix(NA, length(genelist_p2)), genelist_p2)
  genelist_n2 <- cbind(matrix(NA, length(genelist_n2)), genelist_n2)
  
  if(nrow(genelist_p2) != 0){
    write.table(genelist_p2, paste0(savepath, "for_GSEA/cor/", method, "_", cancer, "_",sample_type, "_positive_", num_gene_p, ".gmt"), sep='\t', row.names=T, col.names = F, quote=F)
  }
  
  if(nrow(genelist_n2) != 0){ 
    write.table(genelist_n2, paste0(savepath, "for_GSEA/cor/", method, "_", cancer, "_",sample_type, "_negative_", num_gene_n, ".gmt"), sep='\t', row.names=T, col.names = F, quote=F)
  }
  
  print("End of analysis")
  print("# --------------------------------------------------------------------")
  
}



# jaccard_index <- function(x, y){
#   length(intersect(x,y))/length(union(x,y))
# }

# Otsuka_Ochiai_coefficient <- function(x, y){
#     length(intersect(x, y))/sqrt(length(x) *length(y))
# }

# maxoverlap <- function(x, y){
#     max(length(intersect(x, y))/length(x), length(intersect(x, y))/length(y))
# }



# ITS_reduncancy_detection <- function(cancer,
#                                     tumor_subtype,
#                                     cor_method = "spearman",
#                                     num_gene = 200,
#                                     s = 0.4,
#                                     ITSproof = 2,
#                                     g_vote_threshold = 0.5,
#                                     ITSpath,
#                                     immuneSigInfoPath,
#                                     savepath){

#     cancer = "SKCM"
#     tumor_subtype = "Metastatic"
#     # cancer = "LUAD"
#     # tumor_subtype = "Primary"
#     purity_method = "TUMERIC"
#     cor_method = "spearman"
#     num_gene = 200
#     s=0.4
#     ITSproof = 2
#     g_vote_threshold=0.5
#     ITSpath = paste0("03_drug_immuneSig_enrichment/data/gene_immune_sig_TUMERIC_spearman_r",s , "_ITSproof_",ITSproof,"/for_GSVA/pcor/")
#     immuneSigInfoPath = "06_results_summaryandplotting/data/immene_cell_hieracy/"

#     dir.create(paste0("06_results_summaryandplotting/data/immunesig_jaccard_index_spearman_r",s , "_ITSproof_",ITSproof,"/"))
#     dir.create(paste0("06_results_summaryandplotting/data/immunesig_jaccard_index_spearman_r",s , "_ITSproof_",ITSproof,"/", purity_method))
#     dir.create(paste0("06_results_summaryandplotting/data/immunesig_jaccard_index_spearman_r",s , "_ITSproof_",ITSproof,"/", 
#         purity_method,"/",cancer,"_", tumor_subtype, "/"))
#      dir.create(paste0("06_results_summaryandplotting/data/immunesig_jaccard_index_spearman_r",s , "_ITSproof_",ITSproof,"/", 
#         purity_method,"/",cancer,"_", tumor_subtype, "/TMEcellSig/"))
#      dir.create(paste0("06_results_summaryandplotting/data/immunesig_jaccard_index_spearman_r",s , "_ITSproof_",ITSproof,"/", 
#         purity_method,"/",cancer,"_", tumor_subtype, "/immuneSig/"))

#     dir.create(paste0("06_results_summaryandplotting/data/immunesig_overlapGenes_r",s , "_ITSproof_",ITSproof,"/"))
#     dir.create(paste0("06_results_summaryandplotting/data/immunesig_overlapGenes_r",s , "_ITSproof_",ITSproof,"/", purity_method))
#     dir.create(paste0("06_results_summaryandplotting/data/immunesig_overlapGenes_r",s , "_ITSproof_",ITSproof,"//", 
#         purity_method,"/",cancer,"_", tumor_subtype, "/"))
#     dir.create(paste0("06_results_summaryandplotting/data/immunesig_overlapGenes_r",s , "_ITSproof_",ITSproof,"//", 
#         purity_method,"/",cancer,"_", tumor_subtype, "/TMEcellSig/"))
#     dir.create(paste0("06_results_summaryandplotting/data/immunesig_overlapGenes_r",s , "_ITSproof_",ITSproof,"//", 
#         purity_method,"/",cancer,"_", tumor_subtype, "/immuneSig/"))

#     library(readxl)
#     library(dplyr)
#     load(paste0(ITSpath, cor_method, "_",cancer,"_", tumor_subtype, "_positive_", num_gene, ".Rdata"))
#     # genelist_p = genelist_p[which(lapply(genelist_p, length) == 0)]

#     load(paste0(ITSpath, cor_method, "_",cancer,"_", tumor_subtype, "_negative_", num_gene, ".Rdata"))
#     # genelist_n = genelist_n[which(lapply(genelist_n, length) == 0)]


#     immunesigInfo <- as.data.frame(read_excel(paste0(immuneSigInfoPath, "immune_sig_hieracy_new.xlsx"), sheet = 1))
#     immunesigInfo <- immunesigInfo[c(1,4,5:9)]
#     immunesigInfo$immune_sig = tolower(immunesigInfo$immune_sig)
#     immunesigInfo = immunesig_name_unify(immunesigInfo)
#     immunesigInfoSelected = immunesigInfo[is.element(immunesigInfo$immune_sig, names(genelist_p)), ]

#     TMEcellSig = immunesigInfoSelected[immunesigInfoSelected$Classification == "TME",]
#     immuneSig = immunesigInfoSelected[immunesigInfoSelected$Classification == "Immune signatures",]
#     icb_predictor = immunesigInfoSelected[immunesigInfoSelected$Classification == "ICB predictor",]
#     ICBresisSig = immunesigInfoSelected[immunesigInfoSelected$Classification == "ICB resistant signatures",]
#     TMESig = immunesigInfoSelected[immunesigInfoSelected$Classification == "TMEsignatures",]


#     # LOAD original immune signature gene sets ---------------------------------
#     immunesigGS1 <- as.data.frame(read_excel(paste0(immuneSigGSpath, "immunesig_original_used.xlsx"), sheet = 1))
#     immunesigGS1 <- immunesigGS1[-grep("various", immunesigGS1$Genes), ]
#     immunesigGS1 <- immunesigGS1[-grep("SignatureMatrix", immunesigGS1$Genes), ]
#     immunesigGS1 <- immunesigGS1[-grep("ReferenceCellExpression", immunesigGS1$Genes), ]
#     immunesigGS1 <- immunesigGS1[-grep("mergedSignatures", immunesigGS1$Genes), ]
#     immunesigGS1 <- immunesigGS1[-grep("L22Matrix", immunesigGS1$Genes), ]

#     immunesigGS1 <- immunesigGS1[-c(2:5)]
#     immunesigGS1Name <- as.list(immunesigGS1$immune_sig)
#     immunesigGS1 <- lapply(immunesigGS1Name, function(x){
#       y=immunesigGS1[immunesigGS1$immune_sig == x,][-1]
#       unlist(y[-which(is.na(y))])
#     })
#     names(immunesigGS1) = unlist(immunesigGS1Name)

#     immunesigGS2 <- read_excel(paste0(immuneSigGSpath, "immunesig_original_used.xlsx"), sheet = 2)
#     immunesigGS2 <- as.data.frame(immunesigGS2[immunesigGS2$CancerType == cancer, ])
#     immunesigGS2Name <- as.list(immunesigGS2$Signatures)
#     immunesigGS2 <- lapply(immunesigGS2Name, function(x){
#       y=immunesigGS2[immunesigGS2$Signatures == x,][-1]
#       unlist(y[-which(is.na(y))])
#     })

#     names(immunesigGS2) = paste0(unlist(immunesigGS2Name),"_ConsensusTMEgsva")
#     immunesigGS = c(immunesigGS1, immunesigGS2)
#     immunesigGSname = data.frame(signame = names(immunesigGS), immune_sig = tolower(names(immunesigGS)))
#     immunesigGSname = immunesig_name_unify(immunesigGSname)
#     names(immunesigGS) = immunesigGSname$immune_sig



#     # --------------------------------------------------------------------------
#     # TMEcellSig
#     # --------------------------------------------------------------------------
#     if(length(intersect(names(genelist_p), TMEcellSig$immune_sig)) > 1){

#     ITS_TMEcellSig_p = genelist_p[intersect(names(genelist_p), TMEcellSig$immune_sig)]
#     num_gene_ITS_TMEcellSig_p = lapply(ITS_TMEcellSig_p,length)

#     TMEcellSigIndex = data.frame(table(TMEcellSig$Index))
#     # computing Jaccard Index --------------------------------------------------
#     ITS_TMEcellSig_p_jcindex <- list()
#     for(i in 1:length(ITS_TMEcellSig_p)){
#         x = ITS_TMEcellSig_p[[i]]
#         ITS_TMEcellSig_p_jcindex[[i]] = sapply(ITS_TMEcellSig_p, function(y){
#         jaccard_index(x,y)
#         # length(intersect(x, y))
#         })
#     }
#     names(ITS_TMEcellSig_p_jcindex) = names(ITS_TMEcellSig_p)

#     # plot JC index to see the pattern -----------------------------------------
#     ITS_TMEcellSig_p_jcindex = do.call(cbind, ITS_TMEcellSig_p_jcindex)
#     ITS_TMEcellSig_p_jcindex = ITS_TMEcellSig_p_jcindex[sort(colnames(ITS_TMEcellSig_p_jcindex)), ]
#     ITS_TMEcellSig_p_jcindex = ITS_TMEcellSig_p_jcindex[, sort(colnames(ITS_TMEcellSig_p_jcindex))]
    
#     # draw heatmap to see the pattern before merge
#     rownames(TMEcellSig) <- TMEcellSig$immune_sig
#     pdf(paste0("06_results_summaryandplotting/data/immunesig_jaccard_index_spearman_r",s , "_ITSproof_",ITSproof,"/", 
#         purity_method,"/",cancer,"_", tumor_subtype, "/TMEcellSig/TMEcellSig_jci_p.pdf"), 
#         width=20, height = 18)
#     ng = num_gene_ITS_TMEcellSig_p[colnames(ITS_TMEcellSig_p_jcindex)]
#     # rownames(ITS_TMEcellSig_p_jcindex) = paste0(rownames(ITS_TMEcellSig_p_jcindex), " (", unlist(ng), ")")

#     pheatmap::pheatmap(ITS_TMEcellSig_p_jcindex,
#                         cluster_cols = T,
#                         cluster_rows = T,
#                         display_numbers = T,
#                         number_format = "%.2f",
#                         fontsize = 5, 
#                         annotation_row = TMEcellSig[c(3,4)],
#                         color = colorRampPalette(c("steelblue", "white", "firebrick3"))(20),
#                         border_color = 'grey50' )
#     dev.off()

#     ITS_TMEcellSig_p_jcindex = lapply(as.list(TMEcellSigIndex[TMEcellSigIndex$Freq >= 2,]$Var1), function(x){
#       tmp = TMEcellSig[which(TMEcellSig$Index == x), ]
#       df = ITS_TMEcellSig_p_jcindex[tmp$immune_sig, tmp$immune_sig]
#       ng = num_gene_ITS_TMEcellSig_p[tmp$immune_sig]
#       rownames(df) = paste0(rownames(df), " (", unlist(ng), ")")
#       pdf(paste0("06_results_summaryandplotting/data/immunesig_jaccard_index_spearman_r",s , "_ITSproof_",ITSproof,"/", 
#           purity_method,"/",cancer,"_", tumor_subtype, "/TMEcellSig/", x , "_jci_p.pdf"), 
#           width=8, height =7)
#       pheatmap::pheatmap(df,
#                           cluster_cols = T,
#                           cluster_rows = T,
#                           display_numbers = T,
#                           number_format = "%.2f",
#                           fontsize = 10, 
#                           # annotation_row = data.frame(unlist(ng)),
#                           color = colorRampPalette(c("white", "firebrick3"))(20),
#                           border_color = 'grey50' )
#       dev.off()
#     })


#   # number of overlap genes ----------------------------------------------------
#     ITS_TMEcellSig_overlapGenes <- list()
#     for(i in 1:length(ITS_TMEcellSig_p)){
#         x = ITS_TMEcellSig_p[[i]]
#         ITS_TMEcellSig_overlapGenes[[i]] = sapply(ITS_TMEcellSig_p, function(y){
#         # jaccard_index(x,y)
#         length(intersect(x, y))
#         })
#     }


#     names(ITS_TMEcellSig_overlapGenes) = names(ITS_TMEcellSig_p)
#     ITS_TMEcellSig_overlapGenes = do.call(cbind, ITS_TMEcellSig_overlapGenes)
#     ITS_TMEcellSig_overlapGenes = ITS_TMEcellSig_overlapGenes[sort(colnames(ITS_TMEcellSig_overlapGenes)), ]
#     ITS_TMEcellSig_overlapGenes = ITS_TMEcellSig_overlapGenes[, sort(colnames(ITS_TMEcellSig_overlapGenes))]

#     pdf(paste0("06_results_summaryandplotting/data/immunesig_overlapGenes_r",s , "_ITSproof_",ITSproof,"/", 
#         purity_method,"/",cancer,"_", tumor_subtype, "/TMEcellSig/TMEcellSig_p.pdf"), 
#         width=20, height = 18)
#     # ng = num_gene_ITS_TMEcellSig_p[colnames(ITS_TMEcellSig_p_jcindex)]
#     # rownames(ITS_TMEcellSig_p_jcindex) = paste0(rownames(ITS_TMEcellSig_p_jcindex), " (", unlist(ng), ")")

#     pheatmap::pheatmap(ITS_TMEcellSig_overlapGenes,
#                         cluster_cols = T,
#                         cluster_rows = T,
#                         display_numbers = T,
#                         number_format = "%.2f",
#                         fontsize = 5, 
#                         annotation_row = TMEcellSig[c(3,4)],
#                         color = colorRampPalette(c("steelblue", "white", "firebrick3"))(20),
#                         border_color = 'grey50' )
#     dev.off()


#     ITS_TMEcellSig_p_overlapGenes = lapply(as.list(TMEcellSigIndex[TMEcellSigIndex$Freq >= 2,]$Var1), function(x){
#       tmp = TMEcellSig[which(TMEcellSig$Index == x), ]
#       df = ITS_TMEcellSig_overlapGenes[tmp$immune_sig, tmp$immune_sig]
#       ng = num_gene_ITS_TMEcellSig_p[tmp$immune_sig]
#       rownames(df) = paste0(rownames(df), " (", unlist(ng), ")")
#       pdf(paste0("06_results_summaryandplotting/data/immunesig_overlapGenes_r",s , "_ITSproof_",ITSproof,"/", 
#           purity_method,"/",cancer,"_", tumor_subtype, "/TMEcellSig/", x , "_p.pdf"), 
#           width=8, height =7)
#       pheatmap::pheatmap(df,
#                           cluster_cols = T,
#                           cluster_rows = T,
#                           display_numbers = T,
#                           number_format = "%.2f",
#                           fontsize = 10, 
#                           # annotation_row = data.frame(unlist(ng)),
#                           color = colorRampPalette(c("white", "firebrick3"))(20),
#                           border_color = 'grey50' )
#       dev.off()
#     })



#     # merge ITS in same category and high JC-index -----------------------------

#     ITS_TMEcellSig_p_merged = lapply(as.list(TMEcellSigIndex$Var1), function(x){
#       tmp = TMEcellSig[which(TMEcellSig$Index == x), ]
#       # ITS_TMEcellSig_p_jcindex[tmp$immune_sig, tmp$immune_sig]

#       gs = ITS_TMEcellSig_p[tmp$immune_sig]
#       unique(unlist(gs))
#     })
    


#     ITS_TMEcellSig_p_merged = lapply(as.list(TMEcellSigIndex[TMEcellSigIndex$Freq >= 2,]$Var1), function(x){
#       tmp = TMEcellSig[which(TMEcellSig$Index == x), ]
#       df = genelist_p[tmp$immune_sig]
#       g_all = unique(unlist(df))
#       g_vote = sapply(g_all, function(g){
#         vote = lapply(df, function(y){if(is.element(g, y))return(1)})
#         do.call(sum, vote)/length(vote)
#       })
#       g_f = names(g_vote[g_vote > g_vote_threshold])
#       return(g_f)
#     })
#     names(ITS_TMEcellSig_p_merged) = TMEcellSigIndex[TMEcellSigIndex$Freq >= 2,]$Var1
#     lapply(ITS_TMEcellSig_p_merged, length)

#     ITS_TMEcellSig_p_rest = lapply(as.list(TMEcellSigIndex[TMEcellSigIndex$Freq < 2,]$Var1), function(x){
#       tmp = TMEcellSig[which(TMEcellSig$Index == x), ]
#       df = genelist_p[[tmp$immune_sig]]
#       return(df)
#     })
#     names(ITS_TMEcellSig_p_rest) = TMEcellSigIndex[TMEcellSigIndex$Freq < 2,]$Var1

#     ITS_TMEcellSig_p_new = c(ITS_TMEcellSig_p_merged, ITS_TMEcellSig_p_rest)

#     ITS_TMEcellSig_p_new_jcindex <- list()
#     for(i in 1:length(ITS_TMEcellSig_p_new)){
#         x = ITS_TMEcellSig_p_new[[i]]
#         ITS_TMEcellSig_p_new_jcindex[[i]] = sapply(ITS_TMEcellSig_p_new, function(y){
#           jaccard_index(x,y)
#           # length(intersect(x, y))
#         })
#     }
#     names(ITS_TMEcellSig_p_new_jcindex) = names(ITS_TMEcellSig_p_new)

#     # plot JC index to see the pattern -----------------------------------------
#     ITS_TMEcellSig_p_new_jcindex = do.call(cbind, ITS_TMEcellSig_p_new_jcindex)
#     ITS_TMEcellSig_p_new_jcindex = ITS_TMEcellSig_p_new_jcindex[sort(colnames(ITS_TMEcellSig_p_new_jcindex)), ]
#     ITS_TMEcellSig_p_new_jcindex = ITS_TMEcellSig_p_new_jcindex[, sort(colnames(ITS_TMEcellSig_p_new_jcindex))]
    
#     # draw heatmap to see the pattern
#     dir.create(paste0("06_results_summaryandplotting/data/immunesig_jaccard_index_spearman_r",s , "_ITSproof_",ITSproof,"/", 
#         purity_method,"/",cancer,"_", tumor_subtype, "/TMEcellSig_merged/"))
#     pdf(paste0("06_results_summaryandplotting/data/immunesig_jaccard_index_spearman_r",s , "_ITSproof_",ITSproof,"/", 
#         purity_method,"/",cancer,"_", tumor_subtype, "/TMEcellSig_merged/TMEcellSig_jci_p.pdf"), 
#         width=15, height = 14)

#     pheatmap::pheatmap(ITS_TMEcellSig_p_new_jcindex,
#                         cluster_cols = T,
#                         cluster_rows = T,
#                         display_numbers = T,
#                         number_format = "%.2f",
#                         fontsize = 5, 
#                         # annotation_row = TMEcellSig[c(3,4)],
#                         color = colorRampPalette(c("steelblue", "white", "firebrick3"))(20),
#                         border_color = 'grey50' )
#     dev.off()

#     }




#     if(length(intersect(names(genelist_n), TMEcellSig$immune_sig)) > 1){

#     ITS_TMEcellSig_n = genelist_n[intersect(names(genelist_n), TMEcellSig$immune_sig)]
#     num_gene_ITS_TMEcellSig_n = lapply(ITS_TMEcellSig_n,length)

#     TMEcellSigIndex = data.frame(table(TMEcellSig$Index))
#     # computing Jaccard Index --------------------------------------------------
#     ITS_TMEcellSig_n_jcindex <- list()
#     for(i in 1:length(ITS_TMEcellSig_n)){
#         x = ITS_TMEcellSig_n[[i]]
#         ITS_TMEcellSig_n_jcindex[[i]] = sapply(ITS_TMEcellSig_n, function(y){
#         jaccard_index(x,y)
#         # length(intersect(x, y))
#         })
#     }
#     names(ITS_TMEcellSig_n_jcindex) = names(ITS_TMEcellSig_n)

#     # plot JC index to see the pattern -----------------------------------------
#     ITS_TMEcellSig_n_jcindex = do.call(cbind, ITS_TMEcellSig_n_jcindex)
#     ITS_TMEcellSig_n_jcindex = ITS_TMEcellSig_n_jcindex[sort(colnames(ITS_TMEcellSig_n_jcindex)), ]
#     ITS_TMEcellSig_n_jcindex = ITS_TMEcellSig_n_jcindex[, sort(colnames(ITS_TMEcellSig_n_jcindex))]
    
#     # draw heatmap to see the pattern
#     rownames(TMEcellSig) <- TMEcellSig$immune_sig
#     pdf(paste0("06_results_summaryandplotting/data/immunesig_jaccard_index_spearman_r",s , "_ITSproof_",ITSproof,"/", 
#         purity_method,"/",cancer,"_", tumor_subtype, "/TMEcellSig/TMEcellSig_jci_n.pdf"), 
#         width=20, height = 18)
#     ng = num_gene_ITS_TMEcellSig_n[colnames(ITS_TMEcellSig_n_jcindex)]
#     # rownames(ITS_TMEcellSig_n_jcindex) = paste0(rownames(ITS_TMEcellSig_n_jcindex), " (", unlist(ng), ")")

#     pheatmap::pheatmap(ITS_TMEcellSig_n_jcindex,
#                         cluster_cols = T,
#                         cluster_rows = T,
#                         display_numbers = T,
#                         number_format = "%.2f",
#                         fontsize = 5, 
#                         annotation_row = TMEcellSig[c(3,4)],
#                         color = colorRampPalette(c("steelblue", "white", "firebrick3"))(20),
#                         border_color = 'grey50' )
#     dev.off()

#     ITS_TMEcellSig_n_jcindex = lapply(as.list(TMEcellSigIndex[TMEcellSigIndex$Freq >= 2,]$Var1), function(x){
#       # 
#     #   for(i in 1:length(TMEcellSigIndex[TMEcellSigIndex$Freq >= 2,]$Var1)){

#     #   x = TMEcellSigIndex[TMEcellSigIndex$Freq >= 2,]$Var1[i]
#       tmp = TMEcellSig[which(TMEcellSig$Index == x), ]
#       df = ITS_TMEcellSig_n_jcindex[tmp$immune_sig, tmp$immune_sig]
#       ng = num_gene_ITS_TMEcellSig_n[tmp$immune_sig]
#       if(sum(unlist(ng))>1){
#       rownames(df) = paste0(rownames(df), " (", unlist(ng), ")")
#       pdf(paste0("06_results_summaryandplotting/data/immunesig_jaccard_index_spearman_r",s , "_ITSproof_",ITSproof,"/", 
#           purity_method,"/",cancer,"_", tumor_subtype, "/TMEcellSig/", x , "_jci_n.pdf"), 
#           width=8, height =7)
#       pheatmap::pheatmap(df,
#                           cluster_cols = T,
#                           cluster_rows = T,
#                           display_numbers = T,
#                           number_format = "%.2f",
#                           fontsize = 10, 
#                           # annotation_row = data.frame(unlist(ng)),
#                           color = colorRampPalette(c("white", "firebrick3"))(20),
#                           border_color = 'grey50' )
#       dev.off()

#       }

#     })


#   # number of overlap genes ----------------------------------------------------
#     ITS_TMEcellSig_overlapGenes <- list()
#     for(i in 1:length(ITS_TMEcellSig_n)){
#         x = ITS_TMEcellSig_n[[i]]
#         ITS_TMEcellSig_overlapGenes[[i]] = sapply(ITS_TMEcellSig_n, function(y){
#         # jaccard_index(x,y)
#         length(intersect(x, y))
#         })
#     }


#     names(ITS_TMEcellSig_overlapGenes) = names(ITS_TMEcellSig_n)
#     ITS_TMEcellSig_overlapGenes = do.call(cbind, ITS_TMEcellSig_overlapGenes)
#     ITS_TMEcellSig_overlapGenes = ITS_TMEcellSig_overlapGenes[sort(colnames(ITS_TMEcellSig_overlapGenes)), ]
#     ITS_TMEcellSig_overlapGenes = ITS_TMEcellSig_overlapGenes[, sort(colnames(ITS_TMEcellSig_overlapGenes))]

#     pdf(paste0("06_results_summaryandplotting/data/immunesig_overlapGenes_r",s , "_ITSproof_",ITSproof,"/", 
#         purity_method,"/",cancer,"_", tumor_subtype, "/TMEcellSig/TMEcellSig_n.pdf"), 
#         width=20, height = 18)
#     # ng = num_gene_ITS_TMEcellSig_n[colnames(ITS_TMEcellSig_n_jcindex)]
#     # rownames(ITS_TMEcellSig_n_jcindex) = paste0(rownames(ITS_TMEcellSig_n_jcindex), " (", unlist(ng), ")")

#     pheatmap::pheatmap(ITS_TMEcellSig_overlapGenes,
#                         cluster_cols = T,
#                         cluster_rows = T,
#                         display_numbers = T,
#                         number_format = "%.2f",
#                         fontsize = 5, 
#                         annotation_row = TMEcellSig[c(3,4)],
#                         color = colorRampPalette(c("steelblue", "white", "firebrick3"))(20),
#                         border_color = 'grey50' )
#     dev.off()


#     ITS_TMEcellSig_n_overlapGenes = lapply(as.list(TMEcellSigIndex[TMEcellSigIndex$Freq >= 2,]$Var1), function(x){
#       tmp = TMEcellSig[which(TMEcellSig$Index == x), ]
#       df = ITS_TMEcellSig_overlapGenes[tmp$immune_sig, tmp$immune_sig]
#       ng = num_gene_ITS_TMEcellSig_n[tmp$immune_sig]
#       if(sum(unlist(ng))>1){

#       rownames(df) = paste0(rownames(df), " (", unlist(ng), ")")
#       pdf(paste0("06_results_summaryandplotting/data/immunesig_overlapGenes_r",s , "_ITSproof_",ITSproof,"/", 
#           purity_method,"/",cancer,"_", tumor_subtype, "/TMEcellSig/", x , "_n.pdf"), 
#           width=8, height =7)
#       pheatmap::pheatmap(df,
#                           cluster_cols = T,
#                           cluster_rows = T,
#                           display_numbers = T,
#                           number_format = "%.2f",
#                           fontsize = 10, 
#                           # annotation_row = data.frame(unlist(ng)),
#                           color = colorRampPalette(c("white", "firebrick3"))(20),
#                           border_color = 'grey50' )
#       dev.off()
#       }
#     })



#     # merge ITS in same category and high JC-index -----------------------------

#     ITS_TMEcellSig_n_merged = lapply(as.list(TMEcellSigIndex$Var1), function(x){
#       tmp = TMEcellSig[which(TMEcellSig$Index == x), ]
#       # ITS_TMEcellSig_n_jcindex[tmp$immune_sig, tmp$immune_sig]

#       gs = ITS_TMEcellSig_n[tmp$immune_sig]
#       unique(unlist(gs))
#     })
    


#     ITS_TMEcellSig_n_merged = lapply(as.list(TMEcellSigIndex[TMEcellSigIndex$Freq >= 2,]$Var1), function(x){
#       tmp = TMEcellSig[which(TMEcellSig$Index == x), ]
#       df = genelist_n[tmp$immune_sig]
#       g_all = unique(unlist(df))
#       g_vote = sapply(g_all, function(g){
#         vote = lapply(df, function(y){if(is.element(g, y))return(1)})
#         do.call(sum, vote)/length(vote)
#       })
#       g_f = names(g_vote[g_vote > g_vote_threshold])
#       return(g_f)
#     })
#     names(ITS_TMEcellSig_n_merged) = TMEcellSigIndex[TMEcellSigIndex$Freq >= 2,]$Var1
#     lapply(ITS_TMEcellSig_n_merged, length)

#     ITS_TMEcellSig_n_rest = lapply(as.list(TMEcellSigIndex[TMEcellSigIndex$Freq < 2,]$Var1), function(x){
#       tmp = TMEcellSig[which(TMEcellSig$Index == x), ]
#       df = genelist_n[[tmp$immune_sig]]
#       return(df)
#     })
#     names(ITS_TMEcellSig_n_rest) = TMEcellSigIndex[TMEcellSigIndex$Freq < 2,]$Var1

#     ITS_TMEcellSig_n_new = c(ITS_TMEcellSig_n_merged, ITS_TMEcellSig_n_rest)

#     ITS_TMEcellSig_n_new_jcindex <- list()
#     for(i in 1:length(ITS_TMEcellSig_n_new)){
#         x = ITS_TMEcellSig_n_new[[i]]
#         ITS_TMEcellSig_n_new_jcindex[[i]] = sapply(ITS_TMEcellSig_n_new, function(y){
#           jaccard_index(x,y)
#           # length(intersect(x, y))
#         })
#     }
#     names(ITS_TMEcellSig_n_new_jcindex) = names(ITS_TMEcellSig_n_new)

#     # plot JC index to see the pattern -----------------------------------------
#     ITS_TMEcellSig_n_new_jcindex = do.call(cbind, ITS_TMEcellSig_n_new_jcindex)
#     ITS_TMEcellSig_n_new_jcindex = ITS_TMEcellSig_n_new_jcindex[sort(colnames(ITS_TMEcellSig_n_new_jcindex)), ]
#     ITS_TMEcellSig_n_new_jcindex = ITS_TMEcellSig_n_new_jcindex[, sort(colnames(ITS_TMEcellSig_n_new_jcindex))]
    
#     # draw heatmap to see the pattern
#     dir.create(paste0("06_results_summaryandplotting/data/immunesig_jaccard_index_spearman_r",s , "_ITSproof_",ITSproof,"/", 
#         purity_method,"/",cancer,"_", tumor_subtype, "/TMEcellSig_merged/"))
#     pdf(paste0("06_results_summaryandplotting/data/immunesig_jaccard_index_spearman_r",s , "_ITSproof_",ITSproof,"/", 
#         purity_method,"/",cancer,"_", tumor_subtype, "/TMEcellSig_merged/TMEcellSig_jci_n.pdf"), 
#         width=15, height = 14)

#     pheatmap::pheatmap(ITS_TMEcellSig_n_new_jcindex,
#                         cluster_cols = T,
#                         cluster_rows = T,
#                         display_numbers = T,
#                         number_format = "%.2f",
#                         fontsize = 5, 
#                         # annotation_row = TMEcellSig[c(3,4)],
#                         color = colorRampPalette(c("steelblue", "white", "firebrick3"))(20),
#                         border_color = 'grey50' )
#     dev.off()

#     }


#     # --------------------------------------------------------------------------
#     # --------------------------------------------------------------------------
#     # immuneSig
#     # --------------------------------------------------------------------------
#     # --------------------------------------------------------------------------
#     # detect if the original TME signatures are with similar genes -------------

#     # LOAD original immune signature gene sets ---------------------------------

#     immunesigGS_selected = immunesigGS[intersect(names(immunesigGS), immuneSig$immune_sig)]
#     # lapply(immunesigGS_selected, length)

#     TMESig_p_Ochiai <- list()
#     for(i in 1:length(immunesigGS_selected)){
#       x = immunesigGS_selected[[i]]
#       TMESig_p_Ochiai[[i]] = sapply(immunesigGS_selected, function(y){
#         # Otsuka_Ochiai_coefficient(x,y)
#         # jaccard_index(x, y)
#         maxoverlap(x, y)

#       })
#     }
#     names(TMESig_p_Ochiai) = names(immunesigGS_selected)

#     # plot JC index to see the pattern -----------------------------------------
#     TMESig_p_Ochiai = do.call(cbind, TMESig_p_Ochiai)
#     TMESig_p_Ochiai = TMESig_p_Ochiai[sort(colnames(TMESig_p_Ochiai)), ]
#     TMESig_p_Ochiai = TMESig_p_Ochiai[, sort(colnames(TMESig_p_Ochiai))]
#     TMESig_p_Ochiai[is.na(TMESig_p_Ochiai)] = 0


#     pdf(paste0("06_results_summaryandplotting/data/immunesig_jaccard_index_spearman_r",s , "_ITSproof_",ITSproof,"/oriTMESig_Ochiai.pdf"), 
#         width=18, height = 16)
#     # rownames(ITS_TMEcellSig_p_Ochiai) = paste0(rownames(ITS_TMEcellSig_p_Ochiai), " (", unlist(ng), ")")
#       num_gene_immunesigGS_selected = lapply(immunesigGS_selected,length)
#       ng = num_gene_immunesigGS_selected[colnames(TMESig_p_Ochiai)]
#       colnames(TMESig_p_Ochiai) = paste0(colnames(TMESig_p_Ochiai), " (", unlist(ng), ")")

#       tmp = pheatmap::pheatmap(TMESig_p_Ochiai,
#                         cluster_cols = T,
#                         cluster_rows = T,
#                         display_numbers = T,
#                         number_format = "%.2f",
#                         fontsize = 5, 
#                         # annotation_row = immuneSig[c(3,4)],
#                         color = colorRampPalette(c("steelblue", "white", "firebrick3"))(20),
#                         border_color = 'grey50' )
#       dev.off()


#     # d = dist(TMESig_p_Ochiai, method = 'euclidean')
#     # tree = hclust(d, method = 'complete')
#     tree = tmp$tree_col
#     treedf = as.data.frame(as.matrix(do.call(cbind, tree[c(4,2:3)])))
#     treedf$height = as.numeric(treedf$height)
#     treedf$order = as.numeric(treedf$order)
#     names(treedf) = c("immune_sig", "height","order")
#     treedf$immune_sig = as.character(treedf$immune_sig)
#     v = cutree(tree, 15)[tree$order]
#     gaps = which((v[-1] - v[-length(v)]) != 0)
#     its_cls <- as.data.frame(as.matrix(v))
#     its_cls$immune_sig <- str_split(rownames(its_cls),'[.]',simplify = T)[,1]
#     its_cls = inner_join(treedf, its_cls)
#     its_cls = its_cls[order(its_cls$V1),]
#     its_cls$V1 = as.character(its_cls$V1)
#     rownames(its_cls) = its_cls$immune_sig




    
#     if(length(intersect(names(genelist_p), immuneSig$immune_sig)) > 1){
  
#     ITS_immuneSig_p = genelist_p[intersect(names(genelist_p), immuneSig$immune_sig)]

#     # ITS_immuneSig_p_oriImmuSig_sim = lapply(as.list(intersect(names(ITS_immuneSig_p), names(immunesigGS_selected))),function(x){
#     #   data.frame(rate_on_ori=length(intersect(ITS_immuneSig_p[[x]], immunesigGS_selected[[x]]))/ length(immunesigGS_selected[[x]]),
#     #   num_ITS_p=length(ITS_immuneSig_p[[x]]), 
#     #   num_ori_p=length(immunesigGS_selected[[x]]))
      
#     #   })

#     # names(ITS_immuneSig_p_oriImmuSig_sim) = intersect(names(ITS_immuneSig_p), names(immunesigGS_selected))

#     num_gene_ITS_immuneSig_p = lapply(ITS_immuneSig_p, length)

#     immuneSigIndex = data.frame(table(immuneSig$Index))
#     # computing Jaccard Index --------------------------------------------------
#     ITS_immuneSig_p_jcindex <- list()
#     for(i in 1:length(ITS_immuneSig_p)){
#         x = ITS_immuneSig_p[[i]]
#         ITS_immuneSig_p_jcindex[[i]] = sapply(ITS_immuneSig_p, function(y){
#         jaccard_index(x,y)
#         # length(intersect(x, y))
#         })
#     }
#     names(ITS_immuneSig_p_jcindex) = names(ITS_immuneSig_p)

#     # plot JC index to see the pattern before merge ----------------------------
#     ITS_immuneSig_p_jcindex = do.call(cbind, ITS_immuneSig_p_jcindex)
#     ITS_immuneSig_p_jcindex = ITS_immuneSig_p_jcindex[sort(colnames(ITS_immuneSig_p_jcindex)), ]
#     ITS_immuneSig_p_jcindex = ITS_immuneSig_p_jcindex[, sort(colnames(ITS_immuneSig_p_jcindex))]
    
#     # draw heatmap to see the pattern
#     rownames(immuneSig) <- immuneSig$immune_sig
#     pdf(paste0("06_results_summaryandplotting/data/immunesig_jaccard_index_spearman_r", s, "_ITSproof_",ITSproof,"/", 
#         purity_method,"/",cancer,"_", tumor_subtype, "/immunesig_jci_p.pdf"), 
#         width=20, height = 18)
#     ng = num_gene_ITS_immuneSig_p[colnames(ITS_immuneSig_p_jcindex)]
#     # rownames(ITS_immuneSig_p_jcindex) = paste0(rownames(ITS_immuneSig_p_jcindex), " (", unlist(ng), ")")

#     pheatmap::pheatmap(ITS_immuneSig_p_jcindex,
#                         cluster_cols = T,
#                         cluster_rows = T,
#                         display_numbers = T,
#                         number_format = "%.2f",
#                         fontsize = 5, 
#                         annotation_row = immuneSig[c(3,4)],
#                         color = colorRampPalette(c("steelblue", "white", "firebrick3"))(20),
#                         border_color = 'grey50' )
#     dev.off()



#     # ITS_immuneSig_p_jcindex = lapply(as.list(immuneSigIndex[immuneSigIndex$Freq >= 2,]$Var1[-12]), function(x){
#     #   tmp = immuneSig[which(immuneSig$Index == x), ]
#     #   df = ITS_immuneSig_p_jcindex[tmp$immune_sig, tmp$immune_sig]
#     #   ng = num_gene_ITS_immuneSig_p[tmp$immune_sig]
#     #   rownames(df) = paste0(rownames(df), " (", unlist(ng), ")")
#     #   pdf(paste0("06_results_summaryandplotting/data/immunesig_jaccard_index_spearman_r",s , "_ITSproof_",ITSproof,"/", 
#     #       purity_method,"/",cancer,"_", tumor_subtype, "/immuneSig/", x , "_jci_p.pdf"), 
#     #       width=8, height =7)
#     #   pheatmap::pheatmap(df,
#     #                       cluster_cols = T,
#     #                       cluster_rows = T,
#     #                       display_numbers = T,
#     #                       number_format = "%.2f",
#     #                       fontsize = 10, 
#     #                       # annotation_row = data.frame(unlist(ng)),
#     #                       color = colorRampPalette(c("white", "firebrick3"))(20),
#     #                       border_color = 'grey50' )
#     #   dev.off()
#     # })

#     ITS_immuneSig_p_jcindex = lapply(as.list(unique(its_cls$V1)), function(x){
#       tmp = its_cls[which(its_cls$V1 == x), ]
#       if(nrow(tmp) > 1){

#       df = ITS_immuneSig_p_jcindex[tmp$immune_sig, tmp$immune_sig]
#       ng = num_gene_ITS_immuneSig_p[tmp$immune_sig]
#       rownames(df) = paste0(rownames(df), " (", unlist(ng), ")")
#       pdf(paste0("06_results_summaryandplotting/data/immunesig_jaccard_index_spearman_r",s , "_ITSproof_",ITSproof,"/", 
#           purity_method,"/",cancer,"_", tumor_subtype, "/immuneSig/", x , "_jci_p.pdf"), 
#           width=8, height =7)
#       pheatmap::pheatmap(df,
#                           cluster_cols = T,
#                           cluster_rows = T,
#                           display_numbers = T,
#                           number_format = "%.2f",
#                           fontsize = 10, 
#                           # annotation_row = data.frame(unlist(ng)),
#                           color = colorRampPalette(c("white", "firebrick3"))(20),
#                           border_color = 'grey50' )
#       dev.off()
#       }
#     })


#     ITS_immuneSig_overlapGenes <- list()
#     for(i in 1:length(ITS_immuneSig_p)){
#         x = ITS_immuneSig_p[[i]]
#         ITS_immuneSig_overlapGenes[[i]] = sapply(ITS_immuneSig_p, function(y){
#         # jaccard_index(x,y)
#         length(intersect(x, y))
#         })
#     }


#     names(ITS_immuneSig_overlapGenes) = names(ITS_immuneSig_p)
#     ITS_immuneSig_overlapGenes = do.call(cbind, ITS_immuneSig_overlapGenes)
#     ITS_immuneSig_overlapGenes = ITS_immuneSig_overlapGenes[sort(colnames(ITS_immuneSig_overlapGenes)), ]
#     ITS_immuneSig_overlapGenes = ITS_immuneSig_overlapGenes[, sort(colnames(ITS_immuneSig_overlapGenes))]

#     # ITS_immuneSig_p_overlapGenes = lapply(as.list(immuneSigIndex[immuneSigIndex$Freq >= 2,]$Var1[-12]), function(x){
#     #   tmp = immuneSig[which(immuneSig$Index == x), ]
#     #   df = ITS_immuneSig_overlapGenes[tmp$immune_sig, tmp$immune_sig]
#     #   ng = num_gene_ITS_immuneSig_p[tmp$immune_sig]
#     #   rownames(df) = paste0(rownames(df), " (", unlist(ng), ")")
#     #   pdf(paste0("06_results_summaryandplotting/data/immunesig_overlapGenes_r",s , "_ITSproof_",ITSproof,"/", 
#     #       purity_method,"/",cancer,"_", tumor_subtype, "/immuneSig/", x , "_p.pdf"), 
#     #       width=8, height =7)
#     #   pheatmap::pheatmap(df,
#     #                       cluster_cols = T,
#     #                       cluster_rows = T,
#     #                       display_numbers = T,
#     #                       number_format = "%.2f",
#     #                       fontsize = 10, 
#     #                       # annotation_row = data.frame(unlist(ng)),
#     #                       color = colorRampPalette(c("white", "firebrick3"))(20),
#     #                       border_color = 'grey50' )
#     #   dev.off()
#     # })

#     ITS_immuneSig_p_overlapGenes = lapply(as.list(its_cls$V1), function(x){
#       tmp = its_cls[which(its_cls$Index == x), ]
#       if(nrow(tmp) > 1){
#       df = ITS_immuneSig_overlapGenes[tmp$immune_sig, tmp$immune_sig]
#       ng = num_gene_ITS_immuneSig_p[tmp$immune_sig]
#       rownames(df) = paste0(rownames(df), " (", unlist(ng), ")")
#       pdf(paste0("06_results_summaryandplotting/data/immunesig_overlapGenes_r",s , "_ITSproof_",ITSproof,"/", 
#           purity_method,"/",cancer,"_", tumor_subtype, "/immuneSig/", x , "_p.pdf"), 
#           width=8, height =7)
#       pheatmap::pheatmap(df,
#                           cluster_cols = T,
#                           cluster_rows = T,
#                           display_numbers = T,
#                           number_format = "%.2f",
#                           fontsize = 10, 
#                           # annotation_row = data.frame(unlist(ng)),
#                           color = colorRampPalette(c("white", "firebrick3"))(20),
#                           border_color = 'grey50' )
#       dev.off()
#       }
#     })


#     # merge ITS in same category and high JC-index -----------------------------

#     # ITS_immuneSig_p_merged = lapply(as.list(immuneSigIndex[immuneSigIndex$Freq >= 2,]$Var1), function(x){
#     #   tmp = immuneSig[which(immuneSig$Index == x), ]
#     #   df = genelist_p[tmp$immune_sig]
#     #   g_all = unique(unlist(df))
#     #   g_vote = sapply(g_all, function(g){
#     #     vote = lapply(df, function(y){if(is.element(g, y))return(1)})
#     #     do.call(sum, vote)/length(vote)
#     #   })
#     #   g_f = names(g_vote[g_vote > g_vote_threshold])
#     #   return(g_f)
#     # })
#     # names(ITS_immuneSig_p_merged) = immuneSigIndex[immuneSigIndex$Freq >= 2,]$Var1
#     # lapply(ITS_immuneSig_p_merged, length)

#     # ITS_immuneSig_p_rest = lapply(as.list(immuneSigIndex[immuneSigIndex$Freq < 2,]$Var1), function(x){
#     #   tmp = immuneSig[which(immuneSig$Index == x), ]
#     #   df = genelist_p[[tmp$immune_sig]]
#     #   return(df)
#     # })
#     # names(ITS_immuneSig_p_rest) = immuneSigIndex[immuneSigIndex$Freq < 2,]$Var1

#     # ITS_immuneSig_p_new = c(ITS_immuneSig_p_merged, ITS_immuneSig_p_rest)


#     merge_index = data.frame(table(its_cls$V1))
#     merge_index1 = merge_index[merge_index$Var1!=1 & merge_index$Freq>1,]
#     ITS_immuneSig_p_merged = lapply(as.list(merge_index1$Var1), function(x){
#       tmp = its_cls[which(its_cls$V1 == x), ]
#       df = genelist_p[tmp$immune_sig]
#       g_all = unique(unlist(df))
#       g_vote = sapply(g_all, function(g){
#         vote = lapply(df, function(y){if(is.element(g, y))return(1)})
#         do.call(sum, vote)/length(vote)
#       })
#       g_f = names(g_vote[g_vote > g_vote_threshold])
#       return(g_f)
#     })
#     names(ITS_immuneSig_p_merged) = merge_index1$Var1
#     lapply(ITS_immuneSig_p_merged, length)

#     merge_index2 = merge_index[merge_index$Var1 ==1 | merge_index$Freq == 1,]
#     immunesig_rest = its_cls[is.element(its_cls$V1, merge_index2$Var1), ]
#     ITS_immuneSig_p_rest = lapply(as.list(immunesig_rest$immune_sig), function(x){
#       # tmp = immuneSig[which(immuneSig$Index == x), ]
#       df = genelist_p[[x]]
#       return(df)
#     })
#     names(ITS_immuneSig_p_rest) = immunesig_rest$immune_sig
#     ITS_immuneSig_p_new = c(ITS_immuneSig_p_merged, ITS_immuneSig_p_rest)

#     ITS_immuneSig_p_new_jcindex <- list()
#     for(i in 1:length(ITS_immuneSig_p_new)){
#         x = ITS_immuneSig_p_new[[i]]
#         ITS_immuneSig_p_new_jcindex[[i]] = sapply(ITS_immuneSig_p_new, function(y){
#           # jaccard_index(x,y)
#           maxoverlap(x, y)
#         })
#     }
#     names(ITS_immuneSig_p_new_jcindex) = names(ITS_immuneSig_p_new)

#     # plot JC index to see the pattern -----------------------------------------
#     ITS_immuneSig_p_new_jcindex = do.call(cbind, ITS_immuneSig_p_new_jcindex)
#     ITS_immuneSig_p_new_jcindex = ITS_immuneSig_p_new_jcindex[sort(colnames(ITS_immuneSig_p_new_jcindex)), ]
#     ITS_immuneSig_p_new_jcindex = ITS_immuneSig_p_new_jcindex[, sort(colnames(ITS_immuneSig_p_new_jcindex))]
    
#     # draw heatmap to see the pattern
#     dir.create(paste0("06_results_summaryandplotting/data/immunesig_jaccard_index_spearman_r",s , "_ITSproof_",ITSproof,"/", 
#         purity_method,"/",cancer,"_", tumor_subtype, "/immuneSig_merged/"))
#     pdf(paste0("06_results_summaryandplotting/data/immunesig_jaccard_index_spearman_r",s , "_ITSproof_",ITSproof,"/", 
#         purity_method,"/",cancer,"_", tumor_subtype, "/immuneSig_merged/immuneSig_jci_p_new.pdf"), 
#         width=12, height =12)

#     pheatmap::pheatmap(ITS_immuneSig_p_new_jcindex,
#                         cluster_cols = T,
#                         cluster_rows = T,
#                         display_numbers = T,
#                         number_format = "%.2f",
#                         fontsize = 5, 
#                         # annotation_row = immuneSig[c(3,4)],
#                         color = colorRampPalette(c("steelblue", "white", "firebrick3"))(20),
#                         border_color = 'grey50' )
#     dev.off()


#     # # merge ITS in same category and high JC-index -----------------------------

#     # ITS_immuneSig_p_merged = lapply(as.list(immuneSigIndex$Var1), function(x){
#     #   tmp = immuneSig[which(immuneSig$Index == x), ]
#     #   # ITS_immuneSig_p_jcindex[tmp$immune_sig, tmp$immune_sig]

#     #   gs = ITS_immuneSig_p[tmp$immune_sig]
#     #   unique(unlist(gs))
#     # })
    


#     # ITS_immuneSig_p_merged = lapply(as.list(immuneSigIndex[immuneSigIndex$Freq >= 2,]$Var1), function(x){
#     #   tmp = immuneSig[which(immuneSig$Index == x), ]
#     #   df = genelist_p[tmp$immune_sig]
#     #   g_all = unique(unlist(df))
#     #   g_vote = sapply(g_all, function(g){
#     #     vote = lapply(df, function(y){if(is.element(g, y))return(1)})
#     #     do.call(sum, vote)/length(vote)
#     #   })
#     #   g_f = names(g_vote[g_vote > g_vote_threshold])
#     #   return(g_f)
#     # })
#     # names(ITS_immuneSig_p_merged) = immuneSigIndex[immuneSigIndex$Freq >= 2,]$Var1
#     # lapply(ITS_immuneSig_p_merged, length)

#     # ITS_immuneSig_p_rest = lapply(as.list(immuneSigIndex[immuneSigIndex$Freq < 2,]$Var1), function(x){
#     #   tmp = immuneSig[which(immuneSig$Index == x), ]
#     #   df = genelist_p[[tmp$immune_sig]]
#     #   return(df)
#     # })
#     # names(ITS_immuneSig_p_rest) = immuneSigIndex[immuneSigIndex$Freq < 2,]$Var1

#     # ITS_immuneSig_p_new = c(ITS_immuneSig_p_merged, ITS_immuneSig_p_rest)

#     # ITS_immuneSig_p_new_jcindex <- list()
#     # for(i in 1:length(ITS_immuneSig_p_new)){
#     #     x = ITS_immuneSig_p_new[[i]]
#     #     ITS_immuneSig_p_new_jcindex[[i]] = sapply(ITS_immuneSig_p_new, function(y){
#     #       jaccard_index(x,y)
#     #       # length(intersect(x, y))
#     #     })
#     # }
#     # names(ITS_immuneSig_p_new_jcindex) = names(ITS_immuneSig_p_new)

#     # # plot JC index to see the pattern -----------------------------------------
#     # ITS_immuneSig_p_new_jcindex = do.call(cbind, ITS_immuneSig_p_new_jcindex)
#     # ITS_immuneSig_p_new_jcindex = ITS_immuneSig_p_new_jcindex[sort(colnames(ITS_immuneSig_p_new_jcindex)), ]
#     # ITS_immuneSig_p_new_jcindex = ITS_immuneSig_p_new_jcindex[, sort(colnames(ITS_immuneSig_p_new_jcindex))]
    
#     # # draw heatmap to see the pattern
#     # dir.create(paste0("06_results_summaryandplotting/data/immunesig_jaccard_index_spearman_r",s , "_ITSproof_",ITSproof,"/", 
#     #     purity_method,"/",cancer,"_", tumor_subtype, "/immuneSig_merged/"))
#     # pdf(paste0("06_results_summaryandplotting/data/immunesig_jaccard_index_spearman_r",s , "_ITSproof_",ITSproof,"/", 
#     #     purity_method,"/",cancer,"_", tumor_subtype, "/immuneSig_merged/immuneSig_jci_p.pdf"), 
#     #     width=15, height = 14)

#     # pheatmap::pheatmap(ITS_immuneSig_p_new_jcindex,
#     #                     cluster_cols = T,
#     #                     cluster_rows = T,
#     #                     display_numbers = T,
#     #                     number_format = "%.2f",
#     #                     fontsize = 5, 
#     #                     # annotation_row = immuneSig[c(3,4)],
#     #                     color = colorRampPalette(c("steelblue", "white", "firebrick3"))(20),
#     #                     border_color = 'grey50' )
#     # dev.off()

#     }


#   if(length(intersect(names(genelist_n), immuneSig$immune_sig)) > 1){

#     ITS_immuneSig_n = genelist_n[intersect(names(genelist_n), immuneSig$immune_sig)]

#     ITS_immuneSig_n_oriImmuSig_sim = lapply(as.list(intersect(names(ITS_immuneSig_n), names(immunesigGS_selected))),function(x){
#       length(intersect(ITS_immuneSig_n[[x]], immunesigGS_selected[[x]]))/ length(immunesigGS_selected[[x]])})
#     names(ITS_immuneSig_n_oriImmuSig_sim) = intersect(names(ITS_immuneSig_n), names(immunesigGS_selected))

#     num_gene_ITS_immuneSig_n = lapply(ITS_immuneSig_n,length)

#     immuneSigIndex = data.frame(table(immuneSig$Index))
#     # computing Jaccard Index --------------------------------------------------
#     ITS_immuneSig_n_jcindex <- list()
#     for(i in 1:length(ITS_immuneSig_n)){
#         x = ITS_immuneSig_n[[i]]
#         ITS_immuneSig_n_jcindex[[i]] = sapply(ITS_immuneSig_n, function(y){
#         jaccard_index(x,y)
#         # length(intersect(x, y))
#         })
#     }
#     names(ITS_immuneSig_n_jcindex) = names(ITS_immuneSig_n)

#     # plot JC index to see the pattern -----------------------------------------
#     ITS_immuneSig_n_jcindex = do.call(cbind, ITS_immuneSig_n_jcindex)
#     ITS_immuneSig_n_jcindex = ITS_immuneSig_n_jcindex[sort(colnames(ITS_immuneSig_n_jcindex)), ]
#     ITS_immuneSig_n_jcindex = ITS_immuneSig_n_jcindex[, sort(colnames(ITS_immuneSig_n_jcindex))]
    
#     # draw heatmap to see the pattern
#     rownames(immuneSig) <- immuneSig$immune_sig
#     pdf(paste0("06_results_summaryandplotting/data/immunesig_jaccard_index_spearman_r",s , "_ITSproof_",ITSproof,"/", 
#         purity_method,"/",cancer,"_", tumor_subtype, "/immuneSig/immuneSig_jci_n.pdf"), 
#         width=20, height = 18)
#     ng = num_gene_ITS_immuneSig_n[colnames(ITS_immuneSig_n_jcindex)]
#     # rownames(ITS_immuneSig_n_jcindex) = paste0(rownames(ITS_immuneSig_n_jcindex), " (", unlist(ng), ")")

#     pheatmap::pheatmap(ITS_immuneSig_n_jcindex,
#                         cluster_cols = T,
#                         cluster_rows = T,
#                         display_numbers = T,
#                         number_format = "%.2f",
#                         fontsize = 5, 
#                         annotation_row = immuneSig[c(3,4)],
#                         color = colorRampPalette(c("steelblue", "white", "firebrick3"))(20),
#                         border_color = 'grey50' )
#     dev.off()

#     ITS_immuneSig_n_jcindex = lapply(as.list(immuneSigIndex[immuneSigIndex$Freq >= 2,]$Var1), function(x){
#       # 
#     #   for(i in 1:length(immuneSigIndex[immuneSigIndex$Freq >= 2,]$Var1)){

#     #   x = immuneSigIndex[immuneSigIndex$Freq >= 2,]$Var1[i]
#       tmp = immuneSig[which(immuneSig$Index == x), ]
#       df = ITS_immuneSig_n_jcindex[tmp$immune_sig, tmp$immune_sig]
#       ng = num_gene_ITS_immuneSig_n[tmp$immune_sig]
#       if(sum(unlist(ng))>2){
#       rownames(df) = paste0(rownames(df), " (", unlist(ng), ")")
#       pdf(paste0("06_results_summaryandplotting/data/immunesig_jaccard_index_spearman_r",s , "_ITSproof_",ITSproof,"/", 
#           purity_method,"/",cancer,"_", tumor_subtype, "/immuneSig/", x , "_jci_n.pdf"), 
#           width=8, height =7)
#       pheatmap::pheatmap(df,
#                           cluster_cols = T,
#                           cluster_rows = T,
#                           display_numbers = T,
#                           number_format = "%.2f",
#                           fontsize = 10, 
#                           # annotation_row = data.frame(unlist(ng)),
#                           color = colorRampPalette(c("white", "firebrick3"))(20),
#                           border_color = 'grey50' )
#       dev.off()

#       }

#     })


#   # number of overlap genes ----------------------------------------------------
#     ITS_immuneSig_overlapGenes <- list()
#     for(i in 1:length(ITS_immuneSig_n)){
#         x = ITS_immuneSig_n[[i]]
#         ITS_immuneSig_overlapGenes[[i]] = sapply(ITS_immuneSig_n, function(y){
#         # jaccard_index(x,y)
#         length(intersect(x, y))
#         })
#     }


#     names(ITS_immuneSig_overlapGenes) = names(ITS_immuneSig_n)
#     ITS_immuneSig_overlapGenes = do.call(cbind, ITS_immuneSig_overlapGenes)
#     ITS_immuneSig_overlapGenes = ITS_immuneSig_overlapGenes[sort(colnames(ITS_immuneSig_overlapGenes)), ]
#     ITS_immuneSig_overlapGenes = ITS_immuneSig_overlapGenes[, sort(colnames(ITS_immuneSig_overlapGenes))]

#     pdf(paste0("06_results_summaryandplotting/data/immunesig_overlapGenes_r",s , "_ITSproof_",ITSproof,"/", 
#         purity_method,"/",cancer,"_", tumor_subtype, "/immuneSig/immuneSig_n.pdf"), 
#         width=20, height = 18)
#     # ng = num_gene_ITS_immuneSig_n[colnames(ITS_immuneSig_n_jcindex)]
#     # rownames(ITS_immuneSig_n_jcindex) = paste0(rownames(ITS_immuneSig_n_jcindex), " (", unlist(ng), ")")

#     pheatmap::pheatmap(ITS_immuneSig_overlapGenes,
#                         cluster_cols = T,
#                         cluster_rows = T,
#                         display_numbers = T,
#                         number_format = "%.2f",
#                         fontsize = 5, 
#                         annotation_row = immuneSig[c(3,4)],
#                         color = colorRampPalette(c("steelblue", "white", "firebrick3"))(20),
#                         border_color = 'grey50' )
#     dev.off()


#     ITS_immuneSig_n_overlapGenes = lapply(as.list(immuneSigIndex[immuneSigIndex$Freq >= 2,]$Var1), function(x){
#       tmp = immuneSig[which(immuneSig$Index == x), ]
#       df = ITS_immuneSig_overlapGenes[tmp$immune_sig, tmp$immune_sig]
#       ng = num_gene_ITS_immuneSig_n[tmp$immune_sig]
#       if(sum(unlist(ng))>2){

#       rownames(df) = paste0(rownames(df), " (", unlist(ng), ")")
#       pdf(paste0("06_results_summaryandplotting/data/immunesig_overlapGenes_r",s , "_ITSproof_",ITSproof,"/", 
#           purity_method,"/",cancer,"_", tumor_subtype, "/immuneSig/", x , "_n.pdf"), 
#           width=8, height =7)
#       pheatmap::pheatmap(df,
#                           cluster_cols = T,
#                           cluster_rows = T,
#                           display_numbers = T,
#                           number_format = "%.2f",
#                           fontsize = 10, 
#                           # annotation_row = data.frame(unlist(ng)),
#                           color = colorRampPalette(c("white", "firebrick3"))(20),
#                           border_color = 'grey50' )
#       dev.off()
#       }
#     })



#     # merge ITS in same category and high JC-index -----------------------------

#     ITS_immuneSig_n_merged = lapply(as.list(immuneSigIndex$Var1), function(x){
#       tmp = immuneSig[which(immuneSig$Index == x), ]
#       # ITS_immuneSig_n_jcindex[tmp$immune_sig, tmp$immune_sig]

#       gs = ITS_immuneSig_n[tmp$immune_sig]
#       unique(unlist(gs))
#     })
    


#     ITS_immuneSig_n_merged = lapply(as.list(immuneSigIndex[immuneSigIndex$Freq >= 2,]$Var1), function(x){
#       tmp = immuneSig[which(immuneSig$Index == x), ]
#       df = genelist_n[tmp$immune_sig]
#       g_all = unique(unlist(df))
#       g_vote = sapply(g_all, function(g){
#         vote = lapply(df, function(y){if(is.element(g, y))return(1)})
#         do.call(sum, vote)/length(vote)
#       })
#       g_f = names(g_vote[g_vote > g_vote_threshold])
#       return(g_f)
#     })
#     names(ITS_immuneSig_n_merged) = immuneSigIndex[immuneSigIndex$Freq >= 2,]$Var1
#     lapply(ITS_immuneSig_n_merged, length)

#     ITS_immuneSig_n_rest = lapply(as.list(immuneSigIndex[immuneSigIndex$Freq < 2,]$Var1), function(x){
#       tmp = immuneSig[which(immuneSig$Index == x), ]
#       df = genelist_n[[tmp$immune_sig]]
#       return(df)
#     })
#     names(ITS_immuneSig_n_rest) = immuneSigIndex[immuneSigIndex$Freq < 2,]$Var1

#     ITS_immuneSig_n_new = c(ITS_immuneSig_n_merged, ITS_immuneSig_n_rest)

#     ITS_immuneSig_n_new_jcindex <- list()
#     for(i in 1:length(ITS_immuneSig_n_new)){
#         x = ITS_immuneSig_n_new[[i]]
#         ITS_immuneSig_n_new_jcindex[[i]] = sapply(ITS_immuneSig_n_new, function(y){
#           jaccard_index(x,y)
#           # length(intersect(x, y))
#         })
#     }
#     names(ITS_immuneSig_n_new_jcindex) = names(ITS_immuneSig_n_new)

#     # plot JC index to see the pattern -----------------------------------------
#     ITS_immuneSig_n_new_jcindex = do.call(cbind, ITS_immuneSig_n_new_jcindex)
#     ITS_immuneSig_n_new_jcindex = ITS_immuneSig_n_new_jcindex[sort(colnames(ITS_immuneSig_n_new_jcindex)), ]
#     ITS_immuneSig_n_new_jcindex = ITS_immuneSig_n_new_jcindex[, sort(colnames(ITS_immuneSig_n_new_jcindex))]
    
#     # draw heatmap to see the pattern
#     dir.create(paste0("06_results_summaryandplotting/data/immunesig_jaccard_index_spearman_r",s , "_ITSproof_",ITSproof,"/", 
#         purity_method,"/",cancer,"_", tumor_subtype, "/immuneSig_merged/"))
#     pdf(paste0("06_results_summaryandplotting/data/immunesig_jaccard_index_spearman_r",s , "_ITSproof_",ITSproof,"/", 
#         purity_method,"/",cancer,"_", tumor_subtype, "/immuneSig_merged/immuneSig_jci_n.pdf"), 
#         width=15, height = 14)

#     pheatmap::pheatmap(ITS_immuneSig_n_new_jcindex,
#                         cluster_cols = T,
#                         cluster_rows = T,
#                         display_numbers = T,
#                         number_format = "%.2f",
#                         fontsize = 5, 
#                         # annotation_row = immuneSig[c(3,4)],
#                         color = colorRampPalette(c("steelblue", "white", "firebrick3"))(20),
#                         border_color = 'grey50' )
#     dev.off()

#     }


#  # ICB predictor ---------------------------------------------------------------
#     if(length(intersect(names(genelist_p), icb_predictor$immune_sig)) > 1){

#       ITS_icb_predictor_p = genelist_p[intersect(names(genelist_p), icb_predictor$immune_sig)]
    
#     }

    


#  # ICB resistant signatures ---------------------------------------------------------------
#     if(length(intersect(names(genelist_p), ICBresisSig$immune_sig)) > 1){

#       ITS_ICBresisSig_p = genelist_p[intersect(names(genelist_p), ICBresisSig$immune_sig)]

#     }


#  # ICB resistant signatures ---------------------------------------------------------------
#     if(length(intersect(names(genelist_p), TMESig$immune_sig)) > 1){

#       ITS_ICBresisSig_p = genelist_p[intersect(names(genelist_p), TMESig$immune_sig)]

#     }


#   # generate gmt files for GSEA analysis ---------------------------------------- 
#   print(" save ITS for GSEA analysis.")
  
#   genelist2 <- sapply(genelist, function(x)paste(x, collapse = "\t"))
#   genelist_p2 <- sapply(genelist_p, function(x)paste(x, collapse = "\t"))
#   genelist_n2 <- sapply(genelist_n, function(x)paste(x, collapse = "\t"))
  
#   genelist2 <- cbind(matrix(NA, length(genelist2)), genelist2)
#   genelist_p2 <- cbind(matrix(NA, length(genelist_p2)), genelist_p2)
#   genelist_n2 <- cbind(matrix(NA, length(genelist_n2)), genelist_n2)
  
#   if(nrow(genelist_p2) != 0){
#     write.table(genelist_p2, paste0(savepath, "for_GSEA/cor/", method, "_", cancer, "_",sample_type, "_positive_", num_gene_p, ".gmt"), sep='\t', row.names=T, col.names = F, quote=F)
#   }
  
#   if(nrow(genelist_n2) != 0){ 
#     write.table(genelist_n2, paste0(savepath, "for_GSEA/cor/", method, "_", cancer, "_",sample_type, "_negative_", num_gene_n, ".gmt"), sep='\t', row.names=T, col.names = F, quote=F)
#   }

# }





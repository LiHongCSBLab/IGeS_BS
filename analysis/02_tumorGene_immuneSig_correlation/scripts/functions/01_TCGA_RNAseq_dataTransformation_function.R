# function used for expression precessing


# GeneID2Name -------------------------------------------------------------
# Convert Gene Ensemble ID to HGNC Symbol with biomaRt
GeneID2Name <- function(x,
                        species = c('human', 'mouse')) {
  if (!require(biomaRt)) {
    BiocManager::install("biomaRt")
  }
  require(biomaRt)
  
  if (species == 'human') {
    human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    genelist <-
      getBM(
        attributes = c("ensembl_gene_id", "hgnc_symbol"),
        filters = "ensembl_gene_id",
        values = x,
        mart = human,
        useCache = FALSE
      )
    
    
  } else if (species == 'mouse') {
    mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
    genelist <-
      getBM(
        attributes = c("ensembl_gene_id", "hgnc_symbol"),
        filters = "ensembl_gene_id",
        values = x,
        mart = mouse,
        useCache = FALSE
      )
  }
  return(genelist)
}



# get gene length --------------------------------------------------------------------
# This exon-based gene length: drop overlapped exon, keep nun-redundance gene length
# Output includes: Gene ensemble ID, gene name and gene length

geneLen <- function(v = c('hg19', 'hg38')) {
 if (!require(data.table)) {
    install.packages("data.table")
  }
  if (!require(rtracklayer)) {
    BiocManager::install("rtracklayer")
  }
  if (!require(IRanges)) {
    BiocManager::install("IRanges")
  }
  
  require(data.table)
  require(IRanges)
  require(rtracklayer)
  
  if (v == 'hg19') {
    hg19 <-
      readGFF("/picb/bigdata/project/FengFYM/s7_annotation_files/hg19.gencodev27.gtf")
    anno <- setDT(hg19)
  } else if (v == 'hg38') {
    hg38 <-
      readGFF("/picb/bigdata/project/FengFYM/s7_annotation_files/gencode.v22.annotation.gtf")
    anno <- setDT(hg38)
  }
  
  anno <- anno[type == "exon", ]
  setnames(
    anno,
    c("seqid", "start", "end", "gene_id", "gene_name", "exon_number"),
    c("Chr", "ExonStart", "ExonEnd","Ensemble_id",  "Gene", "Exon_number")
  )
  #mkdir bin and mean by bin
  Exon_region <-
    unique(anno[, .(Chr, ExonStart, ExonEnd, Exon_number, Ensemble_id)])
  Exon_region <-
    Exon_region[, {
      x <-
        IRanges(ExonStart, ExonEnd)
      y <-
        reduce(x)
      list(ExonStart = y@start,
           ExonEnd = y@start + y@width - 1)
    }, by = .(Ensemble_id, Chr)]
  Exon_region[, Exon_num := 1:.N, by = Ensemble_id]
  Exon_region <-
    Exon_region[, .(Chr, ExonStart, ExonEnd, Exon_num, Ensemble_id)]
  Exon_len <-
    Exon_region[, .(ExonLen = ExonEnd - ExonStart + 1), by = .(Exon_num, Ensemble_id)]
  gene_len <- Exon_len[, .(Length = sum(ExonLen)), by = Ensemble_id]
  gene_len <- merge(gene_len, anno[, .(Ensemble_id, Gene)], by = "Ensemble_id", all.x = T)  
  gene_len <- as.data.frame(unique(gene_len))
  gene_len$Ensemble_id <- substr(gene_len$Ensemble_id, 1, 15)
  # write out
  fwrite(
    Exon_region,
    file = paste0(
      "/picb/bigdata/project/FengFYM/s7_annotation_files/All_",
      v,
      "_gene_exon.bed"
    ),
    sep = "\t",
    col.names = T
  )
  fwrite(
    gene_len,
    file = paste0(
      "/picb/bigdata/project/FengFYM/s7_annotation_files/All_",
      v,
      "_gene_len.txt"
    ),
    sep = "\t",
    col.names = T
  )
}


# gene length generated from other R-packages ----------------------------------
# Gene length from GenomicFeatures ---------------------------------------------
#first get the gtf file from, lets say, ensembl
# library(GenomicFeatures)
# # the function above makeTranscriptDbFromGFF is now have different name, 
# # make sure that you type the correct one
# txdb <- makeTxDbFromGFF("/picb/bigdata/project/FengFYM/s7_annotation_files/gencode.v22.annotation.gtf",
#                         format="gtf")
# head(txdb)
# View(txdb)
# # then collect the exons per gene id
# exons.list.per.gene <- exonsBy(txdb,by="gene")
# head(exons.list.per.gene)
# # then for each gene, reduce all the exons to a set of non overlapping exons, 
# # calculate their lengths (widths) and sum then
# exonic.gene.sizes <- lapply(exons.list.per.gene,function(x){sum(width(reduce(x)))})
# head(exonic.gene.sizes)
# #to see them in a table format, you should unlist them
# unlist_geneLength<-unlist(exonic.gene.sizes)
# write.table(unlist_geneLength,
#     "/picb/bigdata/project/FengFYM/s7_annotation_files/hg38_gene_len_GenomicFeatures.txt", 
#     sep = '\t', col.names = T)
#
# # Gene length from IOBR -------------------------------------------------------
# load("/picb/bigdata/project/FengFYM/s7_annotation_files/length_ensembl_IOBR.rda")
# names(length_ensembl) = c("Ensemble_id", "eff_length","symbol")
# tmp = inner_join(length_ensembl, gene_len)  


# gene expression normalization ----------------------------------------
# raw count to FPKM
# CountMat: require raw count table (row: genes(ensemble ID); column: samples)
# v: choose gene length generated from different gtf file
#    'hg38' <- gencode.v22.annotation.gtf
#    'hg19' <- hg19.gencodev27.gtf
#    'mm10' <- gencode.vM27.annotation.gtf (mouse)

Count2FPKM <- function(CountMat,
                       v = c('hg38', 'hg19','mm10')) {
  if (!require(dplyr)) {
    install.packages("dplyr")
  }
  require(dplyr)
  
  if(v %in% c('hg38', 'hg19')){
    gene_length = read.table(
      paste0(
        "/picb/bigdata/project/FengFYM/s7_annotation_files/All_",
        v,
        "_gene_len.txt"
      ),
      sep = '\t',
      header = T
    )
  }else if(v=='mm10'){
    gene_length = read.table(
      paste0(
        "/picb/bigdata/project/FengFYM/s7_annotation_files/All_",
        v,
        "_gene_len_mus.txt"
      ),
      sep = '\t',
      header = T
    )
  }
  CountMat$Gene = rownames(CountMat)
  CountMat_genelength = inner_join(gene_length, CountMat, by = 'Gene')
  genelength_kb = CountMat_genelength$Length / 1000
  rownames(CountMat_genelength) <- CountMat_genelength$Gene
  CountMat = CountMat_genelength[-which(is.element(colnames(CountMat_genelength), 
                                                   names(gene_length)))]
  rpk <- CountMat / genelength_kb
  fpkm <- t(t(rpk) / colSums(CountMat) * 10 ^ 6)
  return(fpkm)
}


# # Count2FPKM_UQ --------------------------------------------------------------
# This function should be reconsidered according to original fomular
#
# Count2FPKM_UQ <- function(CountMat,
#                           v = c('hg38', 'hg19','from_company')) {
#   if (!require(dplyr)) {
#     install.packages("dplyr")
#   }
#   require(dplyr)
#   if(v %in% c('hg38', 'hg19')){
#     gene_length = read.table(
#       paste0(
#         "/picb/bigdata/project/FengFYM/s7_annotation_files/All_",
#         v,
#         "_gene_len.txt"
#       ),
#       sep = '\t',
#       header = T
#     )
#   }else if(v %in% c('mm10','from_company')){
#     gene_length = read.table(
#       paste0(
#         "/picb/bigdata/project/FengFYM/s7_annotation_files/All_",
#         v,
#         "_gene_len_mus.txt"
#       ),
#       sep = '\t',
#       header = T
#     )
#   } 
#   CountMat$Gene = rownames(CountMat)
#   CountMat_genelength = inner_join(gene_length, CountMat, by = 'Gene')
#   genelength_kb = CountMat_genelength$Length / 1000
#   rownames(CountMat_genelength) <- CountMat_genelength$Gene
#   CountMat = CountMat_genelength[-which(is.element(colnames(CountMat_genelength), names(gene_length)))]
#   rpk <- CountMat / genelength_kb
#   fpkm_uq <-
#     t(t(rpk) / apply(CountMat, 2, function(x)
#       summary(x)[5]) * 10 ^ 6)
#   return(fpkm_uq)
# }


# Count2TPM --------------------------------------------------------------------
# CountMat: require raw count table (row: genes(ensemble ID); column: samples)
# v: choose gene length generated from different gtf file
#    'hg38' <- gencode.v22.annotation.gtf
#    'hg19' <- hg19.gencodev27.gtf
#    'mm10' <- gencode.vM27.annotation.gtf (mouse)
#    'from_company' <- Provided by FYS (mouse)

Count2TPM <- function(CountMat,
                      v = c('hg38', 'hg19', 'mm10','from_company')) {
  if (!require(dplyr)) {
    install.packages("dplyr")
  }
  require(dplyr)
  if(v %in% c('hg38', 'hg19')){
    gene_length = read.table(
      paste0(
        "/picb/bigdata/project/FengFYM/s7_annotation_files/All_",
        v,
        "_gene_len.txt"
      ),
      sep = '\t',
      header = T
    )
  }else if(v %in% c('mm10','from_company')){
    gene_length = read.table(
      paste0(
        "/picb/bigdata/project/FengFYM/s7_annotation_files/All_",
        v,
        "_gene_len_mus.txt"
      ),
      sep = '\t',
      header = T
    )
  }
  CountMat$Ensemble_id = rownames(CountMat)
  CountMat_genelength = inner_join(gene_length, CountMat, by = 'Ensemble_id')
  genelength_kb = CountMat_genelength$Length / 1000
  rownames(CountMat_genelength) <- CountMat_genelength$Ensemble_id
  CountMat = CountMat_genelength[-which(is.element(colnames(CountMat_genelength), 
                                                            names(gene_length)))]
  
  rpk <- CountMat / genelength_kb
  tpm <- t(t(rpk) / colSums(rpk) * 10 ^ 6)
  return(tpm)
  
}

# fpkmToTpm --------------------------------------------------------------------
# fpkm: require FPKM table (row: genes; column: samples)

fpkm2TPM <- function(fpkm) {

  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))

}


# TCGA_dataProcessing --------------------------------------------------------------------
# This function is for convert TCGA data for downstreaming analysis, including:
# 1) compute TPM from raw counts table
# 2) convert gene Esemble ID to gene name
# 3) filter out low expression genes (genes not expressed in 20% of samples are removed)
# 3) select samples according to cancer subtype
# 
# cancer: TCGA cancer abbreviation
# sample_type: cancer subtype ('Primary' or "Metastatic")
# rawDataPath: path to TCGA data table downloaded with TCGAbiolinks
# processedDataPath: path to save TCGA processed data for downstreaming analysis

TCGA_dataProcessing <- function(cancer,
                                sample_type = 'Primary',
                                rawDataPath,
                                processedDataPath) {
                                  
  # cancer = 'LIHC'
  cancer_type = paste0('TCGA-', cancer)
  
  # rawDataPath <- "/picb/bigdata/project/FengFYM/s6_TCGAbiolinks_download_normalization/"
  # processedDataPath <- '02_tumorGene_immuneSig_correlation/data/TCGA_processedData/'
  processedDataPath <- paste0(processedDataPath, cancer, '/')
  dir.create(processedDataPath)
  # cibersortDataPath  <- paste0('02_tumorGene_immuneSig_correlation/data/TCGA_cibersortData/', cancer,'/')
  # dir.create(cibersortDataPath)
  # genePatternDataPath  <- paste0('02_tumorGene_immuneSig_correlation/data/TCGA_genePatternData/', cancer,'/')
  # dir.create(genePatternDataPath)
  
  # Prepare TCGA clinical information --------------------------------------------
  
  sampleInfoPath = paste0(
    rawDataPath,
    "TCGA_sampleinfo/",
    cancer,
    ".Rdata"
  )
  
  
  # if (!file.exists(sampleInfoPath)) {
  #   query <- GDCquery(
  #     project = cancer_type,
  #     data.category = "Transcriptome Profiling",
  #     data.type = "Gene Expression Quantification",
  #     workflow.type = "HTSeq - FPKM-UQ"
  #   )
  #   GDCdownload(query, method = "api", files.per.chunk = 200)
  #   data <- GDCprepare(query)
  #   sampleinfo <- colData(data)
  #   save(sampleinfo, file = sampleInfoPath)
  #   
  #   
  # } else{
  #   load(sampleInfoPath)
  #   sampleinfo_cancer <- sampleinfo[c('sample_type')]
  #   sampleinfo_cancer$SampleID <- rownames(sampleinfo_cancer)
  #   print(table(sampleinfo_cancer$sample_type))
  # }
  
  sample_selected = read.csv("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/02_tumorGene_immuneSig_correlation/data/TCGA_patient_ID_tumor_type_selected.csv")
  
  # Prepare TCGA expression matrix --------------------------------------------
  
  # load in full expression matrix
  # drop out low expression genes
  # save to certain path
  
  CountMatPath = paste0(processedDataPath, "mixture_file_TCGA_", cancer,"_",sample_type, "_Count.txt")
  FPKMMatPath = paste0(processedDataPath, "mixture_file_TCGA_", cancer, "_",sample_type,"_FPKM.txt")
  FPKMUQMatPath = paste0(processedDataPath, "mixture_file_TCGA_", cancer,"_",sample_type, "_FPKM_UQ.txt")
  TPMMatPath = paste0(processedDataPath, "mixture_file_TCGA_", cancer,"_",sample_type, "_TPM.txt")
  FPKM2TPMMatPath = paste0(processedDataPath, "mixture_file_TCGA_", cancer,"_",sample_type, "_FPKM2TPM.txt")

  genename = read.csv(
    "02_tumorGene_immuneSig_correlation/data/TCGA_processedData/TCGA_genename.txt",
    sep = '\t',
    header = T
  )
  
  if (!file.exists(CountMatPath)) {
    count <-
      as.data.frame(fread(
        paste0(rawDataPath,
               'Counts/TCGA-',
               cancer,
               '-Counts.csv'
        )
      ))
    rownames(count) <- count$V1
    count = count[-1]
    # TCGA_lihc_genename <- GeneID2Name(rownames(count), species = 'human')
    # write.table(TCGA_lihc_genename, paste0(processedDataPath,"TCGA_",cancer,"_genename.txt"),sep = '\t',row.names = F, quote=F)
    count_filtered = count[rowSums(count > 0) >= 0.2 * ncol(count),]
    
    # count_rowname <- GeneID2Name(rownames(count_filtered), species = 'human')
    count_rowname <- genename
    
    count_filtered$ensembl_gene_id  = rownames(count_filtered)
    count_filtered = inner_join(count_rowname, count_filtered, by = 'ensembl_gene_id')
    count_filtered2 = count_filtered[!duplicated(count_filtered$hgnc_symbol), ]
    rownames(count_filtered2) = count_filtered2$hgnc_symbol
    count_filtered2 = count_filtered2[-which(colnames(count_filtered2) %in% names(genename)), ]

    # count_filtered3 = count_filtered2[sampleinfo_cancer[sampleinfo_cancer$sample_type == sample_type, ]$SampleID]
    
    sample_selected_cancer = sample_selected[sample_selected$TCGA_cancer_type == cancer,]

    count_filtered3 = count_filtered2[is.element(substr(colnames(count_filtered2),1,15),
                                                 sample_selected_cancer[sample_selected_cancer$tumor_type_merged == sample_type,]$TGCA_barcode)]
    
    count_filtered3 %>% write.table(
      paste0(processedDataPath, "mixture_file_TCGA_", cancer,"_",sample_type, "_Count.txt"),
      sep = "\t", 
      quote = F
    )
    # tmp1 <- data.frame(GeneSymbol = rownames(count_filtered3))
    # cbind(tmp1, count_filtered3) %>% write.table(
    #   paste0(
    #     cibersortDataPath,
    #     "mixture_file_TCGA_",
    #     cancer,
    #     "_count_cibersort.txt"
    #   ),
    #   sep = "\t",
    #   quote = F,
    #   row.names = F
    # )
    
  } 

  
  if (!file.exists(FPKMMatPath)) {
    fpkm <-
      as.data.frame(fread(paste0(rawDataPath,
                          'FPKM/TCGA-',cancer,'-FPKM.csv')
      ))
    rownames(fpkm) <- fpkm$V1
    fpkm = fpkm[-1]
    fpkm_filtered = fpkm[rowSums(fpkm > 0) >= 0.2 * ncol(fpkm),]
    
    # fpkm_rowname <- GeneID2Name(rownames(fpkm_filtered), species = 'human')
    fpkm_rowname <- genename
    
    fpkm_filtered$ensembl_gene_id  = rownames(fpkm_filtered)
    fpkm_filtered = inner_join(fpkm_rowname, fpkm_filtered, by = 'ensembl_gene_id')
    fpkm_filtered2 = fpkm_filtered[!duplicated(fpkm_filtered$hgnc_symbol), ]
    rownames(fpkm_filtered2) = fpkm_filtered2$hgnc_symbol
    fpkm_filtered2 = fpkm_filtered2[-c(1, 2)]
    
    # fpkm_filtered3 = fpkm_filtered2[sampleinfo_cancer[sampleinfo_cancer$sample_type == sample_type, ]$SampleID]
    
    sample_selected_cancer = sample_selected[sample_selected$TCGA_cancer_type == cancer,]
    fpkm_filtered3 = fpkm_filtered2[is.element(substr(colnames(fpkm_filtered2),1,15),
                                                 sample_selected_cancer[sample_selected_cancer$tumor_type_merged == sample_type,]$TGCA_barcode)]
    
    
    fpkm_filtered3 %>% write.table(
      paste0(processedDataPath, "mixture_file_TCGA_", cancer,"_",sample_type,  "_FPKM.txt"),
      sep = "\t",
      quote = F
    )
    # tmp3 <- data.frame(GeneSymbol = rownames(fpkm_filtered3))
    # cbind(tmp3, fpkm_filtered3) %>% write.table(
    #   paste0(
    #     cibersortDataPath,
    #     "mixture_file_TCGA_",
    #     cancer,
    #     "_FPKM_cibersort.txt"
    #   ),
    #   sep = "\t",
    #   quote = F,
    #   row.names = F
    # )
    
  } 

  
  
  
  if (!file.exists(FPKMUQMatPath)) {
    fpkmUQ <-
      as.data.frame(fread(paste0(rawDataPath,
                          'FPKM_UQ/TCGA-',cancer,'-FPKM.csv')
      ))
    rownames(fpkmUQ) <- fpkmUQ$V1
    fpkmUQ = fpkmUQ[-1]
    fpkmUQ_filtered = fpkmUQ[rowSums(fpkmUQ > 0) >= 0.2 * ncol(fpkmUQ),]
    
    # fpkmUQ_rowname <- GeneID2Name(rownames(fpkmUQ_filtered), species = 'human')
    fpkmUQ_rowname <- genename
    
    fpkmUQ_filtered$ensembl_gene_id  = rownames(fpkmUQ_filtered)
    fpkmUQ_filtered = inner_join(fpkmUQ_rowname, fpkmUQ_filtered, by = 'ensembl_gene_id')
    fpkmUQ_filtered2 = fpkmUQ_filtered[!duplicated(fpkmUQ_filtered$hgnc_symbol), ]
    rownames(fpkmUQ_filtered2) = fpkmUQ_filtered2$hgnc_symbol
    fpkmUQ_filtered2 = fpkmUQ_filtered2[-c(1, 2)]
    
    # fpkmUQ_filtered3 = fpkmUQ_filtered2[sampleinfo_cancer[sampleinfo_cancer$sample_type == sample_type, ]$SampleID]
    sample_selected_cancer = sample_selected[sample_selected$TCGA_cancer_type == cancer,]
    fpkmUQ_filtered3 = fpkmUQ_filtered2[is.element(substr(colnames(fpkmUQ_filtered2),1,15),
                                               sample_selected_cancer[sample_selected_cancer$tumor_type_merged == sample_type,]$TGCA_barcode)]
    
    fpkmUQ_filtered3 %>% write.table(
      paste0(processedDataPath, "mixture_file_TCGA_", cancer, "_",sample_type, "_FPKM_UQ.txt"),
      sep = "\t",
      quote = F
    )
    # tmp2 <- data.frame(GeneSymbol = rownames(fpkmUQ_filtered3))
    # cbind(tmp2, fpkmUQ_filtered3) %>% write.table(
    #   paste0(
    #     processedDataPath,
    #     "mixture_file_TCGA_",
    #     cancer,
    #     "_FPKMUQ_cibersort.txt"
    #   ),
    #   sep = "\t",
    #   quote = F,
    #   row.names = F
    # )
    
  } 

  
  
  # filtered out normal samples
  
  if (!file.exists(TPMMatPath)) {
    # count_filtered <-
    #   read.csv(
    #     paste0(processedDataPath, "mixture_file_TCGA_", cancer, "_",sample_type, "_Count.txt"),
    #     sep = '\t',
    #     row.names = 1
    #   )
    # colnames(count_filtered) <-
    #   gsub('\\.', '-', colnames(count_filtered))
    
    # TPM_filtered = Count2TPM(count_filtered,
    #                          v = 'hg38')


    count <-
      as.data.frame(fread(
        paste0(rawDataPath,
               'Counts/TCGA-',
               cancer,
               '-Counts.csv'
        ) 
      ))
    rownames(count) <- count$V1
    count = count[-1]

    TPM = Count2TPM(CountMat = count,
                    v = 'hg38')

    write.table(
      TPM,
      paste0("/picb/bigdata/project/FengFYM/s6_TCGAbiolinks_download_normalization/TPM/mixture_file_TCGA_",
             cancer,"_",sample_type, 
             "_TPM_allGeneswithEnsembleID.txt"),
      sep = '\t',
      quote = F
    )

    
    TPM_filtered = as.data.frame(TPM[rowSums(TPM > 0) >= 0.2 * ncol(TPM),])
    TPM_filtered$ensembl_gene_id  = rownames(TPM_filtered)
    
    TPM_filtered_ID2Name = inner_join(genename, TPM_filtered, by = 'ensembl_gene_id')
    TPM_filtered_ID2Name = TPM_filtered_ID2Name[!duplicated(TPM_filtered_ID2Name$hgnc_symbol), ]
    rownames(TPM_filtered_ID2Name) = TPM_filtered_ID2Name$hgnc_symbol
    TPM_filtered_ID2Name = TPM_filtered_ID2Name[-which(colnames(TPM_filtered_ID2Name) %in% names(genename)), ]

    colnames(TPM_filtered_ID2Name) <-
      gsub('\\.', '-', colnames(TPM_filtered_ID2Name))
    sample_selected_cancer = sample_selected[sample_selected$TCGA_cancer_type == cancer,]

    TPM_filtered = TPM_filtered_ID2Name[is.element(substr(colnames(TPM_filtered_ID2Name),1,15),
                                                 sample_selected_cancer[sample_selected_cancer$tumor_type_merged == sample_type,]$TGCA_barcode)]

    write.table(
      TPM_filtered,
      paste0(processedDataPath, "mixture_file_TCGA_",
             cancer,"_",sample_type, 
             "_TPM.txt"),
      sep = '\t',
      quote = F
    )

    # tmp4 <- data.frame(GeneSymbol = rownames(TPM_filtered3))
    # cbind(tmp4, TPM_filtered3) %>% write.table(
    #   paste0(processedDataPath, "mixture_file_TCGA_", cancer, "_TPM_cibersort.txt"),
    #   sep = "\t",
    #   quote = F,
    #   row.names = F
    # )
  } 
  
    if (!file.exists(FPKM2TPMMatPath)) {
      
    fpkm <-
      as.data.frame(fread(paste0(rawDataPath,
                          'FPKM/TCGA-',cancer,'-FPKM.csv')
      ))
    rownames(fpkm) <- fpkm$V1
    fpkm = fpkm[-1]

    TPM2 = as.data.frame(apply(fpkm,2,fpkm2TPM))
    # colSums(TPM)
    TPM_rowname <- genename

    TPM2$ensembl_gene_id  = rownames(TPM2)
    TPM_filtered2 = inner_join(TPM_rowname, TPM2, by = 'ensembl_gene_id')
    TPM_filtered2 = TPM_filtered2[!duplicated(TPM_filtered2$hgnc_symbol), ]
    rownames(TPM_filtered2) = TPM_filtered2$hgnc_symbol
    TPM_filtered2 = TPM_filtered2[-which(colnames(TPM_filtered2) %in% names(genename)), ]

    TPM_filtered3 <- TPM_filtered2[intersect(rownames(TPM_filtered2),rownames(TPM_filtered)),]
    TPM_filtered3 = TPM_filtered3[rowSums(TPM_filtered3 > 0) >= 0.2 * ncol(TPM_filtered3),]

    # count_filtered3 = count_filtered2[sampleinfo_cancer[sampleinfo_cancer$sample_type == sample_type, ]$SampleID]
    
    sample_selected_cancer = sample_selected[sample_selected$TCGA_cancer_type == cancer,]

    TPM_filtered3 = TPM_filtered3[is.element(substr(colnames(TPM_filtered3),1,15),
                                                 sample_selected_cancer[sample_selected_cancer$tumor_type_merged == sample_type,]$TGCA_barcode)]
    
    TPM_filtered3 = TPM_filtered3[rownames(TPM_filtered), colnames(TPM_filtered)]
    
    write.table(
      TPM_filtered3,
      paste0(processedDataPath, "mixture_file_TCGA_",
             cancer,"_",sample_type, 
             "_FPKM2TPM.txt"),
      sep = '\t',
      quote = F
    )


    # tmp4 <- data.frame(GeneSymbol = rownames(TPM_filtered3))
    # cbind(tmp4, TPM_filtered3) %>% write.table(
    #   paste0(processedDataPath, "mixture_file_TCGA_", cancer, "_TPM_cibersort.txt"),
    #   sep = "\t",
    #   quote = F,
    #   row.names = F
    # )
  } 
}
  # prepare file for Gene Pattern
  # dim(fpkmUQ_filtered)
  # dim(fpkm_filtered)
  # dim(TPM_filtered)
  # fpkmUQ_filtered_log = log2(1+fpkmUQ_filtered)
  # fpkm_filtered_log = log2(1+fpkm_filtered)
  # TPM_filtered_log = log2(1+TPM_filtered)
  # annot_fpkmUQ_filtered_log = data.frame(ProbeID = rownames(fpkmUQ_filtered_log),
  #                                         Symbol = rownames(fpkmUQ_filtered_log))
  # annot_fpkm_filtered_log = data.frame(ProbeID = rownames(fpkm_filtered_log),
  #                                       Symbol = rownames(fpkm_filtered_log))
  # annot_TPM_filtered_log = data.frame(ProbeID = rownames(TPM_filtered_log),
  #                                      Symbol = rownames(TPM_filtered_log))
  # fpkmUQ_filtered_log = cbind(annot_fpkmUQ_filtered_log, fpkmUQ_filtered_log)
  # fpkm_filtered_log = cbind(annot_fpkm_filtered_log, fpkm_filtered_log)
  # TPM_filtered_log = cbind(annot_TPM_filtered_log, TPM_filtered_log)
  
  # write.table(fpkmUQ_filtered_log,
  #             paste0(genePatternDataPath,"mixture_file_LIHC_log2fpkmUQ.gct"),
  #             sep = "\t",
  #             quote = F,
  #             row.names = F
  # )
  # write.table(fpkm_filtered_log,
  #             paste0(genePatternDataPath,"mixture_file_LIHC_log2fpkm.gct"),
  #             sep = "\t",
  #             quote = F,
  #             row.names = F
  # )
  # write.table(TPM_filtered_log,
  #             paste0(genePatternDataPath,"mixture_file_LIHC_log2TPM.gct"),
  #             sep = "\t",
  #             quote = F,
  #             row.names = F
  # )
  






# code downloaded from IOBR r-package ------------------------------------------
count2tpm_iobr <- function(countMat, 
                           idType = "Ensembl", 
                           org="hsa",  
                           source = "web", 
                           effLength = NULL, 
                           id = "id", 
                           gene_symbol = "symbol", 
                           length = "eff_length", 
                           remove_redundancy = "mean") {
  # requireNamespace("biomaRt")
  if(class(countMat)!="matrix")  countMat<-as.matrix(countMat)

  if(is.null(effLength) & source == "web"){
    datasets = paste0(c("hsapiens", "mmusculus", "btaurus", "cfamiliaris",
                        "ptroglodytes", "rnorvegicus", "sscrofa"), "_gene_ensembl")
    type = c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol", "start_position", "end_position")
    if(org =="mmu") type[3] = "mgi_symbol"
    # listEnsemblArchives()
    # listMarts()
    # listAttributes()
    ds <- datasets[grepl(org, datasets)]
    mart <- biomaRt::useMart(host = "www.ensembl.org", biomart = 'ENSEMBL_MART_ENSEMBL', dataset = ds)
    ensembl <- biomaRt::getBM(attributes=type, mart = mart)
    ensembl$Length <- abs(ensembl$end_position - ensembl$start_position)

    if(toupper(idType) == "ENSEMBL"){
      len <- ensembl[match(rownames(countMat),ensembl$ensembl_gene_id), "Length"]
      rownames(countMat) = ensembl[match(rownames(countMat),ensembl$ensembl_gene_id), 3]
    }
    else if(toupper(idType) == "SYMBOL")
      len <- ensembl[match(rownames(countMat), ensembl[,3]), "Length"]
    else if(toupper(idType) == "ENTREZ")
      len <- ensembl[match(rownames(countMat), ensembl[,2]), "Length"]
    else
      stop("Please input right type of gene name, such as Ensembl or gene Symbol ...")
  }


  if(source == "default" & tolower(idType) == "ensembl") {

    countMat<-countMat[rownames(countMat)%in%length_ensembl$id,]
    length_ensembl<-length_ensembl[length_ensembl$id%in%rownames(countMat),]

    len<- length_ensembl[match(rownames(countMat), length_ensembl$id), "eff_length"]

    rownames(countMat)<- length_ensembl[match(rownames(countMat),length_ensembl$id), 3]
    countMat<-matrix(as.numeric(countMat), dim(countMat), dimnames = dimnames(countMat))

  }

  if(!is.null(effLength)){

    effLength<-as.data.frame(effLength)
    colnames(effLength)[which(colnames(effLength)==id)]<-"id"
    colnames(effLength)[which(colnames(effLength)==length)]<-"eff_length"
    effLength<-effLength[!duplicated(effLength$id),]


    countMat<-as.matrix(countMat)
    countMat<-countMat[rownames(countMat)%in%effLength$id,]
    effLength<-effLength[effLength$id%in%rownames(countMat),]

    if(id!=gene_symbol){
      # countMat<-as.matrix(countMat)
      colnames(effLength)[which(colnames(effLength)==gene_symbol)]<-"gene_symbol"
      rownames(countMat)<- effLength[match(rownames(countMat),effLength$id), "gene_symbol"]

    }else{
      # countMat<-as.matrix(countMat)
      effLength$gene_symbol<-effLength$id
      # colnames(effLength)[which(colnames(effLength)==gene_symbol)]<-"gene_symbol"
      rownames(countMat)<- effLength[match(rownames(countMat),effLength$id), "gene_symbol"]
    }

    len<- effLength[match(rownames(countMat), effLength[,"gene_symbol"]), "eff_length"]

  }

  na_idx <- which(is.na(len))
  if(length(na_idx)>0){
    warning(paste0("Omit ", length(na_idx), " genes of which length is not available !"))
    countMat <- countMat[!is.na(len),]
    len = len[!is.na(len)]
  }
  tmp <- countMat / len
  TPM <- 1e6 * t(t(tmp) / colSums(tmp))

  if(tolower(remove_redundancy)=="mean"){
    order_index <- apply(TPM,1,function(x) mean(x,na.rm=T))
  }else if(tolower(remove_redundancy)=="sd"){
    order_index <- apply(TPM,1,function(x) sd(x,na.rm=T))
  }else if(tolower(remove_redundancy)=="median"){
    order_index <- apply(TPM,1,function(x) median(x,na.rm=T))
  }

  TPM <-TPM[order(order_index,decreasing=T),]
  TPM <- TPM[!duplicated(rownames(TPM)),]
  TPM <- TPM[!is.na(rownames(TPM)),]
  TPM <- TPM[!rownames(TPM)==" ",]
  TPM <- TPM[,!is.na(colnames(TPM))]
  TPM <- TPM[,!colnames(TPM)==" "]
  return(TPM)
}


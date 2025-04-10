
ITSexpr_purity_detection <- function(genelist_p,
                                   genelist_n,
                                   r_threshold,
                                   cancer = 'LIHC',
                                   sample_type = 'Primary',
                                   tumor_purity_method = "TUMERIC") {
  
  
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
        
        r = cor.test(df$gene, df$purity, method = "pearson")
        
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
        
        r = cor.test(df$gene, df$purity, method = "pearson")
        
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
  genelist_filter <- list(genelist_p, genelist_n)
  return(genelist_filter)
}




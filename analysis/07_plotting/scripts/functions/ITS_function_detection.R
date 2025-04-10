jaccard_index <- function(x, y){
  length(intersect(x,y))/length(union(x,y))
}

Otsuka_Ochiai_coefficient <- function(x, y){
  length(intersect(x, y))/sqrt(length(x) *length(y))
}


ITS_funcSim <- function(ITSenrichRes, 
                        p_threshold = 0.05, 
                        savepath){
  
  ITSenrichRes_sign <- lapply(ITSenrichRes, function(x){
    if(is.null(nrow(x))){
      NULL
    }else if(nrow(x) == 0){
      NULL
    }else {
      x[x$p.adjust < p_threshold, ]$ID
    }
  })

  ITSenrichRes_sign_jcindex <- list()
    for(i in 1:length(ITSenrichRes_sign)){
      x = ITSenrichRes_sign[[i]]
      ITSenrichRes_sign_jcindex[[i]] = sapply(ITSenrichRes_sign, function(y){
        # jaccard_index(x,y)
        gs_similarity(x, y, similarity_method = "jaccard_index")
        # length(intersect(x, y))
      })
    }
    names(ITSenrichRes_sign_jcindex) = names(ITSenrichRes_sign)
    ITSenrichRes_sign_jcindex = do.call(rbind, ITSenrichRes_sign_jcindex)


  ITSenrichRes_sign_Ochiai <- list()
    for(i in 1:length(ITSenrichRes_sign)){
      x = ITSenrichRes_sign[[i]]
      ITSenrichRes_sign_Ochiai[[i]] = sapply(ITSenrichRes_sign, function(y){
        # jaccard_index(x,y)
        gs_similarity(x, y, similarity_method = "Otsuka_Ochiai")
        # length(intersect(x, y))
      })
    }
    names(ITSenrichRes_sign_Ochiai) = names(ITSenrichRes_sign)
    ITSenrichRes_sign_Ochiai = do.call(rbind, ITSenrichRes_sign_Ochiai)


  save(ITSenrichRes_sign_jcindex, ITSenrichRes_sign_Ochiai,
       file= paste0(savepath, "_sim.Rdata"))


}


ITS_function_detection <- function(genelist,
                                 sign = c("p", "n"),
                                 cancer_type,
                                 savepath) {

  dir.create(paste0(savepath, cancer_type,"/"))
  library(org.Hs.eg.db)
  library(topGO)
  library(clusterProfiler)
  library(msigdbr)
  m_df <- msigdbr(species = "Homo sapiens")

  # print(paste0("number of genes = ", length(ITS)))

  genelist_gene2ID <- lapply(genelist, function(ITS){
    gene.df <- bitr(ITS, 
                    fromType = "SYMBOL", 
                    toType = c("ENTREZID", "ENSEMBL", "SYMBOL"), 
                    OrgDb = org.Hs.eg.db)
    return(gene.df)
  })
  names(genelist_gene2ID) <- names(genelist)
  genelist <- genelist_gene2ID

  # cancer hallmark ------------------------------------------------------------
  print("cancer hallmark analysis ############################")
  ITSenrichRes=list()
  for(i in seq(length(genelist))){

    ITSenrichRes[[i]] <- ITSenrichment(ITS = genelist[[i]],
                                       method ="H")
  }
  names(ITSenrichRes) = names(genelist)
  save(ITSenrichRes, 
       file= paste0(savepath, cancer_type,"/", "hallmark_", sign,".Rdata"))

  ITS_funcSim(ITSenrichRes = ITSenrichRes, 
              savepath = paste0(savepath, cancer_type,"/", "hallmark_", sign))

  # GO-BP ---------------------------------------------------------------------
  print("GOBP analysis ############################")
  ITSenrichRes=list()
  for(i in seq(length(genelist))){
    # i=325
    # print(paste0(names(genelist)[i]," number of genes = ", length(genelist[[i]])))
    ITSenrichRes[[i]] <- ITSenrichment(ITS = genelist[[i]],
                                               method ="GOBP")

  }
  names(ITSenrichRes) = names(genelist)
  
  save(ITSenrichRes, 
       file= paste0(savepath, cancer_type,"/", "GOBP_", sign,".Rdata"))
  
  ITS_funcSim(ITSenrichRes = ITSenrichRes, 
              savepath = paste0(savepath, cancer_type,"/", "GOBP_", sign))

    # KEGG ---------------------------------------------------------------------
  print("KEGG analysis ############################")
  ITSenrichRes=list()
  for(i in seq(length(genelist))){
    # i=1
    # print(paste0(names(genelist)[i]," number of genes = ", length(genelist[[i]])))
    ITSenrichRes[[i]] <- ITSenrichment(ITS = genelist[[i]],
                                       method ="KEGG")

  }
  names(ITSenrichRes) = names(genelist)
 
  save(ITSenrichRes, 
       file= paste0(savepath, cancer_type,"/", "KEGG_", sign,".Rdata"))
  
  ITS_funcSim(ITSenrichRes = ITSenrichRes, 
              savepath = paste0(savepath, cancer_type,"/", "KEGG_", sign))

  # Oncogenic signatures ------------------------------------------------------------
  print("Oncogenic signatures analysis ############################")
  
  ITSenrichRes=list()
  for(i in seq(length(genelist))){
    ITSenrichRes[[i]] <-ITSenrichment(ITS = genelist[[i]],
                                        method ="C6")
  }
  names(ITSenrichRes) = names(genelist)

  
  save(ITSenrichRes, 
       file= paste0(savepath, "C6_OncogenicSig_", sign,".Rdata"))
  ITS_funcSim(ITSenrichRes = ITSenrichRes, 
              savepath = paste0(savepath, cancer_type,"/", "C6_OncogenicSig_", sign))

  # # Immunic signatures ------------------------------------------------------------
  # print("Oncogenic signatures analysis ############################")
  
  # ITSenrichRes=list()
  # for(i in seq(length(genelist))){
  #   ITSenrichRes[[i]] <-ITSenrichment(ITS = genelist[[i]],
  #                                       method ="C7")
  # }
  # names(ITSenrichRes) = names(genelist)

  
  # save(ITSenrichRes, 
  #      file= paste0(savepath, cancer_type,"/", "C7_ImmunicSig_", sign,".Rdata"))
  # ITS_funcSim(ITSenrichRes = ITSenrichRes, 
  #             savepath = paste0(savepath, cancer_type,"/", "C7_ImmunicSig_", sign))

}



ITSenrichment <- function(ITS,
                          method = c("KEGG", "GOBP")){
  # functional analysis
  # library(org.Hs.eg.db)
  # library(topGO)
  # library(clusterProfiler)
  # library(msigdbr)
  # m_df <- msigdbr(species = "Homo sapiens")

  # print(paste0("number of genes = ", length(ITS)))
  # gene.df <- bitr(ITS, 
  #                 fromType = "SYMBOL", 
  #                 toType = c("ENTREZID", "ENSEMBL", "SYMBOL"), 
  #                 OrgDb = org.Hs.eg.db)
  gene.df <- ITS
  if(method == "KEGG"){
    res <- enrichKEGG(gene = gene.df$ENTREZID,
                               organism = 'hsa',
                               pvalueCutoff = 1)
    
  }else if(method == "GOBP"){
    res <- enrichGO(gene = gene.df$ENTREZID,
                             'org.Hs.eg.db',
                             ont = "BP",
                             pvalueCutoff = 1)

  }else if(method == "C6"){

    m_t2g <- msigdbr(species = "Homo sapiens", category = "C6") %>% 
        dplyr::select(gs_name, entrez_gene)

    res <- enricher(gene.df$ENTREZID, TERM2GENE=m_t2g)

  }else if(method == "C7"){
    m_t2g <- msigdbr(species = "Homo sapiens", category = "C7") %>% 
        dplyr::select(gs_name, entrez_gene)

    res <- enricher(gene.df$ENTREZID, TERM2GENE=m_t2g)

  }else if(method == "H"){
    m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
        dplyr::select(gs_name, entrez_gene)

    res <- enricher(gene.df$ENTREZID, TERM2GENE=m_t2g)
  }

  # ITSenrichRes <- summary(ITSenrichRes)
  res2 <- as.data.frame(res)
  return(res2)
}

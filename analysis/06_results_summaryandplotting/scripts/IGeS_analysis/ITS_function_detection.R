
 
ITS_function_dection <- function(genelist_p,
                                 genelist_n,
                                 cancer, 
                                 sample_type){
  # 02_tumorGene_immuneSig_correlation/scripts/correlation_dection_plot.R

  # cancer hallmark ------------------------------------------------------------
  print("cancer hallmark analysis ############################")
  
  ITSenrichRes_p=list()
  for(i in seq(length(genelist_p))){
    ITSenrichRes_p[[i]] <-ITSenrichment(ITS = genelist_p[[i]],
                                        method ="H")
  }
  names(ITSenrichRes_p) = names(genelist_p)

  
  ITSenrichRes_n=list()
  for(i in seq(length(genelist_n))){
    ITSenrichRes_n[[i]] <-ITSenrichment(ITS = genelist_n[[i]],
                                        method ="H")
  }
  names(ITSenrichRes_n) = names(genelist_n)
  save(ITSenrichRes_p, ITSenrichRes_n, 
       file= paste0("06_results_summaryandplotting/ITS_function/kegg_", 
                    cancer, "_", sample_type,"_n.Rdata"))

  # Oncogenic signatures ------------------------------------------------------------
  print("Oncogenic signatures analysis ############################")
  
  ITSenrichRes_p=list()
  for(i in seq(length(genelist_p))){
    ITSenrichRes_p[[i]] <-ITSenrichment(ITS = genelist_p[[i]],
                                        method ="C6")
  }
  names(ITSenrichRes_p) = names(genelist_p)

  
  ITSenrichRes_n=list()
  for(i in seq(length(genelist_n))){
    ITSenrichRes_n[[i]] <-ITSenrichment(ITS = genelist_n[[i]],
                                        method ="C6")
  }
  names(ITSenrichRes_n) = names(genelist_n)
  save(ITSenrichRes_p, ITSenrichRes_n, 
       file= paste0("06_results_summaryandplotting/ITS_function/C6_OncogenicSig_", 
                    cancer, "_", sample_type,"_n.Rdata"))

  # Immunic signatures ------------------------------------------------------------
  print("Oncogenic signatures analysis ############################")
  
  ITSenrichRes_p=list()
  for(i in seq(length(genelist_p))){
    ITSenrichRes_p[[i]] <-ITSenrichment(ITS = genelist_p[[i]],
                                        method ="C7")
  }
  names(ITSenrichRes_p) = names(genelist_p)

  
  ITSenrichRes_n=list()
  for(i in seq(length(genelist_n))){
    ITSenrichRes_n[[i]] <-ITSenrichment(ITS = genelist_n[[i]],
                                        method ="C7")
  }
  names(ITSenrichRes_n) = names(genelist_n)
  save(ITSenrichRes_p, ITSenrichRes_n, 
       file= paste0("06_results_summaryandplotting/ITS_function/C7_ImmunicSig_", 
                    cancer, "_", sample_type,"_n.Rdata"))

  # GO-BP ---------------------------------------------------------------------
  ITSenrichRes_p=list()
  for(i in seq(length(genelist_p))){
    # i=325
    # print(paste0(names(genelist_p)[i]," number of genes = ", length(genelist_p[[i]])))
    ITSenrichRes_p[[i]] <-ITSenrichment(ITS = genelist_p[[i]],
                                               method ="gobp")
    # ITSenrichRes_p[[i]] <-ITSenrichment(ITS = genelist_p[[i]],
    #                                     method ="kegg")
    # print(head(ITSenrichRes_p[[i]][ITSenrichRes_p[[i]]$p.adjust < 0.01,][c(2)], 20))
  }
  names(ITSenrichRes_p) = names(genelist_p)
  
  ITSenrichRes_n=list()
  for(i in seq(length(genelist_n))){
    # i=325
    # print(paste0(names(genelist_p)[i]," number of genes = ", length(genelist_p[[i]])))
    ITSenrichRes_n[[i]] <-ITSenrichment(ITS = genelist_n[[i]],
                                               method ="gobp")
    # write.csv(ITSenrichRes_n, paste0("06_results_summaryandplotting/ITS_function/", 
    # print(head(ITSenrichRes[ITSenrichRes$p.adjust < 0.01,][c(1,2)], 10))
  }
  names(ITSenrichRes_n) = names(genelist_n)
  save(ITSenrichRes_p, ITSenrichRes_n, 
       file= paste0("06_results_summaryandplotting/ITS_function/gobp_", 
                    cancer, "_", sample_type,".Rdata"))
  

    # KEGG ---------------------------------------------------------------------
  ITSenrichRes_p=list()
  for(i in seq(length(genelist_p))){
    # i=325
    # print(paste0(names(genelist_p)[i]," number of genes = ", length(genelist_p[[i]])))
    ITSenrichRes_p[[i]] <-ITSenrichment(ITS = genelist_p[[i]],
                                               method ="KEGG")
    # ITSenrichRes_p[[i]] <-ITSenrichment(ITS = genelist_p[[i]],
    #                                     method ="kegg")
    # print(head(ITSenrichRes_p[[i]][ITSenrichRes_p[[i]]$p.adjust < 0.01,][c(2)], 20))
  }
  names(ITSenrichRes_p) = names(genelist_p)
  
  ITSenrichRes_n=list()
  for(i in seq(length(genelist_n))){
    # i=325
    # print(paste0(names(genelist_p)[i]," number of genes = ", length(genelist_p[[i]])))
    ITSenrichRes_n[[i]] <-ITSenrichment(ITS = genelist_n[[i]],
                                               method ="KEGG")
    # write.csv(ITSenrichRes_n, paste0("06_results_summaryandplotting/ITS_function/", 
    # print(head(ITSenrichRes[ITSenrichRes$p.adjust < 0.01,][c(1,2)], 10))
  }
  names(ITSenrichRes_n) = names(genelist_n)
  save(ITSenrichRes_p, ITSenrichRes_n, 
       file= paste0("06_results_summaryandplotting/ITS_function/KEGG_", 
                    cancer, "_", sample_type,".Rdata"))
  

}



ITSenrichment <- function(ITS,
                          method = c("kegg", "gobp")){
  # functional analysis
  library(org.Hs.eg.db)
  library(topGO)
  library(clusterProfiler)
  library(msigdbr)
  m_df <- msigdbr(species = "Homo sapiens")

  print(paste0("number of genes = ", length(ITS)))
  gene.df <- bitr(ITS, 
                  fromType = "SYMBOL", toType = c("ENTREZID", "ENSEMBL", "SYMBOL"), OrgDb = org.Hs.eg.db)

  if(method == "kegg"){
    ITSenrichRes <- enrichKEGG(gene = gene.df$ENTREZID,
                               organism = 'hsa',
                               pvalueCutoff = 1)
    
  }else if(method == "gobp"){
    ITSenrichRes <- enrichGO(gene = gene.df$ENTREZID,
                             'org.Hs.eg.db',
                             ont = "BP",
                             pvalueCutoff = 1)

  }else if(method == "C6"){

    m_t2g <- msigdbr(species = "Homo sapiens", category = "C6") %>% 
        dplyr::select(gs_name, entrez_gene)

    ITSenrichRes <- enricher(gene.df$ENTREZID, TERM2GENE=m_t2g)

  }else if(method == "C7"){
    m_t2g <- msigdbr(species = "Homo sapiens", category = "C7") %>% 
        dplyr::select(gs_name, entrez_gene)

    ITSenrichRes <- enricher(gene.df$ENTREZID, TERM2GENE=m_t2g)
  }else if(method == "H"){
    m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
        dplyr::select(gs_name, entrez_gene)

    ITSenrichRes <- enricher(gene.df$ENTREZID, TERM2GENE=m_t2g)
  }

  ITSenrichRes <- summary(ITSenrichRes)
  return(ITSenrichRes)
}

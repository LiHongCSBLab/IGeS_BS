ITS_GSenrichTest <- function(x,
                             genelist,
                             filename){
  
  geneExprInCCLE <- read.csv("03_drug_immuneSig_enrichment/data/CCLE_geneexpr/geneExprInCCLE.csv")
  geneExprInCCLE_filted <- unique(geneExprInCCLE[geneExprInCCLE$rate < 0.2,]$gene_name)
  genelist = intersect(genelist, geneExprInCCLE_filted)
  
  a = length(intersect(x, genelist))
  b = length(setdiff(genelist, x))
  c = length(setdiff(x, genelist))
  d = length(setdiff(geneExprInCCLE_filted, union(x, genelist)))
  
  N = length(geneExprInCCLE_filted)
  M = length(intersect(genelist, geneExprInCCLE_filted))
  dat = matrix(c(a,c,b,d), ncol=2,nrow=2)
  
  FishTestP = fisher.test(dat)$p.value
  ChisqTestP = chisq.test(dat)$p.value
  # pvalue=1-phyper(a,# 差异基因中，位于通路中基因数量
  #                 a+c, # 差异基因的数量
  #                 N-(a+c), # 全部基因的数量 - 差异数量
  #                 M)  # 全部基因中，位于通路中基因数量
  
  
  # N是数据总库的基因数，
  # M是目标基因集中的基因数，
  # n是获得的差异基因数量，
  # k是在目标基因集中富集的差异基因数
  k = a
  M = length(intersect(genelist, geneExprInCCLE_filted)) # a+b
  n = length(x) # a+c
  N = length(geneExprInCCLE_filted)
  pvalue=phyper(k-1,M, N-M, n, lower.tail=FALSE)                 
  # tmp <- data.frame(gene.not.interest = c(M-k, N-M-n+k), 
  #                   gene.in.interest  = c(k, n-k))
  # row.names(tmp) <- c("In_category", "not_in_category")
  data.frame(filename=filename,
             FishTestP = FishTestP, 
             ChisqTestP = ChisqTestP,
             HypergeometricP = pvalue)
  
  
}

ITS_GSenrichTest_res <- function(x,
                                 savepath){
  res = list()                                
  # ------------------------------------------------------------------------------
  # ITS gene vs drugtarget -------------------------------------------------------
  
  drugtarget <- read.csv('/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/s3_lincs_drug_MoA/drugTargetMerge.txt', sep = '\t')
  
  res[[1]] = ITS_GSenrichTest(x = x, 
                              genelist = unique(drugtarget$target),
                              filename = "ITSp_vs_drugtarget")
  
  
  
  
  # ------------------------------------------------------------------------------
  # ITS gene vs druggable genes --------------------------------------------------
  # ------------------------------------------------------------------------------
  
  drugtarget_OASIS <- read.csv('/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/s3_lincs_drug_MoA/drugtarget_OASIS_TIDE.csv')
  
  res[[2]] = ITS_GSenrichTest(x = x, 
                              genelist = drugtarget_OASIS$gene_name, 
                              filename = "ITSp_vs_druggable_gene")
  
  
  # ------------------------------------------------------------------------------
  # ITS gene vs inhibitory and stimulatory genes ---------------------------------
  # ------------------------------------------------------------------------------
  
  # ------------------------------------------------------------------------------
  # inhibitory
  immuneModulator_PMID29628290_1 <- as.data.frame(read_excel("06_results_summaryandplotting/data/immune_modulator/PMID29628290.xlsx", sheet = 1))
  inhibitory <- unique(immuneModulator_PMID29628290_1[immuneModulator_PMID29628290_1$`Immune Checkpoint` == "Inhibitory",]$`HGNC Symbol`)
  inhibitory <- inhibitory[!is.na(inhibitory)]
  
  res[[3]] = ITS_GSenrichTest(x = x, 
                              genelist = inhibitory, 
                              filename = "ITSp_vs_inhibitor_PMID29628290")
  
  
  # ------------------------------------------------------------------------------
  # stimulatory
  immuneModulator_PMID29628290_1 <- as.data.frame(read_excel("06_results_summaryandplotting/data/immune_modulator/PMID29628290.xlsx", sheet = 1))
  stimulatory <- unique(immuneModulator_PMID29628290_1[immuneModulator_PMID29628290_1$`Immune Checkpoint` != "Inhibitory",]$`HGNC Symbol`)
  stimulatory <- stimulatory[!is.na(stimulatory)]
  
  
  res[[4]] = ITS_GSenrichTest(x = x, 
                              genelist = stimulatory, 
                              filename = "ITSp_vs_stimulatory_PMID29628290")
  
  
  
  # ------------------------------------------------------------------------------
  # co_stim_inhib 
  
  immuneModulator_PMID29628290_2 <- as.data.frame(read_excel("06_results_summaryandplotting/data/immune_modulator/PMID29628290.xlsx", sheet = 2))
  co_stim_inhib <- unique(immuneModulator_PMID29628290_2[immuneModulator_PMID29628290_2$Reason == "co-stim/inhib list", ]$Gene)
  
  res[[5]] = ITS_GSenrichTest(x = x, 
                              genelist = co_stim_inhib, 
                              filename = "ITSp_vs_co_stim_inhib_PMID29628290")
  
  
  # ------------------------------------------------------------------------------
  # ITS gene vs icb_enhancer and icb_suppressor genes ----------------------------
  # ------------------------------------------------------------------------------
  
  immuneModulator_PMID34980132_1 <- as.data.frame(read_excel("06_results_summaryandplotting/data/immune_modulator/PMID34980132/ICB enhancer and suppressor genes.xlsx", sheet = 1))
  # ------------------------------------------------------------------------------
  icb_enhancer <- unique(immuneModulator_PMID34980132_1[immuneModulator_PMID34980132_1$Type == "ICB enhancer gene", ]$Symbol)
  
  res[[6]] = ITS_GSenrichTest(x = x, 
                              genelist = icb_enhancer, 
                              filename = "ITSp_vs_icb_enhancer_PMID34980132")
  
  
  # ------------------------------------------------------------------------------
  immuneModulator_PMID34980132_1 <- as.data.frame(read_excel("06_results_summaryandplotting/data/immune_modulator/PMID34980132/ICB enhancer and suppressor genes.xlsx", sheet = 1))
  # ------------------------------------------------------------------------------
  icb_suppressor <- unique(immuneModulator_PMID34980132_1[immuneModulator_PMID34980132_1$Type == "ICB suppressors gene", ]$Symbol)
  
  res[[7]] = ITS_GSenrichTest(x = x, 
                              genelist = icb_suppressor, 
                              filename = "ITSp_vs_icb_suppressor_PMID34980132")
  
  
  # ------------------------------------------------------------------------------
  # ITS gene vs regulator genes --------------------------------------------------
  # ------------------------------------------------------------------------------
  
  immuneModulator_PMID34980132_2 <- as.data.frame(read_excel("06_results_summaryandplotting/data/immune_modulator/PMID34980132/positive and negative regulators.xlsx", sheet = 1))
  # ------------------------------------------------------------------------------
  # regulator_n
  regulator_p <- unique(immuneModulator_PMID34980132_2[immuneModulator_PMID34980132_2$Type == "Positive regulators", ]$Symbol)
  
  res[[8]] = ITS_GSenrichTest(x = x, 
                              genelist = regulator_p, 
                              filename = "ITSp_vs_regulator_p_PMID34980132")
  
  # ------------------------------------------------------------------------------
  # regulator_n
  immuneModulator_PMID34980132_2 <- as.data.frame(read_excel("06_results_summaryandplotting/data/immune_modulator/PMID34980132/positive and negative regulators.xlsx", sheet = 1))
  regulator_n <- unique(immuneModulator_PMID34980132_2[immuneModulator_PMID34980132_2$Type == "Negative regulators", ]$Symbol)
  
  
  res[[9]] = ITS_GSenrichTest(x = x, 
                              genelist = regulator_n, 
                              filename = "ITSp_vs_regulator_n_PMID34980132")
  
  # ------------------------------------------------------------------------------
  # ITS gene vs resistor and sensitizor  -----------------------------------------
  # ------------------------------------------------------------------------------
  
  immuneModulator_PMID34980132_3 <- as.data.frame(read_excel("06_results_summaryandplotting/data/immune_modulator/PMID34980132/sensitizer_resistor_genes.xlsx", sheet = 1))
  # ------------------------------------------------------------------------------
  # Sensitizer
  Sensitizer <- unique(immuneModulator_PMID34980132_3[immuneModulator_PMID34980132_3$Type == "Sensitizer", ]$Symbol)
  
  res[[10]] = ITS_GSenrichTest(x = x, 
                               genelist = Sensitizer, 
                               filename = "ITSp_vs_Sensitizer_PMID34980132")
  
  # ------------------------------------------------------------------------------
  # Resistor
  Resistor <- unique(immuneModulator_PMID34980132_3[immuneModulator_PMID34980132_3$Type == "Resistor", ]$Symbol)
  
  res[[11]] = ITS_GSenrichTest(x = x, 
                               genelist = Resistor, 
                               filename = "ITSp_vs_Resistor_PMID34980132")
  
  # ------------------------------------------------------------------------------
  # ITS gene vs OG and TSG -------------------------------------------------------
  # ------------------------------------------------------------------------------
  
  immuneModulator_PMID34980132_4 <- as.data.frame(read_excel("06_results_summaryandplotting/data/immune_modulator/PMID34980132/Tumor immunity related OGs TSGs.xlsx", sheet = 1))
  OG <- unique(immuneModulator_PMID34980132_4[immuneModulator_PMID34980132_4$Type == "Tumor immunity-related OGs", ]$Symbol)
  
  
  res[[12]] = ITS_GSenrichTest(x = x, 
                               genelist = OG, 
                               filename = "ITSp_vs_OG_PMID34980132")
  
  TSG <- unique(immuneModulator_PMID34980132_4[immuneModulator_PMID34980132_4$Type == "Tumor immunity-related TSGs", ]$Symbol)
  
  res[[13]] = ITS_GSenrichTest(x = x, 
                               genelist = TSG, 
                               filename = "ITSp_vs_TSG_PMID34980132")
  
  # ------------------------------------------------------------------------------
  # cancer oncogenes -------------------------------------------------------------
  # ------------------------------------------------------------------------------
  canceroncogene <- read.csv("06_results_summaryandplotting/data/canceroncegenes/oncoKB_cancerGeneList.tsv", sep = "\t", header = T)
  oncogene <- canceroncogene[canceroncogene$Is.Oncogene == "Yes",]$Hugo.Symbol
  res[[14]] = ITS_GSenrichTest(x = x, 
                               genelist = oncogene, 
                               filename = "ITSp_vs_OG_oncoKB")
  
  
  suppressorgene <- canceroncogene[canceroncogene$Is.Tumor.Suppressor.Gene == "Yes",]$Hugo.Symbol
  res[[15]] = ITS_GSenrichTest(x = x, 
                               genelist = suppressorgene, 
                               filename = "ITSp_vs_TSG_oncoKB")
  
  res <- do.call(rbind, res)
  # write.csv(res, paste0(savepath))
  return(res)
}




ITS_GS_preparation <- function(){
  gs_all = list()                                
  # ------------------------------------------------------------------------------
  # ITS gene vs drugtarget -------------------------------------------------------
  
  drugtarget <- read.csv('/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/s3_lincs_drug_MoA/drugTargetMerge.txt', sep = '\t')
  gs_all[[1]] <- data.frame(gs_name = "ITSp_vs_drugtarget",
                            SYMBOL = unique(drugtarget$target))
  
  # ------------------------------------------------------------------------------
  # ITS gene vs druggable genes --------------------------------------------------
  # ------------------------------------------------------------------------------
  
  drugtarget_OASIS <- read.csv('/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/s3_lincs_drug_MoA/drugtarget_OASIS_TIDE.csv')
  
  gs_all[[2]] = data.frame(gs_name = "ITSp_vs_druggable_gene",
                           SYMBOL = unique(drugtarget_OASIS$gene_name))
  
  
  # ------------------------------------------------------------------------------
  # ITS gene vs inhibitory and stimulatory genes ---------------------------------
  # ------------------------------------------------------------------------------
  
  # ------------------------------------------------------------------------------
  # inhibitory
  immuneModulator_PMID29628290_1 <- as.data.frame(read_excel("06_results_summaryandplotting/data/immune_modulator/PMID29628290.xlsx", sheet = 1))
  inhibitory <- unique(immuneModulator_PMID29628290_1[immuneModulator_PMID29628290_1$`Immune Checkpoint` == "Inhibitory",]$`HGNC Symbol`)
  inhibitory <- inhibitory[!is.na(inhibitory)]
  
  gs_all[[3]] = data.frame(gs_name = "ITSp_vs_inhibitor_PMID29628290",
                           SYMBOL = inhibitory)
  
  # ------------------------------------------------------------------------------
  # stimulatory
  immuneModulator_PMID29628290_1 <- as.data.frame(read_excel("06_results_summaryandplotting/data/immune_modulator/PMID29628290.xlsx", sheet = 1))
  stimulatory <- unique(immuneModulator_PMID29628290_1[immuneModulator_PMID29628290_1$`Immune Checkpoint` != "Inhibitory",]$`HGNC Symbol`)
  stimulatory <- stimulatory[!is.na(stimulatory)]
  
  gs_all[[4]] = data.frame(gs_name = "ITSp_vs_stimulatory_PMID29628290",
                           SYMBOL = stimulatory)
  # ------------------------------------------------------------------------------
  # co_stim_inhib 
  
  immuneModulator_PMID29628290_2 <- as.data.frame(read_excel("06_results_summaryandplotting/data/immune_modulator/PMID29628290.xlsx", sheet = 2))
  co_stim_inhib <- unique(immuneModulator_PMID29628290_2[immuneModulator_PMID29628290_2$Reason == "co-stim/inhib list", ]$Gene)
  
  gs_all[[5]] = data.frame(gs_name = "ITSp_vs_co_stim_inhib_PMID29628290",
                           SYMBOL = co_stim_inhib)
  
  
  # ------------------------------------------------------------------------------
  # ITS gene vs icb_enhancer and icb_suppressor genes ----------------------------
  # ------------------------------------------------------------------------------
  
  immuneModulator_PMID34980132_1 <- as.data.frame(read_excel("06_results_summaryandplotting/data/immune_modulator/PMID34980132/ICB enhancer and suppressor genes.xlsx", sheet = 1))
  # ------------------------------------------------------------------------------
  icb_enhancer <- unique(immuneModulator_PMID34980132_1[immuneModulator_PMID34980132_1$Type == "ICB enhancer gene", ]$Symbol)
  
  gs_all[[6]] = data.frame(gs_name = "ITSp_vs_icb_enhancer_PMID34980132",
                           SYMBOL = icb_enhancer)
  
  # ------------------------------------------------------------------------------
  immuneModulator_PMID34980132_1 <- as.data.frame(read_excel("06_results_summaryandplotting/data/immune_modulator/PMID34980132/ICB enhancer and suppressor genes.xlsx", sheet = 1))
  # ------------------------------------------------------------------------------
  icb_suppressor <- unique(immuneModulator_PMID34980132_1[immuneModulator_PMID34980132_1$Type == "ICB suppressors gene", ]$Symbol)
  
  gs_all[[7]] = data.frame(gs_name = "ITSp_vs_icb_suppressor_PMID34980132",
                           SYMBOL = icb_suppressor)
  
  
  # ------------------------------------------------------------------------------
  # ITS gene vs regulator genes --------------------------------------------------
  # ------------------------------------------------------------------------------
  
  immuneModulator_PMID34980132_2 <- as.data.frame(read_excel("06_results_summaryandplotting/data/immune_modulator/PMID34980132/positive and negative regulators.xlsx", sheet = 1))
  # ------------------------------------------------------------------------------
  # regulator_n
  regulator_p <- unique(immuneModulator_PMID34980132_2[immuneModulator_PMID34980132_2$Type == "Positive regulators", ]$Symbol)
  
  
  gs_all[[8]] = data.frame(gs_name = "ITSp_vs_regulator_p_PMID34980132",
                           SYMBOL = regulator_p)
  
  # ------------------------------------------------------------------------------
  # regulator_n
  immuneModulator_PMID34980132_2 <- as.data.frame(read_excel("06_results_summaryandplotting/data/immune_modulator/PMID34980132/positive and negative regulators.xlsx", sheet = 1))
  regulator_n <- unique(immuneModulator_PMID34980132_2[immuneModulator_PMID34980132_2$Type == "Negative regulators", ]$Symbol)
  
  
  gs_all[[9]] = data.frame(gs_name = "ITSp_vs_regulator_n_PMID34980132",
                           SYMBOL = regulator_n)
  
  # ------------------------------------------------------------------------------
  # ITS gene vs resistor and sensitizor  -----------------------------------------
  # ------------------------------------------------------------------------------
  
  immuneModulator_PMID34980132_3 <- as.data.frame(read_excel("06_results_summaryandplotting/data/immune_modulator/PMID34980132/sensitizer_resistor_genes.xlsx", sheet = 1))
  # ------------------------------------------------------------------------------
  # Sensitizer
  Sensitizer <- unique(immuneModulator_PMID34980132_3[immuneModulator_PMID34980132_3$Type == "Sensitizer", ]$Symbol)
  
  gs_all[[10]] = data.frame(gs_name = "ITSp_vs_Sensitizer_PMID34980132",
                            SYMBOL = Sensitizer)
  
  # ------------------------------------------------------------------------------
  # Resistor
  Resistor <- unique(immuneModulator_PMID34980132_3[immuneModulator_PMID34980132_3$Type == "Resistor", ]$Symbol)
  
  gs_all[[11]] = data.frame(gs_name = "ITSp_vs_Resistor_PMID34980132",
                            SYMBOL = Resistor)
  
  # ------------------------------------------------------------------------------
  # ITS gene vs OG and TSG -------------------------------------------------------
  # ------------------------------------------------------------------------------
  
  immuneModulator_PMID34980132_4 <- as.data.frame(read_excel("06_results_summaryandplotting/data/immune_modulator/PMID34980132/Tumor immunity related OGs TSGs.xlsx", sheet = 1))
  OG <- unique(immuneModulator_PMID34980132_4[immuneModulator_PMID34980132_4$Type == "Tumor immunity-related OGs", ]$Symbol)
  
  gs_all[[12]] = data.frame(gs_name = "ITSp_vs_OG_PMID34980132",
                            SYMBOL = OG)
  
  TSG <- unique(immuneModulator_PMID34980132_4[immuneModulator_PMID34980132_4$Type == "Tumor immunity-related TSGs", ]$Symbol)
  gs_all[[13]] = data.frame(gs_name = "ITSp_vs_TSG_PMID34980132",
                            SYMBOL = TSG)
  
  # ------------------------------------------------------------------------------
  # cancer oncogenes -------------------------------------------------------------
  # ------------------------------------------------------------------------------
  canceroncogene <- read.csv("06_results_summaryandplotting/data/canceroncegenes/oncoKB_cancerGeneList.tsv", sep = "\t", header = T)
  oncogene <- canceroncogene[canceroncogene$Is.Oncogene == "Yes",]$Hugo.Symbol
  gs_all[[14]] = data.frame(gs_name = "ITSp_vs_OG_oncoKB",
                            SYMBOL = oncogene)
  
  suppressorgene <- canceroncogene[canceroncogene$Is.Tumor.Suppressor.Gene == "Yes",]$Hugo.Symbol
  gs_all[[15]] = data.frame(gs_name = "ITSp_vs_TSG_oncoKB",
                            SYMBOL = suppressorgene)
  
  gs_all <- do.call(rbind, gs_all)
  save(gs_all, file = "07_plotting/02_ITScharacteristic/custom_gs.Rdata")
  
}



ITS_GS_enrich <- function(x,
                          gs=NULL,
                          pvalueCutoff = NULL,
                          type = c("custom", "H", "KEGG", "GOBP"),
                          method = c('ORA', 'GSEA','custom'),
                          savepath = NULL){
  
  
  library(org.Hs.eg.db)
  library(topGO)
  library(clusterProfiler)
  library(msigdbr)
  library(enrichplot)
  m_df <- msigdbr(species = "Homo sapiens")
  print(paste0("number of genes = ", length(x)))
  
  # LOAD GENES IN CCLE -------------------------------------------------------
  geneExprInCCLE <- read.csv("03_drug_immuneSig_enrichment/data/CCLE_geneexpr/geneExprInCCLE.csv")
  
  protein_coding_genes <- read.table("03_drug_immuneSig_enrichment/data/protein_coding_genes_hg38.txt", sep = '\t', header = T)
  geneExprInCCLE_filtered <- inner_join(geneExprInCCLE, protein_coding_genes)
  geneExprInCCLE_filtered <- geneExprInCCLE_filtered[geneExprInCCLE_filtered$rate < 0.2,]
  bg.gene <- bitr(geneExprInCCLE_filtered$gene_name, 
                  fromType = "SYMBOL", 
                  toType = c("ENTREZID", "ENSEMBL", "SYMBOL"), 
                  OrgDb = org.Hs.eg.db)
          
  # load("07_plotting/02_ITScharacteristic/custom_gs.Rdata")
  if(method == "GSEA"){
    
      genename <- bitr(names(x), 
                  fromType = "SYMBOL", 
                  toType = c("ENTREZID",  "SYMBOL"),  #"ENSEMBL",
                  OrgDb = org.Hs.eg.db)
    
    if(is.null(gs)){
      if(type == "KEGG"){
        # names(x) <- genename$SYMBOL
        res <- gseKEGG(geneList     = x,
                       organism     = 'hsa',
                       minGSSize    = 10,
                       pvalueCutoff = pvalueCutoff,
                       verbose      = FALSE)
        
        res2 = data.frame(res)
        n = nrow(res2[res2$p.adjust < 0.1,])
        
        if(n>10){
          pdf(paste0(savepath, ".pdf"), width = 15, height = 12) 
          #ridgeplot(res)       
          gseaplot2(res,geneSetID=1:10,pvalue_table=T)   
          dev.off()
        }else if( n<=10 & n>0 ){
          pdf(paste0(savepath, ".pdf"), width = 15, height = 12) 
          #ridgeplot(res)       
          gseaplot2(res,geneSetID=1:n,pvalue_table=T)   
          dev.off()
        }

      }else if( type == "H"){
        names(x) <- genename$ENTREZID

        m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
          dplyr::select(gs_name, entrez_gene)
        res <- GSEA(x,
                    pvalueCutoff = pvalueCutoff, 
                    TERM2GENE=m_t2g)
        res2 = data.frame(res)
        n = nrow(res2[res2$p.adjust < 0.1,])
        
        if(n>10){
          pdf(paste0(savepath, ".pdf"), width = 15, height = 12) 
          #ridgeplot(res)       
          gseaplot2(res,geneSetID=1:10,pvalue_table=T)   
          dev.off()
        }else if( n<=10 & n>0 ){
          pdf(paste0(savepath, ".pdf"), width = 15, height = 12) 
          #ridgeplot(res)       
          gseaplot2(res,geneSetID=1:n,pvalue_table=T)   
          dev.off()          
        }
        
      }else if( type == "GOBP"){
        names(x) <- genename$ENTREZID

        res <- gseGO(geneList     = x,
                     OrgDb        = org.Hs.eg.db,
                     ont          = "BP",
                     minGSSize    = 0,
                     maxGSSize    = 10000,
                     pvalueCutoff = pvalueCutoff,
                     verbose      = FALSE)
        res2 = data.frame(res)
        n = nrow(res2[res2$p.adjust < 0.1,])
        
        if(n>10){
          pdf(paste0(savepath, ".pdf"), width = 15, height = 12) 
          #ridgeplot(res)       
          gseaplot2(res,geneSetID=1:10,pvalue_table=T)   
          dev.off()
        }else if( n<=10 & n>0 ){
          pdf(paste0(savepath, ".pdf"), width = 15, height = 12) 
          #ridgeplot(res)       x
          gseaplot2(res,geneSetID=1:n,pvalue_table=T)   
          dev.off()          
        }

      }else if(type == "Reactome"){
        names(x) <- genename$ENTREZID

        res <- ReactomePA::gsePathway(x, 
                                      pvalueCutoff = pvalueCutoff,
                                      pAdjustMethod = "BH", 
                                      verbose = FALSE)
        res2 = data.frame(res)
        n = nrow(res2[res2$p.adjust < 0.1,])
                              
        if(n>10){
          pdf(paste0(savepath, ".pdf"), width = 15, height = 12) 
          #ridgeplot(res)       
          gseaplot2(res,geneSetID=1:10,pvalue_table=T)   
          dev.off()
        }else if( n<=10 & n>0 ){
          pdf(paste0(savepath, ".pdf"), width = 15, height = 12) 
          #ridgeplot(res)       
          gseaplot2(res,geneSetID=1:n,pvalue_table=T)   
          dev.off()          
        }
        
      }
      
    }else{
      
      res <- GSEA(geneList = x, 
                  pvalueCutoff = pvalueCutoff, 
                  minGSSize = 0,
                  maxGSSize = 5000,
                  TERM2GENE = gs[c(1,2)])
        # if(n>10){
        #   pdf(paste0(savepath, ".pdf"), width = 15, height = 12) 
        #   #ridgeplot(res)       
        #   gseaplot2(res,geneSetID=1:10,pvalue_table=T)   
        #   dev.off()
        # }else if( n<=10 & n>0 ){
        #   pdf(paste0(savepath, ".pdf"), width = 15, height = 12) 
        #   #ridgeplot(res)       x
        #   gseaplot2(res,geneSetID=1:n,pvalue_table=T)   
        #   dev.off()          
        # }
    }
  }else if(method == "ORA"){
    
    gene.df <- bitr(x, 
                    fromType = "SYMBOL", 
                    toType = c("ENTREZID", "ENSEMBL", "SYMBOL"), 
                    OrgDb = org.Hs.eg.db)
    
    if(is.null(gs)){
      if(type == "KEGG"){
        require(KEGG.db)
        res <- clusterProfiler::enrichKEGG(gene = gene.df$ENTREZID,
                          universe = bg.gene$ENTREZID,
                          organism = 'hsa',
                          keyType = "kegg",
                          # minGSSize    = 0,
                          # maxGSSize    = 10000,
                           use_internal_data = T,
                          pvalueCutoff = pvalueCutoff)

      }else if( type == "H"){
        m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
          dplyr::select(gs_name, gene_symbol, entrez_gene)
        res <- enricher(gene.df$SYMBOL, #$ENTREZID, ,
                        universe = bg.gene$SYMBOL, #ENTREZID,
                        pvalueCutoff =pvalueCutoff, 
                        TERM2GENE=m_t2g)

      }else if( type == "GOBP"){
        res <- enrichGO(gene = gene.df$ENTREZID,
                        universe = bg.gene$ENTREZID,
                        'org.Hs.eg.db',
                        ont = "BP",
                        pvalueCutoff = pvalueCutoff)

      }else if(type == "Reactome"){
        res <- ReactomePA::enrichPathway(gene = gene.df$ENTREZID,
                             universe = bg.gene$ENTREZID,
                             pvalueCutoff = pvalueCutoff,
                             organism='human', 
                             readable=TRUE)
      }
    }else{

      bg.gene2 = bg.gene[c(1,2)]
      bg.gene2$gs_name = "bg"
      bg.gene2 = bg.gene2[c(3,1,2)]

          
      gs2 <- bitr(gs$SYMBOL, 
                      fromType = "SYMBOL", 
                      toType = c("ENTREZID",  "SYMBOL"), 
                      OrgDb = org.Hs.eg.db)
      gs3 <- merge(gs, gs2, by = "SYMBOL")
      gs3 <- gs3[c(2,1,3)]
      gs3 <- gs3[order(gs3$gs_name),]      
      gs3 = rbind(gs3[c(1,3)], bg.gene2[c(1,3)])
      res <- enricher(gene.df$ENTREZID, 
                      universe = bg.gene$ENTREZID,
                      pvalueCutoff = pvalueCutoff, 
                      minGSSize = 0,
                      maxGSSize = 10000,
                      TERM2GENE = gs3)

    }
    
  }else if (method == "custom") {
    res = list()                              
    
    for(i in seq(length(unique(gs$gs_name)))){
      filename = unique(gs$gs_name)[i]
      res[[i]] = ITS_GSenrichTest(x = x, 
                                  genelist = unique(gs[gs$gs_name == filename, ]$SYMBOL),
                                  filename = filename)
    }
    res <- do.call(rbind, res)
    # write.csv(res, paste0(savepath))
    # return(res)
  } 

  if(is.null(savepath)){
    return(res)
  }else{

    res2 = as.data.frame(res)
    write.table(res2, paste0(savepath, ".txt"), quote = F, sep = '\t', row.names = F)
    save(res, file=paste0(savepath, ".Rdata"))
    
    return(res)
  }
  dev.off()
}
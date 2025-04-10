listDrug <- function(cancer = "LIHC",
                     purity_method = "TUMERIC",
                     dataset = "70138",
                     cellline = "HEPG2",
                     filepath){

  load(paste0(filepath, purity_method,"/",cancer,"_",dataset,"_",cellline, "_drugDEG_immunesig.Rdata"))
  drug_index_file = read.csv("r1_drug_immuneSig_correlation_analysis/06_results_summaryandplotting/data/drug_index_LINCS.txt",sep = "\t")
  
  return(drug_index_file[which(is.element(drug_index_file$drug_index, names(drugDEG_immunesig))),])
}



listImmuneSig <- function(cancer = "LIHC",
                          purity_method = "TUMERIC",
                          dataset = "70138",
                          cellline = "HEPG2",
                          immuneKeyword = NULL,
                          filepath ){

  load(paste0(filepath, "/",purity_method,"/",cancer,"_",dataset,"_",cellline, "_drugDEG_immunesig.Rdata"))
  immunesiglist = unique(do.call(rbind, drugDEG_immunesig)$immunesig)
  if(is.null(immuneKeyword)){
    return(immunesiglist)
  }else{
    return(immunesiglist[grep(immuneKeyword, immunesiglist)])
  }
}


# INPUT: 
#   REQUIRE: DRUG, GENE, CANCER TYPE
#   OPTICAL: PURITY_METHOD (DEFAULT: TUMERIC), CORRELATION_METHOD (DEFAULT: PEARSON)
drug_target_DEG_immunesig_match <- function(cancer,
                                            purity_method,
                                            datatype = 'allgenes',
                                            dataset,
                                            drug_keyword, 
                                            r = NULL,
                                            padj = 0.05,
                                            p = NULL,
                                            drug_logFC = NULL,
                                            drug_padj = 0.1,
                                            drug_p = NULL,
                                            drugMergeFCpath,
                                            immusig_path,
                                            resultpath){
  
  
  
  library(reshape2)
  library(stringr)
  library(dplyr)
  
  # cancer = 'LIHC'
  # purity_method = 'TUMERIC'
  # immusig_path <- "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/02_tumorGene_immuneSig_correlation/results/"
  
  if(purity_method == 'IHC'){
    load(paste0(immusig_path, "immunecell_TCGAgenes_pcor_pearson_IHC/pcor_pearson_",cancer,"_Primary_IHC.Rdata"))
    
  }else{
    load(paste0(immusig_path, "immunecell_TCGAgenes_result_pearson_", purity_method,"/pcor_pearson_",cancer,"_Primary_", purity_method,".Rdata"))
    
  }
  
  
  gene_immunesig_sigcor <- lapply(cor_result, function(x){
    if(!is.null(x)){
      # x=cor_result[[1]]
      x = data.frame(r = unlist(x$cor),
                     p.value = unlist(x$p.value),
                     adj.p = unlist(x$p.adj))
      x$gene = row.names(x)
      # x = inner_join(drugTarget['target'],x)
      x <- x[x$adj.p < padj,]
      if(is.null(r)){
        s1 <- summary(x[x$r > 0,]$r)[5]
        s2 <- summary(x[x$r < 0,]$r)[2]
        y <- unique(x[(x$r > s1 | x$r < s2) & x$adj.p < padj,])            
      }else{
        y <- unique(x[(abs(x$r) > r) & x$adj.p < padj,])            
      }
      
      # y <- y[!is.na(y$r),]
      return(y)
    }
  })
  
  
  # datatype = 'allgenes'
  # dataset = '70138'
  drugTargetMoA <- read.table("r1_drug_immuneSig_correlation_analysis/06_results_summaryandplotting/data/drugTargetMoA_merged.txt", sep = '\t', header = T)
  drugindexall <- read.table("r1_drug_immuneSig_correlation_analysis/06_results_summaryandplotting/data/drug_index_LINCS.txt", sep = "\t", header = T)
  drugTargetMoA <- merge(drugTargetMoA, drugindexall, by = c("pert_iname", "drug_index"), all =  T)
  # drugMergeFCpath = "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/s2_drug_inducted_transcriptome_change/mergeFC_r1/"
  drugMergeFCpath = paste0(drugMergeFCpath, "/", datatype, "/", dataset, "/", cancer, "/")
  alldruglist <- gsub("_allgeneFC.csv","", list.files(drugMergeFCpath))

  if(!is.null(drug_keyword)){
    drug_list <- drugTargetMoA[grep( drug_keyword, drugTargetMoA$pert_iname), ]
    
  }else{
    drug_list <- drugTargetMoA
  }
  
  
  if(nrow(drug_list)>0){
    drug_list2 <- intersect(drug_list$drug_index, alldruglist)
    
    drugDEG_immunesig <- lapply(as.list(drug_list2), function(drug){
      allgeneFC <- read.csv(paste0( drugMergeFCpath, drug, "_allgeneFC.csv"))
      names(allgeneFC) <- c('gene','drug_logFC','drug_P.Value','drug_adj.P.Val','drug_index')
      drugGeneFC_immunesig <- lapply(gene_immunesig_sigcor, function(x){
        if(!is.null(x)){
          inner_join(x, allgeneFC)
        }
      })
      
      drugDEG_immunesig1 <- lapply(drugGeneFC_immunesig, function(x){
        if(!is.null(x)){
          # x[abs(x$drug_logFC) > drug_logFC & x$drug_adj.P.Val < drug_padj, ]
          x[x$drug_adj.P.Val < drug_padj, ]
        }
      })
      # lapply(drugDEG_immunesig,nrow)
      drugDEG_immunesig1 <- do.call(rbind, drugDEG_immunesig1)
      drugDEG_immunesig1$immunesig <- gsub("\\.\\d+","",rownames(drugDEG_immunesig1))
      drugDEG_immunesig1$immunesig <- gsub("\\.NA","",drugDEG_immunesig1$immunesig)
      return(drugDEG_immunesig1)
    })
    
    names(drugDEG_immunesig) <- drug_list2
    
    
    save(drugDEG_immunesig, file =paste0(filepath, "/",purity_method,"/",cancer,"_",dataset,"_",cellline, "_drugDEG_immunesig.Rdata"))
    
  }else{
    print('Your drug keyword can not be found in our database.')
  }
  
}



network_visualization <- function(df = NULL,
                                  cancer,
                                  purity_method,
                                  dataset,
                                  cellline,
                                  druglist,
                                  # drugkeywork = NULL, 
                                  immuneKeyword,
                                  PPIresource = c("STRING","inBio"),
                                  PPIconfidentScore = NULL,
                                  r = NULL,
                                  padj = NULL,
                                  p = NULL,
                                  drug_logFC = NULL,
                                  drug_padj = NULL,
                                  drug_p = NULL,
                                  filepath,
                                  immusig_path,
                                  resultpath){
  
    library(igraph)
    library(devtools)
    library(visNetwork)
    library(dplyr)
    
    if(PPIresource == "inBio"){
      load("r1_drug_immuneSig_correlation_analysis/06_results_summaryandplotting/data/PPI/inBio_PPI_net.RData")
      PPI_table = PPI_table[PPI_table$confidence.score >= PPIconfidentScore, ]
    }else{
      load("r1_drug_immuneSig_correlation_analysis/06_results_summaryandplotting/data/PPI/stringPPI.Rdata")
      names(PPI_table) = c("Protein1","Protein2",
                           "neighborhood","fusion",
                           "cooccurence","coexpression",
                           "experimental","database",
                           "textmining","combined_score")
      PPI_table$combined_score = PPI_table$combined_score/1000
      PPI_table = PPI_table[PPI_table$combined_score >= PPIconfidentScore, ]
    }
    
    if(purity_method == 'IHC'){
      load(paste0(immusig_path, "immunecell_TCGAgenes_pcor_pearson_IHC/pcor_pearson_",cancer,"_Primary_IHC.Rdata"))
      
    }else{
      load(paste0(immusig_path, "immunecell_TCGAgenes_result_pearson_", purity_method,"/pcor_pearson_",cancer,"_Primary_", purity_method,".Rdata"))
      
    }
    
    drugTargetMoA <- read.table("r1_drug_immuneSig_correlation_analysis/06_results_summaryandplotting/data/drugTargetMoA_merged.txt", sep = "\t", header = T)
    drugindexall <- read.table("r1_drug_immuneSig_correlation_analysis/06_results_summaryandplotting/data/drug_index_LINCS.txt", sep = "\t", header = T)
    drugTargetMoA <- merge(drugTargetMoA, drugindexall, by = c("pert_iname", "drug_index"), all =  T)

    # dt_lib = read.csv(paste0(load_dir, "/data/drugtarget.csv"), sep = ",",header = T, stringsAsFactors = F)
    if(is.null(df)){
      
      load(paste0(filepath, "/",purity_method,"/",cancer,"_",dataset,"_",cellline, "_drugDEG_immunesig.Rdata"))
      
      for(drug in druglist){

      dDEGimmuneAll = drugDEG_immunesig[[drug]]    
      dDEGimmune = dDEGimmuneAll[grep(immuneKeyword, dDEGimmuneAll$immunesig),]

      if(nrow(dDEGimmune)>0){

        dtMoA <- drugTargetMoA[drugTargetMoA$drug_index == drug, ]
        dDEGimmune_dtMoA = inner_join(dDEGimmune, dtMoA)
        
        dDEG = unique(dDEGimmune_dtMoA[c("pert_iname", "gene", "drug_logFC")])
        names(dDEG) = c("N1","N2","drug_logFC")
        
        
        geneimmune = unique(dDEGimmune_dtMoA[c("gene","immunesig","r")])
        names(geneimmune) = c("N1","N2","pearsonR")
        
        
        dt = unique(dDEGimmune_dtMoA[c("pert_iname","target")])
        names(dt) = c("N1","N2")
        dt$ifTarget = 1
        immune_cor <- cor_result[grep(immuneKeyword, names(cor_result))]
        
        ti <- lapply(immune_cor, function(x){
          if(!is.null(x)){
            # x=cor_result[[1]]
            x = data.frame(r = unlist(x$cor),
                           p.value = unlist(x$p.value),
                           adj.p = unlist(x$p.adj))
            x$N2 = row.names(x)
            if(is.null(r)){
              s1 <- summary(x[x$r > 0,]$r)[5]
              s2 <- summary(x[x$r < 0,]$r)[2]
              y <- unique(x[(x$r > s1 | x$r < s2) & x$adj.p < padj,])            
            }else{
              y <- unique(x[(abs(x$r) > r) & x$adj.p < padj,])            
            }
            y = inner_join(dt['N2'],y)
            
            return(y)
          }
        })
        ti <- do.call(rbind, ti)
        ti$N1 <- gsub("\\.\\d+","",rownames(ti))
        ti$N1 <- gsub("\\.NA","",ti$N1)
        ti <- ti[c('N1','N2','r')]
        names(ti) = c("N1","N2","pearsonR")
        
        dmoa = unique(dDEGimmune_dtMoA[c("pert_iname","MoA")])
        names(dmoa) = c("N1","N2")
        
        genenodes = union(dDEG$N2,dt$N2)
        ggi <- PPI_table[is.element(PPI_table$Protein1,genenodes),]
        ggi <- ggi[is.element(ggi$Protein2,genenodes),]
        
        if(nrow(ggi) > 0){
          if(PPIresource == "inBio"){
            
            names(ggi) = c("N1","N2","confidence.score", "initial.score")
          }else{
            names(ggi) = c("N1","N2",
                           "neighborhood","fusion",
                           "cooccurence","coexpression",
                           "experimental","database",
                           "textmining","combined_score")
          }
          
          df = merge(dDEG, geneimmune, by=c("N1","N2"), all = T)
          df = merge(df, dt, by=c("N1","N2"), all = T) 
          df = merge(df, dmoa, by=c("N1","N2"), all = T)
          df = merge(df, ggi, by=c("N1","N2"), all = T)
          
        }else{
          
          df = merge(dDEG, geneimmune, by=c("N1","N2"), all = T)
          df = merge(df, dt, by=c("N1","N2"), all = T) 
          df = merge(df, dmoa, by=c("N1","N2"), all = T)
        }
        
        if(nrow(ti)>0){
           df = merge(df, ti, by=c("N1","N2","pearsonR"), all = T) 
        }

        
      }else{

        dtMoA <- drugTargetMoA[drugTargetMoA$drug_index == drug, ]

        dt = unique(dtMoA[c("pert_iname","target")])
        names(dt) = c("N1","N2")
        dt$ifTarget = 1
        immune_cor <- cor_result[grep(immuneKeyword, names(cor_result))]
        
        ti <- lapply(immune_cor, function(x){
          if(!is.null(x)){
            # x=cor_result[[1]]
            x = data.frame(r = unlist(x$cor),
                           p.value = unlist(x$p.value),
                           adj.p = unlist(x$p.adj))
            x$N2 = row.names(x)
            if(is.null(r)){
              s1 <- summary(x[x$r > 0,]$r)[5]
              s2 <- summary(x[x$r < 0,]$r)[2]
              y <- unique(x[(x$r > s1 | x$r < s2) & x$adj.p < padj,])            
            }else{
              y <- unique(x[(abs(x$r) > r) & x$adj.p < padj,])            
            }
            y = inner_join(dt['N2'],y)
            
            return(y)
          }
        })
        ti <- do.call(rbind, ti)
        ti$N1 <- gsub("\\.\\d+","",rownames(ti))
        ti$N1 <- gsub("\\.NA","",ti$N1)
        ti <- ti[c('N1','N2','r')]
        names(ti) = c("N1","N2","pearsonR")
        
        dmoa = unique(dtMoA[c("pert_iname","MoA")])
        names(dmoa) = c("N1","N2")
        
        genenodes = unique(dt$N2)
        ggi <- PPI_table[is.element(PPI_table$Protein1,genenodes),]
        ggi <- ggi[is.element(ggi$Protein2,genenodes),]
        
        if(nrow(ggi) > 0){
          if(PPIresource == "inBio"){
            
            names(ggi) = c("N1","N2","confidence.score", "initial.score")
          }else{
            names(ggi) = c("N1","N2",
                           "neighborhood","fusion",
                           "cooccurence","coexpression",
                           "experimental","database",
                           "textmining","combined_score")
          }
          
          df = merge(dt, dmoa, by=c("N1","N2"), all = T)
          df = merge(df, ggi, by=c("N1","N2"), all = T)
          
        }else{
          df = merge(dt, dmoa, by=c("N1","N2"), all = T)
        }

        if(nrow(ti)>0){
           df = merge(df, ti, by=c("N1","N2"), all = T) 
        }


      }

        if(is.null(is.element(NULL, df$N1)))df = df[!is.null(df$N1),]
        if(is.null(is.element(NULL, df$N2)))df = df[!is.null(df$N1),]
        if(is.element(NA, df$N1))df = df[!is.na(df$N1),]
        if(is.element(NA, df$N2))df = df[!is.na(df$N2),]
        if(is.element("", df$N1))df = df[-which(is.element(df$N1, "")),]
        if(is.element("", df$N2))df = df[-which(is.element(df$N2, "")),]

      if(is.null(df) | nrow(df) < 1){
        print(paste0(drug, " can not be plotted for ", immuneKeyword," due to lack of data!"))

      }else{
        drugname = gsub("\\(", "", unique(dmoa$N1))
        drugname = gsub("\\)", "", drugname)


        dir.create(paste0(resultpath, cancer, "_", purity_method, "_", dataset, "_", cellline, "/" ))
        dir.create(paste0(resultpath, cancer, "_", purity_method, "_", dataset, "_", cellline, "/", drug, "_", drugname, "/"))
        dir.create(paste0(resultpath, cancer, "_", purity_method, "_", dataset, "_", cellline, "/", drug, "_", drugname, "/" , immuneKeyword, "/"))
        resultpath_drug = paste0(resultpath, cancer, "_", purity_method, "_", dataset, "_", cellline, "/" , drug, "_", drugname, "/", immuneKeyword, "/")

        write.table(df, paste0(resultpath_drug, "network_df.txt"), sep = "\t", row.names = F, quote = F)
        
        nodes = union(df$N1,df$N2)
        net = igraph::graph_from_data_frame(d=df, vertices = unique(nodes), directed = F)
        
        
        edge.color <- colorRampPalette(c("#D6D6D6","#383838"), alpha=TRUE)
        igraph::E(net)$color <- edge.color(igraph::ecount(net))
        
        
        igraph::V(net)$color = "#E7F3FD"
        igraph::V(net)[which(is.element(igraph::V(net)$name, union(dmoa$N1,dmoa$N2)))]$color = "#99FFFF"
        igraph::V(net)[which(is.element(igraph::V(net)$name, dt$N2))]$color = "#89D0F5"
        igraph::V(net)[grep(immuneKeyword, igraph::V(net)$name)]$color = "#FF9900"
        
        igraph::V(net)$shape = "box"
        igraph::V(net)[which(igraph::V(net)$name %in% union(dmoa$N1,dmoa$N2))]$shape = "ellipse"
        igraph::V(net)$group = "DEG"
        igraph::V(net)[which(igraph::V(net)$name %in%  unique(dt$N2))]$group = "Targets"
        igraph::V(net)[ grep(immuneKeyword, igraph::V(net)$name)]$group = "ImmuneSig"
        igraph::V(net)$font.size = 40
        
        data <- toVisNetworkData(net)
        
        visNetwork(nodes = data$nodes, edges = data$edges, width = "90%", height = "95vh")%>%
          visNodes(size = 40)%>%
          # visHierarchicalLayout(direction = "LR", levelSeparation = 500)%>% 
          visIgraphLayout(layout = "layout_on_sphere",type = "full",randomSeed = 123) %>%
          visOptions(highlightNearest = list(enabled = TRUE,  hideColor = "lightgrey", hover = T),
                     nodesIdSelection =list(enabled = TRUE), selectedBy = "group") %>%
          # visConfigure(enabled = TRUE) %>%
          addFontAwesome() %>%
          visGroups(groupname = "Targets", color = "#89D0F5")%>%
          visGroups(groupname = "ImmuneSig", color = "#FF9900")%>%
          visLegend() %>%
          visInteraction(navigationButtons = TRUE) %>%
          visOptions(manipulation = TRUE) %>% 
          visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
          visSave(file = paste0(resultpath_drug, "networkAnalysis.html"))
        
      }

    }
  }
}

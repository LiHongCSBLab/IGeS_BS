ITS_g_cancer_plot <- function(ITS_set, 
                               type = c("s","r"),
                               genelist, 
                               enrichRes,
                               savepath,
                               filename){

  library(corrplot)
  library(UpSetR)
  library(reshape2)
  library(data.table)
  ITS_dg <- lapply(ITS_set, function(x){
    lapply(x, function(y){intersect(y, genelist)})
  }) 
  
  ITS_dg = lapply(ITS_dg,function(x){x[!is.na(names(x))]})
  
  ITS_dc_list = list()
  for(i in 1:length(ITS_dg)){
    # i = 1
    names(ITS_dg[[i]])
    listinput <- ITS_dg[[i]]
    df <- lapply(as.list(names(listinput)), function(x){
        if(length(listinput[[x]])>0){
        data.frame(immune_sig = x, gene_name = listinput[[x]])
        }
    })
    df <- do.call(rbind, df)
    gene_df <- data.frame(table(df$gene_name))
    names(gene_df) <- c("gene_name", "num_ITS")
    gene_df$num_ITS_rate <- gene_df$num_ITS / length(names(listinput))
    gene_df$cancer_type <- names(ITS_dg)[i]
    ITS_dc_list[[i]] <- gene_df
  }

  names(ITS_dc_list) = names(ITS_dg)
  ITS_dc <- do.call(rbind, ITS_dc_list)

  ITS_dc_df <- dcast(ITS_dc, gene_name~cancer_type, value.var = "num_ITS_rate")
  rownames(ITS_dc_df) <- ITS_dc_df$gene_name
  ITS_dc_df <- as.matrix(ITS_dc_df[-1])


  write.csv(ITS_dc_df, paste0(savepath, filename, ".csv"), quote = F)
  
  library(corrplot)
  if(type == 's'){
    colors <- colorRampPalette(c("blue", "white", "#B33333"))(20)
  }else if(type == 'r'){
    colors <- colorRampPalette(c("blue", "white", "#338080"))(20)
  }

  pdf(paste0(savepath, filename, ".pdf"),
      width = 10,
      height = 10)
  corrplot(ITS_dc_df, method = "pie", 
           col = colors, na.label = " ",
           tl.col = 'black')
  dev.off()
  
  ITS_dc_df2 <- ITS_dc_df
  ITS_dc_df2$gene_name <- row.names(ITS_dc_df2)
  
  enrichRes_geneID <- lapply(as.list(enrichRes$geneID), function(x){
    data.frame(ENTREZID = unlist(strsplit(x, "/")))
  })
  names(enrichRes_geneID) <- enrichRes$ID
  enrichRes_geneID <- do.call(rbind, enrichRes_geneID)
  enrichRes_geneID$ID <-  unlist(lapply(strsplit(rownames(enrichRes_geneID), "\\."), function(x)x[[1]]))

  enrichRes2 <- merge(enrichRes[c('ID', 'Description')], enrichRes_geneID, by = 'ID')

  genename <- bitr(genelist,
                    fromType = "SYMBOL", 
                    toType = c("ENTREZID", "SYMBOL"), 
                    OrgDb = org.Hs.eg.db)
  enrichRes2 <- unique(merge(enrichRes2, genename, by = 'ENTREZID'))
  enrichRes2 <- enrichRes2[order(enrichRes2$ID), ]

  write.csv(enrichRes2, paste0(savepath, filename, "_geneKEGGannot.csv"), quote = F)
  
}

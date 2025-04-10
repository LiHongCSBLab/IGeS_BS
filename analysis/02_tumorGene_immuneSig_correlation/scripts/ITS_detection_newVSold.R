setwd("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/02_tumorGene_immuneSig_correlation/results/")


resultpath = paste0("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/02_tumorGene_immuneSig_correlation/results/immunecell_TCGAgenes_result_pearson_TUMERIC/")
file_name = 'pcor_pearson_COAD_Primary_TUMERIC.Rdata'
method = "pearson"
num_gene = 400
num_gene_p = 200
num_gene_n = 200

file = paste0(resultpath, file_name)
  cancer = unlist(strsplit(file_name, split = "_"))[3]
  sample_type = unlist(strsplit(file_name, split = "_"))[4]
  tumor_purity_method = gsub('.Rdata', '', unlist(strsplit(file_name, split = "_"))[5])
  method = "pearson"
 
  load(file)
  

  genelist_p <- sapply(cor_result, function(x){
    #summary(x$r)
    # x= cor_result[[87]]
    x = data.frame(r = unlist(x$cor),
                   p.value = unlist(x$p.value),
                   adj.p = unlist(x$p.adj))
    x <- x[x$adj.p < 0.05,]
    s1 <- summary(x[x$r > 0,]$r)[5]
    s2 <- summary(x[x$r < 0,]$r)[2]
    y <- unique(x[(x$r > s1 | x$r < s2) & x$adj.p < 0.05,])
    y <- y[!is.na(y$r),]
    y$gene_name = rownames(y)
    
    #length(y)
    if(nrow(y) ==0){
      g_p = NULL
    }else{
      if(nrow(y) <= num_gene_p){
        
        g_p <- y[y$r > 0,]$gene_name
        
      }else{
        y_p <- y[order(y$r, decreasing = T),]
        y_p <- y_p[y_p$r > 0,]
        g_p <- y_p$gene_name[1:num_gene_p]
      }
    }
    
    # g_p <- paste(g_p, collapse = "\t")
    return(g_p)
  })
  
  genelist_n <- sapply(cor_result, function(x){
    # summary(x$r)
    # x= cor_result[[1]]
    x = data.frame(r = unlist(x$cor),
                   p.value = unlist(x$p.value),
                   adj.p = unlist(x$p.adj))
    
    x <- x[x$adj.p < 0.05,]
    s1 <- summary(x[x$r > 0,]$r)[5]
    s2 <- summary(x[x$r < 0,]$r)[2]
    y <- unique(x[(x$r > s1 | x$r < s2) & x$adj.p < 0.05,])
    y <- y[!is.na(y$r),]
    y$gene_name = rownames(y)
    
    # length(y)
    if(nrow(y) ==0){
      g_n = NULL
    }else{
      if(nrow(y) <= num_gene_n) {
        g_n <- y[y$r < 0,]$gene_name
      }else{
        y_n <- y[order(y$r, decreasing = T),]
        y_n <- y_n[y_n$r < 0,]
        g_n <- y_n$gene_name[1:num_gene_n]
      }
    }
    # g_n <- paste(g_n, collapse = "\t")
    return(g_n)
  })

  


genelist_p_new = genelist_p
genelist_n_new = genelist_n




resultpath = paste0("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/02_tumorGene_immuneSig_correlation/results/immunecell_TCGAgenes_result_pearson_TUMERIC_old/")
file_name = 'pcor_pearson_COAD_Primary_TUMERIC.Rdata'
method = "pearson"
num_gene = 400
num_gene_p = 200
num_gene_n = 200


  file = paste0(resultpath, file_name)
  cancer = unlist(strsplit(file_name, split = "_"))[3]
  sample_type = unlist(strsplit(file_name, split = "_"))[4]
  tumor_purity_method = gsub('.Rdata', '', unlist(strsplit(file_name, split = "_"))[5])
  method = "pearson"
 
  load(file)
  

  genelist_p <- sapply(cor_result, function(x){
    #summary(x$r)
    # x= cor_result[[87]]
    x = data.frame(r = unlist(x$cor),
                   p.value = unlist(x$p.value),
                   adj.p = unlist(x$p.adj))
    x <- x[x$adj.p < 0.05,]
    s1 <- summary(x[x$r > 0,]$r)[5]
    s2 <- summary(x[x$r < 0,]$r)[2]
    y <- unique(x[(x$r > s1 | x$r < s2) & x$adj.p < 0.05,])
    y <- y[!is.na(y$r),]
    y$gene_name = rownames(y)
    
    #length(y)
    if(nrow(y) ==0){
      g_p = NULL
    }else{
      if(nrow(y) <= num_gene_p){
        
        g_p <- y[y$r > 0,]$gene_name
        
      }else{
        y_p <- y[order(y$r, decreasing = T),]
        y_p <- y_p[y_p$r > 0,]
        g_p <- y_p$gene_name[1:num_gene_p]
      }
    }
    
    # g_p <- paste(g_p, collapse = "\t")
    return(g_p)
  })
  
  genelist_n <- sapply(cor_result, function(x){
    # summary(x$r)
    # x= cor_result[[1]]
    x = data.frame(r = unlist(x$cor),
                   p.value = unlist(x$p.value),
                   adj.p = unlist(x$p.adj))
    
    x <- x[x$adj.p < 0.05,]
    s1 <- summary(x[x$r > 0,]$r)[5]
    s2 <- summary(x[x$r < 0,]$r)[2]
    y <- unique(x[(x$r > s1 | x$r < s2) & x$adj.p < 0.05,])
    y <- y[!is.na(y$r),]
    y$gene_name = rownames(y)
    
    # length(y)
    if(nrow(y) ==0){
      g_n = NULL
    }else{
      if(nrow(y) <= num_gene_n) {
        g_n <- y[y$r < 0,]$gene_name
      }else{
        y_n <- y[order(y$r, decreasing = T),]
        y_n <- y_n[y_n$r < 0,]
        g_n <- y_n$gene_name[1:num_gene_n]
      }
    }
    # g_n <- paste(g_n, collapse = "\t")
    return(g_n)
  })

  

genelist_p_old = genelist_p
genelist_n_old = genelist_n


df = data.frame(i =  names(genelist_p_new), new_p=0, old_p=0, new_n=0, old_n=0, p=0, n=0)
for(i in seq(length(genelist_p_new))){
    pn = genelist_p_new[[i]]
    po = genelist_p_old[[i]]
    nn = genelist_n_new[[i]]
    no = genelist_n_old[[i]]
    df[df$i==names(genelist_p_new)[i],]$new_p = (length(pn[!is.na(pn)]))
    df[df$i==names(genelist_p_new)[i],]$old_p = (length(po[!is.na(po)]))
    df[df$i==names(genelist_p_new)[i],]$new_n = (length(nn[!is.na(nn)]))
    df[df$i==names(genelist_p_new)[i],]$old_n = (length(no[!is.na(no)]))
    df[df$i==names(genelist_p_new)[i],]$p = (length(intersect(pn[!is.na(pn)], po[!is.na(po)])))
    df[df$i==names(genelist_p_new)[i],]$n = (length(intersect(nn[!is.na(nn)], no[!is.na(no)])))
}

View(df)

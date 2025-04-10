df_for_plot <- function(cancer,
                        tumor_subtype,
                        purity_method = "TUMERIC",
                        dataset,
                        sign,
                        num_gene,
                        datatype,
                        immunesig,
                        ACAT_Pval = 0.05,
                        # work_path, 
                        filepath,
                        savepath) {
  
  dir.create(paste0(savepath, "/"))
  dir.create(paste0(savepath, "/", dataset,"/"))
  dir.create(paste0(savepath, "/", dataset,"/",datatype))
  dir.create(paste0(savepath, "/", dataset,"/",datatype,"/",cancer,"_", tumor_subtype, "/"))
  
  # filepath =  paste0(work_path, "03_drug_immuneSig_enrichment/results_", purity_method,"_meta_oe_original_final/")

  if(datatype == "978genes"){
    resPath <- paste0(filepath, "/drug_immunesig_",purity_method,"_GSEA_result_summarized_978genes/")
    resPath <- paste0(resPath, "/", dataset,"/", cancer, "_", tumor_subtype, "/", sign, "_", num_gene, "/")
  }else{
    resPath <- paste0(filepath, "/drug_immunesig_",purity_method,"_GSEA_result_summarized/")
    resPath <- paste0(resPath, "/", dataset,"/", cancer, "_", tumor_subtype, "/", sign, "_", num_gene, "/")
    
  }
  
  files <- list.files(resPath)
  
  druginfo <- gsub("_result_merged.csv", "", files)
  # memory.limit(size=20000)
  
  df = lapply(as.list(files), function(file){
    # for(file in files){
    if(file.size(paste0(resPath, file)) > 0){
      df <- read.csv(paste0(resPath, file))
      # df <- df[df$mPvalue_fisher < 0.1, ]
      df <- df[df$mPvalue_ACAT < ACAT_Pval, ]
      # df <- df[df$mPvalue_ACAT < 1, ]
      df <- df[c("drug_index","sig","mNES",
                 "mPvalueNaive","mPvalue_fisher","mPvalue_ACAT")]
      return(as.data.frame(df))
    }
  })
  df = do.call(rbind, df)
  
  names(df) = c(c("drug_index","immune_sig","mNES",
                  "mPvalueNaive","mPvalue_fisher","mPvalue_ACAT"))
  
  
  df$immune_sig = tolower(df$immune_sig)
  
  df_immunesig = dplyr::inner_join(df, immunesig)
  immunesig2 = unique(df_immunesig[c("immune_sig", "setindex")])
  immunesig2 <- immunesig2[order(immunesig2$setindex),]
  
  write.csv(immunesig2,  paste0(savepath, "/",dataset,"/",datatype,"/",cancer,"_", tumor_subtype, "/", sign, "_", num_gene,"_immunesig2_",datatype,".csv"), row.names = F, quote = F)
  
  df_immunesig[grep("01" ,df_immunesig$setindex),]$setindex = 1
  df_immunesig[grep("02", df_immunesig$setindex),]$setindex = 2
  
  
  if(nrow(df_immunesig) > 0 ){
    
    # pathway bubble plot
    node_color <- colorRampPalette(colors = c("steelblue", "white", "#DC143C"))
    df_immunesig <- df_immunesig[order(df_immunesig$setindex),]
    df_immunesig$immune_sig <- factor(df_immunesig$immune_sig, levels = unique(immunesig2$immune_sig))
    df_immunesig$log10_p <- -log10(df_immunesig$mPvalue_ACAT)
    if(nrow(df_immunesig[which(df_immunesig$log10_p =="Inf"),])!=0){df_immunesig[which(df_immunesig$log10_p =="Inf"),]$log10_p <- 10}
    
    
    immunesig2 = immunesig2[order(immunesig2$immune_sig),]
    immunesig2 = immunesig2[order(immunesig2$setindex),]
    
    df_immunesig$immune_sig <- factor(df_immunesig$immune_sig, levels = unique(immunesig2$immune_sig))
    
    p1 <- ggplot(df_immunesig, aes(x = immune_sig, y = drug_index)) + 
      geom_point(aes(size = log10_p, color = mNES , )) +
      scale_size(range = c(2, 8)) +
      scale_colour_gradientn(colours = node_color(500), na.value = "grey90", values = seq(0, 1, length.out = 500)) +
      scale_size_area(max_size = 5) 
    p1 <- p1 + theme_bw() +
      theme(axis.text.x  =  element_text(angle=-90, vjust=0.5,hjust = 0, size=8),
            axis.text.y = element_text(size = 15),
            axis.title.x = element_blank()) # +
    # ggtitle("Plot1a: gene positively correlated with immune signatures \n            immune signatures increase risk")
    # p1
    
    plot_savepath <- paste0(savepath, "/",dataset,"/",datatype,"/",cancer,"_", tumor_subtype, "/", sign, "_", num_gene,"_",datatype,".pdf")
    ggsave(plot_savepath, p1, width = 20, height = 35)
    
    
    # df_immunesig3 = df_immunesig[-grep("oe_",df_immunesig$immune_sig),]
    # df_immunesig3 = df_immunesig[!is.element(df_immunesig$immune_sig,
    #                                          c("dysfunction_score_tide",
    #                                            "exclusion_score_tide",
    #                                            #  "mhc1_21978456_sig160",
    #                                            #  "mhc.i_19272155_sig160",
    #                                            #  "mhcii_tcga_caf_immunecls",
    #                                            #  "cd8a_sig160",
    #                                            #  "cd8_cd68_ratio_sig160",
    #                                            #  "cd8.t.cells_sig160",
    #                                            #  "mcd3_cd8_21214954_sig160",
    #                                            "az_ips")),]
    # df_immunesig3 = df_immunesig3[-grep("xcell",df_immunesig3$immune_sig),]
    
    df_immunesig3 = df_immunesig
    p2 <- ggplot(df_immunesig3, aes(x = immune_sig, y = drug_index)) + 
      geom_point(aes(size = log10_p, color = mNES , )) +
      scale_size(range = c(2, 8)) +
      scale_colour_gradientn(colours = node_color(500), na.value = "grey90", values = seq(0, 1, length.out = 500)) +
      scale_size_area(max_size = 5) 
    p2 <- p2 + theme_bw() +
      theme(axis.text.x  =  element_text(angle=-90, vjust=0.5,hjust = 0, size=8),
            axis.text.y = element_text(size = 15),
            axis.title.x = element_blank()) # +
    # ggtitle("Plot1a: gene positively correlated with immune signatures \n            immune signatures increase risk")
    # p2
    plot_savepath <- paste0(savepath, "/",dataset,"/",datatype,"/",cancer,"_", tumor_subtype, "/", sign, "_", num_gene,"_",datatype,"2.pdf")
    ggsave(plot_savepath, p2, width = 20, height = 35)
    
    write.csv(df_immunesig,  paste0(savepath, "/",dataset,"/",datatype,"/",cancer,"_", tumor_subtype, "/", sign, "_", num_gene,"_",datatype,".csv"), row.names = F, quote = F)
    # pdf(file=plot_savepath, width = 20, height = 8)
    # p
    # dev.off()
    # rm("p1"); rm("p2"); rm("p")
    
  }
  return(df)
  
  
}





drugITSRanker <- function(df_immunesig,
                          cancer,
                          tumor_subtype, 
                          weight_type,
                          method = c("glmnet_weight", "simple", "sum_mean"),
                          type = c("unweighted", "weighted", "binary_weighted"),
                          # mergeITS,
                          # mergeITS_type = c("Type", "Index"),
                          d1 = NULL,
                          d2 = NULL,
                          savepath){
  
  require(ggplot2)
  require(GGally)
  library(dplyr)
  library(plot3D)
  library(fmsb)
  library(ggrepel)
  library(fs)
  drugTargetMoA_file=path("Data", "drugTargetMoA_merged.txt")
  drugTargetMoA <- read.table(drugTargetMoA_file, sep = "\t", header = T)
  drugindexall_file=path("Data", "drug_index_LINCS.txt")
  drugindexall <- read.table(drugindexall_file, sep = "\t", header = T)
  
  # res_save <- read.csv(paste0(savepath, "drugITSRank/", weight_type,"_",type, "_",method,"_res.txt"),sep='\t', header = T)
  # savepath = paste0(savepath, type,'/')
  savepath = paste0(savepath, "drugITSRank/")
  dir.create(savepath)

  # plot(df_immunesig[df_immunesig$drug_index == 'drug888',]$mNES_diff, 
  # df_immunesig[df_immunesig$drug_index == 'drug888',]$weight,
  # xlim = c(-2,2), ylim = c(-2,2))
  # abline(h=0)
  # abline(v = 0)
  
  
  if(type == "weighted"){
    
    df_immunesig$mNES_diff_weighted = df_immunesig$mNES_diff * df_immunesig$weight

  }else if(type == "binary_weighted"){
    
    df_immunesig$mNES_diff_binary = 0
    df_immunesig[df_immunesig$mNES_diff > 0, ]$mNES_diff_binary = 1
    df_immunesig[df_immunesig$mNES_diff < 0, ]$mNES_diff_binary = -1
    
    df_immunesig$mNES_diff_weighted = df_immunesig$mNES_diff_binary * df_immunesig$weight
    
    # df_s = df_immunesig[grep("01",df_immunesig$setindex),]
    # df_r = df_immunesig[grep("02",df_immunesig$setindex),]
    # df_s = df_s[c("drug_index","immune_sig","mNES_diff_weighted")]
    # names(df_s) = c("drug_index","immune_sig_s","mNES_diff_s")
    # df_r = df_r[c("drug_index","immune_sig","mNES_diff_weighted")]
    # names(df_r) = c("drug_index","immune_sig_r","mNES_diff_r")
    
  }else if(type == "binary_unweighted"){
    
    df_immunesig$mNES_diff_binary = 0
    df_immunesig[df_immunesig$mNES_diff > 0, ]$mNES_diff_binary = 1
    df_immunesig[df_immunesig$mNES_diff < 0, ]$mNES_diff_binary = -1
    
    df_immunesig$mNES_diff_weighted = df_immunesig$mNES_diff_binary # * df_immunesig$weight
    
    # df_s = df_immunesig[grep("01",df_immunesig$setindex),]
    # df_r = df_immunesig[grep("02",df_immunesig$setindex),]
    # df_s = df_s[c("drug_index","immune_sig","mNES_diff_weighted")]
    # names(df_s) = c("drug_index","immune_sig_s","mNES_diff_s")
    # df_r = df_r[c("drug_index","immune_sig","mNES_diff_weighted")]
    # names(df_r) = c("drug_index","immune_sig_r","mNES_diff_r")
    
  }else if(type == "unweighted"){
    df_immunesig$mNES_diff_weighted = df_immunesig$mNES_diff

    # df_s = df_s[c("drug_index","immune_sig","mNES_diff")]
    # names(df_s) = c(c("drug_index","immune_sig_s","mNES_diff_s"))
    # df_r = df_r[c("drug_index","immune_sig","mNES_diff")]
    # names(df_r) = c(c("drug_index","immune_sig_r","mNES_diff_r"))
    
  }

    
  # if(mergeITS){

  # immuneSigInfoPath=path("Data", "immune_sig_hieracy_new.xlsx")
  # immunesigInfo <- as.data.frame(readxl::read_excel(immuneSigInfoPath, sheet = 1))
  #   immunesigInfo <- as.data.frame(readxl::read_excel(paste0(immuneSigInfoPath, "immune_sig_hieracy_new.xlsx"), sheet = 1))
  #   immunesigInfo <- immunesigInfo[c(1,4,5:9)]
  #   immunesigInfo$immune_sig = tolower(immunesigInfo$immune_sig)
  #   immunesigInfo = immunesig_name_unify(immunesigInfo)

  #   df_s = df_immunesig[grep("01",df_immunesig$setindex),]
  #   df_s = df_s[c("drug_index","immune_sig","mNES_diff_weighted")]
  #   df_s = inner_join(df_s, immunesigInfo)
    
  #   df_s2 = lapply(as.list(unique(df_s$drug_index)), function(d){
  #     df1 = df_s[df_s$drug_index == d, ]
  #     tmp = aggregate(df1['mNES_diff_weighted'], by = list(df1[,mergeITS_type]), mean)
  #     tmp$drug_index = d 
  #     return(tmp)
  #   })
  #   df_s2 = do.call(rbind, df_s2) 
  #   df_s = df_s2[c("drug_index","Group.1","mNES_diff_weighted")]
  #   names(df_s) = c("drug_index","immune_sig_s","mNES_diff_s")


  #   df_r = df_immunesig[grep("02",df_immunesig$setindex),]
  #   df_r = df_r[c("drug_index","immune_sig","mNES_diff_weighted")]
  #   df_r = inner_join(df_r, immunesigInfo)
  #   df_r2 = lapply(as.list(unique(df_r$drug_index)), function(d){
  #     df1 = df_r[df_r$drug_index == d, ]
  #     tmp = aggregate(df1['mNES_diff_weighted'], by = list(df1[,mergeITS_type]), mean)
  #     tmp$drug_index = d 
  #     return(tmp)
  #   })
  #   df_r2 = do.call(rbind, df_r2) 
  #   df_r = df_r2[c("drug_index","Group.1","mNES_diff_weighted")]
  #   names(df_r) = c("drug_index","immune_sig_r","mNES_diff_r")

  # }else{

    df_s = df_immunesig[grep("01",df_immunesig$setindex),]
    df_r = df_immunesig[grep("02",df_immunesig$setindex),]
    
    df_s = df_s[c("drug_index","immune_sig","mNES_diff_weighted")]
    names(df_s) = c("drug_index","immune_sig_s","mNES_diff_s")
    df_r = df_r[c("drug_index","immune_sig","mNES_diff_weighted")]
    names(df_r) = c("drug_index","immune_sig_r","mNES_diff_r")

  # }
  
  
  df = merge(unique(df_s[c("drug_index","mNES_diff_s")]),
             unique(df_r[c("drug_index","mNES_diff_r")]),
             by= "drug_index")
  
  mat_s <-  reshape(df_s,  timevar = "immune_sig_s",
                    idvar = "drug_index",direction = "wide")
  rownames(mat_s) <- mat_s$drug_index
  mat_s <- mat_s[-1]
  
  mat_r <-  reshape(df_r,  timevar = "immune_sig_r",
                    idvar = "drug_index",direction = "wide")
  rownames(mat_r) <- mat_r$drug_index
  mat_r <- mat_r[-1]

  # Rebuild matrices with full row set, filling missing rows with NA
  remake_matrix <- function(mat, all_rows) {
    new_mat <- matrix(NA, nrow = length(all_rows), ncol = ncol(mat),
                      dimnames = list(all_rows, colnames(mat)))
    common_rows <- intersect(rownames(mat), all_rows)
    new_mat[common_rows, ] <- mat[common_rows, ]
    return(new_mat)
  }

  # Apply to both matrices
  all_rows <- union(rownames(mat_s), rownames(mat_r))
  mat_s <- as.data.frame(remake_matrix(as.matrix(mat_s), all_rows))
  mat_r <- as.data.frame(remake_matrix(as.matrix(mat_r), all_rows))

  # mat_s2 <- as.matrix(mat_s[is.element(rownames(mat_s), intersect(rownames(mat_s), rownames(mat_r) )),])
  # mat_r2 <- as.matrix(mat_r[is.element(rownames(mat_r), intersect(rownames(mat_s), rownames(mat_r) )),])
  mat_s2 <- as.matrix(mat_s[intersect(rownames(mat_s), rownames(mat_r)),])
  rownames(mat_s2) = intersect(rownames(mat_s), rownames(mat_r))
  colnames(mat_s2) = names(mat_s)
  mat_r2 <- as.matrix(mat_r[intersect(rownames(mat_s), rownames(mat_r)),])
  rownames(mat_r2) = intersect(rownames(mat_s), rownames(mat_r))
  colnames(mat_r2) = names(mat_r)
  
  mat_s2[is.na(mat_s2)] = 0
  mat_r2[is.na(mat_r2)] = 0
  
  
  
  if(method == "glmnet_weight"){
    library(plotrix)
    # dir.create(paste0(savepath,"/",method,"/"))
    
    CCAdata1_Z <- as.matrix(mat_s2)
    CCAdata2_Z <- as.matrix(mat_r2)
    # CCAdata1_Z <- as.matrix(scale(mat_s2, center = T))
    # CCAdata2_Z <- as.matrix(scale(mat_r2, center = T))
    
    
    if(length(which(is.na(CCAdata1_Z) == T))!=0)CCAdata1_Z[is.na(CCAdata1_Z)] <- 0
    if(length(which(is.na(CCAdata1_Z) == T))!=0)CCAdata2_Z[is.na(CCAdata2_Z)] <- 0
    
    u = as.data.frame(apply(CCAdata1_Z, 1, sum)); names(u) = "u"; u$drug_index = rownames(u)
    v = as.data.frame(apply(CCAdata2_Z, 1, sum)); names(v) = "v"; v$drug_index = rownames(v)
    
    
    res = merge(u, v, by = "drug_index")
    rownames(res) = res$drug_index
    res = res[c("u", "v")]
    res$drug_index <- rownames(res)    
    res = res[c("drug_index","u","v")]
    
    res_save = dplyr::inner_join(drugindexall, res)
    # res_save[res_save$u >0 & res_save$v <0 ,]
    res_save$score = res_save$u + res_save$v
    res_save = res_save[order(res_save$score, decreasing = T), ]
    res_save$rank = 1:nrow(res_save)
    # res_save[res_save$u>0 & res_save$v >0,]  
    
    write.table(res_save, paste0(savepath, weight_type,"_",type, "_",method,"_res.txt"),sep='\t',
                quote = F, row.names = F)    
    
    dat = res_save %>%
      arrange(desc(score)) %>%
      mutate("x" = row_number())


    res1 = head(dat, nrow(dat) * 0.1)  
    res1 = res1[res1$u>0 & res1$v >0,]  
    range = ceiling(max(dat$score))
    
    p_rank <- ggplot(dat, aes(x = x, y = score))+
      geom_point(color="#1762f7", size = 0.5)+
      theme_bw()+
      labs(x = "Ranked drugs", y = "Potential")+
      geom_text_repel(inherit.aes = F, 
                      data = res1, 
                      aes(x = x, y = score, label = pert_iname, color = 'red'),
                      max.overlaps = getOption("ggrepel.max.overlaps", default = 200))+
      scale_y_continuous(breaks = seq(-range,range,2)) +
      ggtitle(paste0(cancer, "_", tumor_subtype, "_", dataset))
    
    
    ggsave(paste0(savepath, weight_type,"_",type, "_",method,"_res.pdf"),p_rank, width = 5, height = 12)
    
    # par(mfrow = c(1,1))
    # plot(dat$x, dat$score, pch=20,
    #      ylab = "Immune impact score")
    # abline(v = nrow(dat) * 0.01, lty = 6)
    # abline(v = nrow(dat) * 0.05, lty = 6)
    # abline(v = nrow(dat) * 0.1, lty = 6)
    # points(res1$x, res1$score, pch=21, col="red", bg="red")
    # text(res1$x, res1$score, paste0(res1$drug_index),
    #      cex = 1, pos = 4, col = "red")
    
    # range =max(abs(res$u),abs(res$v))
    # plot(res$u, res$v, pch=20,
    #      xlim=c(-range, range), ylim=c(-range,range),
    #      xlab ="sensitive score", ylab = "resistance score")
    # abline(v = 0, h = 0, lty = 6)
    # abline(coef = c(0,1), lty = 6)
    # abline(coef = c(0,-1), lty = 6)
    # points(dcca$u, dcca$v, pch=21, col="red", bg="red")
    # arrows(x0 = 0, y0 = 0, x1 = dcca$u, y1 = dcca$v, lwd = 1.5, col = "red", length = 0.1)
    # text(dcca$u, dcca$v, gsub("drug", "D", rownames(dcca)),
    #      cex = 1, pos = 4, col = "red")
    # draw.circle(0,0,0.5, lty = 6)
    # draw.circle(0,0,1, lty = 6)
    
    
  }else if(method == "simple"){
    library(plotrix)
    # dir.create(paste0(savepath,"/",method,"/"))
    
    CCAdata1_Z <- as.matrix(mat_s2)
    CCAdata2_Z <- as.matrix(mat_r2)
    # CCAdata1_Z <- as.matrix(scale(mat_s2, center = T))
    # CCAdata2_Z <- as.matrix(scale(mat_r2, center = T))
    
    
    if(length(which(is.na(CCAdata1_Z) == T))!=0)CCAdata1_Z[is.na(CCAdata1_Z)] <- 0
    if(length(which(is.na(CCAdata1_Z) == T))!=0)CCAdata2_Z[is.na(CCAdata2_Z)] <- 0
    
    u = as.data.frame(apply(CCAdata1_Z, 1, mean)); names(u) = "u"; u$drug_index = rownames(u)
    v = as.data.frame(apply(CCAdata2_Z, 1, mean)); names(v) = "v"; v$drug_index = rownames(v)
    
    
    res = merge(u, v, by = "drug_index")
    rownames(res) = res$drug_index
    res = res[c("u", "v")]
    
    res$drug_index <- rownames(res)    
    res = res[c("drug_index","u","v")]
    
    res_save = inner_join(drugindexall, res)
    # res_save[res_save$u >0 & res_save$v <0 ,]
    res_save$score = (res_save$u - res_save$v) 
    res_save = res_save[order(res_save$score, decreasing = T), ]
    res_save$rank = 1:nrow(res_save)
    
    write.table(res_save, paste0(savepath, weight_type,"_",type, "_",method,"_res.txt"),sep='\t',
                quote = F, row.names = F)    
    
    dat = res_save %>%
      arrange(desc(score)) %>%
      mutate("x" = row_number())
    res1 = head(dat, nrow(dat) * 0.1)  
    res1 = res1[res1$u>0 & res1$v < 0,]  
    range = ceiling(max(dat$score))
    
    p_rank <- ggplot(dat, aes(x = x, y = score))+
      geom_point(color="#1762f7", size = 0.5)+
      theme_bw()+
      labs(x = "Ranked drugs", y = "Potential")+
      geom_text_repel(inherit.aes = F, 
                      data = res1, 
                      aes(x = x, y = score, label = pert_iname, color = 'red'),
                      max.overlaps = getOption("ggrepel.max.overlaps", default = 200))+
      scale_y_continuous(breaks = seq(-range,range,2)) +
      ggtitle(paste0(cancer, "_", tumor_subtype, "_", dataset))
    ggsave(paste0(savepath, weight_type,"_",type, "_",method,"_res.pdf"),p_rank, width = 5, height = 12)
    
    # par(mfrow = c(1,1))
    # plot(dat$x, dat$score, pch=20,
    #      ylab = "Immune impact score")
    # abline(v = nrow(dat) * 0.01, lty = 6)
    # abline(v = nrow(dat) * 0.05, lty = 6)
    # abline(v = nrow(dat) * 0.1, lty = 6)
    # points(res1$x, res1$score, pch=21, col="red", bg="red")
    # text(res1$x, res1$score, paste0(res1$drug_index,"\n", res1$pert_iname),
    #      cex = 1, pos = 4, col = "red")
    
    # range =max(abs(res$u),abs(res$v))
    # plot(res$u, res$v, pch=20,
    #      xlim=c(-range, range), ylim=c(-range,range),
    #      xlab ="sensitive score", ylab = "resistance score")
    # abline(v = 0, h = 0, lty = 6)
    # abline(coef = c(0,1), lty = 6)
    # abline(coef = c(0,-1), lty = 6)
    # points(dcca$u, dcca$v, pch=21, col="red", bg="red")
    # arrows(x0 = 0, y0 = 0, x1 = dcca$u, y1 = dcca$v, lwd = 1.5, col = "red", length = 0.1)
    # text(dcca$u, dcca$v, gsub("drug", "D", rownames(dcca)),
    #      cex = 1, pos = 4, col = "red")
    # draw.circle(0,0,0.5, lty = 6)
    # draw.circle(0,0,1, lty = 6)
    
    
  }else if(method == "sum_mean"){
    
    library(plotrix)
    # dir.create(paste0(savepath,"/",method,"/"))
    
    CCAdata1_Z <- as.matrix(mat_s2)
    CCAdata2_Z <- as.matrix(mat_r2)
    # CCAdata1_Z <- as.matrix(scale(mat_s2, center = T))
    # CCAdata2_Z <- as.matrix(scale(mat_r2, center = T))
    
    num_ITS = ncol(CCAdata1_Z) + ncol(CCAdata2_Z)
    if(length(which(is.na(CCAdata1_Z) == T))!=0)CCAdata1_Z[is.na(CCAdata1_Z)] <- 0
    if(length(which(is.na(CCAdata1_Z) == T))!=0)CCAdata2_Z[is.na(CCAdata2_Z)] <- 0
    
    u = as.data.frame(apply(CCAdata1_Z, 1, sum)); names(u) = "u"; u$drug_index = rownames(u)
    v = as.data.frame(apply(CCAdata2_Z, 1, sum)); names(v) = "v"; v$drug_index = rownames(v)
    
    
    res = merge(u, v, by = "drug_index")
    rownames(res) = res$drug_index
    res = res[c("u", "v")]
    res1 = res[res$u>0 & res$v <0,]  
    res$drug_index <- rownames(res)    
    res = res[c("drug_index","u","v")]
    
    res_save = inner_join(drugindexall, res)
    # res_save[res_save$u >0 & res_save$v <0 ,]
    res_save$score = (res_save$u - res_save$v) / num_ITS
    res_save = res_save[order(res_save$score, decreasing = T), ]
    res_save$rank = 1:nrow(res_save)
    write.table(res_save, paste0(savepath, weight_type,"_",type, "_",method,"_res.txt"),sep='\t',
                quote = F, row.names = F)
    dat = res_save %>%
      arrange(desc(score)) %>%
      mutate("x" = row_number())
    res1 = head(dat, nrow(dat) * 0.1)  
    res1 = res1[res1$u>0 & res1$v < 0,]  
    range = ceiling(max(dat$score))

    p_rank <- ggplot(dat, aes(x = x, y = score))+
      geom_point(color="#1762f7", size = 0.5)+
      theme_bw()+
      labs(x = "Ranked drugs", y = "Potential")+
      geom_text_repel(inherit.aes = F, 
                      data = res1, 
                      aes(x = x, y = score, label = pert_iname, color = 'red'),
                      max.overlaps = getOption("ggrepel.max.overlaps", default = 200))+
      scale_y_continuous(breaks = seq(-range,range,2)) +
      ggtitle(paste0(cancer, "_", tumor_subtype, "_", dataset))
    ggsave(paste0(savepath, weight_type,"_",type, "_",method,"_res.pdf"),p_rank, width = 5, height = 12)
    
    # par(mfrow = c(1,1))
    # plot(dat$x, dat$score, pch=20,
    #      ylab = "Immune impact score")
    # abline(v = nrow(dat) * 0.01, lty = 6)
    # abline(v = nrow(dat) * 0.05, lty = 6)
    # abline(v = nrow(dat) * 0.1, lty = 6)
    # points(dcca$x, dcca$score, pch=21, col="red", bg="red")
    # text(dcca$x, dcca$score, paste0(dcca$drug_index,"\n", dcca$pert_iname),
    #      cex = 1, pos = 4, col = "red")
    # # range =max(abs(res$u),abs(res$v))
    # plot(res$u, res$v, pch=20,
    #      xlim=c(-range, range), ylim=c(-range,range),
    #      xlab ="sensitive score", ylab = "resistance score")
    # abline(v = 0, h = 0, lty = 6)
    # abline(coef = c(0,1), lty = 6)
    # abline(coef = c(0,-1), lty = 6)
    # points(res1$u, res1$v,pch=20,col="red")
    # text(res$u, res$v, gsub("drug","D",rownames(res)),
    #      cex = 0.7, pos = 4)
    # draw.circle(0,0,0.5, lty = 6)
    # draw.circle(0,0,1, lty = 6)
    # dev.off()
    
    
  }else{
    print("Error: ranker method not found! ")
  }
  
  return(res_save)
}




drugITSRanker_meta <- function(df_immunesig,
                               cancer,
                               weight_type,
                               tumor_subtype, 
                               method = c("glmnet_weight", "simple", "sum_mean"),
                               type = c("unweighted", "weighted", "binary_weighted"),
                              #  mergeITS,
                              #  mergeITS_type = c("Type", "Index"),
                               d1 = NULL,
                               d2 = NULL,
                               savepath){
  
  require(ggplot2)
  require(GGally)
  library(dplyr)
  library(plot3D)
  library(fmsb)
  library(ggrepel)

  library(fs)
  drugTargetMoA_file=path("Data", "drugTargetMoA_merged.txt")
  drugTargetMoA <- read.table(drugTargetMoA_file, sep = "\t", header = T)
  drugindexall_file=path("Data", "drug_index_LINCS.txt")
  drugindexall <- read.table(drugindexall_file, sep = "\t", header = T)
  
  # savepath = paste0(savepath, type,'/')
  savepath = paste0(savepath, "drugITSRank/")
  dir.create(savepath)
  df_immunesig$weight <- log2(df_immunesig$OR)


  if(type == "weighted"){

    df_immunesig$mNES_diff_weighted = df_immunesig$mNES_diff * df_immunesig$weight
    # df_s = df_immunesig[grep("01",df_immunesig$setindex),]
    # df_r = df_immunesig[grep("02",df_immunesig$setindex),]   
    # df_s = df_s[c("drug_index","immune_sig","mNES_diff_weighted")]
    # names(df_s) = c(c("drug_index","immune_sig_s","mNES_diff_s"))
    # df_r = df_r[c("drug_index","immune_sig","mNES_diff_weighted")]
    # names(df_r) = c(c("drug_index","immune_sig_r","mNES_diff_r"))
    
  }else if(type == "binary_weighted"){
    
    df_immunesig$mNES_diff_binary = 0
    df_immunesig[df_immunesig$mNES_diff > 0, ]$mNES_diff_binary = 1
    df_immunesig[df_immunesig$mNES_diff < 0, ]$mNES_diff_binary = -1
    
    df_immunesig$mNES_diff_weighted = df_immunesig$mNES_diff_binary * df_immunesig$weight
    
    # df_s = df_immunesig[grep("01",df_immunesig$setindex),]
    # df_r = df_immunesig[grep("02",df_immunesig$setindex),]
    # df_s = df_s[c("drug_index","immune_sig","mNES_diff_weighted")]
    # names(df_s) = c(c("drug_index","immune_sig_s","mNES_diff_s"))
    # df_r = df_r[c("drug_index","immune_sig","mNES_diff_weighted")]
    # names(df_r) = c(c("drug_index","immune_sig_r","mNES_diff_r"))
    
  }else if(type == "binary_unweighted"){
    
    df_immunesig$mNES_diff_binary = 0
    df_immunesig[df_immunesig$mNES_diff > 0, ]$mNES_diff_binary = 1
    df_immunesig[df_immunesig$mNES_diff < 0, ]$mNES_diff_binary = -1
    
    df_immunesig$mNES_diff_weighted = df_immunesig$mNES_diff_binary 
    
    # df_s = df_immunesig[grep("01",df_immunesig$setindex),]
    # df_r = df_immunesig[grep("02",df_immunesig$setindex),]    
    # df_s = df_s[c("drug_index","immune_sig","mNES_diff_weighted")]
    # names(df_s) = c(c("drug_index","immune_sig_s","mNES_diff_s"))
    # df_r = df_r[c("drug_index","immune_sig","mNES_diff_weighted")]
    # names(df_r) = c(c("drug_index","immune_sig_r","mNES_diff_r"))
    
  }else if(type == "unweighted"){
    
    df_immunesig$mNES_diff_weighted = df_immunesig$mNES_diff
    # df_s = df_immunesig[grep("01",df_immunesig$setindex),]
    # df_r = df_immunesig[grep("02",df_immunesig$setindex),]
    
    # df_s = df_s[c("drug_index","immune_sig","mNES_diff")]
    # names(df_s) = c(c("drug_index","immune_sig_s","mNES_diff_s"))
    # df_r = df_r[c("drug_index","immune_sig","mNES_diff")]
    # names(df_r) = c(c("drug_index","immune_sig_r","mNES_diff_r"))
    
  }
  
  
    
  # if(mergeITS){

  #   immuneSigInfoPath = "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/06_results_summaryandplotting/data/immene_cell_hieracy/"
    
  #   immunesigInfo <- as.data.frame(readxl::read_excel(paste0(immuneSigInfoPath, "immune_sig_hieracy_new.xlsx"), sheet = 1))
  #   immunesigInfo <- immunesigInfo[c(1,4,5:9)]
  #   immunesigInfo$immune_sig = tolower(immunesigInfo$immune_sig)
  #   immunesigInfo = immunesig_name_unify(immunesigInfo)

  #   df_s = df_immunesig[grep("01",df_immunesig$setindex),]
  #   df_s = df_s[c("drug_index","immune_sig","mNES_diff_weighted")]
  #   df_s = inner_join(df_s, immunesigInfo)
    
  #   df_s2 = lapply(as.list(unique(df_s$drug_index)), function(d){
  #     df1 = df_s[df_s$drug_index == d, ]
  #     tmp = aggregate(df1['mNES_diff_weighted'], by = list(df1[,mergeITS_type]), mean)
  #     tmp$drug_index = d 
  #     return(tmp)
  #   })
  #   df_s2 = do.call(rbind, df_s2) 
  #   df_s = df_s2[c("drug_index","Group.1","mNES_diff_weighted")]
  #   names(df_s) = c("drug_index","immune_sig_s","mNES_diff_s")


  #   df_r = df_immunesig[grep("02",df_immunesig$setindex),]
  #   df_r = df_r[c("drug_index","immune_sig","mNES_diff_weighted")]
  #   df_r = inner_join(df_r, immunesigInfo)
  #   df_r2 = lapply(as.list(unique(df_r$drug_index)), function(d){
  #     df1 = df_r[df_r$drug_index == d, ]
  #     tmp = aggregate(df1['mNES_diff_weighted'], by = list(df1[,mergeITS_type]), mean)
  #     tmp$drug_index = d 
  #     return(tmp)
  #   })
  #   df_r2 = do.call(rbind, df_r2) 
  #   df_r = df_r2[c("drug_index","Group.1","mNES_diff_weighted")]
  #   names(df_r) = c("drug_index","immune_sig_r","mNES_diff_r")

  # }else{

    df_s = df_immunesig[grep("01",df_immunesig$setindex),]
    df_r = df_immunesig[grep("02",df_immunesig$setindex),]
    
    df_s = df_s[c("drug_index","immune_sig","mNES_diff_weighted")]
    names(df_s) = c("drug_index","immune_sig_s","mNES_diff_s")
    df_r = df_r[c("drug_index","immune_sig","mNES_diff_weighted")]
    names(df_r) = c("drug_index","immune_sig_r","mNES_diff_r")

  # }

  
  df = merge(unique(df_s[c("drug_index","mNES_diff_s")]),
             unique(df_r[c("drug_index","mNES_diff_r")]),
             by= "drug_index")
  
  mat_s <-  reshape(df_s,  timevar = "immune_sig_s",
                    idvar = "drug_index",direction = "wide")
  rownames(mat_s) <- mat_s$drug_index
  mat_s <- mat_s[-1]
  
  mat_r <-  reshape(df_r,  timevar = "immune_sig_r",
                    idvar = "drug_index",direction = "wide")
  rownames(mat_r) <- mat_r$drug_index
  mat_r <- mat_r[-1]
  
  mat_s2 <- as.matrix(mat_s[is.element(rownames(mat_s), intersect(rownames(mat_s), rownames(mat_r) )),])
  mat_r2 <- as.matrix(mat_r[is.element(rownames(mat_r), intersect(rownames(mat_s), rownames(mat_r) )),])
  
  mat_s2[is.na(mat_s2)] = 0
  mat_r2[is.na(mat_r2)] = 0
  
  
  
  if(method == "simple"){
    library(plotrix)
    # dir.create(paste0(savepath,"/",method,"/"))
    
    CCAdata1_Z <- as.matrix(mat_s2)
    CCAdata2_Z <- as.matrix(mat_r2)
    # CCAdata1_Z <- as.matrix(scale(mat_s2, center = T))
    # CCAdata2_Z <- as.matrix(scale(mat_r2, center = T))
    
    
    if(length(which(is.na(CCAdata1_Z) == T))!=0)CCAdata1_Z[is.na(CCAdata1_Z)] <- 0
    if(length(which(is.na(CCAdata1_Z) == T))!=0)CCAdata2_Z[is.na(CCAdata2_Z)] <- 0
    
    u = as.data.frame(apply(CCAdata1_Z, 1, mean)); names(u) = "u"; u$drug_index = rownames(u)
    v = as.data.frame(apply(CCAdata2_Z, 1, mean)); names(v) = "v"; v$drug_index = rownames(v)
    
    
    res = merge(u, v, by = "drug_index")
    rownames(res) = res$drug_index
    res = res[c("u", "v")]
    res1 = res[res$u>0 & res$v <0,]  
    
    res$drug_index <- rownames(res)    
    res = res[c("drug_index","u","v")]
    
    res_save = inner_join(drugindexall, res)
    # res_save[res_save$u >0 & res_save$v <0 ,]
    res_save$score = (res_save$u + res_save$v) 
    res_save = res_save[order(res_save$score, decreasing = T), ]
    res_save$rank = 1:nrow(res_save)
    
    write.table(res_save, paste0(savepath, weight_type,"_",type, "_",method,"_res.txt"),sep='\t',
                quote = F, row.names = F)
    
    dat = res_save %>%
      arrange(desc(score)) %>%
      mutate("x" = row_number())
    res1 = head(dat, nrow(dat) * 0.1)  
    res1 = res1[res1$u>0 & res1$v < 0,]  
    range = ceiling(max(dat$score))

    p_rank <- ggplot(dat, aes(x = x, y = score))+
      geom_point(color="#1762f7", size = 0.5)+
      theme_bw()+
      labs(x = "Ranked drugs", y = "Potential")+
      geom_text_repel(inherit.aes = F, 
                      data = res1, 
                      aes(x = x, y = score, label = pert_iname, color = 'red'),
                      max.overlaps = getOption("ggrepel.max.overlaps", default = 200))+
      scale_y_continuous(breaks = seq(-range,range,2)) +
      ggtitle(paste0(cancer, "_", tumor_subtype, "_", dataset))
    ggsave(paste0(savepath, weight_type,"_",type, "_",method,"_res.pdf"),p_rank, width = 5, height = 12)
    # par(mfrow = c(1,1))
    # range =max(abs(res$u),abs(res$v))
    # plot(res$u, res$v, pch=20,
    #      xlim=c(-range, range), ylim=c(-range,range),
    #      xlab ="sensitive score", ylab = "resistance score")
    # abline(v = 0, h = 0, lty = 6)
    # abline(coef = c(0,1), lty = 6)
    # abline(coef = c(0,-1), lty = 6)
    # points(res1$u, res1$v,pch=20,col="red")
    # text(res$u, res$v, gsub("drug","D",rownames(res)),
    #      cex = 0.7, pos = 4)
    # draw.circle(0,0,0.5, lty = 6)
    # draw.circle(0,0,1, lty = 6)
    # dev.off()
    
  }else if(method == "sum_mean"){
    
    library(plotrix)
    # dir.create(paste0(savepath,"/",method,"/"))
    
    CCAdata1_Z <- as.matrix(mat_s2)
    CCAdata2_Z <- as.matrix(mat_r2)
    # CCAdata1_Z <- as.matrix(scale(mat_s2, center = T))
    # CCAdata2_Z <- as.matrix(scale(mat_r2, center = T))
    
    num_ITS = ncol(CCAdata1_Z) + ncol(CCAdata2_Z)
    if(length(which(is.na(CCAdata1_Z) == T))!=0)CCAdata1_Z[is.na(CCAdata1_Z)] <- 0
    if(length(which(is.na(CCAdata1_Z) == T))!=0)CCAdata2_Z[is.na(CCAdata2_Z)] <- 0
    
    u = as.data.frame(apply(CCAdata1_Z, 1, sum)); names(u) = "u"; u$drug_index = rownames(u)
    v = as.data.frame(apply(CCAdata2_Z, 1, sum)); names(v) = "v"; v$drug_index = rownames(v)
    
    
    res = merge(u, v, by = "drug_index")
    rownames(res) = res$drug_index
    res = res[c("u", "v")]
    res1 = res[res$u>0 & res$v <0,]  
    
    res$drug_index <- rownames(res)    
    res = res[c("drug_index","u","v")]
    
    res_save = inner_join(drugindexall, res)
    # res_save[res_save$u >0 & res_save$v <0 ,]
    res_save$score = (res_save$u + res_save$v) / num_ITS
    res_save = res_save[order(res_save$score, decreasing = T), ]
    res_save$rank = 1:nrow(res_save)
    write.table(res_save, paste0(savepath, weight_type,"_",type, "_",method,"_res.txt"),sep='\t',
                quote = F, row.names = F)
        
    dat = res_save %>%
      arrange(desc(score)) %>%
      mutate("x" = row_number())
    res1 = head(dat, nrow(dat) * 0.1)  
    res1 = res1[res1$u>0 & res1$v < 0,]  
    range = ceiling(max(dat$score))

    p_rank <- ggplot(dat, aes(x = x, y = score))+
      geom_point(color="#1762f7", size = 0.5)+
      theme_bw()+
      labs(x = "Ranked drugs", y = "Potential")+
      geom_text_repel(inherit.aes = F, 
                      data = res1, 
                      aes(x = x, y = score, label = pert_iname, color = 'red'),
                      max.overlaps = getOption("ggrepel.max.overlaps", default = 200))+
      scale_y_continuous(breaks = seq(-range,range,2)) +
      ggtitle(paste0(cancer, "_", tumor_subtype, "_", dataset))

      ggsave(paste0(savepath, weight_type,"_",type, "_",method,"_res.pdf"),p_rank, width = 5, height = 12)
    # ggsave(paste0(savepath, method,"/", "overall_drugranker_analysis.pdf"),p_rank, width = 5, height = 12)
    # par(mfrow = c(1,1))
    # range =max(abs(res$u),abs(res$v))
    # plot(res$u, res$v, pch=20,
    #      xlim=c(-range, range), ylim=c(-range,range),
    #      xlab ="sensitive score", ylab = "resistance score")
    # abline(v = 0, h = 0, lty = 6)
    # abline(coef = c(0,1), lty = 6)
    # abline(coef = c(0,-1), lty = 6)
    # points(res1$u, res1$v,pch=20,col="red")
    # text(res$u, res$v, gsub("drug","D",rownames(res)),
    #      cex = 0.7, pos = 4)
    # draw.circle(0,0,0.5, lty = 6)
    # draw.circle(0,0,1, lty = 6)
    # dev.off()
    
    
  }else{
    print("Error: ranker method not found! ")
  }
  
  return(res_save)
}



drugITSProfile <- function(df_immunesig,
                           cancer,
                           tumor_subtype, 
                           weight_type,
                           type,
                           method = c("glmnet_weight", "simple", "sum_mean"),
                           d1 = NULL,
                           d2 = NULL,
                           savepath){
  
  require(ggplot2)
  require(GGally)
  library(dplyr)
  library(plot3D)
  library(fmsb)
  library(ggrepel)
  library(fs)


  drugTargetMoA_file=path("Data", "drugTargetMoA_merged.txt")
  drugTargetMoA <- read.table(drugTargetMoA_file, sep = "\t", header = T)
  drugindexall_file=path("Data", "drug_index_LINCS.txt")
  drugindexall <- read.table(drugindexall_file, sep = "\t", header = T)
  
  df_s = df_immunesig[grep("01",df_immunesig$setindex),]
  df_r = df_immunesig[grep("02",df_immunesig$setindex),]
  
  df_s = df_s[c("drug_index","immune_sig","mNES_diff")]
  names(df_s) = c("drug_index","immune_sig_s","mNES_diff_s")
  df_r = df_r[c("drug_index","immune_sig","mNES_diff")]
  names(df_r) = c("drug_index","immune_sig_r","mNES_diff_r")
  
    
  # df_s = df_immunesig[grep("01",df_immunesig$setindex),]
  # df_r = df_immunesig[grep("02",df_immunesig$setindex),]
  
  # df_s = df_s[c("drug_index","immune_sig","mNES_diff")]
  # names(df_s) = c(c("drug_index","immune_sig_s","mNES_diff_s"))
  # df_r = df_r[c("drug_index","immune_sig","mNES_diff")]
  # names(df_r) = c(c("drug_index","immune_sig_r","mNES_diff_r"))
  
  
  df = merge(unique(df_s[c("drug_index","mNES_diff_s")]),
             unique(df_r[c("drug_index","mNES_diff_r")]),
             by= "drug_index")
  
  mat_s <-  reshape(df_s,  timevar = "immune_sig_s",
                    idvar = "drug_index",direction = "wide")
  rownames(mat_s) <- mat_s$drug_index
  mat_s <- mat_s[-1]
  
  mat_r <-  reshape(df_r,  timevar = "immune_sig_r",
                    idvar = "drug_index",direction = "wide")
  rownames(mat_r) <- mat_r$drug_index
  mat_r <- mat_r[-1]
  
  # mat_s2 <- as.matrix(mat_s[is.element(rownames(mat_s), intersect(rownames(mat_s), rownames(mat_r) )),])
  # mat_r2 <- as.matrix(mat_r[is.element(rownames(mat_r), intersect(rownames(mat_s), rownames(mat_r) )),])
  mat_s2 <- as.matrix(mat_s[intersect(rownames(mat_s), rownames(mat_r)),])
  rownames(mat_s2) = intersect(rownames(mat_s), rownames(mat_r))
  colnames(mat_s2) = names(mat_s)
  mat_r2 <- as.matrix(mat_r[intersect(rownames(mat_s), rownames(mat_r)),])
  rownames(mat_r2) = intersect(rownames(mat_s), rownames(mat_r))
  colnames(mat_r2) = names(mat_r)
  
  mat_s2[is.na(mat_s2)] = 0
  mat_r2[is.na(mat_r2)] = 0
  
  # IGeS_renaming    
  # iges_rename <- openxlsx::read.xlsx("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/IGeS_BS/IGeS_BS_tool/IGeS_Func_rename.xlsx")
  iges_rename_file=path("Data", "IGeS_Func_rename.xlsx")

  iges_rename <- openxlsx::read.xlsx(iges_rename_file)
  # iges_rename$Functional_summary = gsub(" and ", '&', iges_rename$Functional_summary)
  # iges_rename$Functional_summary = gsub("immune ", 'Imm.', iges_rename$Functional_summary)
  # iges_rename$Functional_summary = gsub("stimulatory", 'sti.', iges_rename$Functional_summary)
  # iges_rename$Functional_summary = gsub("regulatory", 'reg.', iges_rename$Functional_summary)
  # iges_rename$Functional_summary = gsub("antigen processing and presentation", 'antigen present.', iges_rename$Functional_summary)
  # iges_rename$Functional_summary = gsub("ICBresistant", 'ICB resistant', iges_rename$Functional_summary)
  iges_rename$Function_detail = gsub("Interleukin", 'IL', iges_rename$Function_detail)
  iges_rename$Function_detail = gsub("Cytokines", 'CY', iges_rename$Function_detail)
  iges_rename$Function_detail = gsub("Chemokine", 'CK', iges_rename$Function_detail)

  iges_rename$rename = paste0(iges_rename$IGeS.index, '\n', 
                                  iges_rename$Functional_summary, '\n', 
                                  iges_rename$Function_orignin, '\n',
                                  '(',iges_rename$Function_detail, ')')
  iges_rename$rename = gsub("\nNA", '', iges_rename$rename)
  # Reorder the column names 
  colnames(mat_s2) = gsub('mNES_diff_s.','',colnames(mat_s2))
  colnames(mat_r2) = gsub('mNES_diff_r.','',colnames(mat_r2))
  new_colname_s <- iges_rename$rename[match(colnames(mat_s2), iges_rename$IGeS)]
  colnames(mat_s2) <- new_colname_s
  desired_order <- intersect(iges_rename$rename, colnames(mat_s2))
  mat_s2 <- mat_s2[, desired_order]

  new_colname_r <- iges_rename$rename[match(colnames(mat_r2), iges_rename$IGeS)]
  colnames(mat_r2) <- new_colname_r
  desired_order <- intersect(iges_rename$rename, colnames(mat_r2))
  mat_r2 <- mat_r2[, desired_order]

  res_save <- read.csv(paste0(savepath, "drugITSRank/", weight_type,"_",type, "_",method,"_res.txt"),sep='\t', header = T)
  dat = res_save %>%
    arrange(desc(score)) %>%
    mutate("x" = row_number())
  
  dir.create(paste0(savepath, 'drugITSProfile/'))
  # dir.create(paste0(savepath, 'drugITSProfile/',type))
  savepath = paste0(savepath, 'drugITSProfile/')

  
  if(method == "glmnet_weight"){
    library(plotrix)
    savepath = paste0(savepath, weight_type,"_", type, '_', method,"/")
    dir.create(savepath)

    CCAdata1_Z <- as.matrix(mat_s2)
    CCAdata2_Z <- as.matrix(mat_r2)

    # CCAdata1_Z <- as.matrix(scale(mat_s2, center = T))
    # CCAdata2_Z <- as.matrix(scale(mat_r2, center = T))
    
    
    if(length(which(is.na(CCAdata1_Z) == T))!=0)CCAdata1_Z[is.na(CCAdata1_Z)] <- 0
    if(length(which(is.na(CCAdata1_Z) == T))!=0)CCAdata2_Z[is.na(CCAdata2_Z)] <- 0
    
    u = as.data.frame(apply(CCAdata1_Z, 1, sum)); names(u) = "u"; u$drug_index = rownames(u)
    v = as.data.frame(apply(CCAdata2_Z, 1, sum)); names(v) = "v"; v$drug_index = rownames(v)
    
    
    # res = merge(u, v, by = "drug_index")
    # rownames(res) = res$drug_index
    # res = res[c("u", "v")]
    # res1 = res[res$u>0 & res$v <0,]  
    
    # res$drug_index <- rownames(res)    
    # res = res[c("drug_index","u","v")]
    
    # res_save = inner_join(drugindexall, res)
    # # res_save[res_save$u >0 & res_save$v <0 ,]
    # res_save$score = res_save$u + res_save$v
    # res_save = res_save[order(res_save$score, decreasing = T), ]
    # res_save$rank = 1:nrow(res_save)
    
    
    for(drug in sort(unique(unique(res_save$drug_index)))){
      # drug='drug100'
      dtm <- drugTargetMoA[drugTargetMoA$drug_index == drug, ]
      if(is.element("", dtm$MoA))dtm <- dtm[-which(dtm$MoA == ""), ]
      if(nrow(dtm) < 1 ){
        
        dtm <- drugindexall[drugindexall$drug_index == drug, ]
        dtm$MoA = ""
        dtm$target = ""
        
      }
      drugname <- unique(dtm$pert_iname)
      
      dcca <- dat[dat$drug_index == drug, ]
      ds <- CCAdata1_Z[rownames(CCAdata1_Z) == drug, ]
      dr <- CCAdata2_Z[rownames(CCAdata2_Z) == drug, ]


      drugname = gsub("\\(", "", drugname)
      drugname = gsub("\\)", "", drugname)
      drugname = gsub("\\/", ".", drugname)
      
      
      dat_shown <- dat[dat$drug_index == drug, ]
      

      # pdf(paste0(savepath, method,"/",drug,"_",drugname,".pdf"), width = 15, height = 8)
      pdf(paste0(savepath, drug,"_",drugname,".pdf"), width = 15, height = 8)
      
      # par(mfrow = c(1,2), pty = "s")
      mat <- matrix(c(1,1,2,2,2,2), nrow=1, byrow=TRUE)
      # mat
      layout(mat)
      
      plot(dat$x, dat$score, pch=20,
           ylab = "Immune impact score")
      abline(v = nrow(dat) * 0.01, lty = 6)
      abline(v = nrow(dat) * 0.05, lty = 6)
      abline(v = nrow(dat) * 0.1, lty = 6)
      points(dcca$x, dcca$score, pch=21, col="red", bg="red")
      text(dcca$x, dcca$score, paste0(dcca$drug_index,"\n", dcca$pert_iname),
           cex = 2, pos = 4, col = "red")
      
      # range =max(abs(res$u),abs(res$v))
      # plot(res$u, res$v, pch=20,
      #      xlim=c(-range, range), ylim=c(-range,range),
      #      xlab ="sensitive score", ylab = "resistance score")
      # abline(v = 0, h = 0, lty = 6)
      # abline(coef = c(0,1), lty = 6)
      # abline(coef = c(0,-1), lty = 6)
      # points(dcca$u, dcca$v, pch=21, col="red", bg="red")
      # arrows(x0 = 0, y0 = 0, x1 = dcca$u, y1 = dcca$v, lwd = 1.5, col = "red", length = 0.1)
      # text(dcca$u, dcca$v, gsub("drug", "D", rownames(dcca)),
      #      cex = 1, pos = 4, col = "red")
      # draw.circle(0,0,0.5, lty = 6)
      # draw.circle(0,0,1, lty = 6)
      
      ds1 = (matrix(0,length(dr))); rownames(ds1) = names(dr)
      dr1 = (matrix(0,length(ds))); rownames(dr1) = names(ds)
      ds2 = rbind(as.matrix(ds),ds1);colnames(ds2) = "s"
      dr2 = rbind(dr1,as.matrix(dr));colnames(dr2) = "r"
      
      df <- as.data.frame(cbind(ds2, dr2))
      df$max <- ceiling(max(c(max(mat_s2), max(mat_r2))))
      df$min <- floor(min(c(min(mat_s2), min(mat_r2))))
      df$mean <- 0
      df <- df[c("max","min","r","s","mean")]
      # df$V1 <-  (df$V1 - mean(df$V1))/sd(df$V1)
      df <- as.data.frame(t(df))
      colnames(df) <- gsub("mNES_diff_s.", "", colnames(df))
      colnames(df) <- gsub("mNES_diff_r.", "", colnames(df))
      p2 = radarchart(df,
                      axistype = 4,
                      # pfcol = c(rgb(0.2,0.2,0.2,0.9,0.9),NA,NA),
                      pfcol = c(NA,NA, NA),
                      pcol= c(rgb(0.2,0.5,0.5,0.5,0.9),
                              rgb(0.7,0.2,0.2,0.5,0.9),
                              rgb(0.22,0.22,0.22,0.5)),
                      plty = 1, 
                      plwd = 4,
                      # Customize the grid
                      cglcol = "grey70",  
                      cglwd = 0.8,
                      # Customize the axis
                      axislabcol = "black", 
                      # Variable labels
                      cglty = 1,
                      # caxislabels = seq(floor(min(c(min(mat_s3), min(mat_r3)))),
                      #                   ceiling(max(c(max(mat_s3), max(mat_r3)))),
                      #                   3),
                      caxislabels =  c(seq(-4,4,2)),
                      vlcex = 0.6,
                      title = paste0(# "Impact on Immune\n",
                                     # "Drug: ", drug,"_",drugname, "\n",
                                     "Rank: ",dcca$rank,"/", nrow(dat),"\n",
                                     "MoA: ", unique(dtm$MoA), "\n",
                                     "Target: ", paste(unique(dtm$target), collapse="/")))
      # abline(h =0, v = 0, lty = 2)
      dev.off()
      
    }
    
    
  }else if(method == "simple"){
    library(plotrix)
    # dir.create(paste0(savepath,"/",method,"/"))
    savepath = paste0(savepath, weight_type,"_", type, '_', method,"/")
    dir.create(savepath)

    
    CCAdata1_Z <- as.matrix(mat_s2)
    CCAdata2_Z <- as.matrix(mat_r2)
    # CCAdata1_Z <- as.matrix(scale(mat_s2, center = T))
    # CCAdata2_Z <- as.matrix(scale(mat_r2, center = T))
    
    
    if(length(which(is.na(CCAdata1_Z) == T))!=0)CCAdata1_Z[is.na(CCAdata1_Z)] <- 0
    if(length(which(is.na(CCAdata1_Z) == T))!=0)CCAdata2_Z[is.na(CCAdata2_Z)] <- 0
    
    # u = as.data.frame(apply(CCAdata1_Z, 1, mean)); names(u) = "u"; u$drug_index = rownames(u)
    # v = as.data.frame(apply(CCAdata2_Z, 1, mean)); names(v) = "v"; v$drug_index = rownames(v)
    
    
    # res = merge(u, v, by = "drug_index")
    # rownames(res) = res$drug_index
    # res = res[c("u", "v")]
    # res1 = res[res$u>0 & res$v <0,]  
    
    # res$drug_index <- rownames(res)    
    # res = res[c("drug_index","u","v")]
    
    # res_save = inner_join(drugindexall, res)
    # # res_save[res_save$u >0 & res_save$v <0 ,]
    # res_save$score = (res_save$u - res_save$v) 
    # res_save = res_save[order(res_save$score, decreasing = T), ]
    # res_save$rank = 1:nrow(res_save)
    
    for(drug in sort(unique(unique(res_save$drug_index)))){
      
      dtm <- drugTargetMoA[drugTargetMoA$drug_index == drug, ]
      if(is.element("", dtm$MoA))dtm <- dtm[-which(dtm$MoA == ""), ]
      if(nrow(dtm) < 1 ){
        
        dtm <- drugindexall[drugindexall$drug_index == drug, ]
        dtm$MoA = ""
        dtm$target = ""
        
      }
      drugname <- unique(dtm$pert_iname)
      
      dcca <- dat[dat$drug_index == drug, ]
      ds <- CCAdata1_Z[rownames(CCAdata1_Z) == drug, ]
      dr <- CCAdata2_Z[rownames(CCAdata2_Z) == drug, ]
      
      drugname = gsub("\\(", "", drugname)
      drugname = gsub("\\)", "", drugname)
      drugname = gsub("\\/", ".", drugname)
      
      
      
      # pdf(paste0(savepath, method,"/",drug,"_",drugname,".pdf"), width = 15, height = 8)
      pdf(paste0(savepath, drug,"_",drugname,".pdf"), width = 15, height = 8)


      # par(mfrow = c(1,2), pty = "s")
      mat <- matrix(c(1,1,2,2,2,2), nrow=1, byrow=TRUE)
      # mat
      layout(mat)
      
      plot(dat$x, dat$score, pch=20,
           ylab = "Immune impact score")
      abline(v = nrow(dat) * 0.01, lty = 6)
      abline(v = nrow(dat) * 0.05, lty = 6)
      abline(v = nrow(dat) * 0.1, lty = 6)
      points(dcca$x, dcca$score, pch=21, col="red", bg="red")
      text(dcca$x, dcca$score, paste0(dcca$drug_index,"\n", dcca$pert_iname),
           cex = 1, pos = 4, col = "red")
      
      # range =max(abs(res$u),abs(res$v))
      # plot(res$u, res$v, pch=20,
      #      xlim=c(-range, range), ylim=c(-range,range),
      #      xlab ="sensitive score", ylab = "resistance score")
      # abline(v = 0, h = 0, lty = 6)
      # abline(coef = c(0,1), lty = 6)
      # abline(coef = c(0,-1), lty = 6)
      # points(dcca$u, dcca$v, pch=21, col="red", bg="red")
      # arrows(x0 = 0, y0 = 0, x1 = dcca$u, y1 = dcca$v, lwd = 1.5, col = "red", length = 0.1)
      # text(dcca$u, dcca$v, gsub("drug", "D", rownames(dcca)),
      #      cex = 1, pos = 4, col = "red")
      # draw.circle(0,0,0.5, lty = 6)
      # draw.circle(0,0,1, lty = 6)
      
      ds1 = (matrix(0,length(dr))); rownames(ds1) = names(dr)
      dr1 = (matrix(0,length(ds))); rownames(dr1) = names(ds)
      ds2 = rbind(as.matrix(ds),ds1);colnames(ds2) = "s"
      dr2 = rbind(dr1,as.matrix(dr));colnames(dr2) = "r"
      
      df <- as.data.frame(cbind(ds2, dr2))
      df$max <- ceiling(max(c(max(mat_s2), max(mat_r2))))
      df$min <- floor(min(c(min(mat_s2), min(mat_r2))))
      df$mean <- 0
      df <- df[c("max","min","r","s","mean")]
      # df$V1 <-  (df$V1 - mean(df$V1))/sd(df$V1)
      df <- as.data.frame(t(df))
      colnames(df) <- gsub("mNES_diff_s.", "", colnames(df))
      colnames(df) <- gsub("mNES_diff_r.", "", colnames(df))
      p2 = radarchart(df,
                      axistype = 4,
                      # pfcol = c(rgb(0.2,0.2,0.2,0.9,0.9),NA,NA),
                      pfcol = c(NA,NA, NA),
                      pcol= c(rgb(0.2,0.5,0.5,0.5,0.9),
                              rgb(0.7,0.2,0.2,0.5,0.9),
                              rgb(0.22,0.22,0.22,0.5)),
                      plty = 1, 
                      plwd = 4,
                      # Customize the grid
                      cglcol = "grey70",  cglwd = 0.8,
                      # Customize the axis
                      axislabcol = "black", 
                      # Variable labels
                      cglty = 1,
                      # caxislabels = seq(floor(min(c(min(mat_s3), min(mat_r3)))),
                      #                   ceiling(max(c(max(mat_s3), max(mat_r3)))),
                      #                   3),
                      caxislabels =  c(seq(-4,4,2)),
                      vlcex = 0.6,
                      title = paste0(# "Impact on Immune\n",
                                     # "Drug: ", drug,"_",drugname, "\n",
                                     "Rank: ",dcca$rank,"/", nrow(dat),"\n",
                                     "MoA: ", unique(dtm$MoA), "\n",
                                     "Target: ", paste(unique(dtm$target), collapse="/")))
      # abline(h =0, v = 0, lty = 2)
      dev.off()

    }
    
    
  }else if(method == "sum_mean"){
    
    library(plotrix)
    # dir.create(paste0(savepath,"/",method,"/"))
    savepath = paste0(savepath, weight_type,"_", type, '_', method,"/")
    dir.create(savepath)
    
    CCAdata1_Z <- as.matrix(mat_s2)
    CCAdata2_Z <- as.matrix(mat_r2)
    # CCAdata1_Z <- as.matrix(scale(mat_s2, center = T))
    # CCAdata2_Z <- as.matrix(scale(mat_r2, center = T))
    
    num_ITS = ncol(CCAdata1_Z) + ncol(CCAdata2_Z)
    if(length(which(is.na(CCAdata1_Z) == T))!=0)CCAdata1_Z[is.na(CCAdata1_Z)] <- 0
    if(length(which(is.na(CCAdata1_Z) == T))!=0)CCAdata2_Z[is.na(CCAdata2_Z)] <- 0
    
    # u = as.data.frame(apply(CCAdata1_Z, 1, sum)); names(u) = "u"; u$drug_index = rownames(u)
    # v = as.data.frame(apply(CCAdata2_Z, 1, sum)); names(v) = "v"; v$drug_index = rownames(v)
    
    
    # res = merge(u, v, by = "drug_index")
    # rownames(res) = res$drug_index
    # res = res[c("u", "v")]
    # res1 = res[res$u>0 & res$v <0,]  
    
    # res$drug_index <- rownames(res)    
    # res = res[c("drug_index","u","v")]
    
    # res_save = inner_join(drugindexall, res)
    # # res_save[res_save$u >0 & res_save$v <0 ,]
    # res_save$score = (res_save$u - res_save$v) / num_ITS
    # res_save = res_save[order(res_save$score, decreasing = T), ]
    # res_save$rank = 1:nrow(res_save)
    
    for(drug in sort(unique(unique(res_save$drug_index)))){
      
      dtm <- drugTargetMoA[drugTargetMoA$drug_index == drug, ]
      if(is.element("", dtm$MoA))dtm <- dtm[-which(dtm$MoA == ""), ]
      if(nrow(dtm) < 1 ){
        
        dtm <- drugindexall[drugindexall$drug_index == drug, ]
        dtm$MoA = ""
        dtm$target = ""
        
      }
      drugname <- unique(dtm$pert_iname)
      
      dcca <- dat[dat$drug_index == drug, ]
      ds <- CCAdata1_Z[rownames(CCAdata1_Z) == drug, ]
      dr <- CCAdata2_Z[rownames(CCAdata2_Z) == drug, ]
      
      drugname = gsub("\\(", "", drugname)
      drugname = gsub("\\)", "", drugname)
      drugname = gsub("\\/", ".", drugname)
      
      
      # pdf(paste0(savepath, method,"/",drug,"_",drugname,".pdf"), width = 15, height = 8)
      pdf(paste0(savepath, drug,"_",drugname,".pdf"), width = 15, height = 8)

      
      # par(mfrow = c(1,2), pty = "s")
      mat <- matrix(c(1,1,2,2,2,2), nrow=1, byrow=TRUE)
      # mat
      layout(mat)
      
      plot(dat$x, dat$score, pch=20,
           ylab = "Immune impact score")
      abline(v = nrow(dat) * 0.01, lty = 6)
      abline(v = nrow(dat) * 0.05, lty = 6)
      abline(v = nrow(dat) * 0.1, lty = 6)
      points(dcca$x, dcca$score, pch=21, col="red", bg="red")
      text(dcca$x, dcca$score, paste0(dcca$drug_index,"\n", dcca$pert_iname),
           cex = 1, pos = 4, col = "red")
      
      # range =max(abs(res$u),abs(res$v))
      # plot(res$u, res$v, pch=20,
      #      xlim=c(-range, range), ylim=c(-range,range),
      #      xlab ="sensitive score", ylab = "resistance score")
      # abline(v = 0, h = 0, lty = 6)
      # abline(coef = c(0,1), lty = 6)
      # abline(coef = c(0,-1), lty = 6)
      # points(dcca$u, dcca$v, pch=21, col="red", bg="red")
      # arrows(x0 = 0, y0 = 0, x1 = dcca$u, y1 = dcca$v, lwd = 1.5, col = "red", length = 0.1)
      # text(dcca$u, dcca$v, gsub("drug", "D", rownames(dcca)),
      #      cex = 1, pos = 4, col = "red")
      # draw.circle(0,0,0.5, lty = 6)
      # draw.circle(0,0,1, lty = 6)
      
      ds1 = (matrix(0,length(dr))); rownames(ds1) = names(dr)
      dr1 = (matrix(0,length(ds))); rownames(dr1) = names(ds)
      ds2 = rbind(as.matrix(ds),ds1);colnames(ds2) = "s"
      dr2 = rbind(dr1,as.matrix(dr));colnames(dr2) = "r"
      
      df <- as.data.frame(cbind(ds2, dr2))
      df$max <- ceiling(max(c(max(mat_s2), max(mat_r2))))
      df$min <- floor(min(c(min(mat_s2), min(mat_r2))))
      df$mean <- 0
      df <- df[c("max","min","r","s","mean")]
      # df$V1 <-  (df$V1 - mean(df$V1))/sd(df$V1)
      df <- as.data.frame(t(df))
      colnames(df) <- gsub("mNES_diff_s.", "", colnames(df))
      colnames(df) <- gsub("mNES_diff_r.", "", colnames(df))
      p2 = radarchart(df,
                      axistype = 4,
                      # pfcol = c(rgb(0.2,0.2,0.2,0.9,0.9),NA,NA),
                      pfcol = c(NA,NA, NA),
                      pcol= c(rgb(0.2,0.5,0.5,0.5,0.9),
                              rgb(0.7,0.2,0.2,0.5,0.9),
                              rgb(0.22,0.22,0.22,0.5)),
                      plty = 1, 
                      plwd = 4,
                      # Customize the grid
                      cglcol = "grey70",  cglwd = 0.8,
                      # Customize the axis
                      axislabcol = "black", 
                      # Variable labels
                      cglty = 1,
                      # caxislabels = seq(floor(min(c(min(mat_s3), min(mat_r3)))),
                      #                   ceiling(max(c(max(mat_s3), max(mat_r3)))),
                      #                   3),
                      caxislabels =  c(seq(-4,4,2)),
                      vlcex = 0.6,
                      title = paste0(# "Impact on Immune\n",
                                     # "Drug: ", drug,"_",drugname, "\n",
                                     "Rank: ",dcca$rank,"/", nrow(dat),"\n",
                                     "MoA: ", unique(dtm$MoA), "\n",
                                     "Target: ", paste(unique(dtm$target), collapse="/")))
      # abline(h =0, v = 0, lty = 2)
      dev.off()

    }
    
    
  }else {
    print("Error: ranker method not found! ")
  }
  
  return(res_save)
}



drugITSscore <- function(cancer = "LIHC",
                         tumor_subtype = "Primary",
                         purity_method = "TUMERIC",
                         dataset=70138,
                         datatype="allgenes",
                         num_gene = 200,
                         immunesig, 
                         model = c("m1", "m2", "m3"),
                         ACAT_Pval = 0.05,
                        #  work_path = work_path,
                         filepath = filepath,
                         savepath){
  
  # savepath = paste0("06_results_summaryandplotting/results_", purity_method,"_merged/")
  dir.create(savepath)
  
  # immunesig <- immunesig_extraction(cancer = cancer,
  #                                   tumor_subtype = tumor_subtype)
  
  positive_df = df_for_plot(cancer = cancer, 
                            tumor_subtype = tumor_subtype,
                            purity_method = purity_method,
                            dataset=dataset, 
                            sign = "positive", 
                            num_gene = 200, 
                            datatype = datatype, 
                            immunesig = immunesig, 
                            ACAT_Pval = ACAT_Pval,
                            filepath = filepath, 
                            # work_path = work_path,
                            savepath = savepath)

  # negative_df = df_for_plot(cancer = cancer, 
  #                           tumor_subtype = tumor_subtype,
  #                           purity_method = purity_method,
  #                           dataset=dataset, 
  #                           sign = "negative", 
  #                           num_gene = 200,  
  #                           datatype = datatype,  
  #                           immunesig = immunesig,
  #                           ACAT_Pval = ACAT_Pval,
  #                           filepath = filepath, 
  #                           # work_path = work_path,
  #                           savepath = savepath)
  
  
  names(positive_df) = c("drug_index","immune_sig", "mNES_p","mPvalueNaive_p",
                         "mPvalue_fisher_p","mPvalue_ACAT_p")
  # names(negative_df) = c("drug_index","immune_sig", "mNES_n","mPvalueNaive_n",
  #                        "mPvalue_fisher_n","mPvalue_ACAT_n")
  # df = merge(positive_df, negative_df, by = c("drug_index", "immune_sig"), all=T)
  # df[is.na(df$mNES_p),]$mNES_p=0
  # df[is.na(df$mNES_n),]$mNES_n=0
  df=positive_df

  if(model == "m1"){
    
    df$mNES_diff = df$mNES_p - df$mNES_n
    # df = df[which(df$mPvalue_ACAT_p < 0.05 ),]
    # df = df[which(df$mPvalue_ACAT_p < 0.05 & df$mPvalue_ACAT_n < 0.05),]
    # df = df[-grep("oe_",df$immune_sig),]
    
    df_immunesig = merge(df, immunesig, by = "immune_sig")
    if(length(unique(df_immunesig$immune_sig)) < 0.1*length(unique(immunesig$immune_sig))){
      print("analysis only with positive gene set")
      df_immunesig = merge(positive_df, immunesig, by = "immune_sig")
      df_immunesig$mNES_diff = df_immunesig$mNES_p - 0
    }
  }else if(model == "m2"){
    
    df$mNES_diff = df$mNES_p
    df_immunesig = merge(df, immunesig, by = "immune_sig")
    
  }else if(model == "m3"){
    
    df$mNES_diff = df$mNES_n
    df_immunesig = merge(df, immunesig, by = "immune_sig")
    
  }
  
  return(df_immunesig)
}



drugITSscore_plot <- function(cancer = "LIHC",
                              tumor_subtype = "Primary",
                              purity_method = "TUMERIC",
                              immunesig, 
                              dataset="70138",
                              datatype="allgenes",
                              num_gene = 200,
                              model = "m2", # c("m1", "m2"),
                              weight_type = "model_weight", # c("model_weight", "meta_weight"),
                              type = "weighted", #c("unweighted", "weighted", "binary_weighted", "binary_unweighted"),
                              drugrank_method = "glmnet_weight", #c("glmnet_weight", "simple", "sum_mean"),
                              ACAT_Pval = 0.05,
                              ifProfile = FALSE, 
                              filepath, 
                              savepath){
  
  dir.create(paste0(savepath, "/"))
  dir.create(paste0(savepath, "/",dataset,"/"))
  dir.create(paste0(savepath, "/",dataset,"/", datatype,"/"))
  dir.create(paste0(savepath, "/",dataset,"/", datatype,"/", cancer, "_", tumor_subtype))
  
  
  df_immunesig <- drugITSscore(cancer = cancer,
                               tumor_subtype = tumor_subtype,
                               purity_method = purity_method,
                               dataset=dataset,
                               datatype=datatype,
                               num_gene = num_gene,
                               immunesig = immunesig, 
                               model = model,
                               ACAT_Pval = ACAT_Pval,
                               # work_path = work_path,
                               filepath = filepath,
                               savepath = savepath)
  
  
  df_immunesig$log10_p_p <- -log10(df_immunesig$mPvalue_ACAT_p)
  # df_immunesig = df_immunesig[-grep("oe_",df_immunesig$immune_sig),]
  immunesig = immunesig[order(immunesig$immune_sig),]
  immunesig = immunesig[order(immunesig$setindex),]
  # df_immunesig$immune_sig <- factor(df_immunesig$immune_sig, levels = unique(immunesig$immune_sig))
  node_color <- colorRampPalette(colors = c("steelblue", "white", "#DC143C"))
  
  p1 <- ggplot(df_immunesig, aes(x = immune_sig, y = drug_index)) + 
    geom_point(aes(size = log10_p_p, color = mNES_diff , )) +
    scale_size(range = c(2, 8)) +
    scale_colour_gradientn(colours = node_color(500), na.value = "grey90", values = seq(0, 1, length.out = 500)) +
    scale_size_area(max_size = 5) 
  p1 <- p1 + theme_bw() +
    theme(axis.text.x  =  element_text(angle=-90, vjust=0.5,hjust = 0, size=8),
          axis.text.y = element_text(size = 15),
          axis.title.x = element_blank()) # +
  # ggtitle("Plot1a: gene positively correlated with immune signatures \n 
  #           immune signatures increase risk")
  # p1
  
  plot_savepath <- paste0(savepath, "/",dataset,"/", datatype,"/",cancer,"_", tumor_subtype, "/p_minus_n.pdf")
  ggsave(plot_savepath, p1, width = 20, height = 45)
  
  
  
  
  # drug ITS analysis
  # savepath = paste0("06_results_summaryandplotting/results_", purity_method,"/")
  # dir.create(savepath)
  # dir.create(paste0(savepath, "/", dataset, "/"))
  # dir.create(paste0(savepath, "/", dataset, "/", datatype, "/"))
  # dir.create(paste0(savepath, "/", dataset, "/", datatype, "/",cancer, "_", tumor_subtype, "/"))
  savepath = paste0(savepath, "/", dataset, "/", datatype, "/",cancer, "_", tumor_subtype, "/")
  # type='weighted'
  
  
  
  if(weight_type == "model_weight"){
    res_save = drugITSRanker(df_immunesig = df_immunesig,
                             weight_type=weight_type,
                             method = drugrank_method, #  options: "simple","sum_mean", "base" (cca), "scca"
                             type = type,
                            #  mergeITS = mergeITS,
                            #  mergeITS_type = mergeITS_type,
                             cancer = cancer,
                             tumor_subtype = tumor_subtype,
                             d1 = 1,
                             d2 = 2,
                             savepath = savepath)
    
  }else if(weight_type == "meta_weight"){
    
    res_save = drugITSRanker_meta(df_immunesig = df_immunesig,
                                  weight_type=weight_type,
                                  method = drugrank_method, #  options: "simple","sum_mean", "base" (cca), "scca"
                                  # mergeITS = mergeITS,
                                  # mergeITS_type = mergeITS_type,
                                  type = type,
                                  cancer = cancer,
                                  tumor_subtype = tumor_subtype,
                                  d1 = 1,
                                  d2 = 2,
                                  savepath = savepath)
    
  }
  
  if(ifProfile){
    
    drugprofile = drugITSProfile(df_immunesig = df_immunesig,
                                 weight_type=weight_type,
                                 method = drugrank_method, #  options: "simple","sum_mean", "base" (cca), "scca"
                                 type = type,
                                 cancer = cancer,
                                 tumor_subtype = tumor_subtype,
                                 d1 = 1,
                                 d2 = 2,
                                 savepath = savepath)
    
  }
}

# drug_MF(df_immunesig = res_save,
#         res_save = res_save,
#         matrix = "mm", # mm: matrix multiplication
#         savepath = savepath)

# drug_distri_distance(df_immunesig = df_immunesig,
#                      res_save = res_save,
#                      savepath = savepath)


# library(RColorBrewer)
# display.brewer.all()
# head(brewer.pal.info)


# res <- read.csv(paste0(savepath,"/ks_analysis_res.txt"), header = T, sep = "\t")
# # strategy 1
# res <- res[order(res$u_minus_v, decreasing = T), ]
# res$color_1 <- colorRampPalette(c("red", "white","blue"))(nrow(res))

# # strategy 2
# res$weight_1 <- 1
# res[res$ks_p_s_vs_srand >= 0.05 & res$ks_p_r_vs_rrand >= 0.05, ]$weight_1 = -1
# res$score_2 <- res$ks_D_s_vs_r * res$weight_1
# res <- res[order(res$score_2, decreasing = T), ]
# res$color_2 <- colorRampPalette(c("red", "white","blue"))(nrow(res))


# # strategy 3
# res$score_3 <- res$wass_dist_s_vs_r * res$weight_1
# res <- res[order(res$score_3, decreasing = T), ]
# res$color_3 <- colorRampPalette(c("red", "white","blue"))(nrow(res))

# # strategy 2
# res$weight_2_s <- 1
# res$weight_2_r <- -1
# res[res$ks_p_s_vs_srand < 0.05 & res$ks_p_r_vs_rrand < 0.05, ]$weight_2_s = 1
# res[res$ks_p_s_vs_srand < 0.05 & res$ks_p_r_vs_rrand < 0.05, ]$weight_2_r = 1
# res[res$ks_p_s_vs_srand < 0.05 & res$ks_p_r_vs_rrand >= 0.05, ]$weight_2_s = 1
# res[res$ks_p_s_vs_srand < 0.05 & res$ks_p_r_vs_rrand >= 0.05, ]$weight_2_r = -1
# res[res$ks_p_s_vs_srand >= 0.05 & res$ks_p_r_vs_rrand < 0.05, ]$weight_2_s = -1
# res[res$ks_p_s_vs_srand >= 0.05 & res$ks_p_r_vs_rrand < 0.05, ]$weight_2_r = 1
# res[res$ks_p_s_vs_srand >= 0.05 & res$ks_p_r_vs_rrand >= 0.05, ]$weight_2_s = -1
# res[res$ks_p_s_vs_srand >= 0.05 & res$ks_p_r_vs_rrand >= 0.05, ]$weight_2_r = -1

# res$score_4 <- res$ks_D_s_vs_srand * res$weight_2_s + res$ks_D_r_vs_rrand * res$weight_2_r
# res <- res[order(res$score_4, decreasing = T), ]
# res$color_4 <- colorRampPalette(c("red", "white","blue"))(nrow(res))

# res$score_5 <- res$wass_dist_s_vs_srand * res$weight_2_s + res$wass_dist_r_vs_rrand * res$weight_2_r
# res <- res[order(res$score_5, decreasing = T), ]
# res$color_5 <- colorRampPalette(c("red", "white","blue"))(nrow(res))

# write.table(res, paste0(savepath,"/ks_analysis_res_scoring.txt"),
#             sep = "\t",  quote = F, row.names = F)

# res1 = res[res$ks_p_s_vs_r < 0.05, ]
# res2 = res[res$ks_p_s_vs_srand < 0.05 & res$ks_p_r_vs_rrand < 0.05, ]
# res3 = res[res$wass_dist_pval_s_vs_r < 0.05, ]
# res4 = res[res$wass_dist_pval_s_vs_srand < 0.05 & res$wass_dist_pval_r_vs_rrand < 0.05, ]

# pdf(paste0(savepath, "overall_drug_ITS_analysis_ks_p_s_vs_r.pdf"), width = 8, height = 12)
#   par(mfrow = c(3,2))
#   range =max(abs(res$u),abs(res$v))
#   plot(res$u, res$v, pch=20,
#        # col=res$color_1,
#        xlim=c(-range, range), ylim=c(-range,range),
#        xlab ="sensitive score", ylab = "resistance score")
#   abline(v = 0, h = 0, lty = 6)
#   abline(coef = c(0,1), lty = 6)
#   abline(coef = c(0,-1), lty = 6)
#   # points(res1$u, res1$v,pch=20,col="red")
#   # text(res$u, res$v, gsub("drug","D",rownames(res)),
#   #      cex = 0.7, pos = 4)
#   # draw.circle(0,0,0.5, lty = 6)
#   # draw.circle(0,0,1, lty = 6)

#   plot(res$u, res$v, pch=20,
#        col=res$color_1,
#        xlim=c(-range, range), ylim=c(-range,range),
#        xlab ="sensitive score", ylab = "resistance score",
#        main = "mean_s_minus_mean_r")
#   abline(v = 0, h = 0, lty = 6)
#   abline(coef = c(0,1), lty = 6)
#   abline(coef = c(0,-1), lty = 6)
#   for(i in 1:nrow(res1)){
#     points(res1[i,]$u, res1[i,]$v,pch=1, col="black",lwd = 0.2)
#   }

#   plot(res$u, res$v, pch=20,
#        col=res$color_2,
#        xlim=c(-range, range), ylim=c(-range,range),
#        xlab ="sensitive score", ylab = "resistance score",
#        main = "weighted_1_ks_D_s_vs_r")
#   abline(v = 0, h = 0, lty = 6)
#   abline(coef = c(0,1), lty = 6)
#   abline(coef = c(0,-1), lty = 6)
#   for(i in 1:nrow(res1)){
#     points(res1[i,]$u, res1[i,]$v,pch=1,col="black",lwd = 0.2)
#   }
#   plot(res$u, res$v, pch=20,
#        col=res$color_3,
#        xlim=c(-range, range), ylim=c(-range,range),
#        xlab ="sensitive score", ylab = "resistance score",
#        main = "weighted_1_wass_dist_s_vs_r")
#   abline(v = 0, h = 0, lty = 6)
#   abline(coef = c(0,1), lty = 6)
#   abline(coef = c(0,-1), lty = 6)
#   for(i in 1:nrow(res1)){
#     points(res1[i,]$u, res1[i,]$v,pch=1,col="black",lwd = 0.2)
#   }

#   plot(res$u, res$v, pch=20,
#        col=res$color_4,
#        xlim=c(-range, range), ylim=c(-range,range),
#        xlab ="sensitive score", ylab = "resistance score",
#        main = "weighted_2_ks_D_s_vs_r")
#   abline(v = 0, h = 0, lty = 6)
#   abline(coef = c(0,1), lty = 6)
#   abline(coef = c(0,-1), lty = 6)
#   for(i in 1:nrow(res1)){
#     points(res1[i,]$u, res1[i,]$v,pch=1,col="black",lwd = 0.2)
#   }

#   plot(res$u, res$v, pch=20,
#        col=res$color_5,
#        xlim=c(-range, range), ylim=c(-range,range),
#        xlab ="sensitive score", ylab = "resistance score",
#        main = "weighted_2_wass_dist_s_vs_r")
#   abline(v = 0, h = 0, lty = 6)
#   abline(coef = c(0,1), lty = 6)
#   abline(coef = c(0,-1), lty = 6)
#   for(i in 1:nrow(res1)){
#     points(res1[i,]$u, res1[i,]$v,pch=1,col="black",lwd = 0.2)
#   }

#   dev.off()

#   pdf(paste0(savepath, "overall_drug_ITS_analysis_ks_p_s_srand_r_rrand.pdf"), width = 8, height = 12)
#   par(mfrow = c(3,2))
#   range =max(abs(res$u),abs(res$v))
#   plot(res$u, res$v, pch=20,
#        # col=res$color_1,
#        xlim=c(-range, range), ylim=c(-range,range),
#        xlab ="sensitive score", ylab = "resistance score")
#   abline(v = 0, h = 0, lty = 6)
#   abline(coef = c(0,1), lty = 6)
#   abline(coef = c(0,-1), lty = 6)
#   # points(res2$u, res2$v,pch=20,col="red")
#   # text(res$u, res$v, gsub("drug","D",rownames(res)),
#   #      cex = 0.7, pos = 4)
#   # draw.circle(0,0,0.5, lty = 6)
#   # draw.circle(0,0,1, lty = 6)

#   plot(res$u, res$v, pch=20,
#        col=res$color_1,
#        xlim=c(-range, range), ylim=c(-range,range),
#        xlab ="sensitive score", ylab = "resistance score",
#        main = "mean_s_minus_mean_r")
#   abline(v = 0, h = 0, lty = 6)
#   abline(coef = c(0,1), lty = 6)
#   abline(coef = c(0,-1), lty = 6)
#   for(i in 1:nrow(res2)){
#     points(res2[i,]$u, res2[i,]$v,pch=1,col="black",lwd = 0.2)
#   }

#   plot(res$u, res$v, pch=20,
#        col=res$color_2,
#        xlim=c(-range, range), ylim=c(-range,range),
#        xlab ="sensitive score", ylab = "resistance score",
#        main = "weighted_1_ks_D_s_vs_r")
#   abline(v = 0, h = 0, lty = 6)
#   abline(coef = c(0,1), lty = 6)
#   abline(coef = c(0,-1), lty = 6)
#   for(i in 1:nrow(res2)){
#     points(res2[i,]$u, res2[i,]$v,pch=1,col="black",lwd = 0.2)
#   }
#   plot(res$u, res$v, pch=20,
#        col=res$color_3,
#        xlim=c(-range, range), ylim=c(-range,range),
#        xlab ="sensitive score", ylab = "resistance score",
#        main = "weighted_1_wass_dist_s_vs_r")
#   abline(v = 0, h = 0, lty = 6)
#   abline(coef = c(0,1), lty = 6)
#   abline(coef = c(0,-1), lty = 6)
#   for(i in 1:nrow(res2)){
#     points(res2[i,]$u, res2[i,]$v,pch=1,col="black",lwd = 0.2)
#   }

#   plot(res$u, res$v, pch=20,
#        col=res$color_4,
#        xlim=c(-range, range), ylim=c(-range,range),
#        xlab ="sensitive score", ylab = "resistance score",
#        main = "weighted_2_ks_D_s_vs_r")
#   abline(v = 0, h = 0, lty = 6)
#   abline(coef = c(0,1), lty = 6)
#   abline(coef = c(0,-1), lty = 6)
#   for(i in 1:nrow(res2)){
#     points(res2[i,]$u, res2[i,]$v,pch=1,col="black",lwd = 0.2)
#   }

#   plot(res$u, res$v, pch=20,
#        col=res$color_5,
#        xlim=c(-range, range), ylim=c(-range,range),
#        xlab ="sensitive score", ylab = "resistance score",
#        main = "weighted_2_wass_dist_s_vs_r")
#   abline(v = 0, h = 0, lty = 6)
#   abline(coef = c(0,1), lty = 6)
#   abline(coef = c(0,-1), lty = 6)
#   for(i in 1:nrow(res2)){
#     points(res2[i,]$u, res2[i,]$v,pch=1,col="black",lwd = 0.2)
#   }

#   dev.off()




#   pdf(paste0(savepath, "overall_drug_ITS_analysis_wass_dist_pval_s_vs_r.pdf"), width = 8, height = 12)
#   par(mfrow = c(3,2))
#   range =max(abs(res$u),abs(res$v))
#   plot(res$u, res$v, pch=20,
#        # col=res$color_1,
#        xlim=c(-range, range), ylim=c(-range,range),
#        xlab ="sensitive score", ylab = "resistance score")
#   abline(v = 0, h = 0, lty = 6)
#   abline(coef = c(0,1), lty = 6)
#   abline(coef = c(0,-1), lty = 6)
#   # points(res3$u, res3$v,pch=20,col="red")
#   # text(res$u, res$v, gsub("drug","D",rownames(res)),
#   #      cex = 0.7, pos = 4)
#   # draw.circle(0,0,0.5, lty = 6)
#   # draw.circle(0,0,1, lty = 6)

#   plot(res$u, res$v, pch=20,
#        col=res$color_1,
#        xlim=c(-range, range), ylim=c(-range,range),
#        xlab ="sensitive score", ylab = "resistance score",
#        main = "mean_s_minus_mean_r")
#   abline(v = 0, h = 0, lty = 6)
#   abline(coef = c(0,1), lty = 6)
#   abline(coef = c(0,-1), lty = 6)
#   for(i in 1:nrow(res3)){
#     points(res3[i,]$u, res3[i,]$v,pch=1,col="black",lwd = 0.2)
#   }

#   plot(res$u, res$v, pch=20,
#        col=res$color_2,
#        xlim=c(-range, range), ylim=c(-range,range),
#        xlab ="sensitive score", ylab = "resistance score",
#        main = "weighted_1_ks_D_s_vs_r")
#   abline(v = 0, h = 0, lty = 6)
#   abline(coef = c(0,1), lty = 6)
#   abline(coef = c(0,-1), lty = 6)
#   for(i in 1:nrow(res3)){
#     points(res3[i,]$u, res3[i,]$v,pch=1,col="black",lwd = 0.2)
#   }
#   plot(res$u, res$v, pch=20,
#        col=res$color_3,
#        xlim=c(-range, range), ylim=c(-range,range),
#        xlab ="sensitive score", ylab = "resistance score",
#        main = "weighted_1_wass_dist_s_vs_r")
#   abline(v = 0, h = 0, lty = 6)
#   abline(coef = c(0,1), lty = 6)
#   abline(coef = c(0,-1), lty = 6)
#   for(i in 1:nrow(res3)){
#     points(res3[i,]$u, res3[i,]$v,pch=1,col="black",lwd = 0.2)
#   }

#   plot(res$u, res$v, pch=20,
#        col=res$color_4,
#        xlim=c(-range, range), ylim=c(-range,range),
#        xlab ="sensitive score", ylab = "resistance score",
#        main = "weighted_2_ks_D_s_vs_r")
#   abline(v = 0, h = 0, lty = 6)
#   abline(coef = c(0,1), lty = 6)
#   abline(coef = c(0,-1), lty = 6)
#   for(i in 1:nrow(res3)){
#     points(res3[i,]$u, res3[i,]$v,pch=1,col="black",lwd = 0.2)
#   }

#   plot(res$u, res$v, pch=20,
#        col=res$color_5,
#        xlim=c(-range, range), ylim=c(-range,range),
#        xlab ="sensitive score", ylab = "resistance score",
#        main = "weighted_2_wass_dist_s_vs_r")
#   abline(v = 0, h = 0, lty = 6)
#   abline(coef = c(0,1), lty = 6)
#   abline(coef = c(0,-1), lty = 6)
#   for(i in 1:nrow(res3)){
#     points(res3[i,]$u, res3[i,]$v,pch=1,col="black",lwd = 0.2)
#   }

#   dev.off()


#   pdf(paste0(savepath, "overall_drug_ITS_analysis_wass_dist_pval_s_srand_r_rand.pdf"), width = 8, height = 12)
#   par(mfrow = c(3,2))
#   range =max(abs(res$u),abs(res$v))
#   plot(res$u, res$v, pch=20,
#        # col=res$color_1,
#        xlim=c(-range, range), ylim=c(-range,range),
#        xlab ="sensitive score", ylab = "resistance score")
#   abline(v = 0, h = 0, lty = 6)
#   abline(coef = c(0,1), lty = 6)
#   abline(coef = c(0,-1), lty = 6)
#   # points(res4$u, res4$v,pch=20,col="red")
#   # text(res$u, res$v, gsub("drug","D",rownames(res)),
#   #      cex = 0.7, pos = 4)
#   # draw.circle(0,0,0.5, lty = 6)
#   # draw.circle(0,0,1, lty = 6)

#   plot(res$u, res$v, pch=20,
#        col=res$color_1,
#        xlim=c(-range, range), ylim=c(-range,range),
#        xlab ="sensitive score", ylab = "resistance score",
#        main = "mean_s_minus_mean_r")
#   abline(v = 0, h = 0, lty = 6)
#   abline(coef = c(0,1), lty = 6)
#   abline(coef = c(0,-1), lty = 6)
#   for(i in 1:nrow(res4)){
#     points(res4[i,]$u, res4[i,]$v,pch=1,col="black",lwd = 0.2)
#   }

#   plot(res$u, res$v, pch=20,
#        col=res$color_2,
#        xlim=c(-range, range), ylim=c(-range,range),
#        xlab ="sensitive score", ylab = "resistance score",
#        main = "weighted_1_ks_D_s_vs_r")
#   abline(v = 0, h = 0, lty = 6)
#   abline(coef = c(0,1), lty = 6)
#   abline(coef = c(0,-1), lty = 6)
#   for(i in 1:nrow(res4)){
#     points(res4[i,]$u, res4[i,]$v,pch=1,col="black",lwd = 0.2)
#   }
#   plot(res$u, res$v, pch=20,
#        col=res$color_3,
#        xlim=c(-range, range), ylim=c(-range,range),
#        xlab ="sensitive score", ylab = "resistance score",
#        main = "weighted_1_wass_dist_s_vs_r")
#   abline(v = 0, h = 0, lty = 6)
#   abline(coef = c(0,1), lty = 6)
#   abline(coef = c(0,-1), lty = 6)
#   for(i in 1:nrow(res4)){
#     points(res4[i,]$u, res4[i,]$v,pch=1,col="black",lwd = 0.2)
#   }

#   plot(res$u, res$v, pch=20,
#        col=res$color_4,
#        xlim=c(-range, range), ylim=c(-range,range),
#        xlab ="sensitive score", ylab = "resistance score",
#        main = "weighted_2_ks_D_s_vs_r")
#   abline(v = 0, h = 0, lty = 6)
#   abline(coef = c(0,1), lty = 6)
#   abline(coef = c(0,-1), lty = 6)
#   for(i in 1:nrow(res4)){
#     points(res4[i,]$u, res4[i,]$v,pch=1,col="black",lwd = 0.2)
#   }

#   plot(res$u, res$v, pch=20,
#        col=res$color_5,
#        xlim=c(-range, range), ylim=c(-range,range),
#        xlab ="sensitive score", ylab = "resistance score",
#        main = "weighted_2_wass_dist_s_vs_r")
#   abline(v = 0, h = 0, lty = 6)
#   abline(coef = c(0,1), lty = 6)
#   abline(coef = c(0,-1), lty = 6)
#   for(i in 1:nrow(res4)){
#     points(res4[i,]$u, res4[i,]$v,pch=1,col="black",lwd = 0.2)
#   }

#   dev.off()

# res = res[order(res$u_minus_v, decreasing = T), ]
# l1 = head(res, nrow(res) * 0.1)$pert_iname
# res = res[order(res$score_2, decreasing = T), ]
# l2 = head(res, nrow(res) * 0.1)$pert_iname
# res = res[order(res$score_3, decreasing = T), ]
# l3 = head(res, nrow(res) * 0.1)$pert_iname
# res = res[order(res$score_4, decreasing = T), ]
# l4 = head(res, nrow(res) * 0.1)$pert_iname
# res = res[order(res$score_5, decreasing = T), ]
# l5 = head(res, nrow(res) * 0.1)$pert_iname

# length(intersect(l1,l2))
# length(intersect(l1,l3))
# length(intersect(l1,l4))
# length(intersect(l1,l5))

# length(intersect(l2,l3))
# length(intersect(l2,l4))
# length(intersect(l2,l5))

# length(intersect(l3,l4))
# length(intersect(l3,l5))

# length(intersect(l4,l5))

# }



# not in use -------------------------------------------------------------------
drug_distri_distance <- function(df_immunesig,
                                 res_save,
                                 savepath){
  
  
  library(plyr)
  library(philentropy)
  set.seed(1234)
  
  res_save = res_save[order(res_save$u, decreasing = T), ]
  res_save$rank_u = 1:nrow(res_save)
  
  res_save = res_save[order(res_save$v, decreasing = F), ]
  res_save$rank_v = 1:nrow(res_save)
  
  res_save$u_minus_v = res_save$u - res_save$v
  
  res_save = res_save[order(res_save$u_minus_v, decreasing = T), ]
  res_save$rank_u_minus_v = 1:nrow(res_save)
  
  
  res_save = res_save[order(res_save$u_minus_v, decreasing = T), ]
  
  df_1 = df_immunesig # [c("drug_index", "mNES_diff", "flag")]
  # df_1$mNES_diff_max_min = (df_1$mNES_diff - min(df_1$mNES_diff))/(max(df_1$mNES_diff) - min(df_1$mNES_diff))
  res_save$ks_D_s_vs_r = 0
  res_save$ks_p_s_vs_r = 0
  res_save$rand_s_minus_r = 0
  res_save$ks_D_s_vs_srand = 0
  res_save$ks_p_s_vs_srand = 0
  res_save$ks_D_r_vs_rrand = 0
  res_save$ks_p_r_vs_rrand = 0
  # res_save$KL_hist_s_vs_r = 0
  # res_save$KL_hist_s_vs_srand = 0
  # res_save$KL_hist_r_vs_rrand = 0
  res_save$KL_density_s_vs_r = 0
  res_save$KL_density_s_vs_srand = 0
  res_save$KL_density_r_vs_rrand = 0
  res_save$wass_dist_s_vs_r = 0
  res_save$wass_dist_pval_s_vs_r = 0
  res_save$wass_dist_s_vs_srand = 0
  res_save$wass_dist_pval_s_vs_srand = 0
  res_save$wass_dist_r_vs_rrand = 0
  res_save$wass_dist_pval_r_vs_rrand = 0
  
  dist_plot=list()
  
  s_obs_mean = mean(df_1[df_1$flag == "sensitive", ]$mNES_diff)
  s_obs_sd = sd(df_1[df_1$flag == "sensitive", ]$mNES_diff)
  r_obs_mean = mean(df_1[df_1$flag == "resistant", ]$mNES_diff)
  r_obs_sd = sd(df_1[df_1$flag == "resistant", ]$mNES_diff)
  s_minus_r_mean = mean(res_save$u_minus_v)
  s_minus_r_sd = sd(res_save$u_minus_v)
  
  sim.fun <-function (m,f,...) { 
    sample <- 1:m 
    for (i in 1:m) { 
      sample[i] <-f(...) 
    } 
    sample 
  } 
  
  f <- function(n=10,mu=0,sigma=1){ 
    r=rnorm(n,mu,sigma)
    # mean(r)
    (mean(r)-mu)/(sigma/sqrt(n)) 
  }
  
  # rand_s_minus_r=list()
  
  for(i in 1:length(res_save$drug_index)){
    
    d = res_save$drug_index[i]
    print(paste0("Start computing for ", d))
    
    df = df_1[df_1$drug_index == d,]
    s = df[df$flag == "sensitive" ,]
    r = df[df$flag == "resistant" ,]
    
    ks_res_s_vs_r = ks.test(s$mNES_diff,r$mNES_diff, alternative = "less")
    res_save[res_save$drug_index == d,]$ks_D_s_vs_r = ks_res_s_vs_r$statistic
    res_save[res_save$drug_index == d,]$ks_p_s_vs_r = ks_res_s_vs_r$p.value
    
    df_obs <- df[c("mNES_diff", "flag")]
    df_obs$type <- paste0(df_obs$flag, "_obs")
    
    
    rand_s <- sim.fun(1000,f, nrow(s), s_obs_mean, s_obs_sd)   # 模拟1000个样本量为30的N(5,4)随机数
    rand_r <- sim.fun(1000,f, nrow(r), r_obs_mean, r_obs_sd)   # 模拟1000个样本量为30的N(5,4)随机数
    
    res_save[res_save$drug_index == d,]$rand_s_minus_r <- mean(rand_s) - mean(rand_r)
    df_rand <- data.frame(mNES_diff = c(rand_s, rand_r),
                          flag = c(rep("sensitive", length(rand_s)), rep("resistant", length(rand_r))))
    
    ks_res_s_vs_srand = ks.test(s$mNES_diff, rand_s, alternative = "less")
    res_save[res_save$drug_index == d,]$ks_D_s_vs_srand = ks_res_s_vs_srand$statistic
    res_save[res_save$drug_index == d,]$ks_p_s_vs_srand = ks_res_s_vs_srand$p.value
    
    ks_res_r_vs_rrand = ks.test(r$mNES_diff, rand_r, alternative = "greater")
    res_save[res_save$drug_index == d,]$ks_D_r_vs_rrand = ks_res_r_vs_rrand$statistic
    res_save[res_save$drug_index == d,]$ks_p_r_vs_rrand = ks_res_r_vs_rrand$p.value
    
    df_rand$type <- paste0(df_rand$flag, "_rand")
    df2 = rbind(df_obs, df_rand)
    
    
    
    #     # Kulback-Leibler Divergence between P and Q
    
    # library(philentropy)
    
    # p = hist(df_obs[df_obs$flag == "sensitive",]$mNES_diff, xlim=c(-10,10))[c("mids",'density')]
    # p = do.call(cbind, p)
    # names(p) = c("mids", 'p')
    # q = hist(df_obs[df_obs$flag == "resistant",]$mNES_diff, xlim=c(-10,10))[c("mids",'density')]
    # q = do.call(cbind, q)
    # names(q) = c("mids", 'q')
    # tmp = merge(p,q,by='mids',all=T)
    # if(length(which(is.na(tmp$density.y))) >0 )tmp[is.na(tmp$density.y),]$density.y = 0
    # if(length(which(is.na(tmp$density.x))) >0 )tmp[is.na(tmp$density.x),]$density.x = 0
    # tmp = as.data.frame(apply(tmp,2,function(x)x/sum(x)))
    
    # x <- rbind(tmp$density.x, tmp$density.y)
    # res_save[res_save$drug_index == d,]$KL_hist_s_vs_r = KL(x, unit = "log2")
    
    # p_rand_s = hist(rand_s, xlim=c(-10,10))[c("mids",'density')]
    # p_rand_s = do.call(cbind, p_rand_s)
    # names(p_rand_s) = c("mids", 'p_rand_s')
    # tmp = merge(p, p_rand_s, by='mids',all=T)
    # if(length(which(is.na(tmp$density.y))) >0 )tmp[is.na(tmp$density.y),]$density.y = 0
    # if(length(which(is.na(tmp$density.x))) >0 )tmp[is.na(tmp$density.x),]$density.x = 0
    # tmp = as.data.frame(apply(tmp,2,function(x)x/sum(x)))
    # x <- rbind(tmp$density.x, tmp$density.y)
    # res_save[res_save$drug_index == d,]$KL_hist_s_vs_srand = KL(x, unit = "log2")
    
    
    # q_rand_r = hist(rand_r, xlim=c(-10,10))[c("mids",'density')]
    # q_rand_r = do.call(cbind, q_rand_r)
    # names(q_rand_r) = c("mids", 'q_rand_r')
    # tmp = merge(q, q_rand_r, by='mids',all=T)
    # if(length(which(is.na(tmp$density.y))) >0 )tmp[is.na(tmp$density.y),]$density.y = 0
    # if(length(which(is.na(tmp$density.x))) >0 )tmp[is.na(tmp$density.x),]$density.x = 0
    # tmp = as.data.frame(apply(tmp,2,function(x)x/sum(x)))    
    # x <- rbind(tmp$density.x, tmp$density.y)
    # res_save[res_save$drug_index == d,]$KL_hist_r_vs_rrand = KL(x, unit = "log2")
    
    
    
    
    # Kulback-Leibler Divergence between P and Q
    
    # plot(density(df_obs[df_obs$flag == "sensitive",]$mNES_diff, bw = 1))
    # lines(density(df_obs[df_obs$flag == "resistant",]$mNES_diff, bw = 1), col='green')
    
    p = density(df_obs[df_obs$flag == "sensitive",]$mNES_diff)$y 
    p = p/sum(p)
    q = density(df_obs[df_obs$flag == "resistant",]$mNES_diff)$y
    q = q/sum(q)
    x <- rbind(p, q)
    res_save[res_save$drug_index == d,]$KL_density_s_vs_r = KL(x)
    
    p_rand_s = density(rand_s)$y 
    p_rand_s = p_rand_s/sum(p_rand_s)
    x <- rbind(p, p_rand_s)
    res_save[res_save$drug_index == d,]$KL_density_s_vs_srand = KL(x)
    
    q_rand_r = density(rand_r)$y 
    q_rand_r = q_rand_r/sum(q_rand_r)
    x <- rbind(q, q_rand_r)
    res_save[res_save$drug_index == d,]$KL_density_r_vs_rrand = KL(x)
    
    library(waddR)
    
    spec.output<-c("pval","d.wass^2","perc.loc","perc.size","perc.shape")      
    wass_dist_s_vs_r = wasserstein.test(df_obs[df_obs$flag == "sensitive",]$mNES_diff,
                                        df_obs[df_obs$flag == "resistant",]$mNES_diff)[spec.output]
    
    res_save[res_save$drug_index == d,]$wass_dist_s_vs_r = wass_dist_s_vs_r[2]
    res_save[res_save$drug_index == d,]$wass_dist_pval_s_vs_r = wass_dist_s_vs_r[1]
    
    wass_dist_s_vs_srand = wasserstein.test(df_obs[df_obs$flag == "sensitive",]$mNES_diff, rand_s)[spec.output]
    
    res_save[res_save$drug_index == d,]$wass_dist_s_vs_srand = wass_dist_s_vs_srand[2]
    res_save[res_save$drug_index == d,]$wass_dist_pval_s_vs_srand = wass_dist_s_vs_srand[1]
    
    wass_dist_r_vs_rrand = wasserstein.test(df_obs[df_obs$flag == "resistant",]$mNES_diff, rand_r)[spec.output]
    
    res_save[res_save$drug_index == d,]$wass_dist_r_vs_rrand = wass_dist_r_vs_rrand[2]
    res_save[res_save$drug_index == d,]$wass_dist_pval_r_vs_rrand = wass_dist_r_vs_rrand[1]
    
    
    p1 = ggplot(df2, aes(x = mNES_diff, colour = type)) + 
      stat_ecdf(geom = "step", pad = FALSE)+
      theme_bw()+ 
      ggtitle(paste0(d, " ", res_save[res_save$drug_index == d, ]$pert_iname,
                     # "\nu_minus_v = ", res_save[res_save$drug_index == d,]$u_minus_v, 
                     "\nks_D_s_vs_r = ", ks_res_s_vs_r$statistic,
                     "\nks_p_s_vs_r = ", ks_res_s_vs_r$p.value,
                     "\nks_D_s_vs_srand = ", ks_res_s_vs_srand$statistic,
                     "\nks_p_s_vs_srand = ", ks_res_s_vs_srand$p.value,
                     "\nks_D_r_vs_rrand = ", ks_res_r_vs_rrand$statistic,
                     "\nks_p_r_vs_rrand = ", ks_res_r_vs_rrand$p.value))
    
    cdat <- ddply(df, "flag", summarise, rating.mean=mean(mNES_diff))
    names(cdat) <- c("type", "mean")
    cdat$type <- c("resistant_obs", "sensitive_obs")
    p2 = ggplot(df2, aes(x = mNES_diff, colour = type)) + 
      geom_density(alpha = 0.3) +
      geom_vline(data=cdat, aes(xintercept=mean,  colour=type),
                 linetype="dashed", size=1) +
      theme_bw()+ 
      ggtitle(paste0(# d, " ", res_save[res_save$drug_index == d, ]$pert_iname,
        "\nu_minus_v = ", res_save[res_save$drug_index == d,]$u_minus_v, 
        "\nwasserstein_dis_s_vs_r = ", wass_dist_s_vs_r[2],
        "\nwasserstein_dis_p_s_vs_r = ", wass_dist_s_vs_r[1],
        "\nwasserstein_dis_s_vs_srand = ", wass_dist_s_vs_srand[2],
        "\nwasserstein_dis_p_s_vs_srand = ", wass_dist_s_vs_srand[1],
        "\nwasserstein_dis_r_vs_rrand = ", ks_res_r_vs_rrand[2],
        "\nwasserstein_dis_p_r_vs_rrand = ", ks_res_r_vs_rrand[1]))
    
    dist_plot[[i]] <- p2/p1
    
  }
  
  # res_save = res_save[order(res_save$ks_D, decreasing = T), ]
  # res_save$KS_D_RANK = 1:nrow(res_save)
  
  # res_save = res_save[order(res_save$ks_p, decreasing = F), ]
  # res_save$KS_p_RANK = 1:nrow(res_save)
  
  pdf(paste0(savepath, "/distribution_ecdf_plot.pdf"),
      width = 5,
      height = 8, 
      onefile = TRUE)
  for(x in seq(length(p))){
    print(p[[x]])
  }
  dev.off()
  
  # plot(density(res_save$u_minus_v))
  # lines(density(res_save$rand_s_minus_r), col = "green")
  
  write.table(res_save, paste0(savepath,"/ks_analysis_res.txt"),
              sep = "\t",  quote = F, row.names = F)
  
}




drug_MF <- function(df_immunesig,
                    res_save,
                    matrix = c("mm", "svd"), # mm: matrix multiplication
                    type = c( "weighted"),
                    savepath){
  
  df_s = df_immunesig[grep("01",df_immunesig$setindex),]
  df_r = df_immunesig[grep("02",df_immunesig$setindex),]
  
  df_s = df_s[c("drug_index","immune_sig","mNES_diff")]
  names(df_s) = c(c("drug_index","immune_sig_s","mNES_diff_s"))
  df_r = df_r[c("drug_index","immune_sig","mNES_diff")]
  names(df_r) = c(c("drug_index","immune_sig_r","mNES_diff_r"))
  
  
  df = merge(unique(df_s[c("drug_index","mNES_diff_s")]),
             unique(df_r[c("drug_index","mNES_diff_r")]),
             by= "drug_index")
  
  mat_s <-  reshape(df_s,  timevar = "immune_sig_s",
                    idvar = "drug_index",direction = "wide")
  rownames(mat_s) <- mat_s$drug_index
  mat_s <- mat_s[-1]
  
  mat_r <-  reshape(df_r,  timevar = "immune_sig_r",
                    idvar = "drug_index",direction = "wide")
  rownames(mat_r) <- mat_r$drug_index
  mat_r <- mat_r[-1]
  
  mat_s2 <- as.matrix(mat_s[is.element(rownames(mat_s), intersect(rownames(mat_s), rownames(mat_r) )),])
  mat_r2 <- as.matrix(mat_r[is.element(rownames(mat_r), intersect(rownames(mat_s), rownames(mat_r) )),])
  
  mat_s2[is.na(mat_s2)] = 0
  mat_r2[is.na(mat_r2)] = 0
  
  # filter out column with too many zero  ----------------------------------------------------
  
  mat_s3 = mat_s2[,which(colSums(mat_s2==0) < 0.9 * nrow(mat_s2))]
  mat_r3 = mat_r2[,which(colSums(mat_r2==0) < 0.9 * nrow(mat_r2))]
  
  ITSmat <- as.matrix(cbind(mat_s3, mat_r3))
  colnames(ITSmat) <- gsub("mNES_diff_s.", "", colnames(ITSmat))
  colnames(ITSmat) <- gsub("mNES_diff_r.", "", colnames(ITSmat))
  
  
  ITS_type = unique(df_immunesig[c("immune_sig", "flag")])
  ITS_type = table(ITS_type)
  ITS_type_mat = data.frame(sensitive = ITS_type[,2],
                            resistant = ITS_type[,1])
  rownames(ITS_type_mat) = rownames(ITS_type)
  ITS_type_mat = ITS_type_mat[intersect(rownames(ITS_type_mat), colnames(ITSmat)),]
  # ITS_type_mat[ITS_type_mat$resistant == 1, ]$resistant = -1/nrow(ITS_type_mat[ITS_type_mat$resistant == 1, ])
  # ITS_type_mat[ITS_type_mat$sensitive == 1, ]$sensitive = 1/nrow(ITS_type_mat[ITS_type_mat$sensitive == 1, ])
  
  
  # WEIGHTED
  if(type == "weighted"){
    ITS_type_mat_s = unique(df_immunesig[df_immunesig$flag == "sensitive", ][c("immune_sig", "weight")])
    names(ITS_type_mat_s) = c("immune_sig", "sensitive")
    ITS_type_mat_r = unique(df_immunesig[df_immunesig$flag == "resistant", ][c("immune_sig", "weight")])
    names(ITS_type_mat_r) = c("immune_sig", "resistant")
    ITS_type_mat = merge(ITS_type_mat_s, ITS_type_mat_r, by = 'immune_sig', all = T)
    ITS_type_mat[is.na(ITS_type_mat)] = 0
    rownames(ITS_type_mat) = ITS_type_mat$immune_sig
    ITS_type_mat = ITS_type_mat[-1]
  }
  
  
  if(method == "svd"){
    
    A <- ITSmat # matrix(rep(1,18), nrow = 3)
    A = t(A)
    b=as.matrix(ITS_type_mat)
    
    # 欠定方程组，即方程个数少于变量个数的方程组，可以使用SVD法求解。
    sol.svd <- svd(A)
    #获取U D V各个值
    U<-sol.svd$u
    D<-sol.svd$d
    V<-sol.svd$v
    C<-t(U)%*%b
    Y<-C/D
    X<-V%*%Y
    rownames(X) = colnames(A)
    
    
    res = data.frame(drug_index = rownames(X), sensitive = X[,1], resistant = X[,2])
    res$drug_index = rownames(X)
    res$s_plus_r = res$sensitive + res$resistant
    res <- res[order(res$s_plus_r),]
    res$s_plus_r_rank <- nrow(res):1
    
    res <- inner_join(res_save, res)
    res <- res[order(res$s_plus_r_rank),]
    
    
    b=as.matrix(ITS_type_mat[,1])# - ITS_type_mat[,2])
    rownames(b) = rownames(ITS_type_mat)
    
    # 欠定方程组，即方程个数少于变量个数的方程组，可以使用SVD法求解。
    sol.svd <- svd(A)
    #获取U D V各个值
    U<-sol.svd$u
    D<-sol.svd$d
    V<-sol.svd$v
    C<-t(U)%*%b
    Y<-C/D
    X<-V%*%Y
    rownames(X) = colnames(A)
    
    res = data.frame(drug_index = rownames(X), score = X[,1])
    res$drug_index = rownames(X)
    res <- res[order(res$score),]
    res$rank <- nrow(res):1
    res <- inner_join(res_save, res)
    res <- res[order(res$rank),]
    
  }else if(method == "mm"){
    A <- ITSmat # matrix(rep(1,18), nrow = 3)
    A = t(A)
    
    b = as.matrix(ITS_type_mat)
    # rownames(b) = ITS_type_mat$immune_sig
    b = b[intersect(rownames(A), rownames(b)), ]
    
    A = A[intersect(rownames(A), rownames(b)), ]
    
    tmp = t(A)%*%b
    res = data.frame(drug_index = rownames(tmp), sensitive = tmp[,1], resistant = tmp[,2])
    res$drug_index = rownames(X)
    # res$s_minus_r = scale(res$sensitive) - scale(res$resistant)
    res$s_minus_r = res$sensitive - res$resistant
    res <- res[order(res$s_minus_r),]
    res$s_plus_r_rank <- nrow(res):1
    
    res <- inner_join(res_save, res)
    res <- res[order(res$s_plus_r_rank),]
    
  }
  write.table(res, paste0(savepath,"/", method, "_analysis_res.txt"),
              sep = '\t', quote = F, row.names = F)
  
  return(res)
}



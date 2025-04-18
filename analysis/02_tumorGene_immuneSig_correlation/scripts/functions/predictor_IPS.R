
# Immunophenoscore (IPS) --------------------------------------------------

## Read expression data from tab-delimited text file,
## with official human gene symbols (HGNC) in the first columns
## and expression values (i.e. log2(TPM+1)) for each sample in the other columns
## Normalization:
## expMatrix_norm = log2(TPM_matrix+1) 

model_icb_IPS <- function(gene_expression, # non-log TPM or FPKM Expression
                          cancer_type,
                          sample_type,
                          sig_dir,
                          save_path) {
  
  # require(ggplot2)  
  # require(gridExtra)  
  
  if(!file.exists(save_path)){dir.create(save_path)}
  save_path = paste0(save_path,"IPS_result/")
  if(!file.exists(save_path)){dir.create(save_path)}
  
  ## calculate Immunophenoscore
  ipsmap <- function (x) {
    if (x <= 0) {
      ips <- 0
    } else {
      if (x >= 3) {
        ips <- 10
      } else {
        ips <- round(x * 10 / 3, digits = 0)
      }
    }
    return(ips)
  }
  
  ## Assign colors
  #   my_palette <-
  #     colorRampPalette(c("blue", "white", "red"))(n = 1000)
  #   mapcolors <- function (x) {
  #     za <- NULL
  #     if (x >= 3) {
  #       za = 1000
  #     } else {
  #       if (x <= -3) {
  #         za = 1
  #       } else {
  #         za = round(166.5 * x + 500.5, digits = 0)
  #       }
  #     }
  #     return(my_palette[za])
  #   }
  
  #   my_palette2 <- colorRampPalette(c("black", "white"))(n = 1000)
  #   mapbw <- function (x) {
  #     za2 <- NULL
  #     if (x >= 2) {
  #       za2 = 1000
  #     } else {
  #       if (x <= -2) {
  #         za2 = 1
  #       } else {
  #         za2 = round(249.75 * x + 500.5, digits = 0)
  #       }
  #     }
  #     return(my_palette2[za2])
  #   }
  
  
  gene_expression = log2(gene_expression+1) # - apply(log2(gene_expression+1),1,mean)
  
  gene_expression = as.data.frame(gene_expression)
  sample_names <- names(gene_expression)
  ## Read IPS genes and corresponding weights from tab-delimited text file "IPS_genes.txt"
  # For different
  IPSG <-
    read.table(
      paste0(sig_dir, "IPS_genes.txt"),
      header = TRUE,
      sep = "\t",
      dec = ".",
      check.names = FALSE
    )
  unique_ips_genes <- as.vector(unique(IPSG$NAME))
  
  IPS <- NULL
  MHC <- NULL
  CP <- NULL
  EC <- NULL
  SC <- NULL
  AZ <- NULL
  
  # Gene names in expression file
  GVEC <- row.names(gene_expression)
  # Genes names in IPS genes file
  VEC <- as.vector(IPSG$GENE)
  # Match IPS genes with genes in expression file
  ind <- which(is.na(match(VEC, GVEC)))
  # List genes missing or differently named
  MISSING_GENES <- VEC[ind]
  dat <- IPSG[ind, ]
  if (length(MISSING_GENES) > 0) {
    cat("differently named or missing genes for ", cancer_type,": ", MISSING_GENES, "\n")
  }
  # for (x in 1:length(ind)) {
  #   print(IPSG[ind, ])
  # }
  
  for (i in 1:ncol(gene_expression)) {
    # i=3
    GE <- gene_expression[[i]]
    mGE <- mean(GE)
    sGE <- sd(GE)
    Z1 <- (gene_expression[as.vector(IPSG$GENE), i] - mGE) / sGE
    W1 <- IPSG$WEIGHT
    WEIGHT <- NULL
    MIG <- NULL
    k <- 1
    
    for (gen in unique_ips_genes) {
      MIG[k] <-
        mean(Z1[which (as.vector(IPSG$NAME) == gen)], na.rm = TRUE)
      WEIGHT[k] <- mean(W1[which (as.vector(IPSG$NAME) == gen)])
      k <- k + 1
    }
    WG <- MIG * WEIGHT
    MHC[i] <- mean(WG[1:10], na.rm = TRUE)
    CP[i] <- mean(WG[11:20], na.rm = TRUE)
    EC[i] <- mean(WG[21:24], na.rm = TRUE)
    SC[i] <- mean(WG[25:26], na.rm = TRUE)
    AZ[i] <- sum(MHC[i], CP[i], EC[i], SC[i], na.rm = TRUE)
    IPS[i] <- ipsmap(AZ[i])
    
    ## Plot Immunophenogram
    # data_a2 <-
    #   data.frame (
    #     start = c(0, 2.5, 5, 7.5, 10, 15, seq(20, 39), 0, 10, 20, 30),
    #     end = c(2.5, 5, 7.5, 10, 15, seq(20, 40), 10, 20, 30, 40),
    #     y1 = c(rep(2.6, 26), rep(0.4, 4)),
    #     y2 = c(rep(5.6, 26), rep(2.2, 4)),
    #     z = c(MIG[c(21:26, 11:20, 1:10)],
    #           EC[i], SC[i], CP[i], MHC[i]),
    #     vcol = c(unlist(lapply(MIG[c(21:26, 11:20, 1:10)], mapcolors)),
    #              unlist(lapply(
    #                c(EC[i], SC[i], CP[i], MHC[i]), mapbw
    #              ))),
    #     label = c(unique_ips_genes[c(21:26, 11:20, 1:10)],
    #               "EC", "SC", "CP", "MHC")
    #   )
    # data_a2$label <-
    #   factor(data_a2$label, levels = unique(data_a2$label))
    
    # plot_a2 <- ggplot() +
    #   geom_rect(
    #     data = data_a2,
    #     mapping = aes(
    #       xmin = start,
    #       xmax = end,
    #       ymin = y1,
    #       ymax = y2,
    #       fill = label
    #     ),
    #     size = 0.5,
    #     color = "black",
    #     alpha = 1
    #   ) +
    #   coord_polar() +
    #   scale_y_continuous(limits = c(0, 6)) +
    #   scale_fill_manual(values = as.vector(data_a2$vcol), guide = FALSE) +
    #   theme_bw() +
    #   theme(
    #     panel.margin = unit(0, 'mm'),
    #     panel.grid.major = element_blank(),
    #     panel.grid.minor = element_blank(),
    #     panel.border = element_blank(),
    #     panel.background = element_blank(),
    #     axis.line = element_line(colour = "white"),
    #     axis.text = element_blank(),
    #     axis.ticks = element_blank()
    #   ) +
    #   geom_text(aes(x = 5, y = 1.3, label = "EC"), size = 4) +
    #   geom_text(aes(x = 15, y = 1.3, label = "SC"), size = 4) +
    #   geom_text(aes(x = 25, y = 1.3, label = "CP"), size = 4) +
    #   geom_text(aes(x = 35, y = 1.3, label = "MHC"), size = 4)
    
    # plot_a2 <- plot_a2 +
    #   geom_text(aes(x = 1.25, y = 4.1, label = "+ Act CD4"),
    #             angle = 78.75,
    #             size = 4) +
    #   geom_text(aes(x = 3.75, y = 4.1, label = "+ Act CD8"),
    #             angle = 56.25,
    #             size = 4) +
    #   geom_text(aes(x = 6.25, y = 4.1, label = "+ Tem CD4"),
    #             angle = 33.75,
    #             size = 4) +
    #   geom_text(aes(x = 8.75, y = 4.1, label = "+ Tem CD8"),
    #             angle = 11.25,
    #             size = 4) +
    #   geom_text(aes(x = 17.5, y = 4.1, label = "- MDSC"),
    #             angle = -67.5,
    #             size = 4) +
    #   geom_text(aes(x = 12.5, y = 4.1, label = "- Treg"),
    #             angle = -22.5,
    #             size = 4)
    
    # plot_a3 <-
    #   plot_a2 + geom_text(aes(x = 20.5, y = 4.1, label = "PD-1 -"),
    #                       angle = 85.5,
    #                       size = 4) +
    #   geom_text(aes(x = 21.5, y = 4.1, label = "CTLA4 -"),
    #             angle = 76.5,
    #             size = 4) +
    #   geom_text(aes(x = 22.5, y = 4.1, label = "LAG3 -"),
    #             angle = 67.5,
    #             size = 4) +
    #   geom_text(aes(x = 23.5, y = 4.1, label = "TIGIT -"),
    #             angle = 58.5,
    #             size = 4) +
    #   geom_text(aes(x = 24.5, y = 4.1, label = "TIM3 -"),
    #             angle = 49.5,
    #             size = 4) +
    #   geom_text(aes(x = 25.5, y = 4.1, label = "PD-L1 -"),
    #             angle = 40.5,
    #             size = 4) +
    #   geom_text(aes(x = 26.5, y = 4.1, label = "PD-L2 -"),
    #             angle = 31.5,
    #             size = 4) +
    #   geom_text(aes(x = 27.5, y = 4.1, label = "CD27 +"),
    #             angle = 22.5,
    #             size = 4) +
    #   geom_text(aes(x = 28.5, y = 4.1, label = "ICOS +"),
    #             angle = 13.5,
    #             size = 4) +
    #   geom_text(aes(x = 29.5, y = 4.1, label = "IDO1 -"),
    #             angle = 4.5,
    #             size = 4)
    
    # plot_a4 <-
    #   plot_a3 + geom_text(aes(x = 30.5, y = 4.1, label = "B2M +"),
    #                       angle = -4.5,
    #                       size = 4) +
    #   geom_text(aes(x = 31.5, y = 4.1, label = "TAP1 +"),
    #             angle = -13.5,
    #             size = 4) +
    #   geom_text(aes(x = 32.5, y = 4.1, label = "TAP2 +"),
    #             angle = -22.5,
    #             size = 4) +
    #   geom_text(aes(x = 33.5, y = 4.1, label = "HLA-A +"),
    #             angle = -31.5,
    #             size = 4) +
    #   geom_text(aes(x = 34.5, y = 4.1, label = "HLA-B +"),
    #             angle = -40.5,
    #             size = 4) +
    #   geom_text(aes(x = 35.5, y = 4.1, label = "HLA-C +"),
    #             angle = -49.5,
    #             size = 4) +
    #   geom_text(aes(x = 36.5, y = 4.1, label = "HLA-DPA1 +"),
    #             angle = -58.5,
    #             size = 4) +
    #   geom_text(aes(x = 37.5, y = 4.1, label = "HLA-DPB1 +"),
    #             angle = -67.5,
    #             size = 4) +
    #   geom_text(aes(x = 38.5, y = 4.1, label = "HLA-E +"),
    #             angle = -76.5,
    #             size = 4) +
    #   geom_text(aes(x = 39.5, y = 4.1, label = "HLA-F +"),
    #             angle = -85.5,
    #             size = 4)
    
    # plot_a5 <- plot_a4 +
    #   geom_text(
    #     aes(
    #       x = 0,
    #       y = 6,
    #       label = paste("Immunophenoscore: ", IPS[i], sep = "")
    #     ),
    #     angle = 0,
    #     size = 6,
    #     vjust = -0.5
    #   ) +
    #   theme(axis.title = element_blank())
    # plot_a <- plot_a5 +
    #   theme(plot.margin = unit(c(0, 0, 0, 0), "mm")) +
    #   geom_text(
    #     vjust = 1.15,
    #     hjust = 0,
    #     aes(
    #       x = 25.5,
    #       y = 6,
    #       label = "\n\n\n\n   MHC: Antigen Processing                                 EC: Effector Cells\n   CP: Checkpoints | Immunomodulators              SC: Suppressor Cells\n\n",
    #       hjust = 0
    #     ),
    #     size = 4
    #   )
    # #plot_a
    
    # ## Legend sample-wise (averaged) z-scores
    # data_b <- data.frame (
    #   start = rep(0, 23),
    #   end = rep(0.7, 23),
    #   y1 = seq(0, 22, by = 1),
    #   y2 = seq(1, 23, by = 1),
    #   z = seq(-3, 3, by = 6 / 22),
    #   vcol = c(unlist(lapply(
    #     seq(-3, 3, by = 6 / 22), mapcolors
    #   ))),
    #   label = LETTERS[1:23]
    # )
    # data_b_ticks <-
    #   data.frame(
    #     x = rep(1.2, 7),
    #     value = seq(-3, 3, by = 1),
    #     y = seq(0, 6, by = 1) * (22 / 6) + 0.5
    #   )
    # legendtheme <- theme(
    #   plot.margin = unit(c(2, 0, 2, 0), "inch"),
    #   panel.margin = unit(0, "null"),
    #   panel.grid.major = element_blank(),
    #   panel.grid.minor = element_blank(),
    #   panel.border = element_blank(),
    #   panel.background = element_blank(),
    #   axis.line = element_line(colour = "white"),
    #   axis.text = element_blank(),
    #   axis.ticks = element_blank(),
    #   axis.title.x = element_blank()
    # )
    
    # plot_b <- ggplot(hjust = 0) +
    #   geom_rect(
    #     data = data_b,
    #     mapping = aes(
    #       xmin = start,
    #       xmax = end,
    #       ymin = y1,
    #       ymax = y2,
    #       fill = label
    #     ),
    #     size = 0.5,
    #     color = "black",
    #     alpha = 1
    #   ) +
    #   scale_x_continuous(limits = c(0, 1.5), expand = c(0, 0)) +
    #   scale_fill_manual(values = as.vector(data_b$vcol), guide = FALSE) +
    #   geom_text(
    #     data = data_b_ticks,
    #     aes(x = x, y = y, label = value),
    #     hjust = "inward",
    #     size = 4
    #   ) +
    #   theme_bw() +
    #   legendtheme +
    #   ylab("Sample-wise (averaged) z-score")
    # #plot_b
    # ## Legend weighted z-scores
    # data_c <- data.frame (
    #   start = rep(0, 23),
    #   end = rep(0.7, 23),
    #   y1 = seq(0, 22, by = 1),
    #   y2 = seq(1, 23, by = 1),
    #   z = seq(-2, 2, by = 4 / 22),
    #   vcol = c(unlist(lapply(
    #     seq(-2, 2, by = 4 / 22), mapbw
    #   ))),
    #   label = LETTERS[1:23]
    # )
    # data_c_ticks <-
    #   data.frame(
    #     x = rep(1.2, 5),
    #     value = seq(-2, 2, by = 1),
    #     y = seq(0, 4, by = 1) * (22 / 4) + 0.5
    #   )
    # plot_c <- ggplot() +
    #   geom_rect(
    #     data = data_c,
    #     mapping = aes(
    #       xmin = start,
    #       xmax = end,
    #       ymin = y1,
    #       ymax = y2,
    #       fill = label
    #     ),
    #     size = 0.5,
    #     color = "black",
    #     alpha = 1
    #   ) +
    #   scale_x_continuous(limits = c(0, 1.5), expand = c(0, 0)) +
    #   scale_fill_manual(values = as.vector(data_c$vcol), guide = FALSE) +
    #   geom_text(
    #     data = data_c_ticks,
    #     aes(x = x, y = y, label = value),
    #     hjust = "inward",
    #     size = 4
    #   ) +
    #   theme_bw() +
    #   legendtheme +
    #   ylab("Weighted z-score")
    # #plot_c
    # ## Save plot to file (1 pdf file for each sample)
    # file_name <- paste0(save_path, "IPS_", sample_names[i], ".pdf")
    # pdf(file_name, width = 10, height = 8)
    # grid.arrange(plot_a,
    #              plot_b,
    #              plot_c,
    #              ncol = 3,
    #              widths = c(0.8, 0.1, 0.1))
    # dev.off()
  }
  
  DF <-
    data.frame(
      SAMPLE = sample_names,
      MHC = MHC,
      EC = EC,
      SC = SC,
      CP = CP,
      AZ = AZ,
      IPS = IPS
    )
  
  
  ## barplot(DF[[7]],DF[[1]],col=c(rep("yellow",7),rep("green",4)))
  
  
  write.csv(
    DF,
    file = paste0(save_path, "IPS_result_", cancer_type, "_", sample_type, ".csv"),
    quote = FALSE
  )
  # return(DF)
  
}


# Immunophenoscore (IPS) - mouse ------------------------------------------

## Read expression data from tab-delimited text file, with official human gene symbols (HGNC) in the first columns
## and expression values (i.e. log2(TPM+1)) for each sample in the other columns


model_icb_IPS_mouse <- function(gene_expression,
                                cancer_type,
                                sample_type,
                                sig_dir,
                                IPSG,
                                sig_type = c('original','modified_mouse'),
                                save_path) {
  
  
  #FPKM 2 TPM for matching the input of the model
  # fpkmToTpm <- function(fpkm)
  # {
  #   exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
  # }
  # 
  ## calculate Immunophenoscore
  ipsmap <- function (x) {
    if (x <= 0) {
      ips <- 0
    } else {
      if (x >= 3) {
        ips <- 10
      } else {
        ips <- round(x * 10 / 3, digits = 0)
      }
    }
    return(ips)
  }
  
  # gene_expression <- apply(gene_expression, 2, fpkmToTpm)
  gene_expression <- log2(gene_expression + 1)
  
  
  gene_expression = as.data.frame(gene_expression)
  sample_names <- names(gene_expression)
  ## Read IPS genes and corresponding weights from tab-delimited text file "IPS_genes.txt"
  # For different
  
  extra_mat <- data.frame(matrix(0, length(setdiff(IPSG$GENE, rownames(gene_expression))), ncol(gene_expression)))
  rownames(extra_mat) <-
    setdiff(IPSG$GENE, rownames(gene_expression))
  colnames(extra_mat) <- colnames(gene_expression)
  gene_expression <- rbind(gene_expression, extra_mat)
  IPSG <- IPSG[is.element(IPSG$GENE, rownames(gene_expression)), ]
  
  unique_ips_genes <- as.vector(unique(IPSG$NAME))
  
  IPS <- NULL
  MHC <- NULL
  CP <- NULL
  EC <- NULL
  SC <- NULL
  AZ <- NULL
  
  # Gene names in expression file
  GVEC <- row.names(gene_expression)
  # Genes names in IPS genes file
  VEC <- as.vector(IPSG$GENE)
  # Match IPS genes with genes in expression file
  ind <- which(is.na(match(VEC, GVEC)))
  # List genes missing or differently named
  MISSING_GENES <- VEC[ind]
  dat <- IPSG[ind, ]
  if (length(MISSING_GENES) > 0) {
    cat("differently named or missing genes: ", MISSING_GENES, "\n")
  }
  for (x in 1:length(ind)) {
    print(IPSG[ind, ])
  }
  
  # IPSG_expressed <- IPSG[is.element(IPSG$GENE,GVEC), ]
  
  for (i in 1:length(sample_names)) {
    # i=1
    GE <- gene_expression[[i]]
    mGE <- mean(GE)
    sGE <- sd(GE)
    Z1 <- (gene_expression[as.vector(IPSG$GENE), i] - mGE) / sGE
    W1 <- IPSG$WEIGHT
    WEIGHT <- NULL
    MIG <- NULL
    k <- 1
    
    for (gen in unique_ips_genes) {
      MIG[k] <-
        mean(Z1[which (as.vector(IPSG$NAME) == gen)], na.rm = TRUE)
      WEIGHT[k] <- mean(W1[which (as.vector(IPSG$NAME) == gen)])
      k <- k + 1
      
    }
    
    WG <- MIG * WEIGHT
    class_label = unique(IPSG[c('CLASS','NAME')])
    MHC[i] <- mean(WG[which(class_label$CLASS == 'MHC')])
    CP[i] <- mean(WG[which(class_label$CLASS == 'CP')])
    EC[i] <- mean(WG[which(class_label$CLASS == 'EC')])
    SC[i] <- mean(WG[which(class_label$CLASS == 'SC')])
    # MHC[i]<-mean(WG[1:3])
    # CP[i]<-mean(WG[4:13])
    # EC[i]<-mean(WG[14:17])
    # SC[i]<-mean(WG[18:19])
    AZ[i] <- sum(MHC[i], CP[i], EC[i], SC[i])
    IPS[i] <- ipsmap(AZ[i])
  }
  
  DF <-
    data.frame(
      SAMPLE = sample_names,
      MHC = MHC,
      EC = EC,
      SC = SC,
      CP = CP,
      AZ = AZ,
      IPS = IPS
    )
  
  
  # if(!file.exists(save_path)){dir.create(save_path)}
  # save_path = paste0(save_path,"IPS_result_mouse/")
  # if(!file.exists(save_path)){dir.create(save_path)}
  
  if(sig_type == 'original'){
    write.csv(
      DF,
      file = paste0(save_path, "IPS_result_mouse_", sig_type,".csv"),
      quote = FALSE)
  }else if(sig_type == 'modified_mouse'){
    write.csv(
      DF,
      file = paste0(save_path, "IPS_result_mouse_", sig_type,".csv"),
      quote = FALSE)
  }
  
  return(DF)
  
}


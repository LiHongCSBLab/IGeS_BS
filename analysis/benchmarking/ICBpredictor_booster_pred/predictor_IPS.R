
# Immunophenoscore (IPS) --------------------------------------------------

## Read expression data from tab-delimited text file,
## with official human gene symbols (HGNC) in the first columns
## and expression values (i.e. log2(TPM+1)) for each sample in the other columns
## Normalization:
## expMatrix_norm = log2(TPM_matrix+1) 

model_icb_IPS <- function(gene_expression, # non-log TPM or FPKM Expression
                          label,
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
  
  
  # gene_expression = log2(gene_expression+1) # - apply(log2(gene_expression+1),1,mean)
  gene_expression = gene_expression # - apply(log2(gene_expression+1),1,mean)
  
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
  
  DF$SAMPLE = gsub('\\.', ':', DF$SAMPLE)

  write.csv(
    DF,
    file = paste0(save_path, "IPS_result_", cancer_type, "_", sample_type, ".csv"),
    quote = FALSE
  )

  

  treated_group = DF[DF$SAMPLE %in% sample_label[sample_label$treatment == 1, ]$treatment_ID,]
  control_group = DF[DF$SAMPLE %in% sample_label[sample_label$treatment == 0, ]$treatment_ID,]
  res = apply(treated_group[,-1], 2, mean) - apply(control_group[,-1], 2, mean)

  write.csv(res, paste0(save_path, "IPS_pred.csv"),
            quote = F) 


  # return(DF)
  
}


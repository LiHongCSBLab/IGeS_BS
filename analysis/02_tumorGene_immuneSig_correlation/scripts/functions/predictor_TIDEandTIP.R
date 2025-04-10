
# ICB response prediction (TIDE and TIP) ----------------------------------
# ICB response prediction
# including:
# 1) TIDE including TIDE dysfunction score, TIDE exclusion score, expression of
#    PD1/PDL1, CD8A/B, CTL related genes and IFNG related genes.
# 2) TIP (hot/cold genes)
# Input expression matrix need to be FPKM, RPKM or TPM
# Expression matrix will go through log2-transformation and normalization within
# the function:
#   expMatrix_norm = log2(FPKM_matrix+1) - apply(log2(FPKM_matrix+1),1,mean)
#   or:
#   expMatrix_norm = log2(TPM_matrix+1) - apply(log2(TPM_matrix+1),1,mean)
#
# !! Warning: changing data sample size will influence TIDE or TIP score!!

model_icb_tide_TIP <- function(mat,
                               cancername = "SKCM",
                               sample_type,
                               sig_dir,
                               save_path) {
   
  mat = log2(mat+1) - apply(log2(mat+1),1,mean)
  # write.table(mat, paste0(save_path,"mat_log2TPM.txt"),row.names = T,quote = F, sep = '\t')
  mat = t(mat)
   
  TIDE_signature = read.csv(paste0(sig_dir, "data/TIDE_signature/signature"),
                            sep = "\t",
                            header = T)
  TIDE_signature$Gene = rownames(TIDE_signature)
  genename = read.csv(paste0(sig_dir,"data/TIDE_signature/genename.csv"),
                      sep=",",header = T)
  TIDE_signature = merge(TIDE_signature, genename, by = "Gene")
  
  signature.std = read.csv(paste0(sig_dir, "data/TIDE_signature/signature.std"), 
                           sep = "\t")
  
  
  # Calculate patient level dysfunction score
  dys.sig1 = TIDE_signature[!is.na(TIDE_signature$Dysfunction), ]
  mat2 = mat[, colnames(mat) %in% dys.sig1$Symbol]
  mat2 = mat2[, sort(colnames(mat2))]
  dys.sig2 = dys.sig1[is.element(dys.sig1$Symbol, colnames(mat2)), ]
  dys.sig2 = dys.sig2[order(dys.sig2$Symbol), ]
  
  dys.score1 = data.frame(apply(mat2, 1, function(x) {
    cor(x, dys.sig2$Dysfunction,
        use = "complete")
  }))
  
  names(dys.score1) = "dys.score"
  dys.score1$patient = rownames(dys.score1)
  dys.score1 = dys.score1[order(dys.score1$patient, decreasing = F), ]
  
  cancerlist = gsub(
    ".RNASeq.norm_subtract", "",
    rownames(signature.std)[grep(".RNASeq.norm_subtract", rownames(signature.std))]
  )
  
  if (is.element(cancername, cancerlist)) {
    dys.score1$Dysfunction = round(((dys.score1$dys.score) /
                                      signature.std[rownames(signature.std) ==
                                                      paste0(cancername, ".RNASeq.norm_subtract"), ]$Dysfunction), 2)
  } else{
    dys.score1$Dysfunction = round(((dys.score1$dys.score) /
                                      signature.std[rownames(signature.std) ==
                                                      paste0(cancername, ".RNASeq"),]$Dysfunction), 2)
  }
  

  # Calculate patient level exclusion score
  exc.sig = TIDE_signature[!is.na(TIDE_signature$Exclusion), ]
  exc.sig1 = exc.sig[exc.sig$Symbol %in% colnames(mat), ]
  exc.sig1 = exc.sig1[order(exc.sig1$Symbol), ]
  mat3 = mat[, colnames(mat) %in% exc.sig1$Symbol]
  mat3 = mat3[, sort(colnames(mat3))]
  exc.score = data.frame(apply(mat3, 1, function(x) {
    cor(x, exc.sig1$Exclusion,
        use = "complete")
  }))
  names(exc.score) = "exclusion_score"
  exc.score$patient = rownames(exc.score)
  exc.score = exc.score[order(exc.score$patient), ]
  # exc.score$Exclusion  = round(((exc.score$exclusion_score)/
  # signature.std[rownames(signature.std) == paste0( cancername,".RNASeq") ,
  # ]$Exclusion),2)
    
  if (is.element(cancername, cancerlist)) {
    exc.score$Exclusion  = round(((exc.score$exclusion_score) /
                                    signature.std[rownames(signature.std) ==
                                                    paste0(cancername, 
                                                           ".RNASeq.norm_subtract"),]$Exclusion), 2)
  } else{
    exc.score$Exclusion  = round(((exc.score$exclusion_score) /
                                    signature.std[rownames(signature.std) ==
                                                    paste0(cancername, ".RNASeq"), ]$Exclusion), 2)
  }
  TIDE_result = merge(dys.score1, exc.score, by = "patient")
  
  # Generating CTL sign
  mat = t(mat)
  CTL_genes = c("CD8A", "CD8B", "GZMA", "GZMB", "PRF1")
  CTLsig = mat[rownames(mat) %in% CTL_genes, ]
  CTL_class = data.frame(patient = colnames(mat), class = "low")
  for (i in 1:nrow(CTL_class)) {
    if (all(CTLsig[, CTL_class[i, 1]] > 0)) {
      CTL_class[i, 2] = "high"
    }
  }
  
  TIDE_result = merge(TIDE_result, CTL_class, by = "patient")
  TIDE_result$TIDE = TIDE_result$Exclusion
  TIDE_result[TIDE_result$class == "high", ]$TIDE =
    TIDE_result[TIDE_result$class == "high", ]$Dysfunction
  
  CTLsig = data.frame(t(mat[rownames(mat) %in% CTL_genes, ]))
  
  CTLsig$patient = rownames(CTLsig)
  
  CD8_genes = c("CD8A", "CD8B")
  # CD8sig = data.frame(CD8A = mat[is.element(rownames(mat) ,CD8_genes),])
  CD8sig = data.frame(t(mat[rownames(mat) %in% CD8_genes, ]))
  CD8sig$CD8 = apply(CD8sig, 1, mean)
  CD8sig$patient = rownames(CD8sig)
  
  IFNG_genes = c("IFNG", "STAT1", "IDO1", "CXCL10", "CXCL9", "HLA-DRA")
  IFNGsig = data.frame(t(mat[rownames(mat) %in% IFNG_genes, ]))
  IFNGsig$IFNG = apply(IFNGsig, 1, mean)
  IFNGsig$patient = rownames(IFNGsig)
  
  PD1_PDL1 = c("PDCD1", "CD274")
  PD1_PDL1sig = data.frame(t(mat[rownames(mat) %in% PD1_PDL1, ]))
  PD1_PDL1sig$patient = rownames(PD1_PDL1sig)
  
  hot = c(
    "CXCL9",
    "CXCL10",
    "CXCL11",
    "CXCR3",
    "CD3",
    "CD4",
    "CD8A",
    "CD8B",
    "CD274",
    "PDCD1",
    "CXCR4",
    "CCL5"
  )
  cold = c("CXCL1", "CXCL2", "CCL20")
  hot_sig = t(mat[rownames(mat) %in% hot, ])
  #hot_sig$patient = rownames(hot_sig)
  cold_gene = intersect(rownames(mat), cold)
  if(length(cold_gene) == 1){
    cold_sig = t(mat[cold_gene, ])
    rownames(cold_sig) = cold_gene
    cold_sig = t(cold_sig)
  }else{
    cold_sig = t(mat[rownames(mat) %in% cold, ])
  }
  # cold_sig = t(mat[rownames(mat) %in% cold, ])
  #cold_sig$patient = rownames(cold_sig)
  
  hot_sig_binary = sign(hot_sig)
  cold_sig_binary = -1 * sign(cold_sig)
  
  hot_cold_class = data.frame(cbind(hot_sig_binary, cold_sig_binary))
  # hot_cold_class$patient = rownames(hot_cold_class)
  hot_cold_class$TIP_signature = apply(hot_cold_class, 1, sum)
  hot_cold_class$patient = rownames(hot_cold_class)
  hot_cold_sig = cbind(hot_cold_class[c('patient', 'TIP_signature')],
                       hot_sig, cold_sig)
  hot_cold_sig$TIP_averexpr = apply(hot_sig, 1, sum) - apply(cold_sig, 1, sum)
  
  hot_cold_sig$TIP_class = "cold"
  hot_cold_sig[hot_cold_sig$TIP_signature > 0, ]$TIP_class = "hot"
  names(hot_cold_sig) <- c("patient","TIP_signature", 
                           paste0("TIP_",names(hot_cold_sig)[
                             -which(is.element(names(hot_cold_sig),c("patient",
                                                                     "TIP_signature",
                                                                     "TIP_averexpr", 
                                                                     "TIP_class")))
                           ]),
                           "TIP_averexpr", "TIP_class")
  
  TIDE_pre = merge(merge(merge(
    merge(TIDE_result, CD8sig, by = "patient"),
    IFNGsig, by = "patient"
  ), PD1_PDL1sig, by = "patient"),
  hot_cold_sig, by = "patient")
  
  # TIDE_pre$TIDEv2 = TIDE_pre$Exclusion
  # TIDE_pre[TIDE_pre$TIP_class == "hot", ]$TIDEv2 =
  # TIDE_pre[TIDE_pre$TIP_class == "hot", ]$Dysfunction
  
  if(!file.exists(save_path)){dir.create(save_path)}
  save_path = paste0(save_path,"TIDE_TIP_result/")
  if(!file.exists(save_path)){dir.create(save_path)}
  
  write.csv(TIDE_pre, paste0(save_path, "TIDE_TIP_result_", cancername,"_",sample_type, ".csv"),
            quote = F) 
  # return(TIDE_pre) 
}

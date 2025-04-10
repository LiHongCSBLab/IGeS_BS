
# Integrating immune signatures ----------------------------------------------
# Source:
# immune_landscape_result from Thorsson V, Gibbs D L, Brown S D, et al. 
# The immune landscape of cancer[J]. Immunity, 2018, 48(4): 812-830. e14.
# 
# Scores_160_Signatures.tsv was download from GDC
# Scores_160_Signatures contains immune cell infiltration score calculated via 
# Cibersort, other immune signatures from Yasin, Wolf, Attractors, ICR, 
# c7 atoms and Bindea.
#
# 
# TIDE and TIP  --------------------------------------------------------------
# TIDE_sig was downloaded from Jiang P, Gu S, Pan D, et al. Signatures of T 
# cell dysfunction and exclusion predict cancer immunotherapy response[J]. 
# Nature medicine, 2018, 24(10): 1550-1558.
# TIDE_sig includes Dysfunction signature, exclusion signature, MDSC signatures,
# CAF signatures and M2 signatures
# 
# TIDE and TIP local
# compute locally with TIDE Signature files from author (requested) and TIP 
# according to original paper
# TIP reference: Tumor immunological phenotype signature-based high-throughput
# screening for the discovery of combination immunotherapy compounds
#
#
# Immunophenoscore  ----------------------------------------------------------
# IPS_sig is downloaded from website TCIA based on the research Charoentong P,
# Finotello F, Angelova M, et al. Pan-cancer immunogenomic analyses reveal 
# genotype-immunophenotype relationships and predictors of response to 
# checkpoint blockade[J]. Cell reports, 2017, 18(1): 248-262.
# 
#    IPS_sig is downloaded from `https://tcia.at/` web.
#    It contains four scores:
#       IPS -> ips_ctla4_neg_pd1_neg
#       IPS-CTLA4-blocker -> ips_ctla4_pos_pd1_neg
#       IPS-PD1/PDL1/PDL2-blocker -> ips_ctla4_neg_pd1_pos
#       IPS-CTLA4- and PD1/PDL1/PDL2-blocker -> ips_ctla4_pos_pd1_pos
#  
#    IPS local:
#    recalculate with source code 
#    (Github: https://github.com/MayerC-imed/Immunophenogram)
# 
#
# TIL Estimations  -----------------------------------------------------------
# TIMER2.0: Infiltrated immune cell estimation (TIMER2.0 together with other 
# 5 kinds of estimation methods, including CIBERSORT, quanTIseq, xCell, 
# MCP-counter (or mMCP-counter for mouse) and EPIC methods) The result table 
# is downloaded from TIMER2.0 website and required to cite all 6 publications.
# 
# TIL_estimation: Due to the patients are not completely match with other 
# immune signatures, TILhas been recalculated with immunedeconv and other 
# packages, refered to function `immunecell_estiamtor` above.
# 
# ImmunCellAI: Infiltrated immune cell estimation 
# The result table is downloaded from ImmunCellAI website.
#
#
# TIGS: tumor immunogenicity score  ------------------------------------------
# downloaded from: 
# https://xsliulab.github.io/tumor-immunogenicity-score/#immunotherapy-datasets-analyses
# 
# 
# Immune Resistance Program  -------------------------------------------------
# OE and other 47 immune signatures (OE)
# All signatures were merged together via TCGA-SampleID (12 charaters).
# # Note: This signatures is computed from the RNAseq matrix below
#         and the value in the file below is: log2(norm_value+1)
# tcgaRNApath = "/picb/bigdata/project/FengFYM/UCEC_Xena/EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena"
#
# 
# cluster_signatures_cancercell  ---------------------------------------------
# Downloaded from https://github.com/BostonGene
# reference: Conserved pan-cancer microenvironment subtypes predict response to immunotherapy
#
#
#

immunesig_generator <- function(immunesig_path){
    
    library(stringr)
    library(dplyr)
  # Immune signatures ----------------------------------------------------------
  # immune landscape (Cell immunity, 2018) -------------------------------------
  
  if(!file.exists(paste0(immunesig_path, "immuneLandscape_sig160_recompute.txt"))){
    ## recomputed signatures from Wolf et al -------------------------------------
    sigPath <- "sig160_result/zscore_output/"
    a <- list.files(paste0(immunesig_path, sigPath))
    df <- lapply(as.list(a), function(x){
      f <- read.csv(paste0(immunesig_path, sigPath,x), row.names = 1)
      colnames(f) <- str_replace_all(colnames(f), pattern='\\.', replacement='-')
      f = as.data.frame(t(f))
      f$ID = row.names(f)
      f$SampleID <- substr(f$ID,1,12)
      f$type <- gsub(".csv","",x)
      return(f)
    })
    sig160_result <- do.call(rbind, df)
    sig160_result <- sig160_result[c(70,69,71, 1:68)]
    names(sig160_result) <- c("SampleID", "ID", "type", paste0(names(sig160_result)[-c(1:3)],"_sig160"))
    ## recomputed signatures from Senbabaoglu et al, Attractors and ICR ----------
    sigPath <- "sig160ssgsea_result/scaled_output/"
    a <- list.files(paste0(immunesig_path, sigPath))
    df <- lapply(as.list(a), function(x){
      load(paste0(immunesig_path, sigPath,x))
      f <- esOut_scaled
      f = as.data.frame(t(f))
      f$ID = row.names(f)
      f$SampleID <- substr(f$ID,1,12)
      return(f)
    })
    sig160ssgsea_result <- do.call(rbind, df)
    sig160ssgsea_result <- sig160ssgsea_result[c(42,41,1:40)]

    sig160_result_all = inner_join(sig160_result, sig160ssgsea_result)
    write.table(sig160_result_all, paste0(immunesig_path, "immuneLandscape_sig160_recompute.txt"),
                sep = "\t", quote = F, row.names = F)
  }

  sig160_result <- read.table(paste0(immunesig_path, "immuneLandscape_sig160_recompute.txt"),
                         sep = "\t", header = T )


  # consensus immune-CAF clustering (Cell Cancer Cell, 2021) -------------------------------------
  
  if(!file.exists(paste0(immunesig_path, "TCGA_caf_immuneCLS_recompute.txt"))){
    ## recomputed signatures from Wolf et al -------------------------------------
    sigPath <- "TCGA_caf_immuneCLS_result/scaled_output/"
    a <- list.files(paste0(immunesig_path, sigPath))
    df <- lapply(as.list(a), function(x){
      load(paste0(immunesig_path, sigPath,x))
      f <- esOut_scaled
      f = as.data.frame(t(f))
      f$ID = row.names(f)
      f$SampleID <- substr(f$ID,1,12)
      return(f)
    })
    TCGA_caf_immuneCLS_result <- do.call(rbind, df)
    TCGA_caf_immuneCLS_result <- TCGA_caf_immuneCLS_result[c(31,30,1:29)]
    write.table(TCGA_caf_immuneCLS_result, 
                paste0(immunesig_path, "TCGA_caf_immuneCLS_recompute.txt"),
                sep = "\t", quote = F, row.names = F)

  }

  TCGA_caf_immuneCLS_result <- read.table(paste0(immunesig_path, "TCGA_caf_immuneCLS_recompute.txt"),
                                          sep = "\t", header = T )


  # TIDE TCGA name after base --------------------------------------------------
  # TIDE downloaded ------------------------------------------------
  if(!file.exists(paste0(immunesig_path, "TIDE_TCGA_base_signature.txt"))){
    a <- list.files(paste0(immunesig_path, "TIDE/"))
    a <- a[stringr::str_detect(a, "TCGA")]
    a <- a[stringr::str_detect(a, "base")]
    dir <- paste0(immunesig_path,"TIDE/", a)
    f1 <- read.table(dir[1], sep="\t", header= T)
    f1$SampleID <- row.names(f1)
    for(i in 2:length(a)){
      f2 <- read.table(dir[i], sep="\t", header= T)
      f2$SampleID <- row.names(f2)
      f1 <- rbind(f1, f2)
    }
    f1 <- f1[c(6, 1:3)]
    names(f1) <- c("SampleID", paste0(colnames(f1)[-1], "_TIDE"))
    write.table(f1, paste0(immunesig_path, "TIDE_TCGA_base_signature.txt"),
                sep = "\t", quote = F, row.names = F)  
  }
  
  TIDE_sig <- read.table(paste0(immunesig_path, "TIDE_TCGA_base_signature.txt"),
                         sep = "\t", header = T )

  # TIDE and TIP - Calculate locally  ------------------------------------------
  if(!file.exists(paste0(immunesig_path, "TIDE_TIP_signature_TCGA.txt"))){
    a <- list.files(paste0(immunesig_path, "TIDE_TIP_result/"))
    df <- lapply(as.list(a), function(x){
      f <- read.csv(paste0(immunesig_path,"TIDE_TIP_result/", x), row.names = 1)
      f$ID <- f$patient
      f$SampleID <- substr(f$ID, 1, 12)
      f <- f[c( "SampleID", "ID",
                "Dysfunction", # "dys.score", 
                "Exclusion",# "exclusion_score", 
                "TIDE", 
                "CD8", 
                "IFNG", 
                "TIP_signature", 
                "PDCD1", 
                "CD274")]
      return(f)
    })
    TIDE_TIP_sig <- do.call(rbind, df)
    names(TIDE_TIP_sig) <- c("SampleID","ID",
                             "Dysfunction_TIDE","Exclusion_TIDE","TIDE_TIDE",
                             "CD8_TIDE","IFNG_TIDE",
                             "TIP_signature_TIP", 
                             "PDCD1_TIDE","CD274_TIDE")
    write.table(TIDE_TIP_sig, paste0(immunesig_path, "TIDE_TIP_signature_TCGA.txt"),
                sep = "\t", quote = F, row.names = F)
  }
  
  TIDE_TIP_sig <- read.table(paste0(immunesig_path, "TIDE_TIP_signature_TCGA.txt"),
                             sep = "\t", header = T)
  

  # Immunephenoscore -----------------------------------------------------------
  # Immunephenoscore - locally -------------------------------------------------
  if(!file.exists(paste0(immunesig_path, "IPS_signature_locally.txt"))){
    
    a <- list.files(paste0(immunesig_path, "IPS_result/"))
    df <- lapply(as.list(a), function(x){
      f <- read.csv(paste0(immunesig_path,"IPS_result/", x), row.names = 1)
      names(f) <- c("ID", "MHC", "EC", "SC","CP", "AZ", "IPS")
      f$SampleID <- substr(f$ID, 1, 12)
      return(f)
    })
    IPS_sig_locally = do.call(rbind, df)
    IPS_sig_locally = IPS_sig_locally[c(8,1:7)]
    names(IPS_sig_locally) <- c("SampleID", "ID", paste0(colnames(IPS_sig_locally)[-c(1,2)],"_IPS"))
    write.table(IPS_sig_locally, paste0(immunesig_path, "IPS_signature_locally.txt"),
                sep = "\t", quote = F, row.names = F)   
  }
  
  IPS_sig_locally <- read.table(paste0(immunesig_path, "IPS_signature_locally.txt"), sep="\t", header=T)
  

  # Immune Resistance Program  -------------------------------------------------
  # immuneOEsig_path <- "/picb/bigdata/project/FengFYM/immune_targeted_combo/immunecell_DEG_cor/sig_OE_results/"
  # # Note: This signatures is computed from the RNAseq matrix below
  #         and the value in the file below is: log2(norm_value+1)
  # tcgaRNApath = "/picb/bigdata/project/FengFYM/UCEC_Xena/EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena"
  
  if(!file.exists(paste0(immunesig_path, "sig_OE_signature.txt"))){
    
    a <- list.files(paste0(immunesig_path, "sig_OE_results/"))

    df <- lapply(as.list(a), function(x){
      load(paste0(immunesig_path, "sig_OE_results/", x))
      f <- df$res
      f <- as.data.frame(f)
      f$SampleID <- rownames(f)
      return(f)
    })
    sig_OE <- do.call(rbind, df)
    sig_OE <- sig_OE[c('SampleID', colnames(sig_OE)[-ncol(sig_OE)])]
    colnames(sig_OE) <- c('SampleID', paste0(colnames(sig_OE)[-1], "_OE"))
    write.table(sig_OE, paste0(immunesig_path, "sig_OE_signature.txt"), sep="\t", quote=F, row.names=F)
  }
  
  sig_OE <- read.csv(paste0(immunesig_path, "sig_OE_signature.txt"), sep="\t", header = T)
  

  # TIMER2.0 -------------------------------------------------------------
  TIMER <- read.csv(paste0(immunesig_path, "TIMER2.0_infiltration_estimation_for_tcga.csv"))
  colnames(TIMER) <- c("SampleID2",colnames(TIMER)[-1])
  TIMER$SampleID <- substr(TIMER$SampleID2,1,12)
  cibersort_TIMER <- cbind(TIMER[c("SampleID", "SampleID2")], TIMER[grep("CIBERSORT",colnames(TIMER))])
  cibersort_TIMER <- cibersort_TIMER[-grep("CIBERSORT.ABS",colnames(cibersort_TIMER))]
  
    

  # TIL estimation  ------------------------------------------------------------
  if(!file.exists(paste0(immunesig_path, "TILestimation_TCGA.txt"))){
    
    a <- list.files(paste0(immunesig_path, "TIL_estimation/merged/"))
    dir <- paste0(immunesig_path, "TIL_estimation/merged/", a)
    load(dir[1])
    f1 <- as.data.frame(t(res_TIL_all))
    f1$SampleID <- rownames(f1)
    # 7,31,34
    sp_idx = c(grep("DLBC", a), grep("THYM", a), grep("UVM",a))
    for(i in seq(length(a))[-sp_idx]){
      load(dir[i])
      f2 <- as.data.frame(t(res_TIL_all))
      f2$SampleID <- rownames(f2)
      
      f1 <- rbind(f1, f2)
    }
    
    load( paste0(immunesig_path, "TIL_estimation/merged/",a[grep("DLBC", a)]))
    f2 <- as.data.frame(t(res_TIL_all))
    missing_column <- setdiff(colnames(f1), colnames(f2))
    missing_column <- missing_column[-grep("Sample", missing_column)]
    f2_missing_column <- data.frame(matrix(NA, nrow(f2), length(missing_column)))
    colnames(f2_missing_column) <- missing_column
    f2 = cbind(f2, f2_missing_column)
    f2$SampleID <- rownames(f2)
    f2 = f2[colnames(f1)]
    f1 <- rbind(f1, f2)
    
    load( paste0(immunesig_path, "TIL_estimation/merged/",a[grep("THYM", a)]))
    f2 <- as.data.frame(t(res_TIL_all))
    missing_column <- setdiff(colnames(f1), colnames(f2))
    missing_column <- missing_column[-grep("Sample", missing_column)]
    f2_missing_column <- data.frame(matrix(NA, nrow(f2), length(missing_column)))
    colnames(f2_missing_column) <- missing_column
    f2 = cbind(f2, f2_missing_column)
    f2$SampleID <- rownames(f2)
    f2 = f2[colnames(f1)]
    f1 <- rbind(f1, f2)
    
    load( paste0(immunesig_path, "TIL_estimation/merged/",a[grep("UVM", a)]))
    f2 <- as.data.frame(t(res_TIL_all))
    missing_column <- setdiff(colnames(f1), colnames(f2))
    missing_column <- missing_column[-grep("Sample", missing_column)]
    f2_missing_column <- data.frame(matrix(NA, nrow(f2), length(missing_column)))
    colnames(f2_missing_column) <- missing_column
    f2 = cbind(f2, f2_missing_column)
    f2$SampleID <- rownames(f2)
    f2 = f2[colnames(f1)]
    f1 <- rbind(f1, f2)
    
    f1$ID <- f1$SampleID
    f1$SampleID <- substr(f1$ID, 1, 12)
    f1 <- f1[c(122,123, 1:121)]
    colnames(f1) <- gsub("timer.","",colnames(f1))
    colnames(f1) <- gsub("quantiseq.","",colnames(f1))
    colnames(f1) <- gsub("mcp_counter.","",colnames(f1))
    colnames(f1) <- gsub("epic.","",colnames(f1))
    colnames(f1) <- gsub("xcell.","",colnames(f1))

    write.table(f1, paste0(immunesig_path, "TILestimation_TCGA.txt"), sep="\t", quote=F, row.names=F)
  }
  
  TILestimation <- read.csv(paste0(immunesig_path, "TILestimation_TCGA.txt"), sep="\t", header = T)
  

  # ImmuCellAI -----------------------------------------------------------------
  
  if(!file.exists(paste0(immunesig_path, "ImmuCellAI_signature.txt"))){
    
    a <- list.files(paste0(immunesig_path, "ImmuCellAI/"))
    df <- lapply(as.list(a), function(x){
      f <- read.table(paste0(immunesig_path, "ImmuCellAI/",x), sep="\t", header= T)
      f$ID <- rownames(f)
      return(f)
    })
    ImmuCellAI <- do.call(rbind, df)
    ImmuCellAI$ID <- chartr(".","-",ImmuCellAI$ID)
    ImmuCellAI$SampleID <- substr(ImmuCellAI$ID,1,12)
    names(ImmuCellAI) <- c(paste0(names(ImmuCellAI)[-c(26,27)], "_ImmuCellAI"), "ID" ,"SampleID")
    ImmuCellAI <- ImmuCellAI[c("SampleID", "ID", names(ImmuCellAI)[-c(26,27)])]
    write.table(ImmuCellAI, paste0(immunesig_path, "ImmuCellAI_signature.txt"), sep="\t", quote=F, row.names=F)
  }
  
  ImmuCellAI <- read.csv(paste0(immunesig_path, "ImmuCellAI_signature.txt"), sep="\t", header = T)
  
  
    
  # Tumor immunogenicity score (TIGS)  -----------------------------------------
  # downloaded from https://xsliulab.github.io/tumor-immunogenicity-score/#exploration-of-aps-tmb-tigs-at-pan-cancer-level
  # reference: 
  TCGA_ALL_TIGS <- read.csv(paste0(immunesig_path, "TCGA_ALL_TIGS.csv"), 
                            sep=",", header = T)
  TCGA_ALL_TIGS <- TCGA_ALL_TIGS[c("Tumor_Sample_Barcode", "TIGS")]
  colnames(TCGA_ALL_TIGS) <- c("SampleID", "TIGS")
  
  
  # dectect loaded signatures --------------------------------------------------

  dim(sig160_result[!is.element(colnames(sig160_result), c("SampleID","SampleID2", "ID"))])
  dim(TCGA_caf_immuneCLS_result[!is.element(colnames(TCGA_caf_immuneCLS_result), c("SampleID","SampleID2", "ID"))])
  dim(TIDE_sig[!is.element(colnames(TIDE_sig), c("SampleID","SampleID2", "ID"))])
  dim(TIDE_TIP_sig[!is.element(colnames(TIDE_TIP_sig), c("SampleID","SampleID2", "ID"))])
  dim(IPS_sig_locally[!is.element(colnames(IPS_sig_locally), c("SampleID","SampleID2", "ID"))]) 
  dim(sig_OE[!is.element(colnames(sig_OE), c("SampleID","SampleID2", "ID"))])
  dim(cibersort_TIMER[!is.element(colnames(cibersort_TIMER), c("SampleID","SampleID2", "ID"))])
  dim(TILestimation[!is.element(colnames(TILestimation), c("SampleID","SampleID2", "ID"))])
  dim(ImmuCellAI[!is.element(colnames(ImmuCellAI), c("SampleID","SampleID2", "ID"))])
  dim(TCGA_ALL_TIGS[!is.element(colnames(TCGA_ALL_TIGS), c("SampleID","SampleID2", "ID"))])

  head(names(sig160_result[!is.element(colnames(sig160_result), c("SampleID","SampleID2", "ID"))]))
  head(names(TCGA_caf_immuneCLS_result[!is.element(colnames(TCGA_caf_immuneCLS_result), c("SampleID","SampleID2", "ID"))]))
  head(names(TIDE_sig[!is.element(colnames(TIDE_sig), c("SampleID","SampleID2", "ID"))]))
  head(names(TIDE_TIP_sig[!is.element(colnames(TIDE_TIP_sig), c("SampleID","SampleID2", "ID"))]))
  head(names(IPS_sig_locally[!is.element(colnames(IPS_sig_locally), c("SampleID","SampleID2", "ID"))]))
  head(names(sig_OE[!is.element(colnames(sig_OE), c("SampleID","SampleID2", "ID"))]))
  head(names(cibersort_TIMER[!is.element(colnames(cibersort_TIMER), c("SampleID","SampleID2", "ID"))]))
  head(names(TILestimation[!is.element(colnames(TILestimation), c("SampleID","SampleID2", "ID"))]))
  head(names(ImmuCellAI[!is.element(colnames(ImmuCellAI), c("SampleID","SampleID2", "ID"))]))
  head(names(TCGA_ALL_TIGS[!is.element(colnames(TCGA_ALL_TIGS), c("SampleID","SampleID2", "ID"))]))

  table(sig160_result$type)

  # drop patients with multiple samples
  sig160_result_dupID <- unique(sig160_result[duplicated(sig160_result$SampleID), ]$SampleID)
  length(sig160_result_dupID)
  sig160_result_filtered <- sig160_result[!is.element(sig160_result$SampleID, sig160_result_dupID), ]

  TCGA_caf_immuneCLS_result_dupID <- unique(TCGA_caf_immuneCLS_result[duplicated(TCGA_caf_immuneCLS_result$SampleID), ]$SampleID)
  length(TCGA_caf_immuneCLS_result_dupID)
  TCGA_caf_immuneCLS_result_filtered <- TCGA_caf_immuneCLS_result[!is.element(TCGA_caf_immuneCLS_result$SampleID, TCGA_caf_immuneCLS_result_dupID), ]

  TIDE_sig_dupID <- unique(TIDE_sig[duplicated(TIDE_sig$SampleID), ]$SampleID)
  length(TIDE_sig_dupID)
  TIDE_sig_filtered <- TIDE_sig[!is.element(TIDE_sig$SampleID, TIDE_sig_dupID), ]

  TIDE_TIP_sig_dupID <- unique(TIDE_TIP_sig[duplicated(TIDE_TIP_sig$SampleID), ]$SampleID)
  length(TIDE_TIP_sig_dupID)
  TIDE_TIP_sig_filtered <- TIDE_TIP_sig[!is.element(TIDE_TIP_sig$SampleID, TIDE_TIP_sig_dupID), ]

  IPS_sig_locally_dupID <- unique(IPS_sig_locally[duplicated(IPS_sig_locally$SampleID), ]$SampleID)
  length(IPS_sig_locally_dupID)
  IPS_sig_locally_filtered <- IPS_sig_locally[!is.element(IPS_sig_locally$SampleID, IPS_sig_locally_dupID), ]

  sig_OE_dupID <- unique(sig_OE[duplicated(sig_OE$SampleID), ]$SampleID)
  length(sig_OE_dupID)
  sig_OE_filtered <- sig_OE[!is.element(sig_OE$SampleID, sig_OE_dupID), ]

  cibersort_TIMER_dupID <- unique(cibersort_TIMER[duplicated(cibersort_TIMER$SampleID), ]$SampleID)
  length(cibersort_TIMER_dupID)
  cibersort_TIMER_filtered <- cibersort_TIMER[!is.element(cibersort_TIMER$SampleID, cibersort_TIMER_dupID), ]

  TILestimation_dupID <- unique(TILestimation[duplicated(TILestimation$SampleID), ]$SampleID)
  length(TILestimation_dupID)
  TILestimation_filtered <- TILestimation[!is.element(TILestimation$SampleID, TILestimation_dupID), ]

  ImmuCellAI_dupID <- unique(ImmuCellAI[duplicated(ImmuCellAI$SampleID), ]$SampleID)
  length(ImmuCellAI_dupID)
  ImmuCellAI_filtered <- ImmuCellAI[!is.element(ImmuCellAI$SampleID, ImmuCellAI_dupID), ]

  TCGA_ALL_TIGS_dupID <- unique(TCGA_ALL_TIGS[duplicated(TCGA_ALL_TIGS$SampleID), ]$SampleID)
  length(TCGA_ALL_TIGS_dupID)
  TCGA_ALL_TIGS_filtered <- TCGA_ALL_TIGS[!is.element(TCGA_ALL_TIGS$SampleID, TCGA_ALL_TIGS_dupID), ]

  res_immunesig <- full_join(sig160_result_filtered, TCGA_caf_immuneCLS_result_filtered)
  res_ICBpredictor <- left_join(TIDE_TIP_sig_filtered, TIDE_sig_filtered)
  res_ICBpredictor <- full_join(res_ICBpredictor, IPS_sig_locally_filtered)
  res_ICBpredictor <- left_join(res_ICBpredictor, sig_OE_filtered)
  res_ICBpredictor <- left_join(res_ICBpredictor, TCGA_ALL_TIGS_filtered)

  res_TILestimated <- left_join(TILestimation_filtered, cibersort_TIMER_filtered)
  res_TILestimated <- left_join(res_TILestimated, ImmuCellAI_filtered)

  res_immunesig_all <- full_join(res_immunesig, res_ICBpredictor)
  res_immunesig_all <- full_join(res_immunesig_all, res_TILestimated)

  # res_immunesig_all <- res_immunesig_all[-which(is.element(names(res_immunesig_all), "type"))]
  write.table(res_immunesig_all, paste0(immunesig_path, "immune_sigatures.txt"), sep="\t", quote=F, row.names=F)


  ## survival analysis: TCGA sample ID (digit) = 12 ----------------------------
  # sample_selected = read.csv("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/02_tumorGene_immuneSig_correlation/data/TCGA_patient_ID_tumor_type_selected.csv")
  # sample_selected <- sample_selected[c(1,2,5,6)]
  # names(sample_selected) <- c("SampleID", "TCGA_barcode","tumor_type_merged",'TCGA_cancer_type')
  # length(intersect(sample_selected$SampleID, res_immunesig_all$SampleID))
  
  
  return(res_immunesig_all)
}


# Integrating immune signatures ----------------------------------------------
# Source:
# immune_landscape_result from Thorsson V, Gibbs D L, Brown S D, et al. 
# The immune landscape of cancer[J]. Immunity, 2018, 48(4): 812-830. e14.
# 
# Scores_160_Signatures.tsv was download from GDC
# Scores_160_Signatures contains immune cell infiltration score calculated via 
# Cibersort, other immune signatures from Yasin, Wolf, Attractors, ICR, 
# c7 atoms and Bindea.
#
# 
# TIDE and TIP  --------------------------------------------------------------
# TIDE_sig was downloaded from Jiang P, Gu S, Pan D, et al. Signatures of T 
# cell dysfunction and exclusion predict cancer immunotherapy response[J]. 
# Nature medicine, 2018, 24(10): 1550-1558.
# TIDE_sig includes Dysfunction signature, exclusion signature, MDSC signatures,
# CAF signatures and M2 signatures
# 
# TIDE and TIP local
# compute locally with TIDE Signature files from author (requested) and TIP 
# according to original paper
# TIP reference: Tumor immunological phenotype signature-based high-throughput
# screening for the discovery of combination immunotherapy compounds
#
#
# Immunophenoscore  ----------------------------------------------------------
# IPS_sig is downloaded from website TCIA based on the research Charoentong P,
# Finotello F, Angelova M, et al. Pan-cancer immunogenomic analyses reveal 
# genotype-immunophenotype relationships and predictors of response to 
# checkpoint blockade[J]. Cell reports, 2017, 18(1): 248-262.
# 
#    IPS_sig is downloaded from `https://tcia.at/` web.
#    It contains four scores:
#       IPS -> ips_ctla4_neg_pd1_neg
#       IPS-CTLA4-blocker -> ips_ctla4_pos_pd1_neg
#       IPS-PD1/PDL1/PDL2-blocker -> ips_ctla4_neg_pd1_pos
#       IPS-CTLA4- and PD1/PDL1/PDL2-blocker -> ips_ctla4_pos_pd1_pos
#  
#    IPS local:
#    recalculate with source code 
#    (Github: https://github.com/MayerC-imed/Immunophenogram)
# 
#
# TIL Estimations  -----------------------------------------------------------
# TIMER2.0: Infiltrated immune cell estimation (TIMER2.0 together with other 
# 5 kinds of estimation methods, including CIBERSORT, quanTIseq, xCell, 
# MCP-counter (or mMCP-counter for mouse) and EPIC methods) The result table 
# is downloaded from TIMER2.0 website and required to cite all 6 publications.
# 
# TIL_estimation: Due to the patients are not completely match with other 
# immune signatures, TILhas been recalculated with immunedeconv and other 
# packages, refered to function `immunecell_estiamtor` above.
# 
# ImmunCellAI: Infiltrated immune cell estimation 
# The result table is downloaded from ImmunCellAI website.
#
#
# TIGS: tumor immunogenicity score  ------------------------------------------
# downloaded from: 
# https://xsliulab.github.io/tumor-immunogenicity-score/#immunotherapy-datasets-analyses
# 
# 
# Immune Resistance Program  -------------------------------------------------
# OE and other 47 immune signatures (OE)
# All signatures were merged together via TCGA-SampleID (12 charaters).
# # Note: This signatures is computed from the RNAseq matrix below
#         and the value in the file below is: log2(norm_value+1)
# tcgaRNApath = "/picb/bigdata/project/FengFYM/UCEC_Xena/EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena"
#
# 
# cluster_signatures_cancercell  ---------------------------------------------
# Downloaded from https://github.com/BostonGene
# reference: Conserved pan-cancer microenvironment subtypes predict response to immunotherapy
#
#
#

immunesig_generator_old <- function(immunesig_path){
    
  # immune landscape supplement ------------------------------------------------
  immune_landscape_result <- read.csv(paste0(immunesig_path, "immune_landscape_suppl.txt"), sep = "\t")
  immune_landscape_result <- immune_landscape_result[c(1,2,5:64)]
  colnames(immune_landscape_result) <- c("SampleID", paste0(colnames(immune_landscape_result)[-1], "_immunelandscape"))
  
  # remove same signatures from sig160 files.
  immune_landscape_result2 <- immune_landscape_result[- which(is.element(colnames(immune_landscape_result), 
                                                                         c("B.Cells.Memory_immunelandscape" , 
                                                                           "B.Cells.Naive_immunelandscape", 
                                                                           "Dendritic.Cells.Activated_immunelandscape", 
                                                                           "Dendritic.Cells.Resting_immunelandscape", 
                                                                           "Eosinophils_immunelandscape", 
                                                                           "Eosinophils.1_immunelandscape", 
                                                                           "Macrophages_immunelandscape", 
                                                                           "Macrophages.M0_immunelandscape" , 
                                                                           "Macrophages.M1_immunelandscape" , 
                                                                           "Macrophages.M2_immunelandscape" , 
                                                                           "Mast.Cells_immunelandscape" , 
                                                                           "Mast.Cells.Activated_immunelandscape" , 
                                                                           "Mast.Cells.Resting_immunelandscape" , 
                                                                           "Monocytes_immunelandscape", 
                                                                           "Neutrophils_immunelandscape", 
                                                                           "Neutrophils.1_immunelandscape", 
                                                                           "NK.Cells.Activated_immunelandscape" , 
                                                                           "NK.Cells.Resting_immunelandscape" , 
                                                                           "T.Cells.CD4.Memory.Activated_immunelandscape" , 
                                                                           "T.Cells.CD4.Memory.Resting_immunelandscape" , 
                                                                           "T.Cells.CD4.Naive_immunelandscape", 
                                                                           "T.Cells.CD8_immunelandscape", 
                                                                           "T.Cells.Follicular.Helper_immunelandscape", 
                                                                           "T.Cells.gamma.delta_immunelandscape", 
                                                                           "T.Cells.Regulatory.Tregs_immunelandscape" , 
                                                                           "Th1.Cells_immunelandscape", 
                                                                           "Th17.Cells_immunelandscape" , 
                                                                           "Th2.Cells_immunelandscape")))]
  
  # 160 signatures -------------------------------------------------------------
  sig160_result <- read.csv(paste0(immunesig_path, "Scores_160_Signatures.tsv"), sep = "\t")
  colnames(sig160_result) <- chartr(".", "-",colnames(sig160_result))
  sig160_result2 <- as.data.frame(t(sig160_result[-c(1,2)]))
  sig160_result2$SampleID <- row.names(sig160_result2)
  sig160_result2$SampleID <- substr(sig160_result2$SampleID,1,12)
  sig160_result2 <- sig160_result2[c(161,1:160)]
  colnames(sig160_result2) <- c("SampleID", paste0(colnames(sig160_result2)[-1], "_sig160"))
  sig160_result2$ID <- rownames(sig160_result2)
  
  # TIDE - TCGA name after base
  if(!file.exists(paste0(immunesig_path, "TIDE_TCGA_base_signature.txt"))){
    a <- list.files(paste0(immunesig_path, "TIDE/"))
    a <- a[stringr::str_detect(a, "TCGA")]
    a <- a[stringr::str_detect(a, "base")]
    dir <- paste0(immunesig_path,"TIDE/",a)
    f1 <- read.table(dir[1], sep="\t", header= T)
    f1$SampleID <- row.names(f1)
    for(i in 2:length(a)){
      f2 <- read.table(dir[i], sep="\t", header= T)
      f2$SampleID <- row.names(f2)
      f1 <- rbind(f1, f2)
    }
    f1 <- f1[c(6,1:3)]
    write.table(f1, paste0(immunesig_path, "TIDE_TCGA_base_signature.txt"),
                sep = "\t", quote = F, row.names = F
    )
    
  }
  
  TIDE_sig <- read.table(paste0(immunesig_path, "TIDE_TCGA_base_signature.txt"),
                         sep = "\t", header = T
  )
  names(TIDE_sig) <- c("SampleID", paste0(colnames(TIDE_sig)[-1], "_TIDE"))
  
  # TIDE and TIP - Calculate locally  ------------------------------------------
  if(!file.exists(paste0(immunesig_path, "TIDE_TIP_signature_TCGA.txt"))){
    a <- list.files(paste0(immunesig_path, "TIDE_TIP_result/"))
    dir <- paste0(immunesig_path,"TIDE_TIP_result/",a)
    f1 <- read.table(dir[1], sep=",", row.names = 1, header= T)
    f1$SampleID <- f1$patient
    f1 <- f1[c(
      "SampleID", "dys.score", "Dysfunction", "exclusion_score", "Exclusion",
      "TIDE", "CD8", "IFNG", "TIP_signature", "PDCD1", "CD274"
    )]
    for (i in 2:length(a)) {
      f2 <- read.table(dir[i], sep = ",", row.names = 1, header = T)
      f2$SampleID <- f2$patient
      f2 <- f2[c(
        "SampleID", "dys.score", "Dysfunction", "exclusion_score",
        "Exclusion", "TIDE", "CD8", "IFNG", "TIP_signature", "PDCD1", "CD274"
      )]
      f1 <- rbind(f1, f2)
    }
    names(f1) <- c(
      "SampleID", "Dysfunction_score", "Dysfunction", "Exclusion_score", 
      "Exclusion", "TIDE", "CD8", "IFNG", "TIP_signature", "PDCD1", "CD274"
    )
    write.table(f1, paste0(immunesig_path, "TIDE_TIP_signature_TCGA.txt"),
                sep = "\t", quote = F, row.names = F
    )
    
  }
  
  TIDE_TIP_sig <- read.table(paste0(immunesig_path, "TIDE_TIP_signature_TCGA.txt"),
                             sep = "\t", header = T
  )
  names(TIDE_TIP_sig) <- c("SampleID", 
                           paste0(colnames(TIDE_TIP_sig)[2:8], "_TIDE"),
                           paste0(colnames(TIDE_TIP_sig)[9], "_TIP"),
                           paste0(colnames(TIDE_TIP_sig)[10:11], "_TIDE"))
  
  
  # Immunephenoscore -----------------------------------------------------------
  if(!file.exists(paste0(immunesig_path, "IPS_signature.txt"))){
    
    a <- list.files(paste0(immunesig_path, "IPS/"))
    dir <- paste0(immunesig_path, "IPS/",a)
    f1 <- read.table(dir[1], sep = "\t", header = T)
    f1 <- f1[c(
      "barcode", "ips_ctla4_neg_pd1_neg", "ips_ctla4_neg_pd1_pos",
      "ips_ctla4_pos_pd1_neg", "ips_ctla4_pos_pd1_pos"
    )]
    
    for(i in 2:length(a)){
      f2 <- read.csv(dir[i], sep="\t", header= T)
      f2 <- f2[c(
        "barcode", "ips_ctla4_neg_pd1_neg", "ips_ctla4_neg_pd1_pos",
        "ips_ctla4_pos_pd1_neg", "ips_ctla4_pos_pd1_pos"
      )]
      f1 <- rbind(f1, f2)
    }
    
    colnames(f1) <- c("SampleID",colnames(f1)[-1])
    write.table(f1, paste0(immunesig_path, "IPS_signature.txt"),
                sep = "\t",
                quote = F, row.names = F
    )
    
  }
  
  IPS_sig <- read.table(paste0(immunesig_path, "IPS_signature.txt"),
                        sep = "\t", header = T
  )
  names(IPS_sig) <- c("SampleID", paste0(colnames(IPS_sig)[-1],"_IPS"))
  
  # Immunephenoscore - locally -------------------------------------------------
  if(!file.exists(paste0(immunesig_path, "IPS_signature_locally.txt"))){
    
    a <- list.files(paste0(immunesig_path, "IPS_result/"))
    dir <- paste0(immunesig_path, "IPS_result/",a)
    f1 <- read.table(dir[1], sep = ",", row.names = 1, header = T)
    
    for(i in 2:length(a)){
      f2 <- read.csv(dir[i], sep = ",", row.names = 1, header = T)
      f1 <- rbind(f1, f2)
    }
    
    names(f1) <- c("SampleID", "MHC", "EC", "SC","CP", "AZ", "IPS")
    write.table(f1, paste0(immunesig_path, "IPS_signature_locally.txt"),
                sep = "\t", quote = F, row.names = F
    )
    
  }
  
  IPS_sig_locally <- read.table(paste0(immunesig_path, "IPS_signature_locally.txt"), sep="\t", header=T)
  names(IPS_sig_locally) <- c("SampleID", paste0(colnames(IPS_sig_locally)[-1],"_IPS"))
  
  
  # Immune Resistance Program  -------------------------------------------------
  # immuneOEsig_path <- "/picb/bigdata/project/FengFYM/immune_targeted_combo/immunecell_DEG_cor/sig_OE_results/"
  # # Note: This signatures is computed from the RNAseq matrix below
  #         and the value in the file below is: log2(norm_value+1)
  # tcgaRNApath = "/picb/bigdata/project/FengFYM/UCEC_Xena/EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena"
  
  if(!file.exists(paste0(immunesig_path, "sig_OE_signature.txt"))){
    
    a <- list.files(paste0(immunesig_path, "sig_OE_results/"))
    dir <- paste0(immunesig_path, "sig_OE_results/", a)
    load(dir[1])
    f1 <- df$res
    #colnames(f1) <- paste0('OE_', f1)
    #f1$SampleID <- rownames(f1)
    
    for(i in 2:length(a)){
      load(dir[i])
      f2 <- df$res
      f1 <- rbind(f1, f2)
    }
    f1 <- as.data.frame(f1)
    f1$SampleID <- rownames(f1)
    f1 <- f1[c('SampleID', colnames(f1)[-ncol(f1)])]
    colnames(f1) <- c('SampleID', paste0(colnames(f1)[-1], "_OE"))
    
    write.table(f1, paste0(immunesig_path, "sig_OE_signature.txt"), sep="\t", quote=F, row.names=F)
  }
  
  sig_OE <- read.csv(paste0(immunesig_path, "sig_OE_signature.txt"), sep="\t", header = T)
  
  
  
  
  # TIMER2.0 -------------------------------------------------------------
  TIMER <- read.csv(paste0(immunesig_path, "TIMER2.0_infiltration_estimation_for_tcga.csv"))
  colnames(TIMER) <- c("SampleID",colnames(TIMER)[-1])
  cibersort_TIMER <- TIMER[c(1,grep("CIBERSORT",colnames(TIMER)))]
  
  
  
  
  # TIL estimation  ------------------------------------------------------------
  if(!file.exists(paste0(immunesig_path, "TILestimation_TCGA.txt"))){
    
    a <- list.files(paste0(immunesig_path, "TIL_estimation/merged/"))
    dir <- paste0(immunesig_path, "TIL_estimation/merged/", a)
    load(dir[1])
    f1 <- as.data.frame(t(res_TIL_all))
    f1$SampleID <- rownames(f1)
    # 7,31,34
    sp_idx = c(grep("DLBC", a), grep("THYM", a), grep("UVM",a))
    for(i in seq(length(a))[-sp_idx]){
      load(dir[i])
      f2 <- as.data.frame(t(res_TIL_all))
      f2$SampleID <- rownames(f2)
      
      f1 <- rbind(f1, f2)
    }
    
    load( paste0(immunesig_path, "TIL_estimation/merged/",a[grep("DLBC", a)]))
    f2 <- as.data.frame(t(res_TIL_all))
    missing_column <- setdiff(colnames(f1), colnames(f2))
    missing_column <- missing_column[-grep("Sample", missing_column)]
    f2_missing_column <- data.frame(matrix(NA, nrow(f2), length(missing_column)))
    colnames(f2_missing_column) <- missing_column
    f2 = cbind(f2, f2_missing_column)
    f2$SampleID <- rownames(f2)
    f2 = f2[colnames(f1)]
    f1 <- rbind(f1, f2)
    
    load( paste0(immunesig_path, "TIL_estimation/merged/",a[grep("THYM", a)]))
    f2 <- as.data.frame(t(res_TIL_all))
    missing_column <- setdiff(colnames(f1), colnames(f2))
    missing_column <- missing_column[-grep("Sample", missing_column)]
    f2_missing_column <- data.frame(matrix(NA, nrow(f2), length(missing_column)))
    colnames(f2_missing_column) <- missing_column
    f2 = cbind(f2, f2_missing_column)
    f2$SampleID <- rownames(f2)
    f2 = f2[colnames(f1)]
    f1 <- rbind(f1, f2)
    
    load( paste0(immunesig_path, "TIL_estimation/merged/",a[grep("UVM", a)]))
    f2 <- as.data.frame(t(res_TIL_all))
    missing_column <- setdiff(colnames(f1), colnames(f2))
    missing_column <- missing_column[-grep("Sample", missing_column)]
    f2_missing_column <- data.frame(matrix(NA, nrow(f2), length(missing_column)))
    colnames(f2_missing_column) <- missing_column
    f2 = cbind(f2, f2_missing_column)
    f2$SampleID <- rownames(f2)
    f2 = f2[colnames(f1)]
    f1 <- rbind(f1, f2)
    
    f1 <- f1[c("SampleID", colnames(f1)[-ncol(f1)])]
    
    write.table(f1, paste0(immunesig_path, "TILestimation_TCGA.txt"), sep="\t", quote=F, row.names=F)
  }
  
  TILestimation <- read.csv(paste0(immunesig_path, "TILestimation_TCGA.txt"), sep="\t", header = T)
  
  
  
  # ImmuCellAI -----------------------------------------------------------------
  
  if(!file.exists(paste0(immunesig_path, "ImmuCellAI_signature.txt"))){
    
    a <- list.files(paste0(immunesig_path, "ImmuCellAI/"))
    dir <- paste0(immunesig_path, "ImmuCellAI/",a)
    f1 <- read.table(dir[1], sep="\t", header= T)
    f1$ID <- rownames(f1)
    
    for(i in 2:length(a)){
      f2 <- read.csv(dir[i], sep="\t", header= T)
      f2$ID <- rownames(f2)
      f1 <- rbind(f1, f2)
    }
    
    f1$ID <- chartr(".","-",f1$ID)
    f1$SampleID <- substr(f1$ID,1,12)
    f1 <- f1[c("SampleID", names(f1)[-c(26,27)], "ID")]
    write.table(f1, paste0(immunesig_path, "ImmuCellAI_signature.txt"), sep="\t", quote=F, row.names=F)
  }
  
  ImmuCellAI <- read.csv(paste0(immunesig_path, "ImmuCellAI_signature.txt"), sep="\t", header = T)
  colnames(ImmuCellAI) <- c("SampleID", paste0(names(ImmuCellAI)[-which(names(ImmuCellAI) %in% c('SampleID', 'ID'))], "_ImmuCellAI"), "ID")
  
  
  
  # Tumor immunogenicity score (TIGS)  -----------------------------------------
  # downloaded from https://xsliulab.github.io/tumor-immunogenicity-score/#exploration-of-aps-tmb-tigs-at-pan-cancer-level
  # reference: 
  TCGA_ALL_TIGS <- read.csv(paste0(immunesig_path, "TCGA_ALL_TIGS.csv"), 
                            sep=",", header = T)
  TCGA_ALL_TIGS <- TCGA_ALL_TIGS[c("Tumor_Sample_Barcode", "TIGS")]
  colnames(TCGA_ALL_TIGS) <- c("SampleID", "TIGS")
  
  
  
  # cluster_signatures_cancercell ----------------------------------------------
  TCGA_caf_immune_sig <- read.table(paste0(
    immunesig_path,
    "cluster_signatures_cancercell/signatures-tcga.tsv"
  ), sep = "\t", header = T, row.names = 1)
  
  colnames(TCGA_caf_immune_sig) <- gsub("\\." ,"-", colnames(TCGA_caf_immune_sig))
  TCGA_caf_immune_sig <- as.data.frame(t(TCGA_caf_immune_sig))
  TCGA_caf_immune_sig$SampleID = rownames(TCGA_caf_immune_sig)
  TCGA_caf_immune_sig <- TCGA_caf_immune_sig[c("SampleID", names(TCGA_caf_immune_sig)[-which(names(TCGA_caf_immune_sig) == "SampleID")])]
  write.table(TCGA_caf_immune_sig, paste0(immunesig_path, "TCGA_caf_immune_sig.txt"), sep="\t", quote=F, row.names=F)
  
  TCGA_caf_immune_sig <- read.csv(paste0(immunesig_path, "TCGA_caf_immune_sig.txt"), sep="\t", header = T)
  colnames(TCGA_caf_immune_sig) <- c("SampleID", paste0(names(TCGA_caf_immune_sig)[-which(names(TCGA_caf_immune_sig) %in% c('SampleID', 'ID'))], "_TCGA_caf_immuneCLS"))
  
  # survival analysis: TCGA sample ID = 12
  sample_selected = read.csv("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/02_tumorGene_immuneSig_correlation/data/TCGA_patient_ID_tumor_type_selected.csv")
  sample_selected <- sample_selected[c(1,2,5,6)]
  names(sample_selected) <- c("SampleID", "TCGA_barcode","tumor_type_merged",'TCGA_cancer_type')
  
  




  rawDataPath <- "/picb/bigdata/project/FengFYM/s6_TCGAbiolinks_download_normalization/"
  cancerlist <- c("LGG","BRCA","CESC","CHOL","COAD", "ACC","BLCA",
                  "ESCA","GBM","HNSC","KICH","KIRC","KIRP","LIHC","LUAD","LUSC",
                  "DLBC","MESO","OV","PAAD","PCPG","PRAD","READ","SARC","SKCM",
                  "STAD","TGCT","THYM","THCA","UCS","UCEC","LAML") # "UVM", 
  
  cancer = cancerlist[1]
  sampleInfoPath = paste0(
    rawDataPath,
    "TCGA_sampleinfo/",
    cancer,
    ".Rdata"
  )
  load(sampleInfoPath)
  sampleinfo_cancer <- sampleinfo[c('patient','sample_type')]
  names(sampleinfo_cancer) <- c('SampleID','sample_type')
  sampleinfo_cancer$ID <- rownames(sampleinfo_cancer)
  
  for(cancer in cancerlist){
    sampleInfoPath = paste0(
      rawDataPath,
      "TCGA_sampleinfo/",
      cancer,
      ".Rdata"
    )
    load(sampleInfoPath)
    sampleinfo_cancer2 <- sampleinfo[c('patient','sample_type')]
    names(sampleinfo_cancer2) <- c('SampleID','sample_type')
    sampleinfo_cancer2$ID <- rownames(sampleinfo_cancer2)
    
    sampleinfo_cancer <- rbind(sampleinfo_cancer,sampleinfo_cancer2)
  }
  
  # sample_selected = merge(sample_selected, sampleinfo_cancer, by = c())
  
  head(TIDE_TIP_sig)  # SampleID 28; 
  head(IPS_sig_locally) # SampleID 28 
  head(TILestimation) # SampleID 28
  colnames(TILestimation) <- gsub("timer.","",colnames(TILestimation))
  colnames(TILestimation) <- gsub("quantiseq.","",colnames(TILestimation))
  colnames(TILestimation) <- gsub("mcp_counter.","",colnames(TILestimation))
  colnames(TILestimation) <- gsub("epic.","",colnames(TILestimation))
  colnames(TILestimation) <- gsub("xcell.","",colnames(TILestimation))
  immune_sigatures1 <- merge(TIDE_TIP_sig, IPS_sig_locally, by = "SampleID",all=T)
  immune_sigatures1 <- merge(immune_sigatures1, unique(TILestimation), by = "SampleID",all=T)
  immune_sigatures1$ID <- immune_sigatures1$SampleID
  immune_sigatures1$TCGA_barcode <- substr(immune_sigatures1$ID, 1, 15)
  # head(cibersort_TIMER) # SampleID 15 
  # immune_sigatures1 <- merge(immune_sigatures1, unique(cibersort_TIMER), by = "SampleID", all=T)
  immune_sigatures1$SampleID <- substr(immune_sigatures1$ID, 1, 12)
  
  head(sig160_result2)  # SampleID 12; ID 28
  head(ImmuCellAI)  # SampleID 12; ID 28 
  immune_sigatures2 <- merge(unique(sig160_result2), unique(ImmuCellAI), by = c("ID","SampleID"),all=T)
  
  
  head(immune_landscape_result2) # SampleID 12
  head(TIDE_sig)  # SampleID 12 
  head(IPS_sig) # SampleID 12 
  head(sig_OE)  # SampleID 12 
  head(TCGA_ALL_TIGS) # SampleID 12
  head(TCGA_caf_immune_sig) # SampleID 12 
  immune_sigatures3 <- merge(unique(immune_landscape_result2), unique(IPS_sig), by = "SampleID",all=T)
  immune_sigatures3 <- merge(immune_sigatures3, unique(sig_OE), by = "SampleID",all=T)
  immune_sigatures3 <- merge(immune_sigatures3, unique(TCGA_ALL_TIGS), by = "SampleID",all=T)
  immune_sigatures3 <- merge(immune_sigatures3, unique(TCGA_caf_immune_sig), by = "SampleID",all=T)
  TIDE_sig_duplicated = TIDE_sig[duplicated(TIDE_sig$SampleID), ]
  TIDE_sig2 = TIDE_sig[!is.element(TIDE_sig$SampleID, TIDE_sig_duplicated$SampleID), ]
  immune_sigatures3 <- merge(immune_sigatures3, unique(TIDE_sig2[1:4]), by = "SampleID",all=T)
  
  
  
  # select patient with same critria as TCGA data preprocessing -------------------------------
  immune_sigatures1_2 <- merge(sample_selected, immune_sigatures1, by = c("SampleID", "TCGA_barcode"))
  sig1_duplicated <- immune_sigatures1_2[is.element(immune_sigatures1_2$SampleID,
                                                    immune_sigatures1_2[duplicated(immune_sigatures1_2$SampleID),]$SampleID),]
  sig1_duplicated <- sig1_duplicated[order(sig1_duplicated$SampleID),]
  # table(sig1_duplicated$TCGA_cancer_type)
  # BLCA BRCA COAD  GBM KIRC LUAD LUSC PRAD SKCM UCEC 
  #  9   17   35    4   12   26    2    6    4   12
  immune_sigatures1_2 <- immune_sigatures1_2[!is.element(immune_sigatures1_2$SampleID, sig1_duplicated$SampleID),]
  
  sample_selected2 <- unique(immune_sigatures1_2[c("SampleID", "TCGA_barcode","ID","tumor_type_merged","TCGA_cancer_type")])
  immune_sigatures2_2 <- merge(sample_selected2, immune_sigatures2, by = c("ID","SampleID"))
  sig2_duplicated <- immune_sigatures2_2[is.element(immune_sigatures2_2$SampleID,
                                                    immune_sigatures2_2[duplicated(immune_sigatures2_2$SampleID),]$SampleID),]
  sig2_duplicated <- sig2_duplicated[order(sig2_duplicated$SampleID),]
  # table(sig2_duplicated$TCGA_cancer_type)
  # GBM LUSC SKCM 
  #  4    2    4 
  immune_sigatures2_2 <- immune_sigatures2_2[!is.element(immune_sigatures2_2$SampleID, sig2_duplicated$SampleID),]
  
  
  immune_sigatures3_2 <- merge(unique(sample_selected2), immune_sigatures3, by = "SampleID")
  sig3_duplicated <- immune_sigatures3_2[is.element(immune_sigatures3_2$SampleID,
                                                    immune_sigatures3_2[duplicated(immune_sigatures3_2$SampleID),]$SampleID),]
  sig3_duplicated <- sig3_duplicated[order(sig3_duplicated$SampleID),]
  # table(sig3_duplicated$TCGA_cancer_type)
  # BLCA BRCA COAD  GBM KIRC LUAD LUSC PRAD SKCM UCEC 
  #   9   17   35    4   12   26    2    6    4   12 
  immune_sigatures3_2 <- immune_sigatures3_2[!is.element(immune_sigatures2_2$SampleID, sig3_duplicated$SampleID),]
  
  
  head(cibersort_TIMER)
  immune_sigatures4 = cibersort_TIMER
  immune_sigatures4$TCGA_barcode = immune_sigatures4$SampleID
  immune_sigatures4$SampleID = substr(immune_sigatures4$SampleID,1,12)
  immune_sigatures4_2 <- merge(unique(sample_selected2), unique(immune_sigatures4), by = c("TCGA_barcode","SampleID"))
  sig4_duplicated <- immune_sigatures4_2[is.element(immune_sigatures4_2$SampleID,
                                                    immune_sigatures4_2[duplicated(immune_sigatures4_2$SampleID),]$SampleID),]
  sig4_duplicated <- sig4_duplicated[order(sig4_duplicated$SampleID),]
  # table(sig4_duplicated$TCGA_cancer_type)
  # BLCA BRCA COAD  GBM KIRC LUAD LUSC PRAD SKCM UCEC 
  #  9   15   35    4   12   26    2    6    4   12 
  immune_sigatures4_2 <- immune_sigatures4_2[!is.element(immune_sigatures4_2$SampleID, sig4_duplicated$SampleID),]
  
  
  # merge in patient level
  
  immune_sigatures <- merge(immune_sigatures1_2, immune_sigatures2_2, 
                            by = intersect(colnames(immune_sigatures1_2), colnames(immune_sigatures2_2)))
  immune_sigatures <- merge(immune_sigatures, immune_sigatures3_2, 
                            by = intersect(colnames(immune_sigatures), colnames(immune_sigatures3_2)))
  immune_sigatures <- merge(immune_sigatures, immune_sigatures4_2, 
                            by = intersect(colnames(immune_sigatures), colnames(immune_sigatures4_2)))
  
  #   table(immune_sigatures$TCGA_cancer_type)
  
  # ACC BLCA BRCA CESC CHOL COAD DLBC ESCA  GBM HNSC KICH KIRC KIRP  LGG LIHC LUAD 
  #  78  405 1068  303   36  427   48  159  151  495   64  509  283  509  368  498 
  # LUSC MESO   OV PAAD PCPG PRAD READ SARC SKCM STAD TGCT THCA THYM UCEC  UCS  UVM 
  # 486   86  262  176  177  489  156  255  466  371  149  498  119  522   56   80 
  
  immune_sigatures <- immune_sigatures[-which(is.element(colnames(immune_sigatures),
                                                         c("TCGA_barcode","tumor_type_merged", "TCGA_cancer_type")))]
  immune_sigatures <- immune_sigatures[c("SampleID", colnames(immune_sigatures)[-c(1,2)],"ID")]
  
  
  write.table(immune_sigatures, paste0(immunesig_path, "immune_sigatures.txt"), sep="\t", quote=F, row.names=F)
  
  return(immune_sigatures)
}


# using ITS to predict ICB prediction ------------------------------------------

immune_sigAnalysis <- function(data, 
                               gs,
                               method = "gsva",
                               kcdf = "Gaussian",
                               ssgsea.norm) {
  
  ## load R package
  if (!require(GSVA)) {
    BiocManager::install("GSVA")
  }
  
  ## please make sure that the rownames of data
  ## are the same with the gene names in the gene set
  ## that you want to run GSVA
  # rownames(data) <- .....
  
  ## gsva
  esOut <- gsva(data,
                gs = gs, 
                method=method,
                mx.diff=FALSE, kcdf=kcdf,
                verbose=FALSE, parallel.sz=1,
                min.sz=0, max.sz=2000,
                ssgsea.norm=ssgsea.norm)
  return(esOut)
}





immuneSig_ITSgsva_wilcoxTest <- function(ITSp_res,
                                         ITSn_res,
                                         response,
                                         #  cancer = "SKCM",
                                         #  tumor_subtype = "Metastatic",
                                         enrich_method = "gsva",
                                         # savepath_allresult = "02_tumorGene_immuneSig_correlation/results_selected/immunecell_TCGAgenes_result_spearman_TUMERIC/",
                                         savepath){
  
  #   cancer = "SKCM"
  #   tumor_subtype = "Metastatic"
  #   purity_method = "TUMERIC"
  #   datatype="allgenes"
  #   num_gene = 200
  #   workpath = "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immune_sig_correlation_analysis/"
  #   immune_sigGSpath = "02_tumorGene_immune_sig_correlation/data/immune_signatures/originalSignaturesGenes/"
  #   savepath = "r1_drug_immune_sig_correlation_analysis/"
  #     immune_sig_path = "06_results_summaryandplotting/data/immune_sig_selected/"
  
  #   immune_sig <- read.csv(paste0(immune_sig_path,"/", cancer,"_",sample_type,".csv"))
  #   immune_sig <- immune_sig[which(immune_sig$confident >= 0.5 ), ]
  
  dir.create(paste0(savepath, "/immuneSig_ITS_wilcox/"))
  savepath = paste0(savepath, "/immuneSig_ITS_wilcox/")
  
  ## Load ITS gene sets ----------------------------------------------------------
  #   ITSpath = "03_drug_immune_sig_enrichment/data/gene_immune_sig_TUMERIC_spearman/for_GSVA/pcor/"
  #   load(paste0(ITSpath, "spearman_",cancer,"_", tumor_subtype, "_positive_200.Rdata"))
  #   load(paste0(ITSpath, "spearman_",cancer,"_", tumor_subtype, "_negative_200.Rdata"))
  #   savepath_allresult = paste0("02_tumorGene_immuneSig_correlation/results_selected/immunecell_TCGAgenes_result_spearman_TUMERIC/")
  # load(paste0(savepath_allresult, "/pcor_spearman_", cancer, "_",tumor_subtype, "_positive_200.Rdata"))
  # load(paste0(savepath_allresult, "/pcor_spearman_", cancer, "_",tumor_subtype, "_negative_200.Rdata"))
  
  
  
  # ## gsva to generate sample ITS score -------------------------------------------
  # ITSp <- immune_sigAnalysis(data = as.matrix(log2(1+tpm_mat)), 
  #                            gs = genelist_p,
  #                            method = "gsva", 
  #                            kcdf = "Gaussian") 
  # ITSn <- immune_sigAnalysis(data = as.matrix(log2(1+tpm_mat)), 
  #                            gs = genelist_n,
  #                            method = "gsva", 
  #                            kcdf = "Gaussian") 
  
  
  ## wilcox test: R vs NR ---------------------------------------------------------
  ## ITS positive relate to original immune_sig
  
  ITSp_res_wilcox_RgreaterNR <- apply(ITSp_res[-which(is.element(names(ITSp_res), c("patient", "label", "patientID")))], 2, function(x){
    tmp = cbind(ITSp_res[c("patientID", "label")], data.frame(score = x))
    wilcox.test(tmp[tmp$label == "R", ]$score, tmp[tmp$label == "NR", ]$score, alternative = "greater")$p.value
  })
  
  ITSp_res_wilcox_RlessNR <- apply(ITSp_res[-which(is.element(names(ITSp_res), c("patient", "label", "patientID")))], 2, function(x){
    tmp = cbind(ITSp_res[c("patientID", "label")], data.frame(score = x))
    wilcox.test(tmp[tmp$label == "R", ]$score, tmp[tmp$label == "NR", ]$score, alternative = "less")$p.value
  })
  
  ITSp_res_wilcox_RgreaterNR = as.data.frame(ITSp_res_wilcox_RgreaterNR)
  ITSp_res_wilcox_RgreaterNR$immune_sig <- rownames(ITSp_res_wilcox_RgreaterNR)
  ITSp_res_wilcox_RlessNR = as.data.frame(ITSp_res_wilcox_RlessNR)
  ITSp_res_wilcox_RlessNR$immune_sig <- rownames(ITSp_res_wilcox_RlessNR)
  ITSp_res_wilcox = full_join(ITSp_res_wilcox_RgreaterNR, ITSp_res_wilcox_RlessNR)
  ITSp_res_wilcox = ITSp_res_wilcox[c("immune_sig", "ITSp_res_wilcox_RgreaterNR", "ITSp_res_wilcox_RlessNR")]
  # ITSp_res_wilcox = inner_join(immune_sig[c("immune_sig","type1")], ITSp_res_wilcox)
  # View(ITSp_res_wilcox)
  
  
  ## ITS negative relate to original immune_sig
  
  ITSn_res_wilcox_RgreaterNR <- apply(ITSn_res[-which(is.element(names(ITSp_res), c("patient", "label", "patientID")))], 2, function(x){
    tmp = cbind(ITSn_res[c("patientID", "label")], data.frame(score = x))
    wilcox.test(tmp[tmp$label == "R", ]$score, tmp[tmp$label == "NR", ]$score, alternative = "greater")$p.value
  })
  
  ITSn_res_wilcox_RlessNR <- apply(ITSn_res[-which(is.element(names(ITSp_res), c("patient", "label", "patientID")))], 2, function(x){
    tmp = cbind(ITSn_res[c("patientID", "label")], data.frame(score = x))
    wilcox.test(tmp[tmp$label == "R", ]$score, tmp[tmp$label == "NR", ]$score, alternative = "less")$p.value
  })
  
  ITSn_res_wilcox_RgreaterNR = as.data.frame(ITSn_res_wilcox_RgreaterNR)
  ITSn_res_wilcox_RgreaterNR$immune_sig <- rownames(ITSn_res_wilcox_RgreaterNR)
  ITSn_res_wilcox_RlessNR = as.data.frame(ITSn_res_wilcox_RlessNR)
  ITSn_res_wilcox_RlessNR$immune_sig <- rownames(ITSn_res_wilcox_RlessNR)
  ITSn_res_wilcox = full_join(ITSn_res_wilcox_RgreaterNR, ITSn_res_wilcox_RlessNR)
  ITSn_res_wilcox = ITSn_res_wilcox[c("immune_sig", "ITSn_res_wilcox_RgreaterNR", "ITSn_res_wilcox_RlessNR")]
  # ITSn_res_wilcox = inner_join(immune_sig[c("immune_sig","type1")], ITSn_res_wilcox)
  # View(ITSn_res_wilcox)
  
  
  ITS_res_wilcox = full_join(ITSp_res_wilcox, ITSn_res_wilcox)
  
  write.csv(ITS_res_wilcox, paste0(savepath, "ITS_res_wilcox_", enrich_method, ".csv"), 
            row.names =F, quote = F)
  return(ITS_res_wilcox)
  
}


## correlation dot plot ----------------------------------------------------------------------
cor_immunesig_ITSgsva_ICBdata <- function(data,
                                          genelist_p,
                                          genelist_n,
                                          ori_TMEsig,
                                          enrich_method,
                                          type = c("ITS","mergedITS"),
                                          plot = FALSE,
                                          savepath,
                                          savepath_dotplot=NULL){
  
  library(readxl)
  
  immuneSigInfoPath = "06_results_summaryandplotting/data/immene_cell_hieracy/"
  immunesigInfo <- as.data.frame(read_excel(paste0(immuneSigInfoPath, "immune_sig_hieracy_new.xlsx"), sheet = 1))
  immunesigInfo <- immunesigInfo[c(1,4,5:9)]
  immunesigInfo$immune_sig = tolower(immunesigInfo$immune_sig)
  immunesigInfo = immunesig_name_unify(immunesigInfo)
  
  xcell_cellfamily <- read.csv("02_tumorGene_immuneSig_correlation/data/immune_signatures/originalSignaturesGenes/immune_sig_original_files/XCELL/cells_families.txt", sep = '\t')
  xcell_type1 <- xcell_cellfamily[is.element(xcell_cellfamily$Group, 
                                             c("B-cells",
                                               "CD4+ memory T-cells",
                                               "CD4+ T-cells",
                                               "CD8+ T-cells",
                                               "DC",
                                               "Endothelial cells",
                                               "Epithelial cells",
                                               "Fibroblasts",
                                               "Macrophages")), ]
  
  xcell_type2 <- xcell_cellfamily[is.element(xcell_cellfamily$Group, 
                                             "Parent"), ]
  xcell_type3 <- xcell_type2[is.element(xcell_type2$Type, 
                                        c("Lymphoid", 
                                          "Myeloid")), ]
  xcell_type4 <- xcell_type2[is.element(xcell_type2$Cells, 
                                        c("Endothelial cells", 
                                          "Epithelial cells",
                                          "Fibroblasts")), ]
  xcell_select <- rbind(xcell_type1, xcell_type3, xcell_type4)
  xcell_remove <- xcell_cellfamily[!is.element(xcell_cellfamily$Cells, xcell_select$Cells),]
  xcell_remove <- paste0(xcell_remove$Cells, "_xcell")
  xcell_remove <- gsub(" ", ".", xcell_remove)
  
  
  ITSp <- immune_sigAnalysis(data = data, 
                             gs = genelist_p,
                             method = enrich_method, 
                             kcdf = "Gaussian") 
  
  ITSp_res <- as.data.frame(t(ITSp))
  ITSp_res$patientID <- rownames(ITSp_res)
  ITSp_res <- inner_join(response, ITSp_res, by = "patientID")
  write.csv(ITSp_res, paste0(savepath, "/ITSp_", enrich_method, ".csv"), 
            row.names =F, quote = F)
  
  
  ITSn <- immune_sigAnalysis(data = data, 
                             gs = genelist_n,
                             method = enrich_method, 
                             kcdf = "Gaussian") 
  ITSn_res <- as.data.frame(t(ITSn))
  ITSn_res$patientID <- rownames(ITSn_res)
  ITSn_res <- inner_join(response, ITSn_res)
  write.csv(ITSn_res, paste0(savepath, "/ITSn_", enrich_method, ".csv"), 
            row.names =F, quote = F)
  
  # CORRELATION DOT PLOT -----------------------------------------------
  
  # dir.create(paste0("05_immuneSig_ICBresponse/results/01_Hugo_dataset_outcome/immuneSig_ITS_correction"))
  # savepath = paste0("05_immuneSig_ICBresponse/results/01_Hugo_dataset_outcome/immuneSig_ITS_correction")
  
  # ITSgsva_p <- read.csv(paste0("05_immuneSig_ICBresponse/results/01_Hugo_dataset_outcome/ITSasPredictor/ITSp_gsva.csv"))
  # ITSgsva_n <- read.csv(paste0("05_immuneSig_ICBresponse/results/01_Hugo_dataset_outcome/ITSasPredictor/ITSn_gsva.csv"))
  ITSgsva_p = ITSp_res
  ITSgsva_n = ITSn_res
  if (length(which(is.element(names(ITSgsva_p), xcell_remove)))>0) ITSgsva_p <- ITSgsva_p[-which(is.element(names(ITSgsva_p), xcell_remove))]
  if (length(which(is.element(names(ITSgsva_n), xcell_remove)))>0) ITSgsva_n <- ITSgsva_n[-which(is.element(names(ITSgsva_n), xcell_remove))]
  
  ITSgsva_p_name = data.frame(name = names(ITSgsva_p), immune_sig = tolower(names(ITSgsva_p)))
  ITSgsva_p_name = immunesig_name_unify(ITSgsva_p_name)
  names(ITSgsva_p) = ITSgsva_p_name$immune_sig
  
  ITSgsva_n_name = data.frame(name = names(ITSgsva_n), immune_sig = tolower(names(ITSgsva_n)))
  ITSgsva_n_name = immunesig_name_unify(ITSgsva_n_name)
  names(ITSgsva_n) = ITSgsva_n_name$immune_sig
  
  # 1) correlation between ITS and TMEsig ----------------------------------------
  # ori_TMEsig <- read.csv(paste0("05_immuneSig_ICBresponse/results/01_Hugo_dataset_outcome/SKCM_Metastatic_immuneSig_merged.csv"))
  ori_TMEsig_name = data.frame(name = names(ori_TMEsig), immune_sig = tolower(names(ori_TMEsig)))
  ori_TMEsig_name = immunesig_name_unify(ori_TMEsig_name)
  names(ori_TMEsig) = ori_TMEsig_name$immune_sig
  
  if(type == "ITS"){
    
    ori_TMEsig_p = ori_TMEsig[intersect(names(ori_TMEsig), names(ITSgsva_p))]
    rownames(ori_TMEsig_p) = ori_TMEsig$patient
    ITSgsva_p2 = ITSgsva_p[intersect(names(ori_TMEsig), names(ITSgsva_p))]
    rownames(ITSgsva_p2) = ITSgsva_p$patientid
    ITSgsva_p2 = ITSgsva_p2[rownames(ori_TMEsig_p),]
    if(length(grep("NA", rownames(ITSgsva_p2)))){ITSgsva_p2 = ITSgsva_p2[-grep("NA", rownames(ITSgsva_p2)),]}
    ori_TMEsig_p = ori_TMEsig_p[rownames(ITSgsva_p2),]
    cor_p = cor(ori_TMEsig_p, ITSgsva_p2, method = "spearman")
    
    cor_p = data.frame(diag(cor_p))
    cor_p$immune_sig = rownames(cor_p)
    names(cor_p) = c('r', 'immune_sig')
    cor_p = cor_p[order(cor_p$r),]
    write.csv(cor_p, paste0(savepath, "/sig_", enrich_method, ".csv"), row.names = F, quote = F)
    
    if(plot){
      
      low_cor_sig = cor_p[which(cor_p$r < 0.5), ]
      
      p1 = list()
      for(x in seq(nrow(low_cor_sig))){
        signame = low_cor_sig[x,2]
        df = cbind(ITSgsva_p2[signame],
                   ori_TMEsig_p[signame])
        names(df) <- c("ITS", "immune_sig")
        df$patientID = rownames(df)
        df = inner_join(df, response)
        r = cor.test(df$ITS, df$immune_sig, method = "spearman")
        # print(r$estimate)}
        p1[[x]] = ggplot(data = df, mapping = aes(x = ITS, y = immune_sig, , colour = label)) + 
          geom_point() +
          # geom_abline(intercept = 0, slope = 1) +
          theme_bw() +
          ggtitle(paste0(signame,
                         "\n spearman's r = ", round(r$estimate,3), 
                         "\n P value = ", signif(r$p.value,3)))+
          xlab(paste0("ITS_", signame)) + 
          ylab(paste0("immune_sig_", signame))
      }
      
      
      pdf(paste0(savepath, "/low_cor_sig_", enrich_method, ".pdf"),
          width = 5,
          height = 5, 
          onefile = TRUE)
      for(x in seq(length(p1))){
        print(p1[[x]])
      }
      dev.off()
      
      
      high_cor_sig = cor_p[which(cor_p$r >= 0.5), ]
      
      p2 = list()
      for(x in seq(nrow(high_cor_sig))){
        
        signame = high_cor_sig[x,2]
        # signame="mdsc_tcga_caf_immunecls"
        # signame="cytotoxic.cells_sig160"
        # signame="fibroblasts_mcp_counter"
        df = cbind(ITSgsva_p2[signame],
                   ori_TMEsig_p[signame])
        names(df) <- c("ITS", "immune_sig")
        df$patientID = rownames(df)
        df = inner_join(df, response)
        r = cor.test(df$ITS, df$immune_sig, method = "spearman")
        # print(r$estimate)}
        p2[[x]] = ggplot(data = df, mapping = aes(x = ITS, y = immune_sig, , colour = label)) + 
          geom_point() +
          # geom_abline(intercept = 0, slope = 1) +
          theme_bw() +
          ggtitle(paste0(signame,
                         "\n spearman's r = ", round(r$estimate,3), 
                         "\n P value = ", signif(r$p.value,3)))+
          xlab(paste0("ITS_", signame)) + 
          ylab(paste0("immune_sig_", signame))
        
      }
      
      
      pdf(paste0(savepath, "/high_cor_sig_", enrich_method, ".pdf"),
          width = 5,
          height = 5, 
          onefile = TRUE)
      for(x in seq(length(p2))){
        print(p2[[x]])
      }
      dev.off() 
    }
    
  }else if(type == "mergedITS"){
    
    ori_TMEsig_p = ori_TMEsig
    rownames(ori_TMEsig_p) = ori_TMEsig$patient
    ITSgsva_p2 = ITSgsva_p
    rownames(ITSgsva_p2) = ITSgsva_p$patientid
    ITSgsva_p2 = ITSgsva_p2[rownames(ori_TMEsig_p),]
    if(length(grep("NA", rownames(ITSgsva_p2)))){ITSgsva_p2 = ITSgsva_p2[-grep("NA", rownames(ITSgsva_p2)),]}
    ori_TMEsig_p = ori_TMEsig_p[rownames(ITSgsva_p2),]
    cor_p = cor(ori_TMEsig_p[-1], ITSgsva_p2[-c(1,2)], method = "spearman")
    cor_p[upper.tri(cor_p)] = NA
    cor_p = data.frame(melt(cor_p))
    cor_p = cor_p[!is.na(cor_p$value), ]
    names(cor_p) = c('immune_sig',"merged_ITS",'r')
    cor_p = cor_p[order(cor_p$r), ]
    
    write.csv(cor_p, paste0(savepath, "/sig_", enrich_method, ".csv"), row.names = F, quote = F)
    
    if(plot){
      cor_p = inner_join(cor_p, immunesigInfo)
      low_cor_sig = cor_p[which(cor_p$r < 0.5), ]
      
      p1 = list()
      for(x in seq(nrow(low_cor_sig))){
        signame = low_cor_sig[x,2]
        df = cbind(ITSgsva_p2[signame],
                   ori_TMEsig_p[signame])
        names(df) <- c("ITS", "immune_sig")
        df$patientID = rownames(df)
        df = inner_join(df, response)
        r = cor.test(df$ITS, df$immune_sig, method = "spearman")
        # print(r$estimate)}
        p1[[x]] = ggplot(data = df, mapping = aes(x = ITS, y = immune_sig, , colour = label)) + 
          geom_point() +
          # geom_abline(intercept = 0, slope = 1) +
          theme_bw() +
          ggtitle(paste0(signame,
                         "\n spearman's r = ", round(r$estimate,3), 
                         "\n P value = ", signif(r$p.value,3)))+
          xlab(paste0("ITS_", signame)) + 
          ylab(paste0("immune_sig_", signame))
      }
      
      
      pdf(paste0(savepath, "/low_cor_sig_", enrich_method, ".pdf"),
          width = 5,
          height = 5, 
          onefile = TRUE)
      for(x in seq(length(p1))){
        print(p1[[x]])
      }
      dev.off()
      
      
      high_cor_sig = cor_p[which(cor_p$r >= 0.5), ]
      
      p2 = list()
      for(x in seq(nrow(high_cor_sig))){
        
        signame = high_cor_sig[x,2]
        # signame="mdsc_tcga_caf_immunecls"
        # signame="cytotoxic.cells_sig160"
        # signame="fibroblasts_mcp_counter"
        df = cbind(ITSgsva_p2[signame],
                   ori_TMEsig_p[signame])
        names(df) <- c("ITS", "immune_sig")
        df$patientID = rownames(df)
        df = inner_join(df, response)
        r = cor.test(df$ITS, df$immune_sig, method = "spearman")
        # print(r$estimate)}
        p2[[x]] = ggplot(data = df, mapping = aes(x = ITS, y = immune_sig, , colour = label)) + 
          geom_point() +
          # geom_abline(intercept = 0, slope = 1) +
          theme_bw() +
          ggtitle(paste0(signame,
                         "\n spearman's r = ", round(r$estimate,3), 
                         "\n P value = ", signif(r$p.value,3)))+
          xlab(paste0("ITS_", signame)) + 
          ylab(paste0("immune_sig_", signame))
        
      }
      
      
      pdf(paste0(savepath, "/high_cor_sig_", enrich_method, ".pdf"),
          width = 5,
          height = 5, 
          onefile = TRUE)
      for(x in seq(length(p2))){
        print(p2[[x]])
      }
      dev.off()     
    }  
  }
  
  return(list(ITSp_res, ITSn_res))
}
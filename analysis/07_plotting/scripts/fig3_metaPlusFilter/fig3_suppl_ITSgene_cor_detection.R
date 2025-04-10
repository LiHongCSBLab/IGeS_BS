# rm(list=ls())
# Mapping drug DEGs to correlation results
print("start computing")
work_path <- "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/"

setwd(work_path)
options(stringsAsFactors = F)


library(reshape2)
library(stringr)
library(ggplot2)
library(dplyr)
source("03_drug_immuneSig_enrichment/scripts/02_GSEA_geneset_preparation_function.R")
source("06_results_summaryandplotting/scripts/immunesig_selection_function.R")


sg_cor_detection <- function(ITS_set,
                             genelistpath,
                             genelistfile,
                             num_ITS = 50, 
                             cancer= "SKCM",
                             sample_type="Metastatic",
                             tumor_purity_method = "TUMERIC",
                             savepath){
  
  # genelist = df_stat3
  df_stat <- read.csv(genelistpath, row.names = 1)
  df_stat2 <- sort(colSums(sign(df_stat)))
  df_stat3 <- as.data.frame(df_stat2[df_stat2 > num_ITS])
  names(df_stat3) = "num_ITS"
  

  immune_sig <- rownames(df_stat)
  
  glist = rownames(df_stat3)
  glist <- gsub('\\.', '-', glist)

  
  ITS_p_set = ITS_set
  
  cancertype = paste0(cancer,"_", sample_type)
  genelist_p <- ITS_p_set[[cancertype]]
  
  corfilepath <- list.files(paste0("02_tumorGene_immuneSig_correlation/results/immunecell_TCGAgenes_result_spearman_", tumor_purity_method, "/" ))
  
  load(paste0("02_tumorGene_immuneSig_correlation/results/immunecell_TCGAgenes_result_spearman_", tumor_purity_method, "/", corfilepath[grep(cancertype, corfilepath)]))
  imname <- data.frame(immunesig = names(cor_result), immune_sig = tolower(names(cor_result)))
  imname <- immunesig_name_unify(imname)
  names(cor_result) <- imname$immune_sig

  # cancer = unlist(strsplit(cancertype, "_", fixed=TRUE))[1]
  # sample_type = unlist(strsplit(cancertype, "_", fixed=TRUE))[2]
  
  # dotplot (gene expression and Immunesig) ---------------------------------------------------
  immunesig_path <- "02_tumorGene_immuneSig_correlation/data/immune_signatures/"
  immunesig_file <- "02_tumorGene_immuneSig_correlation/data/immune_signatures/immune_sigatures.txt"
  
  patient_ID = read.csv("02_tumorGene_immuneSig_correlation/data/TCGA_patient_ID_tumor_type_selected.csv")
  
  # load in integrated immune signature file
  
  immune_sigatures <-  read.table(immunesig_file, sep = "\t", header = T)
  immune_sigatures_selected <- immune_sigatures[immune_sigatures$type == paste0(cancer, "_", sample_type),]
  immune_sigatures_selected <- immune_sigatures_selected[-which(is.element(names(immune_sigatures_selected), "type"))]
  immune_sigatures_selected$SampleID2 <- substr(immune_sigatures_selected$ID, 1, 15)
  
  imname <- data.frame(immunesig = names(immune_sigatures_selected), immune_sig = tolower(names(immune_sigatures_selected)))
  imname <- immunesig_name_unify(imname)
  names(immune_sigatures_selected) <- imname$immune_sig
  
  processedDataPath <- '02_tumorGene_immuneSig_correlation/data/TCGA_processedData/'
  processedDataPath <- paste0(processedDataPath, cancer,'/')
  TPMMatPath = paste0(processedDataPath, "mixture_file_TCGA_", cancer,"_",
                      sample_type, "_TPM.txt")
  mymatrix_filter <- read.csv(TPMMatPath, sep = '\t', row.names = 1)
  colnames(mymatrix_filter) <- gsub('\\.', '-',   colnames(mymatrix_filter))
  colnames(mymatrix_filter) <- substr(colnames(mymatrix_filter),1,15)
  
  
  for(i in seq(length(glist))) {
      # i= 1
      p1=list()      
      g = glist[i]

  for(j in seq(length(immune_sig))){
      sig = immune_sig[j]
      tmp = cor_result[[sig]]
      if(is.null(tmp)){
        r = NULL
      }else {
        r = tmp[which(rownames(tmp) == g),]
      }
      # g = 'IRF4','PDCD1LG2','IKZF3','LCK'
      gexpr = as.data.frame(t(data.frame(mymatrix_filter[rownames(mymatrix_filter) == g ,])))
      gexpr$sampleid2 = gsub('\\.', '-',   rownames(gexpr))
      df <- inner_join(gexpr, immune_sigatures_selected[c("sampleid2", sig)])
      
      names(df) <- c("gene","sample","immunesig")
      df$gene <- log2(df$gene +1 )
      #  r = cor.test(df[,1], df[,3], method = "spearman")
     
      if(is.null(r)){
      p1[[j]] <- ggplot(data = df, mapping = aes(x = gene, y = immunesig)) + 
        geom_point()+
        # geom_abline(intercept = 0, slope = 1) +
        theme_bw() +
        xlab(g) + ylab(sig)
     }else if(nrow(r) == 0){
      p1[[j]] <- ggplot(data = df, mapping = aes(x = gene, y = immunesig)) + 
        geom_point()+
        # geom_abline(intercept = 0, slope = 1) +
        theme_bw() +
        xlab(g) + ylab(sig)
     }else{
      # print(r$estimate)}
      p1[[j]] <- ggplot(data = df, mapping = aes(x = gene, y = immunesig)) + 
        geom_point()+
        # geom_abline(intercept = 0, slope = 1) +
        theme_bw()+
        ggtitle(paste0(sig, "- gene ", g,
                       "\n spearman's r = ", round(unlist(r$cor), 3), 
                       "\n P value = ", signif(unlist(r$p.adj),3)))+
        xlab(g) + ylab(sig)
      }
      
    }
    # savepath = "07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en"
    # genelistfile = "ITSp_vs_OG_oncoKB"
    dir.create(paste0(savepath,"/ITS_correlation_sgdotplot/"))
    dir.create(paste0(savepath,"/ITS_correlation_sgdotplot/", genelistfile, "/"))
    dir.create(paste0(savepath,"/ITS_correlation_sgdotplot/", genelistfile, "/",cancertype,"/"))

    pdf(paste0(savepath, "/ITS_correlation_sgdotplot/", genelistfile, "/", cancertype, "/", g, "_p.pdf"),
        width = 5, height = 5, onefile = TRUE)
    for(x in seq(length(p1))){
      print(p1[[x]])
    }
    dev.off()

    rm(p1)
   
  }
}


load("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITS_set.Rdata")

genelistfiles <- list.files("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/")
genelistfiles <- genelistfiles[grep(".csv", genelistfiles)]
genelistfiles <- genelistfiles[grep("ITSp", genelistfiles)]
genelistfiles <- gsub(".csv","",genelistfiles)
genelistfiles <- c(genelistfiles[c(7,8,14,15)], genelistfiles[-c(7,8,14,15)])


for(genelistfile in genelistfiles[-c(1,2)] ){
    # genelistfile = "ITSp_vs_OG_oncoKB"
    # genelistfile = genelistfiles[2]
    genelistpath = paste0("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/", genelistfile, ".csv")
    savepath = "07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en"
    dir.create(savepath)
    # sg_cor_detection(ITS_set = ITS_p_set,
    #                 genelistpath = genelistpath,
    #                 genelistfile = genelistfile,
    #                 num_ITS = 50, 
    #                 cancer= "Primary",
    #                 sample_type="Metastatic",
    #                 tumor_purity_method = "CPE",
    #                 savepath = savepath)

    sg_cor_detection(ITS_set = ITS_p_set,
                    genelistpath = genelistpath,
                    genelistfile = genelistfile,
                    num_ITS = 50, 
                    cancer= "SKCM",
                    sample_type="Metastatic",
                    tumor_purity_method = "TUMERIC",
                    savepath = savepath)

    cancerlist <- unique(gsub("spearman_","", 
                            gsub("_Metastatic_positive_200.Rdata","",
                                gsub("_Metastatic_negative_200.Rdata","",
                                        gsub("_Primary_positive_200.Rdata","",
                                            gsub("_Primary_negative_200.Rdata","",
                                                list.files("03_drug_immuneSig_enrichment/data/gene_immune_sig_TUMERIC_spearman/for_GSVA/pcor")))))))
    cancerlist <- cancerlist[-grep("SKCM", cancerlist)]
    sample_type = "Primary"
    for(cancer in cancerlist){ # [14:20]

        sg_cor_detection(ITS_set = ITS_p_set,
                        genelistpath = genelistpath,
                        genelistfile = genelistfile,
                        num_ITS = 50, 
                        cancer= cancer,
                        sample_type=sample_type,
                        tumor_purity_method = "TUMERIC",
                        savepath = savepath)
    }

}

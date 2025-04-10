# Tumor genes - infiltrated immune cell fraction
workpath <- "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/"
setwd(workpath)
library(parallel)
library(psych)
library(ppcor)
library(reshape2)
library(ggplot2)
library(dplyr)
options(stringsAsFactors = F)

source("06_results_summaryandplotting/scripts/immunesig_selection_function.R")

source("02_tumorGene_immuneSig_correlation/scripts/functions/02_immuneSig_generation_function.R")
source("02_tumorGene_immuneSig_correlation/scripts/functions/03_cor_gene_immuneSig_function.R")

source("05_immuneSig_ICBresponse/scripts/functions/predictor_ITS.R")

# tcgaRNApath <- "/picb/bigdata/project/FengFYM/immune_targeted_combo/TCGAbiolinks_download_normalization/FPKM/"
# tcgaRNApath <- "/picb/bigdata/project/FengFYM/UCEC_Xena/EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena"
immunesig_path <- "02_tumorGene_immuneSig_correlation/data/immune_signatures/"
immunesig_file <- "02_tumorGene_immuneSig_correlation/data/immune_signatures/immune_sigatures.txt"
patient_ID = read.csv("02_tumorGene_immuneSig_correlation/data/TCGA_patient_ID_tumor_type_selected.csv")

# load in integrated immune signature file
if(! file.exists(immunesig_file)){
  immune_sigatures <- immunesig_generator(immunesig_path)
}else{
  immune_sigatures <-  read.table(immunesig_file, sep = "\t", header = T)
}

immune_sig <- immune_sigatures


cor_immunesig_ITSgsva <- function(cancer = "SKCM",
                                  sample_type = "Metastatic",
                                  tumor_purity_method = 'TUMERIC'){
  
  
  dir.create("02_tumorGene_immuneSig_correlation/immuneSig_ITS_correction_ssgsea/")
  dir.create(paste0("02_tumorGene_immuneSig_correlation/immuneSig_ITS_correction_ssgsea/", cancer, "_", sample_type))
  savepath = paste0("02_tumorGene_immuneSig_correlation/immuneSig_ITS_correction_ssgsea/", cancer, "_", sample_type)
  
  
  protein_coding_genes <- read.table("03_drug_immuneSig_enrichment/data/protein_coding_genes_hg38.txt", sep = '\t', header = T)
  
  # sample_types <- read.table( paste0(workpath, "data/TCGA_cancertypes_subtypes.txt"), header=T, sep="\t")
  sample_types <- read.csv("02_tumorGene_immuneSig_correlation/data/TCGA_patient_ID_tumor_type_selected.csv")
  cancertypes <- unique(sample_types$TCGA_cancer_type)
  load("01_tumor_purity/data/processedData/TCGA_tumor_purity_merged.Rdata")
  tumor_purity <-purity_merged
  
  resultpath_TUMERIC2 <- "02_tumorGene_immuneSig_correlation/results/immunecell_TCGAgenes_result_spearman_TUMERIC/"
  dir.create(resultpath_TUMERIC2)
  
  immune_sig$SampleID2 <- substr(immune_sig$ID, 1, 15)
  
  tumor_purity = tumor_purity[-which(duplicated(tumor_purity$Sample_ID)),]
  
  # Prepare TCGA expression matrix --------------------------------------------
  
  processedDataPath <- '02_tumorGene_immuneSig_correlation/data/TCGA_processedData/'
  processedDataPath <- paste0(processedDataPath, cancer,'/')
  
  print(paste0("start calculation for ", cancer))
  
  
  
  
  TPMMatPath = paste0(processedDataPath, "mixture_file_TCGA_", cancer,"_",sample_type, "_TPM.txt")
  mymatrix_filter <- read.csv(TPMMatPath, sep = '\t', row.names = 1)
  colnames(mymatrix_filter) <- gsub('\\.', '-',   colnames(mymatrix_filter))
  colnames(mymatrix_filter) <- substr(colnames(mymatrix_filter),1,15)
  mymatrix_filter <- mymatrix_filter[is.element(row.names(mymatrix_filter), protein_coding_genes$gene_name),]
  
  # matching files for correlation
  immune_sig$SampleID2 <- substr(immune_sig$ID, 1, 15)
  c.mat <- immune_sig[immune_sig$SampleID2 %in% unique(colnames(mymatrix_filter)),]
  c.mat <- c.mat[order(c.mat$SampleID2),]
  row.names(c.mat) <- c.mat$SampleID2
  
  t.mat <- mymatrix_filter[colnames(mymatrix_filter) %in% c.mat$SampleID2,]
  t.mat <- mymatrix_filter[, c.mat$SampleID2]
  
  t.mat <- log2(t.mat + 1)
  
  tumor_purity_filter <- tumor_purity[is.element(tumor_purity$Sample_ID,  c.mat$SampleID2),]
  tumor_purity_filter2 <- unique(tumor_purity_filter[c('Sample_ID',tumor_purity_method)])
  names(tumor_purity_filter2) <- c('Sample_ID','purity')
  
  
  
  c.mat <- c.mat[- which(names( c.mat) %in% c("SampleID","SampleID2", "PatientID",
                                              "TCGA.Study_immunelandscape", 
                                              "immunelandscape_TCGA.Study",
                                              "ID"))]
  
  cmat_name <- data.frame(immune_sig = tolower(colnames(c.mat)), signame = colnames(c.mat))
  cmat_name = immunesig_name_unify(cmat_name)
  colnames(c.mat) = cmat_name$immune_sig
  
  # Compute ITS gsva score
  # savepath_allresult = "02_tumorGene_immuneSig_correlation/results_selected/immunecell_TCGAgenes_result_spearman_TUMERIC/"
  savepath_allresult = paste0("02_tumorGene_immuneSig_correlation/results_selected/ITS_oe_original_spearman_", tumor_purity_method, "_r0.4/for_GSVA/pcor/")
  load(paste0(savepath_allresult, "/pcor_spearman_", cancer, "_",sample_type, "_positive_200.Rdata"))
  load(paste0(savepath_allresult, "/pcor_spearman_", cancer, "_",sample_type, "_negative_200.Rdata"))
  
  
  
  # ssgsea to generate sample ITS score -------------------------------------------
  ITSp <- immune_sigAnalysis(data = as.matrix(t.mat), 
                             gs = genelist_p,
                             method = "ssgsea", 
                             kcdf = "Gaussian",
                             ssgsea.norm = TRUE)
  ITSp = data.frame(t(ITSp))
  write.csv(ITSp, paste0(savepath, "/ITSp_ssgsea.csv"), quote = F, row.names = T)
  
  # tmp_p = data.frame(diag(cor(c.mat[colnames(ITSp)], ITSp , method = "spearman")))
  # names(tmp_p) = 'r'
  # tmp_p$immune_sig = rownames(tmp_p)
  # tmp_p = tmp_p[order(tmp_p$r), ]
  
  # write.csv(tmp_p, paste0(savepath, "/cor_sig_ssgsea.csv"), row.names=F, quote = F)
  
  
  ITSp <- read.csv(paste0(savepath, "/ITSp_ssgsea.csv"))
  rownames(ITSp) = ITSp$X
  
  # should add p value
  print("calculating correlation")  
  
  
  cor_result <- lapply(as.list(names(c.mat)[-c(1:3)]), function(x){
    print(x)
    # x = "sc_ips"
    if(is.element(x, names(ITSp))){
    immuneSig = c.mat[x]
    names(immuneSig) = 'immuneSig'
    immuneSig$X = rownames(immuneSig)
    
    tmp = inner_join(immuneSig, ITSp[c('X', x)])
    tmp = tmp[!is.na(tmp$immuneSig),];
    names(tmp) = c("immuneSig","X","ITS")
    cor_result <- tryCatch({a = corr.test(x = tmp$immuneSig, 
                                          y = tmp$ITS,
                                          method = 'spearman');
                            data.frame(r = a$r, p.value = a$p)},
                            error = function(e){
                            return(data.frame(r = NA, p.value = NA,))})
    
    cor_result <- as.data.frame((cor_result));
    names(cor_result) <- c("cor","p.value");
    cor_result$p.adj <- p.adjust(cor_result$p.value, method="fdr")
    cor_result$oriSig = x
    cor_result$ITS = x
    return(cor_result)    
    }
  })  
  
  cor_result = do.call(rbind,cor_result)
  write.csv(cor_result, paste0(savepath, "/cor_p_sig_ssgsea.csv"), row.names=F, quote = F)
  
}

#   # 1) if it is influence by the number of genes in ITS
#   lapply(genelist_p, length)[tmp_p[which(tmp_p$r < 0.5), ]$immune_sig]
#   # PLOT

#   # 1) if THERE ARE MANY ZERO IN EITHER IMMUNESIG SCORE OR ITS
#   ## DOT PLOT
#   low_cor_sig = tmp_p[which(tmp_p$r < 0.5), ]

#   p1=list()
#   for(x in seq(nrow(low_cor_sig))){
#     signame = low_cor_sig[x,2]
#     df = cbind(ITSp[signame],
#                c.mat[signame])
#     names(df) <- c("ITS", "immune_sig")

#     r = cor.test(df$ITS, df$immune_sig, method = "spearman")
#     # print(r$estimate)}
#     p1[[x]] = ggplot(data = df, mapping = aes(x = ITS, y = immune_sig)) + 
#       geom_point() +
#       # geom_abline(intercept = 0, slope = 1) +
#       theme_bw() +
#       ggtitle(paste0(signame,
#                      "\n spearman's r = ", round(r$estimate,3), 
#                      "\n P value = ", signif(r$p.value,3)))+
#       xlab(paste0("ITS_", signame)) + 
#       ylab(paste0("immune_sig_", signame))

#   }


#   dir.create("02_tumorGene_immuneSig_correlation/immuneSig_ITS_correction/")
#   dir.create(paste0("02_tumorGene_immuneSig_correlation/immuneSig_ITS_correction/", cancer, "_", sample_type))
#   savepath = paste0("02_tumorGene_immuneSig_correlation/immuneSig_ITS_correction/", cancer, "_", sample_type)

#   # ggsave(paste0(savepath, "/GRANS_PCA_16704732_sig160.pdf"), p)
#   pdf(paste0(savepath, "/low_cor_sig_ssgsea.pdf"),
#       width = 5,
#       height = 5, 
#       onefile = TRUE)
#   for(x in seq(length(p1))){
#     print(p1[[x]])
#   }
#   dev.off()


#   high_cor_sig = tmp_p[which(tmp_p$r >= 0.5), ]

#   p2=list()
#   for(x in seq(nrow(high_cor_sig))){
#     signame = high_cor_sig[x,2]
#     df = cbind(ITSp[signame],
#                c.mat[signame])
#     names(df) <- c("ITS", "immune_sig")

#     r = cor.test(df$ITS, df$immune_sig, method = "spearman")
#     # print(r$estimate)}
#     p2[[x]] = ggplot(data = df, mapping = aes(x = ITS, y = immune_sig)) + 
#       geom_point() +
#       # geom_abline(intercept = 0, slope = 1) +
#       theme_bw() +
#       ggtitle(paste0(signame,
#                      "\n spearman's r = ", round(r$estimate,3), 
#                      "\n P value = ", signif(r$p.value,3)))+
#       xlab(paste0("ITS_", signame)) + 
#       ylab(paste0("immune_sig_", signame))

#   }

#   pdf(paste0(savepath, "/high_cor_sig_ssgsea.pdf"),
#       width = 5,
#       height = 5, 
#       onefile = TRUE)
#   for(x in seq(length(p2))){
#     print(p2[[x]])
#   }
#   dev.off()

# }


# mSKCM
cor_immunesig_ITSgsva(cancer = "SKCM",
                      sample_type = "Metastatic",
                      tumor_purity_method = 'TUMERIC')

cor_immunesig_ITSgsva(cancer = "SKCM",
                      sample_type = "Primary",
                      tumor_purity_method = 'CPE')


cancerlist <- c("BLCA", "LUAD", "STAD", "KIRC", "GBM",
                 "UCEC","LIHC", "PRAD", "PAAD", "BRCA", "READ","LUSC",
                 "OV",  "COAD", "DLBC", "UCS", "HNSC","CESC","KIRC", "CHOL", 
                 "ESCA", "PCPG", "ACC","LGG", "KICH", "SARC","TGCT","THYM", 
                 "MESO", "KIRP","THCA")
for(cancer in cancerlist){
  cor_immunesig_ITSgsva(cancer = cancer,
                        sample_type = "Primary",
                        tumor_purity_method = 'TUMERIC')
  
}
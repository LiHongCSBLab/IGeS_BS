# rm(list=ls())
work_path <- "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/"

setwd(work_path)
options(stringsAsFactors = F)

library(reshape2)
library(ggplot2)
library(ggplotify)
library(grid)
library(cowplot)
library(parallel)
library(patchwork)
library(dplyr)
library(stringr)
library(data.table)
require(GGally)
library(plot3D)
library(fmsb)
library(readxl)
library(dplyr)

source("06_results_summaryandplotting/scripts/plot_for_drug_function.R")
source("06_results_summaryandplotting/scripts/ITS_analysis/ITS_gene_stat_function.R")
source("06_results_summaryandplotting/scripts/immunesig_selection_function.R")
source("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/scripts/drugITGSscore_lincs.R")


# ----------------------------------------------------------------------------

dir.create("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/")
# ----------------------------------------------------------------------------
load("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITS_set.Rdata")
# ------------------------------------------------------------------------------

ITS_GSenrichTest <- function(ITS_p_set, 
                                genelist,
                                filename){

    dir.create("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/enrichTest/")
    dir.create(paste0("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/enrichTest/", filename))

  for(cancertype in names(ITS_p_set)){
    # cancertype = "SKCM_Metastatic"
    ITSgene = ITS_p_set[[cancertype]]
    
    # filename = "ITSp_vs_OG_oncoKB"
    genelistpath = paste0("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/", filename, ".csv")
    df_stat <- read.csv(genelistpath, row.names = 1)
    df_stat <- df_stat # * 22
    df_stat_col <- rev(sort(colSums(sign(df_stat))))
    df_stat_col <- as.data.frame(df_stat_col)
    names(df_stat_col) = "num_ITS"

    df_stat_col2 <- rev(sort(colSums(df_stat)))
    df_stat_col2 <- as.data.frame(df_stat_col2)
    names(df_stat_col2) = "num_ITS_weighted"

    df_stat_row <- rev(sort(rowSums(sign(df_stat))))
    df_stat_row <- as.data.frame(df_stat_row)
    names(df_stat_row) = "num_Gene"

    df_stat_row2 <- rev(sort(rowSums(df_stat)))
    df_stat_row2 <- as.data.frame(df_stat_row2)
    names(df_stat_row2) = "num_Gene_weighted"

    write.csv(df_stat_row, paste0("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/enrichTest/", filename, "/num_Gene.csv"), quote = F)
    write.csv(df_stat_row2, paste0("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/enrichTest/", filename, "/num_Gene_weighted.csv"), quote = F)
    write.csv(df_stat_col, paste0("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/enrichTest/", filename, "/num_ITS.csv"), quote = F)
    write.csv(df_stat_col2, paste0("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/enrichTest/", filename, "/num_ITS_weighted.csv"), quote = F)
    # canceroncogene <- read.csv("06_results_summaryandplotting/data/canceroncegenes/oncoKB_cancerGeneList.tsv", sep = "\t", header = T) 
    # oncogene <- canceroncogene[canceroncogene$Is.Oncogene == "Yes",]$Hugo.Symbol
    # genelist <- oncogene
    # drugtarget <- read.csv('/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/s3_lincs_drug_MoA/drugTargetMerge.txt', sep = '\t')
    # genelist = drugtarget$target
    geneExprInCCLE <- read.csv("03_drug_immuneSig_enrichment/data/CCLE_geneexpr/geneExprInCCLE.csv")
    geneExprInCCLE_filted <- unique(geneExprInCCLE[geneExprInCCLE$rate < 0.2,]$gene_name)

    genelist = intersect(genelist, geneExprInCCLE_filted)

    ITS_enrichTest <- lapply(ITSgene, function(x){
      
      # x=ITSgene[[1]]

      #           ITS    non-ITS
      # GS          a          b    a+b  =  M  
      # non-GS      c          d    c+d
      #           a+c        b+d    N
      #             n      = N-n

      a = length(intersect(x, genelist))
      b = length(setdiff(genelist, x))
      c = length(setdiff(x, genelist))
      d = length(setdiff(geneExprInCCLE_filted, union(x, genelist)))

      N = length(geneExprInCCLE_filted)
      M = length(intersect(genelist, geneExprInCCLE_filted))
      dat = matrix(c(a,c,b,d), ncol=2,nrow=2)

      FishTestP = fisher.test(dat)$p.value
      ChisqTestP = chisq.test(dat)$p.value
      # pvalue=1-phyper(a,# 差异基因中，位于通路中基因数量
      #                 a+c, # 差异基因的数量
      #                 N-(a+c), # 全部基因的数量 - 差异数量
      #                 M)  # 全部基因中，位于通路中基因数量


    # N是数据总库的基因数，
    # M是目标基因集中的基因数，
    # n是获得的差异基因数量，
    # k是在目标基因集中富集的差异基因数
     k = a
     M = length(intersect(genelist, geneExprInCCLE_filted)) # a+b
     n = length(x) # a+c
     N = length(geneExprInCCLE_filted)
     pvalue=phyper(k-1,M, N-M, n, lower.tail=FALSE)                 
    # tmp <- data.frame(gene.not.interest = c(M-k, N-M-n+k), 
    #                   gene.in.interest  = c(k, n-k))
    # row.names(tmp) <- c("In_category", "not_in_category")
      data.frame(FishTestP = FishTestP, 
                 ChisqTestP = ChisqTestP,
                 HypergeometricP = pvalue)
  })

    ITS_enrichTest <- do.call(rbind, ITS_enrichTest)
    ITS_enrichTest$FishTestFDR = p.adjust(ITS_enrichTest$FishTestP)
    ITS_enrichTest$ChisqTestFDR = p.adjust(ITS_enrichTest$ChisqTestP)
    ITS_enrichTest$HypergeometricFDR = p.adjust(ITS_enrichTest$HypergeometricP)
    # head(ITS_enrichTest[c(1,3)],20)
    write.csv(ITS_enrichTest, paste0("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/enrichTest/", filename, "/", cancertype, ".csv"), quote = F)
    }


  #   Gene_enrichTest <- lapply(as.list(colnames(df_stat)), function(x){
      
  #     # x=colnames(df_stat)[1]
  #     df_stat <- df_stat * 22
  #     df_stat_col2 <- rev(sort(colSums(df_stat)))
  #     df_stat_col2 <- as.data.frame(df_stat_col2)
  #     names(df_stat_col2) = "num_ITS_weighted"


  #     #           g in OG in N selected ITS    g in OG not in selected ITS
  #     # g in GS          a                         b
  #     # g in non-GS      c                         d

  #     a = df_stat_col2[rownames(df_stat_col2) == x,]
  #     b = 22*130 - df_stat_col2[rownames(df_stat_col2) == x,]
  #     c = 22*130 - length(genelist[-grep(x, genelist)])
  #     d = length(setdiff(geneExprInCCLE_filted, union(x, genelist)))
  #     N = length(geneExprInCCLE_filted)
  #     M = length(intersect(genelist, geneExprInCCLE_filted))
  #     dat = matrix((c(a,c,b,d)), ncol=2,nrow=2)

  #     FishTestP = fisher.test(dat)$p.value
  #     ChisqTestP = chisq.test(dat)$p.value
  #     pvalue=1-phyper(a,# 差异基因中，位于通路中基因数量
  #                     a+c, # 差异基因的数量
  #                     N-(a+c), # 全部基因的数量 - 差异数量
  #                     M)  # 全部基因中，位于通路中基因数量

  #     data.frame(FishTestP = FishTestP, 
  #                ChisqTestP = ChisqTestP,
  #                HypergeometricP = pvalue)
  # })

  #   Gene_enrichTest <- do.call(rbind, Gene_enrichTest)
  #   Gene_enrichTest$FishTestFDR = p.adjust(Gene_enrichTest$FishTestP)
  #   Gene_enrichTest$ChisqTestFDR = p.adjust(Gene_enrichTest$ChisqTestP)
  #   Gene_enrichTest$HypergeometricFDR = p.adjust(Gene_enrichTest$HypergeometricP)
    
  #   dir.create("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/enrichTest/")
  #   dir.create(paste0("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/enrichTest/", filename))
  #   write.csv(Gene_enrichTest, paste0("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/enrichTest/", filename, "/", cancertype, ".csv"), quote = F)
    
}


# ------------------------------------------------------------------------------
# ITS gene vs drugtarget -------------------------------------------------------

drugtarget <- read.csv('/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/s3_lincs_drug_MoA/drugTargetMerge.txt', sep = '\t')

ITS_GSenrichTest(ITS_p_set = ITS_p_set, 
                 genelist = unique(drugtarget$target),
                 filename = "ITSp_vs_drugtarget")




# ------------------------------------------------------------------------------
# ITS gene vs druggable genes --------------------------------------------------
# ------------------------------------------------------------------------------

drugtarget_OASIS <- read.csv('/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/s3_lincs_drug_MoA/drugtarget_OASIS_TIDE.csv')

ITS_GSenrichTest(ITS_p_set = ITS_p_set, 
                   genelist = drugtarget_OASIS$gene_name, 
                   filename = "ITSp_vs_druggable_gene")


# ------------------------------------------------------------------------------
# ITS gene vs inhibitory and stimulatory genes ---------------------------------
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# inhibitory
immuneModulator_PMID29628290_1 <- as.data.frame(read_excel("06_results_summaryandplotting/data/immune_modulator/PMID29628290.xlsx", sheet = 1))
inhibitory <- unique(immuneModulator_PMID29628290_1[immuneModulator_PMID29628290_1$`Immune Checkpoint` == "Inhibitory",]$`HGNC Symbol`)
inhibitory <- inhibitory[!is.na(inhibitory)]

ITS_GSenrichTest(ITS_p_set = ITS_p_set, 
                   genelist = inhibitory, 
                   filename = "ITSp_vs_inhibitor_PMID29628290")


# ------------------------------------------------------------------------------
# stimulatory
immuneModulator_PMID29628290_1 <- as.data.frame(read_excel("06_results_summaryandplotting/data/immune_modulator/PMID29628290.xlsx", sheet = 1))
stimulatory <- unique(immuneModulator_PMID29628290_1[immuneModulator_PMID29628290_1$`Immune Checkpoint` != "Inhibitory",]$`HGNC Symbol`)
stimulatory <- stimulatory[!is.na(stimulatory)]


ITS_GSenrichTest(ITS_p_set = ITS_p_set, 
                   genelist = stimulatory, 
                   filename = "ITSp_vs_stimulatory_PMID29628290")



# ------------------------------------------------------------------------------
# co_stim_inhib 

immuneModulator_PMID29628290_2 <- as.data.frame(read_excel("06_results_summaryandplotting/data/immune_modulator/PMID29628290.xlsx", sheet = 2))
co_stim_inhib <- unique(immuneModulator_PMID29628290_2[immuneModulator_PMID29628290_2$Reason == "co-stim/inhib list", ]$Gene)

ITS_GSenrichTest(ITS_p_set = ITS_p_set, 
                   genelist = co_stim_inhib, 
                   filename = "ITSp_vs_co_stim_inhib_PMID29628290")


# ------------------------------------------------------------------------------
# ITS gene vs icb_enhancer and icb_suppressor genes ----------------------------
# ------------------------------------------------------------------------------

immuneModulator_PMID34980132_1 <- as.data.frame(read_excel("06_results_summaryandplotting/data/immune_modulator/PMID34980132/ICB enhancer and suppressor genes.xlsx", sheet = 1))
# ------------------------------------------------------------------------------
icb_enhancer <- unique(immuneModulator_PMID34980132_1[immuneModulator_PMID34980132_1$Type == "ICB enhancer gene", ]$Symbol)

ITS_GSenrichTest(ITS_p_set = ITS_p_set, 
                   genelist = icb_enhancer, 
                   filename = "ITSp_vs_icb_enhancer_PMID34980132")


# ------------------------------------------------------------------------------
immuneModulator_PMID34980132_1 <- as.data.frame(read_excel("06_results_summaryandplotting/data/immune_modulator/PMID34980132/ICB enhancer and suppressor genes.xlsx", sheet = 1))
# ------------------------------------------------------------------------------
icb_suppressor <- unique(immuneModulator_PMID34980132_1[immuneModulator_PMID34980132_1$Type == "ICB suppressors gene", ]$Symbol)

ITS_GSenrichTest(ITS_p_set = ITS_p_set, 
                   genelist = icb_suppressor, 
                   filename = "ITSp_vs_icb_suppressor_PMID34980132")


# ------------------------------------------------------------------------------
# ITS gene vs regulator genes --------------------------------------------------
# ------------------------------------------------------------------------------

immuneModulator_PMID34980132_2 <- as.data.frame(read_excel("06_results_summaryandplotting/data/immune_modulator/PMID34980132/positive and negative regulators.xlsx", sheet = 1))
# ------------------------------------------------------------------------------
# regulator_n
regulator_p <- unique(immuneModulator_PMID34980132_2[immuneModulator_PMID34980132_2$Type == "Positive regulators", ]$Symbol)

ITS_GSenrichTest(ITS_p_set = ITS_p_set, 
                   genelist = regulator_p, 
                   filename = "ITSp_vs_regulator_p_PMID34980132")

# ------------------------------------------------------------------------------
# regulator_n
immuneModulator_PMID34980132_2 <- as.data.frame(read_excel("06_results_summaryandplotting/data/immune_modulator/PMID34980132/positive and negative regulators.xlsx", sheet = 1))
regulator_n <- unique(immuneModulator_PMID34980132_2[immuneModulator_PMID34980132_2$Type == "Negative regulators", ]$Symbol)


ITS_GSenrichTest(ITS_p_set = ITS_p_set, 
                   genelist = regulator_n, 
                   filename = "ITSp_vs_regulator_n_PMID34980132")

# ------------------------------------------------------------------------------
# ITS gene vs resistor and sensitizor  -----------------------------------------
# ------------------------------------------------------------------------------

immuneModulator_PMID34980132_3 <- as.data.frame(read_excel("06_results_summaryandplotting/data/immune_modulator/PMID34980132/sensitizer_resistor_genes.xlsx", sheet = 1))
# ------------------------------------------------------------------------------
# Sensitizer
Sensitizer <- unique(immuneModulator_PMID34980132_3[immuneModulator_PMID34980132_3$Type == "Sensitizer", ]$Symbol)

ITS_GSenrichTest(ITS_p_set = ITS_p_set, 
                   genelist = Sensitizer, 
                   filename = "ITSp_vs_Sensitizer_PMID34980132")

# ------------------------------------------------------------------------------
# Resistor
Resistor <- unique(immuneModulator_PMID34980132_3[immuneModulator_PMID34980132_3$Type == "Resistor", ]$Symbol)

ITS_GSenrichTest(ITS_p_set = ITS_p_set, 
                   genelist = Resistor, 
                   filename = "ITSp_vs_Resistor_PMID34980132")

# ------------------------------------------------------------------------------
# ITS gene vs OG and TSG -------------------------------------------------------
# ------------------------------------------------------------------------------

immuneModulator_PMID34980132_4 <- as.data.frame(read_excel("06_results_summaryandplotting/data/immune_modulator/PMID34980132/Tumor immunity related OGs TSGs.xlsx", sheet = 1))
OG <- unique(immuneModulator_PMID34980132_4[immuneModulator_PMID34980132_4$Type == "Tumor immunity-related OGs", ]$Symbol)


ITS_GSenrichTest(ITS_p_set = ITS_p_set, 
                   genelist = OG, 
                   filename = "ITSp_vs_OG_PMID34980132")

TSG <- unique(immuneModulator_PMID34980132_4[immuneModulator_PMID34980132_4$Type == "Tumor immunity-related TSGs", ]$Symbol)

ITS_GSenrichTest(ITS_p_set = ITS_p_set, 
                   genelist = TSG, 
                   filename = "ITSp_vs_TSG_PMID34980132")

# ------------------------------------------------------------------------------
# cancer oncogenes -------------------------------------------------------------
# ------------------------------------------------------------------------------
canceroncogene <- read.csv("06_results_summaryandplotting/data/canceroncegenes/oncoKB_cancerGeneList.tsv", sep = "\t", header = T)
oncogene <- canceroncogene[canceroncogene$Is.Oncogene == "Yes",]$Hugo.Symbol
ITS_GSenrichTest(ITS_p_set = ITS_p_set, 
                   genelist = oncogene, 
                   filename = "ITSp_vs_OG_oncoKB")

filename = "ITSp_vs_OG_oncoKB"                   
resultpath <- "07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/enrichTest/"
resultfiles <- list.files(paste0(resultpath, filename))
resultfiles <- resultfiles[-grep("num_Gene_weighted.csv", resultfiles)]
resultfiles <- resultfiles[-grep("num_Gene.csv", resultfiles)]
resultfiles <- resultfiles[-grep("num_ITS_weighted.csv", resultfiles)]
resultfiles <- resultfiles[-grep("num_ITS.csv", resultfiles)]

res <- lapply(as.list(resultfiles), function(file){
  df <- read.csv(paste0(resultpath, filename,'/', file))
  df$cancertype <- gsub(".csv", "", file)
  return(df)
})
res <- do.call(rbind, res)
res_sig <- (rowSums(table(res[res$FishTestFDR < 0.05, ][c('X', 'cancertype')])))
res_sig[res_sig>10]
write.csv(res, paste0(resultpath, filename, ".csv"), quote = F)

suppressorgene <- canceroncogene[canceroncogene$Is.Tumor.Suppressor.Gene == "Yes",]$Hugo.Symbol
ITS_GSenrichTest(ITS_p_set = ITS_p_set, 
                   genelist = suppressorgene, 
                   filename = "ITSp_vs_TSG_oncoKB")

filename = "ITSp_vs_TSG_oncoKB"                   
resultpath <- "07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/enrichTest/"
resultfiles <- list.files(paste0(resultpath, filename))
resultfiles <- resultfiles[-grep("num_Gene_weighted.csv", resultfiles)]
resultfiles <- resultfiles[-grep("num_Gene.csv", resultfiles)]
resultfiles <- resultfiles[-grep("num_ITS_weighted.csv", resultfiles)]
resultfiles <- resultfiles[-grep("num_ITS.csv", resultfiles)]

res <- lapply(as.list(resultfiles), function(file){
  df <- read.csv(paste0(resultpath, filename,'/', file))
  df$cancertype <- gsub(".csv", "", file)
  return(df)
})

res <- do.call(rbind, res)
res_sig <- (rowSums(table(res[res$FishTestFDR < 0.05, ][c('X', 'cancertype')])))
res_sig[res_sig>10]
write.csv(res, paste0(resultpath, filename, ".csv"), quote = F)

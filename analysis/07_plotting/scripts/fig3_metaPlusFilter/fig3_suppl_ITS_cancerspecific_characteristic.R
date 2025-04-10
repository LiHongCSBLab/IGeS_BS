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

# source("06_results_summaryandplotting/scripts/plot_for_drug_function.R")
# source("06_results_summaryandplotting/scripts/ITS_analysis/ITS_gene_stat_function.R")
# source("06_results_summaryandplotting/scripts/immunesig_selection_function.R")
# source("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/scripts/drugITGSscore_lincs.R")

load("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITS_set.Rdata")

# ------------------------------------------------------------------------------
# cancer specific ITS genes characteristics
# ------------------------------------------------------------------------------
ITS_cancer_specific <- lapply(as.list(names(ITS_p_set)), function(cancer){

  # cancer = 'SKCM_Metastatic'
  ITS_cancer = ITS_p_set[[cancer]]
  ITS_other = ITS_p_set[-grep(cancer, names(ITS_p_set))]

  ITS_p_cancer <- lapply(as.list(names(ITS_cancer)), function(x){
    # x = names(ITS_cancer)[1]
    gs1 = ITS_cancer[[x]]
    gs2 = lapply(ITS_other, function(y){y[[x]]})
    gs2 = unique(unlist(gs2))
    setdiff(gs1, gs2)
  })
  names(ITS_p_cancer) <- names(ITS_cancer)
  return(ITS_p_cancer)
})
names(ITS_cancer_specific) <- names(ITS_p_set)

tmp = lapply(ITS_cancer_specific, function(x){unique(unlist(x))})
ITS_cancer_specific <- lapply(as.list(names(tmp)), function(cancer){
    # x = names(tmp)[1]
    ITS_cancer = tmp[[cancer]]
    ITS_other = tmp[-grep(cancer, names(tmp))]
    gs1 = unique(unlist(ITS_cancer))
    gs2 = unique(unlist(ITS_other))
    setdiff(gs1, gs2)
})
names(ITS_cancer_specific) <- names(ITS_p_set)

dir.create(("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITS_cancerspecific_gene"))
save( ITS_cancer_specific, file = "07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITS_cancerspecific_gene/ITS_cancer_specific.Rdata")

# ------------------------------------------------------------------------------

ITS_cancerspecific_characteristic <- function(ITS_p_cancer_specific, 
                                                genelist, 
                                                filename){
  library(UpSetR)
  library(scales)
  library(gridExtra)

  ITS_p_dg <- intersect(ITS_p_cancer_specific, genelist)
  ITS_p_dg2 <- sapply(ITS_p_dg, function(x)paste(x, collapse = "\t"))
  
  ITS_p_dg2 <- cbind(matrix(NA, length(ITS_p_dg2)), ITS_p_dg2)
  write.table(ITS_p_dg2, paste0("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/", filename, ".txt"), 
                sep='\t', row.names=T, col.names = F, quote=F)

  glist = unique(unlist(ITS_p_dg))
  listinput_df <- fromList(ITS_p_dg)
  rownames(listinput_df) <- glist
  genestat <- data.frame(apply(listinput_df,1,sum))
  names(genestat) <- "num_ITS"
  genestat$gene = rownames(genestat)
  genestat = genestat[order(genestat$num_ITS, decreasing = T), ]
  genestat$gene = factor(genestat$gene, levels = genestat$gene)
  write.csv(genestat, paste0("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/", filename, "_genestat.csv"), 
            row.names=F, quote=F)

  # genestat_plot = genestat[genestat$num_ITS > 0.5 * length(ITS_p_dg), ]
  genestat_plot = genestat[1:20,]
  p = ggplot(genestat_plot, aes(x=gene, y=num_ITS)) +
    geom_bar(stat="identity") +
    xlab("")+
    theme(axis.text.x = element_text(size = 10,angle = 90))+
    ggtitle(filename)
  ggsave( paste0("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/", filename, "_genestat.pdf"), p, width = 5, height = 4)

}

# ------------------------------------------------------------------------------


lapply(as.list(names(ITS_cancer_specific)), function(cancer){
  gset <- ITS_cancer_specific[[cancer]]
  dir.create(paste0("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITS_cancerspecific_gene/", cancer ))

  # ------------------------------------------------------------------------------
  # ITS gene vs drugtarget -------------------------------------------------------
  # ------------------------------------------------------------------------------

  drugtarget <- read.csv('/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/s3_lincs_drug_MoA/drugTargetMerge.txt', sep = '\t')
  # drugtarget <- data.frame(gene_name = unique(drugtarget$target))
    try({
    ITS_cancerspecific_characteristic(ITS_p_cancer_specific = gset, 
                    genelist = unique(drugtarget$target), 
                    filename = paste0("ITS_cancerspecific_gene/",cancer,"/ITSp_vs_drugtarget"))
    }, silent = T)


  # ------------------------------------------------------------------------------
  # ITS gene vs druggable genes --------------------------------------------------
  # ------------------------------------------------------------------------------

  drugtarget_OASIS <- read.csv('/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/s3_lincs_drug_MoA/drugtarget_OASIS_TIDE.csv')

    try({
  ITS_cancerspecific_characteristic(ITS_p_cancer_specific = gset, 
                    genelist = drugtarget_OASIS$gene_name, 
                    filename = paste0("ITS_cancerspecific_gene/",cancer,"/ITSp_vs_druggable_gene"))
    }, silent = T)


  # ------------------------------------------------------------------------------
  # ITS gene vs inhibitory and stimulatory genes ---------------------------------
  # ------------------------------------------------------------------------------

  # ------------------------------------------------------------------------------
  # inhibitory
  immuneModulator_PMID29628290_1 <- as.data.frame(read_excel("06_results_summaryandplotting/data/immune_modulator/PMID29628290.xlsx", sheet = 1))
  inhibitory <- unique(immuneModulator_PMID29628290_1[immuneModulator_PMID29628290_1$`Immune Checkpoint` == "Inhibitory",]$`HGNC Symbol`)
  inhibitory <- inhibitory[!is.na(inhibitory)]

    try({
  ITS_cancerspecific_characteristic(ITS_p_cancer_specific = gset, 
                    genelist = inhibitory, 
                    filename = paste0("ITS_cancerspecific_gene/",cancer,"/ITSp_vs_inhibitor_PMID29628290"))
    }, silent = T)

  # ------------------------------------------------------------------------------
  # stimulatory
  immuneModulator_PMID29628290_1 <- as.data.frame(read_excel("06_results_summaryandplotting/data/immune_modulator/PMID29628290.xlsx", sheet = 1))
  stimulatory <- unique(immuneModulator_PMID29628290_1[immuneModulator_PMID29628290_1$`Immune Checkpoint` != "Inhibitory",]$`HGNC Symbol`)
  stimulatory <- stimulatory[!is.na(stimulatory)]

    try({
  ITS_cancerspecific_characteristic(ITS_p_cancer_specific = gset, 
                    genelist = stimulatory, 
                    filename = paste0("ITS_cancerspecific_gene/",cancer,"/ITSp_vs_stimulatory_PMID29628290"))
    }, silent = T)



  # ------------------------------------------------------------------------------
  # co_stim_inhib 

  immuneModulator_PMID29628290_2 <- as.data.frame(read_excel("06_results_summaryandplotting/data/immune_modulator/PMID29628290.xlsx", sheet = 2))
  co_stim_inhib <- unique(immuneModulator_PMID29628290_2[immuneModulator_PMID29628290_2$Reason == "co-stim/inhib list", ]$Gene)
    try({
  ITS_cancerspecific_characteristic(ITS_p_cancer_specific = gset, 
                    genelist = co_stim_inhib, 
                    filename = paste0("ITS_cancerspecific_gene/",cancer,"/ITSp_vs_co_stim_inhib_PMID29628290"))
    }, silent = T)


  # ------------------------------------------------------------------------------
  # ITS gene vs icb_enhancer and icb_suppressor genes ----------------------------
  # ------------------------------------------------------------------------------

  immuneModulator_PMID34980132_1 <- as.data.frame(read_excel("06_results_summaryandplotting/data/immune_modulator/PMID34980132/ICB enhancer and suppressor genes.xlsx", sheet = 1))
  # ------------------------------------------------------------------------------
  icb_enhancer <- unique(immuneModulator_PMID34980132_1[immuneModulator_PMID34980132_1$Type == "ICB enhancer gene", ]$Symbol)
    try({
  ITS_cancerspecific_characteristic(ITS_p_cancer_specific = gset, 
                    genelist = icb_enhancer, 
                    filename = paste0("ITS_cancerspecific_gene/",cancer,"/ITSp_vs_icb_enhancer_PMID34980132"))
    }, silent = T)


  # ------------------------------------------------------------------------------
  immuneModulator_PMID34980132_1 <- as.data.frame(read_excel("06_results_summaryandplotting/data/immune_modulator/PMID34980132/ICB enhancer and suppressor genes.xlsx", sheet = 1))
  # ------------------------------------------------------------------------------
  icb_suppressor <- unique(immuneModulator_PMID34980132_1[immuneModulator_PMID34980132_1$Type == "ICB suppressors gene", ]$Symbol)
    try({
  ITS_cancerspecific_characteristic(ITS_p_cancer_specific = gset, 
                    genelist = icb_suppressor, 
                    filename = paste0("ITS_cancerspecific_gene/",cancer,"/ITSp_vs_icb_suppressor_PMID34980132"))
    }, silent = T)


  # ------------------------------------------------------------------------------
  # ITS gene vs regulator genes --------------------------------------------------
  # ------------------------------------------------------------------------------

  immuneModulator_PMID34980132_2 <- as.data.frame(read_excel("06_results_summaryandplotting/data/immune_modulator/PMID34980132/positive and negative regulators.xlsx", sheet = 1))
  # ------------------------------------------------------------------------------
  # regulator_n
  regulator_p <- unique(immuneModulator_PMID34980132_2[immuneModulator_PMID34980132_2$Type == "Positive regulators", ]$Symbol)
    try({
  ITS_cancerspecific_characteristic(ITS_p_cancer_specific = gset, 
                    genelist = regulator_p, 
                    filename = paste0("ITS_cancerspecific_gene/",cancer,"/ITSp_vs_regulator_p_PMID34980132"))
    }, silent = T)

  # ------------------------------------------------------------------------------
  # regulator_n
  immuneModulator_PMID34980132_2 <- as.data.frame(read_excel("06_results_summaryandplotting/data/immune_modulator/PMID34980132/positive and negative regulators.xlsx", sheet = 1))
  regulator_n <- unique(immuneModulator_PMID34980132_2[immuneModulator_PMID34980132_2$Type == "Negative regulators", ]$Symbol)

    try({
  ITS_cancerspecific_characteristic(ITS_p_cancer_specific = gset, 
                    genelist = regulator_n, 
                    filename = paste0("ITS_cancerspecific_gene/",cancer,"/ITSp_vs_regulator_n_PMID34980132"))
    }, silent = T)

  # ------------------------------------------------------------------------------
  # ITS gene vs resistor and sensitizor  -----------------------------------------
  # ------------------------------------------------------------------------------

  immuneModulator_PMID34980132_3 <- as.data.frame(read_excel("06_results_summaryandplotting/data/immune_modulator/PMID34980132/sensitizer_resistor_genes.xlsx", sheet = 1))
  # ------------------------------------------------------------------------------
  # Sensitizer
  Sensitizer <- unique(immuneModulator_PMID34980132_3[immuneModulator_PMID34980132_3$Type == "Sensitizer", ]$Symbol)
    try({
  ITS_cancerspecific_characteristic(ITS_p_cancer_specific = gset, 
                    genelist = Sensitizer, 
                    filename = paste0("ITS_cancerspecific_gene/",cancer,"/ITSp_vs_Sensitizer_PMID34980132"))
    }, silent = T)

  # ------------------------------------------------------------------------------
  # Resistor
  Resistor <- unique(immuneModulator_PMID34980132_3[immuneModulator_PMID34980132_3$Type == "Resistor", ]$Symbol)
    try({
  ITS_cancerspecific_characteristic(ITS_p_cancer_specific = gset, 
                    genelist = Resistor, 
                    filename = paste0("ITS_cancerspecific_gene/",cancer,"/ITSp_vs_Resistor_PMID34980132"))
    }, silent = T)

  # ------------------------------------------------------------------------------
  # ITS gene vs OG and TSG -------------------------------------------------------
  # ------------------------------------------------------------------------------

  immuneModulator_PMID34980132_4 <- as.data.frame(read_excel("06_results_summaryandplotting/data/immune_modulator/PMID34980132/Tumor immunity related OGs TSGs.xlsx", sheet = 1))
  OG <- unique(immuneModulator_PMID34980132_4[immuneModulator_PMID34980132_4$Type == "Tumor immunity-related OGs", ]$Symbol)
    try({
  ITS_cancerspecific_characteristic(ITS_p_cancer_specific = gset, 
                    genelist = OG, 
                    filename = paste0("ITS_cancerspecific_gene/",cancer,"/ITSp_vs_OG_PMID34980132"))
    }, silent = T)

  TSG <- unique(immuneModulator_PMID34980132_4[immuneModulator_PMID34980132_4$Type == "Tumor immunity-related TSGs", ]$Symbol)
    try({
  ITS_cancerspecific_characteristic(ITS_p_cancer_specific = gset, 
                    genelist = TSG, 
                    filename = paste0("ITS_cancerspecific_gene/",cancer,"/ITSp_vs_TSG_PMID34980132"))
    }, silent = T)

  # ------------------------------------------------------------------------------
  # cancer oncogenes -------------------------------------------------------------
  # ------------------------------------------------------------------------------
  canceroncogene <- read.csv("06_results_summaryandplotting/data/canceroncegenes/oncoKB_cancerGeneList.tsv", sep = "\t", header = T)
  oncogene <- canceroncogene[canceroncogene$Is.Oncogene == "Yes",]$Hugo.Symbol
    try({
  ITS_cancerspecific_characteristic(ITS_p_cancer_specific = gset, 
                    genelist = oncogene, 
                    filename = paste0("ITS_cancerspecific_gene/",cancer,"/ITSp_vs_OG_oncoKB"))
    }, silent = T)

  suppressorgene <- canceroncogene[canceroncogene$Is.Tumor.Suppressor.Gene == "Yes",]$Hugo.Symbol
    try({
  ITS_cancerspecific_characteristic(ITS_p_cancer_specific = gset, 
                    genelist = suppressorgene, 
                    filename = paste0("ITS_cancerspecific_gene/",cancer,"/ITSp_vs_TSG_oncoKB"))
    }, silent = T)


})


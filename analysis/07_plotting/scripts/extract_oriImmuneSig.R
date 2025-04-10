# Mapping drug DEGs to correlation results
print("start computing")
work_path <- "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/"

setwd(work_path)
options(stringsAsFactors = F)
library(parallel)
library(reshape2)
library(stringr)
library(dplyr)

source("06_results_summaryandplotting/scripts/immunesig_selection_function.R")

load("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/02_tumorGene_immuneSig_correlation/results_selected/immunecell_TCGAgenes_result_spearman_TUMERIC_r0.4/pcor_spearman_SKCM_Metastatic_positive_200.Rdata")
ITSname = data.frame(signame = names(genelist_p), immune_sig = tolower(names(genelist_p)))
ITSname = immunesig_name_unify(ITSname)

immuneSigGSpath = "02_tumorGene_immuneSig_correlation/data/immune_signatures/originalSignaturesGenes/"
    library(readxl)
    library(dplyr)


immunesigGS1 <- as.data.frame(read_excel(paste0(immuneSigGSpath, "immunesig_original_used_correctedPMID.xlsx"), sheet = 1))

    immunesigGS1 <- immunesigGS1[-c(2:5)]
    immunesigGS1Name <- as.list(immunesigGS1$immune_sig)
    immunesigGS1 <- lapply(immunesigGS1Name, function(x){
     y=immunesigGS1[immunesigGS1$immune_sig == x,][-1]
     unlist(y[-which(is.na(y))])
    })
    names(immunesigGS1) = unlist(immunesigGS1Name)
    immunesigGS2 <- read_excel(paste0(immuneSigGSpath, "immunesig_original_used_correctedPMID.xlsx"), sheet = 2)
    cancer='SKCM'
    immunesigGS2 <- as.data.frame(immunesigGS2[immunesigGS2$CancerType == cancer, ])
    immunesigGS2Name <- as.list(immunesigGS2$Signatures)
    immunesigGS2 <- lapply(immunesigGS2Name, function(x){
     y=immunesigGS2[immunesigGS2$Signatures == x,][-1]
     unlist(y[-which(is.na(y))])
    })
    names(immunesigGS2) = paste0(unlist(immunesigGS2Name),"_ConsensusTMEgsva")
    immunesigGS = c(immunesigGS1, immunesigGS2)
    immunesigGSname = data.frame(signame = names(immunesigGS), immune_sig = tolower(names(immunesigGS)))
    immunesigGSname = immunesig_name_unify(immunesigGSname)
    names(immunesigGS) = immunesigGSname$immune_sig

immunesigGS1 <- as.data.frame(read_excel(paste0(immuneSigGSpath, "immunesig_original_used_correctedPMID.xlsx"), sheet = 1))
names(immunesigGS1)[1] = 'signame'
immunesigGS1 = merge(immunesigGSname, immunesigGS1, by='signame')
immunesigGS1 =immunesigGS1[is.element(immunesigGS1$immune_sig, ITSname$immune_sig),]

immuneSigInfoPath = "06_results_summaryandplotting/data/immene_cell_hieracy/"
immunesigInfo <- as.data.frame(readxl::read_excel(paste0(immuneSigInfoPath, "immune_sig_hieracy_new.xlsx"), sheet = 1))
immunesigInfo <- immunesigInfo[c(1,5:9)]
immunesigInfo$immune_sig = tolower(immunesigInfo$immune_sig)
immunesigInfo = immunesig_name_unify(immunesigInfo)

immunesigGS1 = merge(immunesigInfo, immunesigGS1, by='immune_sig')

write.table(unique(immunesigGS1), paste0(immuneSigGSpath, "immunesig_original_usedinpaper.txt"), sep='\t', quote=F, row.names=F)

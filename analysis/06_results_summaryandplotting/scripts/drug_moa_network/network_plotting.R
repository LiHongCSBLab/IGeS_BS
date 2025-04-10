# setwd("G:/lab/Projects/p2_1_immunotherapy_targetedtherapy/")
setwd("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/")

options(stringsAsFactors = F)
library(igraph)
library(devtools)
library(visNetwork)
library(dplyr)

source("r1_drug_immuneSig_correlation_analysis/06_results_summaryandplotting/scripts/network_plotting_function.R")

cancer = "LIHC"
purity_method = "TUMERIC"
dataset = "70138"
cellline = "HEPG2"
filepath = "r1_drug_immuneSig_correlation_analysis/06_results_summaryandplotting/results/drugDEGfc_immunesig/"


datatype = 'allgenes'
dataset = '70138'
drugMergeFCpath = "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/s2_drug_inducted_transcriptome_change/mergeFC_r1/"
immusig_path <- "r1_drug_immuneSig_correlation_analysis/02_tumorGene_immuneSig_correlation/results/"
resultpath = paste0("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/06_results_summaryandplotting/results/drugDEGfc_immunesig/", purity_method)

# generating interested drug-Target and DEG - immunesig files for network plotting
# drug_keyword: can choose one drug name or part of drug name that matches certain type of drug
# e.g.  drug_keyword = 'nib' to generate file for all drugs name includes 'nib'

drug_target_DEG_immunesig_match(cancer = cancer,
                                purity_method = purity_method,
                                datatype = datatype,
                                dataset = dataset,
                                drug_keyword = NULL, # 'nib',
                                r = NULL,
                                padj = 0.05,
                                p = NULL,
                                drug_logFC = NULL,
                                drug_padj = 0.1,
                                drug_p = NULL,
                                drugMergeFCpath = drugMergeFCpath,
                                immusig_path = immusig_path,
                                resultpath = resultpath)


# checking which drugs are in
druglist = listDrug(cancer = cancer,
         purity_method = purity_method,
         dataset = dataset,
         cellline = cellline,
         filepath = filepath)


# checking which immune signatures are included
listImmuneSig(cancer = cancer,
              purity_method = purity_method,
              dataset = dataset,
              cellline = cellline,
              immuneKeyword = "OE_",
              filepath = filepath)

# checking which immune signatures to plot
immuneKeywordList = c("TIDE_TIDE","TIP_signature_TIP","IPS_IPS", "OE_resu","OE_trt", "OE_Nivolumab.resistant.melanoma.down")
# immuneKeyword = "T.cell.CD4._epic"
# network plotting for one specific interested drug
resultpath = "r1_drug_immuneSig_correlation_analysis/06_results_summaryandplotting/results/networkAnalysis/"
druglist = druglist$drug_index

for(immuneKeyword in immuneKeywordList){

    network_visualization(cancer = cancer,
                            purity_method = purity_method,
                            dataset = dataset,
                            cellline = cellline,
                            druglist = druglist,
                            immuneKeyword = immuneKeyword,
                            PPIresource = "STRING",
                            PPIconfidentScore = 0.7,
                            r = NULL,
                            padj = 0.05,
                            p = NULL,
                            drug_logFC = NULL,
                            drug_padj = NULL,
                            drug_p = NULL,
                            filepath = filepath,
                            immusig_path = immusig_path, 
                            resultpath = resultpath)
}


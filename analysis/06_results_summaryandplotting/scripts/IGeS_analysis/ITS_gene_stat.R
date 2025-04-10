# Mapping drug DEGs to correlation results
# work_path <- "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/"
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
 
source("06_results_summaryandplotting/scripts/plot_for_drug_function.R")
source("06_results_summaryandplotting/scripts/ITS_analysis/ITS_gene_stat_function.R")
source("06_results_summaryandplotting/scripts/immunesig_selection_function.R")
source("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/scripts/drugITGSscore_lincs.R")






# ----------------------------------------------------------------------------
# drug-induced immune signatures NES
purity_method = "TUMERIC"
savepath = paste0("06_results_summaryandplotting/results_", purity_method, "/")
dir.create(savepath)

dir.create("06_results_summaryandplotting/data/immune_sig_rough_filter_temp/")

cancer = "SKCM"
tumor_subtype = "Metastatic"
datatype="allgenes"
num_gene = 200
purity_method = "TUMERIC"


# immunesig <- immune_sig_filter_v2(cancer = cancer,
#                                   tumor_subtype = tumor_subtype,
#                                   # purity_method = purity_method,
#                                   datatype = "allgenes",
#                                   # num_gene = 200,
#                                   auc_threshold = 0.6,
#                                   proof = 1, # how many datasets support the immune sigture to have prediction power
#                                   savepath = savepath)

# # confident = 1: only support by survival analysis
# # confident < 1: only support by response analysis
# # confiden = 1/6: only support by one dataset in response analysis
# immunesig2 = immunesig[immunesig$confident > min(immunesig$confident),][c("immune_sig", "setindex")] #   & immunesig$confident!=1
# immunesig2[immunesig2$setindex == '01', ]$setindex = "01_sensitive"
# immunesig2[immunesig2$setindex == '02', ]$setindex = "02_resistant"
# immunesig2 = unique(immunesig2)
 
res = immuneSigGene_ITS_found_lincs(
                    cancer = "SKCM",
                    tumor_subtype = "Metastatic",
                    purity_method = "TUMERIC",
                    datatype="allgenes",
                    num_gene = 200,
                    r = 0.4,
                    ITSproof = 2,
                    cluster_method = "hclust",
                    similarity_method = "Otsuka_Ochiai",
                    genevote = 0.5,
                    workpath = "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/", 
                    immuneSigGSpath = "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/02_tumorGene_immuneSig_correlation/data/immune_signatures/originalSignaturesGenes/",
                    savepath = "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/"
                    )

# "UCS", "UVM", "DLBC", "ACC","CHOL","KIRC",  "MESO","PCPG","SARC","TGCT","THYM", "KICH",
cancerlist <- unique(gsub("spearman_","", 
                 gsub("_Metastatic_positive_200.Rdata","",
                 gsub("_Metastatic_negative_200.Rdata","",
                 gsub("_Primary_positive_200.Rdata","",
                 gsub("_Primary_negative_200.Rdata","",
                 list.files("03_drug_immuneSig_enrichment/data/gene_immune_sig_TUMERIC_spearman/for_GSVA/pcor")))))))
# cancerlist <- c( "UCEC","LIHC", "PRAD", "PAAD", "BRCA", "READ","LUSC", # "LUAD",
#                   "OV",  "COAD", 
#                   "HNSC", "BLCA","CESC", # "ESCA", 
#                   "LGG", "GBM", 
#                   "STAD", "KIRP","THCA") 
 
tumor_subtype = "Primary"

for(cancer in cancerlist[1:length(cancerlist)]){ # [-18][6:20]

    res = immuneSigGene_ITS_found_lincs(
                        cancer = cancer,
                        tumor_subtype = tumor_subtype,
                        purity_method = purity_method,
                        datatype="allgenes",
                        num_gene = 200,
                        r = 0.4,
                        ITSproof = 2,
                        cluster_method = "hclust",
                        similarity_method = "Otsuka_Ochiai" ,
                        genevote = 0.5,
                        workpath = "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/", 
                        immuneSigGSpath = "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/02_tumorGene_immuneSig_correlation/data/immune_signatures/originalSignaturesGenes/",
                        savepath = "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/"
                        )
}


# for correlation without purity correction

# cancerlist <- unique(gsub("spearman_","", 
#                  gsub("_Metastatic_positive_200.Rdata","",
#                  gsub("_Metastatic_negative_200.Rdata","",
#                  gsub("_Primary_positive_200.Rdata","",
#                  gsub("_Primary_negative_200.Rdata","",
#                  list.files("03_drug_immuneSig_enrichment/data/gene_immune_sig_TUMERIC_spearman/for_GSVA/cor")))))))

# tumor_subtype = "Primary"

# for(cancer in cancerlist){
#     res = immuneSigGene_ITS_found_lincs(
#                         cancer = cancer,
#                         tumor_subtype = tumor_subtype,
#                         purity_method = purity_method,
#                         datatype="allgenes",
#                         num_gene = 200,
#                         r = 0.4,
#                         ITSproof = 2,
#                         cluster_method = "hclust",
#                         similarity_method = "Otsuka_Ochiai" ,
#                         genevote = 0.5,
#                         workpath = "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/", 
#                         immuneSigGSpath = "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/02_tumorGene_immuneSig_correlation/data/immune_signatures/originalSignaturesGenes/",
#                         savepath = "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/"
#                         )
# }

# ----------------------------------------------------------------------------
# drug-induced immune signatures NES
purity_method = "CPE"
savepath = paste0("06_results_summaryandplotting/results_", purity_method, "/")
dir.create(savepath)

dir.create("06_results_summaryandplotting/data/immune_sig_rough_filter_temp/")

cancer = "SKCM"
tumor_subtype = "Primary"
datatype="allgenes"
num_gene = 200


res = immuneSigGene_ITS_found_lincs(
                        cancer = cancer,
                        tumor_subtype = tumor_subtype,
                        purity_method = purity_method,
                        datatype="allgenes",
                        num_gene = 200,
                        r = 0.4,
                        ITSproof = 2,
                        cluster_method = "hclust",
                        similarity_method = "Otsuka_Ochiai" ,
                        genevote = 0.5,
                        workpath = "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/", 
                        immuneSigGSpath = "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/02_tumorGene_immuneSig_correlation/data/immune_signatures/originalSignaturesGenes/",
                        savepath = "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/"
                        )

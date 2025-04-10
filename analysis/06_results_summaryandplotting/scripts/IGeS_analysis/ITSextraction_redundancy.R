
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
# immunesig <- immune_sig_filter_v1(cancer = cancer,
#                               tumor_subtype = tumor_subtype,
#                               purity_method = purity_method,
#                               datatype = "allgenes",
#                               num_gene = 200,
#                               auc_threshold = 0.6,
#                               proof = 1, # how many datasets support the immune sigture to have prediction power
#                               savepath = savepath)

immunesig <- immune_sig_filter_v2(cancer = cancer,
                              tumor_subtype = tumor_subtype,
                              purity_method = purity_method,
                              datatype = "allgenes",
                              num_gene = 200,
                              auc_threshold = 0.6,
                              proof = 1, # how many datasets support the immune sigture to have prediction power
                              savepath = savepath)

# confident = 1: only support by survival analysis
# confident < 1: only support by response analysis
# confiden = 1/6: only support by one dataset in response analysis
immunesig2 = immunesig[immunesig$confident > min(immunesig$confident),][c("immune_sig", "setindex")] #   & immunesig$confident!=1
immunesig2[immunesig2$setindex == '01', ]$setindex = "01_sensitive"
immunesig2[immunesig2$setindex == '02', ]$setindex = "02_resistant"

immunesig_jci <- immune_sig_finefilter(
                        immunesig = immunesig2,
                        cancer = "SKCM",
                        tumor_subtype = "Metastatic",
                        purity_method = "TUMERIC",
                        datatype="allgenes",
                        num_gene = 200,
                        workpath = "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/", 
                        savepath = "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/" )

# immuneSig_hieracy <- read.csv("06_results_summaryandplotting/data/immene_cell_hieracy/immune_sig_rough_filter.csv", row.names=1)
# immuneSig_hieracy[which(immuneSig_hieracy==''),] <- NA
# immuneSig_hieracy$immune_sig = row.names(immuneSig_hieracy)
# immuneSig_hieracy_selected <- merge(immunesig, immuneSig_hieracy, all.x=T, by = 'immune_sig')

keywords = unique(immunesig[,'Type'])
keywords = keywords[keywords != "" ]
keywords = keywords[!is.na(keywords)]

for(keyword in keywords){

  immunesig_jci <- immune_sig_finefilter(
                        immunesig = immunesig2,
                        keyword = keyword, # substr(keyword,1, 6),
                        row_selected = 'Type',
                        cancer = "SKCM",
                        tumor_subtype = "Metastatic",
                        purity_method = "TUMERIC",
                        datatype="allgenes",
                        num_gene = 200,
                        workpath = "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/", 
                        savepath = "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/" )
}



for(cancer in c("LIHC", "BLCA", "BRCA",  "KIRC", "KIRP","COAD",  "LUAD", "LUSC","OV", 
                 "PAAD", "PRAD", "READ",  "STAD", "UCEC", "THCA", 
                "CESC", "GBM", "HNSC", "LGG","ESCA")){

    # cancer = "LIHC"
    tumor_subtype = "Primary"
    purity_method = "TUMERIC"
    datatype = "allgenes"
    num_gene = 200

    savepath = paste0("06_results_summaryandplotting/results_", purity_method,"_V2/")
    dir.create(savepath)

    immunesig <- immune_sig_filter_v1(cancer = cancer,
                                  tumor_subtype = tumor_subtype,
                                  purity_method = purity_method,
                                  datatype = "allgenes",
                                  num_gene = 200,
                                  auc_threshold = 0.6,
                                  proof = 1, # how many datasets support the immune sigture to have prediction power
                                  savepath = savepath)

    immunesig <- immune_sig_filter_v2(cancer = cancer,
                                  tumor_subtype = tumor_subtype,
                                  purity_method = purity_method,
                                  datatype = "allgenes",
                                  num_gene = 200,
                                  auc_threshold = 0.6,
                                  proof = 1, # how many datasets support the immune sigture to have prediction power
                                  savepath = savepath)

    # confident = 1: only support by survival analysis
    # confident < 1: only support by response analysis
    # confiden = 1/6: only support by one dataset in response analysis
    immunesig2 = immunesig[immunesig$confident > min(immunesig$confident),][c("immune_sig", "setindex")] #   & immunesig$confident!=1
    immunesig2[immunesig2$setindex == '01', ]$setindex = "01_sensitive"
    immunesig2[immunesig2$setindex == '02', ]$setindex = "02_resistant"

    immunesig_jci <- immune_sig_finefilter(
                            immunesig = immunesig2,
                            cancer = cancer,
                            tumor_subtype = tumor_subtype,
                            purity_method = purity_method,
                            datatype="allgenes",
                            num_gene = 200,
                            workpath = "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/", 
                            savepath = "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/" )

    # immuneSig_hieracy <- read.csv("06_results_summaryandplotting/data/immene_cell_hieracy/immune_sig_rough_filter.csv", row.names=1)
    # immuneSig_hieracy[which(immuneSig_hieracy==''),] <- NA
    # immuneSig_hieracy$immune_sig = row.names(immuneSig_hieracy)
    # immuneSig_hieracy_selected <- merge(immunesig, immuneSig_hieracy, all.x=T, by = 'immune_sig')

    keywords = unique(immunesig[,'Type'])
    keywords = keywords[keywords != "" ]
    keywords = keywords[!is.na(keywords)]

    for(keyword in keywords){

      immunesig_jci <- immune_sig_finefilter(
                            immunesig = immunesig2,
                            keyword = keyword, # substr(keyword,1, 6),
                            row_selected = 'Type',
                            cancer = cancer,
                            tumor_subtype = tumor_subtype,
                            purity_method = purity_method,
                            datatype="allgenes",
                            num_gene = 200,
                            workpath = "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/", 
                            savepath = "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/" )
    }
    
}


# p="06_results_summaryandplotting/data/immune_sig_selected/TUMERIC/"
# f <- list.files("06_results_summaryandplotting/data/immune_sig_selected/TUMERIC/")
# immunetemp <- read.csv(paste0(p, '/', f[1]))
# for(i in 2:length(f)){
#     immunetemp1 <- read.csv(paste0(p, '/', f[i]))
#     immunetemp = unique(rbind(immunetemp, immunetemp1))
# }
# immunetemp = unique(immunetemp[2])
# write.csv(unique(immunetemp), "06_results_summaryandplotting/data/immune_sig_rough_filter_temp.csv", quote = F, row.names = F)
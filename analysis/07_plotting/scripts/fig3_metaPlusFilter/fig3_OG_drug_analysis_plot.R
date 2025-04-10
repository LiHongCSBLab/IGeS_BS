# rm(list=ls())
work_path <- "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/"
setwd(work_path)
options(stringsAsFactors = F)

# library(reshape2)
# library(ggplotify)
# library(grid)
# library(cowplot)
# library(parallel)
# library(patchwork)
# library(dplyr)
# library(stringr)
# library(data.table)
# require(GGally)
# library(plot3D)
# library(fmsb)
# library(readxl)
library(ggplot2)
library(dplyr)
library(UpSetR)
library(clusterProfiler)
library(enrichplot)
library(ggrepel)
library(patchwork)
library(aplot)
library(reshape2)

# source("06_results_summaryandplotting/scripts/plot_for_drug_function.R")
# source("06_results_summaryandplotting/scripts/ITS_analysis/ITS_gene_stat_function.R")
# source("06_results_summaryandplotting/scripts/immunesig_selection_function.R")
# source("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/scripts/drugITGSscore_lincs.R")
source("07_plotting_v2/scripts/functions/ITS_enrich.R")
source("07_plotting_v2/scripts/functions/ITS_function_detection.R")
source("06_results_summaryandplotting/scripts/immunesig_selection_function.R")
source("07_plotting_v2/scripts/functions/ITS_enrich.R")
source("07_plotting_v2/scripts/functions/ITS_function_detection.R")
source("07_plotting_v2/scripts/functions/ITS_g_cancer_plot_function.R")



# ------------------------------------------------------------------------------

load( "07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITS_set.Rdata")
meta_result <- read.csv("05_immuneSig_ICBresponse/results/meta_results/res_metaRes.csv")
meta_selected <- meta_result[meta_result$Meta_Pval < 0.05, ]
meta_selected$immune_sig = meta_selected$X

load("05_immuneSig_ICBresponse/results_r0.4_ITSssgsea_CombatRef/finalmodel_17datasets_ssgseanorm/metaAnalysis_filter/en/elasticnet_model_pred.Rdata")

en_model_para = m_en$learner.model$glmnet.fit
m_en_imp <- data.frame(en_model_para$beta[,which(en_model_para$lambda == m_en$learner.model$lambda.min)])
names(m_en_imp) <- "weight"
m_en_imp$immune_sig = row.names(m_en_imp)
meta_selected = merge(meta_selected, m_en_imp, by = "immune_sig")
meta_selected = meta_selected[meta_selected$weight != 0, ]


meta_selected$flag = 'sensitive'
meta_selected[meta_selected$OR < 1, ]$flag = 'resistant'
ITS_S = meta_selected[meta_selected$flag == 'sensitive', ]$X
ITS_R = meta_selected[meta_selected$flag == 'resistant', ]$X

ITS_s_selected <- lapply(ITS_p_set, function(x){x[ITS_S]})
ITS_s_selected <- lapply(as.list(names(ITS_s_selected)), function(x){
    tmp = ITS_s_selected[[x]]
    do.call(rbind, lapply(as.list(names(tmp)), function(y){
        if(length( tmp[[y]]) == 0){
            data.frame(cancer = x, immune_sig = y, genename = NA)
        }else{
            data.frame(cancer = x, immune_sig = y, genename = tmp[[y]])
        }
    }))
})
ITS_s_selected <- do.call(rbind,ITS_s_selected)

ITS_r_selected <- lapply(ITS_p_set, function(x){x[ITS_R]})
ITS_r_selected <- lapply(as.list(names(ITS_r_selected)), function(x){
    tmp = ITS_r_selected[[x]]
    do.call(rbind, lapply(as.list(names(tmp)), function(y){
        if(length( tmp[[y]]) == 0){
            data.frame(cancer = x, immune_sig = y, genename = NA)
        }else{
            data.frame(cancer = x, immune_sig = y, genename = tmp[[y]])
        }
    }))
})
ITS_r_selected <- do.call(rbind,ITS_r_selected)


# 2) sITS-OG-Drug -------------------------------------------------------------------
filepath <- ("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/drug_merge3DB/GSEAenrich/ITSp_s_vs_OG_oncoKB")
OG_drug <- read.table(paste0(filepath, ".txt"), sep = '\t',header = T)
OG_drugFDA <- read.table(paste0(filepath, "_FDA.txt"), sep = '\t', header = T)
OG_drugFDAcancer <- OG_drugFDA[which(OG_drugFDA$ifCancer_drugbank == "yes"), ]

OG_drug_ITS_s_selected <- merge(OG_drug, ITS_s_selected, by = 'genename')
res_OG_drug_ITS_s_selected = do.call(rbind, lapply(as.list(unique(OG_drug_ITS_s_selected$drugname)), function(drug){
    dt = OG_drug_ITS_s_selected[OG_drug_ITS_s_selected$drugname == drug, ]
    res = data.frame(colSums(table(unique(dt[c('cancer','immune_sig')]))))
    names(res) = 'numCancer'
    res$immune_sig = rownames(res)
    res$drugname = drug
    return(res)
}))
res_OG_drug_ITS_s_selected = unique(merge(res_OG_drug_ITS_s_selected, 
                                      OG_drug_ITS_s_selected[c('drugname','genename','MoA','immune_sig')], 
                                      by = c('drugname', 'immune_sig')))
write.csv(res_OG_drug_ITS_s_selected, '07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/drug_merge3DB/GSEAenrich/res_OG_drug_ITS_s_selected.csv', row.names = F, quote = F)

OG_drugFDA = unique(OG_drugFDA[-6])
OG_drugFDA_ITS_s_selected <- merge(OG_drugFDA, ITS_s_selected, by = 'genename')
res_OG_drugFDA_ITS_s_selected = do.call(rbind, lapply(as.list(unique(OG_drugFDA_ITS_s_selected$drugname)), function(drug){
    dt = OG_drugFDA_ITS_s_selected[OG_drugFDA_ITS_s_selected$drugname == drug, ]
    res = data.frame(colSums(table(unique(dt[c('cancer','immune_sig')]))))
    names(res) = 'numCancer'
    res$immune_sig = rownames(res)
    res$drugname = drug
    return(res)
}))
res_OG_drugFDA_ITS_s_selected = unique(merge(res_OG_drugFDA_ITS_s_selected, 
                                      OG_drugFDA_ITS_s_selected[c('drugname','genename','MoA','immune_sig')], 
                                      by = c('drugname', 'immune_sig')))
write.csv(res_OG_drugFDA_ITS_s_selected, '07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/drug_merge3DB/GSEAenrich/res_OG_drugFDA_ITS_s_selected.csv', row.names = F, quote = F)

OG_drugFDAcancer = unique(OG_drugFDAcancer[-6])
OG_drugFDAcancer_ITS_s_selected <- merge(OG_drugFDAcancer, ITS_s_selected, by = 'genename')
res_OG_drugFDAcancer_ITS_s_selected = do.call(rbind, lapply(as.list(unique(OG_drugFDAcancer_ITS_s_selected$drugname)), function(drug){
    dt = OG_drugFDAcancer_ITS_s_selected[OG_drugFDAcancer_ITS_s_selected$drugname == drug, ]
    res = data.frame(colSums(table(unique(dt[c('cancer','immune_sig')]))))
    names(res) = 'numCancer'
    res$immune_sig = rownames(res)
    res$drugname = drug
    return(res)
}))
res_OG_drugFDAcancer_ITS_s_selected = unique(merge(res_OG_drugFDAcancer_ITS_s_selected, 
                                      OG_drugFDAcancer_ITS_s_selected[c('drugname','genename','MoA','immune_sig')], 
                                      by = c('drugname', 'immune_sig')))

write.csv(res_OG_drugFDAcancer_ITS_s_selected, '07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/drug_merge3DB/GSEAenrich/res_OG_drugFDAcancer_ITS_s_selected.csv', row.names = F, quote = F)


# 2) rITS-OG-Drug -------------------------------------------------------------------
filepath <- ("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/drug_merge3DB/GSEAenrich/ITSp_r_vs_OG_oncoKB")
OG_drug <- read.table(paste0(filepath, ".txt"), sep = '\t',header = T)
OG_drugFDA <- read.table(paste0(filepath, "_FDA.txt"), sep = '\t', header = T)
OG_drugFDAcancer <- OG_drugFDA[which(OG_drugFDA$ifCancer_drugbank == "yes"), ]

OG_drug_ITS_r_selected <- merge(OG_drug, ITS_r_selected, by = 'genename')
res_OG_drug_ITS_r_selected = do.call(rbind, lapply(as.list(unique(OG_drug_ITS_r_selected$drugname)), function(drug){
    dt = OG_drug_ITS_r_selected[OG_drug_ITS_r_selected$drugname == drug, ]
    res = data.frame(colSums(table(unique(dt[c('cancer','immune_sig')]))))
    names(res) = 'numCancer'
    res$immune_sig = rownames(res)
    res$drugname = drug
    return(res)
}))
res_OG_drug_ITS_r_selected = unique(merge(res_OG_drug_ITS_r_selected, 
                                      OG_drug_ITS_r_selected[c('drugname','genename','MoA','immune_sig')], 
                                      by = c('drugname', 'immune_sig')))
write.csv(res_OG_drug_ITS_r_selected, '07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/drug_merge3DB/GSEAenrich/res_OG_drug_ITS_r_selected.csv', row.names = F, quote = F)


OG_drugFDA = unique(OG_drugFDA[-6])
OG_drugFDA_ITS_r_selected <- merge(OG_drugFDA, ITS_r_selected, by = 'genename')
OG_drugFDAcancer_ITS_r_selected <- merge(OG_drugFDAcancer, ITS_r_selected, by = 'genename')

OG_drugFDA = unique(OG_drugFDA[-6])
OG_drugFDA_ITS_r_selected <- merge(OG_drugFDA, ITS_r_selected, by = 'genename')
res_OG_drugFDA_ITS_r_selected = do.call(rbind, lapply(as.list(unique(OG_drugFDA_ITS_r_selected$drugname)), function(drug){
    dt = OG_drugFDA_ITS_r_selected[OG_drugFDA_ITS_r_selected$drugname == drug, ]
    res = data.frame(colSums(table(unique(dt[c('cancer','immune_sig')]))))
    names(res) = 'numCancer'
    res$immune_sig = rownames(res)
    res$drugname = drug
    return(res)
}))
res_OG_drugFDA_ITS_r_selected = unique(merge(res_OG_drugFDA_ITS_r_selected, 
                                      OG_drugFDA_ITS_r_selected[c('drugname','genename','MoA','immune_sig')], 
                                      by = c('drugname', 'immune_sig')))

write.csv(res_OG_drugFDA_ITS_r_selected, '07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/drug_merge3DB/GSEAenrich/res_OG_drugFDA_ITS_r_selected.csv', row.names = F, quote = F)

OG_drugFDAcancer = unique(OG_drugFDAcancer[-6])
OG_drugFDAcancer_ITS_r_selected <- merge(OG_drugFDAcancer, ITS_r_selected, by = 'genename')
res_OG_drugFDAcancer_ITS_r_selected = do.call(rbind, lapply(as.list(unique(OG_drugFDAcancer_ITS_r_selected$drugname)), function(drug){
    dt = OG_drugFDAcancer_ITS_r_selected[OG_drugFDAcancer_ITS_r_selected$drugname == drug, ]
    res = data.frame(colSums(table(unique(dt[c('cancer','immune_sig')]))))
    names(res) = 'numCancer'
    res$immune_sig = rownames(res)
    res$drugname = drug
    return(res)
}))
res_OG_drugFDAcancer_ITS_r_selected = unique(merge(res_OG_drugFDAcancer_ITS_r_selected, 
                                      OG_drugFDAcancer_ITS_r_selected[c('drugname','genename','MoA','immune_sig')], 
                                      by = c('drugname', 'immune_sig')))

write.csv(res_OG_drugFDAcancer_ITS_r_selected, '07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/drug_merge3DB/GSEAenrich/res_OG_drugFDAcancer_ITS_r_selected.csv', row.names = F, quote = F)

res_OG_drug_ITS_s_selected$flag = 's'
res_OG_drug_ITS_r_selected$flag = 'r'
res_OG_drug_ITS_selected = rbind(res_OG_drug_ITS_s_selected, res_OG_drug_ITS_r_selected)
write.table(res_OG_drug_ITS_selected, '07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/drug_merge3DB/GSEAenrich/res_OG_drug_ITS_selected.txt', sep = '\t', row.names = F, quote = F)

res_OG_drugFDA_ITS_s_selected$flag = 's'
res_OG_drugFDA_ITS_r_selected$flag = 'r'
res_OG_drugFDA_ITS_selected = rbind(res_OG_drugFDA_ITS_s_selected, res_OG_drugFDA_ITS_r_selected)
write.csv(res_OG_drugFDA_ITS_selected, '07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/drug_merge3DB/GSEAenrich/res_OG_drugFDA_ITS_selected.csv', row.names = F, quote = F)

res_OG_drugFDAcancer_ITS_s_selected$flag = 's'
res_OG_drugFDAcancer_ITS_r_selected$flag = 'r'
res_OG_drugFDAcancer_ITS_selected = rbind(res_OG_drugFDAcancer_ITS_s_selected, res_OG_drugFDAcancer_ITS_r_selected)
write.csv(res_OG_drugFDAcancer_ITS_selected, '07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/drug_merge3DB/GSEAenrich/res_OG_drugFDAcancer_ITS_selected.csv', row.names = F, quote = F)

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
library(UpSetR)
library(clusterProfiler)
library(enrichplot)
# source("06_results_summaryandplotting/scripts/plot_for_drug_function.R")
# source("06_results_summaryandplotting/scripts/ITS_analysis/ITS_gene_stat_function.R")
# source("06_results_summaryandplotting/scripts/immunesig_selection_function.R")
# source("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/scripts/drugITGSscore_lincs.R")
source("07_plotting_v2/scripts/functions/ITS_enrich.R")
source("07_plotting_v2/scripts/functions/ITS_function_detection.R")



# ----------------------------------------------------------------------------

dir.create("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/")
dir.create("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSp_sharedgene/drug_merge3DB")
dir.create("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSp_sharedgene/drug_merge3DB/GSEAenrich/")
dir.create("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSn_sharedgene/drug_merge3DB")
dir.create("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSn_sharedgene/drug_merge3DB/GSEAenrich/")

# LOAD GENES IN LINCS ------------------------------------------------------
lincs_gene = read.csv("/picb/bigdata/project/FengFYM/drug_MoA_reanalysis/LINCS_data/GSE70138/file/GSE92742_Broad_LINCS_gene_info.txt", sep = "\t")

# # LOAD drug-target from drugbank ------------------------------------------------------
# drugTarget_drugbank <- as.data.frame(fread("/picb/bigdata/project/FengFYM/s1_drug_databases/drugbank/drug_proteins.tsv", sep='\t',quote=""))
# drugTarget_drugbank <- drugTarget_drugbank[drugTarget_drugbank$category == "target",]
# drugindication <- as.data.frame(fread("/picb/bigdata/project/FengFYM/s1_drug_databases/drugbank/drugbank_indication.tsv", sep='\t',quote=""))
# DT_drugbank <- merge(drugindication, drugTarget_drugbank, by = 'drugbank_id') 


# drugFDA_drugbank <- as.data.frame(fread("/picb/bigdata/project/FengFYM/s1_drug_databases/drugbank/drugbank-approved.tsv", sep='\t',quote=""))
# DTFDA_drugbank <- merge(drugFDA_drugbank, drugTarget_drugbank, by = 'drugbank_id') 
# drugFDAcancer_drugbank <- as.data.frame(fread("/picb/bigdata/project/FengFYM/s1_drug_databases/drugbank/drugbank-approved-neoplastic.tsv", sep='\t',quote=""))
# DTFDAcancer_drugbank <- merge(drugFDAcancer_drugbank, drugTarget_drugbank, by = 'drugbank_id') 

# # LOAD drug-target from TTD ------------------------------------------------------
# dt_ttd <- read_xlsx('/picb/bigdata/project/FengFYM/s1_drug_databases/ttd/P1-07-Drug-TargetMapping.xlsx')
# dt_ttd <- as.data.frame(dt_ttd)

# ttd_synomy <- read.csv( '/picb/bigdata/project/FengFYM/s1_drug_databases/ttd/P1-04-Drug_synonyms.txt', sep = '\t', header = T)
# names(ttd_synomy) <- c("DrugID", "DRUGNAME", "Drug_Name")

# ttd_uniprot <- read.csv("/picb/bigdata/project/FengFYM/s1_drug_databases/ttd/P2-01-TTD_uniprot_all.txt",sep="\t", header = T)
# names(ttd_uniprot) <- c('TargetID','UNIPROID','genename_specis')

# dt_ttd_name <- dplyr::inner_join(dt_ttd, ttd_uniprot)
# dt_ttd_name<- dt_ttd_name[dt_ttd_name$UNIPROID=='UNIPROID',]
# dt_ttd_name$genename_specis <- gsub("; ","\t", dt_ttd_name$genename_specis)
# dt_ttd_name$genename_specis <- gsub("-","\t", dt_ttd_name$genename_specis)

# dt=data.table(dt_ttd_name["genename_specis"])
# dt = tstrsplit(dt_ttd_name$genename_specis, "_", fixed=TRUE)[c(1,2)]
# dt = as.data.frame(do.call(cbind, dt))
# names(dt) = c("genename", "specis")
# dt$genename_specis = dt_ttd_name$genename_specis
# dt = unique(dt)

# dt_ttd <- inner_join(dt_ttd_name, dt)
# dt_ttd_approved <- dt_ttd[dt_ttd$Highest_status == "Approved", ]

# # LOAD drug-target from lincsPortal2 ------------------------------------------------------
# drugTarget_lincs <- as.data.frame(fread("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/s3_lincs_drug_MoA/drugTargetMerge.txt"))

# # merge all drugs ------------------------------------------------------
# DT_drugbank1 <- unique(DT_drugbank[c('name', 'genename','actions')])
# names(DT_drugbank1) <- c("drugname", "genename", "MoA")
# DT_drugbank1$database <- "drugbank"
# DT_drugbank1$drugname_lower <- tolower(DT_drugbank1$drugname)

# dt_ttd1 <- unique(dt_ttd[c('Drug_Name', 'genename','MOA')])
# names(dt_ttd1) <- c("drugname", "genename", "MoA")
# dt_ttd1$database <- "TTD"
# dt_ttd1$drugname_lower <- tolower(dt_ttd1$drugname)

# drugTarget_lincs1 <- drugTarget_lincs[c(1,3)]
# names(drugTarget_lincs1) <- c("drugname", "genename")
# drugTarget_lincs1$MoA <- NA
# drugTarget_lincs1$database <- "lincsPortal2"
# drugTarget_lincs1$drugname_lower <- tolower(drugTarget_lincs1$drugname)

# drugTarget_all <- unique(rbind(DT_drugbank1, dt_ttd1, drugTarget_lincs1))
# drugTarget_all <- drugTarget_all[order(drugTarget_all$drugname),]


# # merge Approved drugs ------------------------------------------------------
# dt_ttd_approved1 <- unique(dt_ttd_approved[c('Drug_Name', 'genename','MOA')])
# names(dt_ttd_approved1) <- c("drugname", "genename", "MoA")
# dt_ttd_approved1$database <- "TTD"
# dt_ttd_approved1$drugname_lower <- tolower(dt_ttd_approved1$drugname)


# DTFDA_drugbank1 <- unique(DTFDA_drugbank[c('name', 'genename','actions')])
# names(DTFDA_drugbank1) <- c("drugname", "genename", "MoA")
# DTFDA_drugbank1$database <- "drugbank"
# DTFDA_drugbank1$drugname_lower <- tolower(DTFDA_drugbank1$drugname)


# D_FDAcancer <- unique(DTFDAcancer_drugbank[c("name", "genename", "actions")])
# names(D_FDAcancer) <- c("drugname", "genename", "MoA")
# D_FDAcancer$ifCancer_drugbank <- 'yes'

# drugTarget_FDA <- unique(rbind(dt_ttd_approved1, DTFDA_drugbank1))
# drugTarget_FDA <- merge(drugTarget_FDA, D_FDAcancer, all.x=T, by = c("drugname", "genename", "MoA"))
# drugTarget_FDA <- drugTarget_FDA[order(drugTarget_FDA$drugname),]

# save(drugTarget_drugbank, DTFDA_drugbank, DTFDAcancer_drugbank, dt_ttd, 
#      drugTarget_lincs, drugTarget_all, drugTarget_FDA,
#      file = "07_plotting_v2/data/drugTarget_merged3DB.Rdata")


# ------------------------------------------------------------------------------
load( "07_plotting_v2/data/drugTarget_merged3DB.Rdata")
# drugTarget_FDA
# drugTarget_all

# ------------------------------------------------------------------------------
load("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITS_set.Rdata")

# ------------------------------------------------------------------------------
# stat: sensitive-related and resistant-related ITS ----------------------------
# ------------------------------------------------------------------------------
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

# ITS_S = ITS_S[ITS_S != "merck18_tide"]

# ------------------------------------------------------------------------------
# sensitive 

# ITS_p_s ------------------------------------------------------------------------------
ITS_p_s_notspecific_df <- read.table("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSp_sharedgene/ITS_p_s_notspecific_df.txt", header = T, sep = '\t')
# head(ITS_p_s_notspecific_df)



# enrich custom-GSEA -----------------------------------------------------------

load("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSp_sharedgene/enrich/custom_GSEA_s_scaled.Rdata")
res2 = data.frame(res)

# 1) drugTarget_drugbank-Drug -------------------------------------------------------------------
filepath <- ("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSp_sharedgene/drug_merge3DB/GSEAenrich/ITSp_s_vs_drugtarget")

glist <- unlist(strsplit(res2[res2$ID == "ITSp_vs_drugtarget", ]$core_enrichment, "/"))

# DT_druglist <- unique(drugTarget_all[is.element(drugTarget_all$genename, glist),]$drugname_lower)
# DT_drug <- drugTarget_all[is.element(drugTarget_all$drugname_lower, DT_druglist),]
DT_drug <- unique(drugTarget_all[is.element(drugTarget_all$genename, glist),])

write.table(DT_drug, paste0(filepath, ".txt"), sep = '\t', quote = F, row.names = F)

# DT_druglistFDA <- unique(drugTarget_FDA[is.element(drugTarget_FDA$genename, glist),]$drugname_lower)
# DT_drugFDA<- drugTarget_FDA[is.element(drugTarget_FDA$drugname_lower, DT_druglist),]
DT_drugFDA <- unique(drugTarget_FDA[is.element(drugTarget_FDA$genename, glist),])
write.table(DT_drugFDA, paste0(filepath, "_FDA.txt"), sep = '\t', quote = F, row.names = F)

DT_drug2 = unique(DT_drug[c('drugname_lower', 'genename')])
DT_drugFDA2 = unique(DT_drugFDA[c('drugname_lower', 'genename', 'ifCancer_drugbank')])

drugstat <- data.frame(group = c("num_drug",  "num_drugFDA", "num_drugFDAcancer"),
                    value = c(length(unique(DT_drug2$drugname_lower)),
                              length(unique(DT_drugFDA2$drugname_lower)), 
                              length(unique(DT_drugFDA2[which(DT_drugFDA2$ifCancer_drugbank == 'yes'), ]$drugname_lower))))

write.table(drugstat, paste0(filepath, "_stat.txt"), sep = '\t', quote = F, row.names = F)

# unique(DT_drugFDAcancer$name)
# sort(table(DT_drugFDAcancer$genename))


# 2) OG-Drug -------------------------------------------------------------------
filepath <- ("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSp_sharedgene/drug_merge3DB/GSEAenrich/ITSp_s_vs_OG_oncoKB")

glist <- unlist(strsplit(res2[res2$ID == "ITSp_vs_OG_oncoKB",]$core_enrichment, "/"))

OG_druglist <- unique(drugTarget_all[is.element(drugTarget_all$genename, glist),]$drugname_lower)
OG_drug <- drugTarget_all[is.element(drugTarget_all$drugname_lower, OG_druglist),]

write.table(OG_drug, paste0(filepath, ".txt"), sep = '\t', quote = F, row.names = F)

OG_druglistFDA <- unique(drugTarget_FDA[is.element(drugTarget_FDA$genename, glist),]$drugname_lower)
OG_drugFDA<- drugTarget_FDA[is.element(drugTarget_FDA$drugname_lower, OG_druglist),]
OG_drugFDA<- unique(drugTarget_FDA[is.element(drugTarget_FDA$genename, glist),])


write.table(OG_drugFDA, paste0(filepath, "_FDA.txt"), sep = '\t', quote = F, row.names = F)

OG_drug2 = unique(OG_drug[c('drugname_lower', 'genename')])
OG_drugFDA2 = unique(OG_drugFDA[c('drugname_lower', 'genename', 'ifCancer_drugbank')])

drugstat <- data.frame(group = c("num_drug",  "num_drugFDA", "num_drugFDAcancer"),
                    value = c(length(unique(OG_drug2$drugname_lower)),
                              length(unique(OG_drugFDA2$drugname_lower)), 
                              length(unique(OG_drugFDA2[which(OG_drugFDA2$ifCancer_drugbank == 'yes'), ]$drugname_lower))))

write.table(drugstat, paste0(filepath, "_stat.txt"), sep = '\t', quote = F, row.names = F)




# 3) TSG-Drug -------------------------------------------------------------------
filepath <- ("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSp_sharedgene/drug_merge3DB/GSEAenrich/ITSp_s_vs_TSG_oncoKB")

glist <- unlist(strsplit(res2[res2$ID == "ITSp_vs_TSG_oncoKB",]$core_enrichment, "/"))


# TSG_druglist <- unique(drugTarget_all[is.element(drugTarget_all$genename, glist),]$drugname_lower)
# TSG_drug <- drugTarget_all[is.element(drugTarget_all$drugname_lower, TSG_druglist),]
TSG_drug <- unique(drugTarget_all[is.element(drugTarget_all$genename, glist),])


write.table(TSG_drug, paste0(filepath, ".txt"), sep = '\t', quote = F, row.names = F)

# TSG_druglistFDA <- unique(drugTarget_FDA[is.element(drugTarget_FDA$genename, glist),]$drugname_lower)
# TSG_drugFDA<- drugTarget_FDA[is.element(drugTarget_FDA$drugname_lower, TSG_druglist),]
TSG_drugFDA <- unique(drugTarget_FDA[is.element(drugTarget_FDA$genename, glist),])

write.table(TSG_drugFDA, paste0(filepath, "_FDA.txt"), sep = '\t', quote = F, row.names = F)

TSG_drug2 = unique(TSG_drug[c('drugname_lower', 'genename')])
TSG_drugFDA2 = unique(TSG_drugFDA[c('drugname_lower', 'genename', 'ifCancer_drugbank')])

drugstat <- data.frame(group = c("num_drug",  "num_drugFDA", "num_drugFDAcancer"),
                    value = c(length(unique(TSG_drug2$drugname_lower)),
                              length(unique(TSG_drugFDA2$drugname_lower)), 
                              length(unique(TSG_drugFDA2[which(TSG_drugFDA2$ifCancer_drugbank == 'yes'), ]$drugname_lower))))

write.table(drugstat, paste0(filepath, "_stat.txt"), sep = '\t', quote = F, row.names = F)




# enrich custom-ORA ---------------------------------------------------------------
dir.create("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSp_sharedgene/")
dir.create("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSp_sharedgene/drug_merge3DB/")
dir.create("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSp_sharedgene/drug_merge3DB/ORAenrich/")
filepath = "07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSp_sharedgene/enrich/"
load(paste0(filepath, "custom_ORA_s.Rdata"))

res2 = data.frame(res)

# 1) drugTarget_drugbank-Drug -------------------------------------------------------------------
filepath <- ("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSp_sharedgene/drug_merge3DB/ORAenrich/ITSp_s_vs_drugtarget")

glist <- unlist(strsplit(res2[res2$ID == "ITSp_vs_drugtarget",]$geneID, "/"))
glist <- ITS_p_s_notspecific_df[is.element(ITS_p_s_notspecific_df$ENTREZID, glist), ]$gene_name

# DT_druglist <- unique(drugTarget_all[is.element(drugTarget_all$genename, glist),]$drugname_lower)
# DT_drug <- drugTarget_all[is.element(drugTarget_all$drugname_lower, DT_druglist),]
DT_drug <- unique(drugTarget_all[is.element(drugTarget_all$genename, glist),])

write.table(DT_drug, paste0(filepath, ".txt"), sep = '\t', quote = F, row.names = F)

# DT_druglistFDA <- unique(drugTarget_FDA[is.element(drugTarget_FDA$genename, glist),]$drugname_lower)
# DT_drugFDA<- drugTarget_FDA[is.element(drugTarget_FDA$drugname_lower, DT_druglist),]
DT_drugFDA <- unique(drugTarget_FDA[is.element(drugTarget_FDA$genename, glist),])
write.table(DT_drugFDA, paste0(filepath, "_FDA.txt"), sep = '\t', quote = F, row.names = F)

DT_drug2 = unique(DT_drug[c('drugname_lower', 'genename')])
DT_drugFDA2 = unique(DT_drugFDA[c('drugname_lower', 'genename', 'ifCancer_drugbank')])

drugstat <- data.frame(group = c("num_drug",  "num_drugFDA", "num_drugFDAcancer"),
                    value = c(length(unique(DT_drug2$drugname_lower)),
                              length(unique(DT_drugFDA2$drugname_lower)), 
                              length(unique(DT_drugFDA2[which(DT_drugFDA2$ifCancer_drugbank == 'yes'), ]$drugname_lower))))

write.table(drugstat, paste0(filepath, "_stat.txt"), sep = '\t', quote = F, row.names = F)


# 2) OG-Drug -------------------------------------------------------------------
filepath <- ("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSp_sharedgene/drug_merge3DB/ORAenrich/ITSp_s_vs_OG_oncoKB")

glist <- unlist(strsplit(res2[res2$ID == "ITSp_vs_OG_oncoKB",]$geneID, "/"))
glist <- ITS_p_s_notspecific_df[is.element(ITS_p_s_notspecific_df$ENTREZID, glist), ]$gene_name

# OG_druglist <- unique(drugTarget_all[is.element(drugTarget_all$genename, glist$$gene_name),]$drugname_lower)
# OG_drug <- drugTarget_all[is.element(drugTarget_all$drugname_lower, OG_druglist),]
OG_drug <- unique(drugTarget_all[is.element(drugTarget_all$genename, glist$$gene_name),])


write.table(OG_drug, paste0(filepath, ".txt"), sep = '\t', quote = F, row.names = F)

OG_druglistFDA <- unique(drugTarget_FDA[is.element(drugTarget_FDA$genename, glist),]$drugname_lower)
OG_drugFDA<- drugTarget_FDA[is.element(drugTarget_FDA$drugname_lower, OG_druglist),]
OG_drugFDA<- unique(drugTarget_FDA[is.element(drugTarget_FDA$genename, glist),])


write.table(OG_drugFDA, paste0(filepath, "_FDA.txt"), sep = '\t', quote = F, row.names = F)

OG_drug2 = unique(OG_drug[c('drugname_lower', 'genename')])
OG_drugFDA2 = unique(OG_drugFDA[c('drugname_lower', 'genename', 'ifCancer_drugbank')])

drugstat <- data.frame(group = c("num_drug",  "num_drugFDA", "num_drugFDAcancer"),
                    value = c(length(unique(OG_drug2$drugname_lower)),
                              length(unique(OG_drugFDA2$drugname_lower)), 
                              length(unique(OG_drugFDA2[which(OG_drugFDA2$ifCancer_drugbank == 'yes'), ]$drugname_lower))))

write.table(drugstat, paste0(filepath, "_stat.txt"), sep = '\t', quote = F, row.names = F)



# 3) TSG-Drug -------------------------------------------------------------------
filepath <- ("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSp_sharedgene/drug_merge3DB/ORAenrich/ITSp_s_vs_TSG_oncoKB")

glist <- unlist(strsplit(res2[res2$ID == "ITSp_vs_TSG_oncoKB",]$geneID, "/"))
glist <- ITS_p_s_notspecific_df[is.element(ITS_p_s_notspecific_df$ENTREZID, glist), ]$gene_name


TSG_druglist <- unique(drugTarget_all[is.element(drugTarget_all$genename, glist$gene_name),]$drugname_lower)
TSG_drug <- drugTarget_all[is.element(drugTarget_all$drugname_lower, TSG_druglist),]

write.table(TSG_drug, paste0(filepath, ".txt"), sep = '\t', quote = F, row.names = F)

# TSG_druglistFDA <- unique(drugTarget_FDA[is.element(drugTarget_FDA$genename, glist),]$drugname_lower)
# TSG_drugFDA<- drugTarget_FDA[is.element(drugTarget_FDA$drugname_lower, TSG_druglist),]
TSG_drugFDA <- unique(drugTarget_FDA[is.element(drugTarget_FDA$genename, glist),])

write.table(TSG_drugFDA, paste0(filepath, "_FDA.txt"), sep = '\t', quote = F, row.names = F)

TSG_drug2 = unique(TSG_drug[c('drugname_lower', 'genename')])
TSG_drugFDA2 = unique(TSG_drugFDA[c('drugname_lower', 'genename', 'ifCancer_drugbank')])

drugstat <- data.frame(group = c("num_drug",  "num_drugFDA", "num_drugFDAcancer"),
                    value = c(length(unique(TSG_drug2$drugname_lower)),
                              length(unique(TSG_drugFDA2$drugname_lower)), 
                              length(unique(TSG_drugFDA2[which(TSG_drugFDA2$ifCancer_drugbank == 'yes'), ]$drugname_lower))))

write.table(drugstat, paste0(filepath, "_stat.txt"), sep = '\t', quote = F, row.names = F)






# ITS_p_r ------------------------------------------------------------------------------
ITS_p_r_notspecific_df <- read.table("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSp_sharedgene/ITS_p_r_notspecific_df.txt", header = T, sep = '\t')
# head(ITS_p_r_notspecific_df)



# enrich custom-GSEA -----------------------------------------------------------

load("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSp_sharedgene/enrich/custom_GSEA_r_scaled.Rdata")
res2 = data.frame(res)

# 1) drugTarget_drugbank-Drug -------------------------------------------------------------------
filepath <- ("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSp_sharedgene/drug_merge3DB/GSEAenrich/ITSp_r_vs_drugtarget")

glist <- unlist(strsplit(res2[res2$ID == "ITSp_vs_drugtarget",]$core_enrichment, "/"))

# DT_druglist <- unique(drugTarget_all[is.element(drugTarget_all$genename, glist),]$drugname_lower)
# DT_drug <- drugTarget_all[is.element(drugTarget_all$drugname_lower, DT_druglist),]

DT_drug <- unique(drugTarget_all[is.element(drugTarget_all$genename, glist),])

write.table(DT_drug, paste0(filepath, ".txt"), sep = '\t', quote = F, row.names = F)

# DT_druglistFDA <- unique(drugTarget_FDA[is.element(drugTarget_FDA$genename, glist),]$drugname_lower)
# DT_drugFDA<- drugTarget_FDA[is.element(drugTarget_FDA$drugname_lower, DT_druglist),]
DT_drugFDA <- unique(drugTarget_FDA[is.element(drugTarget_FDA$genename, glist),])
write.table(DT_drugFDA, paste0(filepath, "_FDA.txt"), sep = '\t', quote = F, row.names = F)

DT_drug2 = unique(DT_drug[c('drugname_lower', 'genename')])
DT_drugFDA2 = unique(DT_drugFDA[c('drugname_lower', 'genename', 'ifCancer_drugbank')])

drugstat <- data.frame(group = c("num_drug",  "num_drugFDA", "num_drugFDAcancer"),
                    value = c(length(unique(DT_drug2$drugname_lower)),
                              length(unique(DT_drugFDA2$drugname_lower)), 
                              length(unique(DT_drugFDA2[which(DT_drugFDA2$ifCancer_drugbank == 'yes'), ]$drugname_lower))))

write.table(drugstat, paste0(filepath, "_stat.txt"), sep = '\t', quote = F, row.names = F)

# unique(DT_drugFDAcancer$name)
# sort(table(DT_drugFDAcancer$genename))


# 2) OG-Drug -------------------------------------------------------------------
filepath <- ("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSp_sharedgene/drug_merge3DB/GSEAenrich/ITSp_r_vs_OG_oncoKB")

glist <- unlist(strsplit(res2[res2$ID == "ITSp_vs_OG_oncoKB",]$core_enrichment, "/"))

OG_druglist <- unique(drugTarget_all[is.element(drugTarget_all$genename, glist),]$drugname_lower)
OG_drug <- drugTarget_all[is.element(drugTarget_all$drugname_lower, OG_druglist),]

write.table(OG_drug, paste0(filepath, ".txt"), sep = '\t', quote = F, row.names = F)

OG_druglistFDA <- unique(drugTarget_FDA[is.element(drugTarget_FDA$genename, glist),]$drugname_lower)
OG_drugFDA<- drugTarget_FDA[is.element(drugTarget_FDA$drugname_lower, OG_druglist),]
OG_drugFDA<- unique(drugTarget_FDA[is.element(drugTarget_FDA$genename, glist),])


write.table(OG_drugFDA, paste0(filepath, "_FDA.txt"), sep = '\t', quote = F, row.names = F)

OG_drug2 = unique(OG_drug[c('drugname_lower', 'genename')])
OG_drugFDA2 = unique(OG_drugFDA[c('drugname_lower', 'genename', 'ifCancer_drugbank')])

drugstat <- data.frame(group = c("num_drug",  "num_drugFDA", "num_drugFDAcancer"),
                    value = c(length(unique(OG_drug2$drugname_lower)),
                              length(unique(OG_drugFDA2$drugname_lower)), 
                              length(unique(OG_drugFDA2[which(OG_drugFDA2$ifCancer_drugbank == 'yes'), ]$drugname_lower))))

write.table(drugstat, paste0(filepath, "_stat.txt"), sep = '\t', quote = F, row.names = F)




# 3) TSG-Drug -------------------------------------------------------------------
filepath <- ("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSp_sharedgene/drug_merge3DB/GSEAenrich/ITSp_r_vs_TSG_oncoKB")

glist <- unlist(strsplit(res2[res2$ID == "ITSp_vs_TSG_oncoKB",]$core_enrichment, "/"))


# TSG_druglist <- unique(drugTarget_all[is.element(drugTarget_all$genename, glist),]$drugname_lower)
# TSG_drug <- drugTarget_all[is.element(drugTarget_all$drugname_lower, TSG_druglist),]
TSG_drug <- unique(drugTarget_all[is.element(drugTarget_all$genename, glist),])


write.table(TSG_drug, paste0(filepath, ".txt"), sep = '\t', quote = F, row.names = F)

# TSG_druglistFDA <- unique(drugTarget_FDA[is.element(drugTarget_FDA$genename, glist),]$drugname_lower)
# TSG_drugFDA<- drugTarget_FDA[is.element(drugTarget_FDA$drugname_lower, TSG_druglist),]
TSG_drugFDA <- unique(drugTarget_FDA[is.element(drugTarget_FDA$genename, glist),])

write.table(TSG_drugFDA, paste0(filepath, "_FDA.txt"), sep = '\t', quote = F, row.names = F)

TSG_drug2 = unique(TSG_drug[c('drugname_lower', 'genename')])
TSG_drugFDA2 = unique(TSG_drugFDA[c('drugname_lower', 'genename', 'ifCancer_drugbank')])

drugstat <- data.frame(group = c("num_drug",  "num_drugFDA", "num_drugFDAcancer"),
                    value = c(length(unique(TSG_drug2$drugname_lower)),
                              length(unique(TSG_drugFDA2$drugname_lower)), 
                              length(unique(TSG_drugFDA2[which(TSG_drugFDA2$ifCancer_drugbank == 'yes'), ]$drugname_lower))))

write.table(drugstat, paste0(filepath, "_stat.txt"), sep = '\t', quote = F, row.names = F)




# enrich custom-ORA ---------------------------------------------------------------
dir.create("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSp_sharedgene/")
dir.create("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSp_sharedgene/drug_merge3DB/")
dir.create("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSp_sharedgene/drug_merge3DB/ORAenrich/")
filepath = "07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSp_sharedgene/enrich/"
load(paste0(filepath, "custom_ORA_s.Rdata"))

res2 = data.frame(res)

# 1) drugTarget_drugbank-Drug -------------------------------------------------------------------
filepath <- ("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSp_sharedgene/drug_merge3DB/ORAenrich/ITS_r_vs_drugtarget")

glist <- unlist(strsplit(res2[res2$ID == "ITSp_vs_drugtarget",]$geneID, "/"))
glist <- ITS_p_r_notspecific_df[is.element(ITS_p_r_notspecific_df$ENTREZID, glist), ]$gene_name

# DT_druglist <- unique(drugTarget_all[is.element(drugTarget_all$genename, glist),]$drugname_lower)
# DT_drug <- drugTarget_all[is.element(drugTarget_all$drugname_lower, DT_druglist),]
DT_drug <- unique(drugTarget_all[is.element(drugTarget_all$genename, glist),])

write.table(DT_drug, paste0(filepath, ".txt"), sep = '\t', quote = F, row.names = F)

# DT_druglistFDA <- unique(drugTarget_FDA[is.element(drugTarget_FDA$genename, glist),]$drugname_lower)
# DT_drugFDA<- drugTarget_FDA[is.element(drugTarget_FDA$drugname_lower, DT_druglist),]
DT_drugFDA <- unique(drugTarget_FDA[is.element(drugTarget_FDA$genename, glist),])
write.table(DT_drugFDA, paste0(filepath, "_FDA.txt"), sep = '\t', quote = F, row.names = F)

DT_drug2 = unique(DT_drug[c('drugname_lower', 'genename')])
DT_drugFDA2 = unique(DT_drugFDA[c('drugname_lower', 'genename', 'ifCancer_drugbank')])

drugstat <- data.frame(group = c("num_drug",  "num_drugFDA", "num_drugFDAcancer"),
                    value = c(length(unique(DT_drug2$drugname_lower)),
                              length(unique(DT_drugFDA2$drugname_lower)), 
                              length(unique(DT_drugFDA2[which(DT_drugFDA2$ifCancer_drugbank == 'yes'), ]$drugname_lower))))

write.table(drugstat, paste0(filepath, "_stat.txt"), sep = '\t', quote = F, row.names = F)


# 2) OG-Drug -------------------------------------------------------------------
filepath <- ("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSp_sharedgene/drug_merge3DB/ORAenrich/ITSp_r_vs_OG_oncoKB")

glist <- unlist(strsplit(res2[res2$ID == "ITSp_vs_OG_oncoKB",]$geneID, "/"))
glist <- ITS_p_r_notspecific_df[is.element(ITS_p_r_notspecific_df$ENTREZID, glist), ]$gene_name

# OG_druglist <- unique(drugTarget_all[is.element(drugTarget_all$genename, glist$$gene_name),]$drugname_lower)
# OG_drug <- drugTarget_all[is.element(drugTarget_all$drugname_lower, OG_druglist),]
OG_drug <- unique(drugTarget_all[is.element(drugTarget_all$genename, glist$$gene_name),])


write.table(OG_drug, paste0(filepath, ".txt"), sep = '\t', quote = F, row.names = F)

OG_druglistFDA <- unique(drugTarget_FDA[is.element(drugTarget_FDA$genename, glist),]$drugname_lower)
OG_drugFDA<- drugTarget_FDA[is.element(drugTarget_FDA$drugname_lower, OG_druglist),]
OG_drugFDA<- unique(drugTarget_FDA[is.element(drugTarget_FDA$genename, glist),])


write.table(OG_drugFDA, paste0(filepath, "_FDA.txt"), sep = '\t', quote = F, row.names = F)

OG_drug2 = unique(OG_drug[c('drugname_lower', 'genename')])
OG_drugFDA2 = unique(OG_drugFDA[c('drugname_lower', 'genename', 'ifCancer_drugbank')])

drugstat <- data.frame(group = c("num_drug",  "num_drugFDA", "num_drugFDAcancer"),
                    value = c(length(unique(OG_drug2$drugname_lower)),
                              length(unique(OG_drugFDA2$drugname_lower)), 
                              length(unique(OG_drugFDA2[which(OG_drugFDA2$ifCancer_drugbank == 'yes'), ]$drugname_lower))))

write.table(drugstat, paste0(filepath, "_stat.txt"), sep = '\t', quote = F, row.names = F)



# 3) TSG-Drug -------------------------------------------------------------------
filepath <- ("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSp_sharedgene/drug_merge3DB/ORAenrich/ITSp_r_vs_TSG_oncoKB")

glist <- unlist(strsplit(res2[res2$ID == "ITSp_vs_TSG_oncoKB",]$geneID, "/"))
glist <- ITS_p_r_notspecific_df[is.element(ITS_p_r_notspecific_df$ENTREZID, glist), ]$gene_name


TSG_druglist <- unique(drugTarget_all[is.element(drugTarget_all$genename, glist$gene_name),]$drugname_lower)
TSG_drug <- drugTarget_all[is.element(drugTarget_all$drugname_lower, TSG_druglist),]

write.table(TSG_drug, paste0(filepath, ".txt"), sep = '\t', quote = F, row.names = F)

# TSG_druglistFDA <- unique(drugTarget_FDA[is.element(drugTarget_FDA$genename, glist),]$drugname_lower)
# TSG_drugFDA<- drugTarget_FDA[is.element(drugTarget_FDA$drugname_lower, TSG_druglist),]
TSG_drugFDA <- unique(drugTarget_FDA[is.element(drugTarget_FDA$genename, glist),])

write.table(TSG_drugFDA, paste0(filepath, "_FDA.txt"), sep = '\t', quote = F, row.names = F)

TSG_drug2 = unique(TSG_drug[c('drugname_lower', 'genename')])
TSG_drugFDA2 = unique(TSG_drugFDA[c('drugname_lower', 'genename', 'ifCancer_drugbank')])

drugstat <- data.frame(group = c("num_drug",  "num_drugFDA", "num_drugFDAcancer"),
                    value = c(length(unique(TSG_drug2$drugname_lower)),
                              length(unique(TSG_drugFDA2$drugname_lower)), 
                              length(unique(TSG_drugFDA2[which(TSG_drugFDA2$ifCancer_drugbank == 'yes'), ]$drugname_lower))))

write.table(drugstat, paste0(filepath, "_stat.txt"), sep = '\t', quote = F, row.names = F)




# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ITS_n   ----------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
dir.create("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSn_sharedgene/drug_merge3DB/")
dir.create("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSn_sharedgene/drug_merge3DB/GSEAenrich/")


# ITS_n_s ------------------------------------------------------------------------------
ITS_n_s_notspecific_df <- read.table("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSn_sharedgene/ITS_n_s_notspecific_df.txt", header = T, sep = '\t')
# head(ITS_n_s_notspecific_df)



# enrich custom-GSEA -----------------------------------------------------------

load("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSn_sharedgene/enrich/custom_GSEA_s_scaled.Rdata")
res2 = data.frame(res)

# 1) drugTarget_drugbank-Drug -------------------------------------------------------------------
filepath <- ("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSn_sharedgene/drug_merge3DB/GSEAenrich/ITSn_s_vs_drugtarget")

glist <- unlist(strsplit(res2[res2$ID == "ITSp_vs_drugtarget",]$core_enrichment, "/"))

# DT_druglist <- unique(drugTarget_all[is.element(drugTarget_all$genename, glist),]$drugname_lower)
# DT_drug <- drugTarget_all[is.element(drugTarget_all$drugname_lower, DT_druglist),]

DT_drug <- unique(drugTarget_all[is.element(drugTarget_all$genename, glist),])

write.table(DT_drug, paste0(filepath, ".txt"), sep = '\t', quote = F, row.names = F)

# DT_druglistFDA <- unique(drugTarget_FDA[is.element(drugTarget_FDA$genename, glist),]$drugname_lower)
# DT_drugFDA<- drugTarget_FDA[is.element(drugTarget_FDA$drugname_lower, DT_druglist),]
DT_drugFDA <- unique(drugTarget_FDA[is.element(drugTarget_FDA$genename, glist),])
write.table(DT_drugFDA, paste0(filepath, "_FDA.txt"), sep = '\t', quote = F, row.names = F)

DT_drug2 = unique(DT_drug[c('drugname_lower', 'genename')])
DT_drugFDA2 = unique(DT_drugFDA[c('drugname_lower', 'genename', 'ifCancer_drugbank')])

drugstat <- data.frame(group = c("num_drug",  "num_drugFDA", "num_drugFDAcancer"),
                    value = c(length(unique(DT_drug2$drugname_lower)),
                              length(unique(DT_drugFDA2$drugname_lower)), 
                              length(unique(DT_drugFDA2[which(DT_drugFDA2$ifCancer_drugbank == 'yes'), ]$drugname_lower))))

write.table(drugstat, paste0(filepath, "_stat.txt"), sep = '\t', quote = F, row.names = F)

# unique(DT_drugFDAcancer$name)
# sort(table(DT_drugFDAcancer$genename))


# 2) OG-Drug -------------------------------------------------------------------
filepath <- ("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSn_sharedgene/drug_merge3DB/GSEAenrich/ITSn_s_vs_OG_oncoKB")

glist <- unlist(strsplit(res2[res2$ID == "ITSp_vs_OG_oncoKB",]$core_enrichment, "/"))

OG_druglist <- unique(drugTarget_all[is.element(drugTarget_all$genename, glist),]$drugname_lower)
OG_drug <- drugTarget_all[is.element(drugTarget_all$drugname_lower, OG_druglist),]

write.table(OG_drug, paste0(filepath, ".txt"), sep = '\t', quote = F, row.names = F)

OG_druglistFDA <- unique(drugTarget_FDA[is.element(drugTarget_FDA$genename, glist),]$drugname_lower)
OG_drugFDA<- drugTarget_FDA[is.element(drugTarget_FDA$drugname_lower, OG_druglist),]
OG_drugFDA<- unique(drugTarget_FDA[is.element(drugTarget_FDA$genename, glist),])


write.table(OG_drugFDA, paste0(filepath, "_FDA.txt"), sep = '\t', quote = F, row.names = F)

OG_drug2 = unique(OG_drug[c('drugname_lower', 'genename')])
OG_drugFDA2 = unique(OG_drugFDA[c('drugname_lower', 'genename', 'ifCancer_drugbank')])

drugstat <- data.frame(group = c("num_drug",  "num_drugFDA", "num_drugFDAcancer"),
                    value = c(length(unique(OG_drug2$drugname_lower)),
                              length(unique(OG_drugFDA2$drugname_lower)), 
                              length(unique(OG_drugFDA2[which(OG_drugFDA2$ifCancer_drugbank == 'yes'), ]$drugname_lower))))

write.table(drugstat, paste0(filepath, "_stat.txt"), sep = '\t', quote = F, row.names = F)




# 3) TSG-Drug -------------------------------------------------------------------
filepath <- ("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSn_sharedgene/drug_merge3DB/GSEAenrich/ITSn_s_vs_TSG_oncoKB")

glist <- unlist(strsplit(res2[res2$ID == "ITSp_vs_TSG_oncoKB",]$core_enrichment, "/"))


# TSG_druglist <- unique(drugTarget_all[is.element(drugTarget_all$genename, glist),]$drugname_lower)
# TSG_drug <- drugTarget_all[is.element(drugTarget_all$drugname_lower, TSG_druglist),]
TSG_drug <- unique(drugTarget_all[is.element(drugTarget_all$genename, glist),])


write.table(TSG_drug, paste0(filepath, ".txt"), sep = '\t', quote = F, row.names = F)

# TSG_druglistFDA <- unique(drugTarget_FDA[is.element(drugTarget_FDA$genename, glist),]$drugname_lower)
# TSG_drugFDA<- drugTarget_FDA[is.element(drugTarget_FDA$drugname_lower, TSG_druglist),]
TSG_drugFDA <- unique(drugTarget_FDA[is.element(drugTarget_FDA$genename, glist),])

write.table(TSG_drugFDA, paste0(filepath, "_FDA.txt"), sep = '\t', quote = F, row.names = F)

TSG_drug2 = unique(TSG_drug[c('drugname_lower', 'genename')])
TSG_drugFDA2 = unique(TSG_drugFDA[c('drugname_lower', 'genename', 'ifCancer_drugbank')])

drugstat <- data.frame(group = c("num_drug",  "num_drugFDA", "num_drugFDAcancer"),
                    value = c(length(unique(TSG_drug2$drugname_lower)),
                              length(unique(TSG_drugFDA2$drugname_lower)), 
                              length(unique(TSG_drugFDA2[which(TSG_drugFDA2$ifCancer_drugbank == 'yes'), ]$drugname_lower))))

write.table(drugstat, paste0(filepath, "_stat.txt"), sep = '\t', quote = F, row.names = F)




# enrich custom-ORA ---------------------------------------------------------------
dir.create("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSn_sharedgene/")
dir.create("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSn_sharedgene/drug_merge3DB/")
dir.create("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSn_sharedgene/drug_merge3DB/ORAenrich/")
filepath = "07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSn_sharedgene/enrich/"
load(paste0(filepath, "custom_ORA_s.Rdata"))

res2 = data.frame(res)

# 1) drugTarget_drugbank-Drug -------------------------------------------------------------------
filepath <- ("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSn_sharedgene/drug_merge3DB/ORAenrich/ITSp_vs_drugtarget")

glist <- unlist(strsplit(res2[res2$ID == "ITSp_vs_drugtarget",]$geneID, "/"))
glist <- ITS_n_s_notspecific_df[is.element(ITS_n_s_notspecific_df$ENTREZID, glist), ]$gene_name

# DT_druglist <- unique(drugTarget_all[is.element(drugTarget_all$genename, glist),]$drugname_lower)
# DT_drug <- drugTarget_all[is.element(drugTarget_all$drugname_lower, DT_druglist),]
DT_drug <- unique(drugTarget_all[is.element(drugTarget_all$genename, glist),])

write.table(DT_drug, paste0(filepath, ".txt"), sep = '\t', quote = F, row.names = F)

# DT_druglistFDA <- unique(drugTarget_FDA[is.element(drugTarget_FDA$genename, glist),]$drugname_lower)
# DT_drugFDA<- drugTarget_FDA[is.element(drugTarget_FDA$drugname_lower, DT_druglist),]
DT_drugFDA <- unique(drugTarget_FDA[is.element(drugTarget_FDA$genename, glist),])
write.table(DT_drugFDA, paste0(filepath, "_FDA.txt"), sep = '\t', quote = F, row.names = F)

DT_drug2 = unique(DT_drug[c('drugname_lower', 'genename')])
DT_drugFDA2 = unique(DT_drugFDA[c('drugname_lower', 'genename', 'ifCancer_drugbank')])

drugstat <- data.frame(group = c("num_drug",  "num_drugFDA", "num_drugFDAcancer"),
                    value = c(length(unique(DT_drug2$drugname_lower)),
                              length(unique(DT_drugFDA2$drugname_lower)), 
                              length(unique(DT_drugFDA2[which(DT_drugFDA2$ifCancer_drugbank == 'yes'), ]$drugname_lower))))

write.table(drugstat, paste0(filepath, "_stat.txt"), sep = '\t', quote = F, row.names = F)


# 2) OG-Drug -------------------------------------------------------------------
filepath <- ("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSn_sharedgene/drug_merge3DB/ORAenrich/ITSn_s_vs_OG_oncoKB")

glist <- unlist(strsplit(res2[res2$ID == "ITSp_vs_OG_oncoKB",]$geneID, "/"))
glist <- ITS_n_s_notspecific_df[is.element(ITS_n_s_notspecific_df$ENTREZID, glist), ]$gene_name

# OG_druglist <- unique(drugTarget_all[is.element(drugTarget_all$genename, glist$$gene_name),]$drugname_lower)
# OG_drug <- drugTarget_all[is.element(drugTarget_all$drugname_lower, OG_druglist),]
OG_drug <- unique(drugTarget_all[is.element(drugTarget_all$genename, glist$$gene_name),])


write.table(OG_drug, paste0(filepath, ".txt"), sep = '\t', quote = F, row.names = F)

OG_druglistFDA <- unique(drugTarget_FDA[is.element(drugTarget_FDA$genename, glist),]$drugname_lower)
OG_drugFDA<- drugTarget_FDA[is.element(drugTarget_FDA$drugname_lower, OG_druglist),]
OG_drugFDA<- unique(drugTarget_FDA[is.element(drugTarget_FDA$genename, glist),])


write.table(OG_drugFDA, paste0(filepath, "_FDA.txt"), sep = '\t', quote = F, row.names = F)

OG_drug2 = unique(OG_drug[c('drugname_lower', 'genename')])
OG_drugFDA2 = unique(OG_drugFDA[c('drugname_lower', 'genename', 'ifCancer_drugbank')])

drugstat <- data.frame(group = c("num_drug",  "num_drugFDA", "num_drugFDAcancer"),
                    value = c(length(unique(OG_drug2$drugname_lower)),
                              length(unique(OG_drugFDA2$drugname_lower)), 
                              length(unique(OG_drugFDA2[which(OG_drugFDA2$ifCancer_drugbank == 'yes'), ]$drugname_lower))))

write.table(drugstat, paste0(filepath, "_stat.txt"), sep = '\t', quote = F, row.names = F)



# 3) TSG-Drug -------------------------------------------------------------------
filepath <- ("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSn_sharedgene/drug_merge3DB/ORAenrich/ITSn_s_vs_TSG_oncoKB")

glist <- unlist(strsplit(res2[res2$ID == "ITSp_vs_TSG_oncoKB",]$geneID, "/"))
glist <- ITS_n_s_notspecific_df[is.element(ITS_n_s_notspecific_df$ENTREZID, glist), ]$gene_name


TSG_druglist <- unique(drugTarget_all[is.element(drugTarget_all$genename, glist$gene_name),]$drugname_lower)
TSG_drug <- drugTarget_all[is.element(drugTarget_all$drugname_lower, TSG_druglist),]

write.table(TSG_drug, paste0(filepath, ".txt"), sep = '\t', quote = F, row.names = F)

# TSG_druglistFDA <- unique(drugTarget_FDA[is.element(drugTarget_FDA$genename, glist),]$drugname_lower)
# TSG_drugFDA<- drugTarget_FDA[is.element(drugTarget_FDA$drugname_lower, TSG_druglist),]
TSG_drugFDA <- unique(drugTarget_FDA[is.element(drugTarget_FDA$genename, glist),])

write.table(TSG_drugFDA, paste0(filepath, "_FDA.txt"), sep = '\t', quote = F, row.names = F)

TSG_drug2 = unique(TSG_drug[c('drugname_lower', 'genename')])
TSG_drugFDA2 = unique(TSG_drugFDA[c('drugname_lower', 'genename', 'ifCancer_drugbank')])

drugstat <- data.frame(group = c("num_drug",  "num_drugFDA", "num_drugFDAcancer"),
                    value = c(length(unique(TSG_drug2$drugname_lower)),
                              length(unique(TSG_drugFDA2$drugname_lower)), 
                              length(unique(TSG_drugFDA2[which(TSG_drugFDA2$ifCancer_drugbank == 'yes'), ]$drugname_lower))))

write.table(drugstat, paste0(filepath, "_stat.txt"), sep = '\t', quote = F, row.names = F)






# ITS_n_r ------------------------------------------------------------------------------
ITS_n_r_notspecific_df <- read.table("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSn_sharedgene/ITS_n_r_notspecific_df.txt", header = T, sep = '\t')
# head(ITS_n_r_notspecific_df)



# enrich custom-GSEA -----------------------------------------------------------

load("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSn_sharedgene/enrich/custom_GSEA_r_scaled.Rdata")
res2 = data.frame(res)

# 1) drugTarget_drugbank-Drug -------------------------------------------------------------------
filepath <- ("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSn_sharedgene/drug_merge3DB/GSEAenrich/ITSn_r_vs_drugtarget")

glist <- unlist(strsplit(res2[res2$ID == "ITSp_vs_drugtarget",]$core_enrichment, "/"))

# DT_druglist <- unique(drugTarget_all[is.element(drugTarget_all$genename, glist),]$drugname_lower)
# DT_drug <- drugTarget_all[is.element(drugTarget_all$drugname_lower, DT_druglist),]

DT_drug <- unique(drugTarget_all[is.element(drugTarget_all$genename, glist),])

write.table(DT_drug, paste0(filepath, ".txt"), sep = '\t', quote = F, row.names = F)

# DT_druglistFDA <- unique(drugTarget_FDA[is.element(drugTarget_FDA$genename, glist),]$drugname_lower)
# DT_drugFDA<- drugTarget_FDA[is.element(drugTarget_FDA$drugname_lower, DT_druglist),]
DT_drugFDA <- unique(drugTarget_FDA[is.element(drugTarget_FDA$genename, glist),])
write.table(DT_drugFDA, paste0(filepath, "_FDA.txt"), sep = '\t', quote = F, row.names = F)

DT_drug2 = unique(DT_drug[c('drugname_lower', 'genename')])
DT_drugFDA2 = unique(DT_drugFDA[c('drugname_lower', 'genename', 'ifCancer_drugbank')])

drugstat <- data.frame(group = c("num_drug",  "num_drugFDA", "num_drugFDAcancer"),
                    value = c(length(unique(DT_drug2$drugname_lower)),
                              length(unique(DT_drugFDA2$drugname_lower)), 
                              length(unique(DT_drugFDA2[which(DT_drugFDA2$ifCancer_drugbank == 'yes'), ]$drugname_lower))))

write.table(drugstat, paste0(filepath, "_stat.txt"), sep = '\t', quote = F, row.names = F)

# unique(DT_drugFDAcancer$name)
# sort(table(DT_drugFDAcancer$genename))


# 2) OG-Drug -------------------------------------------------------------------
filepath <- ("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSn_sharedgene/drug_merge3DB/GSEAenrich/ITSn_r_vs_OG_oncoKB")

glist <- unlist(strsplit(res2[res2$ID == "ITSp_vs_OG_oncoKB",]$core_enrichment, "/"))

OG_druglist <- unique(drugTarget_all[is.element(drugTarget_all$genename, glist),]$drugname_lower)
OG_drug <- drugTarget_all[is.element(drugTarget_all$drugname_lower, OG_druglist),]

write.table(OG_drug, paste0(filepath, ".txt"), sep = '\t', quote = F, row.names = F)

OG_druglistFDA <- unique(drugTarget_FDA[is.element(drugTarget_FDA$genename, glist),]$drugname_lower)
OG_drugFDA<- drugTarget_FDA[is.element(drugTarget_FDA$drugname_lower, OG_druglist),]
OG_drugFDA<- unique(drugTarget_FDA[is.element(drugTarget_FDA$genename, glist),])


write.table(OG_drugFDA, paste0(filepath, "_FDA.txt"), sep = '\t', quote = F, row.names = F)

OG_drug2 = unique(OG_drug[c('drugname_lower', 'genename')])
OG_drugFDA2 = unique(OG_drugFDA[c('drugname_lower', 'genename', 'ifCancer_drugbank')])

drugstat <- data.frame(group = c("num_drug",  "num_drugFDA", "num_drugFDAcancer"),
                    value = c(length(unique(OG_drug2$drugname_lower)),
                              length(unique(OG_drugFDA2$drugname_lower)), 
                              length(unique(OG_drugFDA2[which(OG_drugFDA2$ifCancer_drugbank == 'yes'), ]$drugname_lower))))

write.table(drugstat, paste0(filepath, "_stat.txt"), sep = '\t', quote = F, row.names = F)




# 3) TSG-Drug -------------------------------------------------------------------
filepath <- ("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSn_sharedgene/drug_merge3DB/GSEAenrich/ITSn_r_vs_TSG_oncoKB")

glist <- unlist(strsplit(res2[res2$ID == "ITSp_vs_TSG_oncoKB",]$core_enrichment, "/"))


# TSG_druglist <- unique(drugTarget_all[is.element(drugTarget_all$genename, glist),]$drugname_lower)
# TSG_drug <- drugTarget_all[is.element(drugTarget_all$drugname_lower, TSG_druglist),]
TSG_drug <- unique(drugTarget_all[is.element(drugTarget_all$genename, glist),])


write.table(TSG_drug, paste0(filepath, ".txt"), sep = '\t', quote = F, row.names = F)

# TSG_druglistFDA <- unique(drugTarget_FDA[is.element(drugTarget_FDA$genename, glist),]$drugname_lower)
# TSG_drugFDA<- drugTarget_FDA[is.element(drugTarget_FDA$drugname_lower, TSG_druglist),]
TSG_drugFDA <- unique(drugTarget_FDA[is.element(drugTarget_FDA$genename, glist),])

write.table(TSG_drugFDA, paste0(filepath, "_FDA.txt"), sep = '\t', quote = F, row.names = F)

TSG_drug2 = unique(TSG_drug[c('drugname_lower', 'genename')])
TSG_drugFDA2 = unique(TSG_drugFDA[c('drugname_lower', 'genename', 'ifCancer_drugbank')])

drugstat <- data.frame(group = c("num_drug",  "num_drugFDA", "num_drugFDAcancer"),
                    value = c(length(unique(TSG_drug2$drugname_lower)),
                              length(unique(TSG_drugFDA2$drugname_lower)), 
                              length(unique(TSG_drugFDA2[which(TSG_drugFDA2$ifCancer_drugbank == 'yes'), ]$drugname_lower))))

write.table(drugstat, paste0(filepath, "_stat.txt"), sep = '\t', quote = F, row.names = F)




# enrich custom-ORA ------------------------------------------
# ITSp_vs_drugtarget---------------------
dir.create("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSn_sharedgene/")
dir.create("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSn_sharedgene/drug_merge3DB/")
dir.create("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSn_sharedgene/drug_merge3DB/ORAenrich/")
filepath = "07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSn_sharedgene/enrich/"
load(paste0(filepath, "custom_ORA_r.Rdata"))

res2 = data.frame(res)

# 1) drugTarget_drugbank-Drug -------------------------------------------------------------------
filepath <- ("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSn_sharedgene/drug_merge3DB/ORAenrich/ITSn_r_vs_drugtarget")

glist <- unlist(strsplit(res2[res2$ID == "ITSp_vs_drugtarget",]$geneID, "/"))
glist <- ITS_n_r_notspecific_df[is.element(ITS_n_r_notspecific_df$ENTREZID, glist), ]$gene_name

# DT_druglist <- unique(drugTarget_all[is.element(drugTarget_all$genename, glist),]$drugname_lower)
# DT_drug <- drugTarget_all[is.element(drugTarget_all$drugname_lower, DT_druglist),]
DT_drug <- unique(drugTarget_all[is.element(drugTarget_all$genename, glist),])

write.table(DT_drug, paste0(filepath, ".txt"), sep = '\t', quote = F, row.names = F)

# DT_druglistFDA <- unique(drugTarget_FDA[is.element(drugTarget_FDA$genename, glist),]$drugname_lower)
# DT_drugFDA<- drugTarget_FDA[is.element(drugTarget_FDA$drugname_lower, DT_druglist),]
DT_drugFDA <- unique(drugTarget_FDA[is.element(drugTarget_FDA$genename, glist),])
write.table(DT_drugFDA, paste0(filepath, "_FDA.txt"), sep = '\t', quote = F, row.names = F)

DT_drug2 = unique(DT_drug[c('drugname_lower', 'genename')])
DT_drugFDA2 = unique(DT_drugFDA[c('drugname_lower', 'genename', 'ifCancer_drugbank')])

drugstat <- data.frame(group = c("num_drug",  "num_drugFDA", "num_drugFDAcancer"),
                    value = c(length(unique(DT_drug2$drugname_lower)),
                              length(unique(DT_drugFDA2$drugname_lower)), 
                              length(unique(DT_drugFDA2[which(DT_drugFDA2$ifCancer_drugbank == 'yes'), ]$drugname_lower))))

write.table(drugstat, paste0(filepath, "_stat.txt"), sep = '\t', quote = F, row.names = F)


# 2) OG-Drug -------------------------------------------------------------------
filepath <- ("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSn_sharedgene/drug_merge3DB/ORAenrich/ITSn_r_vs_OG_oncoKB")

glist <- unlist(strsplit(res2[res2$ID == "ITSp_vs_OG_oncoKB",]$geneID, "/"))
glist <- ITS_n_r_notspecific_df[is.element(ITS_n_r_notspecific_df$ENTREZID, glist), ]$gene_name

# OG_druglist <- unique(drugTarget_all[is.element(drugTarget_all$genename, glist$$gene_name),]$drugname_lower)
# OG_drug <- drugTarget_all[is.element(drugTarget_all$drugname_lower, OG_druglist),]
OG_drug <- unique(drugTarget_all[is.element(drugTarget_all$genename, glist$$gene_name),])


write.table(OG_drug, paste0(filepath, ".txt"), sep = '\t', quote = F, row.names = F)

OG_druglistFDA <- unique(drugTarget_FDA[is.element(drugTarget_FDA$genename, glist),]$drugname_lower)
OG_drugFDA<- drugTarget_FDA[is.element(drugTarget_FDA$drugname_lower, OG_druglist),]
OG_drugFDA<- unique(drugTarget_FDA[is.element(drugTarget_FDA$genename, glist),])


write.table(OG_drugFDA, paste0(filepath, "_FDA.txt"), sep = '\t', quote = F, row.names = F)

OG_drug2 = unique(OG_drug[c('drugname_lower', 'genename')])
OG_drugFDA2 = unique(OG_drugFDA[c('drugname_lower', 'genename', 'ifCancer_drugbank')])

drugstat <- data.frame(group = c("num_drug",  "num_drugFDA", "num_drugFDAcancer"),
                    value = c(length(unique(OG_drug2$drugname_lower)),
                              length(unique(OG_drugFDA2$drugname_lower)), 
                              length(unique(OG_drugFDA2[which(OG_drugFDA2$ifCancer_drugbank == 'yes'), ]$drugname_lower))))

write.table(drugstat, paste0(filepath, "_stat.txt"), sep = '\t', quote = F, row.names = F)



# 3) TSG-Drug -------------------------------------------------------------------
filepath <- ("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITScharacteristic_en/ITSn_sharedgene/drug_merge3DB/ORAenrich/ITSn_r_vs_TSG_oncoKB")

glist <- unlist(strsplit(res2[res2$ID == "ITSp_vs_TSG_oncoKB",]$geneID, "/"))
glist <- ITS_n_r_notspecific_df[is.element(ITS_n_r_notspecific_df$ENTREZID, glist), ]$gene_name


TSG_druglist <- unique(drugTarget_all[is.element(drugTarget_all$genename, glist$gene_name),]$drugname_lower)
TSG_drug <- drugTarget_all[is.element(drugTarget_all$drugname_lower, TSG_druglist),]

write.table(TSG_drug, paste0(filepath, ".txt"), sep = '\t', quote = F, row.names = F)

# TSG_druglistFDA <- unique(drugTarget_FDA[is.element(drugTarget_FDA$genename, glist),]$drugname_lower)
# TSG_drugFDA<- drugTarget_FDA[is.element(drugTarget_FDA$drugname_lower, TSG_druglist),]
TSG_drugFDA <- unique(drugTarget_FDA[is.element(drugTarget_FDA$genename, glist),])

write.table(TSG_drugFDA, paste0(filepath, "_FDA.txt"), sep = '\t', quote = F, row.names = F)

TSG_drug2 = unique(TSG_drug[c('drugname_lower', 'genename')])
TSG_drugFDA2 = unique(TSG_drugFDA[c('drugname_lower', 'genename', 'ifCancer_drugbank')])

drugstat <- data.frame(group = c("num_drug",  "num_drugFDA", "num_drugFDAcancer"),
                    value = c(length(unique(TSG_drug2$drugname_lower)),
                              length(unique(TSG_drugFDA2$drugname_lower)), 
                              length(unique(TSG_drugFDA2[which(TSG_drugFDA2$ifCancer_drugbank == 'yes'), ]$drugname_lower))))

write.table(drugstat, paste0(filepath, "_stat.txt"), sep = '\t', quote = F, row.names = F)


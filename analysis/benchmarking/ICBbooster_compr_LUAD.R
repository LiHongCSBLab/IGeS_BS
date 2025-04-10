# rm(list=ls())
work_path <- "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/"
setwd(work_path)
options(stringsAsFactors = F)

library(getopt)

library(dplyr)
library(reshape2)
library(parallel)
library(patchwork)
library(openxlsx)
library(stringr)

# save path setting ------------------------------------------------------------
savepath = "r1_drug_immuneSig_correlation_analysis/07_plotting_v2/ICBbooster_compr/"
dir.create(savepath)


# annotation of drug evidence
drug_proof <- read.table("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/data/FDAapproved_p3Trial_IOtherapy/drug_list_proof_cancertype.txt", header = T, sep = '\t')


# CM-Drug result------------------------------------------------------------------------
res_cmdrug <- read.xlsx("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/07_plotting_v2/CM_drug_res/result_A549/cp/results_with_CM_Drug.xlsx")
res_cmdrug$pert_iname_lowercase = tolower(res_cmdrug$Compound)
res_cmdrug=res_cmdrug[1:20,]
# 70138 ------------------------------------------------------------------------
# IGeS-BS: weighted sum
IGeS_BS1 <- read.csv("r1_drug_immuneSig_correlation_analysis/06_results_summaryandplotting/results_meta_filter_ITSssgsea_CombatRef_17datasets/TUMERIC_m2_meta_en_ACATP0.05/70138/allgenes/LUAD_Primary/drugITSRank/model_weight_weighted_glmnet_weight_res.txt", sep='\t', header=T)
IGeS_BS1$pert_iname_lowercase = tolower(IGeS_BS1$pert_iname)
names(IGeS_BS1)[c(3:6)] <- paste0(names(IGeS_BS1)[c(3:6)], '_IGeS_BS')
IGeS_BS1$dataset = '70138'
# 92742 ------------------------------------------------------------------------
# IGeS-BS: weighted sum
IGeS_BS2 <- read.csv("r1_drug_immuneSig_correlation_analysis/06_results_summaryandplotting/results_meta_filter_ITSssgsea_CombatRef_17datasets/TUMERIC_m2_meta_en_ACATP0.05/92742/allgenes/LUAD_Primary/drugITSRank/model_weight_weighted_glmnet_weight_res.txt", sep='\t', header=T)
IGeS_BS2$pert_iname_lowercase = tolower(IGeS_BS2$pert_iname)
names(IGeS_BS2)[c(3:6)] <- paste0(names(IGeS_BS2)[c(3:6)], '_IGeS_BS')
IGeS_BS2$dataset = '92742'
IGeS_BS = rbind(IGeS_BS1, IGeS_BS2)
res_compr <- merge(IGeS_BS, res_cmdrug, by =  'pert_iname_lowercase' )
res_compr = res_compr[order(res_compr$pert_score_max, decreasing=T),]
res_compr = res_compr[order(res_compr$Times.in.A549, decreasing=T),]
dim(res_compr)
cor(res_compr$score_IGeS_BS, res_compr$pert_score_max)
cor(res_compr$score_IGeS_BS, res_compr$pert_score_max, method='spearman')
head(res_compr[order(res_compr$Times.in.A549, decreasing=T), ],20)
head(res_compr[res_compr$pert_score_max >1, ],20)

res_compr <- merge(res_compr, drug_proof, by=c('drug_index', 'pert_iname'), all.x=T)
write.xlsx(res_compr, paste0(savepath,  "res_compr_LUAD_top20.xlsx"))
res_compr=read.xlsx(paste0(savepath, "res_compr_LUAD.xlsx"))
cor(res_compr$score_IGeS_BS, res_compr$pert_score_max, method='spearman')
cor(res_compr[res_compr$dataset  == 70138, ]$score_IGeS_BS, res_compr[res_compr$dataset  == 70138, ]$pert_score_max, method='spearman')
cor(res_compr[res_compr$dataset  == 92742, ]$score_IGeS_BS, res_compr[res_compr$dataset  == 92742, ]$pert_score_max, method='spearman')

library("ggpubr")
pdf(paste0(savepath, "res_compr_LUAD_CMdrug_top20.pdf"), width=6, height=6)
ggscatter(unique(res_compr[c('score_IGeS_BS', 'pert_score_max')]), 
        x = "score_IGeS_BS", y = "pert_score_max", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "IGeS_BS", ylab = "CM-Drug") +
  ggtitle("all")

ggscatter(unique(res_compr[res_compr$dataset  == 70138, ][c('score_IGeS_BS', 'pert_score_max')]), 
        x = "score_IGeS_BS", y = "pert_score_max", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "IGeS_BS", ylab = "CM-Drug") +
  ggtitle("70138")

ggscatter(unique(res_compr[res_compr$dataset  == 92742, ][c('score_IGeS_BS', 'pert_score_max')]), 
        x = "score_IGeS_BS", y = "pert_score_max", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "IGeS_BS", ylab = "CM-Drug") +
  ggtitle("92742")
dev.off()

df = unique(res_compr[c('score_IGeS_BS', 'pert_score_max')])
nrow(df[df$score_IGeS_BS >0 & df$pert_score_max>1, ])
nrow(df)

tmp=unique(res_compr[which(res_compr$if_with_ICB == 'yes'), ][-c(13:21)])
nrow(tmp[tmp$score_IGeS_BS > 0,])
nrow(tmp[tmp$pert_score_max > 1,])
cor(tmp$score_IGeS_BS, tmp$pert_score_max, method='spearman')

IGeS_BS[IGeS_BS$pert_iname_lowercase == 'camptothecin', ]
tmp=IGeS_BS[IGeS_BS$pert_iname_lowercase %in% tolower(c('mitoxantrone',
'doxorubicin',
'CGP-60474',
'PD-0325901',
'wortmannin',
'withaferin-a',
'dinaciclib',
'daunorubicin',
'idarubicin',
'GSK-429286A',
'AZD-7762',
'KIN001-055',
'crizotinib',
'olaparib',
'TPCA-1',
'AZ-20',
'BRD-K63750851',
'alvocidib',
'GSK-1059615',
'camptothecin')
),]

tmp[order(tmp$pert_iname_lowercase),]
write.xlsx(tmp, paste0(savepath, "res_compr_LUAD_inPaper.xlsx"))



res=read.xlsx(paste0(savepath, "res_compr_LUAD_inPaper.xlsx"))


# rm(list=ls())
workpath <- "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/"
setwd(workpath)

options(stringsAsFactors = F)
library(mlr)
library(dplyr)
library(reshape2)
# rm(list=ls())
work_path <- "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/"
setwd(work_path)
options(stringsAsFactors = F)
library(ggplot2)
library(pheatmap)
library(aplot)
library(reshape2)
library(ggsci)
library(see)



source("03_drug_immuneSig_enrichment/scripts/02_GSEA_geneset_preparation_function.R")
source("03_drug_immuneSig_enrichment/scripts/02_GSEA_geneset_merge_metaAnalysis_function.R")

source("06_results_summaryandplotting/scripts/immunesig_selection_function.R")
source("06_results_summaryandplotting/scripts/ITS_analysis/ITS_function_detection.R")
source("06_results_summaryandplotting/scripts/ITS_analysis/ITSexpr_purity_detection.R")
source("06_results_summaryandplotting/scripts/ITS_analysis/ITSexpr_purity_detection.R")
source("06_results_summaryandplotting/scripts/ITS_analysis/immuneModulatorDetection.R")


meta_result <- read.csv("05_immuneSig_ICBresponse/results/meta_results/res_metaRes.csv")
meta_selected <- meta_result[meta_result$Meta_Pval < 0.05, ]
meta_selected$immune_sig = meta_selected$X


immuneSigInfoPath = "06_results_summaryandplotting/data/immene_cell_hieracy/"
immunesigInfo <- as.data.frame(readxl::read_excel(paste0(immuneSigInfoPath, "immune_sig_hieracy_new.xlsx"), sheet = 1))
immunesigInfo <- immunesigInfo[c(1,4,5:9)]
immunesigInfo$immune_sig = tolower(immunesigInfo$immune_sig)
immunesigInfo = immunesig_name_unify(immunesigInfo)
sigAnnot2 <- immunesigInfo[c("immune_sig","Index","Type","Classification")]
sigAnnot2 <- inner_join(sigAnnot2, meta_selected)
# sigAnnot <- read.table("07_plotting/01_metaAnalysis/immunesig_selected_info.txt",header = T, sep = '\t')
# sigAnnot2 <- sigAnnot[c("immune_sig","upper_OR","Meta_Pval","flag","setindex" ,"Index","Type","Classification")]

# load("05_immuneSig_ICBresponse/results_r0.4_ITSssgsea_Combat/finalmodel_17datasets_ssgseanorm/metaAnalysis/ridge/ridge_model_pred.Rdata")
# ridge_model_para = m_ridge$learner.model$glmnet.fit
# m_ridge_imp <- data.frame(ridge_model_para$beta[,which(ridge_model_para$lambda == m_ridge$learner.model$lambda.min)])
# names(m_ridge_imp) <- 'ridge_weight'
# m_ridge_imp$immune_sig <- rownames(m_ridge_imp)
# m_ridge_imp <- inner_join(m_ridge_imp, sigAnnot2)

# m_ridge_imp_index <- aggregate(m_ridge_imp[c('ridge_weight')], by = list((m_ridge_imp$Index)), mean)
# m_ridge_imp_type <- aggregate(m_ridge_imp[c('ridge_weight')], by = list((m_ridge_imp$Type)), mean)
# m_ridge_imp_index[order(m_ridge_imp_index$ridge_weight),]
# m_ridge_imp_type[order(m_ridge_imp_type$ridge_weight),]



load("05_immuneSig_ICBresponse/results_r0.4_ITSssgsea_Combat/finalmodel_17datasets_ssgseanorm/metaAnalysis/glmnet_ridge/glmnet_ridge_model_pred.Rdata")
glmnet_ridge_imp <- data.frame(ridge.fit$beta[,which(ridge.fit$lambda == cv.ridge$lambda.min)])
names(glmnet_ridge_imp) <- 'ridge_weight'
glmnet_ridge_imp$immune_sig <- rownames(glmnet_ridge_imp)
glmnet_ridge_imp <- inner_join(glmnet_ridge_imp, sigAnnot2)
glmnet_ridge_imp_index <- aggregate(glmnet_ridge_imp[c('ridge_weight')], by = list((glmnet_ridge_imp$Index)), mean)
glmnet_ridge_imp_type <- aggregate(glmnet_ridge_imp[c('ridge_weight')], by = list((glmnet_ridge_imp$Type)), mean)
glmnet_ridge_imp_index[order(glmnet_ridge_imp_index$ridge_weight),]
glmnet_ridge_imp_type[order(glmnet_ridge_imp_type$ridge_weight),]

plot(glmnet_ridge_imp$ridge_weight, glmnet_ridge_imp$OR)

load("05_immuneSig_ICBresponse/results_r0.4_ITSssgsea_Combat/finalmodel_17datasets_ssgseanorm/metaAnalysis/en/elasticnet_model_pred.Rdata")
en_model_para = m_en$learner.model$glmnet.fit
m_en_imp <- data.frame(en_model_para$beta[,which(en_model_para$lambda == m_en$learner.model$lambda.min)])
names(m_en_imp) <- 'en_weight'
m_en_imp$immune_sig <- rownames(m_en_imp)
m_en_imp <- inner_join(m_en_imp, sigAnnot2)
plot(m_en_imp$en_weight, m_en_imp$OR)
abline(v=0,h=1)
m_en_imp <- m_en_imp[m_en_imp$en_weight !=0,]
m_en_imp[m_en_imp$Index == "CD8+ T cell",]
View(m_en_imp[order(m_en_imp$en_weight),])


m_en_imp_index <- aggregate(m_en_imp[c('en_weight')], by = list((m_en_imp$Index)), mean)
m_en_imp_type <- aggregate(m_en_imp[c('en_weight')], by = list((m_en_imp$Type)), mean)
m_en_imp_index[order(m_en_imp_index$en_weight),]
m_en_imp_type[order(m_en_imp_type$en_weight),]



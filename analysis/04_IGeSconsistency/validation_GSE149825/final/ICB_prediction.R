rm(list=ls())
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

# ----------------------------------------------------------------------------
# drug-induced immune signatures NES
cancer = "SKCM"
tumor_subtype = "Metastatic"
purity_method = "TUMERIC"
dataset=70138
datatype="allgenes"
num_gene = 200

savepath = paste0("06_results_summaryandplotting/results_", purity_method,"/")
dir.create(savepath)



# ITSselected <- read.csv("06_results_summaryandplotting/data/immune_sig_selected_new/selected_ITS_0.4.txt", sep = '\t')
meta_res <- read.csv("05_immuneSig_ICBresponse/results/meta_results/res_metaRes.csv")
meta_res_selected <- meta_res[meta_res$Meta_Pval < 0.05, ]
sig_s = meta_res_selected[meta_res_selected$OR > 1,]
sig_r = meta_res_selected[meta_res_selected$OR < 1,]
sig_s = data.frame( immune_sig = sig_s$X, flag = "sensitive", setindex = '01_sensitive')
sig_r = data.frame( immune_sig = sig_r$X, flag = "resistant", setindex = '02_resistant')

load("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/05_immuneSig_ICBresponse/results_r0.4_ITSssgsea_CombatRef/finalmodel_17datasets_ssgseanorm/metaAnalysis_filter/en/elasticnet_model_pred.Rdata")

en_model_para = m_en$learner.model$glmnet.fit
m_en_imp <- data.frame(en_model_para$beta[,which(en_model_para$lambda == m_en$learner.model$lambda.min)])

names(m_en_imp) <- "weight"
m_en_imp$immune_sig = row.names(m_en_imp)


immunesig = rbind(sig_s, sig_r)
immunesig = merge(immunesig, m_en_imp, by = "immune_sig")
immunesig =immunesig[immunesig$weight != 0 ,]


# ------------------------------------------------------------------------------
setwd("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/")
options(stringsAsFactors = F)
genename = read.csv(
  paste0(
    "/picb/bigdata/project/FengFYM/p2_2_liver_cancer_specific_immuno_combo_therapy/",
    "03_mouse_data_analysis/data/mousedata_for_immuClustering_icbResponse/mouse2human_count_manually_genename_merge.txt"
  ),
  sep = "\t",
  header = T
)

# Gene name conversion ---------------------------------------------------------
mat <- read.csv("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/data/immune_combo_dataset/GSE149825/GSE149825_SMAC_ICB_RNASeq_rawcount.csv", row.names = 1)
# metadata <- data.frame(sampleID = names(mat),
#                        # s1795 s2467 s2503 s2521 Vehicle_S1 Vehicle_S2 s2625 s2669
#                        treatment = c("antiPD1","combo","JQ1","antiPD1","vehicle","vehicle","JQ1","combo"))


count_filtered = mat[rowSums(mat > 0) > 0 * ncol(mat), ]


gene_length = read.table(
  paste0(
    "/picb/bigdata/project/FengFYM/s7_annotation_files/All_mm10_gene_len_mus.txt"
  ),
  sep = "\t",
  header = T
)
count_filtered$Gene = rownames(count_filtered)
count_filtered_genelength = inner_join(gene_length, count_filtered, by = "Gene")
genelength_kb = count_filtered_genelength$Length / 1000
rownames(count_filtered_genelength) <- count_filtered_genelength$Gene
CountMat = count_filtered_genelength[-which(is.element(colnames(count_filtered_genelength),
                                                       names(gene_length)))]

rpk <- CountMat / genelength_kb
tpm <- as.data.frame(t(t(rpk) / colSums(rpk) * 10 ^ 6))

tpm$MGI_symbol  = rownames(tpm)
tpm = inner_join(genename, tpm, by = "MGI_symbol")
tpm = tpm[!duplicated(tpm$hgnc_symbol),]
rownames(tpm) = tpm$hgnc_symbol
tpm = tpm[-c(1,2,3)]



# 2) if our model can achieve single sample prediction -------------------------
# Computing simple sample immune signatures
library(GSVA)
load(paste0("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/02_tumorGene_immuneSig_correlation/results_selected/ITS_oe_original_spearman_TUMERIC_r0.4/for_GSVA/pcor/pcor_spearman_",cancer,"_",tumor_subtype,"_positive_200.Rdata"))

# enrichment
ITS_p <- gsva(as.matrix(log2(tpm + 1)),
                genelist_p,
                method = "ssgsea",
                kcdf = "Gaussian",
                min.sz >1,ssgsea.norm = TRUE) 


dir.create("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/result_GSE149825/ICB_prediction/")
write.csv(ITS_p, "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/result_GSE149825/ICB_prediction/ITS_SigScore.csv", quote=F, row.names=F)

# ------------------------------------------------------------------------------
# load model and extract input features
# ------------------------------------------------------------------------------
library(mlr)
library(mlrMBO)

# load EN model ----------------------------------------------------------------
workpath <- "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/"
result_dir <- "05_immuneSig_ICBresponse/results_r0.4_ITSssgsea_CombatRef/"
resultpath = paste0(workpath, result_dir, "finalmodel_17datasets_ssgseanorm/metaAnalysis_filter/")
load(paste0(resultpath, "/en/elasticnet_model_pred.Rdata"))

# extract needed features ------------------------------------------------------
result_dir <- "05_immuneSig_ICBresponse/results_r0.4_ITSssgsea_CombatRef/"
load( paste0(workpath, result_dir, "ITS_afterCombat17datasets_ssgseanorm.Rdata"))

sig_select <- read.csv(paste0(workpath, "07_plotting_v2/01_metaAnalysis/table/metaSelected_immunesigAnnot_subtypedetail_filter.txt"), sep = '\t')
features_selected <- sig_select$immune_sig

# data seperation --------------------------------------------------------------------------

datasets = unique(ITSres_selected$dataset)
# data seperation --------------------------------------------------------------------------

vlSetName <- c("08_Lauss_dataset_outcome", "GSE115821_pre_outcome", "GSE96619_outcome","Liu_dataset_CTLA4Naive_outcome","GSE176307")
dcSetName <- datasets[!is.element(datasets, vlSetName)]

dcSet <- ITSres_selected[is.element(ITSres_selected$dataset, dcSetName),]
vlSet <- ITSres_selected[is.element(ITSres_selected$dataset, vlSetName),]

dcSet_Annot <- dcSet[c('patientID','SampleID', 'cancer','sample_type','dataset')]
dcSet_ITS <- dcSet[-which(is.element(names(dcSet), names(dcSet_Annot)))]

vlSet_Annot <- vlSet[c('patientID','SampleID', 'cancer','sample_type','dataset')]
vlSet_ITS <- vlSet[-which(is.element(names(vlSet), names(vlSet_Annot)))]

# ------------------------------------------------------------------------------
# seperate dataset
# ------------------------------------------------------------------------------
set.seed(1234)
library(sampling)
# library(caret)  
res_df <- dcSet_ITS[!is.na(dcSet_ITS$label), ]
res_df = res_df[order(res_df$label), ]
res_df$label <- factor(res_df$label,levels = c("NR","R"), order=TRUE)
res_df$tag <- 1
res_df[res_df$label == "NR",]$tag = 0
res_df$tag <- factor(res_df$tag,levels = c(0, 1), order=TRUE)

# vlSet_ITS
vl_df <- vlSet_ITS[!is.na(vlSet_ITS$label), ]
vl_df = vl_df[order(vl_df$label), ]
vl_df$label <- factor(vl_df$label,levels = c("NR","R"), order=TRUE)
vl_df$tag <- 1
vl_df[vl_df$label == "NR",]$tag = 0
vl_df$tag <- factor(vl_df$tag,levels = c(0, 1), order=TRUE)
data_train = res_df


data_train_filter <- data_train[c(intersect(features_selected, names(data_train)), 'tag')]
filtered.task = makeClassifTask(id = "ICBresponse", 
                                data = data_train_filter, 
                                target = "tag", 
                                positive = 1)
ITS_p=as.data.frame(t(ITS_p))
ITS_p$sebocytes_xcell=0
ITS_p$tag=1
ITS_p$tag=factor(ITS_p$tag, levels=c(0,1))

data_test <- as.data.frame(ITS_p[c(getTaskFeatureNames(filtered.task), 'tag')])
ICBresponse_testtask = makeClassifTask(id = "ICBresponse", 
                                       data = data_test, 
                                       target = "tag", 
                                       positive = 1)
pred_en = predict(m_en, task = ICBresponse_testtask)
pred_en$data$prob.1
write.csv(pred_en$data, "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/result_GSE149825/ICB_prediction/ITS_prediction_prob.csv", quote=F, row.names=F)


# plot boxplot for the general score -------------------------------------------
# 1) plot for cehicle and 
library(ggplot2)
library(ggpubr)
res=as.data.frame(pred_en$data)
res$sample = rownames(res)
res$group = c(rep('Vehicle', 3), rep('SMAC', 3), rep("ICB", 3), rep('SMAC.ICB', 3))
res$group = factor(res$group, level=c("Vehicle", "SMAC", "ICB", "SMAC.ICB"))
res = res[is.element(res$group, c("Vehicle", "SMAC")),]
p = ggplot(res, aes(group, prob.1, fill = group)) +
        stat_boxplot(geom='errorbar', width=0.15) + 
        geom_boxplot(size=0.5, fill='white', outlier.fill='white', outlier.color='white') +
        geom_jitter(aes(fill=group), width=0.2, shape=21, size=2.5)+
        theme_bw()
ggsave("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/result_GSE149825/ICB_prediction/prediction_2group_boxplot.pdf", p, width=3, height=4)

comparisons <- list(c("Vehicle", "SMAC")) 

p = ggplot(res, aes(group, prob.1, fill = group)) +
        stat_boxplot(geom='errorbar', width=0.15) + 
        geom_boxplot(size=0.5, fill='white', outlier.fill='white', outlier.color='white') +
        geom_jitter(aes(fill=group), width=0.2, shape=21, size=2.5)+
        theme_bw()+
        stat_compare_means(method="t.test",hide.ns = F,
                            comparisons = comparisons,
                            label="p.signif",
                            bracket.size=0.5,
                            size=2.5)
ggsave("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/result_GSE149825/ICB_prediction/prediction_boxplot_t.test_2group.pdf", p, width=5, height=4)


data_test$sample=rownames(data_test)
immuneSig_boxplot = lapply(as.list(names(data_test)[-c(48,58,59)]), function(x){
    # x = "mhc_ips"
    tmp = data_test[c("sample", x)]
    names(tmp) = c("sample", "score")
    tmp$label = "ICB"
    tmp[grep("Vehicle", tmp$sample),]$label = "Vehicle"
    tmp[grep("SMAC_", tmp$sample), ]$label = "SMAC"
    tmp[grep("SMAC.ICB_", tmp$sample), ]$label = "SMAC.ICB"
    tmp$label = factor(tmp$label, c("Vehicle", "SMAC", "ICB", "SMAC.ICB"))
    tmp1 = tmp[grep("Vehicle", tmp$sample),]
    tmp2 = tmp[grep("SMAC_", tmp$sample),]   
    tmp = tmp[is.element(tmp$label, c("Vehicle", "SMAC")), ]
    p = ggplot(tmp, aes(x=label, y=score, fill=label)) +
        stat_boxplot(geom='errorbar', width=0.15) + 
        geom_boxplot(size=0.5, fill='white', outlier.fill='white', outlier.color='white') +
        geom_jitter(aes(fill=label), width=0.2, shape=21, size=2.5)+
        ggtitle(x) + 
        theme_bw()
    # p
    return(p)
})

pdf("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/result_GSE149825/ICB_prediction/singleITS_2group.pdf",
            width =3,
            height = 4, 
            onefile = TRUE)
for(x in seq(length(immuneSig_boxplot))){
        print(immuneSig_boxplot[[x]])
}
dev.off()

immuneSig_boxplot = lapply(as.list(names(data_test)[-c(48,58,59)]), function(x){
    # x = "mhc_ips"
    tmp = data_test[c("sample", x)]
    names(tmp) = c("sample", "score")
    tmp$label = "ICB"
    tmp[grep("Vehicle", tmp$sample),]$label = "Vehicle"
    tmp[grep("SMAC_", tmp$sample), ]$label = "SMAC"
    tmp[grep("SMAC.ICB_", tmp$sample), ]$label = "SMAC.ICB"
    tmp$label = factor(tmp$label, c("Vehicle", "SMAC", "ICB", "SMAC.ICB"))
    tmp1 = tmp[grep("Vehicle", tmp$sample),]
    tmp2 = tmp[grep("SMAC_", tmp$sample),]   
    tmp = tmp[is.element(tmp$label, c("Vehicle", "SMAC")), ]

    p = ggplot(tmp, aes(x=label, y=score, fill=label)) +
        stat_boxplot(geom='errorbar', width=0.15) + 
        geom_boxplot(size=0.5, fill='white', outlier.fill='white', outlier.color='white') +
        geom_jitter(aes(fill=label), width=0.2, shape=21, size=2.5)+
        ggtitle(x) + 
        theme_bw()+
        stat_compare_means(method="t.test",hide.ns = F,
                            comparisons = comparisons,
                            label="p.signif",
                            bracket.size=0.5,
                            size=2.5)
    return(p)
})

pdf("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/result_GSE149825/ICB_prediction/singleITS_t.test_2group.pdf",
            width = 3,
            height = 4, 
            onefile = TRUE)
for(x in seq(length(immuneSig_boxplot))){
        print(immuneSig_boxplot[[x]])
}
dev.off()



# 2) plot for all groups

library(ggplot2)
library(ggpubr)
res=as.data.frame(pred_en$data)
res$sample = rownames(res)
res$group = c(rep('Vehicle', 3), rep('SMAC', 3), rep("ICB", 3), rep('SMAC.ICB', 3))
res$group = factor(res$group, level=c("Vehicle", "SMAC", "ICB", "SMAC.ICB"))
p = ggplot(res, aes(group, prob.1, fill = group)) +
        stat_boxplot(geom='errorbar', width=0.15) + 
        geom_boxplot(size=0.5, fill='white', outlier.fill='white', outlier.color='white') +
        geom_jitter(aes(fill=group), width=0.2, shape=21, size=2.5)+
        theme_bw()
ggsave("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/result_GSE149825/ICB_prediction/prediction_boxplot.pdf", p, width=5, height=4)

my_comparisons1 <- list(c("Vehicle", "SMAC")) 
my_comparisons2 <- list(c("Vehicle", "ICB"))
my_comparisons3 <- list(c("Vehicle", "SMAC.ICB"))
my_comparisons4 <- list(c("SMAC", "ICB"))
my_comparisons5 <- list(c("SMAC", "SMAC.ICB"))
my_comparisons6 <- list(c("ICB", "SMAC.ICB"))
comparisons=c(my_comparisons1,my_comparisons2,my_comparisons3, my_comparisons4,my_comparisons5,my_comparisons6)

p = ggplot(res, aes(group, prob.1, fill = group)) +
        stat_boxplot(geom='errorbar', width=0.15) + 
        geom_boxplot(size=0.5, fill='white', outlier.fill='white', outlier.color='white') +
        geom_jitter(aes(fill=group), width=0.2, shape=21, size=2.5)+
        theme_bw()+
        stat_compare_means(method="t.test",hide.ns = F,
                            comparisons = comparisons,
                            label="p.signif",
                            bracket.size=0.5,
                            size=2.5)
ggsave("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/result_GSE149825/ICB_prediction/prediction_boxplot_t.test.pdf", p, width=5, height=4)


data_test$sample=rownames(data_test)
immuneSig_boxplot = lapply(as.list(names(data_test)[-c(48,58,59)]), function(x){
    # x = "mhc_ips"
    tmp = data_test[c("sample", x)]
    names(tmp) = c("sample", "score")
    tmp$label = "ICB"
    tmp[grep("Vehicle", tmp$sample),]$label = "Vehicle"
    tmp[grep("SMAC_", tmp$sample), ]$label = "SMAC"
    tmp[grep("SMAC.ICB_", tmp$sample), ]$label = "SMAC.ICB"
    tmp$label = factor(tmp$label, c("Vehicle", "SMAC", "ICB", "SMAC.ICB"))
    tmp1 = tmp[grep("Vehicle", tmp$sample),]
    tmp2 = tmp[grep("SMAC_", tmp$sample),]   

    p = ggplot(tmp, aes(x=label, y=score, fill=label)) +
        stat_boxplot(geom='errorbar', width=0.15) + 
        geom_boxplot(size=0.5, fill='white', outlier.fill='white', outlier.color='white') +
        geom_jitter(aes(fill=label), width=0.2, shape=21, size=2.5)+
        ggtitle(x) + 
        theme_bw()
    # p
    return(p)
})

pdf("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/result_GSE149825/ICB_prediction/singleITS.pdf",
            width = 5,
            height = 4, 
            onefile = TRUE)
for(x in seq(length(immuneSig_boxplot))){
        print(immuneSig_boxplot[[x]])
}
dev.off()

immuneSig_boxplot = lapply(as.list(names(data_test)[-c(48,58,59)]), function(x){
    # x = "mhc_ips"
    tmp = data_test[c("sample", x)]
    names(tmp) = c("sample", "score")
    tmp$label = "ICB"
    tmp[grep("Vehicle", tmp$sample),]$label = "Vehicle"
    tmp[grep("SMAC_", tmp$sample), ]$label = "SMAC"
    tmp[grep("SMAC.ICB_", tmp$sample), ]$label = "SMAC.ICB"
    tmp$label = factor(tmp$label, c("Vehicle", "SMAC", "ICB", "SMAC.ICB"))
    tmp1 = tmp[grep("Vehicle", tmp$sample),]
    tmp2 = tmp[grep("SMAC_", tmp$sample),]   

    p = ggplot(tmp, aes(x=label, y=score, fill=label)) +
        stat_boxplot(geom='errorbar', width=0.15) + 
        geom_boxplot(size=0.5, fill='white', outlier.fill='white', outlier.color='white') +
        geom_jitter(aes(fill=label), width=0.2, shape=21, size=2.5)+
        ggtitle(x) + 
        theme_bw()+
        stat_compare_means(method="t.test",hide.ns = F,
                            comparisons = comparisons,
                            label="p.signif",
                            bracket.size=0.5,
                            size=2.5)
    return(p)
})

pdf("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/result_GSE149825/ICB_prediction/singleITS_t.test.pdf",
            width = 5,
            height = 4, 
            onefile = TRUE)
for(x in seq(length(immuneSig_boxplot))){
        print(immuneSig_boxplot[[x]])
}
dev.off()


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
vehicle_tpm = tpm[ grep("Vehicle", colnames(tpm))]
load("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/02_tumorGene_immuneSig_correlation/results_selected/ITS_oe_original_spearman_TUMERIC_r0.4/for_GSVA/pcor/pcor_spearman_SKCM_Metastatic_positive_200.Rdata")

# enrichment
vehicle_p <- gsva(as.matrix(log2(vehicle_tpm + 1)),
                genelist_p,
                method = "ssgsea",
                kcdf = "Gaussian",
                min.sz > 1)
vehicle_p <- as.data.frame(vehicle_p)
vehicle_p$immunesig = rownames(vehicle_p)

# load("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/02_tumorGene_immuneSig_correlation/results_selected/ITS_oe_original_spearman_CPE_r0.4/pcor_spearman_SKCM_Metastatic_negative_200.Rdata")
# # enrichment
# vehicle_n <- gsva(as.matrix(log2(vehicle_tpm + 1)),
#                   genelist_n,
#                   method = "ssgsea",
#                   kcdf = "Gaussian",
#                   min.sz > 1)
# vehicle_n <- as.data.frame(vehicle_n)
# names(vehicle_n) = paste0(names(vehicle_n),"_n")
# vehicle_n$immunesig = rownames(vehicle_n)

vehicle_p <- vehicle_p[c("immunesig", "Vehicle_S1", "Vehicle_S2","Vehicle_S9")]
names(vehicle_p)=c("immune_sig", "Vehicle_S1", "Vehicle_S2","Vehicle_S9")
vehicle_SigScore <- merge(immunesig, vehicle_p, by='immune_sig')
sum(vehicle_SigScore$Vehicle_S1 * vehicle_SigScore$weight)
sum(vehicle_SigScore$Vehicle_S2 * vehicle_SigScore$weight)
sum(vehicle_SigScore$Vehicle_S9 * vehicle_SigScore$weight)
dir.create("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/result_GSE149825/vehicle/")
write.csv(vehicle_SigScore, "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/result_GSE149825/vehicle/vehicle_SigScore.csv", quote=F, row.names=F)




# Single sample drug screening -------------------------------------------------
# drag file to local first
# load("GSE114601_JQ1/LUAD_Primary_drugITS.Rdata")
# mat_s2
# mat_r2
colnames(mat_s2) = gsub("mNES_diff_s.", "", colnames(mat_s2))
mat_s = as.data.frame(t(mat_s2))
vehicle_SigScore_s = vehicle_SigScore[is.element(vehicle_SigScore$immunesig, rownames(mat_s)),]
vehicle_SigScore_s = vehicle_SigScore[vehicle_SigScore$flag == 'sensitive',]


mat_s1 = mat_s[vehicle_SigScore_s$immunesig, ]
rownames(vehicle_SigScore_s) = vehicle_SigScore_s$immunesig
vehicle_SigScore_s = vehicle_SigScore_s[-grep("immunesig",colnames(vehicle_SigScore_s))]

vehicle_pred_s <- lapply(as.list(colnames(vehicle_SigScore_s)), function(s){
  tmp = vehicle_SigScore_s[,s]
  # colSums(tmp * mat_s1)/colSums(mat_s1)
  colSums(tmp + mat_s1)/nrow(mat_s1)
  # tmp_z = (tmp-mean(tmp))/sd(tmp)
  # mat_s1_z = apply(mat_s1, 2, function(x)x-mean(x)/sd(x))
  # colSums(tmp_z + mat_s1_z)/nrow(mat_s1_z)
  # mat_s3 = as.matrix(mat_s1)
  # mat_s3[which(mat_s3 < 0)] = 0
  # colSums(tmp * mat_s3)/apply(mat_s3,2,function(x)length(x[which(x>0)]))

})

names(vehicle_pred_s) = colnames(vehicle_SigScore_s)
vehicle_pred_s = t(do.call(rbind, vehicle_pred_s))
vehicle_pred_s = as.data.frame(vehicle_pred_s)
names(vehicle_pred_s) = paste0(names(vehicle_pred_s), "_s") 
vehicle_pred_s$drug = rownames(vehicle_pred_s)

# vehicle_pred_s[order(vehicle_pred_s$Vehicle_S1, decreasing = T), ]

colnames(mat_r2) = gsub("mNES_diff_r.", "", colnames(mat_r2))
mat_r = as.data.frame(t(mat_r2))
vehicle_SigScore_r = vehicle_SigScore[is.element(vehicle_SigScore$immunesig, rownames(mat_r)),]
mat_r1 = mat_r[vehicle_SigScore_r$immunesig, ]
rownames(vehicle_SigScore_r) = vehicle_SigScore_r$immunesig
vehicle_SigScore_r = vehicle_SigScore_r[-grep("immunesig",colnames(vehicle_SigScore_r))]
vehicle_pred_r <- lapply(as.list(colnames(vehicle_SigScore_r)), function(s){
  tmp = vehicle_SigScore_r[,s]
  # colSums(tmp * mat_r1)/colSums(mat_r1)
  colSums(tmp + mat_r1)/nrow(mat_r1)
  # tmp_z = (tmp-mean(tmp))/sd(tmp)
  # mat_r1_z = apply(mat_r1, 2, function(x)x-mean(x)/sd(x))
  # colSums(tmp_z + mat_r1_z)/nrow(mat_r1_z)

  # mat_r3 = as.matrix(mat_r1)
  # mat_r3[which(mat_r3 > 0)] = 0
  # colSums(tmp * mat_r3)/apply(mat_r3,2,function(x)length(x[which(x <0)]))
  
})
names(vehicle_pred_r) = colnames(vehicle_SigScore_r)
vehicle_pred_r = t(do.call(rbind, vehicle_pred_r))

vehicle_pred_s[rownames(vehicle_pred_s) == "drug1487",]
vehicle_pred_r[rownames(vehicle_pred_r) == "drug1487",]

vehicle_pred_r = as.data.frame(vehicle_pred_r)
# vehicle_pred_r[order(vehicle_pred_r$Vehicle_S1, decreasing = T), ]
names(vehicle_pred_r) = paste0(names(vehicle_pred_r), "_r") 
vehicle_pred_r$drug = rownames(vehicle_pred_r)

# plot
vehicle_SigScore_s$sign = "s"
vehicle_SigScore_r$sign = "r"
vehicle_SigScore_s$x = as.numeric(1:nrow(vehicle_SigScore_s))
vehicle_SigScore_r$x = -1*c(1:nrow(vehicle_SigScore_r))

zeroDot1 = data.frame(Vehicle_S1=0, Vehicle_S2=0, Vehicle_S9=0, sign = "s", x = 0)
rownames(zeroDot1) = 0
zeroDot2 = data.frame(Vehicle_S1=0, Vehicle_S2=0, Vehicle_S9=0, sign = "r", x = 0)
rownames(zeroDot2) = 0

dfplot = rbind(vehicle_SigScore_s, 
               zeroDot1,
               zeroDot2,
               vehicle_SigScore_r)


dfplot$immunesig = rownames(dfplot)
dfplot[dfplot$immunesig=="01",]$immunesig="0"
dfplot$sign = factor(dfplot$sign, levels = c("s","r"))
dfplot$immunesig = factor(dfplot$immunesig)


patient_plot = ggplot(dfplot, aes(x = x, y = Vehicle_S1, color=sign, fill=sign)) + 
  geom_line() + 
  geom_hline(aes(yintercept = 0)) +
  geom_vline(aes(xintercept = 0)) +
  geom_point(size = 3)+
  theme_bw() +
  theme(axis.text.x  =  element_text(angle=-90, vjust=0.5,hjust = 0, size=12),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_blank())+
  scale_x_continuous(breaks = c(seq(min(dfplot$x),0), seq(0,max(dfplot$x))),
                     labels = rev(dfplot$immunesig))
# mat_s1
# mat_r1
drug_s = mat_s1["drug1487"]
drug_r = mat_r1["drug1487"]
drug_s$sign = "s"
drug_r$sign = "r"
drug_s$x = as.numeric(1:nrow(drug_s))
drug_r$x = -1*c(1:nrow(drug_r))

zeroDot1 = data.frame(drug1487=0, sign = "s", x = 0)
rownames(zeroDot1) = 0
zeroDot2 = data.frame(drug1487=0, sign = "r", x = 0)
rownames(zeroDot2) = 0

drug_dfplot = rbind(drug_s, 
               zeroDot1,
               zeroDot2,
               drug_r)


drug_dfplot$immunesig = rownames(drug_dfplot)
drug_dfplot[drug_dfplot$immunesig=="01",]$immunesig="0"
drug_dfplot$sign = factor(drug_dfplot$sign, levels = c("s","r"))
drug_dfplot$immunesig = factor(drug_dfplot$immunesig)


drug_plot = ggplot(drug_dfplot, aes(x = x, y = drug1487, color=sign, fill=sign)) + 
  geom_line() + 
  geom_hline(aes(yintercept = 0)) +
  geom_vline(aes(xintercept = 0)) +
  geom_point(size = 3)+
  theme_bw() +
  theme(axis.text.x  =  element_text(angle=-90, vjust=0.5,hjust = 0, size=12),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_blank())+
  scale_x_continuous(breaks = c(seq(min(drug_dfplot$x),0), seq(0,max(drug_dfplot$x))),
                     labels = rev(drug_dfplot$immunesig))
drug_plot

dfplot2 = dfplot[-c(2,3)]
dfplot2$label = "patient"
names(dfplot2) = c("value","sign","x","immunesig","label")
drug_dfplot$label = "drug1487"
names(drug_dfplot) = c("value","sign","x","immunesig","label")

dftmp = rbind(dfplot2, drug_dfplot)
dftmp$label = factor(dftmp$label)

p = ggplot(dftmp, aes(x = x, y = value, color=label, fill=label)) + 
  geom_line() + 
  geom_hline(aes(yintercept = 0)) +
  geom_vline(aes(xintercept = 0)) +
  geom_point(size = 3)+
  theme_bw() +
  theme(axis.text.x  =  element_text(angle=-90, vjust=0.5,hjust = 0, size=12),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_blank())+
  scale_x_continuous(breaks = c(seq(min(drug_dfplot$x),0), seq(0,max(drug_dfplot$x))),
                     labels = rev(drug_dfplot$immunesig))
p
ggsave("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/data/immune_combo_dataset/GSE149825/Vehicle_S1.pdf",
p, height = 10, width=25 )

dfplot2 = dfplot[-c(1,3)]
dfplot2$label = "patient"
names(dfplot2) = c("value","sign","x","immunesig","label")
drug_dfplot$label = "drug1487"
names(drug_dfplot) = c("value","sign","x","immunesig","label")

dftmp = rbind(dfplot2, drug_dfplot)
dftmp$label = factor(dftmp$label)

p2 = ggplot(dftmp, aes(x = x, y = value, color=label, fill=label)) + 
  geom_line() + 
  geom_hline(aes(yintercept = 0)) +
  geom_vline(aes(xintercept = 0)) +
  geom_point(size = 3)+
  theme_bw() +
  theme(axis.text.x  =  element_text(angle=-90, vjust=0.5,hjust = 0, size=12),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_blank())+
  scale_x_continuous(breaks = c(seq(min(drug_dfplot$x),0), seq(0,max(drug_dfplot$x))),
                     labels = rev(drug_dfplot$immunesig))
p2
ggsave("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/data/immune_combo_dataset/GSE149825/Vehicle_S2.pdf",
        p2, height = 10, width=25 )


dfplot2 = dfplot[-c(1,2)]
dfplot2$label = "patient"
names(dfplot2) = c("value","sign","x","immunesig","label")
drug_dfplot$label = "drug1487"
names(drug_dfplot) = c("value","sign","x","immunesig","label")

dftmp = rbind(dfplot2, drug_dfplot)
dftmp$label = factor(dftmp$label)

p3 = ggplot(dftmp, aes(x = x, y = value, color=label, fill=label)) + 
  geom_line() + 
  geom_hline(aes(yintercept = 0)) +
  geom_vline(aes(xintercept = 0)) +
  geom_point(size = 3)+
  theme_bw() +
  theme(axis.text.x  =  element_text(angle=-90, vjust=0.5,hjust = 0, size=12),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_blank())+
  scale_x_continuous(breaks = c(seq(min(drug_dfplot$x),0), seq(0,max(drug_dfplot$x))),
                     labels = rev(drug_dfplot$immunesig))
p3
ggsave("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/data/immune_combo_dataset/GSE149825/Vehicle_S9.pdf",
        p3, height = 10, width=25 )



vehicle_pred = inner_join(vehicle_pred_s, vehicle_pred_r)
head(vehicle_pred)
vehicle_pred$Vehicle_S1 = vehicle_pred$Vehicle_S1_s - vehicle_pred$Vehicle_S1_r 
vehicle_pred$Vehicle_S2 = vehicle_pred$Vehicle_S2_s - vehicle_pred$Vehicle_S2_r 
vehicle_pred$Vehicle_S9 = vehicle_pred$Vehicle_S9_s - vehicle_pred$Vehicle_S9_r 
vehicle_pred$score = (vehicle_pred$Vehicle_S1 + vehicle_pred$Vehicle_S2 + vehicle_pred$Vehicle_S9)/3
# vehicle_pred = vehicle_pred[c("drug", "Vehicle_S1", "Vehicle_S2")]
head(vehicle_pred[order(vehicle_pred$score, decreasing = T), ],10)
head(vehicle_pred[order(vehicle_pred$Vehicle_S1, decreasing = T), ],10)
head(vehicle_pred[order(vehicle_pred$Vehicle_S2, decreasing = T), ],10)
head(vehicle_pred[order(vehicle_pred$Vehicle_S9, decreasing = T), ],10)

save(vehicle_pred, vehicle_pred_s, vehicle_pred_r, file = "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/data/immune_combo_dataset/GSE149825/prediction_result.Rdata")
# 

vehicle_pred = vehicle_pred[order(vehicle_pred$score, decreasing = T),]
vehicle_pred$rank = seq(nrow(vehicle_pred))
vehicle_pred[vehicle_pred$drug=="drug1487",]

# 
# vehicle_pred$Vehicle_S1 = (vehicle_pred$Vehicle_S1_s - mean(vehicle_pred$Vehicle_S1_s))/sd(vehicle_pred$Vehicle_S1_s) - (vehicle_pred$Vehicle_S1_r - mean(vehicle_pred$Vehicle_S1_r))/sd(vehicle_pred$Vehicle_S1_r)
# vehicle_pred$Vehicle_S2 = (vehicle_pred$Vehicle_S2_s - mean(vehicle_pred$Vehicle_S2_s))/sd(vehicle_pred$Vehicle_S2_s) - (vehicle_pred$Vehicle_S2_r - mean(vehicle_pred$Vehicle_S2_r))/sd(vehicle_pred$Vehicle_S2_r)
# # vehicle_pred = vehicle_pred[c("drug", "Vehicle_S1", "Vehicle_S2")]
# vehicle_pred[vehicle_pred$drug=="drug222",]
# vehicle_pred[vehicle_pred$drug=="drug19",]
# head(vehicle_pred[order(vehicle_pred$Vehicle_S1, decreasing = T), ],10)
# head(vehicle_pred[order(vehicle_pred$Vehicle_S2, decreasing = T), ],10)
# 
# 

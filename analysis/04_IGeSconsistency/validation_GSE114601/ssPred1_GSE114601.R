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
cancer = "LUAD"
tumor_subtype = "Primary"
purity_method = "TUMERIC"
dataset=70138
datatype="allgenes"
num_gene = 200

savepath = paste0("06_results_summaryandplotting/results_", purity_method,"/")
dir.create(savepath)
immunesig <- immune_sig_filter(cancer = cancer,
                              tumor_subtype = tumor_subtype,
                              purity_method = purity_method,
                              datatype="allgenes",
                              num_gene = 200,
                              auc_threshold = 0.6,
                              proof = 1, # how many datasets support the immune sigture to have prediction power
                              savepath = savepath)
# immunesig <- immunesig_extraction(cancer = cancer,
#                                   tumor_subtype = tumor_subtype)

positive_df = df_for_plot(cancer = cancer, 
                          tumor_subtype = tumor_subtype,
                          purity_method = purity_method,
                          dataset=dataset, 
                          sign = "positive", 
                          num_gene = 200, 
                          datatype = datatype, 
                          immunesig = immunesig, 
                          # resPath = resPath, 
                          work_path = work_path,
                          savepath = savepath)
negative_df = df_for_plot(cancer = cancer, 
                          tumor_subtype = tumor_subtype,
                          purity_method = purity_method,
                          dataset=dataset, 
                          sign = "negative", 
                          num_gene = 200,  
                          datatype = datatype,  
                          immunesig = immunesig,
                          # resPath = resPath, 
                          work_path = work_path,
                          savepath = savepath)


names(positive_df) = c("drug_index","immune_sig", "mNES_p","mPvalueNaive_p",
                       "mPvalue_fisher_p","mPvalue_ACAT_p")
names(negative_df) = c("drug_index","immune_sig", "mNES_n","mPvalueNaive_n",
                       "mPvalue_fisher_n","mPvalue_ACAT_n")
df = merge(positive_df, negative_df, by = c("drug_index", "immune_sig"), all=T)
df[is.na(df$mNES_p),]$mNES_p=0
df[is.na(df$mNES_n),]$mNES_n=0

df$mNES_diff = df$mNES_p - df$mNES_n
# df = df[which(df$mPvalue_ACAT_p < 0.05 ),]
# df = df[which(df$mPvalue_ACAT_p < 0.05 & df$mPvalue_ACAT_n < 0.05),]
# df = df[-grep("oe_",df$immune_sig),]

df_immunesig = merge(df, immunesig, by = "immune_sig")
if(length(unique(df_immunesig$immune_sig)) < 0.1*length(unique(immunesig$immune_sig))){
  print("analysis only with positive gene set")
  df_immunesig = merge(positive_df, immunesig, by = "immune_sig")
  df_immunesig$mNES_diff = df_immunesig$mNES_p - 0
}



drugTargetMoA <- read.table("06_results_summaryandplotting/data/drugTargetMoA_merged.txt", sep = "\t", header = T)
drugindexall <- read.table("06_results_summaryandplotting/data/drug_index_LINCS.txt", sep = "\t", header = T)

df_s = df_immunesig[grep("01",df_immunesig$setindex),]
df_r = df_immunesig[grep("02",df_immunesig$setindex),]

df_s = df_s[c("drug_index","immune_sig","mNES_diff")]
names(df_s) = c(c("drug_index","immune_sig_s","mNES_diff_s"))
df_r = df_r[c("drug_index","immune_sig","mNES_diff")]
names(df_r) = c(c("drug_index","immune_sig_r","mNES_diff_r"))


df = merge(unique(df_s[c("drug_index","mNES_diff_s")]),
           unique(df_r[c("drug_index","mNES_diff_r")]),
           by= "drug_index")

mat_s <-  reshape(df_s,  timevar = "immune_sig_s",
                  idvar = "drug_index",direction = "wide")
rownames(mat_s) <- mat_s$drug_index
mat_s <- mat_s[-1]

mat_r <-  reshape(df_r,  timevar = "immune_sig_r",
                  idvar = "drug_index",direction = "wide")
rownames(mat_r) <- mat_r$drug_index
mat_r <- mat_r[-1]

mat_s2 <- as.matrix(mat_s[is.element(rownames(mat_s), intersect(rownames(mat_s), rownames(mat_r) )),])
mat_r2 <- as.matrix(mat_r[is.element(rownames(mat_r), intersect(rownames(mat_s), rownames(mat_r) )),])

mat_s2[is.na(mat_s2)] = 0
mat_r2[is.na(mat_r2)] = 0

save(mat_s2, mat_r2, file = paste0("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/data/drugITS/all_genes/", dataset, "/", cancer, "_", tumor_subtype, "_drugITS.Rdata"))


# ------------------------------------------------------------------------------
setwd("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/")
options(stringsAsFactors = F)
genename = read.csv(
  paste0(
    '/picb/bigdata/project/FengFYM/p2_2_liver_cancer_specific_immuno_combo_therapy/',
    "03_mouse_data_analysis/data/mousedata_for_immuClustering_icbResponse/mouse2human_count_manually_genename_merge.txt"
  ),
  sep = '\t',
  header = T
)

# Gene name conversion ---------------------------------------------------------
mat <- read.csv("data/immune_combo_dataset/GSE114601/counts.raw.csv", row.names = 1)
metadata <- data.frame(sampleID = names(mat),
                       # s1795 s2467 s2503 s2521 s2596 s2617 s2625 s2669
                       treatment = c("antiPD1","combo","JQ1","antiPD1","vehicle","vehicle","JQ1","combo"))


count_filtered = mat[rowSums(mat > 0) > 0 * ncol(mat), ]


gene_length = read.table(
  paste0(
    "/picb/bigdata/project/FengFYM/s7_annotation_files/All_mm10_gene_len_mus.txt"
  ),
  sep = '\t',
  header = T
)
count_filtered$Gene = rownames(count_filtered)
count_filtered_genelength = inner_join(gene_length, count_filtered, by = 'Gene')
genelength_kb = count_filtered_genelength$Length / 1000
rownames(count_filtered_genelength) <- count_filtered_genelength$Gene
CountMat = count_filtered_genelength[-which(is.element(colnames(count_filtered_genelength),
                                                       names(gene_length)))]

rpk <- CountMat / genelength_kb
tpm <- as.data.frame(t(t(rpk) / colSums(rpk) * 10 ^ 6))

tpm$MGI_symbol  = rownames(tpm)
tpm = inner_join(genename, tpm, by = 'MGI_symbol')
tpm = tpm[!duplicated(tpm$hgnc_symbol),]
rownames(tpm) = tpm$hgnc_symbol
tpm = tpm[-c(1,2,3)]

# generating expression files for GSEA

mat_gsea <- cbind(data.frame(NAME = rownames(tpm)), 
                 data.frame(DESCRIPTION = rep(NA, nrow(tpm))),
                 log2(1+tpm[which(is.element(names(tpm), metadata[metadata$treatment == 'JQ1',]$sampleID))]),
                 log2(1+tpm[which(is.element(names(tpm), metadata[metadata$treatment == 'vehicle',]$sampleID))]))
mat_gsea <- mat_gsea[-which(rowSums(mat_gsea[c(3:6)])==0),]

write.table(mat_gsea, "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/data/immune_combo_dataset/GSE114601/GSE114601RNASeq_gsea.txt",
            row.names = F, quote = F, sep = "\t")

# generating class file for GSEA
grp <- list(paste(4, 2, 1, collapse = " "),
            paste("#", "treated", "vehicle", collapse = " "),
            paste(c(rep("treated", 2),
                    rep("vehicle", 2)),
                  collapse = " "))
grp <- do.call(rbind, grp)
write.table(grp, "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/data/immune_combo_dataset/GSE114601/GSE114601RNASeq_gsea.cls", 
            col.names = F, row.names = F, quote = F, sep = "\t")


# 1) if drug effect on cell lines and mice models are consistant -------------------------
gsea_resultpath = "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/result_GSE114601/drug_immunesig_TUMERIC_GSEA_result/70138/LUAD_Primary/"
res_p_path <- list.files(paste0(paste0(gsea_resultpath, "positive_200/")))
res_p_file <- list.files(paste0(gsea_resultpath, "positive_200/", res_p_path))
df_positive_treated <- read.csv(paste0(gsea_resultpath, "positive_200/", res_p_path, "/gsea_report_for_treated_1636621277738.xls"), 
                        sep = "\t", stringsAsFactors = F)
df_positive_treated = df_positive_treated[c(1, 5:9)]
df_positive_treated$immune_sig <- tolower(df_positive_treated$NAME)
df_positive_vehicle <- read.csv(paste0(gsea_resultpath, "positive_200/", res_p_path, "/gsea_report_for_vehicle_1636621277738.xls"), 
                        sep = "\t", stringsAsFactors = F)
df_positive_vehicle = df_positive_vehicle[c(1, 5:9)]
df_positive_vehicle$immune_sig <- tolower(df_positive_vehicle$NAME)
df_positive = rbind(df_positive_treated, df_positive_vehicle)
df_positive = df_positive[df_positive$FDR.q.val < 0.05, ]
df_positive = df_positive[c("immune_sig", "NES", "FDR.q.val")]
names(df_positive) = c("immune_sig", "NES_p", "FDR_p")

res_n_path <- list.files(paste0(paste0(gsea_resultpath, "negative_200/")))
res_n_file <- list.files(paste0(gsea_resultpath, "negative_200/", res_n_path))
df_negative_treated <- read.csv(paste0(gsea_resultpath, "negative_200/", res_n_path, "/gsea_report_for_treated_1636621278104.xls"), 
                        sep = "\t", stringsAsFactors = F)
df_negative_treated = df_negative_treated[c(1, 5:9)]
df_negative_treated$immune_sig <- tolower(df_negative_treated$NAME)
df_negative_vehicle <- read.csv(paste0(gsea_resultpath, "negative_200/", res_n_path, "/gsea_report_for_vehicle_1636621278104.xls"), 
                        sep = "\t", stringsAsFactors = F)
df_negative_vehicle = df_negative_vehicle[c(1, 5:9)]
df_negative_vehicle$immune_sig <- tolower(df_negative_vehicle$NAME)
df_negative = rbind(df_negative_treated, df_negative_vehicle)
df_negative = df_negative[df_negative$FDR.q.val < 0.05, ]
df_negative = df_negative[c("immune_sig", "NES", "FDR.q.val")]
names(df_negative) = c("immune_sig", "NES_n", "FDR_n")

df_gsea = merge(df_positive, df_negative, by = "immune_sig", all = T)
df_gsea[is.na(df_gsea$NES_p),]$NES_p=0
df_gsea[is.na(df_gsea$NES_n),]$NES_n=0

df_gsea$NES_diff = df_gsea$NES_p - df_gsea$NES_n

df_gsea_immunesig = merge(df_gsea, immunesig, by = "immune_sig")
if(length(unique(df_gsea_immunesig$immune_sig)) < 0.1 * length(unique(df_gsea_immunesig$immune_sig))){
  print("analysis only with positive gene set")
  df_gsea_immunesig = merge(df_positive, immunesig, by = "immune_sig")
  df_gsea_immunesig$NES_diff = df_gsea_immunesig$NES_p - 0
}



drug_s = data.frame(drug19 = mat_s2['drug19',])
drug_r = data.frame(drug19 = mat_r2['drug19',])
drug_s$sign = 's'
drug_r$sign = 'r'
drug_s$x = as.numeric(1:nrow(drug_s))
drug_r$x = -1*c(1:nrow(drug_r))

zeroDot1 = data.frame(drug19=0, sign = 's', x = 0)
rownames(zeroDot1) = 0
zeroDot2 = data.frame(drug19=0, sign = 'r', x = 0)
rownames(zeroDot2) = 0

drug_dfplot = rbind(drug_s, 
                    zeroDot1,
                    zeroDot2,
                    drug_r)


drug_dfplot$immunesig = gsub("mNES_diff_s.", "", gsub("mNES_diff_r.", "", rownames(drug_dfplot)))

drug_dfplot[drug_dfplot$immunesig=="01",]$immunesig="0"
drug_dfplot$sign = factor(drug_dfplot$sign, levels = c('s','r'))
drug_dfplot$immunesig = factor(drug_dfplot$immunesig)

df_gsea_immunesig2 = df_gsea_immunesig[c("immune_sig", "NES_diff")]
names(df_gsea_immunesig2) = c("immunesig", "NES_diff")
df_gsea_immunesig2 = inner_join(df_gsea_immunesig2, drug_dfplot)
df_gsea_immunesig2[df_gsea_immunesig2$x <0,]$x = -1*seq(nrow(df_gsea_immunesig2[df_gsea_immunesig2$x < 0,]))
df_gsea_immunesig2[df_gsea_immunesig2$x >0,]$x = seq(nrow(df_gsea_immunesig2[df_gsea_immunesig2$x > 0,]))

dfplot2 = df_gsea_immunesig2[c(1,2,4,5)]
dfplot2$label = 'JQ1_treated_mouse'
names(dfplot2) = c('immunesig','value','sign','x','label')

drug_dfplot = df_gsea_immunesig2[c(1,3,4,5)]
drug_dfplot$label = 'drug19'
names(drug_dfplot) = c('immunesig','value','sign','x','label')

dftmp = rbind(dfplot2, drug_dfplot)
dftmp = dftmp[order(dftmp$x), ]
dftmp$label = factor(dftmp$label)
dftmp[is.na(dftmp$value),]$value = 0

p = ggplot(dftmp, aes(x = x, y = value, color=label, fill=label)) + 
  geom_line() + 
  geom_hline(aes(yintercept = 0)) +
  geom_vline(aes(xintercept = 0)) +
  geom_point(size = 3)+
  theme_bw() +
  theme(axis.text.x  =  element_text(angle=-90, vjust=0.5,hjust = 0, size=12),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_blank())+
  scale_x_continuous(breaks = unique(dftmp$x),
                     labels = unique(dftmp$immunesig))
p
ggsave("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/result_GSE114601/CellvstreatedMouse_JQ1_drug19.pdf",
p, height = 10, width=25 )

write.csv(df_gsea_immunesig2, "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/result_GSE114601/CellvstreatedMouse_JQ1_drug19.csv",
          row.names = F, quote = F)



drug_s = as.data.frame(t(mat_s['drug222',]))
drug_r = as.data.frame(t(mat_r['drug222',]))
drug_s$sign = 's'
drug_r$sign = 'r'
drug_s$x = as.numeric(1:nrow(drug_s))
drug_r$x = -1*c(1:nrow(drug_r))

zeroDot1 = data.frame(drug222=0, sign = 's', x = 0)
rownames(zeroDot1) = 0
zeroDot2 = data.frame(drug222=0, sign = 'r', x = 0)
rownames(zeroDot2) = 0

drug_dfplot = rbind(drug_s, 
               zeroDot1,
               zeroDot2,
               drug_r)


drug_dfplot$immunesig = gsub("mNES_diff_s.", "", gsub("mNES_diff_r.", "", rownames(drug_dfplot)))

drug_dfplot[drug_dfplot$immunesig=="01",]$immunesig="0"
drug_dfplot$sign = factor(drug_dfplot$sign, levels = c('s','r'))
drug_dfplot$immunesig = factor(drug_dfplot$immunesig)

df_gsea_immunesig2 = df_gsea_immunesig[c("immune_sig", "NES_diff")]
names(df_gsea_immunesig2) = c("immunesig", "NES_diff")
df_gsea_immunesig2 = inner_join(df_gsea_immunesig2, drug_dfplot)
df_gsea_immunesig2[df_gsea_immunesig2$x <0,]$x = -1*seq(nrow(df_gsea_immunesig2[df_gsea_immunesig2$x < 0,]))
df_gsea_immunesig2[df_gsea_immunesig2$x >0,]$x = seq(nrow(df_gsea_immunesig2[df_gsea_immunesig2$x > 0,]))

dfplot2 = df_gsea_immunesig2[c(1,2,4,5)]
dfplot2$label = 'JQ1_treated_mouse'
names(dfplot2) = c('immunesig','value','sign','x','label')

drug_dfplot = df_gsea_immunesig2[c(1,3,4,5)]
drug_dfplot$label = 'drug222'
names(drug_dfplot) = c('immunesig','value','sign','x','label')

dftmp = rbind(dfplot2, drug_dfplot)
dftmp = dftmp[order(dftmp$x), ]
dftmp$label = factor(dftmp$label)
dftmp[is.na(dftmp$value),]$value = 0

p = ggplot(dftmp, aes(x = x, y = value, color=label, fill=label)) + 
  geom_line() + 
  geom_hline(aes(yintercept = 0)) +
  geom_vline(aes(xintercept = 0)) +
  geom_point(size = 3)+
  theme_bw() +
  theme(axis.text.x  =  element_text(angle=-90, vjust=0.5,hjust = 0, size=12),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_blank())+
  scale_x_continuous(breaks = unique(dftmp$x),
                     labels = unique(dftmp$immunesig))
p
ggsave("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/result_GSE114601/CellvstreatedMouse_JQ1_drug222.pdf",
p, height = 10, width=25 )

write.csv(df_gsea_immunesig2, "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/result_GSE114601/CellvstreatedMouse_JQ1_drug222.csv",
          row.names = F, quote = F)


# extract JQ1+PD1 --------------------------------------------------------------

mat_gsea <- cbind(data.frame(NAME = rownames(tpm)), 
                 data.frame(DESCRIPTION = rep(NA, nrow(tpm))),
                 log2(1+tpm[which(is.element(names(tpm), metadata[metadata$treatment == 'combo',]$sampleID))]),
                 log2(1+tpm[which(is.element(names(tpm), metadata[metadata$treatment == 'vehicle',]$sampleID))]))
mat_gsea <- mat_gsea[-which(rowSums(mat_gsea[c(3:6)])==0),]

write.table(mat_gsea, "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/data/immune_combo_dataset/GSE114601/GSE114601RNASeq_COMBO_gsea.txt",
            row.names = F, quote = F, sep = "\t")

# generating class file for GSEA
grp <- list(paste(4, 2, 1, collapse = " "),
            paste("#", "treated", "vehicle", collapse = " "),
            paste(c(rep("treated", 2),
                    rep("vehicle", 2)),
                  collapse = " "))
grp <- do.call(rbind, grp)
write.table(grp, "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/data/immune_combo_dataset/GSE114601/GSE114601RNASeq_COMBO_gsea.cls", 
            col.names = F, row.names = F, quote = F, sep = "\t")

# if drug effect on cell lines and mice models are consistant -------------------------
gsea_resultpath = "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/result_GSE114601/drug_immunesig_TUMERIC_GSEA_result_COMBO/70138/LUAD_Primary/"
res_p_path <- list.files(paste0(paste0(gsea_resultpath, "positive_200/")))
res_p_file <- list.files(paste0(gsea_resultpath, "positive_200/", res_p_path))
df_positive_treated <- read.csv(paste0(gsea_resultpath, "positive_200/", res_p_path, "/gsea_report_for_treated_1636630344414.xls"), 
                        sep = "\t", stringsAsFactors = F)
df_positive_treated = df_positive_treated[c(1, 5:9)]
df_positive_treated$immune_sig <- tolower(df_positive_treated$NAME)
df_positive_vehicle <- read.csv(paste0(gsea_resultpath, "positive_200/", res_p_path, "/gsea_report_for_vehicle_1636630344414.xls"), 
                        sep = "\t", stringsAsFactors = F)
df_positive_vehicle = df_positive_vehicle[c(1, 5:9)]
df_positive_vehicle$immune_sig <- tolower(df_positive_vehicle$NAME)
df_positive = rbind(df_positive_treated, df_positive_vehicle)
df_positive = df_positive[df_positive$FDR.q.val < 0.05, ]
df_positive = df_positive[c("immune_sig", "NES", "FDR.q.val")]
names(df_positive) = c("immune_sig", "NES_p", "FDR_p")

res_n_path <- list.files(paste0(paste0(gsea_resultpath, "negative_200/")))
res_n_file <- list.files(paste0(gsea_resultpath, "negative_200/", res_n_path))
df_negative_treated <- read.csv(paste0(gsea_resultpath, "negative_200/", res_n_path, "/gsea_report_for_treated_1636630362380.xls"), 
                        sep = "\t", stringsAsFactors = F)
df_negative_treated = df_negative_treated[c(1, 5:9)]
df_negative_treated$immune_sig <- tolower(df_negative_treated$NAME)
df_negative_vehicle <- read.csv(paste0(gsea_resultpath, "negative_200/", res_n_path, "/gsea_report_for_vehicle_1636630362380.xls"), 
                        sep = "\t", stringsAsFactors = F)
df_negative_vehicle = df_negative_vehicle[c(1, 5:9)]
df_negative_vehicle$immune_sig <- tolower(df_negative_vehicle$NAME)
df_negative = rbind(df_negative_treated, df_negative_vehicle)
df_negative = df_negative[df_negative$FDR.q.val < 0.05, ]
df_negative = df_negative[c("immune_sig", "NES", "FDR.q.val")]
names(df_negative) = c("immune_sig", "NES_n", "FDR_n")

df_gsea = merge(df_positive, df_negative, by = "immune_sig", all = T)
df_gsea[is.na(df_gsea$NES_p),]$NES_p=0
df_gsea[is.na(df_gsea$NES_n),]$NES_n=0

df_gsea$NES_diff = df_gsea$NES_p - df_gsea$NES_n

df_gsea_immunesig = merge(df_gsea, immunesig, by = "immune_sig")
if(length(unique(df_gsea_immunesig$immune_sig)) < 0.1 * length(unique(df_gsea_immunesig$immune_sig))){
  print("analysis only with positive gene set")
  df_gsea_immunesig = merge(df_positive, immunesig, by = "immune_sig")
  df_gsea_immunesig$NES_diff = df_gsea_immunesig$NES_p - 0
}



drug_s = data.frame(drug19 = mat_s2['drug19',])
drug_r = data.frame(drug19 = mat_r2['drug19',])
drug_s$sign = 's'
drug_r$sign = 'r'
drug_s$x = as.numeric(1:nrow(drug_s))
drug_r$x = -1*c(1:nrow(drug_r))

zeroDot1 = data.frame(drug19=0, sign = 's', x = 0)
rownames(zeroDot1) = 0
zeroDot2 = data.frame(drug19=0, sign = 'r', x = 0)
rownames(zeroDot2) = 0

drug_dfplot = rbind(drug_s, 
                    zeroDot1,
                    zeroDot2,
                    drug_r)


drug_dfplot$immunesig = gsub("mNES_diff_s.", "", gsub("mNES_diff_r.", "", rownames(drug_dfplot)))

drug_dfplot[drug_dfplot$immunesig=="01",]$immunesig="0"
drug_dfplot$sign = factor(drug_dfplot$sign, levels = c('s','r'))
drug_dfplot$immunesig = factor(drug_dfplot$immunesig)

df_gsea_immunesig2 = df_gsea_immunesig[c("immune_sig", "NES_diff")]
names(df_gsea_immunesig2) = c("immunesig", "NES_diff")
df_gsea_immunesig2 = inner_join(df_gsea_immunesig2, drug_dfplot)
df_gsea_immunesig2[df_gsea_immunesig2$x <0,]$x = -1*seq(nrow(df_gsea_immunesig2[df_gsea_immunesig2$x < 0,]))
df_gsea_immunesig2[df_gsea_immunesig2$x >0,]$x = seq(nrow(df_gsea_immunesig2[df_gsea_immunesig2$x > 0,]))

dfplot2 = df_gsea_immunesig2[c(1,2,4,5)]
dfplot2$label = 'Combo_treated_mouse'
names(dfplot2) = c('immunesig','value','sign','x','label')

drug_dfplot = df_gsea_immunesig2[c(1,3,4,5)]
drug_dfplot$label = 'drug19'
names(drug_dfplot) = c('immunesig','value','sign','x','label')

dftmp = rbind(dfplot2, drug_dfplot)
dftmp = dftmp[order(dftmp$x), ]
dftmp$label = factor(dftmp$label, levels = c("drug19", "Combo_treated_mouse"))
# dftmp[is.na(dftmp$value),]$value = 0

p = ggplot(dftmp, aes(x = x, y = value, color=label, fill=label)) + 
  geom_line() + 
  geom_hline(aes(yintercept = 0)) +
  geom_vline(aes(xintercept = 0)) +
  geom_point(size = 3)+
  theme_bw() +
  theme(axis.text.x  =  element_text(angle=-90, vjust=0.5,hjust = 0, size=12),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_blank())+
  scale_x_continuous(breaks = unique(dftmp$x),
                     labels = unique(dftmp$immunesig))
p
ggsave("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/result_GSE114601/CellvstreatedMouse_COMBO_drug19.pdf",
p, height = 10, width=25 )

write.csv(df_gsea_immunesig2, "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/result_GSE114601/CellvstreatedMouse_COMBO_drug19.csv",
          row.names = F, quote = F)



drug_s = as.data.frame(t(mat_s['drug222',]))
drug_r = as.data.frame(t(mat_r['drug222',]))
drug_s$sign = 's'
drug_r$sign = 'r'
drug_s$x = as.numeric(1:nrow(drug_s))
drug_r$x = -1*c(1:nrow(drug_r))

zeroDot1 = data.frame(drug222=0, sign = 's', x = 0)
rownames(zeroDot1) = 0
zeroDot2 = data.frame(drug222=0, sign = 'r', x = 0)
rownames(zeroDot2) = 0

drug_dfplot = rbind(drug_s, 
               zeroDot1,
               zeroDot2,
               drug_r)


drug_dfplot$immunesig = gsub("mNES_diff_s.", "", gsub("mNES_diff_r.", "", rownames(drug_dfplot)))

drug_dfplot[drug_dfplot$immunesig=="01",]$immunesig="0"
drug_dfplot$sign = factor(drug_dfplot$sign, levels = c('s','r'))
drug_dfplot$immunesig = factor(drug_dfplot$immunesig)

df_gsea_immunesig2 = df_gsea_immunesig[c("immune_sig", "NES_diff")]
names(df_gsea_immunesig2) = c("immunesig", "NES_diff")
df_gsea_immunesig2 = inner_join(df_gsea_immunesig2, drug_dfplot)
df_gsea_immunesig2[df_gsea_immunesig2$x <0,]$x = -1*seq(nrow(df_gsea_immunesig2[df_gsea_immunesig2$x < 0,]))
df_gsea_immunesig2[df_gsea_immunesig2$x >0,]$x = seq(nrow(df_gsea_immunesig2[df_gsea_immunesig2$x > 0,]))

dfplot2 = df_gsea_immunesig2[c(1,2,4,5)]
dfplot2$label = 'Combo_treated_mouse'
names(dfplot2) = c('immunesig','value','sign','x','label')

drug_dfplot = df_gsea_immunesig2[c(1,3,4,5)]
drug_dfplot$label = 'drug222'
names(drug_dfplot) = c('immunesig','value','sign','x','label')

dftmp = rbind(dfplot2, drug_dfplot)
dftmp = dftmp[order(dftmp$x), ]
dftmp$label = factor(dftmp$label)
dftmp[is.na(dftmp$value),]$value = 0
dftmp$label = factor(dftmp$label, levels = c("drug222", "Combo_treated_mouse"))

p = ggplot(dftmp, aes(x = x, y = value, color=label, fill=label)) + 
  geom_line() + 
  geom_hline(aes(yintercept = 0)) +
  geom_vline(aes(xintercept = 0)) +
  geom_point(size = 3)+
  theme_bw() +
  theme(axis.text.x  =  element_text(angle=-90, vjust=0.5,hjust = 0, size=12),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_blank())+
  scale_x_continuous(breaks = unique(dftmp$x),
                     labels = unique(dftmp$immunesig))
p
ggsave("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/result_GSE114601/CellvstreatedMouse_COMBO_drug222.pdf",
p, height = 10, width=25 )

write.csv(df_gsea_immunesig2, "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/result_GSE114601/CellvstreatedMouse_COMBO_drug222.csv",
          row.names = F, quote = F)





# 2) Computing simple sample immune signatures ------------------------------------
library(GSVA)
vehicle_tpm = tpm[metadata[metadata$treatment == 'vehicle',]$sampleID]
load("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/03_drug_immuneSig_enrichment/data/gene_immune_sig_TUMERIC/for_GSVA/pcor/pearson_LUAD_Primary_positive_200.Rdata")
# enrichment
vehicle_p <- gsva(as.matrix(log2(vehicle_tpm + 1)),
                genelist_p,
                method = "ssgsea",
                kcdf = "Gaussian",
                min.sz > 1)
vehicle_p <- as.data.frame(vehicle_p)
names(vehicle_p) = paste0(names(vehicle_p),"_p")
vehicle_p$immunesig = rownames(vehicle_p)

load("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/03_drug_immuneSig_enrichment/data/gene_immune_sig_TUMERIC/for_GSVA/pcor/pearson_LUAD_Primary_negative_200.Rdata")
# enrichment
vehicle_n <- gsva(as.matrix(log2(vehicle_tpm + 1)),
                  genelist_n,
                  method = "ssgsea",
                  kcdf = "Gaussian",
                  min.sz > 1)
vehicle_n <- as.data.frame(vehicle_n)
names(vehicle_n) = paste0(names(vehicle_n),"_n")
vehicle_n$immunesig = rownames(vehicle_n)


# test using GVA
# library(GSVA)
# load("GSE114601_JQ1/pearson_LUAD_Primary_positive_200.Rdata")
# # enrichment
# vehicle_p <- gsva(as.matrix(log2(tpm + 1)),
#               genelist_p,
#               method = "gsva",
#               kcdf = "Gaussian",
#               min.sz > 1)
# vehicle_p <- as.data.frame(vehicle_p)
# vehicle_p = vehicle_p[metadata[metadata$treatment == 'vehicle',]$sampleID]
# names(vehicle_p) = paste0(names(vehicle_p),"_p")
# vehicle_p$immunesig = rownames(vehicle_p)
# 
# load("GSE114601_JQ1/pearson_LUAD_Primary_negative_200.Rdata")
# # enrichment
# vehicle_n <- gsva(as.matrix(log2(tpm + 1)),
#               genelist_n,
#               method = "gsva",
#               kcdf = "Gaussian",
#               min.sz > 1)
# vehicle_n <- as.data.frame(vehicle_n)
# vehicle_n = vehicle_n[metadata[metadata$treatment == 'vehicle',]$sampleID]
# names(vehicle_n) = paste0(names(vehicle_n),"_n")
# vehicle_n$immunesig = rownames(vehicle_n)

vehicle_SigScore <- inner_join(vehicle_p, vehicle_n)
vehicle_SigScore$s2596 <- vehicle_SigScore$s2596_p - vehicle_SigScore$s2596_n
vehicle_SigScore$s2617 <- vehicle_SigScore$s2617_p - vehicle_SigScore$s2617_n
vehicle_SigScore <- vehicle_SigScore[c('immunesig', c('s2596', 's2617'))]
vehicle_SigScore$immunesig = tolower(vehicle_SigScore$immunesig)

# Single sample drug screening -------------------------------------------------
# drag file to local first
# load("GSE114601_JQ1/LUAD_Primary_drugITS.Rdata")
# mat_s2
# mat_r2
colnames(mat_s2) = gsub("mNES_diff_s.", "", colnames(mat_s2))
mat_s = as.data.frame(t(mat_s2))
vehicle_SigScore_s = vehicle_SigScore[is.element(vehicle_SigScore$immunesig, rownames(mat_s)),]



mat_s1 = mat_s[vehicle_SigScore_s$immunesig, ]
rownames(vehicle_SigScore_s) = vehicle_SigScore_s$immunesig
vehicle_SigScore_s = vehicle_SigScore_s[-1]
vehicle_pred_s <- lapply(as.list(colnames(vehicle_SigScore_s)), function(s){
  tmp = vehicle_SigScore_s[,s]
  # colSums(tmp * mat_s1)/colSums(mat_s1)
  colSums(tmp + mat_s1)/nrow(mat_s1)
  # tmp_z = (tmp-mean(tmp))/sd(tmp)
  # mat_s1_z = apply(mat_s1, 2, function(x)x-mean(x)/sd(x))
  # colSums(tmp_z + mat_s1_z)/nrow(mat_s1_z)
})
names(vehicle_pred_s) = colnames(vehicle_SigScore_s)
vehicle_pred_s = t(do.call(rbind, vehicle_pred_s))
vehicle_pred_s = as.data.frame(vehicle_pred_s)
names(vehicle_pred_s) = paste0(names(vehicle_pred_s), "_s") 
vehicle_pred_s$drug = rownames(vehicle_pred_s)

# vehicle_pred_s[order(vehicle_pred_s$s2596, decreasing = T), ]

colnames(mat_r2) = gsub("mNES_diff_r.", "", colnames(mat_r2))
mat_r = as.data.frame(t(mat_r2))
vehicle_SigScore_r = vehicle_SigScore[is.element(vehicle_SigScore$immunesig, rownames(mat_r)),]
mat_r1 = mat_r[vehicle_SigScore_r$immunesig, ]
rownames(vehicle_SigScore_r) = vehicle_SigScore_r$immunesig
vehicle_SigScore_r = vehicle_SigScore_r[-1]
vehicle_pred_r <- lapply(as.list(colnames(vehicle_SigScore_r)), function(s){
  tmp = vehicle_SigScore_r[,s]
  # colSums(tmp * mat_r1)/colSums(mat_r1)
  colSums(tmp + mat_r1)/nrow(mat_r1)
  # tmp_z = (tmp-mean(tmp))/sd(tmp)
  # mat_r1_z = apply(mat_r1, 2, function(x)x-mean(x)/sd(x))
  # colSums(tmp_z + mat_r1_z)/nrow(mat_r1_z)
  
})
names(vehicle_pred_r) = colnames(vehicle_SigScore_r)
vehicle_pred_r = t(do.call(rbind, vehicle_pred_r))

vehicle_pred_s[rownames(vehicle_pred_s) == "drug19",]
vehicle_pred_r[rownames(vehicle_pred_r) == "drug19",]
vehicle_pred_s[rownames(vehicle_pred_s) == "drug222",]
vehicle_pred_r[rownames(vehicle_pred_r) == "drug222",]


vehicle_pred_r = as.data.frame(vehicle_pred_r)
# vehicle_pred_r[order(vehicle_pred_r$s2596, decreasing = T), ]
names(vehicle_pred_r) = paste0(names(vehicle_pred_r), "_r") 
vehicle_pred_r$drug = rownames(vehicle_pred_r)

# plot
vehicle_SigScore_s$sign = 's'
vehicle_SigScore_r$sign = 'r'
vehicle_SigScore_s$x = as.numeric(1:nrow(vehicle_SigScore_s))
vehicle_SigScore_r$x = -1*c(1:nrow(vehicle_SigScore_r))

zeroDot1 = data.frame(s2596=0, s2617=0, sign = 's', x = 0)
rownames(zeroDot1) = 0
zeroDot2 = data.frame(s2596=0, s2617=0, sign = 'r', x = 0)
rownames(zeroDot2) = 0

dfplot = rbind(vehicle_SigScore_s, 
               zeroDot1,
               zeroDot2,
               vehicle_SigScore_r)


dfplot$immunesig = rownames(dfplot)
dfplot[dfplot$immunesig=="01",]$immunesig="0"
dfplot$sign = factor(dfplot$sign, levels = c('s','r'))
dfplot$immunesig = factor(dfplot$immunesig)


patient_plot = ggplot(dfplot, aes(x = x, y = s2596, color=sign, fill=sign)) + 
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
drug19_s = mat_s1['drug19']
drug19_r = mat_r1['drug19']
drug19_s$sign = 's'
drug19_r$sign = 'r'
drug19_s$x = as.numeric(1:nrow(drug19_s))
drug19_r$x = -1*c(1:nrow(drug19_r))

zeroDot1 = data.frame(drug19=0, sign = 's', x = 0)
rownames(zeroDot1) = 0
zeroDot2 = data.frame(drug19=0, sign = 'r', x = 0)
rownames(zeroDot2) = 0

drug_dfplot = rbind(drug19_s, 
               zeroDot1,
               zeroDot2,
               drug19_r)


drug_dfplot$immunesig = rownames(drug_dfplot)
drug_dfplot[drug_dfplot$immunesig=="01",]$immunesig="0"
drug_dfplot$sign = factor(drug_dfplot$sign, levels = c('s','r'))
drug_dfplot$immunesig = factor(drug_dfplot$immunesig)


drug_plot = ggplot(drug_dfplot, aes(x = x, y = drug19, color=sign, fill=sign)) + 
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

dfplot2 = dfplot[-1]
dfplot2$label = 'patient'
names(dfplot2) = c('value','sign','x','immunesig','label')
drug_dfplot$label = 'drug19'
names(drug_dfplot) = c('value','sign','x','immunesig','label')

dftmp = rbind(dfplot2, drug_dfplot)
dftmp$label = factor(dftmp$label)

p1 = ggplot(dftmp, aes(x = x, y = value, color=label, fill=label)) + 
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
p1

ggsave("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/data/immune_combo_dataset/GSE114601_JQ1/s2596_drug19.pdf",
p1, height = 10, width=25 )


dfplot2 = dfplot[-2]
dfplot2$label = 'patient'
names(dfplot2) = c('value','sign','x','immunesig','label')
drug_dfplot$label = 'drug19'
names(drug_dfplot) = c('value','sign','x','immunesig','label')

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
ggsave("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/data/immune_combo_dataset/GSE114601_JQ1/s2617_drug19.pdf",
p2, height = 10, width=25 )




drug222_s = mat_s1['drug222']
drug222_r = mat_r1['drug222']
drug222_s$sign = 's'
drug222_r$sign = 'r'
drug222_s$x = as.numeric(1:nrow(drug222_s))
drug222_r$x = -1*c(1:nrow(drug222_r))

zeroDot1 = data.frame(drug222=0, sign = 's', x = 0)
rownames(zeroDot1) = 0
zeroDot2 = data.frame(drug222=0, sign = 'r', x = 0)
rownames(zeroDot2) = 0

drug_dfplot = rbind(drug222_s, 
               zeroDot1,
               zeroDot2,
               drug222_r)


drug_dfplot$immunesig = rownames(drug_dfplot)
drug_dfplot[drug_dfplot$immunesig=="01",]$immunesig="0"
drug_dfplot$sign = factor(drug_dfplot$sign, levels = c('s','r'))
drug_dfplot$immunesig = factor(drug_dfplot$immunesig)



dfplot2 = dfplot[-1]
dfplot2$label = 'patient'
names(dfplot2) = c('value','sign','x','immunesig','label')
drug_dfplot$label = 'drug222'
names(drug_dfplot) = c('value','sign','x','immunesig','label')

dftmp = rbind(dfplot2, drug_dfplot)
dftmp$label = factor(dftmp$label)

p1 = ggplot(dftmp, aes(x = x, y = value, color=label, fill=label)) + 
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
p1

ggsave("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/data/immune_combo_dataset/GSE114601_JQ1/s2596_drug222.pdf",
p1, height = 10, width=25 )


dfplot2 = dfplot[-2]
dfplot2$label = 'patient'
names(dfplot2) = c('value','sign','x','immunesig','label')
drug_dfplot$label = 'drug222'
names(drug_dfplot) = c('value','sign','x','immunesig','label')

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

ggsave("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/data/immune_combo_dataset/GSE114601_JQ1/s2617_drug222.pdf",
p2, height = 10, width=25 )



vehicle_pred = inner_join(vehicle_pred_s, vehicle_pred_r)
head(vehicle_pred)
vehicle_pred$s2596 = vehicle_pred$s2596_s - vehicle_pred$s2596_r 
vehicle_pred$s2617 = vehicle_pred$s2617_s - vehicle_pred$s2617_r 
vehicle_pred$score = (vehicle_pred$s2596 + vehicle_pred$s2617)/2
# vehicle_pred = vehicle_pred[c('drug', 's2596', 's2617')]
head(vehicle_pred[order(vehicle_pred$score, decreasing = T), ],10)
head(vehicle_pred[order(vehicle_pred$s2596, decreasing = T), ],10)
head(vehicle_pred[order(vehicle_pred$s2617, decreasing = T), ],10)
save(vehicle_pred, vehicle_pred_s, vehicle_pred_r, file = "GSE114601_JQ1/prediction_result.Rdata")
# 

vehicle_pred = vehicle_pred[order(vehicle_pred$score, decreasing = T),]
vehicle_pred$rank = seq(nrow(vehicle_pred))
vehicle_pred[vehicle_pred$drug=='drug222',]
vehicle_pred[vehicle_pred$drug=='drug19',]
# 
# vehicle_pred$s2596 = (vehicle_pred$s2596_s - mean(vehicle_pred$s2596_s))/sd(vehicle_pred$s2596_s) - (vehicle_pred$s2596_r - mean(vehicle_pred$s2596_r))/sd(vehicle_pred$s2596_r)
# vehicle_pred$s2617 = (vehicle_pred$s2617_s - mean(vehicle_pred$s2617_s))/sd(vehicle_pred$s2617_s) - (vehicle_pred$s2617_r - mean(vehicle_pred$s2617_r))/sd(vehicle_pred$s2617_r)
# # vehicle_pred = vehicle_pred[c('drug', 's2596', 's2617')]
# vehicle_pred[vehicle_pred$drug=='drug222',]
# vehicle_pred[vehicle_pred$drug=='drug19',]
# head(vehicle_pred[order(vehicle_pred$s2596, decreasing = T), ],10)
# head(vehicle_pred[order(vehicle_pred$s2617, decreasing = T), ],10)
# 
# 

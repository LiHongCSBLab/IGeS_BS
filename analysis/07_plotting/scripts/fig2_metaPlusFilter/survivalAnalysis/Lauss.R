# rm(list=ls())
# ------------------------------------------------------------------------------
# source("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/06_results_summaryandplotting/scripts/immunesig_selection_function.R")
workpath <- "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/"
setwd(workpath)
options(stringsAsFactors = F)

# Checking and loading R packages -----------------------------------------
library(dplyr)
library(ggplot2)
library(pheatmap)
library(ROCR)
library(reshape2)
library(survival)
library(survminer) 

# Loading relied functions -----------------------------------------------------
# TIGS, AUC and Wilcox test
# source("05_immuneSig_ICBresponse/scripts/functions/05_immunotherapy_immuneSig_analysis_function.R")
# source("05_immuneSig_ICBresponse/scripts/functions/predictor_ITS.R")
# source("06_results_summaryandplotting/scripts/immunesig_selection_function.R")



result_dir <- "05_immuneSig_ICBresponse/results_r0.4_ITSssgsea_CombatRef/"
resultpath = paste0(result_dir, "finalmodel_17datasets_ssgseanorm/metaAnalysis/")

load(paste0(result_dir, "ITS_afterCombat17datasets_ssgsea.Rdata"))
datasets = unique(ITSres_selected$dataset)
vlSetName <- c("08_Lauss_dataset_outcome", "GSE115821_pre_outcome", "GSE96619_outcome","Liu_dataset_CTLA4Naive_outcome", "GSE176307")
dcSetName <- datasets[!is.element(datasets, vlSetName)]


vlSetName <- c("08_Lauss_dataset_outcome", "GSE115821_pre_outcome", "GSE96619_outcome","Liu_dataset_CTLA4Naive_outcome", "GSE176307")


# ITS-based models -------------------------------------------------------------
m = 'elastic_net'
load(paste0(resultpath, "/en/elasticnet_model_pred.Rdata"))
pred = pred_en$data
names(pred)[grep('truth',names(pred))] = 'tag'

res_pred <- pred
res_pred$SampleID <- row.names(res_pred)
res_pred <- merge(ITSres_selected[c('SampleID', 'label')], res_pred)


# create savepath -----------------------------------------------------------------
dir.create("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/05_immuneSig_ICBresponse/results_r0.4_ITSssgsea_CombatRef/survivalAnalysis/")
dir.create("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/07_plotting_v2/02_ITSprediction_metaPlusFilter/survivalAnalysis/")
plot_savepath = "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/07_plotting_v2/02_ITSprediction_metaPlusFilter/survivalAnalysis/"

# compute ITS and significant test between R and NR ----------------------------
# generating immune signatures -------------------------------------------------

processedDataPath = "/picb/bigdata/project/FengFYM/Immunotherapy_prediction/data/RNAseq/"
dir.create(processedDataPath)

# Gide_dataset -----------------------------------------------------------------
# survival: Days
dset = 'Lauss'
dir.create(paste0(plot_savepath, "/", dset, "/"))

dataset = "GSE100797.Rdata"
load(paste0(processedDataPath, "/Lauss_dataset/processedData/", dataset))
cancer = "SKCM"
sample_type = "Metastatic"




names(survival) = c("patientID", "os", "os_Status", "pfs", "pfs_Status")

# cox regression
dset = 'Lauss'
df_pred = res_pred[grep(dset, res_pred$SampleID),]
df_pred$patientID = gsub("08_Lauss_dataset_outcome", "", df_pred$SampleID)

length(intersect(df_pred$patientID, survival$patientID))
y = merge(df_pred, survival, by = 'patientID')


# OS ---------------------------------------------------------------------------
# real response and survival
fit_os <- survfit(Surv(os, os_Status) ~ label, data = y)
# summary(fit_os)
p_os = ggsurvplot(fit_os, # 创建的拟合对象
           data = y,  # 指定变量数据来源
           conf.int = TRUE, # 显示置信区间
           pval = TRUE, # 添加P值
           # surv.median.line = "hv",  # 添加中位生存时间线
           risk.table = TRUE,
            xlab = "Days")  # 设置x轴刻度间距
# print(p_os)
pdf(paste0(plot_savepath, "/",dset,"/os_dataset_TrueResponse.pdf"))
print(p_os, newpage = FALSE)
dev.off()


# take median as survival cutoff
y$response_pred = 'NR'
y[y$prob.1 > median(y$prob.1), ]$response_pred = 'R'
fit_os <- survfit(Surv(os, os_Status) ~ response_pred, data = y)
# summary(fit_os)
p_os = ggsurvplot(fit_os, # 创建的拟合对象
           data = y,  # 指定变量数据来源
           conf.int = TRUE, # 显示置信区间
           pval = TRUE, # 添加P值
           # surv.median.line = "hv",  # 添加中位生存时间线
           risk.table = TRUE,
            xlab = "Days")  # 设置x轴刻度间距
# print(p_os)
# dev.off()
pdf(paste0(plot_savepath, "/",dset,"/os_dataset_pred_meidan.pdf"))
print(p_os, newpage = FALSE)
dev.off()


# take mean as survival cutoff
y$response_pred = 'NR'
y[y$prob.1 > mean(y$prob.1), ]$response_pred = 'R'
fit_os <- survfit(Surv(os, os_Status) ~ response_pred, data = y)
# summary(fit_os)
p_os = ggsurvplot(fit_os, # 创建的拟合对象
           data = y,  # 指定变量数据来源
           conf.int = TRUE, # 显示置信区间
           pval = TRUE, # 添加P值
           # surv.median.line = "hv",  # 添加中位生存时间线
           risk.table = TRUE,
            xlab = "Days")  # 设置x轴刻度间距
# print(p_os)
pdf(paste0(plot_savepath, "/",dset,"/os_dataset_pred_mean.pdf"))
print(p_os, newpage = FALSE)
dev.off()


# take best surv cutoff as survival cutoff
surv_cutoff <- surv_cutpoint(
   y,
   time = "os",
   event = "os_Status",
   variables = "prob.1"
)
# summary(surv_cutoff)
# plot(surv_cutoff, "prob.1", palette = "npg")
# dev.off()
y$response_pred = 'NR'
y[y$prob.1 > summary(surv_cutoff)[1,1], ]$response_pred = 'R'
fit_os <- survfit(Surv(os, os_Status) ~ response_pred, data = y)
# summary(fit_os)
p_os = ggsurvplot(fit_os, # 创建的拟合对象
           data = y,  # 指定变量数据来源
           conf.int = TRUE, # 显示置信区间
           pval = TRUE, # 添加P值
           # surv.median.line = "hv",  # 添加中位生存时间线
           risk.table = TRUE,
            xlab = "Days")  # 设置x轴刻度间距
# print(p_os)
pdf(paste0(plot_savepath, "/",dset,"/os_dataset_pred_SurvBestCutoff.pdf"))
print(p_os, newpage = FALSE)
dev.off()

# take best response classification cutoff as survival cutoff
library(pROC)
rocobj <- roc(y$label, y$prob.1)
roc_result <- coords(rocobj, "best", best.method= "youden")
y$response_pred = 'NR'
y[y$prob.1 > roc_result[1,1], ]$response_pred = 'R'
fit_os <- survfit(Surv(os, os_Status) ~ response_pred, data = y)
# summary(fit_os)
p_os = ggsurvplot(fit_os, # 创建的拟合对象
           data = y,  # 指定变量数据来源
           conf.int = TRUE, # 显示置信区间
           pval = TRUE, # 添加P值
           # surv.median.line = "hv",  # 添加中位生存时间线
           risk.table = TRUE,
            xlab = "Days")  # 设置x轴刻度间距
# print(p_os)
pdf(paste0(plot_savepath, "/",dset,"/os_dataset_pred_ROCBestCutoff_youden.pdf"))
print(p_os, newpage = FALSE)
dev.off()



rocobj <- roc(y$label, y$prob.1)
roc_result <- coords(rocobj, "best", best.method= "closest.topleft")
roc_result
y$response_pred = 'NR'
y[y$prob.1 > roc_result[1,1], ]$response_pred = 'R'
fit_os <- survfit(Surv(os, os_Status) ~ response_pred, data = y)
# summary(fit_os)
p_os = ggsurvplot(fit_os, # 创建的拟合对象
           data = y,  # 指定变量数据来源
           conf.int = TRUE, # 显示置信区间
           pval = TRUE, # 添加P值
           # surv.median.line = "hv",  # 添加中位生存时间线
           risk.table = TRUE,
            xlab = "Days")  # 设置x轴刻度间距
# print(p_os)
pdf(paste0(plot_savepath, "/",dset,"/os_dataset_pred_ROCBestCutoff_closesttopleft.pdf"))
print(p_os, newpage = FALSE)
dev.off()


# take best response classification cutoff as survival cutoff
library(pROC)
rocobj <- roc(res_pred$label, res_pred$prob.1)
roc_result <- coords(rocobj, "best", best.method= "youden")
roc_result
y$response_pred = 'NR'
y[y$prob.1 > roc_result[1,1], ]$response_pred = 'R'
fit_os <- survfit(Surv(os, os_Status) ~ response_pred, data = y)
# summary(fit_os)
p_os = ggsurvplot(fit_os, # 创建的拟合对象
           data = y,  # 指定变量数据来源
           conf.int = TRUE, # 显示置信区间
           pval = TRUE, # 添加P值
           # surv.median.line = "hv",  # 添加中位生存时间线
           risk.table = TRUE,
            xlab = "Days")  # 设置x轴刻度间距
# print(p_os)
# dev.off()
pdf(paste0(plot_savepath, "/",dset,"/os_global_pred_ROCBestCutoff_youden.pdf"))
print(p_os, newpage = FALSE)
dev.off()

rocobj <- roc(res_pred$label, res_pred$prob.1)
roc_result <- coords(rocobj, "best", best.method= "closest.topleft")
roc_result
y$response_pred = 'NR'
y[y$prob.1 > roc_result[1,1], ]$response_pred = 'R'
fit_os <- survfit(Surv(os, os_Status) ~ response_pred, data = y)
# summary(fit_os)
p_os = ggsurvplot(fit_os, # 创建的拟合对象
           data = y,  # 指定变量数据来源
           conf.int = TRUE, # 显示置信区间
           pval = TRUE, # 添加P值
           # surv.median.line = "hv",  # 添加中位生存时间线
           risk.table = TRUE,
            xlab = "Days")  # 设置x轴刻度间距
# print(p_os)
# dev.off()
pdf(paste0(plot_savepath, "/",dset,"/os_global_pred_ROCBestCutoff_closesttopleft.pdf"))
print(p_os, newpage = FALSE)
dev.off()


# pfs ---------------------------------------------------------------------------
# real response and survival
fit_pfs <- survfit(Surv(pfs, pfs_Status) ~ label, data = y)
# summary(fit_pfs)
p_pfs = ggsurvplot(fit_pfs, # 创建的拟合对象
           data = y,  # 指定变量数据来源
           conf.int = TRUE, # 显示置信区间
           pval = TRUE, # 添加P值
           # surv.median.line = "hv",  # 添加中位生存时间线
           risk.table = TRUE,
            xlab = "Days")  # 设置x轴刻度间距
# print(p_pfs)
pdf(paste0(plot_savepath, "/",dset,"/pfs_dataset_TrueResponse.pdf"))
print(p_pfs, newpage = FALSE)
dev.off()


# take median as survival cutoff
y$response_pred = 'NR'
y[y$prob.1 > median(y$prob.1), ]$response_pred = 'R'
fit_pfs <- survfit(Surv(pfs, pfs_Status) ~ response_pred, data = y)
# summary(fit_pfs)
p_pfs = ggsurvplot(fit_pfs, # 创建的拟合对象
           data = y,  # 指定变量数据来源
           conf.int = TRUE, # 显示置信区间
           pval = TRUE, # 添加P值
           # surv.median.line = "hv",  # 添加中位生存时间线
           risk.table = TRUE,
            xlab = "Days")  # 设置x轴刻度间距
# print(p_pfs)
# dev.off()
pdf(paste0(plot_savepath, "/",dset,"/pfs_dataset_pred_meidan.pdf"))
print(p_pfs, newpage = FALSE)
dev.off()


# take mean as survival cutoff
y$response_pred = 'NR'
y[y$prob.1 > mean(y$prob.1), ]$response_pred = 'R'
fit_pfs <- survfit(Surv(pfs, pfs_Status) ~ response_pred, data = y)
# summary(fit_pfs)
p_pfs = ggsurvplot(fit_pfs, # 创建的拟合对象
           data = y,  # 指定变量数据来源
           conf.int = TRUE, # 显示置信区间
           pval = TRUE, # 添加P值
           # surv.median.line = "hv",  # 添加中位生存时间线
           risk.table = TRUE,
            xlab = "Days")  # 设置x轴刻度间距
# print(p_pfs)
pdf(paste0(plot_savepath, "/",dset,"/pfs_dataset_pred_mean.pdf"))
print(p_pfs, newpage = FALSE)
dev.off()


# take best surv cutoff as survival cutoff
surv_cutoff <- surv_cutpoint(
   y,
   time = "pfs",
   event = "pfs_Status",
   variables = "prob.1"
)
# summary(surv_cutoff)
# plot(surv_cutoff, "prob.1", palette = "npg")
# dev.off()
y$response_pred = 'NR'
y[y$prob.1 > summary(surv_cutoff)[1,1], ]$response_pred = 'R'
fit_pfs <- survfit(Surv(pfs, pfs_Status) ~ response_pred, data = y)
# summary(fit_pfs)
p_pfs = ggsurvplot(fit_pfs, # 创建的拟合对象
           data = y,  # 指定变量数据来源
           conf.int = TRUE, # 显示置信区间
           pval = TRUE, # 添加P值
           # surv.median.line = "hv",  # 添加中位生存时间线
           risk.table = TRUE,
            xlab = "Days")  # 设置x轴刻度间距
# print(p_pfs)
pdf(paste0(plot_savepath, "/",dset,"/pfs_dataset_pred_SurvBestCutoff.pdf"))
print(p_pfs, newpage = FALSE)
dev.off()

# take best response classification cutoff as survival cutoff
library(pROC)
rocobj <- roc(y$label, y$prob.1)
roc_result <- coords(rocobj, "best", best.method= "youden")
y$response_pred = 'NR'
y[y$prob.1 > roc_result[1,1], ]$response_pred = 'R'
fit_pfs <- survfit(Surv(pfs, pfs_Status) ~ response_pred, data = y)
# summary(fit_pfs)
p_pfs = ggsurvplot(fit_pfs, # 创建的拟合对象
           data = y,  # 指定变量数据来源
           conf.int = TRUE, # 显示置信区间
           pval = TRUE, # 添加P值
           # surv.median.line = "hv",  # 添加中位生存时间线
           risk.table = TRUE,
            xlab = "Days")  # 设置x轴刻度间距
# print(p_pfs)
pdf(paste0(plot_savepath, "/",dset,"/pfs_dataset_pred_ROCBestCutoff_youden.pdf"))
print(p_pfs, newpage = FALSE)
dev.off()



rocobj <- roc(y$label, y$prob.1)
roc_result <- coords(rocobj, "best", best.method= "closest.topleft")
roc_result
y$response_pred = 'NR'
y[y$prob.1 > roc_result[1,1], ]$response_pred = 'R'
fit_pfs <- survfit(Surv(pfs, pfs_Status) ~ response_pred, data = y)
# summary(fit_pfs)
p_pfs = ggsurvplot(fit_pfs, # 创建的拟合对象
           data = y,  # 指定变量数据来源
           conf.int = TRUE, # 显示置信区间
           pval = TRUE, # 添加P值
           # surv.median.line = "hv",  # 添加中位生存时间线
           risk.table = TRUE,
            xlab = "Days")  # 设置x轴刻度间距
# print(p_pfs)
pdf(paste0(plot_savepath, "/",dset,"/pfs_dataset_pred_ROCBestCutoff_closesttopleft.pdf"))
print(p_pfs, newpage = FALSE)
dev.off()


# take best response classification cutoff as survival cutoff
library(pROC)
rocobj <- roc(res_pred$label, res_pred$prob.1)
roc_result <- coords(rocobj, "best", best.method= "youden")
roc_result
y$response_pred = 'NR'
y[y$prob.1 > roc_result[1,1], ]$response_pred = 'R'
fit_pfs <- survfit(Surv(pfs, pfs_Status) ~ response_pred, data = y)
# summary(fit_pfs)
p_pfs = ggsurvplot(fit_pfs, # 创建的拟合对象
           data = y,  # 指定变量数据来源
           conf.int = TRUE, # 显示置信区间
           pval = TRUE, # 添加P值
           # surv.median.line = "hv",  # 添加中位生存时间线
           risk.table = TRUE,
            xlab = "Days")  # 设置x轴刻度间距
# print(p_pfs)
# dev.off()
pdf(paste0(plot_savepath, "/",dset,"/pfs_global_pred_ROCBestCutoff_youden.pdf"))
print(p_pfs, newpage = FALSE)
dev.off()

rocobj <- roc(res_pred$label, res_pred$prob.1)
roc_result <- coords(rocobj, "best", best.method= "closest.topleft")
roc_result
y$response_pred = 'NR'
y[y$prob.1 > roc_result[1,1], ]$response_pred = 'R'
fit_pfs <- survfit(Surv(pfs, pfs_Status) ~ response_pred, data = y)
# summary(fit_pfs)
p_pfs = ggsurvplot(fit_pfs, # 创建的拟合对象
           data = y,  # 指定变量数据来源
           conf.int = TRUE, # 显示置信区间
           pval = TRUE, # 添加P值
           # surv.median.line = "hv",  # 添加中位生存时间线
           risk.table = TRUE,
            xlab = "Days")  # 设置x轴刻度间距
# print(p_pfs)
# dev.off()
pdf(paste0(plot_savepath, "/",dset,"/pfs_global_pred_ROCBestCutoff_closesttopleft.pdf"))
print(p_pfs, newpage = FALSE)
dev.off()

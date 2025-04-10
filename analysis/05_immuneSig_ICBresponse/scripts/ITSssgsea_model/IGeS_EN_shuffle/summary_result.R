# rm(list=ls())
workpath <- "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/"
setwd(workpath)

options(stringsAsFactors = F)

# Checking and loading R packages -----------------------------------------

if (!require(readxl)) {  install.packages('readxl')}
if (!require(data.table)) {  install.packages('data.table')}
if (!require(dplyr)) {  install.packages('dplyr')}
library(ROCR)
library(reshape2)
library(getopt)
library(mlr)
library(mlrMBO)
library(openxlsx)

# load inital EN -----------------------------------------


immunesig <- read.csv("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/05_immuneSig_ICBresponse/results/meta_results/res_metaRes.csv")
immunesig <- immunesig[immunesig$Meta_Pval < 0.05, ]
names(immunesig) <- c("immune_sig", names(immunesig)[-1])
immunesig$flag <- "sensitive"
immunesig[immunesig$OR < 1,]$flag <- "resistant"
immunesig$setindex = "01_sensitive"
immunesig[immunesig$flag == "resistant", ]$setindex = "02_resistant"

load("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/05_immuneSig_ICBresponse/results_r0.4_ITSssgsea_CombatRef/finalmodel_17datasets_ssgseanorm/metaAnalysis_filter/en/elasticnet_model_pred.Rdata")

en_model_para = m_en$learner.model$glmnet.fit
m_en_imp <- data.frame(en_model_para$beta[,which(en_model_para$lambda == m_en$learner.model$lambda.min)])

names(m_en_imp) <- "weight"
m_en_imp$immune_sig = row.names(m_en_imp)
immunesig = merge(immunesig, m_en_imp, by = "immune_sig")
immunesig =immunesig[immunesig$weight != 0 ,]


# load rand EN -----------------------------------------
res_auc_rand = as.list(1:10)
res_auc_all_rand = as.list(1:10)
immunesig_rand = as.list(1:10)
for(rand_index in c(1:10)){

    result_dir <- paste0("07_plotting_v2/IGeS_EN_shuffle/rand", rand_index, "/")

    resultpath = paste0(result_dir, "finalmodel_17datasets_ssgseanorm/")
    featureFilterStrategy = "metaAnalysis_filter" 
    resultpath = paste0(resultpath, "/", featureFilterStrategy, "/en")
    load( paste0(resultpath, "/elasticnet_model_pred.Rdata"))
    pred2 <- ROCR::prediction(pred_en$data$prob.1 , pred_en$data$truth)
    aucPerf <- ROCR::performance( pred2, 'auc' )
    AUCValue<-aucPerf@y.values[[1]]
    res_auc_rand[[rand_index]] <- data.frame(rand = rand_index, auc_all = AUCValue)

    df1 <- read.csv(paste0(resultpath, "/res_auc_all.csv"), row.names=1)
    df1$rand = rand_index
    res_auc_all_rand[[rand_index]] <- df1
    df2 <- read.csv(paste0(resultpath, "/immunesig_", rand_index, ".csv"), row.names=1)
    df2$rand = rand_index
    immunesig_rand[[rand_index]] <- df2

}
res_auc_rand_merge = do.call(rbind, res_auc_rand)
res_auc_all_rand_merge <- do.call(rbind, res_auc_all_rand)
summary(res_auc_all_rand_merge$auc)
immunesig_rand_merge <- do.call(rbind, immunesig_rand)
immunesig_count = as.data.frame(table(immunesig_rand_merge$immune_sig))
aggregate(res_auc_all_rand_merge$auc, list(res_auc_all_rand_merge$rand), mean)
length(intersect(immunesig_count[immunesig_count$Freq>3,]$Var1, immunesig$immune_sig))/length(immunesig$immune_sig)
length(intersect(immunesig_count[immunesig_count$Freq>4,]$Var1, immunesig$immune_sig))/length(immunesig$immune_sig)
length(intersect(immunesig_count[immunesig_count$Freq>5,]$Var1, immunesig$immune_sig))/length(immunesig$immune_sig)
length(intersect(immunesig_count[immunesig_count$Freq>6,]$Var1, immunesig$immune_sig))/length(immunesig$immune_sig)
immunesig_count[immunesig_count$Freq>8, ]$Var1
res_auc_rand_merge
lapply(immunesig_rand, function(x){
   length(intersect( x$immune_sig, immunesig$immune_sig))
})

write.xlsx(list(auc_perRand = res_auc_rand_merge,
                auc_per_testst = res_auc_all_rand_merge,
                immunesig_rand = immunesig_rand_merge,
                immunesig_count=immunesig_count),
"07_plotting_v2/IGeS_EN_shuffle/res_auc_all_rand_merge_testset.xlsx")


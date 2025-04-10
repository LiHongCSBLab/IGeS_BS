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

# Loading relied functions -----------------------------------------------------

# mlr --------------------------------------------------------------------------
library(mlr)
library(mlrMBO)

# mlr feature selection ----------------------------------------------------------------
mlr_feature_selection <- function(data,
                                  filtermethod = c("all", "auc","variance", "univariate.model.score",
                                                  "FSelector_gain.ratio","FSelector_chi.squared",
                                                  "FSelector_information.gain","FSelector_relief",
                                                  "FSelector_symmetrical.uncertainty","kruskal.test",
                                                  "party_cforest.importance",
                                                  "praznik_CMIM","praznik_DISR","praznik_JMI","praznik_JMIM",
                                                  "praznik_MIM","praznik_MRMR","praznik_NJMIM",
                                                  "randomForest_importance", 
                                                  "randomForestSRC_importance", "randomForestSRC_var.select",
                                                  "ranger_impurity", "ranger_permutation")
                                  ){

  ICBtask <- makeClassifTask(id = "ICBresponse", 
                         data = data, 
                         target = "tag", 
                         positive = 1)
  if(filtermethod == "all"){
     filtermethod = c("auc","variance","kruskal.test", # "univariate.model.score",
                    "FSelector_gain.ratio","FSelector_chi.squared",
                    "FSelector_information.gain","FSelector_relief",
                    "FSelector_symmetrical.uncertainty",
                    "praznik_CMIM","praznik_DISR","praznik_JMI","praznik_JMIM",
                    "praznik_MIM","praznik_MRMR","praznik_NJMIM",
                    # "party_cforest.importance",
                    "randomForest_importance", 
                    "randomForestSRC_importance", 
                    "randomForestSRC_var.select",
                    "ranger_impurity", "ranger_permutation") 
     # for(f in filtermethod){
     #      print(f)
     fv = generateFilterValuesData(ICBtask, 
                                method = filtermethod)
     # }

  }else{
     fv = generateFilterValuesData(ICBtask, 
                                method = filtermethod)
  }
  
  fv2 = as.data.frame(fv$data)
  return(fv2)
}

# ------------------------------------------------------------------------------
# create directory -------------------------------------------------------------
# ------------------------------------------------------------------------------
resultpath = "05_immuneSig_ICBresponse/results/feature_importance_mlr/"
dir.create(resultpath)

# ------------------------------------------------------------------------------
load("05_immuneSig_ICBresponse/results/oriImmuneSig_res_all.Rdata")


merged_sig_scaled <- lapply(as.list(names(merged_sig)), function(x){
     df = merged_sig[[x]]
     print(x)
     df_scaled <- apply(df[-c(1:5)], 2, scale)
     df_scaled[is.nan(df_scaled)] <- 0
     df_scaled <- cbind(df[1:5], df_scaled)
     df_scaled <- df_scaled[!is.na(df_scaled$label), ]
     df_scaled$tag <- 1
     df_scaled[df_scaled$label == "NR",]$tag = 0
     df_scaled$tag <- factor(df_scaled$tag,levels = c(1,0), order=TRUE)
     return(df_scaled)
})
names(merged_sig_scaled) <- names(merged_sig)
# feature selection within each dataset ----------------------------------------

res <- lapply(as.list(names(merged_sig_scaled)), function(x){
     print(x)
     df = merged_sig_scaled[[x]]
     df = df[-c(1:5)]
     df[is.na(df)] = 0
     feature_importance <- mlr_feature_selection(data = df,
                                                 filtermethod = "all")
     feature_importance = reshape2::dcast(feature_importance, name~filter, value.var = "value")
     feature_importance$dataset <- x
     write.csv(feature_importance, paste0(resultpath, "featureSelection_", x, ".csv"), quote = F, row.names = F)
     return(feature_importance)
})

names(res) <- names(merged_sig)
feature_importance_all <- do.call(rbind,res)

# feature selection in all datasets --------------------------------------------
merged_sig_scaled_all <- do.call(rbind, merged_sig_scaled)
merged_sig_scaled_all[is.na(merged_sig_scaled_all)] = 0

feature_importance <- mlr_feature_selection(data = merged_sig_scaled_all[-c(1:5)],
                                            filtermethod = "all")
feature_importance = reshape2::dcast(feature_importance, name~filter, value.var = "value")
write.csv(feature_importance, paste0(resultpath, "featureSelection_all.csv"), quote = F, row.names = F)


# feature selection in datasets with same cancer type --------------------------
for(cancer in unique(merged_sig_scaled_all$cancer)){
     df = merged_sig_scaled_all[merged_sig_scaled_all$cancer == cancer, ]
     feature_importance <- mlr_feature_selection(data = df[-c(1:5)],
                                                 filtermethod = "all")
     feature_importance = reshape2::dcast(feature_importance, name~filter, value.var = "value")
     # feature_importance$dataset <- x
     write.csv(feature_importance, paste0(resultpath, "featureSelection_", cancer, ".csv"), quote = F, row.names = F)

}

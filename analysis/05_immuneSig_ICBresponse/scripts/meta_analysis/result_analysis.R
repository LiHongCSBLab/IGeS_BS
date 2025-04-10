# rm(list=ls())
work_path <- "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/"
setwd(work_path)
options(stringsAsFactors = F)

dir.create("05_immuneSig_ICBresponse/results/featureSelection_summary/")

# vote -------------------------------------------------------------------------
oriImmuSig_selected <- read.table("06_results_summaryandplotting/data/immune_sig_selected_new/selected_oriImmuSig.txt",
                             sep = "\t", header = T)


list.files("05_immuneSig_ICBresponse/results/meta_results/")

# combined P-value -------------------------------------------------------------
res_merged_p <- read.csv("05_immuneSig_ICBresponse/results/meta_results/res_merged_p.csv")

# meta-analysis ----------------------------------------------------------------
res_metaRes <- read.csv("05_immuneSig_ICBresponse/results/meta_results/res_metaRes.csv")

res_metaRes <- inner_join(res_merged_p, res_metaRes)
names(res_metaRes) <- c('name', names(res_metaRes)[-1])
# mlr feature selection --------------------------------------------------------
list.files("05_immuneSig_ICBresponse/results/feature_importance_mlr/")
featureSelection_all <- read.csv("05_immuneSig_ICBresponse/results/feature_importance_mlr/featureSelection_all.csv")


featureSelection_all <- inner_join(res_metaRes,featureSelection_all)
names(featureSelection_all) <- c('immune_sig', names(featureSelection_all)[-1])
featureSelection_all <- merge(oriImmuSig_selected, featureSelection_all, by = "immune_sig",all = T)

write.table(featureSelection_all, "05_immuneSig_ICBresponse/results/featureSelection_summary/featureSelection_all.txt", 
          sep = '\t', row.names = F, quote = F)


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# cancer specific --------------------------------------------------------------
# cancer = 'SKCM'
for(cancer in c("SKCM", "LUAD", "BLCA", "KIRC", "STAD", "GBM")){

    res_merged_p <- read.csv(paste0("05_immuneSig_ICBresponse/results/meta_results/", cancer, "_res_merged_p.csv"))

    # meta-analysis ----------------------------------------------------------------
    res_metaRes <- read.csv(paste0("05_immuneSig_ICBresponse/results/meta_results/", cancer, "_res_metaRes.csv"))

    res_metaRes <- inner_join(res_merged_p, res_metaRes)
    names(res_metaRes) <- c('name', names(res_metaRes)[-1])
    # mlr feature selection --------------------------------------------------------
    list.files("05_immuneSig_ICBresponse/results/feature_importance_mlr/")
    featureSelection_all <- read.csv(paste0("05_immuneSig_ICBresponse/results/feature_importance_mlr/featureSelection_", cancer, ".csv"))


    featureSelection_all <- inner_join(res_metaRes,featureSelection_all)
    names(featureSelection_all) <- c('immune_sig', names(featureSelection_all)[-1])
    featureSelection_all <- merge(oriImmuSig_selected, featureSelection_all, by = "immune_sig",all = T)

    write.table(featureSelection_all, 
            paste0("05_immuneSig_ICBresponse/results/featureSelection_summary/featureSelection_",cancer,".txt"), 
            sep = '\t', row.names = F, quote = F)

}

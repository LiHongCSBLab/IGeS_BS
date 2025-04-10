# rm(list=ls())
setwd("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/")
options(stringsAsFactors = F)

DS_name <- c("GSE149825",
              "GSE152925",
              "GSE120500",
              "GSE160785",
              "GSE114601",
              "GSE60939",
              "GSE148538")


res_DEGreport <- lapply(as.list(DS_name), function(ds){
    df = read.csv(paste0("result_", ds, "/", ds, "_DEGanalysisReport.csv"))
    df$dataset = ds
    return(df)
})

ds = "GSE33366"
df1 = read.csv(paste0("result_", ds, "/", ds,"_BMS754807_DEGanalysisReport.csv"))
df1$dataset = paste0(ds,"_BMS754807")
names(df1) = names(res_DEGreport[[1]])
df2 = read.csv(paste0("result_", ds, "/", ds, "_Letrozole_DEGanalysisReport.csv"))
df2$dataset = paste0(ds,"_Letrozole")
names(df2) = names(res_DEGreport[[1]])
df3 = read.csv(paste0("result_", ds, "/", ds, "_Tamoxifen_DEGanalysisReport.csv"))
df3$dataset = paste0(ds,"_Tamoxifen")
names(df3) = names(res_DEGreport[[1]])
res_DEGreport <- c(res_DEGreport, list(df1,df2,df3))

ds = "GSE155923"
df1 = read.csv(paste0("result_", ds, "/Everolimus_Vehicle", ds,"_DEGanalysisReport.csv"))
df1$dataset = paste0(ds,"_Everolimus")
df2 = read.csv(paste0("result_", ds, "/NHWD870_Vehicle", ds, "_DEGanalysisReport.csv"))
df2$dataset = paste0(ds,"_NHWD870")
res_DEGreport <- c(res_DEGreport, list(df1,df2))

plot_savepath = "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/07_plotting_v2/04_cellVSmouse_metaPlusFilter/ValDS_DEGreport.Rdata"
save(res_DEGreport, file = plot_savepath)
# rm(list=ls())
work_path <- "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/"
setwd(work_path)
options(stringsAsFactors = F)

library(ggplot2)
library(ggpubr)
library(reshape2)

dir.create("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/topgene_response/")
result_dir <- "05_immuneSig_ICBresponse/results_r0.4_ITSssgsea_CombatRef/"
load( paste0(result_dir, "ITS_afterCombat17datasets_ssgseanorm.Rdata"))
datasets = unique(ITSres_selected$dataset)
response = ITSres_selected[c("SampleID","label")]
load(file = paste0(result_dir, "Exp_afterCombatRef.Rdata"))
vlSetName <- c("08_Lauss_dataset_outcome", 
               "GSE115821_pre_outcome", 
               "GSE96619_outcome",
               "Liu_dataset_CTLA4Naive_outcome", 
               "GSE176307")
datasets = c(datasets, vlSetName)
ITS_p_s_df <- read.csv("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/ITS_p_s_df.txt", sep='\t')
ITS_p_r_df <- read.csv("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/ITS_p_r_df.txt", sep='\t')

topgene_s <- ITS_p_s_df$gene_name[1:20]
topgene_r <- ITS_p_r_df$gene_name[1:20]

data_s = data_ac[, intersect(colnames(data_ac), topgene_s)]
data_s$SampleID = rownames(data_s)
data_s = merge(response, data_s, by='SampleID')
data_r = data_ac[, intersect(colnames(data_ac), topgene_r)]
data_r$SampleID = rownames(data_r)
data_r = merge(response, data_r, by='SampleID')

d_df = melt(data_s)

pdf("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/topgene_response/topgene_s_all.pdf", height=20, width = 10,onefile=T)
p <- ggboxplot(d_df, x = "label", y = "value",
          color = "label", palette = "jco",
          add = "jitter",
          facet.by = "variable", short.panel.labs = FALSE)
# Use only p.format as label. Remove method name.
p = p + stat_compare_means(label = "p.format")
print(p)
dev.off()

pdf("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/topgene_response/topgene_s_all2.pdf", height=6, width = 15,onefile=T)

p <- ggboxplot(d_df, x = "variable", y = "value",
          color = "label", palette = "jco",
          add = "jitter")
p = p + 
        stat_compare_means(aes(group = label), label =  "p.signif") +
        ggtitle(dset)
print(p)
dev.off()

pdf("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/topgene_response/topgene_s.pdf", height=6, width = 15,onefile=T)
for(dset in datasets){
    # dset = datasets[1]
    d1 = data_s[grep(dset, data_s$SampleID), ]
    # boxplot(CD53~label, d1)

    d1_df = melt(d1)
    # p <- ggboxplot(d1_df, x = "label", y = "value",
    #       color = "label", palette = "jco",
    #       add = "jitter",
    #       facet.by = "variable", short.panel.labs = FALSE)
    # # Use only p.format as label. Remove method name.
    # p + stat_compare_means(label = "p.format")

    p <- ggboxplot(d1_df, x = "variable", y = "value",
          color = "label", palette = "jco",
          add = "jitter")
    p = p + 
        stat_compare_means(aes(group = label), label =  "p.signif") +
        ggtitle(dset)
    print(p)
}
dev.off()





data_r = data_ac[, intersect(colnames(data_ac), topgene_r)]
data_r$SampleID = rownames(data_r)
data_r = merge(response, data_r, by='SampleID')
data_r = data_ac[, intersect(colnames(data_ac), topgene_r)]
data_r$SampleID = rownames(data_r)
data_r = merge(response, data_r, by='SampleID')

d_df = melt(data_r)

pdf("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/topgene_response/topgene_r_all.pdf", height=20, width = 10,onefile=T)
p <- ggboxplot(d_df, x = "label", y = "value",
          color = "label", palette = "jco",
          add = "jitter",
          facet.by = "variable", short.panel.labs = FALSE)
# Use only p.format as label. Remove method name.
p = p + stat_compare_means(label = "p.format")
print(p)
dev.off()

pdf("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/topgene_response/topgene_r_all2.pdf", height=6, width = 15,onefile=T)
p <- ggboxplot(d_df, x = "variable", y = "value",
          color = "label", palette = "jco",
          add = "jitter")
p = p + 
        stat_compare_means(aes(group = label), label =  "p.signif") +
        ggtitle(dset)
print(p)
dev.off()

pdf("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/topgene_response/topgene_r.pdf", height=6, width = 15,onefile=T)
for(dset in datasets){
    # dset = datasets[1]
    d1 = data_r[grep(dset, data_r$SampleID), ]
    # boxplot(CD53~label, d1)

    d1_df = melt(d1)
    # p <- ggboxplot(d1_df, x = "label", y = "value",
    #       color = "label", palette = "jco",
    #       add = "jitter",
    #       facet.by = "variable", short.panel.labs = FALSE)
    # # Use only p.format as label. Remove method name.
    # p + stat_compare_means(label = "p.format")

    p <- ggboxplot(d1_df, x = "variable", y = "value",
          color = "label", palette = "jco",
          add = "jitter")
    p = p + 
        stat_compare_means(aes(group = label), label =  "p.signif") +
        ggtitle(dset)
    print(p)
}
dev.off()
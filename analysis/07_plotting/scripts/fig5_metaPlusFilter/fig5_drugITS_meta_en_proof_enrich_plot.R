# rm(list=ls())
work_path <- "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/"

setwd(work_path)
options(stringsAsFactors = F)

library(org.Hs.eg.db)
library(topGO)
library(clusterProfiler)
library(msigdbr)
library(enrichplot)
library(ggplot2)
library(grid)
library(patchwork)

dir.create(paste0("07_plotting_v2/05_drugITSprofiler_CombatRef_metaPlusFilter/res_enrich_metaPlusFilter_plot/"))


drug_proof <- read.table("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/data/FDAapproved_p3Trial_IOtherapy/drug_list_proof_cancertype.txt", header = T, sep = '\t')

for(rankthreshold in rev(seq(0.5,1,0.1))) {

  # rankthreshold = 1

  resfilepath <- list.files("06_results_summaryandplotting/results_meta_filter_ITSssgsea_CombatRef_17datasets/")
  resfilepath <- resfilepath[grep("m2_meta_", resfilepath)]
  resfilepath <- resfilepath[grep("ACATP0.05", resfilepath)]
  resfilepath <- gsub("CPE_","",resfilepath)
  resfilepath <- unique(gsub("TUMERIC_","",resfilepath))

  for(respath in resfilepath){
    dir.create(paste0("07_plotting_v2/05_drugITSprofiler_CombatRef_metaPlusFilter/res_enrich_metaPlusFilter_plot/", respath))
    # respath = resfilepath[1]
    res_summary_all <- read.csv(paste0("07_plotting_v2/05_drugITSprofiler_CombatRef_metaPlusFilter/res_enrich_metaPlusFilter/",respath,"/",rankthreshold,".txt"), sep = '\t',header = T)
    df = res_summary_all[is.element(res_summary_all$ID, c('proof', 'FDA', 'ClinicalTrial_III','ClinicalTrial_beforeIII','animal')), ]
    df$cancer = paste0(df$cancer_type, "_", df$cancer_subtype,"_", df$dataset)

    strategy = gsub("_m2_meta_en_ACATP0.05/", "", unique(df$strategy))
    strategy = gsub("m2_mergeITS_Index/_", "", strategy)
    strategy = gsub("m2_mergeITS_Type/_", "", strategy)
    strategy = gsub("TUMERIC_", "", strategy)
    strategy = unique(gsub("CPE_", "", strategy))

    for(x in strategy){
        # x =strategy[4]

        df1 = df[grep(x, df$strategy), ]
        df1$ID <- factor(df1$ID, levels = c('proof', 'FDA', 'ClinicalTrial_III','ClinicalTrial_beforeIII','animal'))
        df2 = unique(df1[c('ID','cancer','NES', 'pvalue', 'p.adjust')])
        df2_sign = df2[df2$p.adjust < 0.1, ]
        p_p = ggplot(data=df2,aes(x=ID,y=cancer))+
            geom_point(aes(size= -log10(pvalue),color=NES))+
            scale_color_gradient2(low = "#41B0C3", mid="white",high = "#ed3926",
                                breaks=seq(floor(min(df2$NES)),ceiling(max(df2$NES)),0.4))+
            scale_size_continuous(range = c(1,10))+
            guides(size=F)+
            theme_bw()+
            geom_text(data=df2_sign,aes(x=ID,y=cancer,
                                    label=paste0(round(pvalue,4))))+
            theme(legend.key.height = unit(3,'cm'),
                    legend.justification = c(0,0),
                    legend.title = element_blank(),
                    axis.text.x=element_text(angle=90, hjust=1),
                    text = element_text(size=15))

        p_padj = ggplot(data=df2,aes(x=ID,y=cancer))+
            geom_point(aes(size=-log10(p.adjust),color=NES))+
            scale_color_gradient2(low = "#41B0C3", mid="white",high = "#ed3926",
                                breaks=seq(floor(min(df2$NES)),ceiling(max(df2$NES)),0.4))+
            scale_size_continuous(range = c(1,10))+
            guides(size=F)+
            theme_bw()+
            geom_text(data=df2_sign,aes(x=ID,y=cancer,
                                    label=paste0(round(p.adjust,4))))+
            theme(legend.key.height = unit(3,'cm'),
                    legend.justification = c(0,0),
                    legend.title = element_blank(),
                    axis.text.x=element_text(angle=90, hjust=1),
                    text = element_text(size=15))
        p=p_p+p_padj
        ggsave(paste0("07_plotting_v2/05_drugITSprofiler_CombatRef_metaPlusFilter/res_enrich_metaPlusFilter_plot/",respath,"/" ,x,"_",rankthreshold, ".pdf"), p, width = 16, height = 10)
        }
    }
}


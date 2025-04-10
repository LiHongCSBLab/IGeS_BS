# rm(list=ls())
work_path <- "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/"
setwd(work_path)
options(stringsAsFactors = F)
library(ggplot2)
library(pheatmap)
library(aplot)
library(reshape2)
library(ggsci)
library(see)
library(dplyr)
library(UpSetR)
library(clusterProfiler)
library(enrichplot)
source("06_results_summaryandplotting/scripts/immunesig_selection_function.R")
source("07_plotting_v2/scripts/functions/ITS_g_cancer_plot_function.R")
source("07_plotting_v2/scripts/functions/ITS_enrich.R")
source("07_plotting_v2/scripts/functions/ITS_function_detection.R")

dir.create("07_plotting_v2/03_ITScharacteristic_metaPlusFilter")

# # ITS characteristics ----------------------------------------------------------
# summary IGeS gene score

# ----------------------------------------------------------------------------
dir.create("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/")
dir.create("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/")

# ------------------------------------------------------------------------------
load( "07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITS_set.Rdata")

# ------------------------------------------------------------------------------
meta_result <- read.csv("05_immuneSig_ICBresponse/results/meta_results/res_metaRes.csv")
meta_selected <- meta_result[meta_result$Meta_Pval < 0.05, ]
meta_selected$immune_sig = meta_selected$X

load("05_immuneSig_ICBresponse/results_r0.4_ITSssgsea_CombatRef/finalmodel_17datasets_ssgseanorm/metaAnalysis_filter/en/elasticnet_model_pred.Rdata")

en_model_para = m_en$learner.model$glmnet.fit
m_en_imp <- data.frame(en_model_para$beta[,which(en_model_para$lambda == m_en$learner.model$lambda.min)])
names(m_en_imp) <- "weight"
m_en_imp$immune_sig = row.names(m_en_imp)
meta_selected = merge(meta_selected, m_en_imp, by = "immune_sig")
meta_selected = meta_selected[meta_selected$weight != 0, ]


meta_selected$flag = 'sensitive'
meta_selected[meta_selected$OR < 1, ]$flag = 'resistant'
ITS_S = meta_selected[meta_selected$flag == 'sensitive', ]$X
ITS_R = meta_selected[meta_selected$flag == 'resistant', ]$X

# set output file location -----------------------------------------------------
dir.create("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/function_per_IGeS")
savepath = "07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/function_per_IGeS/"

# ------------------------------------------------------------------------------
# 1) Functional enrichment of each ITS in each cancer type
source("07_plotting_v2/scripts/functions/ITS_g_cancer_plot_function.R")
try({
    res = lapply(as.list(names(ITS_p_set)), function(x){
        ITSlist=ITS_p_set[[x]]
        print(x)
        res = lapply(as.list(names(ITSlist)), function(y){
                    print(y)
                    glist=ITSlist[[y]]
                    if(length(glist) < 5){
                        tmp=data.frame( "ID"=NA, "Description"=NA,
                                        "GeneRatio"=NA, "BgRatio"=NA, "pvalue"=NA, 
                                        "p.adjust" =NA, "qvalue"=NA, "geneID"=NA, 
                                        "Count"=NA, ITS=y, cancer=x)
                        return(tmp)
                    }else{
                        # dir.create( paste0(savepath,"hallmark_ORA_", x, '_', y))
                        # enrichHallmark_res = ITS_GS_enrich(glist,
                        #                                 pvalueCutoff = 0.05,
                        #                                 type = "H",
                        #                                 method = 'ORA',
                        #                                 savepath =paste0(savepath,"hallmark_ORA_", x, '_', y))
                        # enrichHallmark_res2 = data.frame(enrichHallmark_res)
                        # return(enrichHallmark_res2[enrichHallmark_res2$p.adjust < 0.05, ])
                        enrichReactome_res = ITS_GS_enrich(glist,
                                                    pvalueCutoff = 0.05,
                                                    type = "Reactome",
                                                    method = 'ORA')
                        enrichReactome_res2 = as.data.frame(enrichReactome_res)
                        if(nrow(enrichReactome_res2) == 0 ){
                            tmp=data.frame( "ID"=NA, "Description"=NA,
                                            "GeneRatio"=NA, "BgRatio"=NA, "pvalue"=NA, 
                                            "p.adjust" =NA, "qvalue"=NA, "geneID"=NA, 
                                            "Count"=NA, ITS=y, cancer=x)
                            return(tmp)
                        }else{
                            tmp=enrichReactome_res2
                            tmp$ITS=y
                            tmp$cancer = x
                            return(tmp)
                            }
                    }
        })
        return(do.call(rbind, res))
    })
}, silent = T)

save(res , file=paste0(savepath, "ITS_Reactome_enrichResult.Rdata"))

# ------------------------------------------------------------------------------
# sensitive-related ITS immune regulatory 
# ------------------------------------------------------------------------------


ITS_p_s = list()
ITS_p_s_genestat <- list()

for(i in 1:length(ITS_S)){
  # i = 1
  ITSname <- ITS_S[i]
  listinput <- lapply(ITS_p_set, function(x){x[[ITSname]]})
  genelist = unique(unlist(listinput))
  listinput_df <- fromList(listinput)
  rownames(listinput_df) <- genelist
  # table(apply(listinput_df,1,sum))
  # apply(listinput_df,2,sum)
  # ITS_p_s[[i]] = unique(names(apply(listinput_df,1,sum)[which(apply(listinput_df,1,sum)>=2)]))
  ITS_p_s_genestat[[i]] <- listinput_df
  
  
  df_stat <- listinput_df # * 22
  df_stat_row <- rev(sort(rowSums(sign(df_stat))))
  df_stat_row <- as.data.frame(df_stat_row)
  names(df_stat_row) = "num_Gene"
  df_stat_row$gene_name = rownames(df_stat_row)
  df_stat_row$immune_sig <- ITSname
  
  
  ITS_p_s[[i]] <- df_stat_row # [df_stat_row$num_Gene > 1,]
  
}


# sIGeS  -----------------------------------------------------------------
load("07_plotting_v2/data/custom_gs.Rdata")

load( "07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/enrich/custom_GSEA_r_scaled.Rdata")

gs1 = gs_all[gs_all$gs_name == 'ITSp_vs_Sensitizer_PMID34980132',]$SYMBOL
gs2 = gs_all[gs_all$gs_name == 'ITSp_vs_Resistor_PMID34980132',]$SYMBOL

# 1) overall GSEA enrichment plot
library(corrplot)
pdf(paste0(savepath, "ITSp_vs_Sensitizer_gseaplot_s.pdf"), width = 8, height = 4)
# gseaplot2(res,geneSetID = "ITSp_vs_OG_oncoKB")
gseaplot(res, geneSetID = "ITSp_vs_Sensitizer_PMID34980132", by = "runningScore", title = "ITSp_vs_Sensitizer_PMID34980132")
dev.off()



# ------------------------------------------------------------------------------
res_s <- read.table("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/enrich/custom_GSEA_s_scaled.txt", sep = '\t', header = T)
res_r <- read.table("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/enrich/custom_GSEA_r_scaled.txt", sep = '\t', header = T)

GS_selected <- c("ITSp_vs_Sensitizer_PMID34980132",
                        "ITSp_vs_Resistor_PMID34980132",
                        "ITSp_vs_icb_enhancer_PMID34980132",
                        "ITSp_vs_icb_suppressor_PMID34980132",
                        "ITSp_vs_stimulatory_PMID29628290",
                        "ITSp_vs_inhibitor_PMID29628290",
                        "ITSp_vs_co_stim_inhib_PMID29628290",
                        "ITSp_vs_regulator_p_PMID34980132",
                        "ITSp_vs_regulator_n_PMID34980132")

res_s2 <- res_s[is.element(res_s$ID, GS_selected), ]
res_s2$type <- 's'
res_r2 <- res_r[is.element(res_r$ID, GS_selected), ]
res_r2$type <- 'r'

res = rbind(res_s2, res_r2)
res$signif <- ""
res[res$pvalue < 0.05, ]$signif =  "*"
# res[res$p.adjust < 0.2, ]$signif =  "*"
res[res$p.adjust < 0.1, ]$signif =  "**"
res[res$p.adjust < 0.01, ]$signif =  "***"

res$type <- factor(res$type, levels = c('s','r'))
res$ID <- factor(res$ID, levels = rev(GS_selected))


library(stringr)
res$ID <- gsub("pmid","PMID", gsub("Icb","ICB", str_to_title(gsub("ITSp_vs_", "", res$ID))))
res = res[res$ID == 'Sensitizer_PMID34980132']
p_heatmap = ggplot(res, aes(type, ID)) + 
  geom_tile(aes(fill = NES), colour = "grey", 
            size = 1,
            lwd = 2,
            linetype = 1)+
  scale_fill_gradient2(low = "#5C5DAF",mid = "white",high = "#EA2E2D", midpoint = 0) + 
  geom_text(aes(label=signif),col ="black",size = 8) +
  theme_bw() + 
  theme(axis.title.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.title.y=element_blank(), 
        axis.text.x = element_text(angle = 90, hjust = 1, size = 16), 
        axis.text.y = element_text(size = 12)) +
  labs(color = " * p < 0.05\n ** padj < 0.1\n *** padj < 0.01") 

ggsave("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/GSenrichment.pdf", p_heatmap, height = 5, width = 5)

p_bubble = ggplot(res, aes(type, ID)) + 
  geom_point(aes(color=NES, size=-log10(p.adjust)))+ # 
  scale_colour_gradient2(low="#5C5DAF",mid = 'white',high="#EA2E2D")+
  geom_text(aes(label=signif),col ="black",size = 8) +
  theme_bw() + 
  theme(axis.title.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.title.y=element_blank(), 
        axis.text.x = element_text(angle = 90, hjust = 1, size = 16), 
        axis.text.y = element_text(size = 12)) +
  labs(color =" * p < 0.05\n ** padj < 0.1\n *** padj < 0.01") 

ggsave("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/GSenrichment_bubble.pdf", p_bubble, height = 5, width = 5)


# rm(list=ls())
work_path <- "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/"
setwd(work_path)
options(stringsAsFactors = F)

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



ITS_p_s_df <- read.csv("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/ITS_p_s_df.txt", sep='\t')
ITS_p_r_df <- read.csv("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/ITS_p_r_df.txt", sep='\t')

length(union(ITS_p_s_df$gene_name, ITS_p_r_df$gene_name))
nrow(ITS_p_s_df)
nrow(ITS_p_s_df[ITS_p_s_df$weighted_rate > 0.01,])

ITS_p_s_df$num_cancer_s = ITS_p_s_df$weighted_numITS * length(ITS_S)

nrow(ITS_p_r_df)
nrow(ITS_p_r_df[ITS_p_r_df$weighted_rate > 0.01,])

names(ITS_p_s_df)=c(names(ITS_p_s_df)[1], paste0(names(ITS_p_s_df)[-1],'_s'))
names(ITS_p_r_df)=c(names(ITS_p_r_df)[1], paste0(names(ITS_p_r_df)[-1],'_r'))

ITS_p_srShare = merge(ITS_p_s_df, ITS_p_r_df, by = 'gene_name')
nrow(ITS_p_srShare)
plot(ITS_p_srShare$weighted_rate_s, ITS_p_srShare$weighted_rate_r)

load("07_plotting_v2/data/custom_gs.Rdata")
length(intersect(ITS_p_s_df$gene_name, gs_all[gs_all$gs_name == 'ITSp_vs_Sensitizer_PMID34980132',]$SYMBOL))
gs_sensitizer = gs_all[gs_all$gs_name == 'ITSp_vs_Sensitizer_PMID34980132',]
names(gs_sensitizer) <- c("gs_name","gene_name")
ITS_p_s_sensitizer = merge(ITS_p_s_df, gs_sensitizer, by = 'gene_name', all.x = T)
ITS_p_s_sensitizer[!is.element(ITS_p_s_sensitizer$gs_name,'ITSp_vs_Sensitizer_PMID34980132'),]$gs_name = 'none'
wilcox.test(ITS_p_s_sensitizer[is.element(ITS_p_s_sensitizer$gs_name,'ITSp_vs_Sensitizer_PMID34980132'),]$weighted_rate_s,
            ITS_p_s_sensitizer[!is.element(ITS_p_s_sensitizer$gs_name,'ITSp_vs_Sensitizer_PMID34980132'),]$weighted_rate_s)
boxplot(weighted_rate_s~gs_name, ITS_p_s_sensitizer)


length(intersect(ITS_p_r_df$gene_name, gs_all[gs_all$gs_name == 'ITSp_vs_Sensitizer_PMID34980132',]$SYMBOL))
ITS_p_r_sensitizer = merge(ITS_p_r_df, gs_sensitizer, by = 'gene_name', all.x = T)
ITS_p_r_sensitizer[!is.element(ITS_p_r_sensitizer$gs_name,'ITSp_vs_Sensitizer_PMID34980132'),]$gs_name = 'none'
wilcox.test(ITS_p_r_sensitizer[is.element(ITS_p_r_sensitizer$gs_name,'ITSp_vs_Sensitizer_PMID34980132'),]$weighted_rate_r,
            ITS_p_r_sensitizer[!is.element(ITS_p_r_sensitizer$gs_name,'ITSp_vs_Sensitizer_PMID34980132'),]$weighted_rate_r,
            alternative='less')
boxplot(weighted_rate_r~gs_name, ITS_p_r_sensitizer)


length(intersect(ITS_p_s_df$gene_name, gs_all[gs_all$gs_name == 'ITSp_vs_Resistor_PMID34980132',]$SYMBOL))
gs_resistor = gs_all[gs_all$gs_name == 'ITSp_vs_Resistor_PMID34980132',]
names(gs_resistor) <- c("gs_name","gene_name")
ITS_p_s_resistor = merge(ITS_p_s_df, gs_resistor, by = 'gene_name', all.x = T)
ITS_p_s_resistor[!is.element(ITS_p_s_resistor$gs_name,'ITSp_vs_Resistor_PMID34980132'),]$gs_name = 'none'
wilcox.test(ITS_p_s_resistor[is.element(ITS_p_s_resistor$gs_name,'ITSp_vs_Resistor_PMID34980132'),]$weighted_rate_s,
            ITS_p_s_resistor[!is.element(ITS_p_s_resistor$gs_name,'ITSp_vs_Resistor_PMID34980132'),]$weighted_rate_s,
            alternative='less')

boxplot(weighted_rate_s~gs_name, ITS_p_s_resistor)


length(intersect(ITS_p_r_df$gene_name, gs_all[gs_all$gs_name == 'ITSp_vs_Resistor_PMID34980132',]$SYMBOL))
gs_resistor = gs_all[gs_all$gs_name == 'ITSp_vs_Resistor_PMID34980132',]
names(gs_resistor) <- c("gs_name","gene_name")
ITS_p_r_resistor = merge(ITS_p_r_df, gs_resistor, by = 'gene_name', all.x = T)
ITS_p_r_resistor[!is.element(ITS_p_r_resistor$gs_name,'ITSp_vs_Resistor_PMID34980132'),]$gs_name = 'none'
wilcox.test(ITS_p_r_resistor[is.element(ITS_p_r_resistor$gs_name,'ITSp_vs_Resistor_PMID34980132'),]$weighted_rate_r,
            ITS_p_r_resistor[!is.element(ITS_p_r_resistor$gs_name,'ITSp_vs_Resistor_PMID34980132'),]$weighted_rate_r)
boxplot(weighted_rate_r~gs_name, ITS_p_r_resistor)


# rm(list=ls())
work_path <- "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/"
setwd(work_path)
options(stringsAsFactors = F)
library(ggplot2)
library(gg.gap)
library(pheatmap)
library(aplot)
library(reshape2)
library(ggsci)
library(see)
library(tidyverse)
library(patchwork)


source("03_drug_immuneSig_enrichment/scripts/02_GSEA_geneset_preparation_function.R")
source("03_drug_immuneSig_enrichment/scripts/02_GSEA_geneset_merge_metaAnalysis_function.R")

source("06_results_summaryandplotting/scripts/immunesig_selection_function.R")
source("06_results_summaryandplotting/scripts/ITS_analysis/ITS_function_detection.R")
source("06_results_summaryandplotting/scripts/ITS_analysis/ITSexpr_purity_detection.R")
source("06_results_summaryandplotting/scripts/ITS_analysis/ITSexpr_purity_detection.R")
source("06_results_summaryandplotting/scripts/ITS_analysis/immuneModulatorDetection.R")


dir.create("07_plotting_v2/01_metaAnalysis/")
meta_result <- read.csv("05_immuneSig_ICBresponse/results/meta_results/res_metaRes.csv")
meta_selected <- meta_result[meta_result$Meta_Pval < 0.05, ]
meta_selected$immune_sig = meta_selected$X
# meta_selected <- meta_selected[-c(grep("caf_tide", meta_selected$immune_sig),
#                                   grep("m2_tide", meta_selected$immune_sig)), ]

immunesigMeta = unique(meta_selected[c('immune_sig', 'OR', "lower_OR","upper_OR", "Meta_Pval")]) 
immunesigMeta$flag = 'sensitive'
immunesigMeta[immunesigMeta$OR < 1, ]$flag = 'resistant'

# Fig1c: mean number of genes in the signatures  --------------------------------------------
# number of genes in the signatures expressed in CCLE data --------------------------------------------
files <- list.files("06_results_summaryandplotting/immuneSigGene_ITS_found_lincs_20220414/TUMERIC_spearman_r0.4_ITSproof_2/")
files <- c(files, "SKCM_Primary.xlsx")

genestat_selected <- lapply(as.list(files), function(f){
  if(f == "SKCM_Primary.xlsx"){
    df =  as.data.frame(readxl::read_excel(paste0("06_results_summaryandplotting/immuneSigGene_ITS_found_lincs_20220414/CPE_spearman_r0.4_ITSproof_2/", f), sheet = 1))
    
  }else{
    df =  as.data.frame(readxl::read_excel(paste0("06_results_summaryandplotting/immuneSigGene_ITS_found_lincs_20220414/TUMERIC_spearman_r0.4_ITSproof_2/", f), sheet = 1))
  }
  df = df[c("immune_sig", 
            "gene_in_OriginalSig",
            "gene_in_OriginalSig_expressedinCCLE")]
  #  "gene_in_OriginalSig_lincs_overlap", 
  #  "gene_in_OriginalSig_lincs_all_overlap"
  df_selected = merge(df, immunesigMeta[c('immune_sig','flag')], by = 'immune_sig')
  df_selected$gene_in_OriginalSig_not_expressedinCCLE = df_selected$gene_in_OriginalSig - df_selected$gene_in_OriginalSig_expressedinCCLE
  
  df_selected$gene_in_OriginalSig_expressedinCCLE_rate = df_selected$gene_in_OriginalSig_expressedinCCLE / df_selected$gene_in_OriginalSig
  df_selected$gene_in_OriginalSig_not_expressedinCCLE_rate = df_selected$gene_in_OriginalSig_not_expressedinCCLE / df_selected$gene_in_OriginalSig

  df_selected = melt(df_selected[c("immune_sig", 
                                   "gene_in_OriginalSig",
                                   "gene_in_OriginalSig_not_expressedinCCLE",
                                   "gene_in_OriginalSig_expressedinCCLE",
                                   "gene_in_OriginalSig_expressedinCCLE_rate",
                                   "gene_in_OriginalSig_not_expressedinCCLE_rate")])
  df_selected$cancer = gsub(".xlsx","", f)
  return(df_selected)
})

genestat_selected2 = do.call(rbind, genestat_selected)

genestat_selected2 %>% na.omit() %>% 
  group_by(variable,immune_sig) %>% 
  summarise(mean_value=mean(value),
            sd_value=sd(value),
            min_value=min(value),
            max_value=max(value)) -> df1


df1 %>% 
  group_by(immune_sig) %>% 
  mutate(new_col=cumsum(mean_value)) -> df_ORI

df_ORI$variable<-factor(df_ORI$variable, levels = unique(df_ORI$variable))


# ITS varies across cancer types -----------------------------------------------------

files <- list.files("06_results_summaryandplotting/immuneSigGene_ITS_found_lincs_20220414/TUMERIC_spearman_r0.4_ITSproof_2/")
files <- c(files, "SKCM_Primary.xlsx")

ITSgenestat_selected <- lapply(as.list(files), function(f){
  if(f == "SKCM_Primary.xlsx"){
    df =  as.data.frame(readxl::read_excel(paste0("06_results_summaryandplotting/immuneSigGene_ITS_found_lincs_20220414/CPE_spearman_r0.4_ITSproof_2/", f), sheet = 1))
    
  }else{
    df =  as.data.frame(readxl::read_excel(paste0("06_results_summaryandplotting/immuneSigGene_ITS_found_lincs_20220414/TUMERIC_spearman_r0.4_ITSproof_2/", f), sheet = 1))
  }
  
  df = df[c("immune_sig", 
            # "gene_in_OriginalSig",
            # "gene_in_OriginalSig_expressedinCCLE",
            "gene_in_ITS_p_expressedinCCLE",
            "gene_in_OriginalSig_and_ITS_p")]
  #  "gene_in_OriginalSig_lincs_overlap", 
  #  "gene_in_OriginalSig_lincs_all_overlap"
  df_selected = merge(df, immunesigMeta[c('immune_sig','flag')], by = 'immune_sig')
  df_selected$gene_in_ITS_not_OriginalSig = df_selected$gene_in_ITS_p - df_selected$gene_in_OriginalSig_and_ITS_p
  
  df_selected = melt(df_selected[c("immune_sig", 
                                    'gene_in_ITS_p_expressedinCCLE',
                                   "gene_in_ITS_not_OriginalSig",
                                   "gene_in_OriginalSig_and_ITS_p")])
  df_selected$cancer = gsub(".xlsx","", f)
  return(df_selected)
})

ITSgenestat_selected2 = do.call(rbind, ITSgenestat_selected)

ITSgenestat_selected2 %>% na.omit() %>% 
  group_by(variable,immune_sig) %>% 
  summarise(mean_value=mean(value),
            sd_value=sd(value),
            min_value=min(value),
            max_value=max(value)) -> df1

# df1

df1 %>% 
  group_by(immune_sig) %>% 
  mutate(new_col=cumsum(mean_value)) -> df_ITS

df_ITS$variable<-factor(df_ITS$variable,
                     levels = unique(df_ITS$variable))


# density ---------------------------------------------------------------------- 

df_ORI_density = density(df_ORI[df_ORI$variable == 'gene_in_OriginalSig',]$mean_value, from = 0, to =600)
df_ORI_density=data.frame(Name = 'df_ORI_density', x = df_ORI_density$x, y = df_ORI_density$y)
df_ORI_density_upper = density(df_ORI[df_ORI$variable == 'gene_in_OriginalSig',]$mean_value + df_ORI[df_ORI$variable == 'gene_in_OriginalSig',]$sd_value, from = 0, to =600)
df_ORI_density$upper = df_ORI_density_upper$y
df_ORI_density_bottom = density(df_ORI[df_ORI$variable == 'gene_in_OriginalSig',]$mean_value - df_ORI[df_ORI$variable == 'gene_in_OriginalSig',]$sd_value, from = 0, to =600)
df_ORI_density$bottom = df_ORI_density_bottom$y


df_ORI_NOTexpinCCLE_density = density(df_ORI[df_ORI$variable == 'gene_in_OriginalSig_not_expressedinCCLE',]$mean_value, from = 0, to =600)
df_ORI_NOTexpinCCLE_density=data.frame(Name = 'df_ORI_NOTexpinCCLE_density', x = df_ORI_NOTexpinCCLE_density$x, y = df_ORI_NOTexpinCCLE_density$y)
df_ORI_NOTexpinCCLE_density_upper = density(df_ORI[df_ORI$variable == 'gene_in_OriginalSig_not_expressedinCCLE',]$mean_value + df_ORI[df_ORI$variable == 'gene_in_OriginalSig_not_expressedinCCLE',]$sd_value, from = 0, to =600)
df_ORI_NOTexpinCCLE_density$upper = df_ORI_NOTexpinCCLE_density_upper$y
df_ORI_NOTexpinCCLE_density_bottom = density(df_ORI[df_ORI$variable == 'gene_in_OriginalSig_not_expressedinCCLE',]$mean_value - df_ORI[df_ORI$variable == 'gene_in_OriginalSig_not_expressedinCCLE',]$sd_value, from = 0, to =600)
df_ORI_NOTexpinCCLE_density$bottom = df_ORI_NOTexpinCCLE_density_bottom$y


df_ORI_expinCCLE_density = density(df_ORI[df_ORI$variable == 'gene_in_OriginalSig_expressedinCCLE',]$mean_value, from = 0, to =600)
df_ORI_expinCCLE_density=data.frame(Name = 'df_ORI_expinCCLE_density', x = df_ORI_expinCCLE_density$x, y = df_ORI_expinCCLE_density$y)
df_ORI_expinCCLE_density_upper = density(df_ORI[df_ORI$variable == 'gene_in_OriginalSig_expressedinCCLE',]$mean_value + df_ORI[df_ORI$variable == 'gene_in_OriginalSig_expressedinCCLE',]$sd_value, from = 0, to =600)
df_ORI_expinCCLE_density$upper = df_ORI_expinCCLE_density_upper$y
df_ORI_expinCCLE_density_bottom = density(df_ORI[df_ORI$variable == 'gene_in_OriginalSig_expressedinCCLE',]$mean_value - df_ORI[df_ORI$variable == 'gene_in_OriginalSig_expressedinCCLE',]$sd_value, from = 0, to =600)
df_ORI_expinCCLE_density$bottom = df_ORI_expinCCLE_density_bottom$y


df_ITS_density = density(df_ITS[df_ITS$variable == 'gene_in_ITS_p_expressedinCCLE',]$mean_value, from = 0, to =600)
df_ITS_density=data.frame(Name = 'df_ITS_density', x = df_ITS_density$x, y = df_ITS_density$y)
df_ITS_density_upper = density(df_ITS[df_ITS$variable == 'gene_in_ITS_p_expressedinCCLE',]$mean_value + df_ITS[df_ITS$variable == 'gene_in_ITS_p_expressedinCCLE',]$sd_value, from = 0, to =600)
df_ITS_density$upper = df_ITS_density_upper$y
df_ITS_density_bottom = density(df_ITS[df_ITS$variable == 'gene_in_ITS_p_expressedinCCLE',]$mean_value - df_ITS[df_ITS$variable == 'gene_in_ITS_p_expressedinCCLE',]$sd_value, from = 0, to =600)
df_ITS_density$bottom = df_ITS_density_bottom$y

df_ITS_NOTinORI_density = density(df_ITS[df_ITS$variable == 'gene_in_ITS_not_OriginalSig',]$mean_value, from = 0, to =600)
df_ITS_NOTinORI_density=data.frame(Name = 'df_ITS_NOTinORI_density', x = df_ITS_NOTinORI_density$x, y = df_ITS_NOTinORI_density$y)
df_ITS_NOTinORI_density_upper = density(df_ITS[df_ITS$variable == 'gene_in_ITS_not_OriginalSig',]$mean_value + df_ITS[df_ITS$variable == 'gene_in_ITS_not_OriginalSig',]$sd_value, from = 0, to =600)
df_ITS_NOTinORI_density$upper = df_ITS_NOTinORI_density_upper$y
df_ITS_NOTinORI_density_bottom = density(df_ITS[df_ITS$variable == 'gene_in_ITS_not_OriginalSig',]$mean_value - df_ITS[df_ITS$variable == 'gene_in_ITS_not_OriginalSig',]$sd_value, from = 0, to =600)
df_ITS_NOTinORI_density$bottom = df_ITS_NOTinORI_density_bottom$y

df_ITS_inORI_density = density(df_ITS[df_ITS$variable == 'gene_in_OriginalSig_and_ITS_p',]$mean_value, from = 0, to =600)
df_ITS_inORI_density=data.frame(Name = 'df_ITS_inORI_density', x = df_ITS_inORI_density$x, y = df_ITS_inORI_density$y)
df_ITS_inORI_density_upper = density(df_ITS[df_ITS$variable == 'gene_in_OriginalSig_and_ITS_p',]$mean_value + df_ITS[df_ITS$variable == 'gene_in_OriginalSig_and_ITS_p',]$sd_value, from = 0, to =600)
df_ITS_inORI_density$upper = df_ITS_inORI_density_upper$y
df_ITS_inORI_density_bottom = density(df_ITS[df_ITS$variable == 'gene_in_OriginalSig_and_ITS_p',]$mean_value - df_ITS[df_ITS$variable == 'gene_in_OriginalSig_and_ITS_p',]$sd_value, from = 0, to =600)
df_ITS_inORI_density$bottom = df_ITS_inORI_density_bottom$y


# Fig1e: mean number of genes in original signature OR in ITS: expressed in CCLE or not -------------------------------------

df_ORI %>% filter(variable == "gene_in_OriginalSig_expressedinCCLE") %>% select(-new_col, -variable) %>% as.data.frame -> df1
# df1$min_value = df1$mean_value - df1$sd_value
# df1$max_value = df1$mean_value + df1$sd_value
names(df1) = c(names(df1)[1], gsub('value', "IS_expressedinCCLE", names(df1)[2:5]))

df_ORI %>% filter(variable == "gene_in_OriginalSig_not_expressedinCCLE") %>% select(-new_col,-variable) %>% as.data.frame -> df2
# df2$min_value = df2$mean_value - df2$sd_value
# df2$max_value = df2$mean_value + df2$sd_value
names(df2) = c(names(df2)[1], gsub('value', "not_IS_expressedinCCLE", names(df2)[2:5]))

df_ITS %>% filter(variable == "gene_in_ITS_p_expressedinCCLE") %>% select(-new_col,-variable) %>% as.data.frame -> df3
# df3$min_value = df3$mean_value - df3$sd_value
# df3$max_value = df3$mean_value + df3$sd_value
names(df3) = c(names(df3)[1], gsub('value', "ITS_expressedinCCLE", names(df3)[2:5]))

df = inner_join(inner_join(df1, df2), df3)
names(df)

p_meanNumISGene = ggplot(data = df,aes(x = mean_IS_expressedinCCLE, y = mean_not_IS_expressedinCCLE)) + 
                geom_point(size = 0.5) + 
                geom_errorbar(aes(xmin = min_IS_expressedinCCLE,xmax = max_IS_expressedinCCLE), width = .01) + 
                geom_errorbar(aes(ymin = min_not_IS_expressedinCCLE,ymax = max_not_IS_expressedinCCLE), width = .01) +
                geom_abline(slope = 1, intercept = 0, color = 'grey') +
                # geom_abline(slope = -1, intercept = 600, color = 'grey') + 
                geom_hline(aes(yintercept = 10), colour = "grey", linetype = "dashed") +
                geom_vline(aes(xintercept = 10), colour = "grey", linetype = "dashed") +
                scale_x_continuous(limits = c(0, 600), expand = c(0,0)) +
                scale_y_continuous(limits = c(0, 600), expand = c(0,0)) +
                theme_bw()
# p_meanNumISGene

p_density_ISexpCCLE = ggplot(df_ORI_expinCCLE_density, aes(x=x, color = Name, fill = Name)) +
  geom_line(aes(y = y)) +
  geom_vline(aes(xintercept = 10), colour = "grey", linetype = "dashed") +
  scale_x_continuous(limits = c(0, 600)) +
  theme_bw()

p_density_ISnotexpCCLE = ggplot(df_ORI_NOTexpinCCLE_density, aes(x=x, color = Name, fill = Name)) +
  geom_line(aes(y = y)) +
  geom_vline(aes(xintercept = 10), colour = "grey", linetype = "dashed") +
  scale_x_continuous(limits = c(0, 600)) +
  coord_flip() +
  theme_bw()

p_meanNumISGene %>%
  insert_top(p_density_ISexpCCLE, height  = .2) %>%
  insert_right(p_density_ISnotexpCCLE, width  = .2) -> p_meanNumISGene

# (p_density_ISexpCCLE / p_meanNumISGene) + p_density_ISnotexpCCLE

ggsave("07_plotting_v2/01_metaAnalysis/fig1e_meanNumISgene_dot.pdf", p_meanNumISGene, height = 6.5, width = 9)


p_meanNumISGene = ggplot(data = df,aes(x = mean_IS_expressedinCCLE, y = mean_not_IS_expressedinCCLE)) + 
                geom_point(size = 0.5) + 
                # geom_errorbar(aes(xmin = min_IS_expressedinCCLE,xmax = max_IS_expressedinCCLE), width = .01) + 
                # geom_errorbar(aes(ymin = min_not_IS_expressedinCCLE,ymax = max_not_IS_expressedinCCLE), width = .01) +
                geom_abline(slope = 1, intercept = 0, color = 'grey') +
                # geom_abline(slope = -1, intercept = 600, color = 'grey') + 
                geom_hline(aes(yintercept = 10), colour = "grey", linetype = "dashed") +
                geom_vline(aes(xintercept = 10), colour = "grey", linetype = "dashed") +
                scale_x_continuous(limits = c(0, 600), expand = c(0,0)) +
                scale_y_continuous(limits = c(0, 600), expand = c(0,0)) +
                theme_bw()
p_meanNumISGene %>%
  insert_top(p_density_ISexpCCLE, height  = .2) %>%
  insert_right(p_density_ISnotexpCCLE, width  = .2) -> p_meanNumISGene

# (p_density_ISexpCCLE / p_meanNumISGene) + p_density_ISnotexpCCLE

ggsave("07_plotting_v2/01_metaAnalysis/fig1e_meanNumISgene_dot_v2.pdf", p_meanNumISGene, height = 6.5, width = 9)

                
# p_IS_ITSgene = ggplot(data = df, aes(x = mean_IS_expressedinCCLE, y = mean_ITS_expressedinCCLE)) + 
#                 geom_point(size = 1) + 
#                 geom_errorbar(aes(ymin = min_ITS_expressedinCCLE, ymax = max_ITS_expressedinCCLE), width = .01) + 
#                 geom_errorbar(aes(xmin = min_IS_expressedinCCLE, xmax = max_IS_expressedinCCLE), width = .01) +
#                 # geom_abline(slope = 1, intercept = 0, color = 'grey') +
#                 # geom_abline(slope = -1, intercept = 600, color = 'grey') + 
#                 # geom_hline(aes(yintercept = 10), colour = "grey", linetype = "dashed") +
#                 geom_vline(aes(xintercept = 10), colour = "grey", linetype = "dashed") +
#                 scale_x_continuous(limits = c(0, 550)) +
#                 scale_y_continuous(limits = c(0, 550)) +
#                 theme_bw()

# p_IS_ITSgene
# gg.gap(plot = p_IS_ITSgene, ylim = c(0,550), segments = c(250, 500))
# ggsave("07_plotting_v2/01_metaAnalysis/fig1e_IS_ITSgene_dot.pdf", p_IS_ITSgene, height = 5, width = 5)


df_ORI %>% filter(variable == "gene_in_OriginalSig_expressedinCCLE_rate") %>% select(-new_col, -variable) %>% as.data.frame -> df1
# df1$min_value = df1$mean_value - df1$sd_value
# df1$max_value = df1$mean_value + df1$sd_value
names(df1) = c(names(df1)[1], gsub('value', "IS_expressedinCCLE_rate", names(df1)[2:5]))

df_ORI %>% filter(variable == "gene_in_OriginalSig_not_expressedinCCLE_rate") %>% select(-new_col,-variable) %>% as.data.frame -> df2
# df2$min_value = df2$mean_value - df2$sd_value
# df2$max_value = df2$mean_value + df2$sd_value
names(df2) = c(names(df2)[1], gsub('value', "not_IS_expressedinCCLE_rate", names(df2)[2:5]))

df_ITS %>% filter(variable == "gene_in_ITS_p_expressedinCCLE") %>% select(-new_col,-variable) %>% as.data.frame -> df3
# df3$min_value = df3$mean_value - df3$sd_value
# df3$max_value = df3$mean_value + df3$sd_value
names(df3) = c(names(df3)[1], gsub('value', "ITS_expressedinCCLE", names(df3)[2:5]))

df = inner_join(inner_join(df1, df2), df3)
names(df)

p_IS_ITSgeneRatio = ggplot(data = df,aes(x = mean_ITS_expressedinCCLE, y = mean_not_IS_expressedinCCLE_rate)) + 
                geom_point(size=1) + 
                geom_errorbar(aes(xmin = min_ITS_expressedinCCLE, xmax = max_ITS_expressedinCCLE), width = .01) + 
                geom_errorbar(aes(ymin = min_not_IS_expressedinCCLE_rate, ymax = max_not_IS_expressedinCCLE_rate), width = .01) +
                # geom_abline(slope = 1, intercept = 0, color = 'grey') +
                # geom_abline(slope = -1, intercept = 600, color = 'grey') + 
                # geom_hline(aes(yintercept = 10), colour = "grey", linetype = "dashed") +
                # geom_vline(aes(xintercept = 10), colour = "grey", linetype = "dashed") +
                scale_y_continuous(limits = c(0, 1)) +
                scale_x_continuous(limits = c(0, 550)) +
                theme_bw()

# p_IS_ITSgeneRatio = gg.gap(plot = p_IS_ITSgene, ylim = c(0,550), segments = c(250, 500))
p_IS_ITSgeneRatio

ggsave("07_plotting_v2/01_metaAnalysis/fig1e_IS_ITSgeneRatio_dot.pdf", p_IS_ITSgeneRatio, height = 5, width = 5)



# Fig1f: mean number of genes in original signature OR in ITS: IF GENE IN ITS ALSO IN ORIGINAL SIGNATURES -------------------------------------
df_ORIvsITS = rbind(df_ORI_density,
                    df_ITS_NOTinORI_density, 
                    df_ITS_inORI_density)

p_meanNumGene_ORIvsITS = ggplot(df_ORIvsITS, aes(x=x, color = Name, fill = Name)) +
  geom_line(aes(y=y)) +
#   geom_ribbon(aes(ymin=bottom, ymax = upper), alpha = 0.3) +
  theme_bw()
ggsave("07_plotting_v2/01_metaAnalysis/fig1d_meanNumGene_ORIvsITS.pdf", p_meanNumGene_ORIvsITS, height = 5, width = 7, limitsize = FALSE)


df_ITS %>% filter(variable == "gene_in_ITS_not_OriginalSig") %>% select(-new_col, -variable) %>% as.data.frame -> df1
# df1$min_value = df1$mean_value - df1$sd_value
# df1$max_value = df1$mean_value + df1$sd_value
names(df1) = c(names(df1)[1], gsub('value', "gene_in_ITS_not_OriginalSig", names(df1)[2:5]))

df_ITS %>% filter(variable == "gene_in_OriginalSig_and_ITS_p") %>% select(-new_col,-variable) %>% as.data.frame -> df2
# df2$min_value = df2$mean_value - df2$sd_value
# df2$max_value = df2$mean_value + df2$sd_value
names(df2) = c(names(df2)[1], gsub('value', "gene_in_OriginalSig_and_ITS_p", names(df2)[2:5]))


df = inner_join(df1, df2)
names(df)

p_IS_ITSoverlap = ggplot(data = df,aes(x = mean_gene_in_OriginalSig_and_ITS_p, y = mean_gene_in_ITS_not_OriginalSig)) + 
                geom_point(size=0.5) + 
                geom_errorbar(aes(xmin = min_gene_in_OriginalSig_and_ITS_p, xmax = max_gene_in_OriginalSig_and_ITS_p), width = .1) + 
                geom_errorbar(aes(ymin = min_gene_in_ITS_not_OriginalSig, ymax = max_gene_in_ITS_not_OriginalSig), width = .1) +
                # geom_errorbar(aes(xmin = mean_gene_in_OriginalSig_and_ITS_p - sd_gene_in_OriginalSig_and_ITS_p, xmax = mean_gene_in_OriginalSig_and_ITS_p+sd_gene_in_OriginalSig_and_ITS_p), width = .1) + 
                # geom_errorbar(aes(ymin = mean_gene_in_ITS_not_OriginalSig - sd_gene_in_ITS_not_OriginalSig, ymax = mean_gene_in_ITS_not_OriginalSig + sd_gene_in_ITS_not_OriginalSig), width = .1) +
                scale_y_continuous(limits = c(0, 250), expand = c(0,0)) +
                scale_x_continuous(limits = c(0, 550), expand = c(0,0)) +
                theme_bw()

# p_IS_ITSoverlap

p_density_ISinITS = ggplot(df_ITS_inORI_density, aes(x=x, color = Name, fill = Name)) +
  geom_line(aes(y = y)) +
  geom_vline(aes(xintercept = 10), colour = "grey", linetype = "dashed") +
  scale_x_continuous(limits = c(0, 600)) +
  theme_bw()

p_density_ISnotinITS = ggplot(df_ITS_NOTinORI_density, aes(x=x, color = Name, fill = Name)) +
  geom_line(aes(y = y)) +
  geom_vline(aes(xintercept = 10), colour = "grey", linetype = "dashed") +
  scale_x_continuous(limits = c(0, 600)) +
  coord_flip() +
  theme_bw()

p_IS_ITSoverlap %>%
  insert_top(p_density_ISinITS, height  = .2) %>%
  insert_right(p_density_ISnotinITS, width  = .2) -> p_IS_ITSoverlap
p_IS_ITSoverlap

ggsave("07_plotting_v2/01_metaAnalysis/fig1f_IS_ITSoverlap.pdf", p_IS_ITSoverlap, height = 6.5, width = 9)

p_IS_ITSoverlap = ggplot(data = df,aes(x = mean_gene_in_OriginalSig_and_ITS_p, y = mean_gene_in_ITS_not_OriginalSig)) + 
                geom_point(size=0.5) + 
                scale_y_continuous(limits = c(0, 250), expand = c(0,0)) +
                scale_x_continuous(limits = c(0, 550), expand = c(0,0)) +
                theme_bw()

p_IS_ITSoverlap %>%
  insert_top(p_density_ISinITS, height  = .2) %>%
  insert_right(p_density_ISnotinITS, width  = .2) -> p_IS_ITSoverlap
p_IS_ITSoverlap

ggsave("07_plotting_v2/01_metaAnalysis/fig1f_IS_ITSoverlap_v2.pdf", p_IS_ITSoverlap, height = 6.5, width = 9)



# fig 1g: Correlation between ITS and original signatures ------------------------------
files <- list.files("02_tumorGene_immuneSig_correlation/immuneSig_ITS_correction/")

ITScor_selected <- lapply(as.list(files), function(f){
  if(file.exists(paste0("02_tumorGene_immuneSig_correlation/immuneSig_ITS_correction/", f, "/cor_sig_ssgsea.csv"))){
    df =  read.csv(paste0("02_tumorGene_immuneSig_correlation/immuneSig_ITS_correction/", f, "/cor_sig_ssgsea.csv"))
    df_selected = merge(df, immunesigMeta[c('immune_sig','flag')], by = 'immune_sig')
    df_selected$cancer = f
    return(df_selected)
  }
})

ITScor_selected2 = do.call(rbind, ITScor_selected)

ITScor_selected2[c('immune_sig','r','flag','cancer')] %>% na.omit() %>% 
  group_by(immune_sig) %>% 
  summarise(mean_value=mean(r),
            sd_value=sd(r)) -> df_cor_ORIvsITS

df_cor_ORIvsITS$upper = df_cor_ORIvsITS$mean_value + df_cor_ORIvsITS$sd_value
df_cor_ORIvsITS$bottom = df_cor_ORIvsITS$mean_value - df_cor_ORIvsITS$sd_value
df_cor_ORIvsITS <- df_cor_ORIvsITS[order(df_cor_ORIvsITS$mean_value, decreasing = T), ]
df_cor_ORIvsITS$immune_sig<-factor(df_cor_ORIvsITS$immune_sig, levels = unique(df_cor_ORIvsITS$immune_sig))
p_cor_ORIvsITS = ggplot(df_cor_ORIvsITS, aes(x=immune_sig)) +
  geom_point(aes(x=immune_sig, y=mean_value)) +
  geom_errorbar(aes(x=immune_sig, ymin=bottom, ymax = upper), width = 0.5) +
  # geom_ribbon(aes(x=immune_sig, ymin=bottom, ymax = upper), alpha = 0.3) +
  geom_hline(aes(yintercept = 0.5), colour = "#636363", linetype = "dashed") +
  geom_hline(aes(yintercept = 0.75), colour = "#b4b4b4", linetype = "dashed") +
  theme_bw() +
  theme(axis.text.x = element_text(hjust = 1, size = 7, angle = 90)) 
ggsave("07_plotting_v2/01_metaAnalysis/fig1e_cor_ORIvsITS.pdf", p_cor_ORIvsITS, height = 7, width = 5, limitsize = FALSE)
ggsave("07_plotting_v2/01_metaAnalysis/fig1e_cor_ORIvsITS2.pdf", p_cor_ORIvsITS, height = 5, width = 4, limitsize = FALSE)


# p_ITScor <- ggplot(data=df_cor_ORIvsITS,aes(x=immune_sig,
#                                 y=mean_value))+
#   geom_bar(position = "stack",stat="identity", fill = "#9292fd", alpha = .5, width = 0.5) +
#   geom_errorbar(aes(ymin=mean_value-sd_value,
#                     ymax=mean_value+sd_value),
#                 width=0.1)+
#   geom_hline(aes(yintercept = 0.5), colour = "#636363", linetype = "dashed") +
#   geom_hline(aes(yintercept = 0.75), colour = "#b4b4b4", linetype = "dashed") +
#   scale_y_continuous(limits = c(0,1),
#                      breaks = seq(0,1,0.25))+
#   # scale_fill_material_d()+
#   theme_bw()+coord_flip()+
#   labs(x=NULL,y=NULL)



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
library(tidyverse)



source("03_drug_immuneSig_enrichment/scripts/02_GSEA_geneset_preparation_function.R")
source("03_drug_immuneSig_enrichment/scripts/02_GSEA_geneset_merge_metaAnalysis_function.R")

source("06_results_summaryandplotting/scripts/immunesig_selection_function.R")
source("06_results_summaryandplotting/scripts/ITS_analysis/ITS_function_detection.R")
source("06_results_summaryandplotting/scripts/ITS_analysis/ITSexpr_purity_detection.R")
source("06_results_summaryandplotting/scripts/ITS_analysis/ITSexpr_purity_detection.R")
source("06_results_summaryandplotting/scripts/ITS_analysis/immuneModulatorDetection.R")


dir.create("07_plotting_v2/01_metaAnalysis/")


immuneSigInfoPath = "06_results_summaryandplotting/data/immene_cell_hieracy/"

immunesigInfo <- as.data.frame(readxl::read_excel(paste0(immuneSigInfoPath, "immune_sig_hieracy_new.xlsx"), sheet = 1))
immunesigInfo <- immunesigInfo[c(1,4,5:9)]
immunesigInfo$immune_sig = tolower(immunesigInfo$immune_sig)
immunesigInfo = immunesig_name_unify(immunesigInfo)



meta_result <- read.csv("05_immuneSig_ICBresponse/results/meta_results/res_metaRes.csv")
meta_selected <- meta_result[meta_result$Meta_Pval < 0.05, ]
meta_selected$immune_sig = meta_selected$X
# meta_selected <- meta_selected[-c(grep("caf_tide", meta_selected$immune_sig),
#                                   grep("m2_tide", meta_selected$immune_sig)), ]


# stat -------------------------------------------------------------------------
immunesigInfo_used = immunesigInfo[is.element(immunesigInfo$immune_sig, meta_result$X), ]
length(intersect(meta_result$X, immunesigInfo_used$immune_sig))
table(immunesigInfo_used$Classification)


# ------------------------------------------------------------------------------
# 2. meta analysis [P < 0.05  & (OR > 1.5 | OR < 1)] --------------------------------------------------
# ------------------------------------------------------------------------------
immunesig_icbResponse_all_results <- read.csv("05_immuneSig_ICBresponse/results/meta_results/immunesig_icbResponse_all_results.csv", row.names = 1)

immunesig_icbResponse_all_results <- merge(immunesig_icbResponse_all_results, meta_selected, by = "immune_sig")
immunesig_icbResponse_s <- immunesig_icbResponse_all_results[is.element(immunesig_icbResponse_all_results$immune_sig, meta_selected[meta_selected$OR > 1, ]$X), ]
immunesig_icbResponse_s$flag <- 'sensitive'
immunesig_icbResponse_s <- immunesig_icbResponse_s[order(immunesig_icbResponse_s$OR, decreasing = T), ]
immunesig_icbResponse_r <- immunesig_icbResponse_all_results[is.element(immunesig_icbResponse_all_results$immune_sig, meta_selected[meta_selected$OR < 1, ]$X), ]
immunesig_icbResponse_r$flag <- 'resistant'
immunesig_icbResponse_r <- immunesig_icbResponse_r[order(immunesig_icbResponse_r$OR, decreasing = T), ]

immunesig_icbResponse_selected <- rbind(immunesig_icbResponse_s, immunesig_icbResponse_r)
# write.table(immunesig_icbResponse_selected, "07_plotting_v2/01_metaAnalysis/table/immunesig_icbResponse_selected.txt", sep = '\t', row.names = F, quote = F)


data <- immunesig_icbResponse_selected[immunesig_icbResponse_selected$seTE < 10,]
data <- data[!is.na(data$immune_sig), ]
data <- merge(data, immunesigInfo, by = "immune_sig", all.x = T)

# immunesigCL <- read.csv("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/03_drug_immuneSig_enrichment/data/ITS_hclust_Otsuka_Ochiai_TUMERIC_spearman_r_0.4_pan_metaP_0.05_genevote_0.5/SKCM_immunesig_ITS_index.txt", sep = '\t')
immunesigAnnot <- merge(unique(meta_selected), immunesigInfo, by='immune_sig')
write.table(immunesigAnnot, "07_plotting_v2/01_metaAnalysis/metaSelected_immunesigAnnot.txt", sep = '\t', row.names = F, quote = F)

immunesigAnnot <- merge(unique(meta_selected[c("immune_sig")]), immunesigInfo, by='immune_sig')
immunesigCL <- immunesigAnnot[c("immune_sig", "Type")]
names(immunesigCL) <- c("immune_sig", "Index_CL")


data <- merge(data, immunesigCL, by = "immune_sig")
# select signatures for plotting
# data_selected <- lapply(as.list(unique(data$Index_CL)), function(x){
#   tmp = data[data$Index_CL == x, ]
#   tmp = tmp[tmp$OR == max(tmp$OR), ]
#   return(tmp)
# })

# data <- unique(do.call(rbind, data_selected))
names(data)
# data <- data[data$OR >1.5 | data$OR < 1,]
data$TE_plot <- data$TE
if(nrow(data[data$TE_plot >= 3,]) > 0) data[data$TE_plot > 3,]$TE_plot = 5
if(nrow(data[data$TE_plot < -3,]) > 0) data[data$TE_plot < -3,]$TE_plot = -5

data <- data %>%
  mutate(text = case_when( 
    (auc>0.6 & TE_Pval > 0.05) ~ paste(" "), 
    (auc<0.6 & TE_Pval > 0.05) ~ paste(" "), 
    (auc>0.6 & TE_Pval < 0.05 & TE_Pval > 0.01) ~ paste("*"), 
    (auc<0.6 & TE_Pval < 0.05 & TE_Pval > 0.01) ~ paste("*"), 
    (auc>0.6 & TE_Pval < 0.01) ~ paste("**"),
    (auc<0.6 & TE_Pval < 0.01) ~ paste("**")))
    
data <- data[order(data$cancer), ]
data$dataset <- factor(data$dataset, levels = unique(data$dataset)) 
data_s <- data[data$flag == "sensitive", ]
data_s <- data_s[order(data_s$OR, decreasing = T), ]
data_s <- data_s[order(data_s$SubType), ]
data_s <- data_s[order(data_s$Type), ]
data_r <- data[data$flag == "resistant", ]
data_r <- data_r[order(data_r$OR, decreasing = T), ]
data_r <- data_r[order(data_r$SubType), ]
data_r <- data_r[order(data_r$Type), ]
data <- rbind(data_s, data_r)
data$immune_sig <- factor(data$immune_sig, levels = rev(unique(data$immune_sig))) 
data$immune_sig <- factor(data$immune_sig, levels = (unique(data$immune_sig))) 

data$TE_color = cut(data$TE_plot, 2)

p1 = ggplot(data, aes(immune_sig, dataset,  fill = TE_plot)) + 
  geom_tile(aes(fill = TE_plot), colour = "white", 
            size = 1, lwd = 2, linetype = 1)+
  scale_fill_gradient2(low = "#228383",mid = "white",high = "#B33333", midpoint = 0) + 
  geom_text(aes(label=text),col ="black",size = 3) +
  theme_bw() + 
  theme(axis.title.x=element_blank(), # remove title
        axis.ticks.x=element_blank(), # remove x axis
        axis.title.y=element_blank(), # remove y axis
        # axis.text.x = element_text(angle = 90, hjust = 1, size = 14, face = "bold"), 
        axis.text.x = element_blank(),
        axis.text.y = element_text(hjust = 1, size = 14, face = "bold")) + # remove text on y axis
  labs(fill =paste0("+ auc > 0.6", " * 0.01 < p < 0.05","\n\n","** p < 0.01","\n\n","TE","\n\n")) +   # 修改 legend 内容
  scale_x_discrete(position = "bottom") + # set x axis on top
  scale_y_discrete(position = "left") # set x axis on right
# p1


p2 = ggplot(data, aes(immune_sig, dataset,)) + 
  # geom_point(aes(size=abs(TE_plot), color = factor(TE_color)), shape=15) +
  geom_tile(aes(fill = factor(TE_color, levels=unique(TE_color))), colour = "white", size = 1, lwd = 2, linetype = 1) +
  scale_color_manual(values = c("#33808099","#ea2d2d99")) +  
  geom_text(aes(label=text),col ="black",size = 3) +
  theme_bw() + 
  theme(axis.title.x=element_blank(), # remove title
        axis.ticks.x=element_blank(), # remove x axis
        axis.title.y=element_blank(), # remove y axis
        axis.text.x = element_blank(),
        axis.text.y = element_text(hjust = 1, size = 14, face = "bold")) + # remove text on y axis
  labs(fill =paste0(" * 0.01 < p < 0.05","\n\n","** p < 0.01","\n\n","TE","\n\n")) +   # 修改 legend 内容
  scale_x_discrete(position = "bottom") + # set x axis on top
  scale_y_discrete(position = "left") # set x axis on right
# p2

# dataset annotation --------------------------------------------------
dataset_cancer = unique(data[c('dataset', 'cancer')]) %>% 
  # mutate(group=cancer) %>%
  mutate(p="cancer")
p_cancer = ggplot(dataset_cancer, aes(dataset, y=p, fill=cancer))+
  geom_tile() + 
  scale_y_discrete(position="right") +
  # theme_minimal() +
  theme_bw() + # 不要背景
  xlab(NULL) + ylab(NULL) +
  theme(axis.text.x = element_text(angle = 90),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  labs(fill = "cancer")+ 
  scale_fill_jama() +
  coord_flip() + 
  scale_x_discrete(position = "bottom")

dataset_treatment = unique(data[c('dataset', 'treatment_type')]) %>% 
  # mutate(group=cancer) %>%
  mutate(p="treatment_type")
p_treatment = ggplot(dataset_treatment, aes(dataset, y=p, fill=treatment_type))+
  geom_tile() + 
  scale_y_discrete(position="right") +
  theme_bw() + 
  xlab(NULL) + ylab(NULL) +
  theme(axis.text.x = element_text(angle = 90),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())+
  labs(fill = "cancer") +
  scale_fill_jco() +
  coord_flip() + 
  scale_x_discrete(position = "bottom")
  

# immune signature annotation --------------------------------------------------
immunesigAnnot <- merge(unique(data[c("immune_sig")]), immunesigInfo, by='immune_sig')
immunesigType = unique(immunesigAnnot[c('immune_sig', 'Type')]) %>% 
  mutate(p="Type")
p_Type = ggplot(immunesigType, aes(x=immune_sig, y=p, fill=Type))+
  geom_tile() + 
  scale_y_discrete(position="right") +
  theme_bw() + 
  xlab(NULL) + ylab(NULL) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10),
        axis.ticks.x=element_blank())+
  labs(fill = "Type") +
  scale_fill_d3(palette = "category20") 

immunesigSubType = unique(immunesigAnnot[c('immune_sig', 'SubType')]) %>% 
  mutate(p="SubType")
p_SubType = ggplot(immunesigSubType, aes(x=immune_sig, y=p, fill=SubType))+
  geom_tile() + # show.legend = F
  scale_y_discrete(position="right") +
  theme_bw() + 
  xlab(NULL) + ylab(NULL) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x=element_blank())+
  labs(fill = "SubType")

# immunesigClassification = unique(immunesigAnnot[c('immune_sig', 'Classification')]) %>% 
#   mutate(p="Classification")

# immunesigFlag = unique(data[c('immune_sig', 'flag')]) %>% 
#   mutate(p="flag")
# p_flag = ggplot(immunesigFlag, aes(x=p, y= immune_sig,  fill=flag))+
#   geom_tile() + 
#   scale_y_discrete(position="right") +
#   theme_bw() + 
#   xlab(NULL) + ylab(NULL) +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10),
#         axis.text.y = element_blank())+
#   labs(fill = "flag") +
#   scale_fill_manual(values = c("#338080","#B33333"))


# meta analysis forest plot ----------------------------------------------------
immunesigMeta = unique(data[c('immune_sig', 'OR', "lower_OR","upper_OR", "Meta_Pval", "flag")]) 
dotCOLS = c("#33808099","#e79999")
barCOLS = c("#338080","#B33333")

p_metaForest <- ggplot(immunesigMeta, aes(x=immune_sig, y=OR, ymin=lower_OR, ymax=upper_OR,col=flag,fill=flag)) + 
  geom_linerange(size=2,position=position_dodge(width = 0.1)) +
  geom_hline(yintercept=1, lty=2) +
  geom_point(size=2, shape=21, colour="white", stroke = 0.5, position=position_dodge(width = 0.1)) +
  scale_fill_manual(values=barCOLS)+
  scale_color_manual(values=dotCOLS)+
  scale_x_discrete(name="immune_sig") +
  scale_y_continuous(name="Odds ratio", limits = c(0,3))  + 
  theme_bw() + 
  xlab(NULL) + ylab(NULL)+
  theme(axis.title.x=element_blank(), 
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank()) # +
  # coord_flip()
# p_metaForest

p_all1 = p1 %>%
  insert_right(p_treatment, width = .01)%>%
  insert_right(p_cancer, width = .01)%>%
  insert_bottom(p_SubType, height  = .05) %>%
  insert_bottom(p_Type, height  = .05) %>%
  # insert_bottom(p_Classification, height  = .05) %>%
  # insert_bottom(p_flag, height = .05)
  insert_top(p_metaForest, height  = .3)
# p_all1
ggsave("07_plotting_v2/01_metaAnalysis/fig1b_v1.pdf", p_all1, height = 9, width = 17, limitsize = FALSE)
ggsave("07_plotting_v2/01_metaAnalysis/fig1b_v1_2.pdf", p_all1, height = 12, width = 24, limitsize = FALSE)

p_all2 = p2 %>%
  insert_right(p_treatment, width = .01)%>%
  insert_right(p_cancer, width = .01)%>%
  insert_bottom(p_SubType, height  = .05) %>%
  insert_bottom(p_Type, height  = .05) %>%
  # insert_bottom(p_Classification, height  = .05) %>%
  # insert_bottom(p_flag, height = .05)
  insert_top(p_metaForest, height  = .3)
# p_all2
ggsave("07_plotting_v2/01_metaAnalysis/fig1b_v2.pdf", p_all2, height = 9, width = 17, limitsize = FALSE)
ggsave("07_plotting_v2/01_metaAnalysis/fig1b_v2_2.pdf", p_all2, height = 12, width = 24, limitsize = FALSE)


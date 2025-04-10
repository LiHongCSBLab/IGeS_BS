# rm(list=ls())
work_path <- "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/"
setwd(work_path)
options(stringsAsFactors = F)

# library(reshape2)
# library(ggplotify)
# library(grid)
# library(cowplot)
# library(parallel)
# library(patchwork)
# library(dplyr)
# library(stringr)
# library(data.table)
# require(GGally)
# library(plot3D)
# library(fmsb)
# library(readxl)
library(ggplot2)
library(dplyr)
library(UpSetR)
library(clusterProfiler)
library(enrichplot)
# source("06_results_summaryandplotting/scripts/plot_for_drug_function.R")
# source("06_results_summaryandplotting/scripts/ITS_analysis/ITS_gene_stat_function.R")
# source("06_results_summaryandplotting/scripts/immunesig_selection_function.R")
# source("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/scripts/drugITGSscore_lincs.R")
source("07_plotting_v2/scripts/functions/ITS_enrich.R")
source("07_plotting_v2/scripts/functions/ITS_function_detection.R")
source("06_results_summaryandplotting/scripts/immunesig_selection_function.R")
source("07_plotting_v2/scripts/functions/ITS_enrich.R")
source("07_plotting_v2/scripts/functions/ITS_function_detection.R")
source("07_plotting_v2/scripts/functions/ITS_g_cancer_plot_function.R")


# ------------------------------------------------------------------------------
# Drug target analysis
# ------------------------------------------------------------------------------
load( "07_plotting_v2/data/drugTarget_merged3DB.Rdata")
# drugTarget_all
# drugTarget_FDA
drugTarget_FDAcancer <- drugTarget_FDA[which(drugTarget_FDA$ifCancer_drugbank == 'yes'),]

DTall <- unique(drugTarget_all[c('drugname', 'drugname_lower', 'genename')])

# ------------------------------------------------------------------------------
ITS_p_s_df <- read.table("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/ITS_p_s_df.txt", sep = '\t', header = T)
ITS_p_r_df <- read.table("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/ITS_p_r_df.txt", sep = '\t', header = T)

# ------------------------------------------------------------------------------
# drug enrichment according to target drugs 
library(org.Hs.eg.db)
library(topGO)
library(clusterProfiler)
library(msigdbr)
library(enrichplot)
DTall2 <- bitr(DTall$genename, 
               fromType = "SYMBOL", 
               toType = c("ENTREZID",  "SYMBOL"), 
                OrgDb = org.Hs.eg.db)
DTall$SYMBOL <- DTall$genename
DTall2 <- merge(DTall, DTall2, by = "SYMBOL")


D_Enrich_ITS_s <- lapply(as.list(unique(DTall$drugname_lower)), function(drug){
    print(drug)
    dt <- unique(DTall[DTall$drugname_lower == drug, ])
    dt <- dt[c('drugname_lower','genename','SYMBOL')]

    if(length(intersect(ITS_p_s_df$gene_name, dt$genename))>0){

        dt_its <- ITS_p_s_df[is.element(ITS_p_s_df$gene_name, dt$genename), ]

        dt_score <- data.frame(drugname = drug,
                                naive_score = length(intersect(unique(dt$genename), unique(ITS_p_s_df$gene_name)))/ length(unique(dt$genename)),
                                weighted_score = sum(dt_its$weighted_rate)/length(unique(dt$genename)))

        gs <- dt[c(1,3)]
        names(gs) <- c("gs_name", "SYMBOL")
        tmp2 = as.vector(scale(ITS_p_s_df$weighted_rate))
        names(tmp2) = ITS_p_s_df$gene_name
        res <- GSEA(tmp2, pvalueCutoff = 1, TERM2GENE=gs)
        res2 = data.frame(res)
        if(nrow(res2) > 0){
            dt_score$weighted_NES <- res2$NES
            dt_score$weighted_pvalue <- res2$pvalue
            dt_score$weighted_p.adjust <- res2$p.adjust
        }else{
            dt_score$weighted_NES <- NA
            dt_score$weighted_pvalue <- NA
            dt_score$weighted_p.adjust <- NA
        }
        
          
      gs2 <- bitr(gs$SYMBOL, 
                      fromType = "SYMBOL", 
                      toType = c("ENTREZID",  "SYMBOL"), 
                      OrgDb = org.Hs.eg.db)
      gs3 <- merge(gs, gs2, by = "SYMBOL")
      gs3 <- gs3[c(2,1,3)]
      gs3 <- gs3[order(gs3$gs_name),]      

      res <- enricher(DTall2$ENTREZID, 
                      pvalueCutoff = pvalueCutoff, 
                      minGSSize = 0,
                      maxGSSize = 10000,
                      TERM2GENE = gs3)
        if(is.null(res)){
            dt_score$ORA_Count <- NA
            dt_score$ORA_pvalue <- NA
            dt_score$ORA_p.adjust <- NA

        }else{
            dt_score$ORA_Count <- res$Count
            dt_score$ORA_pvalue <- res$pvalue
            dt_score$ORA_p.adjust <- res$p.adjust

        }
        

    }else{
        dt_score <- data.frame(drugname = drug,
                               naive_score = NA,
                               weighted_score = NA,
                               weighted_NES = NA,
                               weighted_pvalue = NA,
                               weighted_p.adjust = NA,
                               ORA_Count = NA,
                               ORA_pvalue = NA,
                               ORA_p.adjust = NA)
    }
    return(dt_score)           
})


D_Enrich_ITS_s_score <- do.call(rbind, D_Enrich_ITS_s)
D_Enrich_ITS_s_score <- D_Enrich_ITS_s_score[order(D_Enrich_ITS_s_score$naive_score,decreasing = T),]
D_Enrich_ITS_s_score$naive_rank <- 1:nrow(D_Enrich_ITS_s_score)
D_Enrich_ITS_s_score$naive_rankTile <- D_Enrich_ITS_s_score$naive_rank/nrow(D_Enrich_ITS_s_score)
D_Enrich_ITS_s_score <- D_Enrich_ITS_s_score[order(D_Enrich_ITS_s_score$weighted_score,decreasing = T),]
D_Enrich_ITS_s_score$weighted_rank <- 1:nrow(D_Enrich_ITS_s_score)
D_Enrich_ITS_s_score$weighted_rankTile <- D_Enrich_ITS_s_score$weighted_rank/nrow(D_Enrich_ITS_s_score)

D_Enrich_ITS_s_score$ifFDA <- 'no'
D_Enrich_ITS_s_score[is.element(D_Enrich_ITS_s_score$drugname, drugTarget_FDA$drugname_lower),]$ifFDA <- 'FDA'
D_Enrich_ITS_s_score[is.element(D_Enrich_ITS_s_score$drugname, drugTarget_FDAcancer$drugname_lower),]$ifFDA <- 'FDAcancer'
head(D_Enrich_ITS_s_score[D_Enrich_ITS_s_score$ifFDA == 'FDAcancer',])

write.table(D_Enrich_ITS_s_score, "07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/drug_merge3DB/D_Enrich_ITS_s_score.txt", sep = '\t', quote = F, row.names = F)

# ITS_r
D_Enrich_ITS_r <- lapply(as.list(unique(DTall$drugname_lower)), function(drug){
    print(drug)
    dt <- unique(DTall[DTall$drugname_lower == drug, ])
    dt <- dt[c('drugname_lower','genename','SYMBOL')]
    if(length(intersect(ITS_p_r_df$gene_name, dt$genename))>0){

        dt_its <- ITS_p_r_df[is.element(ITS_p_r_df$gene_name, dt$genename), ]
        dt_score <- data.frame(drugname = drug,
                                naive_score = length(intersect(unique(dt$genename), unique(ITS_p_r_df$gene_name)))/ length(unique(dt$genename)),
                                weighted_score = sum(dt_its$weighted_rate)/length(unique(dt$genename)))
        gs <- dt[c(1,3)]
        names(gs) <- c("gs_name", "SYMBOL")
        tmp2 = as.vector(scale(ITS_p_r_df$weighted_rate))
        names(tmp2) = ITS_p_r_df$gene_name
        res <- GSEA(tmp2, pvalueCutoff = 1, TERM2GENE=gs)
        res2 = data.frame(res)
        if(nrow(res2) > 0){
            dt_score$weighted_NES <- res2$NES
            dt_score$weighted_pvalue <- res2$pvalue
            dt_score$weighted_p.adjust <- res2$p.adjust
        }else{
            dt_score$weighted_NES <- NA
            dt_score$weighted_pvalue <- NA
            dt_score$weighted_p.adjust <- NA
        }

          gs2 <- bitr(gs$SYMBOL, 
                      fromType = "SYMBOL", 
                      toType = c("ENTREZID",  "SYMBOL"), 
                      OrgDb = org.Hs.eg.db)
      gs3 <- merge(gs, gs2, by = "SYMBOL")
      gs3 <- gs3[c(2,1,3)]
      gs3 <- gs3[order(gs3$gs_name),]      

      res <- enricher(DTall2$ENTREZID, 
                      pvalueCutoff = pvalueCutoff, 
                      minGSSize = 0,
                      maxGSSize = 10000,
                      TERM2GENE = gs3)
          if(is.null(res)){
            dt_score$ORA_Count <- NA
            dt_score$ORA_pvalue <- NA
            dt_score$ORA_p.adjust <- NA

        }else{
            dt_score$ORA_Count <- res$Count
            dt_score$ORA_pvalue <- res$pvalue
            dt_score$ORA_p.adjust <- res$p.adjust

        }
        
    }else{
        dt_score <- data.frame(drugname = drug,
                               naive_score = NA,
                               weighted_score = NA,
                               weighted_NES = NA,
                               weighted_pvalue = NA,
                               weighted_p.adjust = NA,
                               ORA_Count = NA,
                               ORA_pvalue = NA,
                               ORA_p.adjust = NA)
    }
    return(dt_score)           
})
D_Enrich_ITS_r_score <- do.call(rbind, D_Enrich_ITS_r)


D_Enrich_ITS_r_score <- D_Enrich_ITS_r_score[order(D_Enrich_ITS_r_score$naive_score,decreasing = T),]
D_Enrich_ITS_r_score$naive_rank <- 1:nrow(D_Enrich_ITS_r_score)
D_Enrich_ITS_r_score$naive_rankTile <- D_Enrich_ITS_r_score$naive_rank/nrow(D_Enrich_ITS_r_score)
D_Enrich_ITS_r_score <- D_Enrich_ITS_r_score[order(D_Enrich_ITS_r_score$weighted_score,decreasing = T),]
D_Enrich_ITS_r_score$weighted_rank <- 1:nrow(D_Enrich_ITS_r_score)
D_Enrich_ITS_r_score$weighted_rankTile <- D_Enrich_ITS_r_score$weighted_rank/nrow(D_Enrich_ITS_r_score)

D_Enrich_ITS_r_score$ifFDA <- 'no'
D_Enrich_ITS_r_score[is.element(D_Enrich_ITS_r_score$drugname, drugTarget_FDA$drugname_lower),]$ifFDA <- 'FDA'
D_Enrich_ITS_r_score[is.element(D_Enrich_ITS_r_score$drugname, drugTarget_FDAcancer$drugname_lower),]$ifFDA <- 'FDAcancer'
head(D_Enrich_ITS_r_score[D_Enrich_ITS_r_score$ifFDA == 'FDAcancer',])
write.table(D_Enrich_ITS_r_score, "07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/drug_merge3DB/D_Enrich_ITS_r_score.txt", sep = '\t', quote = F, row.names = F)


# ------------------------------------------------------------------------------
# dot rank plot, 
# marking top 0.05%, FDA approved anti-cancer and other diseases treatment seperately
library(gg.gap)
library(ggplot2)
library(patchwork)

library(ggrepel)
# dat_s
D_Enrich_ITS_s_score <- read.csv("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/drug_merge3DB/D_Enrich_ITS_s_score.txt", sep = '\t', header = T)
D_Enrich_ITS_r_score <- read.csv("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/drug_merge3DB/D_Enrich_ITS_r_score.txt", sep = '\t', header = T)

dat_s <- D_Enrich_ITS_s_score
dat_s$hits <- 'no'
dat_s[dat_s$ifFDA == 'FDA' & dat_s$weighted_rankTile < 0.01, ]$hits <- 'FDA'
dat_s[dat_s$ifFDA == 'FDAcancer' & dat_s$weighted_rankTile < 0.01, ]$hits <- 'FDAcancer'
# dat_s[which(dat_s$weighted_pvalue < 0.05 & dat_s$weighted_NES > 0),]$hits <- 'Target enriched'

dat_s = dat_s %>%
  arrange(desc(weighted_score)) %>%
  mutate("x" = row_number())
dat_s = dat_s[!is.na(dat_s$weighted_score),]
dat_s$x_prop <- dat_s$x / max(dat_s$x)
dat_s_shown <- dat_s[-grep("no", dat_s$hits), ]
# dat_s_shown$hits <- factor(dat_s_shown$hits, levels = c('FDA','FDAcancer', 'Target enriched'))
dat_s_shown$hits <- factor(dat_s_shown$hits, levels = c('FDA','FDAcancer'))
# dat_s_shown <- dat_s[which(dat_s$weighted_pvalue < 0.05),]
# dat_s_shown$hits <- dat_s_shown$ifFDA

p_s_rank <- ggplot(dat_s, aes(x = x_prop, y = weighted_score, color=hits))+
  geom_point(color="#B33333", size = 0.5)+
  geom_vline(xintercept = 0, color = 'grey30')+
  geom_hline(yintercept = 0, color = 'grey30')+
  geom_vline(xintercept = 0.01, color = 'grey', linetype = 'dashed')+
  geom_vline(xintercept = 0.05, color = 'grey', linetype = 'dashed')+
  geom_vline(xintercept = 0.1, color = 'grey', linetype = 'dashed')+
  geom_rug(data = subset(dat_s, ifFDA == 'FDA'), color = "#ee766f", sides = "b",length = unit(3,"mm")) +
  geom_rug(data = subset(dat_s, ifFDA == 'FDAcancer'), color = "#2eb7be", sides = "b",length = unit(3,"mm")) +
  # geom_rug(data = subset(dat_s, ifFDA == 'no'), color = "white", sides = "b") +
  theme_bw()+
  labs(x = "Ranked drugs", y = "Potential")+
  geom_text_repel(inherit.aes = F, 
                  data = dat_s_shown, 
                  aes(x = x_prop, y = weighted_score, label = drugname, color = hits),
                  max.overlaps = getOption("ggrepel.max.overlaps", default = 200))+
# scale_y_continuous(breaks = seq(0,1,0.2)) +
scale_x_continuous(limits = c(0, 1)) # length(unique(dat_s$drugname))
p_s_rank

# dat_s[dat_s$drugname == 'Birinapant',]
# dat_r
dat_r <- D_Enrich_ITS_r_score
dat_r$hits <- 'no'
dat_r[dat_r$ifFDA == 'FDA' & dat_r$weighted_rankTile < 0.01, ]$hits <- 'FDA'
dat_r[dat_r$ifFDA == 'FDAcancer' & dat_r$weighted_rankTile < 0.01, ]$hits <- 'FDAcancer'
# dat_r[which(dat_r$weighted_pvalue < 0.05 & dat_r$weighted_NES > 0),]$hits <- 'Target enriched'
dat_r = dat_r %>%
  arrange(desc(weighted_score)) %>%
  mutate("x" = row_number())
dat_r = dat_r[!is.na(dat_r$weighted_score),]
# dat_r_shown <- dat_r[grep("FDA", dat_r$hits), ]
dat_r$x_prop <- dat_r$x / max(dat_r$x)
dat_r_shown <- dat_r[-grep("no", dat_r$hits), ]
# dat_r_shown$hits <- factor(dat_r_shown$hits, levels = c('FDA','FDAcancer', 'Target enriched'))
dat_r_shown$hits <- factor(dat_r_shown$hits, levels = c('FDA','FDAcancer'))

p_r_rank <- ggplot(dat_r, aes(x = x_prop, y = weighted_score, color=hits))+
  geom_point(color="#338080", size = 0.5)+
  geom_vline(xintercept = 0, color = 'grey30')+
  geom_hline(yintercept = 0, color = 'grey30')+
  geom_vline(xintercept = 0.01, color = 'grey', linetype = 'dashed')+
  geom_vline(xintercept = 0.05, color = 'grey', linetype = 'dashed')+
  geom_vline(xintercept = 0.1, color = 'grey', linetype = 'dashed')+
  geom_rug(data = subset(dat_r, ifFDA == 'FDA'), color = "#ee766f", sides = "b",length = unit(3,"mm")) +
  geom_rug(data = subset(dat_r, ifFDA == 'FDAcancer'), color = "#2eb7be", sides = "b",length = unit(3,"mm")) +
  theme_bw()+
  labs(x = "Ranked drugs", y = "Potential")+
  geom_text_repel(inherit.aes = F, 
                  data = dat_r_shown, 
                  aes(x = x_prop, y = weighted_score, label = drugname, color = hits),
                  max.overlaps = getOption("ggrepel.max.overlaps", default = 200))+
# scale_y_continuous(breaks = seq(0,1,0.2))+
scale_x_continuous(limits = c(0, 1))
p_r_rank

# dat_r[dat_r$drugname == 'birinapant',]
ggsave("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/drug_merge3DB/D_Enrich_ITS_s_rankdot.pdf", p_s_rank, height = 6, width = 6)
ggsave("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/drug_merge3DB/D_Enrich_ITS_r_rankdot.pdf", p_r_rank, height = 6, width = 6)


# only plot for FDA approved drugs -----------------------------------------------------------
dat_s <- D_Enrich_ITS_s_score[D_Enrich_ITS_s_score$ifFDA != "no", ]

dat_s$hits <- 'no'
# dat_s[dat_s$ifFDA == 'FDA' & dat_s$weighted_rankTile < 0.01, ]$hits <- 'FDA'
dat_s[dat_s$ifFDA == 'FDAcancer' & dat_s$weighted_rankTile < 0.02, ]$hits <- 'FDAcancer'
# dat_s[which(dat_s$weighted_pvalue < 0.05 & dat_s$weighted_NES > 0),]$hits <- 'Target enriched'

dat_s = dat_s %>%
  arrange(desc(weighted_score)) %>%
  mutate("x" = row_number())
# dat_s = dat_s[!is.na(dat_s$weighted_score),]
dat_s$x_prop <- dat_s$x / max(dat_s$x)
dat_s_shown <- dat_s[-grep("no", dat_s$hits), ]
# dat_s_shown$hits <- factor(dat_s_shown$hits, levels = c('FDA','FDAcancer', 'Target enriched'))
dat_s_shown$hits <- factor(dat_s_shown$hits, levels = c('FDA','FDAcancer'))

# dat_s = dat_s[!is.na(dat_s$weighted_score),]
# dat_s_shown <- dat_s[grep("FDA", dat_s$hits), ]

p_s_rank <- ggplot(dat_s, aes(x = x_prop, y = weighted_score, color=hits))+
  geom_point(color="#B33333", size = 0.5)+
  geom_vline(xintercept = 0, color = 'grey30')+
  geom_hline(yintercept = 0, color = 'grey30')+
  geom_vline(xintercept = 0.02, color = 'grey', linetype = 'dashed')+
  # geom_vline(xintercept = 0.05, color = 'grey', linetype = 'dashed')+
  geom_vline(xintercept = 0.1, color = 'grey', linetype = 'dashed')+
  geom_rug(data = subset(dat_s, ifFDA == 'FDA'), color = "#ee766f", sides = "b",length = unit(3,"mm")) +
  geom_rug(data = subset(dat_s, ifFDA == 'FDAcancer'), color = "#2eb7be", sides = "b",length = unit(3,"mm")) +
  theme_bw()+
  labs(x = "Ranked drugs", y = "Potential")+
  geom_text_repel(inherit.aes = F, 
                  data = dat_s_shown, 
                  aes(x = x_prop, y = weighted_score, label = drugname), color = "#2eb7be", 
                  max.overlaps = getOption("ggrepel.max.overlaps", default = 200))+
# scale_y_continuous(breaks = seq(0,1,0.6)) +
scale_x_continuous(limits = c(0, 1))
p_s_rank


# dat_r
dat_r <- D_Enrich_ITS_r_score[D_Enrich_ITS_r_score$ifFDA != "no", ]

dat_r$hits <- 'no'
# dat_r[dat_r$ifFDA == 'FDA' & dat_r$weighted_rankTile < 0.01, ]$hits <- 'FDA'
dat_r[dat_r$ifFDA == 'FDAcancer' & dat_r$weighted_rankTile < 0.02, ]$hits <- 'FDAcancer'
# dat_r[which(dat_r$weighted_pvalue < 0.05 & dat_r$weighted_NES > 0),]$hits <- 'Target enriched'

dat_r = dat_r %>%
  arrange(desc(weighted_score)) %>%
  mutate("x" = row_number())
# dat_r = dat_r[!is.na(dat_r$weighted_score),]
# dat_r_shown <- dat_r[grep("FDA", dat_r$hits), ]
dat_r$x_prop <- dat_r$x / max(dat_r$x)
dat_r_shown <- dat_r[-grep("no", dat_r$hits), ]
# dat_r_shown$hits <- factor(dat_r_shown$hits, levels = c('FDA','FDAcancer', 'Target enriched'))
dat_r_shown$hits <- factor(dat_r_shown$hits, levels = c('FDA','FDAcancer'))


p_r_rank <- ggplot(dat_r, aes(x = x_prop, y = weighted_score, color=hits))+
  geom_point(color="#338080", size = 0.5)+
  geom_vline(xintercept = 0, color = 'grey30')+
  geom_hline(yintercept = 0, color = 'grey30')+
  geom_vline(xintercept = 0.02, color = 'grey', linetype = 'dashed')+
  # geom_vline(xintercept = 0.05, color = 'grey', linetype = 'dashed')+
  geom_vline(xintercept = 0.1, color = 'grey', linetype = 'dashed')+
  geom_rug(data = subset(dat_r, ifFDA == 'FDA'), color = "#ee766f", sides = "b", length = unit(3,"mm")) +
  geom_rug(data = subset(dat_r, ifFDA == 'FDAcancer'), color = "#2eb7be", sides = "b", length = unit(3,"mm")) +
  theme_bw()+
  labs(x = "Ranked drugs", y = "Potential")+
  geom_text_repel(inherit.aes = F, 
                  data = dat_r_shown, 
                  aes(x = x_prop, y = weighted_score, label = drugname),color = "#2eb7be", 
                  max.overlaps = getOption("ggrepel.max.overlaps", default = 200))+
# scale_y_continuous(breaks = seq(0,1, 0.2))+
scale_x_continuous(limits = c(0, 1))
p_r_rank

p_s_rank + p_r_rank
ggsave("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/drug_merge3DB/Dfda_Enrich_ITS_s_rankdot.pdf", p_s_rank, height = 6, width = 4)
ggsave("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/drug_merge3DB/Dfda_Enrich_ITS_r_rankdot.pdf", p_r_rank, height = 6, width = 4)
ggsave("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/drug_merge3DB/Dfda_Enrich_ITS_rankdot.pdf", p_s_rank+p_r_rank, height = 6, width = 8)

write.table(dat_s, "07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/drug_merge3DB/Dfda_Enrich_ITS_s_score.txt", sep = '\t', quote = F, row.names = F)
write.table(dat_r, "07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/drug_merge3DB/Dfda_Enrich_ITS_r_score.txt", sep = '\t', quote = F, row.names = F)


dat_s$rank2 = dat_s$x/nrow(dat_s)
nrow(dat_s[dat_s$rank2 < 0.1 & dat_s$ifFDA == "FDAcancer", ])
nrow(dat_s[dat_s$rank2 < 0.1 & dat_s$ifFDA == "FDAcancer", ])/nrow(dat_s[dat_s$ifFDA == "FDAcancer", ])

nrow(dat_s[dat_s$rank2 < 0.05 & dat_s$ifFDA == "FDAcancer", ])
nrow(dat_s[dat_s$rank2 < 0.05 & dat_s$ifFDA == "FDAcancer", ])/nrow(dat_s[dat_s$ifFDA == "FDAcancer", ])


dat_r$rank2 = dat_r$x/nrow(dat_r)
nrow(dat_r[dat_r$rank2 < 0.1 & dat_r$ifFDA == "FDAcancer", ])
nrow(dat_r[dat_r$rank2 < 0.1 & dat_r$ifFDA == "FDAcancer", ])/nrow(dat_r[dat_r$ifFDA == "FDAcancer", ])

nrow(dat_r[dat_r$rank2 < 0.05 & dat_r$ifFDA == "FDAcancer", ])
nrow(dat_r[dat_r$rank2 < 0.05 & dat_r$ifFDA == "FDAcancer", ])/nrow(dat_r[dat_r$ifFDA == "FDAcancer", ])

dat_s[dat_s$rank2 < 0.02 & dat_s$ifFDA == "FDAcancer", ]$drugname
dat_r[dat_r$rank2 < 0.05 & dat_r$ifFDA == "FDAcancer", ]$drugname
intersect(dat_s[dat_s$rank2 < 0.1 & dat_s$ifFDA == "FDAcancer", ]$drugname, dat_r[dat_r$rank2 < 0.1 & dat_r$ifFDA == "FDAcancer", ]$drugname)
intersect(dat_s[dat_s$rank2 < 0.05 & dat_s$ifFDA == "FDAcancer", ]$drugname, dat_r[dat_r$rank2 < 0.05 & dat_r$ifFDA == "FDAcancer", ]$drugname)

tmp_s = DTall[is.element(DTall$drugname_lower, 
                         dat_s[dat_s$rank2 < 0.05 & dat_s$ifFDA == "FDAcancer", ]$drugname),]
tmp_r = DTall[is.element(DTall$drugname_lower, 
                        dat_r[dat_r$rank2 < 0.05 & dat_r$ifFDA == "FDAcancer", ]$drugname),]
tmp = DTall[is.element(DTall$drugname_lower, 
intersect(dat_s[dat_s$rank2 < 0.05 & dat_s$ifFDA == "FDAcancer", ]$drugname, dat_r[dat_r$rank2 < 0.05 & dat_r$ifFDA == "FDAcancer", ]$drugname)),]

sort(table(tmp_s$genename))
sort(table(tmp_r$genename))
sort(table(tmp$genename))

# rm(list=ls())
work_path <- "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/"

setwd(work_path)
options(stringsAsFactors = F)
library(ggplot2)
library(aplot)
library(tidyr)
library(colorspace)
library(RColorBrewer)


dir.create("07_plotting_v2/05_drugITSprofiler_CombatRef_metaPlusFilter/sharedDrug_metaPlusFilter_alldrugFDA_all/")


drug_proof <- read.table("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/data/FDAapproved_p3Trial_IOtherapy/drug_list_proof_cancertype.txt", header = T, sep = '\t')
druglist <- unique(drug_proof[drug_proof$if_with_ICB == 'yes', ])
druglist <- druglist[!is.na(druglist$drug_index), ]
druglist <- druglist[druglist$proof != 'animal',]
drugAnnot_proof <- unique(druglist[c('pert_iname', 'proof')])
drugAnnot_cancer <- unique(druglist[c('pert_iname', 'Cancer_TCGA')])

    

# inst_GSE70138_is_gold <- read.csv("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/06_results_summaryandplotting/results_meta_oe_original/data/LINCS_is_gold/inst_GSE70138_is_gold.csv", header = T, sep = '\t')
filepath <- "06_results_summaryandplotting/results_meta_filter_ITSssgsea_CombatRef_17datasets/TUMERIC_m2_meta_en_ACATP0.05/70138/allgenes/"
resfilepath <- list.files(filepath)

filenames <- list.files(paste0(filepath, "SKCM_Metastatic/drugITSRank/"))
filenames <- filenames[grep('.txt', filenames)]


drugMoA <- read.csv('/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/s3_lincs_drug_MoA/lincs_portal2.0/drugMoA_lincsPortal2_2.csv',header = T)
drugTarget <- read.csv('/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/s3_lincs_drug_MoA/lincs_portal2.0/drugtarget_lincsPortal2.csv', header = T)
drug_maxfda <- read.csv('/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/s3_lincs_drug_MoA/lincs_portal2.0/drug_max_fda_phase_lincsPortal2.csv', header = T)
drug_maxfda34 <- unique(drug_maxfda[is.element(drug_maxfda$max_fda_phase,c(3,4)), ])
drug_maxfda34$drugname_lower <- tolower(drug_maxfda34$pert_iname)
load( "07_plotting_v2/data/drugTarget_merged3DB.Rdata")
drug_FDA <- unique(drugTarget_FDA[c(1,5,6)])
drug_FDA$drug_index = NA
drug_FDA <- unique(drug_FDA[c(1,4,3,2)])
names(drug_FDA) <- names(drug_maxfda34)
drug_FDA$max_fda_phase = 'approved'
drug_maxfda34 <- rbind(drug_maxfda34, drug_FDA)

for(f in filenames){
    # f = filenames[1]

    res <- read.csv(paste0("07_plotting_v2/05_drugITSprofiler_CombatRef_metaPlusFilter/drugITSres_metaPlusFilter/meta_en_ACATP0.05_",f), sep = '\t', header = T)
    res$drugname_lower <- tolower(res$pert_iname)
    res <- res[is.element(res$drugname_lower, unique(drug_maxfda34$drugname_lower)), ]
    res$cancer_dataset = paste0(res$cancer, "_", res$dataset)
    table(res$cancer_dataset)
    res2 <- lapply(as.list(unique(res$cancer_dataset)), function(x){
        tmp = res[is.element(res$cancer_dataset, x), ]
        tmp <- tmp[order(tmp$score, decreasing = T),]
        tmp$rank <- 1:nrow(tmp)
        tmp$rank_tile <- tmp$rank/nrow(tmp)
        return(tmp)
    })
    res = do.call(rbind, res2)
    

    table(res$cancer_dataset)
    
    res70138 <- res[res$dataset == '70138', ]
    res70138_stat <- data.frame(table(res70138$drug_index))
    dim(res70138_stat[res70138_stat$Freq >= length(unique(res70138$cancer_dataset)), ])
    res70138_stat <- data.frame(table(res70138$drug_index))
    dim(res70138_stat[res70138_stat$Freq >= length(unique(res70138$cancer_dataset)), ])
    dlist70138 = res70138_stat[res70138_stat$Freq == length(unique(res70138$cancer_dataset)), ]$Var1
    res70138_shared = res70138[is.element(res70138$drug_index, dlist70138), ]
    porder70138 = do.call(rbind, lapply(as.list(unique(res70138_shared$pert_iname)), function(d){
        tmp = res70138_shared[res70138_shared$pert_iname == d, ]
        return(data.frame(pert_iname = d, num = nrow(tmp[tmp$rank_tile < 0.1, ])))
    }))
    porder70138 = porder70138[order(porder70138$num, decreasing = T),]
    # porder70138 = porder70138[porder70138$num != 0, ]
    res70138_shared = res70138_shared[is.element(res70138_shared$pert_iname, porder70138$pert_iname), ]
    res70138_shared = res70138_shared[!is.na(res70138_shared$pert_iname), ]


    dAproof70138 <- merge(res70138_shared, drugAnnot_proof, by = 'pert_iname', all.x = T)
    dAproof70138 <- as.data.frame(table(unique(dAproof70138[c('pert_iname', 'proof')])))
    p_dAproof70138 <- ggplot(dAproof70138, aes(x = proof, y = pert_iname, color=Freq))+
                        geom_point(size = 5, shape = 15) +
                        scale_color_gradient(high = "#ff79136b", low = '#ffffff00') + 
                        scale_y_discrete(position="right") +
                        theme_minimal() + xlab(NULL) + ylab(NULL) +
                        theme(axis.text.x = element_text(size = 10,angle=90, hjust=1),
                              axis.text.y = element_blank()) + 
                        labs(fill = "proof")

    dAcancer70138 <- merge(res70138_shared, drugAnnot_cancer, by = 'pert_iname', all.x = T)
    dAcancer70138 <- as.data.frame(table(unique(dAcancer70138[c('pert_iname', 'Cancer_TCGA')])))
    p_dAcancer70138 <- ggplot(dAcancer70138, aes(x = Cancer_TCGA, y = pert_iname,color=Freq))+
                        geom_point(size = 5, shape = 15) +
                        scale_color_gradient(high = "#ff79136b", low = '#ffffff00') + 
                        scale_y_discrete(position="right") +
                        theme_minimal() + xlab(NULL) + ylab(NULL) +
                        theme(axis.text.x = element_text(size = 10,angle=90, hjust=1),
                              axis.text.y = element_blank()) + 
                        labs(fill = "Cancer_TCGA")
    # p_dAcancer70138

    res70138_shared$pert_iname = factor(res70138_shared$pert_iname, levels = unique(porder70138$pert_iname ))
    res70138_top10 = res70138_shared[res70138_shared$rank_tile < 0.1, ]

    p_70138 = ggplot(data=res70138_shared, aes(x=cancer_dataset, y = pert_iname))+
            geom_point(aes(color = cut(rank_tile, 10))) +
            scale_color_manual(values = brewer.pal(10, "RdYlBu")) +
            # scale_color_binned_sequential("Terrain 2", rev = F, begin = 0, end = 1) +
            # scale_color_gradient2(high = "#41B0C3",mid = 'white', low = "#ed3926", midpoint = 0.5,
            #                     breaks=seq(floor(min(res70138_shared$rank_tile)),ceiling(max(res70138_shared$rank_tile)),0.1))+
            scale_size_continuous(range = c(1,10))+
            guides(size=F)+
            theme_bw()+ xlab(NULL) + ylab(NULL) +
            geom_text(data=res70138_top10,
                      aes(x=cancer_dataset, y = pert_iname, label=paste0(round(rank_tile*100, 2), "%"))) +
            theme(legend.key.height = unit(3,'cm'),
                    legend.justification = c(0,0),
                    legend.title = element_blank(),
                    axis.text.x=element_text(angle=90, hjust=1),
                    text = element_text(size=15))
    # p_70138

    p_70138 %>% 
    insert_left(p_dAproof70138, width = .3) %>% 
    insert_right(p_dAcancer70138, width = .5) %>%
    ggsave(filename = paste0("07_plotting_v2/05_drugITSprofiler_CombatRef_metaPlusFilter/sharedDrug_metaPlusFilter_alldrugFDA_all/", gsub("_res.txt", "",f), "_70138.pdf"), width = 15, height = 15)





    # res <- read.csv(paste0("07_plotting_v2/05_drugITSprofiler_CombatRef_metaPlusFilter/drugITSres_metaPlusFilter/meta_en_ACATP0.05_",f), sep = '\t', header = T)
    # res$cancer_dataset = paste0(res$cancer, "_", res$dataset)
    res92742 <- res[res$dataset == '92742', ]
    res92742_stat <- data.frame(table(res92742$drug_index))
    dim(res92742_stat[res92742_stat$Freq >= length(unique(res92742$cancer_dataset)), ])
    res92742_stat <- data.frame(table(res92742$drug_index))
    dim(res92742_stat[res92742_stat$Freq >= length(unique(res92742$cancer_dataset)), ])
    dlist92742 = res92742_stat[res92742_stat$Freq == length(unique(res92742$cancer_dataset)), ]$Var1
    res92742_shared = res92742[is.element(res92742$drug_index, dlist92742), ]
    porder92742 = do.call(rbind, lapply(as.list(unique(res92742_shared$pert_iname)), function(d){
        tmp = res92742_shared[res92742_shared$pert_iname == d, ]
        return(data.frame(pert_iname = d, num = nrow(tmp[tmp$rank_tile < 0.1, ])))
    }))
    porder92742 = porder92742[order(porder92742$num, decreasing = T),]
    # porder92742 = porder92742[porder92742$num != 0, ]
    res92742_shared = res92742_shared[is.element(res92742_shared$pert_iname, porder92742$pert_iname), ]
    res92742_shared = res92742_shared[!is.na(res92742_shared$pert_iname), ]



    dAproof92742 <- merge(res92742_shared, drugAnnot_proof, by = 'pert_iname', all.x = T)
    dAproof92742 <- as.data.frame(table(unique(dAproof92742[c('pert_iname', 'proof')])))
    p_dAproof92742 <- ggplot(dAproof92742, aes(x = proof, y = pert_iname, color=Freq))+
                        geom_point(size = 5, shape = 15) +
                        scale_color_gradient(high = "#ff79136b", low = '#ffffff00') + 
                        scale_y_discrete(position="right") +
                        theme_minimal() + xlab(NULL) + ylab(NULL) +
                        theme(# panel.grid = element_blank(),
                              # panel.border = element_rect(color="grey"),
                              # axis.ticks = element_blank(),
                              axis.text.x = element_text(size = 10,angle=90, hjust=1),
                              axis.text.y = element_blank()) + 
                        labs(fill = "proof")

    dAcancer92742 <- merge(res92742_shared, drugAnnot_cancer, by = 'pert_iname', all.x = T)
    dAcancer92742 <- as.data.frame(table(unique(dAcancer92742[c('pert_iname', 'Cancer_TCGA')])))
    p_dAcancer92742 <- ggplot(dAcancer92742, aes(x = Cancer_TCGA, y = pert_iname,color=Freq))+
                        geom_point(size = 5, shape = 15) +
                        scale_color_gradient(high = "#ff79136b", low = '#ffffff00') + 
                        scale_y_discrete(position="right") +
                        theme_minimal() + xlab(NULL) + ylab(NULL) +
                        theme(#panel.grid = element_blank(),
                              #panel.border = element_rect(color="grey"),
                              #axis.ticks = element_blank(),
                              axis.text.x = element_text(size = 10,angle=90, hjust=1),
                              axis.text.y = element_blank()) + 
                        labs(fill = "Cancer_TCGA")
    # p_dAcancer92742
    res92742_shared$pert_iname = factor(res92742_shared$pert_iname, levels = porder92742$pert_iname )
    res92742_top10 = res92742_shared[res92742_shared$rank_tile < 0.1, ]

      p_92742 = ggplot(data=res92742_shared, aes(x=cancer_dataset, y = pert_iname))+
            geom_point(aes(color = cut(rank_tile, 10))) +
            scale_color_manual(values = brewer.pal(10, "RdYlBu")) +
            # scale_color_binned_sequential("Terrain 2", rev = F, begin = 0, end = 1) +
            # scale_color_gradient2(high = "#41B0C3",mid = 'white', low = "#ed3926", midpoint = 0.5,
            #                     breaks=seq(floor(min(res70138_shared$rank_tile)),ceiling(max(res70138_shared$rank_tile)),0.1))+
            scale_size_continuous(range = c(1,10))+
            guides(size=F)+
            theme_bw()+ xlab(NULL) + ylab(NULL) +
            geom_text(data=res92742_top10,
                      aes(x=cancer_dataset, y = pert_iname, label=paste0(round(rank_tile*100, 2), "%"))) +
            theme(legend.key.height = unit(3,'cm'),
                    legend.justification = c(0,0),
                    legend.title = element_blank(),
                    axis.text.x=element_text(angle=90, hjust=1),
                    text = element_text(size=15))
    # p_92742

    p_92742 %>% 
    insert_left(p_dAproof92742, width = .3) %>% 
    insert_right(p_dAcancer92742, width = .5) %>%
    ggsave(filename = paste0("07_plotting_v2/05_drugITSprofiler_CombatRef_metaPlusFilter/sharedDrug_metaPlusFilter_alldrugFDA_all/", gsub("_res.txt", "",f), "_92742.pdf"), width = 20, height = 30)



    # res <- read.csv(paste0("07_plotting_v2/05_drugITSprofiler_CombatRef_metaPlusFilter/drugITSres_metaPlusFilter/meta_en_ACATP0.05_",f), sep = '\t', header = T)
    # res$cancer_dataset = paste0(res$cancer, "_", res$dataset)
    # table(res$cancer_dataset)

    res_stat <- data.frame(table(res$drug_index))
    dim(res_stat[res_stat$Freq >= length(unique(res$cancer_dataset)), ])
    
    # dlist = res_stat[res_stat$Freq >= 15, ]$Var1
    dlist = union(dlist70138, dlist92742)
    res_shared = res[is.element(res$drug_index, dlist), ]
    porder = do.call(rbind, lapply(as.list(unique(res_shared$pert_iname)), function(d){
        tmp = res_shared[res_shared$pert_iname == d, ]
        return(data.frame(pert_iname = d, num = nrow(tmp[tmp$rank_tile < 0.1, ])))
    }))
    porder = porder[order(porder$num, decreasing = T),]
    # porder = porder[porder$num != 0, ]
    res_shared = res_shared[is.element(res_shared$pert_iname, porder$pert_iname), ]
    res_shared = res_shared[!is.na(res_shared$pert_iname), ]

    dAproof <- merge(res_shared, drugAnnot_proof, by = 'pert_iname', all.x = T)
    dAproof <- as.data.frame(table(unique(dAproof[c('pert_iname', 'proof')])))
    p_dAproof <- ggplot(dAproof, aes(x = proof, y = pert_iname, color=Freq))+
                        geom_point(size = 5, shape = 15) +
                        scale_color_gradient(high = "#ff79136b", low = '#ffffff00') + 
                        scale_y_discrete(position="right") +
                        theme_minimal() + xlab(NULL) + ylab(NULL) +
                        theme(# panel.grid = element_blank(),
                              # panel.border = element_rect(color="grey"),
                              # axis.ticks = element_blank(),
                              axis.text.x = element_text(size = 10,angle=90, hjust=1),
                              axis.text.y = element_blank()) + 
                        labs(fill = "proof")

    dAcancer <- merge(res_shared, drugAnnot_cancer, by = 'pert_iname', all.x = T)
    dAcancer <- as.data.frame(table(unique(dAcancer[c('pert_iname', 'Cancer_TCGA')])))
    p_dAcancer <- ggplot(dAcancer, aes(x = Cancer_TCGA, y = pert_iname,color=Freq))+
                        geom_point(size = 5, shape = 15) +
                        scale_color_gradient(high = "#ff79136b", low = '#ffffff00') + 
                        scale_y_discrete(position="right") +
                        theme_minimal() + xlab(NULL) + ylab(NULL) +
                        theme(#panel.grid = element_blank(),
                              #panel.border = element_rect(color="grey"),
                              #axis.ticks = element_blank(),
                              axis.text.x = element_text(size = 10,angle=90, hjust=1),
                              axis.text.y = element_blank()) + 
                        labs(fill = "Cancer_TCGA")
    # p_dAcancer
    res_shared$pert_iname = factor(res_shared$pert_iname, levels = porder$pert_iname )
    res_top10 = res_shared[res_shared$rank_tile < 0.1, ]
    p = ggplot(data=res_shared, aes(x=cancer_dataset, y = pert_iname))+
            geom_point(aes(color = cut(rank_tile, 10))) +
            scale_color_manual(values = brewer.pal(10, "RdYlBu")) +
            # scale_color_binned_sequential("Terrain 2", rev = F, begin = 0, end = 1) +
            # scale_color_gradient2(high = "#41B0C3",mid = 'white', low = "#ed3926", midpoint = 0.5,
            #                     breaks=seq(floor(min(res70138_shared$rank_tile)),ceiling(max(res70138_shared$rank_tile)),0.1))+
            scale_size_continuous(range = c(1,10))+
            guides(size=F)+
            theme_bw()+ xlab(NULL) + ylab(NULL) +
            geom_text(data=res_top10,
                      aes(x=cancer_dataset, y = pert_iname, label=paste0(round(rank_tile*100, 2), "%"))) +
            theme(legend.key.height = unit(3,'cm'),
                    legend.justification = c(0,0),
                    legend.title = element_blank(),
                    axis.text.x=element_text(angle=90, hjust=1),
                    text = element_text(size=15))
    # p



    p %>% 
    insert_left(p_dAproof, width = .3) %>% 
    insert_right(p_dAcancer, width = .5) %>%
    ggsave(filename = paste0("07_plotting_v2/05_drugITSprofiler_CombatRef_metaPlusFilter/sharedDrug_metaPlusFilter_alldrugFDA_all/", gsub("_res.txt", "",f), ".pdf"), width = 20, height = 30)


}




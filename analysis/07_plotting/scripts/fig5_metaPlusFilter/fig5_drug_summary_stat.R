# rm(list=ls())
work_path <- "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/"
setwd(work_path)
options(stringsAsFactors = F)
library(dplyr)
library(ggrepel)
library(pheatmap)
library(reshape2)
dir.create("07_plotting_v2/05_drugITSprofiler_CombatRef_metaPlusFilter/drugITSres_en_dotplot/")


# ------------------------------------------------------------------------------
# run for dataset GSE70138

dataset = '70138'
filepath <- "06_results_summaryandplotting/results_meta_filter_ITSssgsea_CombatRef_17datasets/TUMERIC_m2_meta_en_ACATP0.05/70138/allgenes/"
cancerlist <- list.files(filepath)

filenames <- list.files(paste0(filepath, "SKCM_Metastatic/drugITSRank/"))
filenames <- filenames[grep('.txt', filenames)]
if(length(grep("unweight",  filenames))>0)filenames <- filenames[-grep('unweight',filenames)]
f = "model_weight_weighted_glmnet_weight_res.txt"

dat_70138 <- lapply(as.list(cancerlist), function(cancer){
    
    # cancer = 'SKCM_Metastatic' 
    
    cancer_type = strsplit(cancer,'_')[[1]][1]
    cancer_subtype = strsplit(cancer,'_')[[1]][2]    
    drugrankFile = paste0("06_results_summaryandplotting/results_meta_filter_ITSssgsea_CombatRef_17datasets/TUMERIC_m2_meta_en_ACATP0.05/",dataset,"/allgenes/",cancer,"/drugITSRank/",f)
    drugrank <- read.csv(drugrankFile, sep = '\t', header = T)
    drugrank$rank_tile = drugrank$rank / nrow(drugrank)
    drugrank$cancer = cancer
    return(drugrank)

})
dat_70138 <- do.call(rbind, dat_70138)


cancer = 'SKCM_Primary' 
cancer_type = strsplit(cancer,'_')[[1]][1]
cancer_subtype = strsplit(cancer,'_')[[1]][2]
drugrankFile = paste0("06_results_summaryandplotting/results_meta_filter_ITSssgsea_CombatRef_17datasets/CPE_m2_meta_en_ACATP0.05/",dataset,"/allgenes/",cancer,"/drugITSRank/", f)
drugrank <- read.csv(drugrankFile, sep = '\t', header = T)
drugrank$rank_tile = drugrank$rank / nrow(drugrank)
drugrank$cancer = cancer

dat_70138 <- rbind(dat_70138, drugrank)
dat_70138$dataset = '70138'
# names(dat_70138) = c(cancerlist, 'SKCM_Primary')



# ------------------------------------------------------------------------------
# run for dataset GSE92742

dataset = '92742'
filepath <- "06_results_summaryandplotting/results_meta_filter_ITSssgsea_CombatRef_17datasets/TUMERIC_m2_meta_en_ACATP0.05/92742/allgenes/"
cancerlist <- list.files(filepath)

filenames <- list.files(paste0(filepath, "SKCM_Metastatic/drugITSRank/"))
filenames <- filenames[grep('.txt', filenames)]
if(length(grep("unweight",  filenames))>0)filenames <- filenames[-grep('unweight',filenames)]

dat_92742 <- lapply(as.list(cancerlist), function(cancer){
    
    cancer_type = strsplit(cancer,'_')[[1]][1]
    cancer_subtype = strsplit(cancer,'_')[[1]][2]
    drugrankFile = paste0("06_results_summaryandplotting/results_meta_filter_ITSssgsea_CombatRef_17datasets/TUMERIC_m2_meta_en_ACATP0.05/",dataset,"/allgenes/",cancer,"/drugITSRank/",f)
    drugrank <- read.csv(drugrankFile, sep = '\t', header = T)
    drugrank$rank_tile = drugrank$rank / nrow(drugrank)
    drugrank$cancer = cancer
    return(drugrank)

})

dat_92742 <- do.call(rbind, dat_92742)


cancer = 'SKCM_Primary' 
cancer_type = strsplit(cancer,'_')[[1]][1]
cancer_subtype = strsplit(cancer,'_')[[1]][2]
  
drugrankFile = paste0("06_results_summaryandplotting/results_meta_filter_ITSssgsea_CombatRef_17datasets/CPE_m2_meta_en_ACATP0.05/",dataset,"/allgenes/",cancer,"/drugITSRank/", f)
drugrank <- read.csv(drugrankFile, sep = '\t', header = T)
drugrank$rank_tile = drugrank$rank / nrow(drugrank)
drugrank$cancer = cancer

dat_92742 <- rbind(dat_92742, drugrank)
dat_92742$dataset = dataset

dat_all = rbind(dat_70138, dat_92742)
dat_all$cancer_dataset = paste0(dat_all$cancer, "_", dat_all$dataset)

save(dat_70138, dat_92742, dat_all, file = paste0("06_results_summaryandplotting/results_meta_filter_ITSssgsea_CombatRef_17datasets/",f,".Rdata"))

# ------------------------------------------------------------------------------
# overall plot
# ------------------------------------------------------------------------------

drug_proof <- read.table("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/data/FDAapproved_p3Trial_IOtherapy/drug_list_proof_cancertype.txt", header = T, sep = '\t')
druglist <- unique(drug_proof[drug_proof$if_with_ICB == 'yes', ])
druglist <- druglist[!is.na(druglist$drug_index), ]
druglist <- druglist[druglist$proof != 'animal',]


drugMoA <- read.csv('/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/s3_lincs_drug_MoA/lincs_portal2.0/drugMoA_merged_2.txt',sep='\t', header = T)

drugTarget <- read.csv('/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/s3_lincs_drug_MoA/lincs_portal2.0/drugtarget_lincsPortal2.csv', header = T)
drug_maxfda <- read.csv('/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/s3_lincs_drug_MoA/lincs_portal2.0/drug_max_fda_phase_lincsPortal2.csv', header = T)
drug_maxfda34 <- unique(drug_maxfda[is.element(drug_maxfda$max_fda_phase, c(3,4)), ])
drug_maxfda34$drugname_lower <- tolower(drug_maxfda34$pert_iname)


load( "07_plotting_v2/data/drugTarget_merged3DB.Rdata")
drug_FDA <- unique(drugTarget_FDA[c(1,5,6)])
drug_FDA$drug_index = NA
drug_FDA <- unique(drug_FDA[c(1,4,3,2)])
names(drug_FDA) <- names(drug_maxfda34)
drug_FDA$max_fda_phase = 'approved'
drug_maxfda34 <- rbind(drug_maxfda34, drug_FDA)

drug_FDAcancer <- unique(drugTarget_FDA[which(drugTarget_FDA$ifCancer_drugbank == 'yes'),])



load(file = paste0("06_results_summaryandplotting/results_meta_filter_ITSssgsea_CombatRef_17datasets/",f,".Rdata"))
# dat_70138, dat_92742, dat_all
# ------------------------------------------------------------------------------

plot_savepath = '/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/07_plotting_v2/05_drugITSprofiler_CombatRef_metaPlusFilter/'
dir.create(paste0(plot_savepath, 'drug_overall_stat'))

save(dat_70138, dat_92742, dat_all, file = paste0(plot_savepath, "drug_overall_stat/glmnet_weight_res_meta.Rdata"))

drugproof_annot = unique(dat_all[c(1,2)])
drugproof_annot$proof = 'no'
drugproof_annot[is.element(drugproof_annot$pert_iname, unique(druglist$pert_iname)), ]$proof = 'yes'

drugMoA2 = unique(drugMoA[c('drug_index', "MoA")])
names(drugMoA2) = c('drug_index', "MoA_category")
drugproof_annot = merge(drugproof_annot, drugMoA2, by = 'drug_index', all=T)
drugproof_annot = drugproof_annot[!is.na(drugproof_annot$pert_iname),]
rownames(drugproof_annot) = drugproof_annot$pert_iname
mat_70138_2 <- dcast(dat_70138, drug_index ~ cancer , value.var = 'score' )
rownames(mat_70138_2)= mat_70138_2$drug_index
mat_70138_2 <- merge(drugproof_annot, mat_70138_2, by = 'drug_index', all.y=T)
write.table(mat_70138_2, paste0(plot_savepath, 'drug_overall_stat/drugAll_70138_annot.txt'), sep='\t', quote=F, row.names = F)

mat_70138 <- dcast(dat_70138, pert_iname ~ cancer , value.var = 'score' )
rownames(mat_70138)= mat_70138$pert_iname
mat_70138[is.na(mat_70138)] = -666
pdf(paste0(plot_savepath, 'drug_overall_stat/drugAll_70138_heatmap.pdf' ), width = 10, height = 20)
rangevalue = max(abs(min(dat_70138$score)), max(dat_70138$score))
bk <- c(seq(-rangevalue,-0.5,by=0.5),seq(0,rangevalue, by=0.5))

pheatmap(mat_70138[-1],
         scale = "none", 
         cluster_row = T, cluster_col = T,
         annotation_row = drugproof_annot[c('proof', 'MoA_category')],
         color = c('grey', colorRampPalette(colors = c("blue","white"))(length(bk)/2),colorRampPalette(colors = c("white","red"))(length(bk)/2)),
         legend_breaks=c(-666, seq(-rangevalue, rangevalue,2)),
         breaks=bk)
dev.off()



dat_70138_fda <- dat_70138[is.element(dat_70138$pert_iname, drug_maxfda34$drugname_lower), ]
mat_70138_fda  <- dcast(dat_70138_fda, pert_iname ~ cancer , value.var = 'score' )
rownames(mat_70138_fda)= mat_70138_fda$pert_iname
mat_70138_fda[is.na(mat_70138_fda)] = -666
pdf(paste0(plot_savepath, 'drug_overall_stat/drugFDA_70138_heatmap.pdf'), width = 10, height = 20)

rangevalue = max(abs(min(dat_70138_fda$score)), max(dat_70138_fda$score))
bk <- c(seq(-rangevalue,-0.5,by=0.5),seq(0,rangevalue, by=0.5))
pheatmap(mat_70138_fda[-1],
         scale = "none", 
         cluster_row = T, cluster_col = T,
         annotation_row = drugproof_annot[c('proof', 'MoA_category')],
         color = c('grey', colorRampPalette(colors = c("blue","white"))(length(bk)/2),colorRampPalette(colors = c("white","red"))(length(bk)/2)),
         legend_breaks=c(-666, seq(-rangevalue, rangevalue,2)),
         breaks=bk)
dev.off()

dat_70138_fdacancer <- dat_70138[is.element(dat_70138$pert_iname, drug_FDAcancer$drugname_lower), ]
mat_70138_fdacancer  <- dcast(dat_70138_fdacancer, pert_iname ~ cancer , value.var = 'score' )
rownames(mat_70138_fdacancer)= mat_70138_fdacancer$pert_iname
mat_70138_fdacancer[is.na(mat_70138_fdacancer)] = -666
pdf(paste0(plot_savepath, 'drug_overall_stat/drugFDAcancer_70138_heatmap.pdf'), width = 10, height = 20)
rangevalue = max(abs(min(dat_70138_fdacancer$score)), max(dat_70138_fdacancer$score))
bk <- c(seq(-rangevalue,-0.5,by=0.5),seq(0,rangevalue, by=0.5))
pheatmap(mat_70138_fdacancer[-1],
         scale = "none", 
         cluster_row = T, cluster_col = T,
         annotation_row = drugproof_annot[c('proof', 'MoA_category')],
         color = c('grey', colorRampPalette(colors = c("blue","white"))(length(bk)/2),colorRampPalette(colors = c("white","red"))(length(bk)/2)),
         legend_breaks=c(-666, seq(-rangevalue, rangevalue,2)),
         breaks=bk)
dev.off()

pheatmap(mat_70138_fdacancer[-c(1,4,5)],
         scale = "none", 
         cluster_row = T, cluster_col = T,
         annotation_row = drugproof_annot[c('proof', 'MoA_category')],
         color = c('grey', colorRampPalette(colors = c("blue","white"))(length(bk)/2),colorRampPalette(colors = c("white","red"))(length(bk)/2)),
         legend_breaks=c(-666, seq(-rangevalue, rangevalue,2)),
         breaks=bk)
dev.off()

# ------------------------------------------------------------------------------
mat_92742_2 <- dcast(dat_92742, drug_index ~ cancer , value.var = 'score' )
rownames(mat_92742_2)= mat_92742_2$drug_index
mat_92742_2 <- merge(drugproof_annot, mat_92742_2, by = 'drug_index', all.y=T)
write.table(mat_92742_2, paste0(plot_savepath, 'drug_overall_stat/drugAll_92742_annot.txt'),sep='\t', quote=F, row.names = F)

mat_92742 <- dcast(dat_92742, pert_iname ~ cancer , value.var = 'score' )
rownames(mat_92742)= mat_92742$pert_iname
mat_92742[is.na(mat_92742)] = -666
pdf(paste0(plot_savepath, 'drug_overall_stat/drugAll_92742_heatmap.pdf' ), width = 10, height = 20)

rangevalue = max(abs(min(dat_92742$score)), max(dat_92742$score))
bk <- c(seq(-rangevalue,-0.5,by=0.5),seq(0,rangevalue, by=0.5))
pheatmap(mat_92742[-1],
         scale = "none", 
         cluster_row = T, cluster_col = T,
         annotation_row = drugproof_annot[c('proof', 'MoA_category')],
         color = c('grey', colorRampPalette(colors = c("blue","white"))(length(bk)/2),colorRampPalette(colors = c("white","red"))(length(bk)/2)),
         legend_breaks=c(-666, seq(-rangevalue, rangevalue,2)),
         breaks=bk)
dev.off()

dat_92742_fda <- dat_92742[is.element(dat_92742$pert_iname, drug_maxfda34$drugname_lower), ]
mat_92742_fda  <- dcast(dat_92742_fda, pert_iname ~ cancer , value.var = 'score' )
rownames(mat_92742_fda)= mat_92742_fda$pert_iname
mat_92742_fda[is.na(mat_92742_fda)] = -666
pdf(paste0(plot_savepath, 'drug_overall_stat/drugFDA_92742_heatmap.pdf' ), width = 10, height = 20)

rangevalue = max(abs(min(dat_92742_fda$score)), max(dat_92742_fda$score))
bk <- c(seq(-rangevalue,-0.5,by=0.5),seq(0,rangevalue, by=0.5))
pheatmap(mat_92742_fda[-1],
         scale = "none", 
         cluster_row = T, cluster_col = T,
         annotation_row = drugproof_annot[c('proof', 'MoA_category')],
         color = c('grey', colorRampPalette(colors = c("blue","white"))(length(bk)/2),colorRampPalette(colors = c("white","red"))(length(bk)/2)),
         legend_breaks=c(-666, seq(-rangevalue, rangevalue,2)),
         breaks=bk)
dev.off()


dat_92742_fdacancer <- dat_92742[is.element(dat_92742$pert_iname, drug_FDAcancer$drugname_lower), ]
mat_92742_fdacancer  <- dcast(dat_92742_fdacancer, pert_iname ~ cancer , value.var = 'score' )
rownames(mat_92742_fdacancer)= mat_92742_fdacancer$pert_iname
mat_92742_fdacancer[is.na(mat_92742_fdacancer)] = -666
pdf(paste0(plot_savepath, 'drug_overall_stat/drugFDAcancer_92742_heatmap.pdf'), width = 10, height = 20)
rangevalue = max(abs(min(dat_92742_fdacancer$score)), max(dat_92742_fdacancer$score))
bk <- c(seq(-rangevalue,-0.5,by=0.5),seq(0,rangevalue, by=0.5))
pheatmap(mat_92742_fdacancer[-1],
         scale = "none", 
         cluster_row = T, cluster_col = T,
         annotation_row = drugproof_annot[c('proof', 'MoA_category')],
         color = c('grey', colorRampPalette(colors = c("blue","white"))(length(bk)/2),colorRampPalette(colors = c("white","red"))(length(bk)/2)),
         legend_breaks=c(-666, seq(-rangevalue, rangevalue,2)),
         breaks=bk)
dev.off()

# ------------------------------------------------------------------------------
dat_all_2 <- dcast(dat_all, drug_index ~ cancer , value.var = 'score' )
rownames(dat_all_2)= dat_all_2$drug_index
dat_all_2 <- merge(drugproof_annot, dat_all_2, by = 'drug_index', all.y=T)
write.table(dat_all_2, paste0(plot_savepath, 'drug_overall_stat/drugAll_all_annot.txt'),sep='\t', quote=F, row.names = F)

mat_all <- dcast(dat_all, pert_iname ~ cancer_dataset , value.var = 'score' )
rownames(mat_all)= mat_all$pert_iname
mat_all[is.na(mat_all)] = -666
pdf(paste0(plot_savepath, 'drug_overall_stat/drugAll_all_heatmap.pdf' ), width = 15, height = 20)
rangevalue = max(abs(min(dat_all$score)), max(dat_all$score))
bk <- c(seq(-rangevalue,-0.5,by=0.5),seq(0,rangevalue, by=0.5))
pheatmap(mat_all[-1],
         scale = "none", 
         cluster_row = T, cluster_col = T,
         annotation_row = drugproof_annot[c('proof', 'MoA_category')],
         color = c('grey', colorRampPalette(colors = c("blue","white"))(length(bk)/2),colorRampPalette(colors = c("white","red"))(length(bk)/2)),
         legend_breaks=c(-666, seq(-rangevalue, rangevalue,2)),
         breaks=bk)
dev.off()


dat_all_fda <- dat_all[is.element(dat_all$pert_iname, drug_maxfda34$drugname_lower), ]
mat_all_fda  <- dcast(dat_all_fda, pert_iname ~ cancer_dataset , value.var = 'score' )
rownames(mat_all_fda)= mat_all_fda$pert_iname
mat_all_fda[is.na(mat_all_fda)] = -666
pdf(paste0(plot_savepath, 'drug_overall_stat/drugFDA_all_heatmap.pdf' ), width = 15, height = 20)
rangevalue = max(abs(min(dat_all_fda$score)), max(dat_all_fda$score))
bk <- c(seq(-rangevalue,-0.5,by=0.5),seq(0,rangevalue, by=0.5))
pheatmap(mat_all_fda[-1],
         scale = "none", 
         cluster_row = T, cluster_col = T,
         annotation_row = drugproof_annot[c('proof', 'MoA_category')],
         color = c('grey', colorRampPalette(colors = c("blue","white"))(length(bk)/2),colorRampPalette(colors = c("white","red"))(length(bk)/2)),
         legend_breaks=c(-666, seq(-rangevalue, rangevalue,2)),
         breaks=bk)
dev.off()


dat_all_fdacancer <- dat_all[is.element(dat_all$pert_iname, drug_FDAcancer$drugname_lower), ]
mat_all_fdacancer  <- dcast(dat_all_fdacancer, pert_iname ~ cancer_dataset , value.var = 'score' )
rownames(mat_all_fdacancer)= mat_all_fdacancer$pert_iname
mat_all_fdacancer[is.na(mat_all_fdacancer)] = -666

pdf(paste0(plot_savepath, 'drug_overall_stat/drugFDAcancer_all_heatmap.pdf'), width = 10, height = 20)
rangevalue = max(abs(min(dat_all_fdacancer$score)), max(dat_all_fdacancer$score))
bk <- c(seq(-rangevalue,-0.5,by=0.5),seq(0,rangevalue, by=0.5))
pheatmap(mat_all_fdacancer[-1],
         scale = "none", 
         cluster_row = T, cluster_col = T,
         annotation_row = drugproof_annot[c('proof', 'MoA_category')],
         color = c('grey', colorRampPalette(colors = c("blue","white"))(length(bk)/2),colorRampPalette(colors = c("white","red"))(length(bk)/2)),
         legend_breaks=c(-666, seq(-rangevalue, rangevalue,2)),
         breaks=bk)
dev.off()


tmp=drugproof_annot[is.element(tolower(drugproof_annot$pert_iname), mat_all_fdacancer$pert_iname), ]


# rm(list=ls())
work_path <- "D:/veronica_files/lab/projects/p2_1_immunotherapy_targetedtherapy/"
setwd(work_path)
options(stringsAsFactors = F)
library(reshape2)
library(pheatmap)

load("07_plotting_v2/05_drugITSprofiler_CombatRef_metaPlusFilter_corr2/drug_overall_stat/glmnet_weight_res_meta.Rdata")


mat_70138 <- dcast(dat_70138, pert_iname ~ cancer , value.var = 'score' )
rownames(mat_70138)= mat_70138$pert_iname
mat_70138[is.na(mat_70138)] = -666

cls <- (cutree(hclust(dist(mat_70138[-c(4,5)])), k=50))
cls = data.frame(pert_iname = names(cls), cls=cls)
cls=cls[order(cls$cls),]
cls$cls = paste0('cls', cls$cls)
mat_70138 = mat_70138[cls$pert_iname,]

rangevalue = max(abs(min(dat_70138$score)), max(dat_70138$score))
bk <- c(seq(-rangevalue,-0.5,by=0.5),seq(0,rangevalue, by=0.5))

pheatmap(t(mat_70138[-c(1)]),
         scale = "none", 
         cluster_row = T, cluster_col = T,
         border="white",
         cellwidth = 0.5,
         cellheight = 20,
         fontsize = 4,
         annotation_col = cls[2],# drugproof_annot['proof'],
         cutree_cols = 30, treeheight_col=30,
         # cutree_cols = 5,
         color = c('grey95', colorRampPalette(colors = c("purple","white"))(length(bk)/2),colorRampPalette(colors = c("white","red"))(length(bk)/2)),
         legend_breaks=c(-666, seq(-rangevalue, rangevalue,2)),
         breaks=bk)
dev.off()

drugmoa <- read.csv('07_plotting_v2/data/drugMoA_mergeNew.txt', sep='\t', header = T)
drugmoa <- unique(drugmoa[-c(1,3,5)])
mat_70138 <- dcast(dat_70138, pert_iname ~ cancer , value.var = 'score' )
rownames(mat_70138)= mat_70138$pert_iname
# mat_70138[is.na(mat_70138)] = -666

mat_70138$average = apply(mat_70138[-1], 1, function(x)mean(x[-which(is.na(x))]))
drugmoa_70138 <- merge(mat_70138, drugmoa, all.x=T, by='pert_iname')
drugmoa_cls_70138 <- merge(drugmoa_70138, cls, by='pert_iname')
drugmoa_cls_70138 <- drugmoa_cls_70138[order(drugmoa_cls_70138$cls), ]
drugmoa_cls_70138 <- drugmoa_cls_70138[order(drugmoa_cls_70138$average), ]
summary(drugmoa_cls_70138$average)
sort(table(drugmoa_cls_70138[drugmoa_cls_70138$average>3, ]$cls))
sort(table(drugmoa_cls_70138[drugmoa_cls_70138$average>3, ]$MoA_merge))
write.table(drugmoa_cls_70138, 
            "07_plotting_v2/05_drugITSprofiler_CombatRef_metaPlusFilter_corr2/drug_overall_stat/drugmoa_cls_70138.txt",
            sep = '\t', quote = F, row.names=F)
drugmoa_annot = unique(drugmoa_cls_70138[c(1,12)])

mat_70138[is.na(mat_70138)] = -666
cls=cls[order(cls$cls),]
mat_70138 = mat_70138[cls$pert_iname,]
rangevalue = max(abs(min(dat_70138$score)), max(dat_70138$score))
bk <- c(seq(-rangevalue,-0.5,by=0.5),seq(0,rangevalue, by=0.5))

cls$clsShow = "other"
cls[cls$cls == 'cls6',]$clsShow='cls6'
cls[cls$cls == 'cls15',]$clsShow='cls15'
cls[cls$cls == 'cls24',]$clsShow='cls24'
cls[cls$cls == 'cls12',]$clsShow='cls12'
cls[cls$cls == 'cls4',]$clsShow='cls4'
mat_70138 = mat_70138[order(mat_70138$average, decreasing = T),]
cls=cls[order(cls$cls),]
mat_70138 = mat_70138[cls$pert_iname,]

pdf('07_plotting_v2/05_drugITSprofiler_CombatRef_metaPlusFilter_corr2/drug_overall_stat/moa_cls_70138_heatmap.pdf', width = 15, height = 7)
pheatmap(t(mat_70138[-c(1,11)]),
         scale = "none", 
         cluster_row = F, cluster_col = F,
         border="white",
         cellwidth = 0.5,
         cellheight = 20,
         fontsize = 4,
         annotation_col = cls[-c(1,2)], # drugproof_annot['proof'],
         # cutree_cols = 30, 
         treeheight_col=30,
         # cutree_cols = 5,
         color = c('grey95', colorRampPalette(colors = c("#228383","white"))(length(bk)/2),colorRampPalette(colors = c("white","#B33333"))(length(bk)/2)),
         legend_breaks=c(-666, seq(-rangevalue, rangevalue,2)),
         breaks=bk)
dev.off()

mat_70138 = mat_70138[order(mat_70138$average, decreasing = T),]
mat_70138 = mat_70138[c("pert_iname","LIHC_Primary","LUAD_Primary",
                        "PAAD_Primary","SKCM_Primary",
                        "READ_Primary","PRAD_Primary",
                        "SKCM_Metastatic","BRCA_Primary",
                        "COAD_Primary", "average")]

pdf('07_plotting_v2/05_drugITSprofiler_CombatRef_metaPlusFilter_corr2/drug_overall_stat/moa_topcls_70138_heatmap.pdf', width = 15, height = 7)
clsShow = cls[cls$cls == 'cls24',]
mat_70138_clsShow = mat_70138[is.element(rownames(mat_70138),clsShow$pert_iname),]
mat_70138_clsShow = mat_70138_clsShow[order(mat_70138_clsShow$average, decreasing = T),]
drugmoa_cls24 = drugmoa_annot[is.element(drugmoa_annot$pert_iname,clsShow$pert_iname), ]
pheatmap(t(mat_70138_clsShow[-c(1,11)]),
         scale = "none", 
         cluster_row = F, cluster_col = F,
         border="white",
         cellwidth = 8,
         cellheight = 20,
         fontsize = 7,
         # annotation_col = drugproof_annot['proof'],
         # cutree_cols = 30, 
         treeheight_col=30,
         # cutree_cols = 5,
         color = c('grey95', colorRampPalette(colors = c("#228383","white"))(length(bk)/2),colorRampPalette(colors = c("white","#B33333"))(length(bk)/2)),
         legend_breaks=c(-666, seq(-rangevalue, rangevalue,2)),
         breaks=bk,
         main='cls24')
# dev.off()
sort(table(drugmoa_cls_70138[drugmoa_cls_70138$cls == 'cls24', ]$MoA_merge))
summary(drugmoa_cls_70138[drugmoa_cls_70138$cls == 'cls24', ]$average)


clsShow = cls[cls$cls == 'cls4',]
mat_70138_clsShow = mat_70138[is.element(rownames(mat_70138),clsShow$pert_iname),]
mat_70138_clsShow = mat_70138_clsShow[order(mat_70138_clsShow$average, decreasing = T),]
drugmoa_cls4 = drugmoa_annot[is.element(drugmoa_annot$pert_iname,clsShow$pert_iname), ]
# pdf('07_plotting_v2/05_drugITSprofiler_CombatRef_metaPlusFilter_corr2/drug_overall_stat/moa_cls4_70138_heatmap.pdf', width = 15, height = 7)
pheatmap(t(mat_70138_clsShow[-c(1,11)]),
         scale = "none", 
         cluster_row = F, cluster_col = F,
         border="white",
         cellwidth = 8,
         cellheight = 20,
         fontsize = 7,
         # annotation_col = drugproof_annot['proof'],
         # cutree_cols = 30, 
         treeheight_col=30,
         # cutree_cols = 5,
         color = c('grey95', colorRampPalette(colors = c("#228383","white"))(length(bk)/2),colorRampPalette(colors = c("white","#B33333"))(length(bk)/2)),
         legend_breaks=c(-666, seq(-rangevalue, rangevalue,2)),
         breaks=bk,
         main='cls4')
# dev.off()
sort(table(drugmoa_cls_70138[drugmoa_cls_70138$cls == 'cls4', ]$MoA_merge))
summary(drugmoa_cls_70138[drugmoa_cls_70138$cls == 'cls4', ]$average)


clsShow = cls[cls$cls == 'cls12',]
mat_70138_clsShow = mat_70138[is.element(rownames(mat_70138),clsShow$pert_iname),]
mat_70138_clsShow = mat_70138_clsShow[order(mat_70138_clsShow$average, decreasing = T),]
drugmoa_cls12 = drugmoa_annot[is.element(drugmoa_annot$pert_iname,clsShow$pert_iname), ]
# pdf('07_plotting_v2/05_drugITSprofiler_CombatRef_metaPlusFilter_corr2/drug_overall_stat/moa_cls12_70138_heatmap.pdf', width = 15, height = 7)
pheatmap(t(mat_70138_clsShow[-c(1,11)]),
         scale = "none", 
         cluster_row = F, cluster_col = F,
         border="white",
         cellwidth = 8,
         cellheight = 20,
         fontsize = 7,
         # annotation_col = drugproof_annot['proof'],
         # cutree_cols = 30, 
         treeheight_col=30,
         # cutree_cols = 5,
         color = c('grey95', colorRampPalette(colors = c("#228383","white"))(length(bk)/2),colorRampPalette(colors = c("white","#B33333"))(length(bk)/2)),
         legend_breaks=c(-666, seq(-rangevalue, rangevalue,2)),
         breaks=bk,
         main='cls12')
# dev.off()
sort(table(drugmoa_cls_70138[drugmoa_cls_70138$cls == 'cls12', ]$MoA_merge))
summary(drugmoa_cls_70138[drugmoa_cls_70138$cls == 'cls12', ]$average)

clsShow = cls[cls$cls == 'cls6',]
mat_70138_clsShow = mat_70138[is.element(rownames(mat_70138),clsShow$pert_iname),]
mat_70138_clsShow = mat_70138_clsShow[order(mat_70138_clsShow$average, decreasing = T),]
drugmoa_cls6 = drugmoa_annot[is.element(drugmoa_annot$pert_iname,clsShow$pert_iname), ]
# pdf('07_plotting_v2/05_drugITSprofiler_CombatRef_metaPlusFilter_corr2/drug_overall_stat/moa_cls6_70138_heatmap.pdf', width = 15, height = 7)
pheatmap(t(mat_70138_clsShow[-c(1,11)]),
         scale = "none", 
         cluster_row = F, cluster_col = F,
         border="white",
         cellwidth = 4,
         cellheight = 20,
         fontsize = 7,
         # annotation_col = drugproof_annot['proof'],
         # cutree_cols = 30, 
         treeheight_col=30,
         # cutree_cols = 5,
         color = c('grey95', colorRampPalette(colors = c("#228383","white"))(length(bk)/2),colorRampPalette(colors = c("white","#B33333"))(length(bk)/2)),
         legend_breaks=c(-666, seq(-rangevalue, rangevalue,2)),
         breaks=bk,
         main='cls6')
# dev.off()
sort(table(drugmoa_cls_70138[drugmoa_cls_70138$cls == 'cls6', ]$MoA_merge))
summary(drugmoa_cls_70138[drugmoa_cls_70138$cls == 'cls6', ]$average)

clsShow = cls[cls$cls == 'cls15',]
mat_70138_clsShow = mat_70138[is.element(rownames(mat_70138),clsShow$pert_iname),]
mat_70138_clsShow = mat_70138_clsShow[order(mat_70138_clsShow$average, decreasing = T),]
drugmoa_cls15 = drugmoa_annot[is.element(drugmoa_annot$pert_iname,clsShow$pert_iname), ]
# pdf('07_plotting_v2/05_drugITSprofiler_CombatRef_metaPlusFilter_corr2/drug_overall_stat/moa_cls15_70138_heatmap.pdf', width = 15, height = 7)
pheatmap(t(mat_70138_clsShow[-c(1,11)]),
         scale = "none", 
         cluster_row = F, cluster_col = F,
         border="white",
         cellwidth = 8,
         cellheight = 20,
         fontsize = 7,
         # annotation_col = drugproof_annot['proof'],
         # cutree_cols = 30, 
         treeheight_col=30,
         # cutree_cols = 5,
         color = c('grey95', colorRampPalette(colors = c("#228383","white"))(length(bk)/2),colorRampPalette(colors = c("white","#B33333"))(length(bk)/2)),
         legend_breaks=c(-666, seq(-rangevalue, rangevalue,2)),
         breaks=bk,
         main='cls15')
dev.off()
sort(table(drugmoa_cls_70138[drugmoa_cls_70138$cls == 'cls15', ]$MoA_merge))
summary(drugmoa_cls_70138[drugmoa_cls_70138$cls == 'cls15', ]$average)


sort(table(drugmoa_cls4$MoA_merge))
sort(table(drugmoa_cls24$MoA_merge))
sort(table(drugmoa_cls12$MoA_merge))
sort(table(drugmoa_cls6$MoA_merge))
sort(table(drugmoa_cls15$MoA_merge))
# 92742 ------------------------------------------------------------------------
mat_92742 <- dcast(dat_92742, pert_iname ~ cancer , value.var = 'score' )
rownames(mat_92742)= mat_92742$pert_iname
mat_92742[is.na(mat_92742)] = -666

cls <- (cutree(hclust(dist(mat_92742[-c(4,5)])), k=50))
cls = data.frame(pert_iname = names(cls), cls=cls)
cls=cls[order(cls$cls),]
cls$cls = paste0('cls', cls$cls)
mat_92742 = mat_92742[cls$pert_iname,]

rangevalue = max(abs(min(dat_92742$score)), max(dat_92742$score))
bk <- c(seq(-rangevalue,-0.5,by=0.5),seq(0,rangevalue, by=0.5))

pheatmap(mat_92742[-c(1)],
         scale = "none", 
         cluster_row = T, cluster_col = F,
         border="white",
         cellwidth = 0.5,
         cellheight = 20,
         fontsize = 4,
         annotation_col = cls[-c(1,2)], # drugproof_annot['proof'],
         # cutree_cols = 30, 
         treeheight_col=30,
         # cutree_cols = 5,
         color = c('grey95', colorRampPalette(colors = c("purple","white"))(length(bk)/2),colorRampPalette(colors = c("white","#ff4500"))(length(bk)/2)),
         legend_breaks=c(-666, seq(-rangevalue, rangevalue,2)),
         breaks=bk)
dev.off()

drugmoa <- read.csv('07_plotting_v2/data/drugMoA_mergeNew.txt', sep='\t', header = T)
drugmoa <- unique(drugmoa[-c(1,3,5)])
mat_92742 <- dcast(dat_92742, pert_iname ~ cancer , value.var = 'score' )
rownames(mat_92742)= mat_92742$pert_iname
# mat_92742[is.na(mat_92742)] = -666

mat_92742$average = apply(mat_92742[-1], 1, function(x)mean(x[-which(is.na(x))]))
drugmoa_92742 <- merge(mat_92742, drugmoa, all.x=T, by='pert_iname')
drugmoa_cls_92742 <- merge(drugmoa_92742, cls, by='pert_iname')
drugmoa_cls_92742 <- drugmoa_cls_92742[order(drugmoa_cls_92742$cls), ]
drugmoa_cls_92742 <- drugmoa_cls_92742[order(drugmoa_cls_92742$average), ]
summary(drugmoa_cls_92742$average)
sort(table(drugmoa_cls_92742[drugmoa_cls_92742$average>5, ]$cls))
sort(table(drugmoa_cls_92742[drugmoa_cls_92742$average>3, ]$MoA_merge))


write.table(drugmoa_cls_92742, 
            "07_plotting_v2/05_drugITSprofiler_CombatRef_metaPlusFilter_corr2/drug_overall_stat/drugmoa_cls_92742.txt",
            sep = '\t', quote = F, row.names=F)


# -------------------------------------------------------------------------



mat_92742[is.na(mat_92742)] = -666
cls=cls[order(cls$cls),]
mat_92742 = mat_92742[cls$pert_iname,]
rangevalue = max(abs(min(dat_92742$score)), max(dat_92742$score))
bk <- c(seq(-rangevalue, -0.5, by=0.5),seq(0, rangevalue, by=0.5))

cls$clsShow = "other"
cls[cls$cls == 'cls14',]$clsShow='cls14'
cls[cls$cls == 'cls1',]$clsShow='cls1'
cls[cls$cls == 'cls29',]$clsShow='cls29'
cls[cls$cls == 'cls19',]$clsShow='cls19'
cls[cls$cls == 'cls3',]$clsShow='cls3'
mat_92742 = mat_92742[order(mat_92742$average, decreasing = T),]
cls=cls[order(cls$cls),]
mat_92742 = mat_92742[cls$pert_iname,]

pdf('07_plotting_v2/05_drugITSprofiler_CombatRef_metaPlusFilter_corr2/drug_overall_stat/moa_cls_92742_heatmap.pdf', width = 28, height = 10)
pheatmap(t(mat_92742[-c(1,14)]),
         scale = "none", 
         cluster_row = T, cluster_col = F,
         border="white",
         cellwidth = 0.1,
         cellheight = 20,
         fontsize = 4,
         annotation_col = cls[-c(1,2)], # drugproof_annot['proof'],
         # cutree_cols = 30, 
         treeheight_col=30,
         # cutree_cols = 5,
         color = c('grey95', colorRampPalette(colors = c("#228383","white"))(length(bk)/2),colorRampPalette(colors = c("white","#B33333"))(length(bk)/2)),
         legend_breaks=c(-666, seq(-rangevalue, rangevalue,2)),
         breaks=bk)
dev.off()


cls$clsShow = "other"
cls[cls$cls == 'cls14',]$clsShow='cls14'
pheatmap(t(mat_92742[-c(1,14)]),
         scale = "none", 
         cluster_row = T, cluster_col = F,
         border="white",
         cellwidth = 0.5,
         cellheight = 20,
         fontsize = 4,
         annotation_col = cls[-c(1,2)], # drugproof_annot['proof'],
         # cutree_cols = 30, 
         treeheight_col=30,
         # cutree_cols = 5,
         color = c('grey95', colorRampPalette(colors = c("purple","white"))(length(bk)/2),colorRampPalette(colors = c("white","#ff4500"))(length(bk)/2)),
         legend_breaks=c(-666, seq(-rangevalue, rangevalue,2)),
         breaks=bk)
sort(table(drugmoa_cls_92742[drugmoa_cls_92742$cls == 'cls14', ]$MoA_merge))
summary(drugmoa_cls_92742[drugmoa_cls_92742$cls == 'cls14', ]$average)


cls$clsShow = "other"
cls[cls$cls == 'cls1',]$clsShow='cls1'
pheatmap(t(mat_92742[-c(1,14)]),
         scale = "none", 
         cluster_row = T, cluster_col = F,
         border="white",
         cellwidth = 0.5,
         cellheight = 20,
         fontsize = 4,
         annotation_col = cls[-c(1,2)], # drugproof_annot['proof'],
         # cutree_cols = 30, 
         treeheight_col=30,
         # cutree_cols = 5,
         color = c('grey95', colorRampPalette(colors = c("purple","white"))(length(bk)/2),colorRampPalette(colors = c("white","#ff4500"))(length(bk)/2)),
         legend_breaks=c(-666, seq(-rangevalue, rangevalue,2)),
         breaks=bk)
sort(table(drugmoa_cls_92742[drugmoa_cls_92742$cls == 'cls1', ]$MoA_merge))
summary(drugmoa_cls_92742[drugmoa_cls_92742$cls == 'cls1', ]$average)


cls$clsShow = "other"
cls[cls$cls == 'cls29',]$clsShow='cls29'
pheatmap(t(mat_92742[-c(1,14)]),
         scale = "none", 
         cluster_row = T, cluster_col = F,
         border="white",
         cellwidth = 0.5,
         cellheight = 20,
         fontsize = 4,
         annotation_col = cls[-c(1,2)], # drugproof_annot['proof'],
         # cutree_cols = 30, 
         treeheight_col=30,
         # cutree_cols = 5,
         color = c('grey95', colorRampPalette(colors = c("purple","white"))(length(bk)/2),colorRampPalette(colors = c("white","#ff4500"))(length(bk)/2)),
         legend_breaks=c(-666, seq(-rangevalue, rangevalue,2)),
         breaks=bk)
sort(table(drugmoa_cls_92742[drugmoa_cls_92742$cls == 'cls29', ]$MoA_merge))
summary(drugmoa_cls_92742[drugmoa_cls_92742$cls == 'cls29', ]$average)
# cls29


cls$clsShow = "other"
cls[cls$cls == 'cls19',]$clsShow='cls19'
pheatmap(t(mat_92742[-c(1,14)]),
         scale = "none", 
         cluster_row = T, cluster_col = F,
         border="white",
         cellwidth = 0.5,
         cellheight = 20,
         fontsize = 4,
         annotation_col = cls[-c(1,2)], # drugproof_annot['proof'],
         # cutree_cols = 30, 
         treeheight_col=30,
         # cutree_cols = 5,
         color = c('grey95', colorRampPalette(colors = c("purple","white"))(length(bk)/2),colorRampPalette(colors = c("white","#ff4500"))(length(bk)/2)),
         legend_breaks=c(-666, seq(-rangevalue, rangevalue,2)),
         breaks=bk)
sort(table(drugmoa_cls_92742[drugmoa_cls_92742$cls == 'cls19', ]$MoA_merge))
summary(drugmoa_cls_92742[drugmoa_cls_92742$cls == 'cls19', ]$average)


cls$clsShow = "other"
cls[cls$cls == 'cls3',]$clsShow='cls3'
pheatmap(t(mat_92742[-c(1,14)]),
         scale = "none", 
         cluster_row = T, cluster_col = F,
         border="white",
         cellwidth = 0.5,
         cellheight = 20,
         fontsize = 4,
         annotation_col = cls[-c(1,2)], # drugproof_annot['proof'],
         # cutree_cols = 30, 
         treeheight_col=30,
         # cutree_cols = 5,
         color = c('grey95', colorRampPalette(colors = c("purple","white"))(length(bk)/2),colorRampPalette(colors = c("white","#ff4500"))(length(bk)/2)),
         legend_breaks=c(-666, seq(-rangevalue, rangevalue,2)),
         breaks=bk)
sort(table(drugmoa_cls_92742[drugmoa_cls_92742$cls == 'cls3', ]$MoA_merge))
summary(drugmoa_cls_92742[drugmoa_cls_92742$cls == 'cls3', ]$average)



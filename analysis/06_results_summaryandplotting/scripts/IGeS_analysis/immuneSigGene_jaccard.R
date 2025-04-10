

cancer = "SKCM"
sample_type = "Metastatic"
tumor_purity_method = "TUMERIC"
datatype="allgenes"
num_gene = 200
workpath = work_path
immuneSigGSpath = "02_tumorGene_immuneSig_correlation/data/immune_signatures/originalSignaturesGenes/"
savepath = "r1_drug_immuneSig_correlation_analysis/"

# LOAD original immune signature gene sets ---------------------------------
immunesigGS1 <- as.data.frame(read_excel(paste0(immuneSigGSpath, "immunesig_original_used.xlsx"), sheet = 1))
immunesigGS1 <- immunesigGS1[-grep("various", immunesigGS1$Genes), ]
immunesigGS1 <- immunesigGS1[-grep("SignatureMatrix", immunesigGS1$Genes), ]
immunesigGS1 <- immunesigGS1[-grep("ReferenceCellExpression", immunesigGS1$Genes), ]
immunesigGS1 <- immunesigGS1[-grep("mergedSignatures", immunesigGS1$Genes), ]
immunesigGS1 <- immunesigGS1[-grep("L22Matrix", immunesigGS1$Genes), ]

immunesigGS1 <- immunesigGS1[-c(2:5)]
immunesigGS1Name <- as.list(immunesigGS1$immune_sig)
immunesigGS1 <- lapply(immunesigGS1Name, function(x){
  y=immunesigGS1[immunesigGS1$immune_sig == x,][-1]
  unlist(y[-which(is.na(y))])
})
names(immunesigGS1) = unlist(immunesigGS1Name)

immunesigGS2 <- read_excel(paste0(immuneSigGSpath, "immunesig_original_used.xlsx"), sheet = 2)
immunesigGS2 <- as.data.frame(immunesigGS2[immunesigGS2$CancerType == cancer, ])
immunesigGS2Name <- as.list(immunesigGS2$Signatures)
immunesigGS2 <- lapply(immunesigGS2Name, function(x){
  y=immunesigGS2[immunesigGS2$Signatures == x,][-1]
  unlist(y[-which(is.na(y))])
})
names(immunesigGS2) = paste0(unlist(immunesigGS2Name),"_ConsensusTMEgsva")
immunesigGS = c(immunesigGS1, immunesigGS2)
immunesigGSname = data.frame(signame = names(immunesigGS), immune_sig = tolower(names(immunesigGS)))
immunesigGSname = immunesig_name_unify(immunesigGSname)
names(immunesigGS) = immunesigGSname$immune_sig


ITSselected2 <- read.table("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/06_results_summaryandplotting/data/immune_sig_selected_new/selected_ITS.txt",
                           sep = "\t", header = T)

immunesigInfo <- as.data.frame(read_excel(paste0(immuneSigInfoPath, "immune_sig_hieracy_new.xlsx"), sheet = 1))
immunesigInfo <- immunesigInfo[c(1,4,5:9)]
immunesigInfo$immune_sig = tolower(immunesigInfo$immune_sig)
immunesigInfo = immunesig_name_unify(immunesigInfo)
immunesigInfoSelected = immunesigInfo[is.element(immunesigInfo$immune_sig, names(immunesigGS)), ]

TMEcellSig = immunesigInfoSelected[immunesigInfoSelected$Classification == "TME",]
immuneSig = immunesigInfoSelected[immunesigInfoSelected$Classification == "Immune signatures",]
icb_predictor = immunesigInfoSelected[immunesigInfoSelected$Classification == "ICB predictor",]
ICBresisSig = immunesigInfoSelected[immunesigInfoSelected$Classification == "ICB resistant signatures",]
TMESig = immunesigInfoSelected[immunesigInfoSelected$Classification == "TMEsignatures",]



immunesigGS_selected = immunesigGS[intersect(names(immunesigGS), immuneSig$immune_sig)]

jaccard_index <- function(x, y){
  length(intersect(x,y))/length(union(x,y))
}

Otsuka_Ochiai_coefficient <- function(x, y){
    length(intersect(x, y))/sqrt(length(x) *length(y))
}

maxoverlap <- function(x, y){
    max(length(intersect(x, y))/length(x), length(intersect(x, y))/length(y))
}
TMESig_p_jcindex <- list()
for(i in 1:length(immunesigGS_selected)){
  x = immunesigGS_selected[[i]]
  TMESig_p_jcindex[[i]] = sapply(immunesigGS_selected, function(y){
    jaccard_index(x,y)
    # length(intersect(x, y))
  })
}
names(TMESig_p_jcindex) = names(immunesigGS_selected)

# plot JC index to see the pattern -----------------------------------------
TMESig_p_jcindex = do.call(cbind, TMESig_p_jcindex)
TMESig_p_jcindex = TMESig_p_jcindex[sort(colnames(TMESig_p_jcindex)), ]
TMESig_p_jcindex = TMESig_p_jcindex[, sort(colnames(TMESig_p_jcindex))]


pdf(paste0("06_results_summaryandplotting/data/immunesig_jaccard_index_spearman_r",s , "_ITSproof_",ITSproof,"/oriTMESig_jci.pdf"), 
    width=18, height = 16)
# rownames(ITS_TMEcellSig_p_jcindex) = paste0(rownames(ITS_TMEcellSig_p_jcindex), " (", unlist(ng), ")")
    num_gene_immunesigGS_selected = lapply(immunesigGS_selected,length)
    ng = num_gene_immunesigGS_selected[colnames(TMESig_p_jcindex)]
    colnames(TMESig_p_jcindex) = paste0(colnames(TMESig_p_jcindex), " (", unlist(ng), ")")

pheatmap::pheatmap(TMESig_p_jcindex,
                   cluster_cols = T,
                   cluster_rows = T,
                   display_numbers = T,
                   number_format = "%.2f",
                   fontsize = 5, 
                   annotation_row = immuneSig[c(3,4)],
                   color = colorRampPalette(c("steelblue", "white", "firebrick3"))(20),
                   border_color = 'grey50' )
dev.off()


TMESig_p_Ochiai <- list()
for(i in 1:length(immunesigGS_selected)){
  x = immunesigGS_selected[[i]]
  TMESig_p_Ochiai[[i]] = sapply(immunesigGS_selected, function(y){
    Otsuka_Ochiai_coefficient(x,y)
    # length(intersect(x, y))
  })
}
names(TMESig_p_Ochiai) = names(immunesigGS_selected)

# plot JC index to see the pattern -----------------------------------------
TMESig_p_Ochiai = do.call(cbind, TMESig_p_Ochiai)
TMESig_p_Ochiai = TMESig_p_Ochiai[sort(colnames(TMESig_p_Ochiai)), ]
TMESig_p_Ochiai = TMESig_p_Ochiai[, sort(colnames(TMESig_p_Ochiai))]
TMESig_p_Ochiai[is.na(TMESig_p_Ochiai)] = 0


pdf(paste0("06_results_summaryandplotting/data/immunesig_jaccard_index_spearman_r",s , "_ITSproof_",ITSproof,"/oriTMESig_Ochiai.pdf"), 
    width=18, height = 16)
# rownames(ITS_TMEcellSig_p_Ochiai) = paste0(rownames(ITS_TMEcellSig_p_Ochiai), " (", unlist(ng), ")")
    num_gene_immunesigGS_selected = lapply(immunesigGS_selected,length)
    ng = num_gene_immunesigGS_selected[colnames(TMESig_p_Ochiai)]
    colnames(TMESig_p_Ochiai) = paste0(colnames(TMESig_p_Ochiai), " (", unlist(ng), ")")

tmp = pheatmap::pheatmap(TMESig_p_Ochiai,
                   cluster_cols = T,
                   cluster_rows = T,
                   display_numbers = T,
                   number_format = "%.2f",
                   fontsize = 5, 
                   annotation_row = immuneSig[c(3,4)],
                   color = colorRampPalette(c("steelblue", "white", "firebrick3"))(20),
                   border_color = 'grey50' )
dev.off()


# d = dist(TMESig_p_Ochiai, method = 'euclidean')
# tree = hclust(d, method = 'complete')
tree = tmp$tree_col
treedf = as.data.frame(as.matrix(do.call(cbind, tree[c(4,2:3)])))
treedf$height = as.numeric(treedf$height)
treedf$order = as.numeric(treedf$order)
names(treedf) = c("immune_sig", "height","order")
treedf$immune_sig = as.character(treedf$immune_sig)
v = cutree(tree, 16)[tree$order]
gaps = which((v[-1] - v[-length(v)]) != 0)
its_cls <- as.data.frame(as.matrix(v))
its_cls$immune_sig <- str_split(rownames(its_cls),'[.]',simplify = T)[,1]
its_cls = inner_join(treedf, its_cls)
its_cls = its_cls[order(its_cls$V1),]
its_cls$V1 = as.character(its_cls$V1)
rownames(its_cls) = its_cls$immune_sig

pdf(paste0("06_results_summaryandplotting/data/immunesig_jaccard_index_spearman_r",s , "_ITSproof_",ITSproof,"/oriTMESig_Ochiai_3.pdf"), 
    width=18, height = 16)
pheatmap::pheatmap(TMESig_p_Ochiai[rownames(its_cls),rownames(its_cls)],
                   cluster_cols = F,
                   cluster_rows = F,
                   display_numbers = T,
                   number_format = "%.2f",
                   fontsize = 5, 
                   annotation_row = its_cls[4],# immuneSig[c(3,4)],
                   annotation_col = its_cls[4],
                   color = colorRampPalette(c("steelblue", "white", "firebrick3"))(20),
                   border_color = 'grey50' )
dev.off()




pdf(paste0("06_results_summaryandplotting/data/immunesig_jaccard_index_spearman_r",s , "_ITSproof_",ITSproof,"/ITSTMESig_TEST_n.pdf"), 
    width=18, height = 16)
pheatmap::pheatmap(ITS_immuneSig_n_jcindex[rownames(its_cls),rownames(its_cls)],
                   cluster_cols = F,
                   cluster_rows = F,
                   display_numbers = T,
                   number_format = "%.2f",
                   fontsize = 5, 
                   annotation_row = its_cls[4],# immuneSig[c(3,4)],
                   annotation_col = its_cls[4],
                   color = colorRampPalette(c("steelblue", "white", "firebrick3"))(20),
                   border_color = 'grey50' )
dev.off()


TMESig_geneoverlap <- list()
for(i in 1:length(immunesigGS_selected)){
  x = immunesigGS_selected[[i]]
  TMESig_geneoverlap[[i]] = sapply(immunesigGS_selected, function(y){
    # jaccard_index(x,y)
   maxoverlap(x, y)
  })
}
names(TMESig_geneoverlap) = names(immunesigGS_selected)

# plot JC index to see the pattern -----------------------------------------
TMESig_geneoverlap = do.call(cbind, TMESig_geneoverlap)
TMESig_geneoverlap = TMESig_geneoverlap[sort(colnames(TMESig_geneoverlap)), ]
TMESig_geneoverlap = TMESig_geneoverlap[, sort(colnames(TMESig_geneoverlap))]
TMESig_geneoverlap[is.na(TMESig_geneoverlap)] = 0

pdf(paste0("06_results_summaryandplotting/data/immunesig_jaccard_index_spearman_r",s , "_ITSproof_",ITSproof,"/oriTMESig_geneoverlap.pdf"), 
    width=18, height = 16)
    num_gene_immunesigGS_selected = lapply(immunesigGS_selected,length)
    ng = num_gene_immunesigGS_selected[colnames(TMESig_geneoverlap)]
    colnames(TMESig_geneoverlap) = paste0(colnames(TMESig_geneoverlap), " (", unlist(ng), ")")

    pheatmap::pheatmap(TMESig_geneoverlap,
                    cluster_cols = T,
                    cluster_rows = T,
                    display_numbers = T,
                    number_format = "%.2f",
                    fontsize = 5, 
                    annotation_row = immuneSig[c(3,4)],
                    color = colorRampPalette(c("steelblue", "white", "firebrick3"))(20),
                    border_color = 'grey50' )
dev.off()


# ------------------------------------------------------------------------------
TMEcellSigGS_selected = immunesigGS[intersect(names(immunesigGS), TMEcellSig$immune_sig)]

jaccard_index <- function(x, y){
  length(intersect(x,y))/length(union(x,y))
}


TMEcellSig_p_jcindex <- list()
for(i in 1:length(TMEcellSigGS_selected)){
  x = TMEcellSigGS_selected[[i]]
  TMEcellSig_p_jcindex[[i]] = sapply(TMEcellSigGS_selected, function(y){
    jaccard_index(x,y)
    # length(intersect(x, y))
  })
}
names(TMEcellSig_p_jcindex) = names(TMEcellSigGS_selected)

# plot JC index to see the pattern -----------------------------------------
TMEcellSig_p_jcindex = do.call(cbind, TMEcellSig_p_jcindex)
TMEcellSig_p_jcindex = TMEcellSig_p_jcindex[sort(colnames(TMEcellSig_p_jcindex)), ]
TMEcellSig_p_jcindex = TMEcellSig_p_jcindex[, sort(colnames(TMEcellSig_p_jcindex))]


pdf(paste0("06_results_summaryandplotting/data/immunesig_jaccard_index_spearman_r",s , "_ITSproof_",ITSproof,"/oriTMEcellSig_jci.pdf"), 
    width=18, height = 16)
# rownames(ITS_TMEcellSig_p_jcindex) = paste0(rownames(ITS_TMEcellSig_p_jcindex), " (", unlist(ng), ")")

pheatmap::pheatmap(TMEcellSig_p_jcindex,
                   cluster_cols = T,
                   cluster_rows = T,
                   display_numbers = T,
                   number_format = "%.2f",
                   fontsize = 5, 
                   annotation_row = TMEcellSig[c(3,4)],
                   color = colorRampPalette(c("steelblue", "white", "firebrick3"))(20),
                   border_color = 'grey50' )
dev.off()

TMEcellSig_p_Ochiai <- list()
for(i in 1:length(TMEcellSigGS_selected)){
  x = TMEcellSigGS_selected[[i]]
  TMEcellSig_p_Ochiai[[i]] = sapply(TMEcellSigGS_selected, function(y){
    jaccard_index(x,y)
    # length(intersect(x, y))
  })
}
names(TMEcellSig_p_Ochiai) = names(TMEcellSigGS_selected)

# plot JC index to see the pattern -----------------------------------------
TMEcellSig_p_Ochiai = do.call(cbind, TMEcellSig_p_Ochiai)
TMEcellSig_p_Ochiai = TMEcellSig_p_Ochiai[sort(colnames(TMEcellSig_p_Ochiai)), ]
TMEcellSig_p_Ochiai = TMEcellSig_p_Ochiai[, sort(colnames(TMEcellSig_p_Ochiai))]
TMEcellSig_p_Ochiai[is.na(TMEcellSig_p_Ochiai)] = 0

pdf(paste0("06_results_summaryandplotting/data/immunesig_jaccard_index_spearman_r",s , "_ITSproof_",ITSproof,"/oriTMEcellSig_Ochiai.pdf"), 
    width=18, height = 16)
# rownames(ITS_TMEcellSig_p_Ochiai) = paste0(rownames(ITS_TMEcellSig_p_Ochiai), " (", unlist(ng), ")")

pheatmap::pheatmap(TMEcellSig_p_Ochiai,
                   cluster_cols = T,
                   cluster_rows = T,
                   display_numbers = T,
                   number_format = "%.2f",
                   fontsize = 5, 
                   annotation_row = TMEcellSig[c(3,4)],
                   color = colorRampPalette(c("steelblue", "white", "firebrick3"))(20),
                   border_color = 'grey50' )
dev.off()




TMEcellSig_geneoverlap <- list()
for(i in 1:length(TMEcellSigGS_selected)){
  x = TMEcellSigGS_selected[[i]]
  TMEcellSig_geneoverlap[[i]] = sapply(TMEcellSigGS_selected, function(y){
    # jaccard_index(x,y)
   maxoverlap(x, y)
  })
}
names(TMEcellSig_geneoverlap) = names(TMEcellSigGS_selected)

# plot JC index to see the pattern -----------------------------------------
TMEcellSig_geneoverlap = do.call(cbind, TMEcellSig_geneoverlap)
TMEcellSig_geneoverlap = TMEcellSig_geneoverlap[sort(colnames(TMEcellSig_geneoverlap)), ]
TMEcellSig_geneoverlap = TMEcellSig_geneoverlap[, sort(colnames(TMEcellSig_geneoverlap))]
TMEcellSig_geneoverlap[is.na(TMEcellSig_geneoverlap)] = 0

pdf(paste0("06_results_summaryandplotting/data/immunesig_jaccard_index_spearman_r",s , "_ITSproof_",ITSproof,"/oriTMEcellSig_geneoverlap.pdf"), 
    width=18, height = 16)
    num_gene_TMEcellSigGS_selected = lapply(TMEcellSigGS_selected,length)
    ng = num_gene_TMEcellSigGS_selected[colnames(TMEcellSig_geneoverlap)]
    colnames(TMEcellSig_geneoverlap) = paste0(colnames(TMEcellSig_geneoverlap), " (", unlist(ng), ")")

    pheatmap::pheatmap(TMEcellSig_geneoverlap,
                    cluster_cols = T,
                    cluster_rows = T,
                    display_numbers = T,
                    number_format = "%.2f",
                    fontsize = 5, 
                    annotation_row = TMEcellSig[c(3,4)],
                   color = colorRampPalette(c("steelblue", "white", "firebrick3"))(20),
                    border_color = 'grey50' )
dev.off()

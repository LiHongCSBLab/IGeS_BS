cancer='SKCM'
sample_type = "Metastatic"
cancer='LUAD'
sample_type = "Primary"

immunesig_path <- "02_tumorGene_immuneSig_correlation/data/immune_signatures/"
immunesig_file <- "02_tumorGene_immuneSig_correlation/data/immune_signatures/immune_sigatures.txt"

patient_ID = read.csv("02_tumorGene_immuneSig_correlation/data/TCGA_patient_ID_tumor_type_selected.csv")

# load in integrated immune signature file

immune_sigatures <-  read.table(immunesig_file, sep = "\t", header = T)
immune_sigatures_selected <- immune_sigatures[immune_sigatures$type == paste0(cancer, "_", sample_type),]
immune_sigatures_selected <- immune_sigatures_selected[-which(is.element(names(immune_sigatures_selected), "type"))]
immune_sigatures_selected$SampleID2 <- substr(immune_sigatures_selected$ID, 1, 15)
df_names  <- data.frame(immunesig = names(immune_sigatures_selected), immune_sig = tolower(names(immune_sigatures_selected) ))
df_names <- immunesig_name_unify(df_names)
names(immune_sigatures_selected) <- df_names$immune_sig
rownames(immune_sigatures_selected) <- immune_sigatures_selected$id
# df = immune_sigatures_selected[,intersect(names(immune_sigatures_selected),TMEcellSig$immune_sig)]
df = immune_sigatures_selected[,intersect(names(immune_sigatures_selected),immuneSig$immune_sig)]
df = df[,intersect(names(df), ITSselected2[ITSselected2$flag=='sensitive', ]$immune_sig)]
df[is.na(df)] = 0
dfcor = cor(df, df)
dfcor[dfcor < 0.8] = -1
tmp = pheatmap::pheatmap(dfcor,
                        cluster_cols = T,
                        cluster_rows = T,
                        display_numbers = T,
                        number_format = "%.2f",
                        fontsize = 5, 
                        # annotation_row = immuneSig[c(3,4)],
                        color = colorRampPalette(c("steelblue", "white", "firebrick3"))(20),
                        border_color = 'grey50' )
    tree = tmp$tree_col
    treedf = as.data.frame(as.matrix(do.call(cbind, tree[c(4,2:3)])))
    treedf$height = as.numeric(treedf$height)
    treedf$order = as.numeric(treedf$order)
    names(treedf) = c("immune_sig", "height","order")
    treedf$immune_sig = as.character(treedf$immune_sig)
    v = cutree(tree, k = 10)[tree$order] # , h = 1.5
    # gaps = which((v[-1] - v[-length(v)]) != 0)
    its_cls <- as.data.frame(as.matrix(v))
    its_cls$immune_sig <- str_split(rownames(its_cls),'[.]',simplify = T)[,1]
    its_cls = inner_join(treedf, its_cls)
    its_cls = its_cls[order(its_cls$V1),]
    its_cls$V1 = as.character(its_cls$V1)
    rownames(its_cls) = its_cls$immune_sig

    its_cls

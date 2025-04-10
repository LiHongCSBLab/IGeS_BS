# rm(list=ls())
# Mapping drug DEGs to correlation results
print("start computing")
work_path <- "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/"

setwd(work_path)
options(stringsAsFactors = F)
library(parallel)
library(reshape2)
library(stringr)
library(ggplot2)
source("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/03_drug_immuneSig_enrichment/scripts/02_GSEA_geneset_preparation_function.R")
source("06_results_summaryandplotting/scripts/immunesig_selection_function.R")

load("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/03_drug_immuneSig_enrichment/data/gene_immune_sig_IHC_spearman/for_GSVA/pcor/spearman_ACC_Primary_positive.Rdata")
cancer= "LIHC"
sample_type="Primary"
dir.create( paste0("06_results_summaryandplotting/ITS_correlation/", 
                cancer, "_", sample_type,"/"))
for(i in which(lapply(genelist_p,length) > 1000)){
    # print(paste0(names(genelist_p)[i]," number of genes = ", length(genelist_p[[i]])))
    ITSenrichRes_p <-ITSenrichment(ITS = genelist_p[[i]],
                                        method ="gobp")
    write.csv(ITSenrichRes_p, paste0("06_results_summaryandplotting/ITS_correlation/", 
                cancer, "_", sample_type,"/", names(genelist_p)[i],length(genelist_p[[i]]),"_p.csv"))
    print(head(ITSenrichRes_p[ITSenrichRes_p$p.adjust < 0.01,][c(2)], 20))
}

for(i in which(lapply(genelist_n,length) > 1000)){
    # print(paste0(names(genelist_p)[i]," number of genes = ", length(genelist_p[[i]])))
    ITSenrichRes_n <-ITSenrichment(ITS = genelist_n[[i]],
                                        method ="gobp")
    write.csv(ITSenrichRes_n, paste0("06_results_summaryandplotting/ITS_correlation/", 
                cancer, "_", sample_type,"/", names(genelist_n)[i],length(genelist_n[[i]]),".csv"))
    # print(head(ITSenrichRes[ITSenrichRes$p.adjust < 0.01,][c(1,2)], 10))
}


# if gene expression has high negative correlation witg purity--------------------------------
cancer = 'LIHC'
sample_type = 'Primary'
tumor_purity_method = "TUMERIC"
processedDataPath <- '02_tumorGene_immuneSig_correlation/data/TCGA_processedData/'
processedDataPath <- paste0(processedDataPath, cancer,'/')
TPMMatPath = paste0(processedDataPath, "mixture_file_TCGA_", cancer,"_",
                    sample_type, "_TPM.txt")
mymatrix_filter <- read.csv(TPMMatPath, sep = '\t', row.names = 1)
colnames(mymatrix_filter) <- gsub('\\.', '-',   colnames(mymatrix_filter))
colnames(mymatrix_filter) <- substr(colnames(mymatrix_filter),1,15)

load("01_tumor_purity/data/processedData/TCGA_tumor_purity_merged.Rdata")
tumor_purity <- purity_merged
tumor_purity = tumor_purity[-which(duplicated(tumor_purity$Sample_ID)),]

tumor_purity_filter <- tumor_purity[is.element(tumor_purity$Sample_ID,  colnames(mymatrix_filter)),]
tumor_purity_filter2 <- unique(tumor_purity_filter[c('Sample_ID',tumor_purity_method)])
names(tumor_purity_filter2) <- c('Sample_ID','purity')


# for(j in seq(length(genelist_p))){
#   gs = genelist_p[[j]]

 genelist_p_filtered <- lapply(genelist_p, function(gs){
    # for(i in seq(length(gs))){
    # gs = genelist_p[[22]]
    gs_filtered <- lapply(as.list(gs), function(g){
        # g=gs[i]
        if(is.element(g, rownames(mymatrix_filter))){
            gexpr = as.data.frame(t(data.frame(mymatrix_filter[rownames(mymatrix_filter) == g ,])))
            gexpr$Sample_ID = gsub('\\.', '-',   rownames(gexpr))
            df <- inner_join(gexpr, tumor_purity_filter2)
            names(df) <- c("gene","sample","purity")
            df <- df[!is.na(df$purity), ]

            r = cor.test(df$gene, df$purity, method = "pearson")
                
            # ggplot(data = df, mapping = aes(x = gene, y = purity)) + 
            #         geom_point() +
            #         # geom_abline(intercept = 0, slope = 1) +
            #         theme_bw() +
            #         ggtitle(paste0(names(genelist_p)[j],"- gene ", g,
            #                     "\n spearman's r = ", round(r$estimate,3), 
            #                     "\n P value = ", signif(r$p.value,3)))+
            #         xlab(g) + ylab("purity")

                if(r$p.value < 0.05 & r$estimate < -0.4) {
                    g=NULL
                }else {
                    g=g
                }
        }
        return(g)
    })
   gs_filtered <- unique(unlist(gs_filtered))
   return(gs_filtered)
})

lapply(genelist_p_filtered, length)


for(i in which(lapply(genelist_p_filtered,length) > 1000)){
    # print(paste0(names(genelist_p)[i]," number of genes = ", length(genelist_p[[i]])))
    
    ITSenrichRes_p_filtered <-ITSenrichment(ITS = genelist_p_filtered[[i]],
                                        method ="gobp") # 
    write.csv(ITSenrichRes_p_filtered, paste0("06_results_summaryandplotting/ITS_correlation/", 
                cancer, "_", sample_type,"/", names(genelist_p_filtered)[i],length(genelist_p_filtered[[i]]),"_p_filtered.csv"))
    print(head(ITSenrichRes_p_filtered[ITSenrichRes_p_filtered$p.adjust < 0.01,][c(2)], 10))
}



 genelist_n_filtered <- lapply(genelist_n, function(gs){
    # for(i in seq(length(gs))){
    
    gs_filtered <- lapply(as.list(gs), function(g){
        # g=gs[i]
        if(is.element(g, rownames(mymatrix_filter))){
            gexpr = as.data.frame(t(data.frame(mymatrix_filter[rownames(mymatrix_filter) == g ,])))
            gexpr$Sample_ID = gsub('\\.', '-',   rownames(gexpr))
            df <- inner_join(gexpr, tumor_purity_filter2)
            names(df) <- c("gene","sample","purity")
            df <- df[!is.na(df$purity), ]

            r = cor.test(df$gene, df$purity, method = "pearson")
                
            # ggplot(data = df, mapping = aes(x = gene, y = purity)) + 
            #         geom_point() +
            #         # geom_abline(intercept = 0, slope = 1) +
            #         theme_bw() +
            #         ggtitle(paste0(names(genelist_p)[j],"- gene ", g,
            #                     "\n spearman's r = ", round(r$estimate,3), 
            #                     "\n P value = ", signif(r$p.value,3)))+
            #         xlab(g) + ylab("purity")

                if(r$p.value < 0.05 & r$estimate < -0.2) {
                    g=NULL
                }else {
                    g=g
                }
        }
        return(g)
    })
   gs_filtered <- unique(unlist(gs_filtered))
   return(gs_filtered)
})

lapply(genelist_n_filtered, length)


for(i in which(lapply(genelist_n_filtered,length) > 1000)){
    # print(paste0(names(genelist_p)[i]," number of genes = ", length(genelist_p[[i]])))
    ITSenrichRes_n_filtered <-ITSenrichment(ITS = genelist_n_filtered[[i]],
                                        method ="gobp")
    write.csv(ITSenrichRes_n_filtered, paste0("06_results_summaryandplotting/ITS_correlation/", 
                cancer, "_", sample_type,"/", names(genelist_n_filtered)[i],length(genelist_n_filtered[[i]]),"_n_filtered.csv"))
    print(head(ITSenrichRes_n_filtered[ITSenrichRes_n_filtered$p.adjust < 0.01,][c(1,2)], 10))
}

# dotplot (gene expression and Immunesig) ---------------------------------------------------
immunesig_path <- "02_tumorGene_immuneSig_correlation/data/immune_signatures/"
immunesig_file <- "02_tumorGene_immuneSig_correlation/data/immune_signatures/immune_sigatures.txt"

patient_ID = read.csv("02_tumorGene_immuneSig_correlation/data/TCGA_patient_ID_tumor_type_selected.csv")

# load in integrated immune signature file
names(cor_result) <- tolower(names(cor_result))

immune_sigatures <-  read.table(immunesig_file, sep = "\t", header = T)
immune_sigatures_selected <- immune_sigatures[immune_sigatures$type == paste0(cancer, "_", sample_type),]
immune_sigatures_selected <- immune_sigatures_selected[-which(is.element(names(immune_sigatures_selected), "type"))]
immune_sigatures_selected$SampleID2 <- substr(immune_sigatures_selected$ID, 1, 15)
names(immune_sigatures_selected) <- tolower(names(immune_sigatures_selected) )

processedDataPath <- '02_tumorGene_immuneSig_correlation/data/TCGA_processedData/'
processedDataPath <- paste0(processedDataPath, cancer,'/')
TPMMatPath = paste0(processedDataPath, "mixture_file_TCGA_", cancer,"_",
                    sample_type, "_TPM.txt")
mymatrix_filter <- read.csv(TPMMatPath, sep = '\t', row.names = 1)
colnames(mymatrix_filter) <- gsub('\\.', '-',   colnames(mymatrix_filter))
colnames(mymatrix_filter) <- substr(colnames(mymatrix_filter),1,15)

dir.create(paste0("06_results_summaryandplotting/ITS_correlation/", 
                cancer, "_", sample_type,"/"))
for(j in which(lapply(genelist_p,length) > 1000)[1:8]){
    set.seed(123)                

    p1=list()
    a = sample(seq(length(genelist_p[[j]])), 20)
    tmp = cor_result[[names(genelist_p)[j]]]

    for(x in seq(length(a))){
        i=a[x]
        g=genelist_p[[j]][i]
        gexpr = as.data.frame(t(data.frame(mymatrix_filter[rownames(mymatrix_filter) == g ,])))
        gexpr$sampleid2 = gsub('\\.', '-',   rownames(gexpr))
        df <- inner_join(gexpr, immune_sigatures_selected[c("sampleid2", names(genelist_p)[j])])
        names(df) <- c("gene","sample","immunesig")
        # df$gene <- log2(df$gene +1 )
        # r = cor.test(df[,1], df[,3], method = "spearman")
       r = tmp[rownames(tmp) == g,]
       # print(r$estimate)}
       p1[[x]] <- ggplot(data = df, mapping = aes(x = gene, y = immunesig)) + 
                geom_point()+
                # geom_abline(intercept = 0, slope = 1) +
                theme_bw()+
                ggtitle(paste0(names(genelist_p)[j],"- gene ", g,
                            "\n spearman's r = ", round(unlist(r$cor), 3), 
                            "\n P value = ", signif(unlist(r$p.adj),3)))+
                xlab(g) + ylab(names(genelist_p)[j])
    }

    pdf(paste0("06_results_summaryandplotting/ITS_correlation/", 
                cancer, "_", sample_type,"/", names(genelist_p)[j],"_p.pdf"),
            width = 5,
            height = 5, 
            onefile = TRUE)
    for(x in seq(length(p1))){
        print(p1[[x]])
    }
    dev.off()
    rm(p1)
}


for(j in which(lapply(genelist_n,length) > 1000)){

    p1=list()
    a = sample(seq(length(genelist_n[[j]])), 20)
    tmp = cor_result[[names(genelist_p)[j]]]

    for(x in seq(length(a))){
        i=a[x]
        g=genelist_n[[j]][i]
        gexpr = as.data.frame(t(data.frame(mymatrix_filter[rownames(mymatrix_filter) == g ,])))
        gexpr$sampleid2 = gsub('\\.', '-',   rownames(gexpr))
        df <- inner_join(gexpr, immune_sigatures_selected[c("sampleid2", names(genelist_n)[j])])
        names(df) <- c("gene","sample","immunesig")
        df$gene <- log2(df$gene +1 )
       r = tmp[rownames(tmp) == g,]
        # print(r$estimate)}
       p1[[x]] <- ggplot(data = df, mapping = aes(x = gene, y = immunesig)) + 
                geom_point()+
                # geom_abline(intercept = 0, slope = 1) +
                theme_bw()+
                ggtitle(paste0(names(genelist_n)[j],"- gene ", g,
                            "\n spearman's r = ", round(unlist(r$cor), 3), 
                            "\n P value = ", signif(unlist(r$p.adj),3)))+
                xlab(g) + ylab(names(genelist_n)[j])
    }

    pdf(paste0("06_results_summaryandplotting/ITS_correlation/", 
                cancer, "_", sample_type,"/", names(genelist_n)[j],"_n.pdf"),
            width = 5,
            height = 5, 
            onefile = TRUE)
    for(x in seq(length(p1))){
        print(p1[[x]])
    }
    dev.off()
    rm(p1)
}

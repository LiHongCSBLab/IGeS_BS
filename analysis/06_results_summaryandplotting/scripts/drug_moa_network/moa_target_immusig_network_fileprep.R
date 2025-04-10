
# 4) network plotting
print("prepare files for ploting")
work_path <- "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/"
setwd(work_path)
options(stringsAsFactors = F)


# INPUT: 
#   REQUIRE: DRUG, GENE, CANCER TYPE
#   OPTICAL: PURITY_METHOD (DEFAULT: TUMERIC), CORRELATION_METHOD (DEFAULT: PEARSON)
drug_target_DEG_immunesig_match <- function(cancer,
                                            purity_method,
                                            immusig_path,
                                            datatype,
                                            dataset,
                                            drug_keyword, 
                                            drugMergeFCpath,
                                            resultpath){

work_path <- "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/"

options(stringsAsFactors = F)
library(parallel)
library(reshape2)
library(stringr)
library(dplyr)

# cancer = 'LIHC'
# purity_method = 'TUMERIC'
# immusig_path <- "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/02_tumorGene_immuneSig_correlation/results/"

if(purity_method == 'IHC'){
    load(paste0(immusig_path, "immunecell_TCGAgenes_pcor_pearson_IHC/pcor_pearson_",cancer,"_Primary_IHC.Rdata"))

}else{
    load(paste0(immusig_path, "immunecell_TCGAgenes_result_pearson_", purity_method,"/pcor_pearson_",cancer,"_Primary_", purity_method,".Rdata"))

}


gene_immunesig_sigcor <- lapply(cor_result, function(x){
    if(!is.null(x)){
      # x=cor_result[[1]]
      x = data.frame(r = unlist(x$cor),
                     p.value = unlist(x$p.value),
                     adj.p = unlist(x$p.adj))
      x$gene = row.names(x)
      # x = inner_join(drugTarget['target'],x)
      x <- x[x$adj.p < 0.05,]
      s1 <- summary(x[x$r > 0,]$r)[5]
      s2 <- summary(x[x$r < 0,]$r)[2]
      y <- unique(x[(x$r > s1 | x$r < s2) & x$adj.p < 0.05,])
      # y <- y[!is.na(y$r),]
      return(y)
    }
})


# datatype = 'allgenes'
# dataset = '70138'
drugTargetMoA <- read.table(paste0(work_path, "06_results_summaryandplotting/data/drugTargetMoA_merged.txt"), sep = '\t', header = T)

# drugMergeFCpath = "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/s2_drug_inducted_transcriptome_change/mergeFC_r1/"
drugMergeFCpath = paste0(drugMergeFCpath, "/", datatype, "/", dataset, "/", cancer, "/")
alldruglist <- gsub("_allgeneFC.csv","", list.files(drugMergeFCpath))

drug_list <- drugTargetMoA[grep( drug_keyword, drugTargetMoA$pert_iname), ]
drug_list2 <- intersect(drug_list$drug_index, alldruglist)

drugDEG_immunesig <- lapply(as.list(drug_list2), function(drug){
    allgeneFC <- read.csv(paste0( drugMergeFCpath, drug, "_allgeneFC.csv"))
    names(allgeneFC) <- c('gene','drug_logFC','drug_P.Value','drug_adj.P.Val','drug_index')
    drugGeneFC_immunesig <- lapply(gene_immunesig_sigcor, function(x){
        if(!is.null(x)){
            inner_join(x, allgeneFC)
        }
    })

    drugDEG_immunesig <- lapply(drugGeneFC_immunesig, function(x){
        if(!is.null(x)){
            # x[abs(x$drug_logFC) > 0.5 & x$drug_adj.P.Val < 0.1, ]
            x[x$drug_adj.P.Val < 0.1, ]
        }
    })
    # lapply(drugDEG_immunesig,nrow)
    drugDEG_immunesig <- do.call(rbind, drugDEG_immunesig)
    drugDEG_immunesig$immunesig <- gsub("\\.\\d+","",rownames(drugDEG_immunesig))
    drugDEG_immunesig$immunesig <- gsub("\\.NA","",drugDEG_immunesig$immunesig)
    return(drugDEG_immunesig)
})

names(drugDEG_immunesig) <- drug_list2
save(drugDEG_immunesig, file = paste0(resultpath,"/", cancer, "_drugDEG_immunesig.Rdata"))
    
}



cancer = 'LIHC'
purity_method = 'TUMERIC'
immusig_path <- "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/02_tumorGene_immuneSig_correlation/results/"
datatype = 'allgenes'
dataset = '70138'
drugMergeFCpath = "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/s2_drug_inducted_transcriptome_change/mergeFC_r1/"
resultpath = paste0("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/06_results_summaryandplotting/results/drugDEGfc_immunesig/", purity_method)

drugDEG_immunesig <- drug_target_DEG_immunesig_match(cancer = cancer,
                                            purity_method = purity_method,
                                            immusig_pathimmusig_path,
                                            datatype = datatype,
                                            dataset = dataset,
                                            drug_keyword = 'nib',
                                            drugMergeFCpath = drugMergeFCpath,
                                            resultpath = resultpath)

sorafenib_res <- drugDEG_immunesig[['drug101']]
sorafenib_CD4_res <- sorafenib_res[grep("CD4", sorafenib_res$immunesig), ]
sorafenib_CD4_res_strict <- sorafenib_CD4_res[abs(sorafenib_CD4_res$drug_logFC) > 1, ]
length(unique(sorafenib_CD4_res_strict[sorafenib_CD4_res_strict$drug_logFC > 0 & sorafenib_CD4_res_strict$r >0 ,]$gene))
length(unique(sorafenib_CD4_res_strict[sorafenib_CD4_res_strict$drug_logFC > 0 & sorafenib_CD4_res_strict$r <0 ,]$gene))
length(unique(sorafenib_CD4_res_strict[sorafenib_CD4_res_strict$drug_logFC < 0 & sorafenib_CD4_res_strict$r >0,]$gene))
length(unique(sorafenib_CD4_res_strict[sorafenib_CD4_res_strict$drug_logFC < 0 & sorafenib_CD4_res_strict$r <0,]$gene))


sorafenib_CD8_res <- sorafenib_res[grep("CD8", sorafenib_res$immunesig), ]
sorafenib_CD8_res_strict <- sorafenib_CD8_res[abs(sorafenib_CD8_res$drug_logFC) > 1, ]
length(unique(sorafenib_CD8_res_strict[sorafenib_CD8_res_strict$drug_logFC > 0 & sorafenib_CD8_res_strict$r >0 ,]$gene))
length(unique(sorafenib_CD8_res_strict[sorafenib_CD8_res_strict$drug_logFC > 0 & sorafenib_CD8_res_strict$r <0 ,]$gene))
length(unique(sorafenib_CD8_res_strict[sorafenib_CD8_res_strict$drug_logFC < 0 & sorafenib_CD8_res_strict$r >0,]$gene))
length(unique(sorafenib_CD8_res_strict[sorafenib_CD8_res_strict$drug_logFC < 0 & sorafenib_CD8_res_strict$r <0,]$gene))

sunitinib_res <- drugDEG_immunesig[['drug187']]
sunitinib_CD4_res <- sorafenib_res[grep("CD4", sunitinib_res$immunesig), ]
sunitinib_CD4_res_strict <- sunitinib_CD4_res[abs(sunitinib_CD4_res$drug_logFC) > 0.5, ]
length(unique(sunitinib_CD4_res_strict[sunitinib_CD4_res_strict$drug_logFC > 0 & sunitinib_CD4_res_strict$r >0 ,]$gene))
length(unique(sunitinib_CD4_res_strict[sunitinib_CD4_res_strict$drug_logFC > 0 & sunitinib_CD4_res_strict$r <0 ,]$gene))
length(unique(sunitinib_CD4_res_strict[sunitinib_CD4_res_strict$drug_logFC < 0 & sunitinib_CD4_res_strict$r >0,]$gene))
length(unique(sunitinib_CD4_res_strict[sunitinib_CD4_res_strict$drug_logFC < 0 & sunitinib_CD4_res_strict$r <0,]$gene))


sunitinib_CD8_res <- sorafenib_res[grep("CD8", sunitinib_res$immunesig), ]
sunitinib_CD8_res_strict <- sunitinib_CD8_res[abs(sunitinib_CD8_res$drug_logFC) > 0.5, ]
length(unique(sunitinib_CD8_res_strict[sunitinib_CD8_res_strict$drug_logFC > 0 & sunitinib_CD8_res_strict$r >0 ,]$gene))
length(unique(sunitinib_CD8_res_strict[sunitinib_CD8_res_strict$drug_logFC > 0 & sunitinib_CD8_res_strict$r <0 ,]$gene))
length(unique(sunitinib_CD8_res_strict[sunitinib_CD8_res_strict$drug_logFC < 0 & sunitinib_CD8_res_strict$r >0,]$gene))
length(unique(sunitinib_CD8_res_strict[sunitinib_CD8_res_strict$drug_logFC < 0 & sunitinib_CD8_res_strict$r <0,]$gene))

lapply(drugDEG_immunesig,function(x)nrow(x[x$gene == 'LCK',]))


drug = 'drug128'
# linifanib
drug_res <- drugDEG_immunesig[[drug]]
drug_CD4_res <- drug_res[grep("CD4", drug_res$immunesig), ]
drug_CD4_res_strict <- drug_CD4_res[abs(drug_CD4_res$drug_logFC) > 1.5 & drug_CD4_res$drug_adj.P.Val < 0.01, ]
length(unique(drug_CD4_res_strict[drug_CD4_res_strict$drug_logFC > 0 & drug_CD4_res_strict$r >0 ,]$gene))
length(unique(drug_CD4_res_strict[drug_CD4_res_strict$drug_logFC > 0 & drug_CD4_res_strict$r <0 ,]$gene))
length(unique(drug_CD4_res_strict[drug_CD4_res_strict$drug_logFC < 0 & drug_CD4_res_strict$r >0,]$gene))
length(unique(drug_CD4_res_strict[drug_CD4_res_strict$drug_logFC < 0 & drug_CD4_res_strict$r <0,]$gene))
table(drug_CD4_res_strict$immunesig)
(intersect(drug_CD4_res_strict[drug_CD4_res_strict$immunesig == "T.cell.CD4._epic",]$gene,
                 drug_CD4_res_strict[drug_CD4_res_strict$immunesig == "T.cell.CD4._timer",]$gene))

drug_CD8_res <- drug_res[grep("CD8", drug_res$immunesig), ]
drug_CD8_res_strict <- drug_CD8_res[abs(drug_CD8_res$drug_logFC) > 1.5 & drug_CD8_res$drug_adj.P.Val < 0.01, ]
length(unique(drug_CD8_res_strict[drug_CD8_res_strict$drug_logFC > 0 & drug_CD8_res_strict$r >0 ,]$gene))
length(unique(drug_CD8_res_strict[drug_CD8_res_strict$drug_logFC > 0 & drug_CD8_res_strict$r <0 ,]$gene))
length(unique(drug_CD8_res_strict[drug_CD8_res_strict$drug_logFC < 0 & drug_CD8_res_strict$r >0,]$gene))
length(unique(drug_CD8_res_strict[drug_CD8_res_strict$drug_logFC < 0 & drug_CD8_res_strict$r <0,]$gene))
table(drug_CD8_res_strict$immunesig)

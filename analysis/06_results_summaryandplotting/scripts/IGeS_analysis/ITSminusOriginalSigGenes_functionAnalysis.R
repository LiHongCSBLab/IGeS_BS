load(paste0(immuneSigGSpath, "/ITSGene_purity_correlation/",
                    cancer, "_", sample_type,"_p.Rdata"))
#   genelist_p, genelist_p_filtered, res_p, enrichRes_p_gobp, enrichRes_p_kegg, 

load( paste0(immuneSigGSpath, "/immuneSigGene_purity_correlation/",
                     cancer, "_", sample_type,".Rdata"))
  
# Mapping drug DEGs to correlation results
# work_path <- "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/"
work_path <- "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/"

setwd(work_path)
options(stringsAsFactors = F)
 
library(reshape2)
library(ggplot2)
library(ggplotify)
library(grid)
library(cowplot)
library(parallel)
library(patchwork)
library(dplyr)
library(stringr)
library(data.table)
require(GGally)
library(plot3D)
library(fmsb)
 
source("06_results_summaryandplotting/scripts/plot_for_drug_function.R")
source("06_results_summaryandplotting/scripts/immunesig_selection_function.R")
source("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/scripts/drugITGSscore_lincs.R")
source("06_results_summaryandplotting/scripts/ITS_analysis/ITSexpr_purity_detection.R")
source("06_results_summaryandplotting/scripts/ITS_analysis/ITS_function_detection.R")

# c


immuneSigGene_found_lincs <- function(cancer = "SKCM",
                                      tumor_subtype = "Metastatic",
                                      purity_method = "TUMERIC",
                                      datatype="allgenes",
                                      num_gene = 200,
                                      workpath = "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/", 
                                      immuneSigGSpath = "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/02_tumorGene_immuneSig_correlation/data/immune_signatures/originalSignaturesGenes/",
                                      ITSpath = "03_drug_immuneSig_enrichment/data/gene_immune_sig_TUMERIC_spearman/for_GSVA/pcor/",
                                      savepath = "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/"
                                      ){

    cancer = "SKCM"
    tumor_subtype = "Metastatic"
    purity_method = "TUMERIC"
    datatype="allgenes"
    num_gene = 200
    workpath = work_path
    immuneSigGSpath = "02_tumorGene_immuneSig_correlation/data/immune_signatures/originalSignaturesGenes/"
    savepath = "r1_drug_immuneSig_correlation_analysis/"

    library(readxl)
    library(dplyr)

    dir.create("06_results_summaryandplotting/immuneSigGene_found_lincs")
    dir.create(paste0("06_results_summaryandplotting/immuneSigGene_found_lincs/", purity_method))

    # LOAD GENES IN CCLE -------------------------------------------------------
    geneExprInCCLE <- read.csv("03_drug_immuneSig_enrichment/data/CCLE_geneexpr/geneExprInCCLE.csv")

    protein_coding_genes <- read.table("03_drug_immuneSig_enrichment/data/protein_coding_genes_hg38.txt", sep = '\t', header = T)
    geneExprInCCLE_filtered <- inner_join(geneExprInCCLE, protein_coding_genes)
    geneExprInCCLE_filtered <- geneExprInCCLE_filtered[geneExprInCCLE_filtered$rate < 0.2,]

    # LOAD GENES IN LINCS ------------------------------------------------------
    lincs_gene = read.csv("/picb/bigdata/project/FengFYM/drug_MoA_reanalysis/LINCS_data/GSE70138/file/GSE92742_Broad_LINCS_gene_info.txt", sep = "\t")
    
    # LOAD drug target GENES ---------------------------------------------------
    drugtarget <- read.csv('/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/s3_lincs_drug_MoA/drugTargetMerge.txt', sep = '\t')
    drugtarget <- data.frame(gene_name = unique(drugtarget$target))
    drugtarget_OASIS <- read.csv('/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/s3_lincs_drug_MoA/drugtarget_OASIS_TIDE.csv')

    # LOAD original immune signature gene sets ---------------------------------
    immunesigGS1 <- as.data.frame(read_excel(paste0(immuneSigGSpath, "immunesig_original_used.xlsx"), sheet = 1))
    # immunesigGS1 <- immunesigGS1[-grep("Matrix", immunesigGS1$Genes), ]
    # immunesigGS1 <- immunesigGS1[-grep("ReferenceCellExpression", immunesigGS1$Genes), ]
    # immunesigGS1 <- immunesigGS1[-grep("SignatureMatrix", immunesigGS1$Genes), ]
    immunesigGS1 <- immunesigGS1[-grep("various", immunesigGS1$Genes), ]
    # immunesigGS1 <- immunesigGS1[- which(is.na(immunesigGS1$Genes)), ]
    # immunesigGS1 <- immunesigGS1[which(immunesigGS1$Genes != "NA"), ]

    immunesigGS1 <- immunesigGS1[-c(2:5)]
    immunesigGS1Name <- as.list(immunesigGS1$immune_sig)
    immunesigGS1 <- lapply(immunesigGS1Name, function(x){
     y=immunesigGS1[immunesigGS1$immune_sig == x,][-1]
     unlist(y[-which(is.na(y))])
    })
    names(immunesigGS1) = unlist(immunesigGS1Name)
    # immunesigGS1 = immunesigGS1[-grep("_immunelandscape", names(immunesigGS1))]
    # immunesigGS1 = immunesigGS1[-grep("CIBERSORT", names(immunesigGS1))]
    # immunesigGS1 = immunesigGS1[-grep("C7_Atom", names(immunesigGS1))]
    # immunesigGS1 = immunesigGS1[-grep("timer", names(immunesigGS1))]
    # immunesigGS1 = immunesigGS1[-grep("epic", names(immunesigGS1))]
    # immunesigGS1 = immunesigGS1[-grep("quantiseq", names(immunesigGS1))]
    # immunesigGS1 = immunesigGS1[-grep("_ConsensusTMEgsva", names(immunesigGS1))]
    # immunesigGS1 = immunesigGS1[-grep("IPS", names(immunesigGS1))]

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



    
    # LOAD ITS gene sets
    # workpath = "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/"
    ITSpath = "03_drug_immuneSig_enrichment/data/gene_immune_sig_TUMERIC_spearman/for_GSVA/pcor/"
    # savepath = "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/"

    # dir.create("06_results_summaryandplotting/ITSgene_found_lincs")
    # dir.create(paste0("06_results_summaryandplotting/ITSgene_found_lincs/", purity_method))


    load(paste0(ITSpath, "spearman_",cancer,"_", tumor_subtype, "_positive_200.Rdata"))
    load(paste0(ITSpath, "spearman_",cancer,"_", tumor_subtype, "_negative_200.Rdata"))



    immunesigGS_selected = immunesigGS[intersect(names(immunesigGS), names(genelist_p))]
    immunesigGS_filtered <- lapply(immunesigGS_selected, function(x){
      intersect(x, geneExprInCCLE_filtered$gene_name)
    })
    names(immunesigGS_filtered) = names(immunesigGS_selected)

    genelist_p = genelist_p[names(immunesigGS_selected)]
    genelist_n = genelist_n[names(immunesigGS_selected)]

      print("6. conduct function analysis for each ITS.")
  # 02_tumorGene_immuneSig_correlation/scripts/correlation_dection_plot.R

  # ITSminusOriginalSigGenes ----------------------------
  ITSminusOriginalSigGenes_p <- lapply(as.list(names(genelist_p)), function(x){
      setdiff(genelist_p[[x]], immunesigGS_selected[[x]])
  })
  names(ITSminusOriginalSigGenes_p) = names(genelist_p)
   ITSminusOriginalSigGenes_n <- lapply(as.list(names(genelist_n)), function(x){
      setdiff(genelist_n[[x]], immunesigGS_selected[[x]])
  })
  names(ITSminusOriginalSigGenes_n) = names(genelist_n)


  print("KEGG analysis ############################")
  
  enrichRes_p_kegg=list()
  for(i in seq(length(ITSminusOriginalSigGenes_p))){
    enrichRes_p_kegg[[i]] <-ITSenrichment(ITS = ITSminusOriginalSigGenes_p[[i]],
                                                 method ="kegg")
  }
  
  names(enrichRes_p_kegg) = names(ITSminusOriginalSigGenes_p)
  
  print("GO-BP analysis ############################")
  enrichRes_p_gobp=list()
  for(i in seq(length(ITSminusOriginalSigGenes_p))){
    enrichRes_p_gobp[[i]] <-ITSenrichment(ITS = ITSminusOriginalSigGenes_p[[i]],
                                                 method ="gobp")
  }
  
  names(enrichRes_p_gobp) = names(ITSminusOriginalSigGenes_p)

    dir.create(paste0(immuneSigGSpath, "/ITSminusOriginalSigGenes/"))
    save(   enrichRes_p_kegg, enrichRes_p_gobp, 
         file= paste0(immuneSigGSpath, "/ITSminusOriginalSigGenes/",
                      cancer, "_", sample_type,"_p.Rdata"))
    enrichRes_n_kegg=list()
    for(i in seq(length(ITSminusOriginalSigGenes_n))){
      enrichRes_n_kegg[[i]] <- ITSenrichment(ITS = ITSminusOriginalSigGenes_n[[i]],
                                                    method ="kegg")
    }
    names(enrichRes_n_kegg) = names(ITSminusOriginalSigGenes_n)
    
    
    print("GO-BP analysis ############################")
    
    # 02_tumorGene_immuneSig_correlation/scripts/correlation_dection_plot.R
    enrichRes_n_gobp=list()
    for(i in seq(length(ITSminusOriginalSigGenes_n))){
      enrichRes_n_gobp[[i]] <-ITSenrichment(ITS = ITSminusOriginalSigGenes_n[[i]],
                                                   method ="gobp")
    }
    names(enrichRes_n_gobp) = names(ITSminusOriginalSigGenes_n)

    save(enrichRes_n_gobp, enrichRes_n_gobp, 
         file= paste0(immuneSigGSpath, "/ITSminusOriginalSigGenes/",
                      cancer, "_", sample_type,"_n.Rdata"))
                      

}



immuneSigGene_ITS_found_lincs <- function(cancer = "SKCM",
                                          tumor_subtype = "Metastatic",
                                          purity_method = "TUMERIC",
                                          datatype="allgenes",
                                          num_gene = 200,
                                          r = 0.4,
                                          ITSproof = 2,
                                          cluster_method = c("mclust","hclust"),
                                          similarity_method =  c("jaccard_index", "Otsuka_Ochiai" ),
                                          genevote = 0.5,
                                          workpath = "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/", 
                                          immuneSigGSpath = "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/02_tumorGene_immuneSig_correlation/data/immune_signatures/originalSignaturesGenes/",
                                          # ITSpath = "03_drug_immuneSig_enrichment/data/gene_immune_sig_TUMERIC_spearman/for_GSVA/pcor/",
                                          savepath = "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/"
                                          ){


    # cancer = "SKCM"
    # tumor_subtype = "Metastatic"
    # purity_method = "TUMERIC"
    # datatype="allgenes"
    # num_gene = 200
    # workpath = work_path
    # immuneSigGSpath = "02_tumorGene_immuneSig_correlation/data/immune_signatures/originalSignaturesGenes/"
    # savepath = "r1_drug_immuneSig_correlation_analysis/"

    library(readxl)
    library(dplyr)

    dir.create("06_results_summaryandplotting/immuneSigGene_ITS_found_lincs_20220414")
    dir.create(paste0(savepath, "/06_results_summaryandplotting/immuneSigGene_ITS_found_lincs_20220414/",purity_method, "_spearman_r",r,"_ITSproof_",ITSproof, "/"))
    dir.create(paste0(savepath, "/06_results_summaryandplotting/ITS_function_analysis_20220414/"))
    dir.create(paste0(savepath, "/06_results_summaryandplotting/ITS_function_analysis_20220414/", purity_method, "_spearman_r_", r, "_ITSproof_",ITSproof,"_", cluster_method, "_", similarity_method,"_genevote_", genevote,"/"))

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

    # cancer oncogenes ---------------------------------------------------------
    canceroncogene <- read.csv("06_results_summaryandplotting/data/canceroncegenes/oncoKB_cancerGeneList.tsv", sep = "\t", header = T)
    oncogene <- canceroncogene[canceroncogene$Is.Oncogene == "Yes",]$Hugo.Symbol
    suppressorgene <- canceroncogene[canceroncogene$Is.Tumor.Suppressor.Gene == "Yes",]$Hugo.Symbol


    # immune modulatory genes --------------------------------------------------

    load("06_results_summaryandplotting/data/immune_modulator/GSVA_ICG.Rdata")
    # GSVA_ICG

    immuneModulator_PMID29628290_1 <- as.data.frame(read_excel("06_results_summaryandplotting/data/immune_modulator/PMID29628290.xlsx", sheet = 1))
    inhibitory <- unique(immuneModulator_PMID29628290_1[immuneModulator_PMID29628290_1$`Immune Checkpoint` == "Inhibitory",]$`HGNC Symbol`)
    inhibitory <- inhibitory[!is.na(inhibitory)]
    stimulatory <- unique(immuneModulator_PMID29628290_1[immuneModulator_PMID29628290_1$`Immune Checkpoint` != "Inhibitory",]$`HGNC Symbol`)
    stimulatory <- stimulatory[!is.na(stimulatory)]


    immuneModulator_PMID29628290_2 <- as.data.frame(read_excel("06_results_summaryandplotting/data/immune_modulator/PMID29628290.xlsx", sheet = 2))
    co_stim_inhib <- unique(immuneModulator_PMID29628290_2[immuneModulator_PMID29628290_2$Reason == "co-stim/inhib list", ]$Gene)

    immuneModulator_PMID34980132_1 <- as.data.frame(read_excel("06_results_summaryandplotting/data/immune_modulator/PMID34980132/ICB enhancer and suppressor genes.xlsx", sheet = 1))
    icb_enhancer <- unique(immuneModulator_PMID34980132_1[immuneModulator_PMID34980132_1$Type == "ICB enhancer gene", ]$Symbol)
    icb_suppressor <- unique(immuneModulator_PMID34980132_1[immuneModulator_PMID34980132_1$Type == "ICB suppressors gene", ]$Symbol)


    immuneModulator_PMID34980132_2 <- as.data.frame(read_excel("06_results_summaryandplotting/data/immune_modulator/PMID34980132/positive and negative regulators.xlsx", sheet = 1))
    regulator_p <- unique(immuneModulator_PMID34980132_2[immuneModulator_PMID34980132_2$Type == "Positive regulators", ]$Symbol)
    regulator_n <- unique(immuneModulator_PMID34980132_2[immuneModulator_PMID34980132_2$Type == "Negative regulators", ]$Symbol)

    immuneModulator_PMID34980132_3 <- as.data.frame(read_excel("06_results_summaryandplotting/data/immune_modulator/PMID34980132/sensitizer_resistor_genes.xlsx", sheet = 1))
    Sensitizer <- unique(immuneModulator_PMID34980132_3[immuneModulator_PMID34980132_3$Type == "Sensitizer", ]$Symbol)
    Resistor <- unique(immuneModulator_PMID34980132_3[immuneModulator_PMID34980132_3$Type == "Resistor", ]$Symbol)

    immuneModulator_PMID34980132_4 <- as.data.frame(read_excel("06_results_summaryandplotting/data/immune_modulator/PMID34980132/Tumor immunity related OGs TSGs.xlsx", sheet = 1))
    OG <- unique(immuneModulator_PMID34980132_4[immuneModulator_PMID34980132_4$Type == "Tumor immunity-related OGs", ]$Symbol)
    TSG <- unique(immuneModulator_PMID34980132_4[immuneModulator_PMID34980132_4$Type == "Tumor immunity-related TSGs", ]$Symbol)


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
    # ITSpath = "03_drug_immuneSig_enrichment/data/gene_immune_sig_TUMERIC_spearman/for_GSVA/pcor/"
    # savepath = "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/"

    # dir.create("06_results_summaryandplotting/ITSgene_found_lincs")
    # dir.create(paste0("06_results_summaryandplotting/ITSgene_found_lincs/", purity_method))

    # ITSpath = paste0("03_drug_immuneSig_enrichment/data/gene_immune_sig_TUMERIC_spearman_r",r,"_ITSproof_",ITSproof,"/for_GSVA/pcor/")
    ITSpath = paste0("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/02_tumorGene_immuneSig_correlation/results_selected/ITS_oe_original_spearman_",purity_method,"_r",r,"/for_GSVA/pcor/")
    load(paste0(ITSpath, "pcor_spearman_",cancer,"_", tumor_subtype, "_positive_200.Rdata"))
    load(paste0(ITSpath, "pcor_spearman_",cancer,"_", tumor_subtype, "_negative_200.Rdata"))



    immunesigGS_selected = immunesigGS[intersect(  names(immunesigGS), names(genelist_p))]
    immunesigGS_filtered <- lapply(immunesigGS_selected, function(x){
      intersect(x, geneExprInCCLE_filtered$gene_name)
    })
    names(immunesigGS_filtered) = names(immunesigGS_selected)

    genelist_p_ori = genelist_p[names(immunesigGS_selected)]
    genelist_n_ori = genelist_n[names(immunesigGS_selected)]

    genelist_p <- lapply(genelist_p_ori, function(x){
      intersect(x, geneExprInCCLE_filtered$gene_name)
    })
    names(genelist_p) = names(genelist_p_ori)

    genelist_n <- lapply(genelist_n_ori, function(x){
      intersect(x, geneExprInCCLE_filtered$gene_name)
    })
    names(genelist_n) = names(genelist_n_ori)


    res = data.frame(immune_sig = names(genelist_p),
                     gene_in_OriginalSig = 0,
                     gene_in_OriginalSig_expressedinCCLE = 0, 
                     gene_in_OriginalSig_lincs_overlap = 0,
                     gene_in_OriginalSig_lincs_all_overlap = 0,

                     gene_in_OriginalSig_drugtarget_overlap = 0,
                     gene_in_OriginalSig_drugtargetOASIS_overlap = 0,
                     gene_in_OriginalSig_oncogene_overlap = 0,
                     gene_in_OriginalSig_suppressorgene_overlap = 0,
                     gene_in_OriginalSig_inhibitory_overlap = 0,
                     gene_in_OriginalSig_stimulatory_overlap = 0,
                     gene_in_OriginalSig_co_stim_inhib_overlap = 0,
                     gene_in_OriginalSig_icb_enhancer_overlap = 0,
                     gene_in_OriginalSig_icb_suppressor_overlap = 0,
                     gene_in_OriginalSig_regulator_p_overlap = 0,
                     gene_in_OriginalSig_regulator_n_overlap = 0,
                     gene_in_OriginalSig_Sensitizer_overlap = 0,
                     gene_in_OriginalSig_Resistor_overlap = 0,
                     gene_in_OriginalSig_OG_overlap = 0,
                     gene_in_OriginalSig_TSG_overlap = 0,

                     gene_in_ITS_p = 0,
                     gene_in_ITS_p_expressedinCCLE = 0,
                     gene_in_ITS_p_lincs_overlap = 0,
                     gene_in_ITS_p_lincs_all_overlap = 0,
                     gene_in_OriginalSig_and_ITS_p = 0,

                     gene_in_ITS_p_drugtarget_overlap = 0,
                     gene_in_ITS_p_drugtargetOASIS_overlap = 0,
                     gene_in_ITS_p_oncogene_overlap = 0,
                     gene_in_ITS_p_suppressorgene_overlap = 0,
                     gene_in_ITS_p_inhibitory_overlap = 0,
                     gene_in_ITS_p_stimulatory_overlap = 0,
                     gene_in_ITS_p_co_stim_inhib_overlap = 0,
                     gene_in_ITS_p_icb_enhancer_overlap = 0,
                     gene_in_ITS_p_icb_suppressor_overlap = 0,
                     gene_in_ITS_p_regulator_p_overlap = 0,
                     gene_in_ITS_p_regulator_n_overlap = 0,
                     gene_in_ITS_p_Sensitizer_overlap = 0,
                     gene_in_ITS_p_Resistor_overlap = 0,
                     gene_in_ITS_p_OG_overlap = 0,
                     gene_in_ITS_p_TSG_overlap = 0,


                     gene_in_ITS_n = 0,
                     gene_in_ITS_n_expressedinCCLE = 0,
                     gene_in_ITS_n_lincs_overlap = 0,
                     gene_in_ITS_n_lincs_all_overlap = 0,
                     gene_in_OriginalSig_and_ITS_n = 0,

                     gene_in_ITS_n_drugtarget_overlap = 0,
                     gene_in_ITS_n_drugtargetOASIS_overlap = 0,
                     gene_in_ITS_n_drugtarget_overlap = 0,
                     gene_in_ITS_n_drugtargetOASIS_overlap = 0,
                     gene_in_ITS_n_oncogene_overlap = 0,
                     gene_in_ITS_n_suppressorgene_overlap = 0,
                     gene_in_ITS_n_inhibitory_overlap = 0,
                     gene_in_ITS_n_stimulatory_overlap = 0,
                     gene_in_ITS_n_co_stim_inhib_overlap = 0,
                     gene_in_ITS_n_icb_enhancer_overlap = 0,
                     gene_in_ITS_n_icb_suppressor_overlap = 0,
                     gene_in_ITS_n_regulator_n_overlap = 0,
                     gene_in_ITS_n_regulator_n_overlap = 0,
                     gene_in_ITS_n_Sensitizer_overlap = 0,
                     gene_in_ITS_n_Resistor_overlap = 0,
                     gene_in_ITS_n_OG_overlap = 0,
                     gene_in_ITS_n_TSG_overlap = 0,

                     
                     gene_in_CCLE = 0, 
                     gene_in_CCLEfiltered = 0,
                     gene_in_lincs = 0,
                     gene_in_lincs_all = 0
                     )


    for(i in seq(length(immunesigGS_selected))){

        res[which(res$immune_sig == names(immunesigGS_selected)[i]),]$gene_in_OriginalSig = length(immunesigGS_selected[[i]])
        res[which(res$immune_sig == names(immunesigGS_selected)[i]),]$gene_in_OriginalSig_expressedinCCLE = length(immunesigGS_filtered[[i]])
        res[which(res$immune_sig == names(immunesigGS_selected)[i]),]$gene_in_OriginalSig_lincs_overlap = length(intersect(immunesigGS_filtered[[i]], lincs_gene[lincs_gene$pr_is_lm ==1,]$pr_gene_symbol))
        res[which(res$immune_sig == names(immunesigGS_selected)[i]),]$gene_in_OriginalSig_lincs_all_overlap = length(intersect(immunesigGS_filtered[[i]], lincs_gene$pr_gene_symbol))


        res[which(res$immune_sig == names(immunesigGS_selected)[i]),]$gene_in_OriginalSig_drugtarget_overlap = length(intersect(immunesigGS_filtered[[i]], drugtarget$gene_name))
        res[which(res$immune_sig == names(immunesigGS_selected)[i]),]$gene_in_OriginalSig_drugtargetOASIS_overlap = length(intersect(immunesigGS_filtered[[i]], drugtarget_OASIS$gene_name))
        res[which(res$immune_sig == names(immunesigGS_selected)[i]),]$gene_in_OriginalSig_oncogene_overlap = length(intersect(immunesigGS_filtered[[i]], oncogene))
        res[which(res$immune_sig == names(immunesigGS_selected)[i]),]$gene_in_OriginalSig_suppressorgene_overlap = length(intersect(immunesigGS_filtered[[i]], suppressorgene))
        res[which(res$immune_sig == names(immunesigGS_selected)[i]),]$gene_in_OriginalSig_inhibitory_overlap = length(intersect(immunesigGS_filtered[[i]], inhibitory))
        res[which(res$immune_sig == names(immunesigGS_selected)[i]),]$gene_in_OriginalSig_stimulatory_overlap = length(intersect(immunesigGS_filtered[[i]], stimulatory))
        res[which(res$immune_sig == names(immunesigGS_selected)[i]),]$gene_in_OriginalSig_co_stim_inhib_overlap = length(intersect(immunesigGS_filtered[[i]], co_stim_inhib))
        res[which(res$immune_sig == names(immunesigGS_selected)[i]),]$gene_in_OriginalSig_icb_enhancer_overlap = length(intersect(immunesigGS_filtered[[i]], icb_enhancer))
        res[which(res$immune_sig == names(immunesigGS_selected)[i]),]$gene_in_OriginalSig_icb_suppressor_overlap = length(intersect(immunesigGS_filtered[[i]], icb_suppressor))
        res[which(res$immune_sig == names(immunesigGS_selected)[i]),]$gene_in_OriginalSig_regulator_p_overlap = length(intersect(immunesigGS_filtered[[i]], regulator_p))
        res[which(res$immune_sig == names(immunesigGS_selected)[i]),]$gene_in_OriginalSig_regulator_n_overlap = length(intersect(immunesigGS_filtered[[i]], regulator_n))
        res[which(res$immune_sig == names(immunesigGS_selected)[i]),]$gene_in_OriginalSig_Sensitizer_overlap = length(intersect(immunesigGS_filtered[[i]], Sensitizer))
        res[which(res$immune_sig == names(immunesigGS_selected)[i]),]$gene_in_OriginalSig_Resistor_overlap = length(intersect(immunesigGS_filtered[[i]], Resistor))
        res[which(res$immune_sig == names(immunesigGS_selected)[i]),]$gene_in_OriginalSig_OG_overlap = length(intersect(immunesigGS_filtered[[i]], OG))
        res[which(res$immune_sig == names(immunesigGS_selected)[i]),]$gene_in_OriginalSig_TSG_overlap = length(intersect(immunesigGS_filtered[[i]], TSG))


        res[which(res$immune_sig == names(genelist_p)[i]),]$gene_in_ITS_p = length(genelist_p[[i]])
        res[which(res$immune_sig == names(genelist_p)[i]),]$gene_in_ITS_p_expressedinCCLE = length(intersect(genelist_p[[i]], geneExprInCCLE_filtered$gene_name))

        res[which(res$immune_sig == names(genelist_p)[i]),]$gene_in_ITS_p_lincs_overlap = length(intersect(genelist_p[[i]], lincs_gene[lincs_gene$pr_is_lm ==1,]$pr_gene_symbol))
        res[which(res$immune_sig == names(genelist_p)[i]),]$gene_in_ITS_p_lincs_all_overlap = length(intersect(genelist_p[[i]], lincs_gene$pr_gene_symbol))
        res[which(res$immune_sig == names(genelist_p)[i]),]$gene_in_OriginalSig_and_ITS_p = length(intersect(genelist_p[[i]], immunesigGS_filtered[[i]]))

        res[which(res$immune_sig == names(genelist_p)[i]),]$gene_in_ITS_p_drugtarget_overlap = length(intersect(genelist_p[[i]], drugtarget$gene_name))
        res[which(res$immune_sig == names(genelist_p)[i]),]$gene_in_ITS_p_drugtargetOASIS_overlap = length(intersect(genelist_p[[i]], drugtarget_OASIS$gene_name))

        res[which(res$immune_sig == names(genelist_p)[i]),]$gene_in_ITS_p_oncogene_overlap = length(intersect(genelist_p[[i]], oncogene))
        res[which(res$immune_sig == names(genelist_p)[i]),]$gene_in_ITS_p_suppressorgene_overlap = length(intersect(genelist_p[[i]], suppressorgene))
        res[which(res$immune_sig == names(genelist_p)[i]),]$gene_in_ITS_p_inhibitory_overlap = length(intersect(genelist_p[[i]], inhibitory))
        res[which(res$immune_sig == names(genelist_p)[i]),]$gene_in_ITS_p_stimulatory_overlap = length(intersect(genelist_p[[i]], stimulatory))
        res[which(res$immune_sig == names(genelist_p)[i]),]$gene_in_ITS_p_co_stim_inhib_overlap = length(intersect(genelist_p[[i]], co_stim_inhib))
        res[which(res$immune_sig == names(genelist_p)[i]),]$gene_in_ITS_p_icb_enhancer_overlap = length(intersect(genelist_p[[i]], icb_enhancer))
        res[which(res$immune_sig == names(genelist_p)[i]),]$gene_in_ITS_p_icb_suppressor_overlap = length(intersect(genelist_p[[i]], icb_suppressor))
        res[which(res$immune_sig == names(genelist_p)[i]),]$gene_in_ITS_p_regulator_p_overlap = length(intersect(genelist_p[[i]], regulator_p))
        res[which(res$immune_sig == names(genelist_p)[i]),]$gene_in_ITS_p_regulator_n_overlap = length(intersect(genelist_p[[i]], regulator_n))
        res[which(res$immune_sig == names(genelist_p)[i]),]$gene_in_ITS_p_Sensitizer_overlap = length(intersect(genelist_p[[i]], Sensitizer))
        res[which(res$immune_sig == names(genelist_p)[i]),]$gene_in_ITS_p_Resistor_overlap = length(intersect(genelist_p[[i]], Resistor))
        res[which(res$immune_sig == names(genelist_p)[i]),]$gene_in_ITS_p_OG_overlap = length(intersect(genelist_p[[i]], OG))
        res[which(res$immune_sig == names(genelist_p)[i]),]$gene_in_ITS_p_TSG_overlap = length(intersect(genelist_p[[i]], TSG))






        res[which(res$immune_sig == names(genelist_n)[i]),]$gene_in_ITS_n = length(genelist_n[[i]])
        res[which(res$immune_sig == names(genelist_p)[i]),]$gene_in_ITS_n_expressedinCCLE = length(intersect(genelist_n[[i]], geneExprInCCLE_filtered$gene_name))

        res[which(res$immune_sig == names(genelist_n)[i]),]$gene_in_ITS_n_lincs_overlap = length(intersect(genelist_n[[i]], lincs_gene[lincs_gene$pr_is_lm ==1,]$pr_gene_symbol))
        res[which(res$immune_sig == names(genelist_n)[i]),]$gene_in_ITS_n_lincs_all_overlap = length(intersect(genelist_n[[i]], lincs_gene$pr_gene_symbol))
        res[which(res$immune_sig == names(genelist_n)[i]),]$gene_in_OriginalSig_and_ITS_n = length(intersect(genelist_n[[i]], immunesigGS_filtered[[i]]))


        res[which(res$immune_sig == names(genelist_n)[i]),]$gene_in_ITS_n_drugtarget_overlap = length(intersect(genelist_n[[i]], drugtarget$gene_name))
        res[which(res$immune_sig == names(genelist_n)[i]),]$gene_in_ITS_n_drugtargetOASIS_overlap = length(intersect(genelist_n[[i]], drugtarget_OASIS$gene_name))
        res[which(res$immune_sig == names(genelist_n)[i]),]$gene_in_ITS_n_oncogene_overlap = length(intersect(genelist_n[[i]], oncogene))
        res[which(res$immune_sig == names(genelist_n)[i]),]$gene_in_ITS_n_suppressorgene_overlap = length(intersect(genelist_n[[i]], suppressorgene))
        res[which(res$immune_sig == names(genelist_n)[i]),]$gene_in_ITS_n_inhibitory_overlap = length(intersect(genelist_n[[i]], inhibitory))
        res[which(res$immune_sig == names(genelist_n)[i]),]$gene_in_ITS_n_stimulatory_overlap = length(intersect(genelist_n[[i]], stimulatory))
        res[which(res$immune_sig == names(genelist_n)[i]),]$gene_in_ITS_n_co_stim_inhib_overlap = length(intersect(genelist_n[[i]], co_stim_inhib))
        res[which(res$immune_sig == names(genelist_n)[i]),]$gene_in_ITS_n_icb_enhancer_overlap = length(intersect(genelist_n[[i]], icb_enhancer))
        res[which(res$immune_sig == names(genelist_n)[i]),]$gene_in_ITS_n_icb_suppressor_overlap = length(intersect(genelist_n[[i]], icb_suppressor))
        res[which(res$immune_sig == names(genelist_n)[i]),]$gene_in_ITS_n_regulator_n_overlap = length(intersect(genelist_n[[i]], regulator_p))
        res[which(res$immune_sig == names(genelist_n)[i]),]$gene_in_ITS_n_regulator_n_overlap = length(intersect(genelist_n[[i]], regulator_n))
        res[which(res$immune_sig == names(genelist_n)[i]),]$gene_in_ITS_n_Sensitizer_overlap = length(intersect(genelist_n[[i]], Sensitizer))
        res[which(res$immune_sig == names(genelist_n)[i]),]$gene_in_ITS_n_Resistor_overlap = length(intersect(genelist_n[[i]], Resistor))
        res[which(res$immune_sig == names(genelist_n)[i]),]$gene_in_ITS_n_OG_overlap = length(intersect(genelist_n[[i]], OG))
        res[which(res$immune_sig == names(genelist_n)[i]),]$gene_in_ITS_n_TSG_overlap = length(intersect(genelist_n[[i]], TSG))



        res[which(res$immune_sig == names(immunesigGS_selected)[i]),]$gene_in_CCLE = nrow(geneExprInCCLE)
        res[which(res$immune_sig == names(immunesigGS_selected)[i]),]$gene_in_CCLEfiltered = nrow(geneExprInCCLE_filtered)
        res[which(res$immune_sig == names(immunesigGS_selected)[i]),]$gene_in_lincs = nrow(lincs_gene[lincs_gene$pr_is_lm ==1,])
        res[which(res$immune_sig == names(immunesigGS_selected)[i]),]$gene_in_lincs_all = nrow(lincs_gene)
    }


    immunesigInfo <- as.data.frame(read_excel( "06_results_summaryandplotting/data/immene_cell_hieracy/immune_sig_hieracy_new.xlsx", sheet = 1))
    immunesigInfo <- immunesigInfo[c(1, 4:9)]
    immunesigInfo$immune_sig <- tolower(immunesigInfo$immune_sig)
    immunesigInfo <- immunesig_name_unify(immunesigInfo)

    res <- inner_join(immunesigInfo, res)

    library(openxlsx)
    write.xlsx(list("Sheet1" = res), file = paste0(savepath, "/06_results_summaryandplotting/immuneSigGene_ITS_found_lincs_20220414/",purity_method, "_spearman_r",r,"_ITSproof_",ITSproof, "/", cancer,"_",tumor_subtype,".xlsx"))



  #   ITSpath = paste0("03_drug_immuneSig_enrichment/data/ITS_", cluster_method, "_", similarity_method, "_",purity_method, "_spearman_r_",r,"_pan_metaP_0.05_genevote_", genevote, "/for_GSVA/pcor/")
  #   load(paste0(ITSpath, "spearman_",cancer,"_", tumor_subtype, "_positive_200.Rdata"))
  #   load(paste0(ITSpath, "spearman_",cancer,"_", tumor_subtype, "_negative_200.Rdata"))

  genelist_p_merged = genelist_p
  genelist_n_merged = genelist_n

    res_merged = data.frame(immune_sig = names(genelist_p_merged),

                            gene_in_ITS_p = 0,
                            gene_in_ITS_p_lincs_overlap = 0,
                            gene_in_ITS_p_lincs_all_overlap = 0,
                            gene_in_OriginalSig_and_ITS_p = 0,

                            gene_in_ITS_p_drugtarget_overlap = 0,
                            gene_in_ITS_p_drugtargetOASIS_overlap = 0,
                            gene_in_ITS_p_oncogene_overlap = 0,
                            gene_in_ITS_p_suppressorgene_overlap = 0,
                            gene_in_ITS_p_inhibitory_overlap = 0,
                            gene_in_ITS_p_stimulatory_overlap = 0,
                            gene_in_ITS_p_co_stim_inhib_overlap = 0,
                            gene_in_ITS_p_icb_enhancer_overlap = 0,
                            gene_in_ITS_p_icb_suppressor_overlap = 0,
                            gene_in_ITS_p_regulator_p_overlap = 0,
                            gene_in_ITS_p_regulator_n_overlap = 0,
                            gene_in_ITS_p_Sensitizer_overlap = 0,
                            gene_in_ITS_p_Resistor_overlap = 0,
                            gene_in_ITS_p_OG_overlap = 0,
                            gene_in_ITS_p_TSG_overlap = 0,


                            gene_in_ITS_n = 0,
                            gene_in_ITS_n_lincs_overlap = 0,
                            gene_in_ITS_n_lincs_all_overlap = 0,
                            gene_in_OriginalSig_and_ITS_n = 0,

                            gene_in_ITS_n_drugtarget_overlap = 0,
                            gene_in_ITS_n_drugtargetOASIS_overlap = 0,
                            gene_in_ITS_n_drugtarget_overlap = 0,
                            gene_in_ITS_n_drugtargetOASIS_overlap = 0,
                            gene_in_ITS_n_oncogene_overlap = 0,
                            gene_in_ITS_n_suppressorgene_overlap = 0,
                            gene_in_ITS_n_inhibitory_overlap = 0,
                            gene_in_ITS_n_stimulatory_overlap = 0,
                            gene_in_ITS_n_co_stim_inhib_overlap = 0,
                            gene_in_ITS_n_icb_enhancer_overlap = 0,
                            gene_in_ITS_n_icb_suppressor_overlap = 0,
                            gene_in_ITS_n_regulator_n_overlap = 0,
                            gene_in_ITS_n_regulator_n_overlap = 0,
                            gene_in_ITS_n_Sensitizer_overlap = 0,
                            gene_in_ITS_n_Resistor_overlap = 0,
                            gene_in_ITS_n_OG_overlap = 0,
                            gene_in_ITS_n_TSG_overlap = 0
                            )


    for(i in seq(length(genelist_p_merged))){

        res_merged[which(res_merged$immune_sig == names(genelist_p_merged)[i]),]$gene_in_ITS_p = length(genelist_p_merged[[i]])

        res_merged[which(res_merged$immune_sig == names(genelist_p_merged)[i]),]$gene_in_ITS_p_lincs_overlap = length(intersect(genelist_p_merged[[i]], lincs_gene[lincs_gene$pr_is_lm ==1,]$pr_gene_symbol))
        res_merged[which(res_merged$immune_sig == names(genelist_p_merged)[i]),]$gene_in_ITS_p_lincs_all_overlap = length(intersect(genelist_p_merged[[i]], lincs_gene$pr_gene_symbol))
        res_merged[which(res_merged$immune_sig == names(genelist_p_merged)[i]),]$gene_in_OriginalSig_and_ITS_p = length(intersect(genelist_p_merged[[i]], immunesigGS_filtered[[i]]))

        res_merged[which(res_merged$immune_sig == names(genelist_p_merged)[i]),]$gene_in_ITS_p_drugtarget_overlap = length(intersect(genelist_p_merged[[i]], drugtarget$gene_name))
        res_merged[which(res_merged$immune_sig == names(genelist_p_merged)[i]),]$gene_in_ITS_p_drugtargetOASIS_overlap = length(intersect(genelist_p_merged[[i]], drugtarget_OASIS$gene_name))

        res_merged[which(res_merged$immune_sig == names(genelist_p_merged)[i]),]$gene_in_ITS_p_oncogene_overlap = length(intersect(genelist_p_merged[[i]], oncogene))
        res_merged[which(res_merged$immune_sig == names(genelist_p_merged)[i]),]$gene_in_ITS_p_suppressorgene_overlap = length(intersect(genelist_p_merged[[i]], suppressorgene))
        res_merged[which(res_merged$immune_sig == names(genelist_p_merged)[i]),]$gene_in_ITS_p_inhibitory_overlap = length(intersect(genelist_p_merged[[i]], inhibitory))
        res_merged[which(res_merged$immune_sig == names(genelist_p_merged)[i]),]$gene_in_ITS_p_stimulatory_overlap = length(intersect(genelist_p_merged[[i]], stimulatory))
        res_merged[which(res_merged$immune_sig == names(genelist_p_merged)[i]),]$gene_in_ITS_p_co_stim_inhib_overlap = length(intersect(genelist_p_merged[[i]], co_stim_inhib))
        res_merged[which(res_merged$immune_sig == names(genelist_p_merged)[i]),]$gene_in_ITS_p_icb_enhancer_overlap = length(intersect(genelist_p_merged[[i]], icb_enhancer))
        res_merged[which(res_merged$immune_sig == names(genelist_p_merged)[i]),]$gene_in_ITS_p_icb_suppressor_overlap = length(intersect(genelist_p_merged[[i]], icb_suppressor))
        res_merged[which(res_merged$immune_sig == names(genelist_p_merged)[i]),]$gene_in_ITS_p_regulator_p_overlap = length(intersect(genelist_p_merged[[i]], regulator_p))
        res_merged[which(res_merged$immune_sig == names(genelist_p_merged)[i]),]$gene_in_ITS_p_regulator_n_overlap = length(intersect(genelist_p_merged[[i]], regulator_n))
        res_merged[which(res_merged$immune_sig == names(genelist_p_merged)[i]),]$gene_in_ITS_p_Sensitizer_overlap = length(intersect(genelist_p_merged[[i]], Sensitizer))
        res_merged[which(res_merged$immune_sig == names(genelist_p_merged)[i]),]$gene_in_ITS_p_Resistor_overlap = length(intersect(genelist_p_merged[[i]], Resistor))
        res_merged[which(res_merged$immune_sig == names(genelist_p_merged)[i]),]$gene_in_ITS_p_OG_overlap = length(intersect(genelist_p_merged[[i]], OG))
        res_merged[which(res_merged$immune_sig == names(genelist_p_merged)[i]),]$gene_in_ITS_p_TSG_overlap = length(intersect(genelist_p_merged[[i]], TSG))



        res_merged[which(res_merged$immune_sig == names(genelist_n_merged)[i]),]$gene_in_ITS_n = length(genelist_n_merged[[i]])

        res_merged[which(res_merged$immune_sig == names(genelist_n_merged)[i]),]$gene_in_ITS_n_lincs_overlap = length(intersect(genelist_n_merged[[i]], lincs_gene[lincs_gene$pr_is_lm ==1,]$pr_gene_symbol))
        res_merged[which(res_merged$immune_sig == names(genelist_n_merged)[i]),]$gene_in_ITS_n_lincs_all_overlap = length(intersect(genelist_n_merged[[i]], lincs_gene$pr_gene_symbol))
        res_merged[which(res_merged$immune_sig == names(genelist_n_merged)[i]),]$gene_in_OriginalSig_and_ITS_n = length(intersect(genelist_n_merged[[i]], immunesigGS_filtered[[i]]))


        res_merged[which(res_merged$immune_sig == names(genelist_n_merged)[i]),]$gene_in_ITS_n_drugtarget_overlap = length(intersect(genelist_n_merged[[i]], drugtarget$gene_name))
        res_merged[which(res_merged$immune_sig == names(genelist_n_merged)[i]),]$gene_in_ITS_n_drugtargetOASIS_overlap = length(intersect(genelist_n_merged[[i]], drugtarget_OASIS$gene_name))
        res_merged[which(res_merged$immune_sig == names(genelist_n_merged)[i]),]$gene_in_ITS_n_oncogene_overlap = length(intersect(genelist_n_merged[[i]], oncogene))
        res_merged[which(res_merged$immune_sig == names(genelist_n_merged)[i]),]$gene_in_ITS_n_suppressorgene_overlap = length(intersect(genelist_n_merged[[i]], suppressorgene))
        res_merged[which(res_merged$immune_sig == names(genelist_n_merged)[i]),]$gene_in_ITS_n_inhibitory_overlap = length(intersect(genelist_n_merged[[i]], inhibitory))
        res_merged[which(res_merged$immune_sig == names(genelist_n_merged)[i]),]$gene_in_ITS_n_stimulatory_overlap = length(intersect(genelist_n_merged[[i]], stimulatory))
        res_merged[which(res_merged$immune_sig == names(genelist_n_merged)[i]),]$gene_in_ITS_n_co_stim_inhib_overlap = length(intersect(genelist_n_merged[[i]], co_stim_inhib))
        res_merged[which(res_merged$immune_sig == names(genelist_n_merged)[i]),]$gene_in_ITS_n_icb_enhancer_overlap = length(intersect(genelist_n_merged[[i]], icb_enhancer))
        res_merged[which(res_merged$immune_sig == names(genelist_n_merged)[i]),]$gene_in_ITS_n_icb_suppressor_overlap = length(intersect(genelist_n_merged[[i]], icb_suppressor))
        res_merged[which(res_merged$immune_sig == names(genelist_n_merged)[i]),]$gene_in_ITS_n_regulator_n_overlap = length(intersect(genelist_n_merged[[i]], regulator_p))
        res_merged[which(res_merged$immune_sig == names(genelist_n_merged)[i]),]$gene_in_ITS_n_regulator_n_overlap = length(intersect(genelist_n_merged[[i]], regulator_n))
        res_merged[which(res_merged$immune_sig == names(genelist_n_merged)[i]),]$gene_in_ITS_n_Sensitizer_overlap = length(intersect(genelist_n_merged[[i]], Sensitizer))
        res_merged[which(res_merged$immune_sig == names(genelist_n_merged)[i]),]$gene_in_ITS_n_Resistor_overlap = length(intersect(genelist_n_merged[[i]], Resistor))
        res_merged[which(res_merged$immune_sig == names(genelist_n_merged)[i]),]$gene_in_ITS_n_OG_overlap = length(intersect(genelist_n_merged[[i]], OG))
        res_merged[which(res_merged$immune_sig == names(genelist_n_merged)[i]),]$gene_in_ITS_n_TSG_overlap = length(intersect(genelist_n_merged[[i]], TSG))



    }


    write.xlsx(res_merged, file = paste0(savepath, "/06_results_summaryandplotting/ITS_function_analysis_20220414/", purity_method, "_spearman_r_", r, "_ITSproof_",ITSproof,"_", cluster_method, "_", similarity_method, "_genevote_", genevote,"/", cancer,"_",tumor_subtype,"_merged.xlsx"))

  # return(res)
}

# apply(res[-1], 2, summary)
# apply(res_merged[-1], 2, summary)


# # length(genelist_p[["cd4..memory.t.cells_xcell"]])
# # length(immunesigGS_selected[["cd4..memory.t.cells_xcell"]])
# # length(intersect(genelist_p[["cd4..memory.t.cells_xcell"]],immunesigGS_selected[["cd4..memory.t.cells_xcell"]] ))
# # intersect(genelist_p[["cd4..memory.t.cells_xcell"]],immunesigGS_selected[["cd4..memory.t.cells_xcell"]] )

# plot(density(res$Sheet1$gene_in_OriginalSig_drugtarget_overlap))
# lines(density(res$Sheet1$gene_in_ITS_p_drugtarget_overlap), col = "red")
# lines(density(res$Sheet1$gene_in_ITS_n_drugtarget_overlap), col = "green")

# # plot(density(res$gene_in_OriginalSig_drugtargetOASIS_overlap))
# # lines(density(res$gene_in_ITS_p_drugtargetOASIS_overlap), col = "red")
# # lines(density(res$gene_in_ITS_n_drugtargetOASIS_overlap), col = "green")
# # devtools::install_github("clauswilke/ggjoy")
# library(ggjoy)
# library(ggplot2)

# dtf = res[c("immune_sig", "gene_in_OriginalSig_expressedinCCLE", "gene_in_ITS_p_lincs_all_overlap")]
# dtf = melt(dtf)

# ggplot(dtf, aes(x = value, y = as.factor(variable), fill = variable)) +
#   geom_joy(scale = 2,alpha = .5,rel_min_height = 0.01) + 
#   theme_joy() 

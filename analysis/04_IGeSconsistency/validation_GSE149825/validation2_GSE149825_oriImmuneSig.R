setwd("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/")
options(stringsAsFactors = F)
# Checking and loading R packages ----------------------------------------------
if (!require(readxl)) {install.packages('readxl')}
if (!require(data.table)) {install.packages('data.table')}
if (!require(dplyr)) {install.packages('dplyr')}
# if (!require(TCGAbiolinks)) {BiocManager::install("TCGAbiolinks") }
if (!require(SummarizedExperiment)) {BiocManager::install("SummarizedExperiment")}
if (!require(GSVA)) {BiocManager::install("GSVA")}
if (!require(limma)) {BiocManager::install("limma")}
if (!require(ConsensusClusterPlus)) {BiocManager::install("ConsensusClusterPlus")}
if (!require(immunedeconv)) {devtools::install_github("grst/immunedeconv")}

library(ROCR)
library(reshape2)

# set result saving path -------------------------------------------------------

GEO_ID = "GSE149825"
drugindex = "drug1487" 

result_dir = paste0("result_", GEO_ID, "/original_immuneSig/")
dir.create(result_dir)
cancer = "SKCM"
sample_type = "Metastatic"


# Data preparation -------------------------------------------------------------
# Count to TPM -----------------------------------------------------------------

genename = read.csv(
  paste0(
    "/picb/bigdata/project/FengFYM/p2_2_liver_cancer_specific_immuno_combo_therapy/",
    "03_mouse_data_analysis/data/mousedata_for_immuClustering_icbResponse/mouse2human_count_manually_genename_merge.txt"
  ),
  sep = "\t",
  header = T
)

# Gene name conversion ---------------------------------------------------------
count <- read.csv("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/data/immune_combo_dataset/GSE149825/GSE149825_SMAC_ICB_RNASeq_rawcount.csv", row.names = 1)

gene_length = read.table(paste0("/picb/bigdata/project/FengFYM/s7_annotation_files/All_mm10_gene_len_mus.txt"),
  sep = "\t", header = T)

count$Gene = rownames(count)
count_genelength = inner_join(gene_length, count, by = "Gene")
genelength_kb = count_genelength$Length / 1000
rownames(count_genelength) <- count_genelength$Gene
CountMat = count_genelength[-which(is.element(colnames(count_genelength),
                                   names(gene_length)))]

rpk <- CountMat / genelength_kb
tpm_mat <- as.data.frame(t(t(rpk) / colSums(rpk) * 10 ^ 6))

tpm_mat$MGI_symbol  = rownames(tpm_mat)
tpm_mat = inner_join(genename, tpm_mat, by = "MGI_symbol")
tpm_mat = tpm_mat[!duplicated(tpm_mat$hgnc_symbol),]
rownames(tpm_mat) = tpm_mat$hgnc_symbol
tpm_mat = tpm_mat[-c(1,2,3)]

dir.create(paste0(result_dir, "/processedData_TPM/"))
write.csv(tpm_mat, paste0(result_dir, "/processedData_TPM/tpm_mat.csv"),
          row.names = T, quote = F)






# Original immune signatures analysis ------------------------------------------
workpath <- "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/"
# setwd(workpath)

# Loading relied functions -----------------------------------------------------
# TIL estimation, TIDE, TIP, IPS and Immune Resistance Pragram
source(paste0(workpath, "02_tumorGene_immuneSig_correlation/scripts/functions/02_immuneSig_generation_function.R"))
source(paste0(workpath, "02_tumorGene_immuneSig_correlation/scripts/functions/TILestimator.R"))
source(paste0(workpath, "02_tumorGene_immuneSig_correlation/scripts/functions/predictor_TIDEandTIP.R"))
source(paste0(workpath, "02_tumorGene_immuneSig_correlation/scripts/functions/predictor_IPS.R"))
source(paste0(workpath, "02_tumorGene_immuneSig_correlation/scripts/functions/predictor_icb_ICI_resistance_program.R"))
source(paste0(workpath, "02_tumorGene_immuneSig_correlation/scripts/functions/immuneSig_sig160.R"))
source(paste0(workpath, "02_tumorGene_immuneSig_correlation/scripts/functions/immuneSig_TCGA_caf_immuneCLSssgsea.R"))
# TIGS, AUC and Wilcox test
source(paste0(workpath, "05_immuneSig_ICBresponse/scripts/functions/05_immunotherapy_immuneSig_analysis_function.R"))
source(paste0(workpath, "05_immuneSig_ICBresponse/scripts/functions/predictor_ITS.R"))

source(paste0(workpath, "06_results_summaryandplotting/scripts/immunesig_selection_function.R"))



# generate expr matrix for cibersort webserver ---------------------------------
dir.create(paste0(result_dir, "/processedData_forCIBERSORT/"))
dir.create(paste0(result_dir, "/processedData_forImmunCellAI/"))
dir.create(paste0(result_dir, "/processedData_forTIDE/"))
dir.create(paste0(result_dir, "/cibersortx/"))
dir.create(paste0(result_dir, "/ImmuCellAI/"))
dir.create(paste0(result_dir, "/TIDEweb/"))

expr4cibersort(expr = tpm_mat,
               dataset = "GSE149825",
               savepath = paste0(result_dir, "/processedData_forCIBERSORT/"))
expr4immuncelai(expr = tpm_mat,
          dataset = "GSE149825",
          savepath = paste0(result_dir, "/processedData_forImmunCellAI/"))

expr4tide(expr = tpm_mat,
          dataset = "GSE149825",
          savepath = paste0(result_dir, "/processedData_forTIDE/"))

# Immunesig --------------------------------------------------------------------
sourceCodePath <- '/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/02_tumorGene_immuneSig_correlation/scripts/softwares/'
immuneSig_sig160(data = tpm_mat,
                 cancer = cancer,
                 sample_type = sample_type,
                 sourceCodePath = sourceCodePath,
                 savePath = result_dir)

sig_dir <- "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/02_tumorGene_immuneSig_correlation/data/immune_signatures/"
immuneSig_sig160_othersigs(data = tpm_mat,
                           cancer = cancer,
                           sample_type = sample_type,
                           method = "ssgsea",
                           kcdf =  "Gaussian", 
                           min.sz=1, 
                           max.sz=1000,
                           sig_dir = sig_dir,
                           savePath = result_dir)

sig_dir <- "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/02_tumorGene_immuneSig_correlation/data/immune_signatures/"
TCGA_caf_immuneCLS(data = tpm_mat,
                   cancer = cancer,
                   sample_type = sample_type, 
                   method = "ssgsea",
                   kcdf =  "Gaussian", 
                   min.sz=0, 
                   max.sz=1000,
                   sig_dir = sig_dir,
                   savePath = result_dir)


# ICB Predictor ----------------------------------------------------------------
## TIDE and TIP


model_icb_tide_TIP(mat = tpm_mat,
                   cancername = cancer,
                   sample_type = sample_type,
                   sig_dir = "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r2_drug_sensitize_icb_prediction/",
                   save_path = result_dir)


model_icb_IPS(gene_expression = tpm_mat,
              cancer_type = cancer,
              sample_type = sample_type,
              sig_dir = paste0(workpath, "02_tumorGene_immuneSig_correlation/scripts/softwares/Immunophenogram-master/"),
              save_path = result_dir)

model_icb_ICI_resistance_program(gene_expression = tpm_mat,
                                 cancer_type = cancer,
                                 sample_type = sample_type,
                                 sig_dir= paste0(workpath, "02_tumorGene_immuneSig_correlation/scripts/softwares/"),
                                 save_path = result_dir)

# TIL estimation --------------------------------------------------------------
res_TIL_all1 <- immunecell_estiamtor(exprmat = tpm_mat, cancer = cancer,
  sample_type = sample_type, estimator_tool = "immunedeconv",  method = "all",
  savepath = result_dir)
res_TIL_all2 <- immunecell_estiamtor(exprmat = as.matrix(tpm_mat), cancer = cancer,
  sample_type = sample_type, estimator_tool = "ConsensusTME", method = "gsva",
  savepath = result_dir)

rownames(res_TIL_all2) <- paste0(rownames(res_TIL_all2),"_ConsensusTMEgsva")
res_TIL_all = rbind(do.call(rbind, res_TIL_all1),res_TIL_all2)
dir.create(paste0(result_dir, "TIL_estimation/merged/"))
save(res_TIL_all, file = paste0(result_dir,"TIL_estimation/merged/",  
                                cancer, "_", sample_type, "_TILestimation.Rdata"))



# merge result files -----------------------------------------------------------
immuneSig_sig160_res <- read.csv(paste0(result_dir, "sig160_result/zscore_output/", cancer, "_", sample_type, ".csv"),row.names=1)
immuneSig_sig160_res <- as.data.frame(t(immuneSig_sig160_res))
colnames(immuneSig_sig160_res) <- paste0(colnames(immuneSig_sig160_res), "_sig160")
immuneSig_sig160_res$patient <- rownames(immuneSig_sig160_res)

load(paste0(result_dir, "sig160ssgsea_result/scaled_output/ssgsea_", cancer, "_", sample_type, "_immueSig.Rdata"))
immuneSig_sig160_othersigs_res <- esOut_scaled
immuneSig_sig160_othersigs_res <- as.data.frame(t(immuneSig_sig160_othersigs_res))
immuneSig_sig160_othersigs_res$patient <- rownames(immuneSig_sig160_othersigs_res)

load(paste0(result_dir, "TCGA_caf_immuneCLS_result/scaled_output/ssgsea_", cancer, "_", sample_type, "_immueSig.Rdata"))
TCGA_caf_immuneCLS_res <- esOut_scaled
TCGA_caf_immuneCLS_res <- as.data.frame(t(TCGA_caf_immuneCLS_res))
TCGA_caf_immuneCLS_res$patient <- rownames(TCGA_caf_immuneCLS_res)

TIDE_TIP_res <- read.csv(paste0(result_dir, "TIDE_TIP_result/TIDE_TIP_result_", cancer, "_", sample_type, ".csv"),row.names=1)
TIDE_TIP_res <- TIDE_TIP_res[-which(is.element(colnames(TIDE_TIP_res), c("class", "TIP_class")))]
# TIDE_TIP_res <- TIDE_TIP_res[c(1,3,5,6,9,11,18)]
TIDE_TIP_res <- TIDE_TIP_res[c("patient","Dysfunction","Exclusion",
                         "TIDE", "CD8","IFNG","TIP_signature")]
names(TIDE_TIP_res) <- c("patient","Dysfunction_TIDE","Exclusion_TIDE",
                         "TIDE_TIDE", "CD8_TIDE","IFNG_TIDE","TIP_signature_TIP")

IPS_res <- read.csv(paste0(result_dir, "IPS_result/IPS_result_", cancer, "_", sample_type, ".csv"),row.names=1)
names(IPS_res) <- c("patient","MHC_IPS","EC_IPS","SC_IPS","CP_IPS","AZ_IPS","IPS_IPS")

load(paste0(result_dir, "sig_OE_results/sig_OE_", cancer, "_", sample_type, ".Rdata"))
ICIRP_res <- df$res
names(ICIRP_res) <- paste0("OE_", names(ICIRP_res))
ICIRP_res$patient <- rownames(ICIRP_res)

load(paste0(result_dir,"TIL_estimation/merged/", cancer, "_", sample_type, "_TILestimation.Rdata"))
TIL_res <- as.data.frame(t(res_TIL_all))
TIL_res$patient <- rownames(TIL_res)
names(TIL_res) <- gsub("timer.", "", names(TIL_res))
names(TIL_res) <- gsub("quantiseq.", "", names(TIL_res))
names(TIL_res) <- gsub("mcp_counter.", "", names(TIL_res))
names(TIL_res) <- gsub("epic.", "", names(TIL_res))
names(TIL_res) <- gsub("xcell.", "", names(TIL_res))

TIDEweb_res <- read.csv(paste0(result_dir, "TIDEweb/GSE149825_yestoPD1treatment.csv"), row.names=1)
colnames(TIDEweb_res) <- paste0(colnames(TIDEweb_res),"_TIDEweb")
TIDEweb_res$patient = row.names(TIDEweb_res)
TIDEweb_res <- TIDEweb_res[c("patient", paste0(c("MSI.Expr.Sig","Merck18","MDSC", "CAF", "TAM.M2"),"_TIDEweb"))]

immuncellai_res <- read.csv(paste0(result_dir, "ImmuCellAI/ImmuCellAI_abundance_result.txt"), sep = '\t', row.names=1)
colnames(immuncellai_res) <- paste0(colnames(immuncellai_res),"_immuncellai")
immuncellai_res$patient = row.names(immuncellai_res)

cibersortx_res <- read.csv(paste0(result_dir, "/cibersortx/CIBERSORTx_Results.csv"), row.names=1)
cibersortx_res <- cibersortx_res[-c(23:25)]
colnames(cibersortx_res) <- paste0(colnames(cibersortx_res),"_CIBERSORT")
cibersortx_res$patient = row.names(cibersortx_res)

# cibersortxABS_res <- read.csv(paste0(result_dir, "/cibersortx/CIBERSORTxABS_Results.csv"), row.names=1)
# cibersortxABS_res <- cibersortxABS_res[-c(23:25)]
# colnames(cibersortxABS_res) <- paste0(colnames(cibersortxABS_res),"_cibersortxABS")
# cibersortxABS_res$patient = row.names(cibersortxABS_res)

merged_res <- merge(TIDE_TIP_res,
                merge(IPS_res, 
                merge(ICIRP_res, 
                merge(TIL_res, 
                merge(TIDEweb_res,
                merge(immuncellai_res,
                      cibersortx_res,
                by = "patient"), 
                by = "patient"), 
                by = "patient"), 
                by = "patient"), 
                by = "patient"), 
                by = "patient")

merged_res <- merge(merged_res,
              merge(immuneSig_sig160_res,
              merge(immuneSig_sig160_othersigs_res,
                    TCGA_caf_immuneCLS_res,
                    by = "patient"),
                by = "patient"), 
                by = "patient")


write.csv(merged_res, paste0(result_dir, "/", GEO_ID, "_", cancer, "_", sample_type, "_immuneSig_merged.csv"), quote = F, row.names = F)



immuneSig_boxplot = lapply(as.list(names(merged_res)[-1]), function(x){
    # x = "CD8 T cells_mcp_counter"
    tmp = merged_res[c("patient", x)]
    names(tmp) = c("patient", "score")
    tmp$label = "ICB"
    tmp[grep("Vehicle", tmp$patient),]$label = "vehicle"
    tmp[grep("SMAC_", tmp$patient), ]$label = "SMAC"
    tmp[grep("SMAC.ICB_", tmp$patient), ]$label = "SMAC+ICB"
    tmp$label = factor(tmp$label, c("vehicle", "SMAC", "ICB", "SMAC+ICB"))
    tmp1 = tmp[grep("Vehicle", tmp$patient),]
    tmp2 = tmp[grep("SMAC_", tmp$patient),]   
    pval = wilcox.test(tmp1$score, tmp2$score)$p.value
    p <- ggplot(tmp, aes(x=label, y=score, fill=label)) + 
        geom_boxplot() + 
        geom_jitter(shape=16, position=position_jitter(0.2))+
        ggtitle(paste0(x, "\nSMACvsVehicle p = ", round(pval, 4))) + 
        theme_bw()
    # p
    return(p)
})

pdf(paste0(result_dir, "/", GEO_ID, "_", cancer, "_", sample_type, "_immuneSig_merged.pdf"),
            width = 5,
            height = 5, 
            onefile = TRUE)
for(x in seq(length(immuneSig_boxplot))){
        print(immuneSig_boxplot[[x]])
}
dev.off()


merged_res_colnames = data.frame(names = names(merged_res), immune_sig = tolower(names(merged_res)))
merged_res_colnames = immunesig_name_unify(merged_res_colnames)
names(merged_res) = merged_res_colnames$immune_sig
# Mapping drug DEGs to correlation results -------------------------------------
source(paste0(workpath, "06_results_summaryandplotting/scripts/immunesig_selection_function.R"))

# ----------------------------------------------------------------------------
# drug-induced immune signatures NES
cancer = "SKCM"
tumor_subtype = "Metastatic"
purity_method = "TUMERIC"
dataset=70138
datatype="allgenes"
num_gene = 200



ITSproof=2
immunesig <- immune_sig_filter_v2(auc_threshold = 0.6,
                                    immuSigproof = 2, # how many datasets support the immune sigture to have prediction power
                                    ITSproof = 2, # how many datasets support the immune sigture to have prediction power
                                    workpath = "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/",
                                    savepath  = "06_results_summaryandplotting/data/immune_sig_selected_new/")

immunesig <- immunesig[immunesig$Freq >= ITSproof, ]

merged_res_selected <- merged_res[c("patient", immunesig$immune_sig)]


immuneSig_selected_boxplot = lapply(as.list(names(merged_res_selected)[-1]), function(x){
    # x = "CD8 T cells_mcp_counter"
    tmp = merged_res_selected[c("patient", x)]
    names(tmp) = c("patient", "score")
    tmp$label = "ICB"
    tmp[grep("Vehicle", tmp$patient),]$label = "vehicle"
    tmp[grep("SMAC_", tmp$patient), ]$label = "SMAC"
    tmp[grep("SMAC.ICB_", tmp$patient), ]$label = "SMAC+ICB"
    tmp$label = factor(tmp$label, c("vehicle", "SMAC", "ICB", "SMAC+ICB"))
    tmp1 = tmp[grep("Vehicle", tmp$patient),]
    tmp2 = tmp[grep("SMAC_", tmp$patient),]   
    pval = wilcox.test(tmp1$score, tmp2$score)$p.value
    p <- ggplot(tmp, aes(x=label, y=score, fill=label)) + 
        geom_boxplot() + 
        geom_jitter(shape=16, position=position_jitter(0.2))+
        ggtitle(paste0(x, "\nSMACvsVehicle p = ", round(pval, 4))) + 
        theme_bw()
    # p
    return(p)
})

pdf(paste0(result_dir, "/", GEO_ID, "_", cancer, "_", sample_type, "_immuneSig_merged_selected.pdf"),
            width = 5,
            height = 5, 
            onefile = TRUE)
for(x in seq(length(immuneSig_selected_boxplot))){
        print(immuneSig_selected_boxplot[[x]])
}
dev.off()

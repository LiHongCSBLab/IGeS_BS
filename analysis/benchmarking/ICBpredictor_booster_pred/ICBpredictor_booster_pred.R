# rm(list=ls())
# Tumor genes - infiltrated immune cell fraction
workpath <- "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/"
setwd(workpath)

options(stringsAsFactors = F)

# Checking and loading R packages -----------------------------------------
 
if (!require(readxl)) {
  install.packages('readxl')
}
if (!require(data.table)) {
  install.packages('data.table')
}
if (!require(dplyr)) {
  install.packages('dplyr')
}

# if (!require(TCGAbiolinks)) {
#   BiocManager::install("TCGAbiolinks")
# }
if (!require(SummarizedExperiment)) {
  BiocManager::install("SummarizedExperiment")
}
if (!require(GSVA)) {
  BiocManager::install("GSVA")
}
if (!require(limma)) {
  BiocManager::install("limma")
}
if (!require(ConsensusClusterPlus)) {
  BiocManager::install("ConsensusClusterPlus")
}

if (!require(immunedeconv)) {
  devtools::install_github("grst/immunedeconv")
}

library(ROCR)
library(reshape2)

# Loading relied functions -----------------------------------------------------
# # TIL estimation, TIDE, TIP, IPS and Immune Resistance Pragram
# source("02_tumorGene_immuneSig_correlation/scripts/functions/02_immuneSig_generation_function.R")
# source("02_tumorGene_immuneSig_correlation/scripts/functions/TILestimator.R")
# source("02_tumorGene_immuneSig_correlation/scripts/functions/predictor_IPS.R")
# source("02_tumorGene_immuneSig_correlation/scripts/functions/predictor_icb_ICI_resistance_program.R")
# source("02_tumorGene_immuneSig_correlation/scripts/functions/immuneSig_sig160.R")
# source("02_tumorGene_immuneSig_correlation/scripts/functions/immuneSig_TCGA_caf_immuneCLSssgsea.R")
# # TIGS, AUC and Wilcox test
# source("05_immuneSig_ICBresponse/scripts/functions/05_immunotherapy_immuneSig_analysis_function.R")
# source("05_immuneSig_ICBresponse/scripts/functions/05_immunotherapy_immuneSig_analysis_function.R")
# source("05_immuneSig_ICBresponse/scripts/functions/predictor_ITS.R")

# source("06_results_summaryandplotting/scripts/immunesig_selection_function.R")


source("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/07_plotting_v2/scripts/revision/ICBpredictor_booster_pred/predictor_TIDEandTIP.R")
source("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/07_plotting_v2/scripts/revision/ICBpredictor_booster_pred/predictor_IPS.R")
source("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/07_plotting_v2/scripts/revision/ICBpredictor_booster_pred/predictor_icb_ICI_resistance_program.R")

# set result saving path -------------------------------------------------------
savepath="07_plotting_v2/ICBpredictor_booster_pred"


# generating immune signatures -------------------------------------------------


# dataset = '70138'
# cancer = 'LIHC'
# sample_type = "Primary"
# dir.create(paste0(savepath, '/', dataset,'/'))
# dir.create(paste0(savepath, '/', dataset,'/', cancer, '/'))

# processedDataPath = "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/03_drug_immuneSig_enrichment/data/drug_treated_expression_matrix/"
# drug_id = list.files(paste0(processedDataPath, '/', dataset,'/', cancer))

# for(drug in drug_id){
#     # drug='drug82'
#     dir.create(paste0(savepath, '/', dataset, '/', cancer, '/'))
#     dir.create(paste0(savepath, '/', dataset, '/', cancer, '/', drug, "/"))
#     drug_files = list.files(paste0(processedDataPath, '/', dataset, '/', cancer, '/', drug))

#     for(filename in drug_files){
#         # filename=drug_files[1]       
#         dir.create(paste0(savepath, '/', dataset, '/', cancer, '/', drug, "/", gsub('.txt', '', filename)))
#         result_dir = paste0(savepath, '/', dataset, '/', cancer, '/', drug, "/", gsub('.txt', '', filename), '/')
#         exp_mat = read.table(paste0(processedDataPath, '/', dataset, '/', cancer, '/', drug, "/", filename), sep='\t', row.names=1, header=T)
#         exp_mat = exp_mat[-1]
#         labelPath = "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/03_drug_immuneSig_enrichment/data/drug_treated_design/for_GSVA"
#         sample_label = read.table(paste0(labelPath, '/', dataset, '/', cancer, '/', drug, "/", filename), sep='\t', row.names=1, header=T)
        
#         # ICB Predictor ----------------------------------------------------------------
#         model_icb_tide_TIP(mat = exp_mat,
#                         label=sample_label,
#                         cancername = cancer,
#                         sample_type = sample_type,
#                         sig_dir = "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r2_drug_sensitize_icb_prediction/",
#                         save_path = result_dir)

#         model_icb_IPS(gene_expression = exp_mat,
#                       label=sample_label,
#                       cancer_type = cancer,
#                       sample_type = sample_type,
#                       sig_dir = "02_tumorGene_immuneSig_correlation/scripts/softwares/Immunophenogram-master/",
#                       save_path = result_dir)

#         model_icb_ICI_resistance_program(gene_expression = exp_mat,
#                                          label=sample_label,
#                                         cancer_type = cancer,
#                                         sample_type = sample_type,
#                                         sig_dir="02_tumorGene_immuneSig_correlation/scripts/softwares/",
#                                         save_path = result_dir)
#     }
#     print(paste0("Finished for", drug))
# }


dataset = '92742'
cancer = 'LIHC'
sample_type = "Primary"
dir.create(paste0(savepath, '/', dataset,'/'))
dir.create(paste0(savepath, '/', dataset,'/', cancer, '/'))

processedDataPath = "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/03_drug_immuneSig_enrichment/data/drug_treated_expression_matrix/"
drug_id = list.files(paste0(processedDataPath, '/', dataset,'/', cancer))

for(drug in drug_id){
    # drug='drug82'
    dir.create(paste0(savepath, '/', dataset, '/', cancer, '/'))
    dir.create(paste0(savepath, '/', dataset, '/', cancer, '/', drug, "/"))

    drug_files = list.files(paste0(processedDataPath, '/', dataset, '/', cancer, '/', drug))

    for(filename in drug_files){
        # filename=drug_files[1]
        
        dir.create(paste0(savepath, '/', dataset, '/', cancer, '/', drug, "/", gsub('.txt', '', filename)))
        result_dir = paste0(savepath, '/', dataset, '/', cancer, '/', drug, "/", gsub('.txt', '', filename), '/')

        processedDataPath = "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/03_drug_immuneSig_enrichment/data/drug_treated_expression_matrix/"
        exp_mat = read.table(paste0(processedDataPath, '/', dataset, '/', cancer, '/', drug, "/", filename), sep='\t', row.names=1, header=T)
        exp_mat = exp_mat[-1]

        labelPath = "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/03_drug_immuneSig_enrichment/data/drug_treated_design/for_GSVA"
        sample_label = read.table(paste0(labelPath, '/', dataset, '/', cancer, '/', drug, "/", filename), sep='\t', row.names=1, header=T)
        
        # ICB Predictor ----------------------------------------------------------------

        model_icb_tide_TIP(mat = exp_mat,
                        label=sample_label,
                        cancername = cancer,
                        sample_type = sample_type,
                        sig_dir = "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r2_drug_sensitize_icb_prediction/",
                        save_path = result_dir)


        model_icb_IPS(gene_expression = exp_mat,
                      label=sample_label,
                      cancer_type = cancer,
                      sample_type = sample_type,
                      sig_dir = "02_tumorGene_immuneSig_correlation/scripts/softwares/Immunophenogram-master/",
                      save_path = result_dir)

        # model_icb_ICI_resistance_program(gene_expression = exp_mat,
        #                                  label=sample_label,
        #                                 cancer_type = cancer,
        #                                 sample_type = sample_type,
        #                                 sig_dir="02_tumorGene_immuneSig_correlation/scripts/softwares/",
        #                                 save_path = result_dir)

    }
}


# rm(list=ls())
work_path <- "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/"

setwd(work_path)
options(stringsAsFactors = F)
set.seed(1234)

dir.create('07_plotting_v2/05_drugITSprofiler_CombatRef_metaPlusFilter/')
dir.create('07_plotting_v2/05_drugITSprofiler_CombatRef_metaPlusFilter/drugrank_diff_lincsPortal2_metaPlusFilter/')
dir.create('07_plotting_v2/05_drugITSprofiler_CombatRef_metaPlusFilter/drugrank_diff_lincsPortal2_metaPlusFilter/70138/')
dir.create('07_plotting_v2/05_drugITSprofiler_CombatRef_metaPlusFilter/drugrank_diff_lincsPortal2_metaPlusFilter/92742/')

drug_proof <- read.table("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/data/FDAapproved_p3Trial_IOtherapy/drug_list_proof_cancertype.txt", header = T, sep = '\t')
druglist <- unique(drug_proof[drug_proof$if_with_ICB == 'yes', ])
druglist <- druglist[!is.na(druglist$drug_index), ]
druglist <- druglist[druglist$proof != 'animal',]


# load( "07_plotting_v2/data/drugTarget_merged3DB.Rdata")
# drug_FDAcancer <- unique(drugTarget_FDA[which(drugTarget_FDA$ifCancer_drugbank == 'yes'),])
# # drug_FDAcancer <- data.frame(drugname = unique(DTFDA_drugbank$name), drugname_lower = tolower(unique(DTFDA_drugbank$name)))
# # drug_FDAcancer <- data.frame(drugname = unique(DTFDAcancer_drugbank$name), drugname_lower = tolower(unique(DTFDAcancer_drugbank$name)))



# drugMoA <- read.csv('/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/03_drug_immuneSig_enrichment/data/drug_moa_merged_cleaned.txt', sep = '\t', header = T)
# drugMoA <- read.csv('/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/s3_lincs_drug_MoA/lincs_portal2.0/drugMoA_lincsPortal2.csv',header = F)
# drugMoA <- do.call(rbind, 
#                   lapply(as.list(drugMoA$V2), function(x){
#                   tmp = drugMoA[drugMoA$V2 == x,]
#                   df = data.frame(drug_index = x, MoA = unique(unlist(tmp[-c(1,2)])))
#                   return(df) }))
# drugMoA <- unique(drugMoA[drugMoA$MoA != "", ])
# write.csv(drugMoA, '/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/s3_lincs_drug_MoA/lincs_portal2.0/drugMoA_lincsPortal2_2.csv', quote = F, row.names = F)
drugMoA <- read.csv('/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/s3_lincs_drug_MoA/lincs_portal2.0/drugMoA_lincsPortal2_2.csv',header = T)
drugTarget <- read.csv('/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/s3_lincs_drug_MoA/lincs_portal2.0/drugtarget_lincsPortal2.csv', header = T)
drug_maxfda <- read.csv('/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/s3_lincs_drug_MoA/lincs_portal2.0/drug_max_fda_phase_lincsPortal2.csv', header = T)
drug_maxfda34 <- unique(drug_maxfda[is.element(drug_maxfda$max_fda_phase,c(3,4)), ])

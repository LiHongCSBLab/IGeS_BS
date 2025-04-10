# This function is for TCGA patient selection
# The basic rules are:
# 1) For all TCGA cancer types,
#    patient should seperate based on the original tissue status,
#    more specifically, primary or metastasis.
#    Note that "normal" tissue should also be excluded.
# 2) For BRCA, patients will be further grouped accouding to tumor subtypes
#    (Basel, Her2, Lum A, Lum B and TN).
# 3) For COAD/READ, Microsatellite Stability influenced both clinical outcome
#    and immunotherapy response (Reference), therefore patients will be assigned
#    to MSS, MSI-L and MSI-H subtypes.
# 4)(Optional) Some researches bring out the effect of viral-induced immune
#    reaction which influence immune microenvironment (immune) in tumor tissue,
#    therefore, cancer types like LIHC, HNSC and CESC can be further grouped
#    based on virus infection status.
# Note that these criteria should be consistent in later analysis.

# work_path <- "/picb/bigdata/project/FengFYM/immune_targeted_combo/immunecell_DEG_cor/"
workpath <- "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/"

setwd(work_path)
library(data.table)
options(stringsAsFactors = F)


# 1)
exprpath = "/picb/bigdata/project/FengFYM/UCEC_Xena/EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena"
tcgaRNAmatrix = as.data.frame(data.table::fread(exprpath))

# Sample Type Codes
# 01    Primary solid tumor
# 02    Recurrent Solid Tumor
# 03    Primary Blood Derived Cancer - Peripheral Blood
# 05    Additional - New Primary
# 06    Metastatic
# 07    Additional Metastatic
# 11    Solid Tissue Normal
patient_ID <- data.frame(SampleID = substr(colnames(tcgaRNAmatrix)[-1], 1, 12),
                         TGCA_barcode = colnames(tcgaRNAmatrix)[-1],
                         tumor_type_index = substring(colnames(tcgaRNAmatrix)[-1],14,15),
                         tumor_type = "Primary_solid_tumor")
patient_ID[patient_ID$tumor_type_index %in% "02",]$tumor_type <- "Recurrent_Solid_Tumor"
patient_ID[patient_ID$tumor_type_index %in% "03",]$tumor_type <- "Primary_Blood_Derived_Cancer"
patient_ID[patient_ID$tumor_type_index %in% "05",]$tumor_type <- "Additional_New_Primary"
patient_ID[patient_ID$tumor_type_index %in% "06",]$tumor_type <- "Metastatic"
patient_ID[patient_ID$tumor_type_index %in% "07",]$tumor_type <- "Additional_Metastaticy"
patient_ID[patient_ID$tumor_type_index %in% "11",]$tumor_type <- "Solid_Tissue_Normal"

patient_ID$tumor_type_merged <- "Others"
patient_ID[patient_ID$tumor_type_index %in% c("01", "03"),]$tumor_type_merged <- "Primary"
patient_ID[patient_ID$tumor_type_index %in% "06",]$tumor_type_merged <- "Metastatic"
patient_ID[patient_ID$tumor_type_index %in% "11",]$tumor_type_merged <- "Normal"

# 2) BRCA
sample_types <- read.table( paste0(work_path, "data/TCGA_cancertypes_subtypes.txt"), header=T, sep="\t")
names(sample_types) <- c("SampleID", "TCGA_cancer_type", "Immune.Subtype", "TCGA.Subtype")
patient_ID <- merge(patient_ID, sample_types, by = "SampleID", all.x = T)

# 3) COAD/READ

COAD_MSI <- read.table("/picb/bigdata/project/FengFYM/s6_TCGAbiolinks_download_normalization/clinical_data/TCGA_COAD_clinical.msi.tsv", sep = "\t", header = T)[c(1, 3)]
READ_MSI <- read.table("/picb/bigdata/project/FengFYM/s6_TCGAbiolinks_download_normalization/clinical_data/TCGA_READ_clinical.msi.tsv", sep = "\t", header = T)[c(1, 3)]
names(COAD_MSI) = c("SampleID", "MSI_status")
names(READ_MSI) = c("SampleID", "MSI_status")

patient_ID = merge(patient_ID, rbind(COAD_MSI, READ_MSI), by = "SampleID", all.x = T)
write.csv(patient_ID, "data/TCGA_patient_ID_tumor_type_all.csv", quote = F, row.names = F)


patient_ID_SKCM = unique(patient_ID[patient_ID$TCGA_cancer_type == "SKCM" & patient_ID$tumor_type_merged %in% c("Primary", "Metastatic"), ])
patient_ID_others = unique(patient_ID[patient_ID$TCGA_cancer_type != "SKCM" & patient_ID$tumor_type_merged %in% "Primary", ])
patient_ID_SKCM = patient_ID_SKCM[!is.na(patient_ID_SKCM$SampleID),]
patient_ID_others = patient_ID_others[!is.na(patient_ID_others$SampleID),]

patient_ID_selected = rbind(patient_ID_SKCM, patient_ID_others)
patient_ID_selected = patient_ID_selected[order(patient_ID_selected$TCGA_cancer_type, decreasing = F),]
write.csv(patient_ID_selected, "data/TCGA_patient_ID_tumor_type_selected.csv", quote = F, row.names = F)


cancerlist = unique(patient_ID$TCGA_cancer_type)
unique(patient_ID$TCGA.Subtype)
table(patient_ID$TCGA.Subtype)
for(i in 1:length(cancerlist)){
    print(cancerlist[i])
    print(table(patient_ID[patient_ID$TCGA_cancer_type == cancerlist[i],]$tumor_type_merged))
}

table(sample_types$TCGA.Study)
table(patient_ID$TCGA_cancer_type)

table(patient_ID[c("TCGA_cancer_type", "tumor_type_merged")])
samples = intersect(colnames(tcgaRNAmatrix),immune_sig$SampleID)

sample_types2 = sample_types[is.element(samples, sample_types$TCGA.Participant.Barcode), ]
patient_ID2 = patient_ID[is.element(samples, patient_ID$SampleID), ]
table(sample_types2$TCGA.Study)
table(patient_ID2$TCGA_cancer_type)


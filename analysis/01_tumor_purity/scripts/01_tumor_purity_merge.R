setwd("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/")
options(stringsAsFactors = F)
library(dplyr)

filePath = "01_tumor_purity/data/rawData/"
savePath = "01_tumor_purity/data/processedData/"

# 1) IHC-derived tumor purity  from  Aran D, Sirota M, Butte AJ. Systematic 
# pan-cancer analysis of tumour purity. 
# Nat Commun. 2015;6:8971.

# 2) CPE - consensus tumor purity derived from IHC-based, ESTIMATE, ABSOLUTE and 
# LUMP methods from Aran D, Sirota M, Butte AJ. Systematic pan-cancer analysis 
# of tumour purity. Nat Commun. 2015;6:8971.
purity_IHCbased = read.csv(paste0(filePath, "IHC_estimated_purity.csv"))
purity_IHC_CPE = purity_IHCbased[c("Sample.ID", "Cancer.type","IHC","CPE")]

# 3) TUMERIC - consensus tumor purity derived from AbsCNseq, PurBayes, Ascat and 
# ESTIMATE. from Pan-cancer analysis of ligand-receptor crosstalk in the tumor 
# microenvironment
purity_TUMERIC = read.csv(paste0(filePath, "purity_separate.csv"))
purity_TUMERIC = purity_TUMERIC[c(1,2)]

# 4) tumor purity downloaded from PanCanAtlas Publications 
purity_PanCanAtlas = read.csv(paste0(filePath, 
                       "TCGA_mastercalls.abs_tables_JSedit.fixed_ABSOLUTE.txt"),
                       sep = '\t')
purity_PanCanAtlas = purity_PanCanAtlas[c(1,2,4)]

# column rename
names(purity_IHC_CPE) = c("Sample_ID2", "Cancer_type", "IHC", "CPE")
purity_IHC_CPE$Sample_ID = substr(purity_IHC_CPE$Sample_ID2,1,15)
names(purity_TUMERIC) = c("Sample_ID", "TUMERIC")
names(purity_PanCanAtlas) = c("Sample_ID", "Sample", "ABSOLUTE")

purity_merged = full_join(purity_IHC_CPE, purity_TUMERIC, 
                        by = "Sample_ID")
purity_merged = full_join(purity_PanCanAtlas, 
                          purity_merged,
                          by = "Sample_ID")
purity_merged = purity_merged[c("Sample",
                                "Sample_ID",
                                "Sample_ID2",
                                "Cancer_type",
                                "IHC",
                                "CPE",
                                "TUMERIC",
                                "ABSOLUTE")]

save(purity_merged, file = paste0(savePath, "TCGA_tumor_purity_merged.Rdata"))

cancers = unique(purity_merged$Cancer_type)
for(cancer in cancers[!is.na(cancers)]){
    tmp = purity_merged[is.element(purity_merged$Cancer_type,cancer),]
    print(nrow(tmp))
    print(summary(tmp[c("IHC",
                        "CPE",
                        "TUMERIC",
                        "ABSOLUTE")]))
}

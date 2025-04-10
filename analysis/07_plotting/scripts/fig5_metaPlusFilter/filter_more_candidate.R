res <- read.csv(paste0("07_plotting_v2/05_drugITSprofiler_CombatRef_metaPlusFilter/drugITSres_metaPlusFilter/meta_en_ACATP0.05_",f), sep = '\t', header = T)
res$drugname_lower <- tolower(res$pert_iname)
res <- merge(res, drug_maxfda34[c('drug_index', 'max_fda_phase')], by='drug_index', all.x=T)
res$cancer_dataset = paste0(res$cancer, "_", res$dataset)


res_annotated <- merge(res, drugAnnot_proof, by = 'pert_iname', all.x = T)
write.xlsx(res_annotated,
           '07_plotting_v2/05_drugITSprofiler_CombatRef_metaPlusFilter/res_annotatedFDA.xlsx')

dAproof70138
write.xlsx(dAproof70138,
           '07_plotting_v2/05_drugITSprofiler_CombatRef_metaPlusFilter/sharedDrug_metaPlusFilter_alldrugFDA_all/dAproof70138.xlsx')
write.xlsx(res70138_shared,
           '07_plotting_v2/05_drugITSprofiler_CombatRef_metaPlusFilter/sharedDrug_metaPlusFilter_alldrugFDA_all/res70138_shared.xlsx')

write.xlsx(dAproof92742,
           '07_plotting_v2/05_drugITSprofiler_CombatRef_metaPlusFilter/sharedDrug_metaPlusFilter_alldrugFDA_all/dAproof92742.xlsx')

write.xlsx(res92742_shared,
           '07_plotting_v2/05_drugITSprofiler_CombatRef_metaPlusFilter/sharedDrug_metaPlusFilter_alldrugFDA_all/res92742_shared.xlsx')

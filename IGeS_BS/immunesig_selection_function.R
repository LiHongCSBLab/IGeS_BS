immunesig_name_unify <- function(immunesig){
  if(nrow(immunesig[grep("_oe", immunesig$immune_sig),]) > 0){
    immunesig[grep("_oe", immunesig$immune_sig),]$immune_sig = paste0("oe_", gsub("_oe","",immunesig[grep("_oe", immunesig$immune_sig),]$immune_sig))
  }
  immunesig$immune_sig = gsub("c7_", "",immunesig$immune_sig)
  immunesig$immune_sig = gsub("xcell.", "",immunesig$immune_sig)
  immunesig$immune_sig = gsub("mcp_counter.", "",immunesig$immune_sig)
  immunesig$immune_sig = gsub("epic.", "",immunesig$immune_sig)
  immunesig$immune_sig = gsub("quantiseq.", "",immunesig$immune_sig)
  immunesig$immune_sig = gsub("timer.", "",immunesig$immune_sig)
  immunesig$immune_sig = gsub("_immuncellai", "_immucellai",immunesig$immune_sig)
  immunesig$immune_sig = gsub("-",".", immunesig$immune_sig)
  immunesig$immune_sig = gsub(" ",".", immunesig$immune_sig)
  immunesig$immune_sig = gsub("\\+",".", immunesig$immune_sig)
  immunesig$immune_sig = gsub("\\(",".", immunesig$immune_sig)
  immunesig$immune_sig = gsub("\\)",".", immunesig$immune_sig) 
  immunesig$immune_sig = gsub("macrophages","macrophage", immunesig$immune_sig)
  immunesig$immune_sig = gsub("neutrophils","neutrophil", immunesig$immune_sig)
  immunesig$immune_sig = gsub("nk.cells","nk.cell", immunesig$immune_sig)
  immunesig$immune_sig = gsub("t.cells.cd4","t.cell.cd4", immunesig$immune_sig)
  immunesig$immune_sig = gsub("b.cells.memory","b.cell.memory", immunesig$immune_sig)
  immunesig$immune_sig = gsub("b.cells.naive","b.cell.naive", immunesig$immune_sig)
  immunesig$immune_sig = gsub("t.cell.cd4.memory.activated_cibersort","t.cell.cd4..memory.activated_cibersort", immunesig$immune_sig)
  immunesig$immune_sig = gsub("t.cell.cd4.memory.resting_cibersort","t.cell.cd4..memory.resting_cibersort", immunesig$immune_sig)
  immunesig$immune_sig = gsub("t.cells.cd8_cibersort","t.cell.cd8._cibersort", immunesig$immune_sig)
  immunesig$immune_sig = gsub("t.cells.follicular.helper_cibersort","t.cell.follicular.helper_cibersort", immunesig$immune_sig)
  immunesig$immune_sig = gsub("cibersortx","cibersort", immunesig$immune_sig)
  immunesig$immune_sig = gsub("tip_signature_tip","tip_signature", immunesig$immune_sig)
  immunesig$immune_sig = gsub("tip_signature","tip_signature_tip", immunesig$immune_sig)
  immunesig$immune_sig = gsub("macrophage.m1_xcell","macrophages.m1_xcell", immunesig$immune_sig)
  immunesig$immune_sig = gsub("macrophage.m2_xcell","macrophages.m2_xcell", immunesig$immune_sig)
  immunesig$immune_sig = gsub("macrophage_xcell","macrophages_xcell", immunesig$immune_sig)
  immunesig$immune_sig = gsub("xcell.", "", immunesig$immune_sig)
  immunesig$immune_sig = gsub("mcp_counter.","", immunesig$immune_sig)
  immunesig$immune_sig = gsub("quantiseq.","", immunesig$immune_sig)
  immunesig$immune_sig = gsub("nk.cell_tcga_caf_immunecls","nk_cells_tcga_caf_immunecls", immunesig$immune_sig)
  immunesig$immune_sig = gsub("macrophage_tcga_caf_immunecls", "macrophages_tcga_caf_immunecls", immunesig$immune_sig)
  immunesig$immune_sig = gsub("t.cell.cd4_consensustmegsva", "t_cells_cd4_consensustmegsva", immunesig$immune_sig) # new
  immunesig$immune_sig = gsub("bcell_immucellai", "b_cell_immucellai", immunesig$immune_sig)
  immunesig$immune_sig = gsub("nk.cell_mcp_counter", "nk.cells_mcp_counter", immunesig$immune_sig)
  immunesig$immune_sig = gsub("neutrophil_xcell","neutrophils_xcell", immunesig$immune_sig)
  immunesig$immune_sig = gsub("nk.cell_xcell","nk.cells_xcell", immunesig$immune_sig)
  immunesig$immune_sig = gsub("mast.cells.resting_cibersort","mast.cell.resting_cibersort", immunesig$immune_sig)
  immunesig$immune_sig = gsub("_cibersort.abs","_cibersortabs", immunesig$immune_sig)
  
  immunesig$immune_sig = gsub("_tideweb","_tide", immunesig$immune_sig)
  immunesig$immune_sig = gsub("tam.m2_tide","m2_tide", immunesig$immune_sig)
  immunesig$immune_sig = gsub("eosinophils_cibersort","eosinophil_cibersort", immunesig$immune_sig)
  immunesig$immune_sig = gsub("mast.cells.activated_cibersort","mast.cell.activated_cibersort", immunesig$immune_sig)
  immunesig$immune_sig = gsub("monocytes_cibersort","monocyte_cibersort", immunesig$immune_sig)
  immunesig$immune_sig = gsub("dendritic.cells.activated_cibersort","myeloid.dendritic.cell.activated_cibersort", immunesig$immune_sig)
  immunesig$immune_sig = gsub("dendritic.cells.resting_cibersort","myeloid.dendritic.cell.resting_cibersort", immunesig$immune_sig)
  immunesig$immune_sig = gsub("t.cell.cd4..naive_cibersort","t.cell.cd4.naive_cibersort", immunesig$immune_sig)
  immunesig$immune_sig = gsub("t.cell.gamma.delta_cibersort","t.cells.gamma.delta_cibersort", immunesig$immune_sig)
  immunesig$immune_sig = gsub("t.cell.regulatory..tregs._cibersort","t.cells.regulatory..tregs._cibersort", immunesig$immune_sig)
  immunesig$immune_sig = gsub("plasma.cells_cibersort","b.cell.plasma_cibersort", immunesig$immune_sig)
  
  immunesig$immune_sig = gsub("tc_immucellai","cytotoxic_immucellai", immunesig$immune_sig)
  immunesig$immune_sig = gsub("tcm_immucellai","central_memory_immucellai", immunesig$immune_sig)
  immunesig$immune_sig = gsub("tem_immucellai","effector_memory_immucellai", immunesig$immune_sig)
  immunesig$immune_sig = gsub("tex_immucellai","exhausted_immucellai", immunesig$immune_sig)
  immunesig$immune_sig = gsub("tgd_immucellai","gamma_delta_immucellai", immunesig$immune_sig)
  
  immunesig$immune_sig = gsub("oe_tme.stromal","oe_tme.stroma", immunesig$immune_sig)
  immunesig$immune_sig = gsub("stromalscore_xcell","stromascore_xcell", immunesig$immune_sig)
  return(immunesig)
}




xcell_sig_remover <- function(){
  # prefilter out xcell unrelated cell type
  xcell_cellfamily <- read.csv(paste0(workpath, "02_tumorGene_immuneSig_correlation/data/immune_signatures/originalSignaturesGenes/immune_sig_original_files/XCELL/cells_families.txt"), sep = '\t')
  xcell_type1 <- xcell_cellfamily[is.element(xcell_cellfamily$Group, 
                                             c("B-cells",
                                               "CD4+ memory T-cells",
                                               "CD4+ T-cells",
                                               "CD8+ T-cells",
                                               "DC",
                                               "Endothelial cells",
                                               "Epithelial cells",
                                               "Fibroblasts",
                                               "Macrophages")), ]
  
  xcell_type2 <- xcell_cellfamily[is.element(xcell_cellfamily$Group, 
                                             "Parent"), ]
  xcell_type3 <- xcell_type2[is.element(xcell_type2$Type, 
                                        c("Lymphoid", 
                                          "Myeloid")), ]
  xcell_type4 <- xcell_type2[is.element(xcell_type2$Cells, 
                                        c("Endothelial cells", 
                                          "Epithelial cells",
                                          "Fibroblasts")), ]
  xcell_select <- rbind(xcell_type1, xcell_type3, xcell_type4)
  xcell_remove <- xcell_cellfamily[!is.element(xcell_cellfamily$Cells, xcell_select$Cells),]
  xcell_remove <- paste0(xcell_remove$Cells, "_xcell")
  xcell_remove <- gsub(" ", ".", xcell_remove)

  return(xcell_remove)
}




immune_sig_filter_v2 <- function(auc_threshold = 0.6,
                                 immuSigproof = 2, # how many datasets support the immune sigture to have prediction power
                                 ITSproof = 2, # how many datasets support the immune sigture to have prediction power
                                 workpath = "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/",
                                 savepath  = "06_results_summaryandplotting/data/immune_sig_selected_new/"){


  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(readxl)
  # prefilter out xcell unrelated cell type
  xcell_cellfamily <- read.csv(paste0(workpath, "02_tumorGene_immuneSig_correlation/data/immune_signatures/originalSignaturesGenes/immune_sig_original_files/XCELL/cells_families.txt"), sep = '\t')
  xcell_type1 <- xcell_cellfamily[is.element(xcell_cellfamily$Group, 
                                             c("B-cells",
                                               "CD4+ memory T-cells",
                                               "CD4+ T-cells",
                                               "CD8+ T-cells",
                                               "DC",
                                               "Endothelial cells",
                                               "Epithelial cells",
                                               "Fibroblasts",
                                               "Macrophages")), ]
  
  xcell_type2 <- xcell_cellfamily[is.element(xcell_cellfamily$Group, 
                                             "Parent"), ]
  xcell_type3 <- xcell_type2[is.element(xcell_type2$Type, 
                                        c("Lymphoid", 
                                          "Myeloid")), ]
  xcell_type4 <- xcell_type2[is.element(xcell_type2$Cells, 
                                        c("Endothelial cells", 
                                          "Epithelial cells",
                                          "Fibroblasts")), ]
  xcell_select <- rbind(xcell_type1, xcell_type3, xcell_type4)
  xcell_remove <- xcell_cellfamily[!is.element(xcell_cellfamily$Cells, xcell_select$Cells),]
  xcell_remove <- paste0(xcell_remove$Cells, "_xcell")
  xcell_remove <- gsub(" ", ".", xcell_remove)
  
  # select signatures with known impact direction on immunotherapy
  # ---------------------------------------------------------------------
  immuneSigInfoPath = "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/06_results_summaryandplotting/data/immene_cell_hieracy/"
  
  immunesigInfo <- as.data.frame(read_excel(paste0(immuneSigInfoPath, "immune_sig_hieracy_new.xlsx"), sheet = 1))
  immunesigInfo <- immunesigInfo[c(1,4,5:9)]
  immunesigInfo$immune_sig = tolower(immunesigInfo$immune_sig)
  immunesigInfo = immunesig_name_unify(immunesigInfo)
  
  # ---------------------------------------------------------------------
  icb_res_path = "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/05_immuneSig_ICBresponse/results_ITS_r_0.4_ITSproof_2/"
  # icb_res_path = "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/05_immuneSig_ICBresponse/results/"
  files = list.files(paste0(icb_res_path, "immunesig_icbResponse_all_results/"))
  files = files[grep("gsva.csv", files)]
  files = files[-grep("Lauss", files)]

  resall = lapply(as.list(files), function(x){
    df = read.csv(paste0(icb_res_path, "immunesig_icbResponse_all_results/",x))
    df$dataset = gsub("_test_result.csv","",x)
    return(df)
  })
  
  resall = unique(do.call(rbind,resall))
  resall <- resall[-which(is.element(resall$immune_sig, tolower(xcell_remove))),]
  write.csv(resall, paste0(icb_res_path, "immunesig_icbResponse_all_results.csv"), quote = F)
  
  # auc_threshold = 0.6
  resall2 = resall[which(resall$immuneSig_auc >= auc_threshold), ]
  tmp = as.data.frame(table(resall2$immune_sig))
  tmp <- tmp[order(tmp$Freq, decreasing = T), ]
  head(tmp)

  # immuSigproof = 0
  immuneSig_s = unique(resall2[which(resall2$immuneSig_R_greater_NR_P_value < 0.05),])
  immuneSig_s_selected <- as.data.frame(table(immuneSig_s$immune_sig))
  names(immuneSig_s_selected) <- c("immune_sig", "Freq")
  immuneSig_s_selected <- immuneSig_s_selected[order(immuneSig_s_selected$Freq, decreasing = T), ]
  table(immuneSig_s_selected$Freq)
  # sort(immuneSig_s_selected[immuneSig_s_selected$Freq >= immuSigproof,]$immune_sig)


  immuneSig_r = unique(resall2[which(resall2$immuneSig_R_less_NR_P_value < 0.05),])
  immuneSig_r_selected <- as.data.frame(table(immuneSig_r$immune_sig))
  names(immuneSig_r_selected) <- c("immune_sig", "Freq")
  immuneSig_r_selected <- immuneSig_r_selected[order(immuneSig_r_selected$Freq, decreasing = T), ]
  # sort(immuneSig_r_selected[immuneSig_r_selected$Freq >= immuSigproof,]$immune_sig)
  table(immuneSig_r_selected$Freq)


  res_s <- unique(immuneSig_s_selected[immuneSig_s_selected$Freq >= immuSigproof,]$immune_sig)
  res_r <- unique(immuneSig_r_selected[immuneSig_r_selected$Freq >= immuSigproof,]$immune_sig)
  if(length(which(is.element(res_s, intersect(res_s,res_r))))){
      res_s <- res_s[-which(is.element(res_s, intersect(res_s, res_r)))]}
  if(length(which(is.element(res_r, intersect(res_s,res_r))))){
      res_r <- res_r[-which(is.element(res_r, intersect(res_s, res_r)))]}
  
  immuneSig_s2 <- immuneSig_s[is.element(immuneSig_s$immune_sig, res_s),]
   table(table(immuneSig_s2$immune_sig))
   length(unique(immuneSig_s2$immune_sig))

  immuneSig_r2 <- immuneSig_r[is.element(immuneSig_r$immune_sig, res_r),]
  table(table(immuneSig_r2$immune_sig))
  length(unique(immuneSig_r2$immune_sig))

  
  res_sITS <- immuneSig_s2[which(immuneSig_s2$ITSp_res_wilcox_RgreaterNR < 0.05),]
  length(unique(res_sITS$immune_sig))

  res_rITS <- immuneSig_r2[which(immuneSig_r2$ITSp_res_wilcox_RlessNR < 0.05),]
  length(unique(res_rITS$immune_sig))


  # ITSproof = 2
  res_sITS_selected <- as.data.frame(table(res_sITS$immune_sig))
  names(res_sITS_selected) <- c("immune_sig", "Freq")
  res_sITS_selected <- res_sITS_selected[which(res_sITS_selected$Freq >= ITSproof), ]
  res_sITS_selected <- res_sITS_selected[order(res_sITS_selected$Freq, decreasing = T), ]
  res_sITS_selected_class <- inner_join(res_sITS_selected, immunesigInfo)
  table(res_sITS_selected_class$Freq)
  length(unique(res_sITS_selected_class$immune_sig))
#   sort(table(res_sITS_selected_class[res_sITS_selected_class$Freq > ITSproof,]$Index))


  res_rITS_selected <- as.data.frame(table(res_rITS$immune_sig))
  names(res_rITS_selected) <- c("immune_sig", "Freq")
  res_rITS_selected <- res_rITS_selected[which(res_rITS_selected$Freq >= ITSproof), ]
  res_rITS_selected <- res_rITS_selected[order(res_rITS_selected$Freq, decreasing = T), ]
  res_rITS_selected_class <- inner_join(res_rITS_selected, immunesigInfo)
  table(res_rITS_selected_class$Freq)
  length(unique(res_rITS_selected_class$immune_sig))
#   sort(table(res_rITS_selected_class[res_rITS_selected_class$Freq > ITSproof,]$Index))
  
  res_sITS_selected_class$flag = "sensitive"
  res_rITS_selected_class$flag = "resistant"

  ITSselected_class = rbind(res_sITS_selected_class, res_rITS_selected_class)
  # need to be polished and merged
  ITSselected2 = ITSselected_class[-c(grep("mesa", (ITSselected_class$immune_sig)),
                                      grep("mela", (ITSselected_class$immune_sig)),
                                      grep("astro", (ITSselected_class$immune_sig)),
                                      grep("hsc", (ITSselected_class$immune_sig)),
                                      grep("mep", (ITSselected_class$immune_sig)),
                                      grep("mus", (ITSselected_class$immune_sig)), 
                                      # grep("up", (ITSselected_class$immune_sig)), 
                                      # grep("down", (ITSselected_class$immune_sig)), 
                                      grep("az_", (ITSselected_class$immune_sig)),
                                      grep("ips_ips", (ITSselected_class$immune_sig)),
                                      grep("tide_tide", (ITSselected_class$immune_sig))),]
  
  
  # sort(table(ITSselected2$immune_sig))
  
  
  
  dir.create(paste0(workpath, "/", savepath))
  write.table(ITSselected2, 
            file = paste0(workpath, "/", savepath, "/selected_ITS.txt"),
            sep = "\t", 
            row.names = F, quote = F)
  
  return(ITSselected2)
  
  #   immuneSig_hieracy <- read.csv("06_results_summaryandplotting/data/immene_cell_hieracy/immune_sig_rough_filter.csv", row.names=1)
  #   immuneSig_hieracy[which(immuneSig_hieracy==''),] <- NA
  #   immuneSig_hieracy$immune_sig = row.names(immuneSig_hieracy)
  #   immuneSig_hieracy_selected <- merge(immunesig, immuneSig_hieracy, all.x=T, by = 'immune_sig')
  #   # View(immuneSig_hieracy_selected)
  #   # table(immuneSig_hieracy_selected$confident)
  #   dir.create("06_results_summaryandplotting/data/immune_sig_hieracy_selected/")
  #   dir.create(paste0("06_results_summaryandplotting/data/immune_sig_hieracy_selected/", purity_method))
  #   write.csv(immuneSig_hieracy_selected, 
  #             file = paste0("06_results_summaryandplotting/data/immune_sig_hieracy_selected/", purity_method,"/",cancer,"_",tumor_subtype,".csv"),
  #             row.names = F, quote = F)
  #   return(immuneSig_hieracy_selected)
}




jaccard_index <- function(x, y){
  length(intersect(x,y))/length(union(x,y))
}


# immune_sig_finefilter <- function(immunesig,
#                                   keyword = NULL,
#                                   row_selected = NULL,
#                                   cancer = "SKCM",
#                                   tumor_subtype = "Metastatic",
#                                   purity_method = "TUMERIC",
#                                   datatype="allgenes",
#                                   num_gene = 200,
#                                   workpath = "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/", 
#                                   savepath = "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/"
# ){
  
#   # library(org.Hs.eg.db)
#   # library(topGO)
#   # library(clusterProfiler)
  
#   dir.create("06_results_summaryandplotting/data/immunesig_jaccard_index")
#   dir.create(paste0("06_results_summaryandplotting/data/immunesig_jaccard_index/", purity_method))
#   dir.create(paste0("06_results_summaryandplotting/data/immunesig_jaccard_index/", 
#                     purity_method,"/",cancer,"_", tumor_subtype, "/"))
#   dir.create(paste0("06_results_summaryandplotting/data/immune_sig_integration/"))
#   dir.create(paste0("06_results_summaryandplotting/data/immune_sig_integration/",purity_method,"/"))
  
#   immuneSig_hieracy <- read.csv("06_results_summaryandplotting/data/immene_cell_hieracy/immune_sig_rough_filter.csv", row.names=1)
#   immuneSig_hieracy[which(immuneSig_hieracy==''),] <- NA
#   immuneSig_hieracy$immune_sig = row.names(immuneSig_hieracy)
#   immuneSig_hieracy_selected <- merge(immunesig, immuneSig_hieracy, all.x=T, by = 'immune_sig')
  
#   if(!is.null(keyword)){
#     immunesig = immuneSig_hieracy_selected[grep(keyword, immuneSig_hieracy_selected[, row_selected]), ][c("setindex", "immune_sig" )]
#   }
  
  
#   if(nrow(immunesig) > 1){
    
    
#     load(paste0("03_drug_immuneSig_enrichment/data/gene_immune_sig_", purity_method, "/for_GSVA/pcor/pearson_",cancer,"_", tumor_subtype, "_positive_200.Rdata"))
#     genelist_p_TUMERIC = genelist_p
#     names(genelist_p_TUMERIC) = tolower(names(genelist_p_TUMERIC))
    
#     if(length(names(genelist_p_TUMERIC)[grep("_oe", names(genelist_p_TUMERIC))]) > 0){
#       names(genelist_p_TUMERIC)[grep("_oe", names(genelist_p_TUMERIC))] = 
#         paste0("oe_", gsub("_oe","",names(genelist_p_TUMERIC)[grep("_oe", names(genelist_p_TUMERIC))]))
#     }
    
#     names(genelist_p_TUMERIC) = gsub("c7_", "",names(genelist_p_TUMERIC))
#     names(genelist_p_TUMERIC) = gsub("xcell.", "",names(genelist_p_TUMERIC))
#     names(genelist_p_TUMERIC) = gsub("mcp_counter.", "",names(genelist_p_TUMERIC))
#     names(genelist_p_TUMERIC) = gsub("epic.", "",names(genelist_p_TUMERIC))
#     names(genelist_p_TUMERIC) = gsub("quantiseq.", "",names(genelist_p_TUMERIC))
#     names(genelist_p_TUMERIC) = gsub("timer.", "",names(genelist_p_TUMERIC))
#     names(genelist_p_TUMERIC) = gsub("_immuncellai", "_immucellai",names(genelist_p_TUMERIC))
#     names(genelist_p_TUMERIC) = gsub("-",".", names(genelist_p_TUMERIC))
#     names(genelist_p_TUMERIC) = gsub(" ",".", names(genelist_p_TUMERIC))
#     names(genelist_p_TUMERIC) = gsub("\\+",".", names(genelist_p_TUMERIC))
#     names(genelist_p_TUMERIC) = gsub("\\(",".", names(genelist_p_TUMERIC))
#     names(genelist_p_TUMERIC) = gsub("\\)",".", names(genelist_p_TUMERIC))
#     names(genelist_p_TUMERIC) = gsub("macrophages","macrophage", names(genelist_p_TUMERIC))
#     names(genelist_p_TUMERIC) = gsub("neutrophils","neutrophil", names(genelist_p_TUMERIC))
#     names(genelist_p_TUMERIC) = gsub("nk.cells","nk.cell", names(genelist_p_TUMERIC))
#     names(genelist_p_TUMERIC) = gsub("b.cells.memory","b.cell.memory", names(genelist_p_TUMERIC))
#     names(genelist_p_TUMERIC) = gsub("b.cells.naive","b.cell.naive", names(genelist_p_TUMERIC))
#     names(genelist_p_TUMERIC) = gsub("t.cell.cd4..memory.activated_cibersort","t.cell.cd4.memory.activated_cibersort", names(genelist_p_TUMERIC))
#     names(genelist_p_TUMERIC) = gsub("t.cell.cd4..memory.resting_cibersort","t.cell.cd4.memory.resting_cibersort", names(genelist_p_TUMERIC))
#     names(genelist_p_TUMERIC) = gsub("t.cells.cd8_cibersort","t.cell.cd8._cibersort", names(genelist_p_TUMERIC))
#     names(genelist_p_TUMERIC) = gsub("t.cells.follicular.helper_cibersort","t.cell.follicular.helper_cibersort", names(genelist_p_TUMERIC))
#     names(genelist_p_TUMERIC) = gsub("cibersortx","cibersort", names(genelist_p_TUMERIC))
#     # names(genelist_p_TUMERIC) = gsub("tip_signature_tip","tip_signature", names(genelist_p_TUMERIC))
#     names(genelist_p_TUMERIC) = gsub("macrophage.m1_xcell","macrophages.m1_xcell", names(genelist_p_TUMERIC))
#     names(genelist_p_TUMERIC) = gsub("macrophage.m2_xcell","macrophages.m2_xcell", names(genelist_p_TUMERIC))
#     names(genelist_p_TUMERIC) = gsub("macrophage_xcell","macrophages_xcell", names(genelist_p_TUMERIC))
#     names(genelist_p_TUMERIC) = gsub("xcell.", "", names(genelist_p_TUMERIC))
#     names(genelist_p_TUMERIC) = gsub("mcp_counter.","", names(genelist_p_TUMERIC))
#     names(genelist_p_TUMERIC) = gsub("quantiseq.","", names(genelist_p_TUMERIC))
#     names(genelist_p_TUMERIC) = gsub("b_cell_immucellai","b_cell_immucellai", names(genelist_p_TUMERIC))
#     names(genelist_p_TUMERIC) = gsub("mast.cell.resting_cibersort","mast.cells.resting_cibersort", names(genelist_p_TUMERIC))
    
#     names(genelist_p_TUMERIC) = gsub("macrophage_tcga_caf_immunecls", "macrophages_tcga_caf_immunecls",names(genelist_p_TUMERIC))
#     names(genelist_p_TUMERIC) = gsub("mast.cells.resting_cibersort", "mast.cell.resting_cibersort",names(genelist_p_TUMERIC))
#     names(genelist_p_TUMERIC) = gsub("neutrophil_xcell", "neutrophils_xcell",names(genelist_p_TUMERIC))
#     names(genelist_p_TUMERIC) = gsub("nk.cell_tcga_caf_immunecls", "nk_cells_tcga_caf_immunecls",names(genelist_p_TUMERIC))
#     names(genelist_p_TUMERIC) = gsub("nk.cell_mcp_counter", "nk.cells_mcp_counter",names(genelist_p_TUMERIC))
#     names(genelist_p_TUMERIC) = gsub("nk.cell_xcell", "nk.cells_xcell",names(genelist_p_TUMERIC))
#     names(genelist_p_TUMERIC) = gsub("t.cell.cd4.memory.activated_cibersort", "t.cell.cd4..memory.activated_cibersort",names(genelist_p_TUMERIC))
#     names(genelist_p_TUMERIC) = gsub("t.cell.cd4.memory.resting_cibersort", "t.cell.cd4..memory.resting_cibersort",names(genelist_p_TUMERIC))
    
#     genelist_p_TUMERIC = genelist_p_TUMERIC[intersect(immunesig$immune_sig, names(genelist_p_TUMERIC))]
    
    
#     load(paste0("03_drug_immuneSig_enrichment/data/gene_immune_sig_", purity_method, "/for_GSVA/pcor/pearson_",cancer,"_", tumor_subtype, "_negative_200.Rdata"))
#     genelist_n_TUMERIC = genelist_n
#     names(genelist_n_TUMERIC)=tolower(names(genelist_n_TUMERIC))
    
#     if(length(names(genelist_n_TUMERIC)[grep("_oe", names(genelist_n_TUMERIC))]) > 0){
#       names(genelist_n_TUMERIC)[grep("_oe", names(genelist_n_TUMERIC))] = 
#         paste0("oe_", gsub("_oe","",names(genelist_n_TUMERIC)[grep("_oe", names(genelist_n_TUMERIC))]))
#     }
    
#     names(genelist_n_TUMERIC) = gsub("c7_", "",names(genelist_n_TUMERIC))
#     names(genelist_n_TUMERIC) = gsub("xcell.", "",names(genelist_n_TUMERIC))
#     names(genelist_n_TUMERIC) = gsub("mcp_counter.", "",names(genelist_n_TUMERIC))
#     names(genelist_n_TUMERIC) = gsub("epic.", "",names(genelist_n_TUMERIC))
#     names(genelist_n_TUMERIC) = gsub("quantiseq.", "",names(genelist_n_TUMERIC))
#     names(genelist_n_TUMERIC) = gsub("timer.", "",names(genelist_n_TUMERIC))
#     names(genelist_n_TUMERIC) = gsub("_immuncellai", "_immucellai",names(genelist_n_TUMERIC))
#     names(genelist_n_TUMERIC) = gsub("-",".", names(genelist_n_TUMERIC))
#     names(genelist_n_TUMERIC) = gsub(" ",".", names(genelist_n_TUMERIC))
#     names(genelist_n_TUMERIC) = gsub("\\+",".", names(genelist_n_TUMERIC))
#     names(genelist_n_TUMERIC) = gsub("\\(",".", names(genelist_n_TUMERIC))
#     names(genelist_n_TUMERIC) = gsub("\\)",".", names(genelist_n_TUMERIC))
#     names(genelist_n_TUMERIC) = gsub("macrophages","macrophage", names(genelist_n_TUMERIC))
#     names(genelist_n_TUMERIC) = gsub("neutrophils","neutrophil", names(genelist_n_TUMERIC))
#     names(genelist_n_TUMERIC) = gsub("nk.cells","nk.cell", names(genelist_n_TUMERIC))
#     names(genelist_n_TUMERIC) = gsub("b.cells.memory","b.cell.memory", names(genelist_n_TUMERIC))
#     names(genelist_n_TUMERIC) = gsub("b.cells.naive","b.cell.naive", names(genelist_n_TUMERIC))
#     names(genelist_n_TUMERIC) = gsub("t.cell.cd4..memory.activated_cibersort","t.cell.cd4.memory.activated_cibersort", names(genelist_n_TUMERIC))
#     names(genelist_n_TUMERIC) = gsub("t.cell.cd4..memory.resting_cibersort","t.cell.cd4.memory.resting_cibersort", names(genelist_n_TUMERIC))
#     names(genelist_n_TUMERIC) = gsub("t.cells.cd8_cibersort","t.cell.cd8._cibersort", names(genelist_n_TUMERIC))
#     names(genelist_n_TUMERIC) = gsub("t.cells.follicular.helper_cibersort","t.cell.follicular.helper_cibersort", names(genelist_n_TUMERIC))
#     names(genelist_n_TUMERIC) = gsub("cibersortx","cibersort", names(genelist_n_TUMERIC))
#     # names(genelist_n_TUMERIC) = gsub("tip_signature_tip","tip_signature", names(genelist_n_TUMERIC))
#     names(genelist_n_TUMERIC) = gsub("macrophage.m1_xcell","macrophages.m1_xcell", names(genelist_n_TUMERIC))
#     names(genelist_n_TUMERIC) = gsub("macrophage.m2_xcell","macrophages.m2_xcell", names(genelist_n_TUMERIC))
#     names(genelist_n_TUMERIC) = gsub("macrophage_xcell","macrophages_xcell", names(genelist_n_TUMERIC))
#     names(genelist_n_TUMERIC) = gsub("xcell.", "", names(genelist_n_TUMERIC))
#     names(genelist_n_TUMERIC) = gsub("mcp_counter.","", names(genelist_n_TUMERIC))
#     names(genelist_n_TUMERIC) = gsub("quantiseq.","", names(genelist_n_TUMERIC))
#     names(genelist_n_TUMERIC) = gsub("b_cell_immucellai","b_cell_immucellai", names(genelist_n_TUMERIC))
#     names(genelist_n_TUMERIC) = gsub("mast.cell.resting_cibersort","mast.cells.resting_cibersort", names(genelist_n_TUMERIC))
#     names(genelist_n_TUMERIC) = gsub("macrophage_tcga_caf_immunecls", "macrophages_tcga_caf_immunecls",names(genelist_n_TUMERIC))
#     names(genelist_n_TUMERIC) = gsub("mast.cells.resting_cibersort", "mast.cell.resting_cibersort",names(genelist_n_TUMERIC))
#     names(genelist_n_TUMERIC) = gsub("neutrophil_xcell", "neutrophils_xcell",names(genelist_n_TUMERIC))
#     names(genelist_n_TUMERIC) = gsub("nk.cell_tcga_caf_immunecls", "nk_cells_tcga_caf_immunecls",names(genelist_n_TUMERIC))
#     names(genelist_n_TUMERIC) = gsub("nk.cell_mcp_counter", "nk.cells_mcp_counter",names(genelist_n_TUMERIC))
#     names(genelist_n_TUMERIC) = gsub("nk.cell_xcell", "nk.cells_xcell",names(genelist_n_TUMERIC))
#     names(genelist_n_TUMERIC) = gsub("t.cell.cd4.memory.activated_cibersort", "t.cell.cd4..memory.activated_cibersort",names(genelist_n_TUMERIC))
#     names(genelist_n_TUMERIC) = gsub("t.cell.cd4.memory.resting_cibersort", "t.cell.cd4..memory.resting_cibersort",names(genelist_n_TUMERIC))
    
#     genelist_n_TUMERIC = genelist_n_TUMERIC[intersect(immunesig$immune_sig, names(genelist_n_TUMERIC))]
    
    
#     rm(genelist_p)
#     rm(genelist_n)
    
    
#     # filter out ITSs with too less genes
#     lincs_gene = read.csv("/picb/bigdata/project/FengFYM/drug_MoA_reanalysis/LINCS_data/GSE70138/file/GSE92742_Broad_LINCS_gene_info.txt", sep = "\t")
    
#     res = data.frame(immunesig = names(genelist_p_TUMERIC),
#                      gene_in_sig = 0,
#                      gene_in_lincs = 0,
#                      gene_in_lincs_all = 0,
#                      sig_lincs_overlap = 0,
#                      sig_lincs_all_overlap = 0
#     )
    
#     for(i in seq(length(genelist_p_TUMERIC))){
      
#       res[res$immunesig == names(genelist_p_TUMERIC)[i],]$gene_in_sig = length(genelist_p_TUMERIC[[i]])
#       res[res$immunesig == names(genelist_p_TUMERIC)[i],]$gene_in_lincs = nrow(lincs_gene[lincs_gene$pr_is_lm ==1,])
#       res[res$immunesig == names(genelist_p_TUMERIC)[i],]$gene_in_lincs_all = nrow(lincs_gene)
#       res[res$immunesig == names(genelist_p_TUMERIC)[i],]$sig_lincs_overlap = length(intersect(genelist_p_TUMERIC[[i]],
#                                                                                                lincs_gene[lincs_gene$pr_is_lm ==1,]$pr_gene_symbol))
#       res[res$immunesig == names(genelist_p_TUMERIC)[i],]$sig_lincs_all_overlap = length(intersect(genelist_p_TUMERIC[[i]],
#                                                                                                    lincs_gene$pr_gene_symbol))
      
#     }
    
#     res = res[res$immunesig %in% intersect(immunesig$immune_sig, res$immunesig), ]
    
#     pdf(paste0(savepath, "/06_results_summaryandplotting/gene_num_stat_new/",purity_method, "/", cancer,"_",tumor_subtype,".pdf"),
#         width = 10, height = 6)
#     par(mfrow = c(1,2))
#     plot(density(res$sig_lincs_all_overlap))
#     # hist(res$sig_lincs_all_overlap,freq=T)
#     m<-seq(0, 200, by=10) #设置一个区间范围
#     a = table(cut(res$sig_lincs_all_overlap,m)) #计算各个区间的频数
#     barplot(a)
#     dev.off()
    
#     write.csv(res, 
#               paste0(savepath, "/06_results_summaryandplotting/gene_num_stat_new/",purity_method, "/", cancer,"_",tumor_subtype,".csv"),
#               row.names = F, quote = F)
    
    
#     a = data.frame(a)
#     its_remove = res[res$sig_lincs_all_overlap < 200*0.05,]$immunesig
#     immunesig2 = immunesig[!is.element(immunesig$immune_sig, its_remove), ]
    
#     genelist_p_TUMERIC = genelist_p_TUMERIC[intersect(immunesig2$immune_sig, names(genelist_p_TUMERIC))]
#     genelist_n_TUMERIC = genelist_n_TUMERIC[intersect(immunesig2$immune_sig, names(genelist_n_TUMERIC))]
    
    
#     # create file to integrate highly dependent signatures
#     if(length(genelist_p_TUMERIC) >1 & length(genelist_n_TUMERIC) >1){
      
#       # 1) 
#       p_jcindex <- list()
#       for(i in 1:length(genelist_p_TUMERIC)){
#         x = genelist_p_TUMERIC[[i]]
#         p_jcindex[[i]] = sapply(genelist_p_TUMERIC, function(y){
#           jaccard_index(x,y)
#         })
#       }
#       names(p_jcindex) = names(genelist_p_TUMERIC)
#       # p_jcindex2 = do.call(rbind, p_jcindex)
#       # p_jcindex2[lower.tri(p_jcindex2)] <- NA
#       # diag(p_jcindex2) <- NA
#       # p_jcindex3 = reshape2::melt(p_jcindex2)
#       # p_jcindex3 = p_jcindex3[!is.na(p_jcindex3$value),]
      
#       p_mat = do.call(cbind, p_jcindex)
#       p_mat = p_mat[sort(colnames(p_mat)), ]
#       p_mat = p_mat[, sort(colnames(p_mat))]
      
#       # draw heatmap to see the pattern
#       pdf(paste0("06_results_summaryandplotting/data/immunesig_jaccard_index/", 
#                  purity_method,"/",cancer,"_", tumor_subtype, "/jci_p.pdf"), 
#           width=20, height = 18)
      
#       pheatmap::pheatmap(p_mat,
#                          cluster_cols = T,
#                          cluster_rows = T,
#                          display_numbers = T,
#                          number_format = "%.2f",
#                          fontsize = 5, 
#                          annotation_row = immuneSig_hieracy[1:2],
#                          color = colorRampPalette(c("steelblue", "white", "firebrick3"))(20),
#                          border_color = 'grey50' )
#       dev.off()
      
#       # cut distance
      
#       d = dist(p_mat, method = 'euclidean')
#       tree = hclust(d, method = 'complete')
      
#       treedf = as.data.frame(as.matrix(do.call(cbind, tree[c(4,2:3)])))
#       treedf$height = as.numeric(treedf$height)
#       treedf$order = as.numeric(treedf$order)
#       names(treedf) = c("immune_sig", "height","order")
#       treedf$immune_sig = as.character(treedf$immune_sig)
#       v = cutree(tree, 10)[tree$order]
#       gaps = which((v[-1] - v[-length(v)]) != 0)
#       its_cls <- as.data.frame(as.matrix(v))
#       its_cls$immune_sig <- str_split(rownames(its_cls),'[.]',simplify = T)[,1]
#       its_cls = inner_join(treedf, its_cls)
      
#       immunesig3 = merge(immunesig2, unique(its_cls), all.x = T, by = "immune_sig")
#       # set threshold
#       p_mat_upper = p_mat
#       p_mat_upper[lower.tri(p_mat_upper, diag = T)]=-666
#       p_mat_JCI = melt(p_mat_upper)
#       p_mat_JCI = p_mat_JCI[p_mat_JCI$value != -666, ]
#       p_mat_JCI = p_mat_JCI[order(p_mat_JCI$value, decreasing = T),]
      
#       # plot jaccard index distribution
#       # par(mfrow = c(1,2))
#       # plot(density(p_mat_JCI$value))
#       # hist(p_mat_JCI$value,freq=T)
#       # dev.off()
      
      
#       tmp = as.data.frame(as.matrix(p_mat_JCI[p_mat_JCI$value > 0.7,]))
#       # tmp = p_mat_JCI[1:round(nrow(p_mat_JCI)*0.05),]
#       # min(p_mat_JCI[1:round(nrow(p_mat_JCI)*0.05),]$value)
#       names(its_cls) = c("Var1","Var1_h","Var1_o","Var1_v")
#       tmp2 = merge(tmp, its_cls, by = "Var1")
#       names(its_cls) = c("Var2","Var2_h","Var2_o","Var2_v")
#       tmp2 = merge(tmp2, its_cls, by = "Var2")
      
#       x = tmp2[1,]
#       if(x$Var1_v == x$Var2_v){
#         y = data.frame(immune_sig = c(x$Var1, x$Var2), label = x$Var1_v)
#       }
#       for(i in 2:nrow(tmp2)){      
#         x = tmp2[i,]
#         if(x$Var1_v == x$Var2_v){
#           y2 = data.frame(immune_sig = c(x$Var1, x$Var2), label = x$Var1_v)
#         }
#         y = rbind(y, y2)
#       }
#       y = unique(y)
#       immunesig3 = merge(immunesig3, y, by = "immune_sig", all.x = T)
      
#       write.csv(immunesig3, paste0("06_results_summaryandplotting/data/immune_sig_integration/",purity_method,"/" ,  cancer,"_", tumor_subtype,".csv"),
#                 row.names = F, quote = F)
#       return(immunesig3)
      
#       immuneSig_hieracy <- read.csv("06_results_summaryandplotting/data/immene_cell_hieracy/immune_sig_rough_filter.csv", row.names=1)
#       immuneSig_hieracy[which(immuneSig_hieracy==''),] <- NA
      
#       immuneSig_hieracy <- immuneSig_hieracy[is.element(row.names(immuneSig_hieracy),row.names(p_mat)),]
#       if(!is.null(keyword)){
        
#         pdf(paste0("06_results_summaryandplotting/data/immunesig_jaccard_index/", 
#                    purity_method,"/",cancer,"_", tumor_subtype, "/", keyword, "_jci_p.pdf"), 
#             width=15, height = 10)
#         pheatmap::pheatmap(p_mat,
#                            cluster_cols = T,
#                            cluster_rows = T,
#                            display_numbers = T,
#                            number_format = "%.2f",
#                            # fontsize = 5, 
#                            annotation_row = immuneSig_hieracy[1:2],
#                            color = colorRampPalette(c("steelblue", "white", "firebrick3"))(20),
#                            border_color = 'grey50' )
#         dev.off()
        
#       }
      
      
#       n_jcindex <- list()
#       for(i in 1:length(genelist_n_TUMERIC)){
#         x = genelist_n_TUMERIC[[i]]
#         n_jcindex[[i]] = sapply(genelist_n_TUMERIC,function(y){
#           jaccard_index(x,y)
#         })
#       }
#       names(n_jcindex) = names(genelist_n_TUMERIC)
#       n_jcindex2 = do.call(rbind, n_jcindex)
#       n_jcindex2[lower.tri(n_jcindex2)] <- NA
#       diag(n_jcindex2) <- NA
#       n_jcindex3 = reshape2::melt(n_jcindex2)
#       n_jcindex3 = n_jcindex3[!is.na(n_jcindex3$value),]
      
#       n_mat = do.call(cbind, n_jcindex)
#       n_mat = n_mat[sort(colnames(n_mat)), ]
#       n_mat = n_mat[, sort(colnames(n_mat))]
      
#       if(!is.null(keyword)){
#         pdf(paste0("06_results_summaryandplotting/data/immunesig_jaccard_index/", purity_method,"/",
#                    cancer,"_", tumor_subtype, "/", keyword, "_jci_n.pdf"), 
#             width=15, height = 10)
#         pheatmap::pheatmap(n_mat,
#                            cluster_cols = T,
#                            cluster_rows = T,
#                            display_numbers = T,
#                            # fontsize = 5,
#                            number_format = "%.2f",
#                            color = colorRampPalette(c("steelblue", "white", "firebrick3"))(20),
#                            border_color = 'grey50' )
#         dev.off()
        
#       }else{
#         pdf(paste0("06_results_summaryandplotting/data/immunesig_jaccard_index/", purity_method,"/",
#                    cancer,"_", tumor_subtype, "/jci_n.pdf"), 
#             width=20, height = 18)
#         pheatmap::pheatmap(n_mat,
#                            cluster_cols = T,
#                            cluster_rows = T,
#                            display_numbers = T,
#                            fontsize = 5,
#                            number_format = "%.2f",
#                            color = colorRampPalette(c("steelblue", "white", "firebrick3"))(20),
#                            border_color = 'grey50' )
#         dev.off()
        
        
#       }
      
#       names(p_jcindex3) = c("Var1", "Var2", "jaccard_p")
#       names(n_jcindex3) = c("Var1", "Var2", "jaccard_n")
#       res_jcindex = inner_join(p_jcindex3, n_jcindex3)
      
#       write.csv(res_jcindex, 
#                 paste0("06_results_summaryandplotting/data/immunesig_jaccard_index/", 
#                        purity_method,"/",cancer,"_", tumor_subtype,  "/", keyword, "_jci.csv"), quote=F, row.names=F)
#       # tmp <- res_jcindex[grep('cd8',res_jcindex$Var1),];tmp <- tmp[grep('cd8',tmp$Var2),];tmp[tmp$jaccard_p>0.5,]
#       # tmp <- res_jcindex[grep('cd4',res_jcindex$Var1),];tmp <- tmp[grep('cd4',tmp$Var2),];tmp[tmp$jaccard_p>0.5,]
#       # tmp <- res_jcindex[grep('macrophage',res_jcindex$Var1),];tmp <- tmp[grep('macrophage',tmp$Var2),];tmp[tmp$jaccard_p>0.5,]
#       # tmp <- res_jcindex[grep('b.cell',res_jcindex$Var1),];tmp <- tmp[grep('b.cell',tmp$Var2),];tmp[tmp$jaccard_p>0.5,]
      
#       # tmp = dcast(res_jcindex, Var1~Var2,value.var='jaccard_p')
#       # rownames(tmp) = tmp$Var1
#       # tmp = tmp[-1]
#       # tmp =  tmp[sort(colnames(tmp))]
#       # tmp =  tmp[sort(colnames(tmp)), ]
#     }
#   }
# }



# Aanlysis immune signature with metrics
# 1) correlation between the rest signatures
# immunesigname <- list()
# immuneSig_hieracy <- list()
# immuneSig_hieracy_file <- file("06_results_summaryandplotting/data/immene_cell_hieracy/immune_sig_rough_filter.csv",open="r")
# n=1
# while ( TRUE ) {
#   line = readLines(immuneSig_hieracy_file, n = 1)
#   if ( length(line) == 0 ) {
#     break
#   }
#   immuneSig_hieracy[[n]] <- unlist(strsplit(line, split = ",")) # cat(n, line,sep = ",")
#   immunesigname[[n]] <- unlist(strsplit(line, split = ","))[1]
#   n = n+1
# }

# names(immuneSig_hieracy) <- unlist(immunesigname)
# immuneSig_hieracy <- immuneSig_hieracy[-1]
# immuneSig_hieracy <- do.call(rbind, lapply(immuneSig_hieracy, data.frame))
# immuneSig_hieracy$immunesig <- gsub("\\.1","",row.names(immuneSig_hieracy))
# immuneSig_hieracy$immunesig <- gsub("\\.2","",immuneSig_hieracy$immunesig)
# immuneSig_hieracy$immunesig <- gsub("\\.3","",immuneSig_hieracy$immunesig)
# immuneSig_hieracy$immunesig <- gsub("\\.4","",immuneSig_hieracy$immunesig)
# immuneSig_hieracy$immunesig <- gsub("\\.5","",immuneSig_hieracy$immunesig)
# immuneSig_hieracy$immunesig <- gsub("\\.6","",immuneSig_hieracy$immunesig)
# immuneSig_hieracy$immunesig <- gsub("\\.7","",immuneSig_hieracy$immunesig)
# immuneSig_hieracy$immunesig <- gsub("\\.8","",immuneSig_hieracy$immunesig)
# immuneSig_hieracy$immunesig <- gsub("\\.9","",immuneSig_hieracy$immunesig)
# names(immuneSig_hieracy) <- c('Var1', 'Var2')

# if(length(which(immuneSig_hieracy$Var1 == ""))>0){
# immuneSig_hieracy = immuneSig_hieracy[-which(immuneSig_hieracy$Var1 == ""), ]}
# if(length(which(immuneSig_hieracy$Var2 == ""))>0){
# immuneSig_hieracy = immuneSig_hieracy[-which(immuneSig_hieracy$Var2 == ""), ]}
# immuneSig_hieracy = immuneSig_hieracy[-which(immuneSig_hieracy$Var2 == "NA1"), ]
# immuneSig_hieracy = immuneSig_hieracy[-which(immuneSig_hieracy$Var2 == "NA"), ]
# immuneSig_hieracy$jaccard_p = 1
# immuneSig_hieracy$jaccard_n = 1


# immuneSig_hieracy_selected1 <- unique(inner_join(immuneSig_hieracy, immunesig_jci['Var1']))
# immuneSig_hieracy_selected2 <-  unique(inner_join(immuneSig_hieracy, immunesig_jci['Var2']))

# immunesig_NET <- rbind(immuneSig_hieracy_selected1, immuneSig_hieracy_selected2, immunesig_jci)

# library(igraph)
# library(visNetwork)

# net = igraph::graph_from_data_frame(d=immunesig_NET, 
#                                     vertices = union(immunesig_NET$Var1, immunesig_NET$Var2), 
#                                     directed = F)


# edge.color <- colorRampPalette(c("#D6D6D6","#383838"), alpha=TRUE)
# igraph::E(net)$color <- edge.color(igraph::ecount(net))
# igraph::V(net)$color = "#E7F3FD"
# data <- toVisNetworkData(net)

# visNetwork(nodes = data$nodes, edges = data$edges, width = "90%", height = "95vh") %>%
#   visNodes(size = 40) %>%
#   visHierarchicalLayout() %>% #direction = "LR", levelSeparation = 500
#   visIgraphLayout(layout = "layout_on_sphere",type = "full",randomSeed = 123) %>%
#   visOptions(highlightNearest = list(enabled = TRUE,  hideColor = "lightgrey", hover = T),
#              nodesIdSelection =list(enabled = TRUE), selectedBy = "group") %>%
#   # visConfigure(enabled = TRUE) %>%
#   addFontAwesome() %>%
#   visGroups(groupname = "Targets", color = "#89D0F5")%>%
#   visGroups(groupname = "ImmuneSig", color = "#FF9900")%>%
#   visLegend() %>%
#   visInteraction(navigationButtons = TRUE) %>%
#   visOptions(manipulation = TRUE) %>% 
#   visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
#   visSave(file = paste0(resultpath_drug, "networkAnalysis.html"))


# 2) function detection

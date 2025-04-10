
work_path <- "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/"

setwd(work_path)
options(stringsAsFactors = F)
  
library(reshape2)
library(ggplot2)
library(ggplotify)
library(grid)
library(cowplot)
library(patchwork)
library(dplyr)
source("06_results_summaryandplotting/scripts/plot_for_drug_function.R")
source("06_results_summaryandplotting/scripts/immunesig_selection_function.R")
source("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/scripts/drugITGSscore_lincs.R")
source("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/scripts/drugITGSscore_GEOdataset.R")

xcell_cellfamily <- read.csv("02_tumorGene_immuneSig_correlation/data/immune_signatures/originalSignaturesGenes/immune_sig_original_files/XCELL/cells_families.txt")

  # select signatures with known impact direction on immunotherapy
  immune_cell_icbresponse_hugo = read.csv("05_immuneSig_ICBresponse/results/immunesig_icbResponse_test_results/01_Hugo_dataset_test_result.csv", row.names = 1)
  immune_cell_icbresponse_riaz = read.csv("05_immuneSig_ICBresponse/results/immunesig_icbResponse_test_results/02_Riaz_dataset_test_result.csv", row.names = 1)
  immune_cell_icbresponse_gide = read.csv("05_immuneSig_ICBresponse/results/immunesig_icbResponse_test_results/05_Gide_dataset_test_result.csv", row.names = 1)
  immune_cell_icbresponse_kim = read.csv("05_immuneSig_ICBresponse/results/immunesig_icbResponse_test_results/06_Kim_dataset_test_result.csv", row.names = 1)
  immune_cell_icbresponse_Lauss = read.csv("05_immuneSig_ICBresponse/results/immunesig_icbResponse_test_results/08_Lauss_dataset_test_result.csv", row.names = 1)
  immune_cell_icbresponse_IMvigor210 = read.csv("05_immuneSig_ICBresponse/results/immunesig_icbResponse_test_results/07_IMvigor210_dataset_test_result.csv", row.names = 1)
  
  
  immunesig_response_predictor <- list(immune_cell_icbresponse_hugo, 
                                       immune_cell_icbresponse_riaz,
                                       immune_cell_icbresponse_gide,
                                       immune_cell_icbresponse_kim,
                                       immune_cell_icbresponse_Lauss,
                                       immune_cell_icbresponse_IMvigor210)
   immunesig_response_predictor[[1]]$dataset = "hugo_dataset"                              
  immunesig_response_predictor[[2]]$dataset = "riaz_dataset"                              
  immunesig_response_predictor[[3]]$dataset = "gide_dataset"                              
  immunesig_response_predictor[[4]]$dataset = "kim_dataset"                              
  immunesig_response_predictor[[5]]$dataset = "lauss_dataset"                              
  immunesig_response_predictor[[6]]$dataset = "IMvigor210_dataset"                              

  immunesig_response_predictor <- do.call( rbind,immunesig_response_predictor)
  for(cl in xcell_remove){
    immunesig_response_predictor = immunesig_response_predictor[-grep(cl, immunesig_response_predictor$immunesig),]
  }
  immunesig_response_predictor$immunesig <- tolower(immunesig_response_predictor$immunesig)
  immunesig_response_predictor = immunesig_response_predictor[-
    c(grep("az_ips",immunesig_response_predictor$immunesig), # duplicated with IPS score
      # immune landscape or sig160
      grep("hematopoietic.stem.cell", immunesig_response_predictor$immunesig), # explanation ???      
      # oe
      grep("melanoma.cells..tirosh.",immunesig_response_predictor$immunesig),# ICB predictor|Melanoma cell state
      grep("mapki.resis", immunesig_response_predictor$immunesig),  # MAPKi resistant
      # cell types 
      grep("uncharacterized", immunesig_response_predictor$immunesig)
    ), ]
  
  names(immunesig_response_predictor) = c("wilcox_P_value","R_greater_NR_P_value",
                                        "R_less_NR_P_value","auc","immune_sig","dataset")
  immunesig_response_predictor <- immunesig_name_unify(immunesig_response_predictor)
  write.csv(immunesig_response_predictor, "05_immuneSig_ICBresponse/results/immunesig_response_predictor_summary.csv",
            row.names = F, quote = F)



  dfs = df[df$auc > auc_threshold & df$R_greater_NR_P_value < 0.05, ]
  dfsm = dfs[c("immune_sig", "dataset", "auc")]

  p1=ggplot(dfsm, aes(x=dataset,y=immune_sig))+ 
    geom_tile(aes(fill=auc), color = 'white', alpha=0.8)+ 
    scale_fill_gradient2(low = "#6A82FB",mid="white", high = "#FC5C7D")+
    theme_bw()+
    ggtitle("immunesig_response_predictor-sensitive")+
    theme(axis.text.x=element_text(angle=90, size=9, hjust = 0,vjust=0.05))
  
  dfr = df[df$auc > auc_threshold & df$R_less_NR_P_value < 0.05, ]
  dfrm = dfr[c("immune_sig", "dataset", "auc")]
  p2=ggplot(dfrm, aes(x=dataset,y=immune_sig))+ 
    geom_tile(aes(fill=auc), color = 'white', alpha=0.8)+ 
    scale_fill_gradient2(low = "#FC5C7D",mid="white", high = "#6A82FB")+
    theme_bw()+
    ggtitle("immunesig_response_predictor-resistant")+
    theme(axis.text.x=element_text(angle=90, size=9, hjust = 0,vjust=0.05))
  p =   p1 + p2 
  ggsave("05_immuneSig_ICBresponse/results/immunesig_response_predictor_summary.pdf", p, 
        width = 13, height = 16)


# SURVIVAL --------------------------------------------------------------------

cancerlist <- c( "UCEC","LIHC", "PRAD", "PAAD", "BRCA", "READ","LUAD","LUSC",
                  "OV",  "SKCM","COAD", "DLBC", 
                  "UCS", "HNSC", "BLCA","CESC","KIRC", "CHOL", "ESCA", "PCPG",
                  "ACC","LGG", "GBM", "KICH",
                  "SARC","STAD","TGCT","THYM", "MESO", "KIRP","THCA")
tumor_subtype = "Primary"

cancer = "SKCM"
tumor_subtype = "Metastatic"
for(cancer in cancerlist){
      
  if(cancer %in%  c("COAD","READ")){
    if(tumor_subtype == "Primary"){
      load(paste0("04_immuneSig_survival/results/04_TCGA_survival_result/coxresult_", cancer, ".Rdata"))
      
    }else{
      load(paste0("04_immuneSig_survival/results/04_TCGA_survival_result/coxresult_", cancer,"_", tumor_subtype, ".Rdata"))
      
    }
    
  }else if(cancer == "BRCA"){
    
    if(tumor_subtype == "Primary"){
      load(paste0("04_immuneSig_survival/results/04_TCGA_survival_result/coxresult_", cancer, ".Rdata"))
      
    }else{
      load(paste0("04_immuneSig_survival/results/04_TCGA_survival_result/coxresult_", cancer,".", tumor_subtype, ".Rdata"))
      
    }
    
  }else if(cancer == "SKCM"){
    
    load(paste0("04_immuneSig_survival/results/04_TCGA_survival_result/coxresult_", cancer,"_", tumor_subtype, ".Rdata"))
    
  }else{
    load(paste0("04_immuneSig_survival/results/04_TCGA_survival_result/coxresult_", cancer, ".Rdata"))
  }
  
  #survival_result = result
  survival_result$immune_sig = rownames(survival_result)

  for(cl in xcell_remove){
    survival_result = survival_result[-grep(cl, survival_result$immune_sig),]
  }
  survival_result$immune_sig <- tolower(survival_result$immune_sig)
  survival_result = survival_result[-
    c(grep("az_ips",survival_result$immune_sig), # duplicated with IPS score
      # immune landscape or sig160
      grep("hematopoietic.stem.cell", survival_result$immune_sig), # explanation ???      
      # oe
      grep("melanoma.cells..tirosh.",survival_result$immune_sig),# ICB predictor|Melanoma cell state
      grep("mapki.resis", survival_result$immune_sig),  # MAPKi resistant
      # cell types 
      grep("uncharacterized", survival_result$immune_sig)
    ), ]
  
  survival_result <- immunesig_name_unify(survival_result)

  df = survival_result[which(survival_result$p < 0.05),]
  write.csv(df, 
            paste0("04_immuneSig_survival/results/result_select/", cancer,"_", tumor_subtype, ".csv"),
            row.names = F, quote = F)


}


df=list()
cancer = "SKCM"
tumor_subtype = "Metastatic"

df[[1]] =  read.csv(paste0("04_immuneSig_survival/results/result_select/", cancer,"_", tumor_subtype, ".csv"))
df[[1]]$cancertype = paste0(cancer,"_",tumor_subtype)

tumor_subtype = "Primary"
# cancerlist <- c( "UCEC","LIHC", "PRAD", "PAAD", "BRCA", "READ","LUAD","LUSC",
#                   "OV",  "SKCM","COAD")
cancerlist <- c( "UCEC","LIHC", "PRAD", "PAAD", "BRCA", "READ","LUAD","LUSC",
                  "OV",  "SKCM","COAD", "DLBC", 
                  "UCS", "HNSC", "BLCA","CESC","KIRC", "CHOL", "ESCA", "PCPG",
                  "ACC","LGG", "GBM", "KICH",
                  "SARC","STAD","TGCT","THYM", "MESO", "KIRP","THCA")


for(i in 2: length(cancerlist)){
    cancer = cancerlist[i]
    df[[i]] =  read.csv(paste0("04_immuneSig_survival/results/result_select/", cancer,"_", tumor_subtype, ".csv"))
    df[[i]]$cancertype = paste0(cancer,"_",tumor_subtype)
}
df = do.call(rbind, df)
df$immune_sig = tolower(df$immune_sig)
df = immunesig_name_unify(df)



df_increaserisk = df[df$p < 0.05 & df$z > 0, ]
df_increaserisk = df_increaserisk[order(df_increaserisk$immune_sig), ]
tmp1 = as.data.frame(as.matrix(table(df_increaserisk$immune_sig)))
tmp1$immune_sig = row.names(tmp1)
df_increaserisk = inner_join(df_increaserisk, tmp1)
df_increaserisk = df_increaserisk[order(df_increaserisk$V1, decreasing = T), ]
tmp1 = tmp1[order(tmp1$V1, decreasing = T),]
tmp1$immune_sig = factor(tmp1$immune_sig, levels = tmp1$immune_sig)
df_increaserisk$immune_sig = factor(df_increaserisk$immune_sig, levels = (levels(tmp1$immune_sig)))

df_decreaserisk = df[df$p < 0.05 & df$z < 0, ]
df_decreaserisk = df_decreaserisk[order(df_decreaserisk$immune_sig), ]
tmp2 = as.data.frame(as.matrix(table(df_decreaserisk$immune_sig)))
tmp2$immune_sig = row.names(tmp2)
df_decreaserisk = inner_join(df_decreaserisk, tmp2)
df_decreaserisk = df_decreaserisk[order(df_decreaserisk$V1, decreasing = T), ]
tmp2 = tmp2[order(tmp2$V1, decreasing = T),]
tmp2$immune_sig = factor(tmp2$immune_sig, levels = tmp2$immune_sig)
df_decreaserisk$immune_sig = factor(df_decreaserisk$immune_sig, levels = (levels(tmp2$immune_sig)))

  p3=ggplot(df_decreaserisk, aes(x=cancertype,y=immune_sig))+ 
    geom_tile(aes(fill=z), color = 'white', alpha=0.8)+ 
    scale_fill_gradient2(low = "#FC5C7D",mid="white", high = "#6A82FB")+
    theme_bw()+
    ggtitle("survival analysis - decrease risk")+
    theme(axis.text.x=element_text(angle=90, size=9, hjust = 0,vjust=0.05))

  p4=ggplot(df_increaserisk, aes(x=cancertype,y=immune_sig))+ 
    geom_tile(aes(fill=z), color = 'white', alpha=0.8)+ 
    scale_fill_gradient2(low = "#FC5C7D",mid="white", high = "#6A82FB")+
    theme_bw()+
    ggtitle("survival analysis - increase risk")+
    theme(axis.text.x=element_text(angle=90, size=9, hjust = 0,vjust=0.05))
   p=   p3 + p4
  ggsave("04_immuneSig_survival/results/survival_summary.pdf", p, 
        width = 13, height = 20)


# ------------------------------------------------------------------------------


library(getopt)
options(stringsAsFactors = F)
options(quote = NULL)
options(stringsAsFactors=F)

command=matrix(c( 
  "cancer",        "c", 1, "character", "the path of gene and fold-change",
  "dataset",       "d", 1, "character", "dataset ID",
  "work_path",     "w", 1, "character", "the work path for reading and saving data",
  "result_path",   "r", 1, "character", "the path for saving result",
  "help",          "h", 0 , "logical",   "help file"),
  byrow=T,ncol=5)

args=getopt(command)

work_path = args$work_path 
result_path = args$result_path 
cancer = args$cancer
dataset = args$dataset
result_path=args$result_path

# cancer = "BRCA" 
# dataset = "70138"
# work_path <- "/picb/bigdata/project/FengFYM/drug_MoA_reanalysis/"

setwd(work_path)
library(ggplot2)
library(prada)
library(rhdf5)
library(limma)
library(dplyr)
library(fs)

# source("l1ktools-master/R/cmap/io.R")
func1=path("dependency/l1ktools-master/R/cmap/", "io.R")
source(func1)



print("start to prepare expression profile")
# result_path <- "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/03_drug_immuneSig_enrichment/data/"
result_path <- paste0(result_path, '/data/')
dir.create(paste0(result_path, "drug_treated_expression_matrix/"))
dir.create(paste0(result_path, "drug_treated_expression_matrix/", dataset))
dir.create(paste0(result_path, "drug_treated_expression_matrix/", dataset, "/", cancer))
# dir.create(paste0(result_path, "drug_treated_expression_978genes_matrix/"))
# dir.create(paste0(result_path, "drug_treated_expression_978genes_matrix/", dataset))
# dir.create(paste0(result_path, "drug_treated_expression_978genes_matrix/", dataset, "/", cancer))

dir.create(paste0(result_path, "drug_treated_design/"))
dir.create(paste0(result_path, "drug_treated_design/for_GSEA/"))
# dir.create(paste0(result_path, "drug_treated_design/for_GSVA/"))
dir.create(paste0(result_path, "drug_treated_design/for_GSEA/", dataset))
# dir.create(paste0(result_path, "drug_treated_design/for_GSVA/", dataset))
dir.create(paste0(result_path, "drug_treated_design/for_GSEA/", dataset, "/", cancer, "/"))
# dir.create(paste0(result_path, "drug_treated_design/for_GSVA/", dataset, "/", cancer, "/"))

# drug_diff_res <- result_path  # "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/s2_drug_inducted_transcriptome_change/"
# dir.create(paste0(drug_diff_res, "/allgenes/"))
# dir.create(paste0(drug_diff_res, "/allgenes/", dataset))
# dir.create(paste0(drug_diff_res, "/allgenes/", dataset, "/", cancer))
# dir.create(paste0(drug_diff_res, "/978genes/"))
# dir.create(paste0(drug_diff_res, "/978genes/", dataset))
# dir.create(paste0(drug_diff_res, "/978genes/", dataset, "/", cancer))


# inst_GSE70138 <- read.csv("LINCS_data/GSE70138/file/GSE70138_Broad_LINCS_inst_info.txt", header = T, sep = "\t")
# #   length(unique(inst_GSE70138$pert_iname))
# inst_GSE92742 <- read.csv("LINCS_data/GSE92742/file/GSE92742_Broad_LINCS_inst_info.txt", header=T, sep="\t")
# #   length(unique(inst_GSE92742$pert_iname))

# drug_index <- unique(rbind(
#     unique(inst_GSE70138["pert_iname"]),
#     unique(inst_GSE92742["pert_iname"])
# ))

# drug_index$drug_index <- paste0("drug", seq(nrow(drug_index)))
# write.table(drug_index,  paste0(result_path, 'drug_index_LINCS.txt'), sep = '\t', quote = F, row.names = F)
drug_index <- read.table(paste0(result_path, 'drug_index_LINCS.txt'), sep = '\t', header = T)

if(dataset == '70138'){
  
  # dataset <- "70138"
  inst_infopath=path("Data/LINCS_data/", "GSE70138_Broad_LINCS_inst_info.txt")
  inst_GSE70138 <- read.csv(inst_infopath, header = T, sep = "\t")
  # inst_GSE70138 <- read.csv("LINCS_data/GSE70138/file/GSE70138_Broad_LINCS_inst_info.txt", header = T, sep = "\t")
  inst_GSE70138 <- inst_GSE70138[which(is.element(inst_GSE70138$pert_type, c("ctl_vehicle", "trt_cp"))),]
  # assign index to drugs
  inst_GSE70138 <- inner_join(drug_index, inst_GSE70138)
  
  # update the cellline_cancertype table
  # cell_annot <-  read.csv("LINCS_data/GSE70138/file/GSE92742_Broad_LINCS_cell_info.txt", header = T, sep = "\t")
  # cellline_cancertype <- read.csv("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/s1_cellline_selection/TCGA_110CL/TCGA_CL_name.csv")
  # names(cellline_cancertype) <- c('cell_id','TCGA_cancer_type')
  # # filter out cell line that is not in cellline_cancertype
  # cell_annot_selected <- inner_join(cellline_cancertype, cell_annot)
  # cell_annot_selected <- cell_annot_selected[c('cell_id','TCGA_cancer_type')]
  # write.table(cell_annot_selected, paste0(result_path, 'cell_annot_selected_70138.txt'), sep = '\t', quote = F, row.names = F)
  cell_annot_selected_path=path("Data/LINCS_data/", "cell_annot_selected_70138.txt")
  cell_annot_selected <- read.table(cell_annot_selected_path, sep = '\t',header = T)
  # cell_annot_selected <- read.table(paste0(result_path, 'cell_annot_selected_70138.txt'), sep = '\t',header = T)

  inst_GSE70138 <- inner_join(cell_annot_selected, inst_GSE70138)
  inst_GSE70138_cancer <- inst_GSE70138[inst_GSE70138$TCGA_cancer_type == cancer, ]
  
  geneinfo_path=path("Data/LINCS_data/", "GSE92742_Broad_LINCS_gene_info.txt")
  geneinfo <- read.csv(geneinfo_path, sep = "\t", header = T)
  geneinfo_sig <- geneinfo[which(geneinfo$pr_is_lm == 1), ]
  
  for (drug in unique(inst_GSE70138_cancer[-which(inst_GSE70138_cancer$pert_iname=='DMSO'),]$drug_index)) {

    dir.create(paste0(result_path, "drug_treated_expression_matrix/", dataset, "/", cancer, "/", drug))
    # dir.create(paste0(result_path, "drug_treated_expression_978genes_matrix/", dataset, "/", cancer, "/", drug))
    dir.create(paste0(result_path, "drug_treated_design/for_GSEA/",   dataset, "/", cancer, "/", drug))
    # dir.create(paste0(result_path, "drug_treated_design/for_GSVA/",   dataset, "/", cancer, "/", drug))


    dir.create(paste0(drug_diff_res, "/allgenes/", dataset, "/", cancer, "/", drug))
    # dir.create(paste0(drug_diff_res, "/978genes/", dataset, "/", cancer, "/", drug))
    
    inst_GSE70138_selected <- inst_GSE70138_cancer[inst_GSE70138_cancer$drug_index == drug, ]
    
    for(cell in unique(inst_GSE70138_selected$cell_id)){
      for (treatment_time in unique(inst_GSE70138_selected$pert_time)) { 
        
        inst_GSE70138_selected2 <- inst_GSE70138_selected[inst_GSE70138_selected$pert_time == treatment_time & inst_GSE70138_selected$cell_id == cell, ]
        
        if (nrow(inst_GSE70138_selected2) != 0) {
          
          doses = unique(inst_GSE70138_selected2$pert_dose)
          # doses = doses[-grep(-666, doses)]
          for (treatment_dose in doses) { 
            
            inst_GSE70138_selected3 <- inst_GSE70138_selected2[inst_GSE70138_selected2$pert_dose == treatment_dose, ]
            
            if (nrow(inst_GSE70138_selected3) != 0) {
              
              DMSO <- inst_GSE70138_cancer[inst_GSE70138_cancer$pert_time == treatment_time &
                                             inst_GSE70138_cancer$cell_id == cell & 
                                             inst_GSE70138_cancer$pert_dose == -666, ]
              
              if (drug != "DMSO") {
                # drug <- unique(inst_GSE70138_selected2$pert_iname)[1]
                print(paste(cell, treatment_time, drug, treatment_dose))
                
                # d <- inst_GSE70138_selected2[inst_GSE70138_selected2$pert_iname == drug, ]
                d <- inst_GSE70138_selected3
                
                DMSO_match <- DMSO[DMSO$det_plate %in% d$det_plate, ]
                
                treatment_design <- data.frame(
                  treatment = c(
                    rep(1, length(d$inst_id)),
                    rep(0, length(DMSO_match$inst_id))
                  ),
                  treatment_ID = c(
                    d$inst_id,
                    DMSO_match$inst_id
                  ),
                  treatment_dose = c(
                    d$pert_dose,
                    DMSO_match$pert_dose
                  )
                )
                row.names(treatment_design) <- treatment_design$treatment_ID

                # write.table(treatment_design, paste0(result_path, "drug_treated_design/for_GSVA/",dataset, "/", cancer, "/", drug, "/", cell, "_", treatment_time,"_", treatment_dose, ".txt"),
                #             row.names = T, quote = F, sep = "\t")
                
                grp <- list(paste(nrow(treatment_design), 2, 1, collapse = " "),
                            paste("#", "treated", "DMSO", collapse = " "),
                            paste(c(rep("treated", nrow(treatment_design[treatment_design$treatment == 1,])),
                                    rep("DMSO", nrow(treatment_design[treatment_design$treatment == 0,]))),
                                  collapse = " "))
                grp <- do.call(rbind, grp)
                write.table(grp, paste0(result_path, "drug_treated_design/for_GSEA/", dataset, "/", cancer, "/", drug, "/", cell, "_", treatment_time,"_", treatment_dose, ".cls"),
                            col.names = F, row.names = F, quote = F, sep = "\t")

                gex_path=path("Data/LINCS_data/", "GSE70138_Broad_LINCS_Level3_INF_mlr12k_n345976x12328.gctx")
                # drug_70138_3h_L3 <- parse.gctx("LINCS_data/GSE70138/GSE70138_Broad_LINCS_Level3_INF_mlr12k_n345976x12328.gctx",
                #                                cid = treatment_design$treatment_ID)
                drug_70138_3h_L3 <- parse.gctx(gex_path,
                                               cid = treatment_design$treatment_ID)
                mat <- as.data.frame(drug_70138_3h_L3@mat)
                mat$pr_gene_id <- rownames(mat)

                mat_all <- merge(mat, geneinfo, by = "pr_gene_id")
                rownames(mat_all) <- mat_all$pr_gene_symbol
                mat_all <- mat_all[c(treatment_design$treatment_ID)]
                mat_all2 <- cbind(data.frame(NAME = rownames(mat_all)),
                              data.frame(DESCRIPTION = rep(NA, nrow(mat_all))),
                              mat_all)
                write.table(mat_all2, paste0(result_path, "drug_treated_expression_matrix/", dataset, "/", cancer, "/", drug, "/", cell, "_", treatment_time,"_", treatment_dose, ".txt"),
                            row.names = F, quote = F, sep = "\t")
                
                # # differential expression analysis
                # design = cbind(data.frame(drug = rep(1, nrow(treatment_design))),treatment_design[1] )
                # fit <- lmFit(mat_all, design)
                # fit <- eBayes(fit)
                # # topTable(fit, coef="treatment")
                # # DEgenes <- topTable(fit, coef="treatment", number=Inf,
                # #                     p.value=0.1, adjust="fdr")
                # allgenes <- topTable(fit, coef = "treatment", number = Inf)
                
                # write.csv(allgenes, paste0(drug_diff_res, "/allgenes/",dataset, "/", cancer, "/", drug, "/", cell, "_", treatment_time,"_", treatment_dose, ".csv"),
                #           quote = F)
                
                # mat_sig <- merge(mat, geneinfo_sig, by = "pr_gene_id")
                # rownames(mat_sig) <- mat_sig$pr_gene_symbol
                # mat_sig <- mat_sig[c(treatment_design$treatment_ID)]
                # mat_sig2 <- cbind(data.frame(NAME = rownames(mat_sig)),
                #               data.frame(DESCRIPTION = rep(NA, nrow(mat_sig))),
                #               mat_sig)
                # write.table(mat_sig2, paste0(result_path, "drug_treated_expression_978genes_matrix/", dataset, "/", cancer, "/", drug, "/", cell, "_", treatment_time,"_", treatment_dose, ".txt"),
                #             row.names = F, quote = F, sep = "\t")

                # # differential expression analysis
                # design = cbind(data.frame(drug = rep(1, nrow(treatment_design))),treatment_design[1] )
                # fit <- lmFit(mat_sig, design)
                # fit <- eBayes(fit)
                # # topTable(fit, coef="treatment")
                # # DEgenes <- topTable(fit, coef="treatment", number=Inf,
                # #                     p.value=0.1, adjust="fdr")
                # allgenes <- topTable(fit, coef = "treatment", number = Inf)
                
                # write.csv(allgenes, paste0(drug_diff_res, "/978genes/",dataset, "/", cancer, "/", drug, "/", cell, "_", treatment_time,"_", treatment_dose, ".csv"),
                #           quote = F)

              }
            }
          }
        }
      }
    }
  }
}else if(dataset == '92742'){
  # dataset = '92742'
  inst_infopath=path("Data/LINCS_data/", "GSE92742_Broad_LINCS_inst_info.txt")
  inst_GSE70138 <- read.csv(inst_infopath, header = T, sep = "\t")
  # inst_GSE92742 <- read.csv("LINCS_data/GSE92742/file/GSE92742_Broad_LINCS_inst_info.txt", header=T, sep="\t")

  inst_GSE92742[inst_GSE92742$pert_time == "24|4",]$pert_time <- "24"
  inst_GSE92742[inst_GSE92742$pert_time == "6|4",]$pert_time <- "6"
  inst_GSE92742 = inst_GSE92742[which(is.element(inst_GSE92742$pert_type, c("ctl_vehicle", "trt_cp"))), ]
  inst_GSE92742 = inst_GSE92742[-which(is.element(inst_GSE92742$pert_iname, c("PBS","H2O"))), ]

  inst_GSE92742 <- inner_join(drug_index, inst_GSE92742)
  
  # geneinfo <- read.csv("LINCS_data/GSE70138/file/GSE92742_Broad_LINCS_gene_info.txt", sep="\t", header = T)
  geneinfo_path=path("Data/LINCS_data/", "GSE92742_Broad_LINCS_gene_info.txt")
  geneinfo <- read.csv(geneinfo_path, sep = "\t", header = T)
  geneinfo_sig <- geneinfo[which(geneinfo$pr_is_lm == 1), ]

  # filter out cell line that is not in cellline_cancertype
  #   cellline_cancertype <- read.csv("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/s1_cellline_selection/TCGA_110CL/TCGA_CL_name.csv")
  # cell_annot_selected <- inner_join(cellline_cancertype, cell_annot)
  # cell_annot_selected <- cell_annot_selected[c('cell_id','TCGA_cancer_type')]
  # write.table(cell_annot_selected, paste0(result_path, 'cell_annot_selected_92742.txt'), sep = '\t', quote = F, row.names = F)
  cell_annot_selected_path=path("Data/LINCS_data/", "cell_annot_selected_92742.txt")
  cell_annot_selected <- read.table(cell_annot_selected_path, sep = '\t',header = T)
  # cell_annot_selected <- read.table(paste0(result_path, 'cell_annot_selected_92742.txt'), sep = '\t',header = T)
  inst_GSE92742 <- inner_join(cell_annot_selected, inst_GSE92742)
  inst_GSE92742_cancer <- inst_GSE92742[inst_GSE92742$TCGA_cancer_type == cancer, ]
  
  for (drug in unique(inst_GSE92742_cancer[-which(inst_GSE92742_cancer$pert_iname=='DMSO'),]$drug_index)) {

    dir.create(paste0(result_path, "drug_treated_expression_matrix/", dataset, "/", cancer, "/", drug))
    # dir.create(paste0(result_path, "drug_treated_expression_978genes_matrix/", dataset, "/", cancer, "/", drug))
    dir.create(paste0(result_path, "drug_treated_design/for_GSEA/",   dataset, "/", cancer, "/", drug))
    # dir.create(paste0(result_path, "drug_treated_design/for_GSVA/",   dataset, "/", cancer, "/", drug))
    
    dir.create(paste0(drug_diff_res, "/allgenes/", dataset, "/", cancer, "/", drug))
    # dir.create(paste0(drug_diff_res, "/978genes/", dataset, "/", cancer, "/", drug))

    inst_GSE92742_selected <- inst_GSE92742_cancer[inst_GSE92742_cancer$drug_index == drug, ]
    
    for(cell in unique(inst_GSE92742_selected$cell_id)){
      for (treatment_time in unique(inst_GSE92742$pert_time)) { 
        
        inst_GSE92742_selected2 <- inst_GSE92742_selected[inst_GSE92742_selected$pert_time == treatment_time &
                                                            inst_GSE92742_selected$cell_id == cell, ]
        inst_GSE92742_selected2$pert_dose <- as.numeric(inst_GSE92742_selected2$pert_dose)
        
        if (nrow(inst_GSE92742_selected2) != 0) {
          
          doses = unique(inst_GSE92742_selected2$pert_dose)
          #   doses = doses[-grep(-666, doses)]
          for (treatment_dose in doses) { # c(6, 24)
            
            inst_GSE92742_selected3 <- inst_GSE92742_selected2[inst_GSE92742_selected2$pert_dose == treatment_dose ,]
            
            if (nrow(inst_GSE92742_selected3) != 0) {
              
              DMSO <- inst_GSE92742_cancer[inst_GSE92742_cancer$pert_time == treatment_time &
                                             inst_GSE92742_cancer$cell_id == cell & 
                                             inst_GSE92742_cancer$pert_dose == -666, ]
              
              print(paste(cell, treatment_time, length(unique(inst_GSE92742_selected3$pert_iname))) )
              
              if (drug != "DMSO") {
                print(paste(cell, treatment_time, drug))
                
                d <- inst_GSE92742_selected3
                DMSO_match <- DMSO[DMSO$rna_plate %in% d$rna_plate, ]
                treatment_design <- unique(data.frame(
                  treatment = c(rep(1, length(d$distil_id)), rep(0, length(DMSO_match$distil_id))),
                  treatment_ID = c(d$distil_id, DMSO_match$distil_id),
                  treatment_dose = c(d$pert_dose, DMSO_match$pert_dose)
                ))
                
                row.names(treatment_design) <- treatment_design$treatment_ID
                
                if (length(unique(treatment_design$treatment)) > 1) {
                  gex_path=path("Data/LINCS_data/", "GSE92742_Broad_LINCS_Level3_INF_mlr12k_n1319138x12328.gctx")
                  drug_92742_3h_L3 <- parse.gctx(gex_path,
                                                 cid = treatment_design$treatment_ID )
                  
                  mat <- as.data.frame(drug_92742_3h_L3@mat)
                  mat$pr_gene_id <- rownames(mat)
                  mat_all <- merge(mat, geneinfo, by = "pr_gene_id")
                  rownames(mat_all) <- mat_all$pr_gene_symbol
                  mat_all <- mat_all[c(treatment_design$treatment_ID)]
                  
                  mat_all2 <- cbind(
                    data.frame(NAME = rownames(mat_all)),
                    data.frame(DESCRIPTION = rep(NA, nrow(mat_all))),
                    mat_all
                  )
                  
                  write.table(mat_all2, paste0(result_path, "drug_treated_expression_matrix/", dataset, "/", cancer, "/", drug, "/", cell, "_", treatment_time, "_", treatment_dose, ".txt"),
                              row.names = F, quote = F, sep = "\t")

                  # # differential expression analysis
                  # design = cbind(data.frame(drug = rep(1, nrow(treatment_design))),treatment_design[1] )
                  # fit <- lmFit(mat_all, design)
                  # fit <- eBayes(fit)
                  # # topTable(fit, coef="treatment")
                  # # DEgenes <- topTable(fit, coef="treatment", number=Inf,
                  # #                     p.value=0.1, adjust="fdr")
                  # allgenes <- topTable(fit, coef = "treatment", number = Inf)
                  
                  # write.csv(allgenes, paste0(drug_diff_res, "/allgenes/",dataset, "/", cancer, "/", drug, "/", cell, "_", treatment_time,"_", treatment_dose, ".csv"),
                  #           quote = F)

                  # mat_sig <- merge(mat, geneinfo_sig, by = "pr_gene_id")
                  # rownames(mat_sig) <- mat_sig$pr_gene_symbol
                  # mat_sig <- mat_sig[c(treatment_design$treatment_ID)]
                  
                  # mat_sig2 <- cbind(data.frame(NAME = rownames(mat_sig)),
                  #   data.frame(DESCRIPTION = rep(NA, nrow(mat_sig))),
                  #   mat_sig)
                  

                  
                  # write.table(mat_sig2, paste0(result_path, "drug_treated_expression_978genes_matrix/", dataset, "/", cancer, "/", drug, "/", cell, "_", treatment_time, "_", treatment_dose, ".txt"),
                  #             row.names = F, quote = F, sep = "\t")
                  # # differential expression analysis
                  # design = cbind(data.frame(drug = rep(1, nrow(treatment_design))),treatment_design[1] )
                  # fit <- lmFit(mat_sig, design)
                  # fit <- eBayes(fit)
                  # # topTable(fit, coef="treatment")
                  # # DEgenes <- topTable(fit, coef="treatment", number=Inf,
                  # #                     p.value=0.1, adjust="fdr")
                  # allgenes <- topTable(fit, coef = "treatment", number = Inf)
                  
                  # write.csv(allgenes, paste0(drug_diff_res, "/978genes/",dataset, "/", cancer, "/", drug, "/", cell, "_", treatment_time,"_", treatment_dose, ".csv"),
                  #           quote = F)


                  # write.table(treatment_design, paste0(result_path, "drug_treated_design/for_GSVA/", dataset, "/", cancer, "/", drug, "/", cell, "_", treatment_time, "_", treatment_dose, ".txt"),
                  #   row.names = T, quote = F, sep = "\t")
                  

                  # grp <- list(paste(nrow(treatment_design), 2, 1 , collapse = " "),
                  #             paste("#", "treated", "DMSO", collapse = " "),
                  #             paste(c(rep("treated", nrow(treatment_design[treatment_design$treatment == 1,])),
                  #                     rep("DMSO", nrow(treatment_design[treatment_design$treatment == 0,]))),
                  #                   collapse = " "))
                  # grp <- do.call(rbind, grp)
                  # write.table(grp, paste0(result_path, "drug_treated_design/for_GSEA/", dataset, "/", cancer, "/", drug, "/", cell, "_", treatment_time, "_", treatment_dose, ".cls"),
                  #             col.names = F, row.names = F, quote = F, sep = "\t")
                  
  
                  
                }
              }
            }
          }
        }
      }
    }
  }
}

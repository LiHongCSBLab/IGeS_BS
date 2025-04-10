print("start to prepare expression profile")
work_path <- "/picb/bigdata/project/FengFYM/DComboNetV2/drug_MoA_reanalysis/"
setwd(work_path)
result_path <- "/picb/bigdata/project/FengFYM/immune_targeted_combo/immunecell_DEG_cor/GSEA_GSVA/"

dir.create(paste0(result_path, "drug_treated_expression_matrix/"))
dir.create(paste0(result_path, "drug_treated_design/"))
dir.create(paste0(result_path, "drug_treated_design/for_GSEA/"))
dir.create(paste0(result_path, "drug_treated_design/for_GSVA/"))

library(prada)
library(rhdf5)
library(limma)

source("l1ktools-master/R/cmap/io.R")
options(stringsAsFactors = F)
options(quote = NULL)

dataset <- "92742"
# cell <- "HEPG2"
# treatment_time <- 6
inst_GSE92742 <- read.csv("LINCS_data/GSE92742/file/GSE92742_Broad_LINCS_inst_info.txt", header=T, sep="\t")

geneinfo <- read.csv("LINCS_data/GSE70138/file/GSE92742_Broad_LINCS_gene_info.txt", sep="\t", header = T)

# update the cellline_cancertype table
# update the cellline_cancertype table
cellline_cancertype <- read.csv("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/s1_cellline_selection/TCGA_110CL/TCGA_CL_name.csv")

# filter out cell line that is not in cellline_cancertype
cell_annot_selected <- inner_join(cellline_cancertype, cell_annot)
cell_annot_selected <- cell_annot_selected[c('cell_id','TCGA_cancer_type')]
write.table(cell_annot_selected, paste0(result_path, 'cell_annot_selected_92742.txt'), sep = '\t', quote = F, row.names = F)

# tmp = cell_annot[is.element(cell_annot$cell_id,setdiff(unique(gsub('.101','',gsub('.311','',union(inst_GSE70138$cell_id, inst_GSE92742$cell_id)))), cellline_cancertype$cellline)),]
# tmp[which(tmp$sample_type !='normal'),]
# filter out cell line that is not in cellline_cancertype


inst_GSE92742[inst_GSE92742$pert_time == "24|4",]$pert_time <- "24"
inst_GSE92742[inst_GSE92742$pert_time == "6|4",]$pert_time <- "6"

for (cell in unique(inst_GSE92742$cell_id)) {
  
  cancer <- unique(cell_annot_selected[cell_annot_selected$cell_id == cell, ]$TCGA_cancer_type)
  
  if (length(cancer) == 1) {
    
    dir.create(paste0(result_path, "drug_treated_expression_matrix/", cancer))
    dir.create(paste0(result_path, "drug_treated_design/for_GSEA/", cancer))
    dir.create(paste0(result_path, "drug_treated_design/for_GSVA/", cancer))
    
    # for (treatment_time in c("24", "6")) {
    for (treatment_time in unique(inst_GSE92742$pert_time)) { # c(6, 24)
      
      inst_GSE92742_selected <- inst_GSE92742[inst_GSE92742$pert_time == treatment_time & inst_GSE92742$cell_id == cell,]
      inst_GSE92742_selected$pert_dose <- as.numeric(inst_GSE92742_selected$pert_dose)
      
      if (nrow(inst_GSE70138_selected) != 0) {
        for (treatment_dose in unique(inst_GSE70138_selected$pert_dose)) { # c(6, 24)
          
          # inst_GSE92742_selected2 <- inst_GSE92742_selected[inst_GSE92742_selected$pert_dose >= 5 ,]
          inst_GSE92742_selected2 <- inst_GSE92742_selected[inst_GSE92742_selected$pert_dose == treatment_dose ,]
          
          if (nrow(inst_GSE92742_selected2) != 0) {
            
            DMSO <- inst_GSE92742[inst_GSE92742$pert_time == treatment_time & inst_GSE92742$cell_id == cell & inst_GSE92742$pert_dose==-666,]
            print(paste(cell, treatment_time, length(unique(inst_GSE92742_selected2$pert_iname))) )
            
            for (drug in unique(inst_GSE92742_selected2$pert_iname)) {
              if (drug != "DMSO") {
                print(paste(cell, treatment_time, drug))
                d <- inst_GSE92742_selected2[inst_GSE92742_selected2$pert_iname == drug, ]
                # d <- d[d$pert_dose >= 5, ]
                DMSO_match <- DMSO[DMSO$rna_plate %in% d$rna_plate, ]
                treatment_design <- data.frame(treatment = c(rep(1, length(d$distil_id)), rep(0, length(DMSO_match$distil_id))),
                                               treatment_ID =  c(d$distil_id, DMSO_match$distil_id ),
                                               treatment_dose = c(d$pert_dose, DMSO_match$pert_dose))
                
                row.names(treatment_design) <- treatment_design$treatment_ID
                
                if (length(unique(treatment_design$treatment)) > 1) {
                  drug_92742_3h_L3 <- parse.gctx("LINCS_data/GSE92742/GSE92742_Broad_LINCS_Level3_INF_mlr12k_n1319138x12328.gctx",
                                                 cid = treatment_design$treatment_ID )
                  
                  mat <- as.data.frame(drug_92742_3h_L3@mat)
                  mat$pr_gene_id <- rownames(mat)
                  mat <- merge(mat, geneinfo, by = "pr_gene_id")
                  rownames(mat) <- mat$pr_gene_symbol
                  mat <- mat[c(treatment_design$treatment_ID)]
                  
                  mat2 <- cbind(
                    data.frame(NAME = rownames(mat)),
                    data.frame(DESCRIPTION = rep(NA, nrow(mat))),
                    mat
                  )
                  
                  write.table(mat2, paste0(result_path, "drug_treated_expression_matrix/", cancer, "/", drug, "_", dataset, "_", cell, "_", treatment_time, "_", treatment_dose, ".txt"),
                              row.names = F, quote = F, sep = "\t")
                  write.table(treatment_design, paste0(result_path, "drug_treated_design/for_GSVA/", cancer, "/", drug, "_", dataset, "_", cell, "_", treatment_time, "_", treatment_dose, ".txt"),
                              row.names = T, quote = F, sep = "\t")
                  grp <- list(paste(nrow(treatment_design), 2, 1 , collapse = " "),
                              paste("#", "treated", "DMSO", collapse = " "),
                              paste(c(rep("treated", nrow(treatment_design[treatment_design$treatment == 1,])),
                                      rep("DMSO", nrow(treatment_design[treatment_design$treatment == 0,]))),
                                    collapse = " "))
                  grp <- do.call(rbind, grp)
                  write.table(grp, paste0(result_path, "drug_treated_design/for_GSEA/", cancer, "/", drug, "_", dataset, "_", cell, "_", treatment_time, "_", treatment_dose, ".cls"),
                              col.names = F, row.names = F, quote = F, sep = "\t"
                  )
                }
              }
            }
          }
        }
      }
    }
  }else if(length(cancer) > 1) {
    
    for(c in cancer) {
      dir.create(paste0(result_path, "drug_treated_expression_matrix/", c))
      dir.create(paste0(result_path, "drug_treated_design/for_GSEA/", c))
      dir.create(paste0(result_path, "drug_treated_design/for_GSVA/", c))
      
      # for (treatment_time in c("24", "6")) {
      for (treatment_time in unique(inst_GSE92742$pert_time)) { # c(6, 24)
        
        inst_GSE92742_selected <- inst_GSE92742[inst_GSE92742$pert_time == treatment_time & inst_GSE92742$cell_id == cell,]
        inst_GSE92742_selected$pert_dose <- as.numeric(inst_GSE92742_selected$pert_dose)
        
        if (nrow(inst_GSE70138_selected) != 0) {
          for (treatment_dose in unique(inst_GSE70138_selected$pert_dose)) { # c(6, 24)
            
            # inst_GSE92742_selected2 <- inst_GSE92742_selected[inst_GSE92742_selected$pert_dose >= 5 ,]
            inst_GSE92742_selected2 <- inst_GSE92742_selected[inst_GSE92742_selected$pert_dose == treatment_dose ,]
            
            if (nrow(inst_GSE92742_selected2) != 0) {
              
              
              DMSO <- inst_GSE92742[inst_GSE92742$pert_time == treatment_time & inst_GSE92742$cell_id == cell & inst_GSE92742$pert_dose==-666,]
              print(paste(cell, treatment_time, length(unique(inst_GSE92742_selected2$pert_iname))) )
              
              for(drug in unique(inst_GSE92742_selected2$pert_iname)) {
                
                if(drug != "DMSO") {
                  print(paste(cell, treatment_time, drug))
                  d <- inst_GSE92742_selected2[inst_GSE92742_selected2$pert_iname == drug,]
                  # d <- d[d$pert_dose >= 5, ]
                  DMSO_match <- DMSO[DMSO$rna_plate %in% d$rna_plate, ]
                  treatment_design <- data.frame(treatment = c(rep(1,length(d$distil_id)),rep(0,length(DMSO_match$distil_id))),
                                                 treatment_ID =  c(d$distil_id, DMSO_match$distil_id ),
                                                 treatment_dose = c(d$pert_dose, DMSO_match$pert_dose))
                  
                  row.names(treatment_design) <- treatment_design$treatment_ID
                  
                  if(length(unique(treatment_design$treatment))>1) {
                    drug_92742_3h_L3 <- parse.gctx("LINCS_data/GSE92742/GSE92742_Broad_LINCS_Level3_INF_mlr12k_n1319138x12328.gctx",
                                                   cid = treatment_design$treatment_ID )
                    
                    mat <- as.data.frame(drug_92742_3h_L3@mat)
                    mat$pr_gene_id <- rownames(mat)
                    mat <- merge(mat, geneinfo, by = "pr_gene_id")
                    rownames(mat) <- mat$pr_gene_symbol
                    mat <- mat[c(treatment_design$treatment_ID)]
                    
                    mat2 <- cbind(
                      data.frame(NAME = rownames(mat)),
                      data.frame(DESCRIPTION = rep(NA, nrow(mat))),
                      mat
                    )
                    
                    write.table(mat2, paste0(result_path, "drug_treated_expression_matrix/", c, "/", drug, "_", dataset, "_", cell, "_", treatment_time, "_", treatment_dose, ".txt"),
                                row.names = F, quote = F, sep = "\t")
                    write.table(treatment_design, paste0(result_path, "drug_treated_design/for_GSVA/", c, "/", drug, "_", dataset, "_", cell, "_", treatment_time, "_", treatment_dose, ".txt"),
                                row.names = T, quote = F, sep = "\t")
                    grp <- list(paste(nrow(treatment_design), 2, 1 , collapse = " "),
                                paste("#", "treated", "DMSO", collapse = " "),
                                paste(c(rep("treated", nrow(treatment_design[treatment_design$treatment == 1,])),
                                        rep("DMSO", nrow(treatment_design[treatment_design$treatment == 0,]))),
                                      collapse = " "))
                    grp <- do.call(rbind, grp)
                    write.table(grp, paste0(result_path, "drug_treated_design/for_GSEA/", c, "/", drug, "_", dataset, "_", cell, "_", treatment_time, "_", treatment_dose, ".cls"),
                                col.names = F, row.names = F, quote = F, sep = "\t")
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}

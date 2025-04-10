# TIL estimation with Immunedeconv ----------------------------------------
# TPM-normalized
# not log-transformed.

immunecell_estiamtor <- function(exprmat,
                                 cancer,
                                 sample_type,
                                 estimator_tool = c("ConsensusTME",
                                                    "immunedeconv",
                                                    "immucellai",
                                                    "immucellai_mouse"),
                                 method = c(
                                   "all",
                                   "xcell",
                                   "epic",
                                   "mcp_counter",
                                   "quantiseq",
                                   "timer",
                                   "ciberbort",
                                   "ciberbortABS",
                                   "ssgsea",
                                   "gsva",
                                   "plage",
                                   "zscore",
                                   "singScore"),
                                 savepath) {
  if (!require(ConsensusTME)) {
    devtools::install_github("cansysbio/ConsensusTME")
  }
  if (!require(immunedeconv)) {
    install.packages("remotes")
    remotes::install_github("grst/immunedeconv")
    # devtools::install_github("grst/immunedeconv")
  }
  if (!require(MCPcounter)) {
    devtools::install_github("ebecht/MCPcounter", ref = "master", subdir = "Source")
  }
  if (!require(xCell)) {
    devtools::install_github('dviraran/xCell')
  }
  
  
  res_ConsensusTME <- list()
  res_immunedeconv <- list()
  
  if (!is.element(
    estimator_tool,
    c(
      "ConsensusTME",
      "immunedeconv",
      "immucellai",
      "immucellai_mouse"
    )
  )) {
    print("ERROR: other immune cell estimators should be generated elsewhere! ")
  } else if (estimator_tool == "immucellai" |
             estimator_tool == "immucellai_mouse") {
    print("Warning: immucellai only provides webserver, ")
    print("PLZ use their website http://bioinfo.life.hust.edu.cn/ImmuCellAI#!/!")
    
  } else if (estimator_tool == "ConsensusTME") {
    if (is.null(cancer)) {
      print("Error: Cancer abbreviation must be provided for ConsensusTME")
      print("       Supported cancer types: ")
      print(ConsensusTME::cancerAll)
      
    } else{
      if (method == "all") {
        res_ConsensusTME[[1]] <-
          ConsensusTME::consensusTMEAnalysis(as.matrix(exprmat),
                                             cancer = cancer,
                                             statMethod = "ssgsea")
        res_ConsensusTME[[2]] <-
          ConsensusTME::consensusTMEAnalysis(as.matrix(exprmat),
                                             cancer = cancer,
                                             statMethod = "gsva")
        res_ConsensusTME[[3]] <-
          ConsensusTME::consensusTMEAnalysis(as.matrix(exprmat),
                                             cancer = cancer,
                                             statMethod = "plage")
        res_ConsensusTME[[4]] <-
          ConsensusTME::consensusTMEAnalysis(as.matrix(exprmat),
                                             cancer = cancer,
                                             statMethod = "zscore")
        res_ConsensusTME[[5]] <-
          ConsensusTME::consensusTMEAnalysis(as.matrix(exprmat),
                                             cancer = cancer,
                                             statMethod = "singScore")[[1]]
        
        rownames(res_ConsensusTME[[1]]) <- paste0(rownames(res_ConsensusTME[[1]]), '_ssgsea')
        rownames(res_ConsensusTME[[2]]) <- paste0(rownames(res_ConsensusTME[[2]]), '_gsva')
        rownames(res_ConsensusTME[[3]]) <- paste0(rownames(res_ConsensusTME[[3]]), '_plage')
        rownames(res_ConsensusTME[[4]]) <- paste0(rownames(res_ConsensusTME[[4]]), '_zscore')
        rownames(res_ConsensusTME[[5]]) <- paste0(rownames(res_ConsensusTME[[5]]), '_singScore')
        names(res_ConsensusTME) <-
          c("ssgsea", "gsva", "plage", "zscore", "singScore")
        save(
          res_ConsensusTME,
          file = paste0(
            savepath,'TIL_estimation/TCGA_',
            cancer,"_",sample_type,
            '_ConsensusTME_',
            method,
            '.Rdata'
          )
        )
        return(res_ConsensusTME)
      } else {
        res_ConsensusTME <- ConsensusTME::consensusTMEAnalysis(exprmat,
                                                               cancer = cancer,
                                                               statMethod = method)
        save(
          res_ConsensusTME,
          file = paste0(
            savepath, 'TIL_estimation/TCGA_',
            cancer,"_",sample_type,
            '_ConsensusTME_',
            method,
            '.Rdata'
          )
        )
        return(res_ConsensusTME)
        
      }
      # else if (is.element(method,
      #                 c("all", "ssgsea", "gsva", "plage", "zscore", "singScore"))) {
      #   print("Error: 'method' should be one of the following methods for ConsensusTME")
      #   print(c("all", "ssgsea", "gsva", "plage", "zscore", "singScore"))
      #   
      # } 
    }
    
  } else if (estimator_tool ==  "immunedeconv") {
    if (!is.element(
      method,
      c(
        "all",
        "ciberbort",
        "ciberbortABS",
        "xcell",
        "epic",
        "mcp_counter",
        "quantiseq"
      )
    )) {
      print("Error: 'method' should be one of the following methods for immunedeconv")
      print(
        c(
          "all",
          "ciberbort",
          "ciberbortABS",
          "xcell",
          "epic",
          "mcp_counter",
          "quantiseq"
        )
      )
      
    } else{
      if (method == "all") {
        # immunedeconv only calculate 36 cell types while xcell its own package can compute 64 types
        # instead of using immunedeconv, we choose to use xcell package
        
        # res_immunedeconv_xcell <-
        # as.data.frame(immunedeconv::deconvolute(exprmat, "xcell"))
        source("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/02_tumorGene_immuneSig_correlation/scripts/softwares/xCell-master/R/xCell.R")
        res_immunedeconv_xcell <- as.data.frame(xCellAnalysis(exprmat))
        
        res_immunedeconv_epic <-
          as.data.frame(immunedeconv::deconvolute(exprmat, "epic", tumor = T))
        # res_mcp_counter <- as.data.frame(immunedeconv::deconvolute(exprmat, "mcp_counter"))
        
        probesets = read.table(
          "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/02_tumorGene_immuneSig_correlation/scripts/softwares/mcp_count_need_files/probesets.txt",
          sep = "\t",
          stringsAsFactors = FALSE,
          colClasses = "character"
        )
        genes = read.table(
          "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/02_tumorGene_immuneSig_correlation/scripts/softwares/mcp_count_need_files/genes.txt",
          sep = "\t",
          stringsAsFactors = FALSE,
          header = TRUE,
          colClasses = "character",
          check.names = FALSE
        )
        
        res_immunedeconv_mcp_counter =   as.data.frame(
          MCPcounter::MCPcounter.estimate(
            exprmat,
            probesets = probesets,
            genes = genes,
            featuresType = "HUGO_symbols"
          )
        )
        
        res_immunedeconv_quantiseq <-
          as.data.frame(immunedeconv::deconvolute(exprmat, "quantiseq"))
        
        if(cancer == "LAML"){
          res_immunedeconv_timer <- matrix(NA, 6, ncol(exprmat))
          colnames(res_immunedeconv_timer) <- colnames(exprmat)
          res_immunedeconv_timer <- cbind(data.frame(cell_type = c("B cell",
                                                                   "T cell CD4+",
                                                                   "T cell CD8+",
                                                                   "Neutrophil",
                                                                   "Macrophage",
                                                                   "Myeloid dendritic cell"
                                                                   )),
                                            res_immunedeconv_timer)
        }else{
          res_immunedeconv_timer <-
           as.data.frame(immunedeconv::deconvolute(exprmat, "timer",
                                                   indications =
                                                   rep(cancer, ncol(exprmat))))

        }
        
        # rownames(res_immunedeconv_xcell) = paste0(res_immunedeconv_xcell$cell_type,'_xcell')
        rownames(res_immunedeconv_xcell) = paste0(rownames(res_immunedeconv_xcell),'_xcell')
        rownames(res_immunedeconv_epic) = paste0(res_immunedeconv_epic$cell_type,'_epic')
        # rownames(res_immunedeconv_epic) = paste0(rownames(res_immunedeconv_epic),'_epic')
        rownames(res_immunedeconv_mcp_counter) = paste0(rownames(res_immunedeconv_mcp_counter),'_mcp_counter')
        rownames(res_immunedeconv_quantiseq) = paste0(res_immunedeconv_quantiseq$cell_type,'_quantiseq')
        # rownames(res_immunedeconv_quantiseq) = paste0(rownames(res_immunedeconv_quantiseq),'_quantiseq')
        rownames(res_immunedeconv_timer) = paste0(res_immunedeconv_timer$cell_type ,'_timer')
        # rownames(res_immunedeconv_timer) = paste0(rownames(res_immunedeconv_timer) ,'_timer')

        # res_immunedeconv_xcell = res_immunedeconv_xcell[-which(names(res_immunedeconv_xcell) == 'cell_type')]
        res_immunedeconv_epic = res_immunedeconv_epic[-which(names(res_immunedeconv_epic) == 'cell_type')]
        res_immunedeconv_quantiseq = res_immunedeconv_quantiseq[-which(names(res_immunedeconv_quantiseq) == 'cell_type')]
        res_immunedeconv_timer = res_immunedeconv_timer[-which(names(res_immunedeconv_timer) == 'cell_type')]
        
        res_immunedeconv <- list()
        res_immunedeconv[[1]] <- res_immunedeconv_xcell
        res_immunedeconv[[2]] <- res_immunedeconv_epic
        res_immunedeconv[[3]] <- res_immunedeconv_mcp_counter
        res_immunedeconv[[4]] <- res_immunedeconv_quantiseq
        res_immunedeconv[[5]] <- res_immunedeconv_timer
        names(res_immunedeconv) <-
          c("xcell", "epic", "mcp_counter", "quantiseq", "timer")
        dir.create(paste0(
          savepath, "TIL_estimation/"
        ))

        save(
          res_immunedeconv,
          file = paste0(
            savepath,'TIL_estimation/TCGA_',
            cancer,"_",sample_type,
            '_immunedeconv_',
            method,
            '.Rdata'
          )
        )
        return(res_immunedeconv)
        
      } else if (method == "cibersort" | method == "ciberbortABS") {
        print("Warning: unable to download CIBERSORT R script, ")
        print("PLZ use their webserver on https://cibersort.stanford.edu/index.php!")
        # cibersort
        # source('G:/lab/Projects/p2_immunotherapy_targetedtherapy/immune_clustering/scripts/tidybulk-master/R/CIBERSORT.R')
        # results <- CIBERSORT('G:/lab/Projects/p2_immunotherapy_targetedtherapy/immune_clustering/scripts/CIBERSORT-master/CIBERSORT_data/LM22.csv',
        #                      "mixture_file.txt",
        #                      perm=100, QN=TRUE)
        
      } else if (method == 'mcp_counter') {
        # res_mcp_counter <- as.data.frame(immunedeconv::deconvolute(exprmat, "mcp_counter"))
        
        probesets = read.table(
          "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/02_tumorGene_immuneSig_correlation/scripts/softwares/mcp_count_need_files/probesets.txt",
          sep = "\t",
          stringsAsFactors = FALSE,
          colClasses = "character"
        )
        genes = read.table(
          "/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/02_tumorGene_immuneSig_correlation/scripts/softwares/mcp_count_need_files/genes.txt",
          sep = "\t",
          stringsAsFactors = FALSE,
          header = TRUE,
          colClasses = "character",
          check.names = FALSE
        )
        
        res_immunedeconv =   as.data.frame(
          MCPcounter::MCPcounter.estimate(
            exprmat,
            probesets = probesets,
            genes = genes,
            featuresType = "HUGO_symbols"
          )
        )
        save(
          res_immunedeconv,
          file = paste0(
            savepath,'TIL_estimation/TCGA_',
            cancer, "_",sample_type,
            '_immunedeconv_',
            method,
            '.Rdata'
          )
        )
        return(res_immunedeconv)
        
      } else if (method == 'timer') {
        res_immunedeconv <-
          as.data.frame(immunedeconv::deconvolute(
            exprmat,
            method = method,
            indications =
              rep(cancer, ncol(exprmat))
          ))
        save(
          res_immunedeconv,
          paste0(
            savepath,'TIL_estimation/TCGA_',
            cancer,"_",sample_type,
            '_immunedeconv_',
            method,
            '.Rdata'
          )
        )
        return(res_immunedeconv)
        
      } else{
        res_immunedeconv[[1]] <-
          as.data.frame(immunedeconv::deconvolute(exprmat, method = method))
        names(res_immunedeconv) = method
        save(
          res_immunedeconv,
          file = paste0(
            savepath,'TIL_estimation/TCGA_',
            cancer,"_",sample_type,
            '_immunedeconv_',
            method,
            '.Rdata'
          )
        )
        return(res_immunedeconv)
        
      }
    }
  }
}
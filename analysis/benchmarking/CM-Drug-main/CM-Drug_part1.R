#### If you have any questions or are interested in collaborating, please feel free to contact us at hsiayun@foxmail.com. Let's work together to advance the field of immunotherapy.
#### Any use of this code should refer to the article: https://doi.org/10.1038/s41590-024-01789-x (Nature Immunology)
#### Code of CM-Drug

#Core & Minor gene sets
# setwd("workingDirectory/")
setwd("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/07_plotting_v2/CM_drug_res/")


list_CM_gene_set=list()

list_CM_gene_set[[1]] <- 
  c("CD8A","HLA-DRA","HLA-DPA1","HLA-DQB1","CD74","IFNG","HLA-DRB5","HLA-DPB1","HLA-DQA1","HLA-DMA","HLA-DRB1","KLRC1","HLA-DQA2",
    "CTSS","KLRD1","HLA-DOA","CIITA","CD8B","KIR2DL4","KLRC3","CD1B","HLA-DMB","CTSB","CD209","CD4","HLA-B",   
    "MICB","HLA-C","B2M","HLA-F","HLA-E","HLA-A","KLRC2")
names(list_CM_gene_set)[1] <- c("Antigen_Processing_and_Presentation")

list_CM_gene_set[[2]] <- 
  c("SH2D1A","IFNG","GZMB","LCP2","PTPN6","PIK3CD","ITGB2","LCK","FASLG","KLRD1","PRF1","HLA-E","CD247","CD48",   
    "CD244","KIR2DL4","KLRC3","IFNB1","KLRC1","MICB","SH2D1B","FCER1G","HCST","HLA-B","HLA-C","HLA-A","TYROBP","KLRC2")
names(list_CM_gene_set)[2] <- c("NaturalKiller_Cell_Cytotoxicity")

list_CM_gene_set[[3]] <- 
  c("K","CD3G","CD3D","ICOS","GRAP2","IFNG","LCP2","CD3E","PTPN6","PIK3CD","NFKBIA","CARD11","CD247","ITK","PTPRC","PDCD1",
    "CTLA4","PIK3R5","IL2","CD8B","CD8A","CD28","CD4","NFKBIE")
names(list_CM_gene_set)[3] <- c("TCR_Signaling_Pathway")

list_CM_gene_set[[4]] <- 
  c("SIRPG","GPR171","CRTAM","GZMA","LAG3","CTSW","PRF1","NKG7","CCR5","C1QB","GZMB","GZMH","LY9","CD7","LAX1","IL7R",
    "ITK","IL2RB","LCP2","KLRG1","SELL","CD8B","CD8A","GNLY")
names(list_CM_gene_set)[4] <- c("Cytotoxiclty_of_ImmuCellAI")

list_CM_gene_set[[5]] <- unique(c(list_CM_gene_set[[1]],list_CM_gene_set[[2]],list_CM_gene_set[[3]],list_CM_gene_set[[4]]))

names(list_CM_gene_set)[5] <- c("Core_gene")

list_CM_gene_set[[6]] <- 
  c("CXCL11","APOBEC3A","CXCL9","IL2","ORM2","ISG20","TCHHL1","OASL","MX1","IL27","TNF","OAS1","CCL8",
    "ISG15","BST2","MX2","CCL1","TLR7","MUC4","ORM1","DEFA3","CCR3","CHIT1","REG3G","C8G","CD40LG",
    "CCR7","MPO","IDO1","IL22","CXCL10","CCL5","PAEP","GNLY","CD8A","CXCL13","CCL4","CCL7","IFNG",
    "PDCD1","CCL4L2","CCL3L3","CCL3","FASLG","FGR","APOBEC3H","HLA-B","MMP12","TLR8","APOBEC3G","IRF7","SLC29A3",
    "HAMP","CD40","NOD2","CXCL2","CTSS","B2M","IRF5","IL10","MARCO","BPI","HCK","CYBB","IL6", 
    "CXCR1","DEFA4","DEFA1","S100A12","CAMP","PGLYRP1","CCL13","AZU1")
names(list_CM_gene_set)[6] <- c("Antimicrobials")

list_CM_gene_set[[7]] <- 
  c("CD79A","CD19","CD79B","CARD11","CD72","INPP5D","PLCG2","PIK3CG","PTPN6","PRKCB","RAC2","VAV1","BTK","SYK","PIK3CD","PIK3R5",
    "NFKBIA","LYN","CR2","LILRB3","NFKBIE","JUN")
names(list_CM_gene_set)[7] <- c("BCR_Signaling_Pathway")

list_CM_gene_set[[8]] <- unique(c(list_CM_gene_set[[6]],list_CM_gene_set[[7]]))

names(list_CM_gene_set)[8] <- c("Minor_gene")

#we set the name of CM gene sets as the "super" initially. So in the code, for consistency, we keep the name 
super.third.human.pd1.all <- list_CM_gene_set

saveRDS(super.third.human.pd1.all,"./super.third.human.pd1.all.RDS")

#Function definition
.systemInfo = function(eachThreadDir, spid, type, idx) {
  tempInfo = sprintf('%s%s.%s.', eachThreadDir, spid, type)
  return(paste0(tempInfo, idx, '.trds'))
}

function.remove = function(paths, check = FALSE) {
  stopifnot(is.character(paths))
  stopifnot(is.logical(check))
  
  if (check) {
    existsVector = file.exists(paths) | dir.exists(paths)
    invisible(suppressWarnings(file.remove(paths[existsVector])))
  } else {
    invisible(suppressWarnings(file.remove(paths)))
  }
}

function.doGC = function() {
  # silence garbage collection (gc)
  invisible(gc())
}

function.getPID = function() {
  # Get the process
  NODE = unlist(strsplit(Sys.info()['nodename'], '\\.'))[1]
  PID = Sys.getpid()
  tempString = sprintf('%s-%s', NODE, PID)
  return(tempString)
}


#temporary directory
function.getTempDir = function(useRAM = TRUE, usecache = TRUE) {
  
  currentNode = unlist(strsplit(Sys.info()['nodename'], '\\.'), use.names = FALSE)[1]
  
  
  if (usecache) {
    tempDirPath = '/'
    return(tempDirPath)
  }
  
  if (useRAM) {
    tempDirPath = '/the_temp_path/'
    return(tempDirPath)
    
  } else {
    
    tempDirPath = '/'
    return(tempDirPath)
  }
}

function.freadXZ = function(file, sep = '\t', header = TRUE) {
  stopifnot(is.character(file))
  require(data.table)
  
  tempOutDir = sprintf('%sfunction.freadXZ_%s/', function.getTempDir(useRAM = TRUE, usecache = TRUE), function.getPID())
  dir.create(tempOutDir, showWarnings = FALSE, recursive = TRUE)
  
  secondStamp = format(Sys.time(), format = '%s')
  fileID = .FileNameScramble(file)
  tempFilePath = sprintf('%s%s_%s.TMP', tempOutDir, secondStamp, fileID)
  xzCommand = sprintf('xz -d -c %s > %s', file, tempFilePath)
  system(command = xzCommand, ignore.stdout = FALSE, ignore.stderr = TRUE)

  tempDT = fread(file = tempFilePath, sep = sep, header = header, showProgress = FALSE)
  
  file.remove(tempFilePath)
  
  return(tempDT)
}

.FileNameScramble = function(charVec) {
  require(stringi)
  
  charVec = stri_trans_tolower(charVec)
  
  charVec = tail(unlist(strsplit(charVec, split = '/'), use.names = FALSE), n = 1)
  
  charVec = unlist(strsplit(charVec, split = ''), use.names = FALSE)
  
  charVec = paste0(na.omit(match(charVec, letters)), collapse = '')
  
  return(charVec)
}


function.XZSaveRDS = function(obj, file, threads = 32, compression = 6) {
  
  stopifnot(is.character(file))
  
  function.remove(paths = file, check = TRUE)
  
  xzCommand = sprintf('xz -z -T %s -%s > %s', threads, compression, file)
  
  xzConnection = pipe(description = xzCommand, open = 'wb')
  
  saveRDS(object = obj, file = xzConnection)
  
  close(xzConnection)
}

################################################################################
################################################################################


require(data.table)
require(stringi)
require(parallel)

# Declaring Global Variables
inDir = c('Data'='Data')

options(
  stringsAsFactors = FALSE,
  warn = 1
)

# Global setting
# Output
outDir = './Data/Compound_Data/'
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)

# Custom Functions
SelectionResolver = function(conditionDF, controlDF) {
  # Output Resolving Function
  finalSampleSize = min(conditionDF$sampleCount, controlDF$sampleCount, 20)
  tempVector = c(
    'source_dataset' = conditionDF$source_dataset,
    'rna_centre' = conditionDF$rna_centre,
    'cell_id' = conditionDF$cell_id,
    
    'trt_type' = conditionDF$pert_type,
    'trt_cval' = conditionDF$pert_cval,
    'trt_iname' = conditionDF$pert_iname,
    'trt_cdose' = conditionDF$pert_cdose,
    'trt_ctime' = conditionDF$pert_ctime,
    'trt_orig_sampleCount' = conditionDF$sampleCount,
    'trt_final_sampleCount' = finalSampleSize,
    
    'ctl_type' = controlDF$pert_type,
    'ctl_cval' = controlDF$pert_cval,
    'ctl_iname' = controlDF$pert_iname,
    'ctl_cdose' = controlDF$pert_cdose,
    'ctl_ctime' = controlDF$pert_ctime,
    'ctl_orig_sampleCount' = controlDF$sampleCount,
    'ctl_final_sampleCount' = finalSampleSize
  )
  return(tempVector)
}

#Matched the Sample and Control
finalSampleMetadata = NULL
finalMatchedMetadata = NULL

for (indext_iii in 1) {
  case_jjj = names(inDir)[indext_iii]
  dir_kkk = inDir[indext_iii]
  
  cat(sprintf('\nLoading Metadata - %s.\n', case_jjj))
  tempPath = sprintf('pertInfo.%s', dir_kkk)
  metadata = fread(input = tempPath, sep = '\t', header = TRUE, data.table = FALSE)
  
  metadata$pert_cdose = sprintf(
    '%s%s',
    metadata$pert_dose,
    metadata$pert_dose_unit
  )
  metadata$pert_cdose[grepl('^NA', metadata$pert_cdose)] = ''
  
  metadata$pert_ctime = sprintf(
    '%s%s',
    metadata$pert_time,
    metadata$pert_time_unit
  )
  
  metadata$pert_cval = sprintf(
    '%s %s for %s in %s',
    metadata$pert_cdose,
    metadata$pert_iname,
    metadata$pert_ctime,
    metadata$cell_id
  )
  metadata$pert_cval = stri_trim_left(metadata$pert_cval)
  
  metadata$rna_centre = sapply(strsplit(x = as.character(metadata$rna_plate), split = '_'), 
                               USE.NAMES = FALSE, FUN = function(i) {return(i[[1]])})
  
  selectionVector = c(
    'Num',
    'rna_plate',
    'rna_centre',
    'pert_iname',
    'pert_type',
    'pert_dose',
    'pert_dose_unit',
    'pert_time',
    'pert_cdose',
    'pert_ctime',
    'pert_cval',
    'cell_id'
  )
  sampleMetadata = metadata[, selectionVector]
  sampleMetadata$source_dataset = rep(case_jjj, nrow(metadata))
  
  # Variable Cleanup
  rm(selectionVector, metadata)
  
  # ----- Generate Baseline-Matched Condition Metadata -----
  
  cat('Generating Baseline-Matched Condition Metadata.\n')
  
  # Match Contrast-Baseline
  selectionVector = c(
    'source_dataset',
    'pert_cval',
    'rna_centre',
    'pert_type',
    'cell_id',
    'pert_iname',
    'pert_ctime',
    'pert_cdose'
  )
  useDF = unique(sampleMetadata[, selectionVector])
  
  #Sample Counts
  sampleMetadata = as.data.table(sampleMetadata)
  setkey(sampleMetadata, rna_centre, pert_type, pert_cval)
  useDF$sampleCount = unlist(mclapply(1:nrow(useDF), mc.preschedule = TRUE, mc.cores = 160, mc.cleanup = TRUE, FUN = function(i) {
    tempCondition = useDF[i, ]
    return(sampleMetadata[.(tempCondition$rna_centre, tempCondition$pert_type, tempCondition$pert_cval), .N, nomatch = 0])
  }), recursive = FALSE, use.names = FALSE)
  sampleMetadata = as.data.frame(sampleMetadata)
  
  # Filter out Sample Size == 1
  useDF = subset(useDF, sampleCount > 1)
  
  # Prepare Subsets
  casesDF = subset(useDF, grepl('^trt', useDF$pert_type))
  controlsDF = subset(useDF, grepl('^ctl', useDF$pert_type))
  
  tempList = mclapply(1:nrow(casesDF), mc.preschedule = TRUE, mc.cores = 32, mc.cleanup = TRUE, FUN = function(currentRowIndex) {
    condition_aaa = casesDF[currentRowIndex, ]
    centre_bbb = condition_aaa$rna_centre
    type_ccc = condition_aaa$pert_type
    cell_ddd = condition_aaa$cell_id
    time_eee = condition_aaa$pert_ctime
    dose_fff = condition_aaa$pert_cdose
    
    controlSubset = subset(controlsDF, rna_centre == centre_bbb & cell_id == cell_ddd & pert_ctime == time_eee)
    
    # Compound Data Resolving
    if (type_ccc == 'trt_cp') {
      controlSubset = subset(controlSubset, pert_type %in% c('ctl_vehicle', 'ctl_untrt'))
      
      # SKIP if No Matching: 
      if (nrow(controlSubset) == 0) {
        return(NULL)
      }
      
      if ('DMSO' %in% controlSubset$pert_iname) {
     
        controlSubset = subset(controlSubset, pert_iname == 'DMSO')
        finalOutput = controlSubset[order(controlSubset$sampleCount, decreasing = TRUE), ][1, ]
        return(SelectionResolver(conditionDF = condition_aaa, controlDF = finalOutput))
      }
      
      if (any(c('PBS', 'H2O', 'UnTrt') %in% controlSubset$pert_iname)) {
        
        finalOutput = controlSubset[order(controlSubset$sampleCount, decreasing = TRUE), ][1, ]
        return(SelectionResolver(conditionDF = condition_aaa, controlDF = finalOutput))
        
      } else {
        finalOutput = controlSubset[order(controlSubset$sampleCount, decreasing = TRUE), ][1, ]
        return(SelectionResolver(conditionDF = condition_aaa, controlDF = finalOutput))
      }
    }
    
  })
  
  # Condition Metadata
  cat('Generateing Condition Metadata.\n')
  matchedMetadata = as.data.frame(do.call('rbind', tempList))
  
  # Final Metadata
  finalSampleMetadata = rbind(finalSampleMetadata, sampleMetadata)
  finalMatchedMetadata = rbind(finalMatchedMetadata, matchedMetadata)
  
  # Variable Cleanup
  rm(case_jjj, dir_kkk)
  rm(selectionVector, useDF, casesDF, controlsDF, tempList)
  rm(sampleMetadata, matchedMetadata)
}

# Appending Condition ID
finalMatchedMetadata = data.frame(
  case_ID = paste0('CM.', 1:nrow(finalMatchedMetadata)),
  finalMatchedMetadata
)

#Switching type
finalMatchedMetadata$trt_orig_sampleCount = as.numeric(finalMatchedMetadata$trt_orig_sampleCount)
finalMatchedMetadata$trt_final_sampleCount = as.numeric(finalMatchedMetadata$trt_final_sampleCount)

finalMatchedMetadata$ctl_orig_sampleCount = as.numeric(finalMatchedMetadata$ctl_orig_sampleCount)
finalMatchedMetadata$ctl_final_sampleCount = as.numeric(finalMatchedMetadata$ctl_final_sampleCount)

#Saving
cat('Now saving...\n')

#Sample Metadata
tempPath = sprintf('%sCM_Drug.sample', outDir)
function.XZSaveRDS(obj = finalSampleMetadata, file = tempPath)

#Matched Metadata
tempPath = sprintf('%sCM_Drug.condition', outDir)
function.XZSaveRDS(obj = finalMatchedMetadata, file = tempPath)

# cat(sprintf('START TIME: %s\n', .startTime))
# START TIME: Fri Apr  1 10:08:28 2022
# cat(sprintf('END TIME: %s\n\n', date()))
# END TIME: Fri Apr  1 20:23:49 2022

# Load Libraries
require(cmapR)
require(data.table)
require(foreach)
require(doParallel)

# Declaring Global Variables
inPath = c('Data'='Data_CM_Drug.gctx'
)
metadataDir = './Data/Compound_Data/'

options(
  stringsAsFactors = FALSE,
  warn = 1
)

#Preparation 
################################################################################

# Creation of Output Directory
outDir = 'Data/Compound_Data/'
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)

# To increase the read and write speed of the results, this portion is output to the memory disk
ramDir = 'ramdisk/'
dir.create(ramDir, recursive = TRUE, showWarnings = FALSE)

################################################################################
cat('Loading Metadata.\n')
tempPath = sprintf('%sCM_Drug.sample', metadataDir)
sampleDF = readRDS(tempPath)

tempPath = sprintf('%sCM_Drug.condition', metadataDir)
matchedDF = readRDS(tempPath)

# Sample Metadata Trimming
uniqueCVals = unique(c(matchedDF$trt_cval, matchedDF$ctl_cval))
sampleDF = subset(sampleDF, pert_cval %in% uniqueCVals)
rm(uniqueCVals)

# Enable data.table optimization
sampleDF = as.data.table(sampleDF)
setkey(sampleDF, rna_centre, pert_type, pert_cval)

cat('Loading Data.\n')
referenceRow = cmapR::read_gctx_ids(gctx_path = inPath['Data'],
                                    #dimension = 'row'
)
dataList = lapply(inPath, function(tempPath) {
  tempMatrix = parse_gctx(fname = tempPath)@mat[referenceRow, ]
  function.doGC()
  return(tempMatrix)
})

gctxMatrix = dataList[[1]]
rm(dataList)
function.doGC()

gctxMatrix = gctxMatrix[, sampleDF$Num]
function.doGC()

# Calculating
cat('Calculating...\n')

# Initiate Multi-Thread
.startTime = date()
SPID_mmm = function.getPID()
referenceColumn = matchedDF$case_ID
caseCount = length(matchedDF$case_ID)
geneCount = length(referenceRow)
registerDoParallel(cores = 40)

# Cleaning Cache
allPath = c(
  .systemInfo(eachThreadDir = ramDir, spid = SPID_mmm, type = 'CM_Drug.data', idx = 1:caseCount)
)
function.remove(paths = allPath, check = TRUE)

# Parallel Processing
cat('Begin Calculating Job: ')
tempOutput = foreach(i = 1:caseCount, .inorder = TRUE) %dopar% {
  # Status Update
  if (i %% 500 == 0) {
    cat(i, ' . ', sep = '')
  }
  
  # Obtain Matching-Sample Mapping
  task_nnn = matchedDF[i, ]
  trtSamples = sampleDF[.(task_nnn$rna_centre, task_nnn$trt_type, task_nnn$trt_cval), Num]
  ctlSamples = sampleDF[.(task_nnn$rna_centre, task_nnn$ctl_type, task_nnn$ctl_cval), Num]
  
  # Sample-Size Control
  if (task_nnn$trt_orig_sampleCount != task_nnn$trt_final_sampleCount) {
    tempLogical = sample(x = 1:task_nnn$trt_orig_sampleCount, size = task_nnn$trt_orig_sampleCount, replace = FALSE)
    tempLogical = (tempLogical %in% 1:task_nnn$trt_final_sampleCount)
    trtSamples = trtSamples[tempLogical]
    rm(tempLogical)
  }
  
  if (task_nnn$ctl_orig_sampleCount != task_nnn$ctl_final_sampleCount) {
    tempLogical = sample(x = 1:task_nnn$ctl_orig_sampleCount, size = task_nnn$ctl_orig_sampleCount, replace = FALSE)
    tempLogical = (tempLogical %in% 1:task_nnn$ctl_final_sampleCount)
    ctlSamples = ctlSamples[tempLogical]
    rm(tempLogical)
  }
  
  trtData = gctxMatrix[, trtSamples]
  ctlData = gctxMatrix[, ctlSamples]
  
  tempFC = rowMeans(trtData) - rowMeans(ctlData)
  
  # Save File
  tempPath = .systemInfo(eachThreadDir = ramDir, spid = SPID_mmm, type = 'CM_Drug.data', idx = i)
  function.XZSaveRDS(obj = tempFC, file = tempPath)
  
  # Variable Cleanup
  rm(task_nnn, trtSamples, ctlSamples)
  rm(trtData, ctlData)
  rm(tempResult, tempFC)
  function.doGC()
  return(NULL)
  #it takes time to run
}

# stop the multi-threads
registerDoSEQ()

#clean-up the variable 
rm(matchedDF, sampleDF, gctxMatrix, allPath, tempOutput)
function.doGC()
cat('clean-up have been done.\n')

tempMatrix = foreach(i = 1:caseCount, .inorder = TRUE, .combine = cbind, .maxcombine = 1000) %do% {
  tempPath = .systemInfo(eachThreadDir = ramDir, spid = SPID_mmm, type = 'CM_Drug.data', idx = i)
  return(readRDS(tempPath))
}
colnames(tempMatrix) = referenceColumn
rownames(tempMatrix) = referenceRow

tempPath = sprintf('%sCM_Drug.data', outDir)
function.XZSaveRDS(obj = tempMatrix, file = tempPath)
rm(tempMatrix)
function.doGC()
function.remove(paths = .systemInfo(eachThreadDir = ramDir, spid = SPID_mmm, type = 'CM_Drug.data', idx = 1:caseCount))

# Directory Cleanup
function.remove(ramDir)

# Print Timestamp
cat(sprintf('START TIME: %s\n', .startTime))
cat(sprintf('END TIME: %s\n\n', date()))

################################################################################
################################################################################
library(dplyr)
library(magrittr)
CM_Drug.data <- readRDS("./Data/Compound_Data/CM_Drug.data")
CM_Drug.data.gtc<- new("GCT", mat=CM_Drug.data)
write_gctx(CM_Drug.data.gtc, 
           compression_level = 9,
           "./Data/Compound_Data/CM_Drug.data.gctx",
           appenddim = FALSE)

rm(list = ls())
# setwd("~/tools/repository/CM-Drug")
CM_Drug.condition <- readRDS("./Data/Compound_Data/CM_Drug.condition")

data_trt_cp <- list()
data_trt_cp.df <- list()
fgsea.sam.trt_cp <- list()
fgsea.res.trt_cp <- list()
trt_cp_number.group <- list()

trt_cp_number <- which(CM_Drug.condition$trt_type == "trt_cp")

for(i in 1:4){trt_cp_number.group[[i]] <- trt_cp_number[((i-1)*10000+1):(i*10000)]}

trt_cp_number.group[[5]] <- trt_cp_number[40001:length(trt_cp_number)]

template <- parse_gctx("Data/Compound_Data/CM_Drug.data.gctx",
                       rid=1:12328, cid=1:10)@mat %>% as.data.frame()

template <- template %>% dplyr::mutate(gene_id=rownames(template ))

gene_df <- read.delim("./gene_df")

template$gene_id <- as.integer(template$gene_id)

template <- dplyr::left_join(template,gene_df,by= "gene_id")

saveRDS(template,"template.RDS")
save.image("template.RData")

################################################################################
################################################################################

# #Use fgsea to perform GSEA 
# options(
#   stringsAsFactors = FALSE,
#   warn = 0
# )
# library(fgsea)
# library(cmapR)
# library(tidyverse)
# setwd("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r1_drug_immuneSig_correlation_analysis/07_plotting_v2/CM_drug_res/")

# # setwd("~/tools/repository/CM-Drug")
# resultDir="./result/cp/"
# dir.create(resultDir, recursive = TRUE, showWarnings = FALSE)
# load("./template.RData")

# #library(tictoc)
# library(furrr)
# library(future)
# plan(multisession, workers = 40)

# .startTime = date()

# for(j in c(1:5)){
#   # setwd("~/tools/repository/CM-Drug")
#   load("./template.RData")
  
#   fun_cmap_fgsea <- function(x){
#     x <- x %>% unlist() 
#     names(x) <-template$gene_symbol
#     fgsea.res <- fgsea(pathways = genesets, stats = x,eps= 0.0, minSize  = 5, maxSize  = 500)
#     return(fgsea.res)
#   }
#   genesets <- readRDS("super.third.human.pd1.all.RDS")
#   data_trt_cp[[j]] <- parse_gctx("./Data/Compound_Data/CM_Drug.data.gctx", 
#                                  rid=1:12328, cid=trt_cp_number.group[[j]])
#   data_trt_cp.df[[j]] <- as.data.frame(data_trt_cp[[j]]@mat)
#   fgsea.res.trt_cp[[j]] <- furrr::future_map(data_trt_cp.df[[j]], ~ fun_cmap_fgsea(.x))
#   setwd("./result/cp/");saveRDS(fgsea.res.trt_cp[[j]],paste("c",j,"_fgsea.res.trt_cp_super.RDS",sep = ""));rm(list=ls());
# }

# cat(sprintf('START TIME: %s\n', .startTime))
# cat(sprintf('END TIME: %s\n\n', date()))

################################################################################
################################################################################

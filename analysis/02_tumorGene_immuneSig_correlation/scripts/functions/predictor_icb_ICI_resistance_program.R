
# ICI resistance program --------------------------------------------------

## Read expression data from tab-delimited text file, with official human
## gene symbols (HGNC) in the first columns
## Note: Input gene expression in TPM/FPKM/RPKM, and the code below will do log2
## transformation (log2(norm_value+1)) automatically.

model_icb_ICI_resistance_program <- function(gene_expression,
                                             cancer_type,
                                             sample_type,
                                             sig_dir,
                                             save_path) {
  # setwd('/picb/bigdata/project/FengFYM/immune_targeted_combo/immunecell_DEG_cor/')
  
  #source('Code/ImmRes_source.R')
  
  gene_expression = log2(gene_expression+1)
  
  compute.samples.res.scores <- function(r,
                                         res.sig,
                                         cell.sig,
                                         residu.flag = F,
                                         cc.sig = NULL,
                                         num.rounds = 1000) {
    r$res.ori <- NULL
    r$res <- NULL
    r$tme <- NULL
    r$X <- NULL
    names(res.sig) <- gsub(" ", ".", names(res.sig))
    two.sided <-
      unique(gsub(".up", "", gsub(".down", "", names(res.sig))))
    b <- is.element(paste0(two.sided, ".up"), names(res.sig)) &
      is.element(paste0(two.sided, ".down"), names(res.sig))
    two.sided <- two.sided[b]
    
    # r: TCGA RNAseq-expression profile (log2(TPM + 1))
    # res.sig: gene signatures, lists, list names are the name of signautres,
    #          every list is a gene list with gens names corresponded to signatures
    
    r$res <-
      get.OE.bulk(r,
                  res.sig,
                  num.rounds = num.rounds)
    
    res2 <-
      r$res[, paste0(two.sided, ".up")] - r$res[, paste0(two.sided, ".down")]
    if (!is.matrix(res2)) {
      res2 <- as.matrix(res2)
    }
    colnames(res2) <- two.sided
    r$res <- cbind(res2, r$res)
    rownames(r$res) <- colnames(r$tpm)
    
    
    r$tme <- get.OE.bulk(r,
                         cell.sig,
                         num.rounds = num.rounds)
    
    r$res <- cmb.res.scores(r)
    
    if (residu.flag) {
      print("Computing residuales")
      r$res <- t(get.residuals(t(r$res), colSums(r$tpm > 0)))
      r$tme <- t(get.residuals(t(r$tme), colSums(r$tpm > 0)))
    }
    
    if (!is.null(cc.sig)) {
      print("Filtering cell-cycle effects")
      if (is.null(r$cc.scores)) {
        r$cc.scores <- get.OE.bulk(r, cc.sig, num.rounds = num.rounds)
      }
      r$res.ori <- r$res
      r$res <- t(get.residuals(t(r$res), r$cc.scores))
    }
    if (is.element("resF", colnames(r$res)) &
        is.element("T.CD8", colnames(r$tme))) {
      r$res <-
        cbind.data.frame(r$res, resF.minus.TIL = r$res[, "resF"] - r$tme[, "T.CD8"])
      if (!is.null(r$res.ori)) {
        r$res.ori <-
          cbind.data.frame(r$res.ori, resF.minus.TIL = r$res.ori[, "resF"] - r$tme[, "T.CD8"])
      }
    }
    return(r)
  }
  
  
  cmb.res.scores <- function(r,
                             res.sig = NULL,
                             bulk.flag = T) {
    res <- r$res
    if (!is.null(res.sig)) {
      res <- get.OE(r = r,
                    sig = res.sig,
                    bulk.flag = bulk.flag)
    }
    
    if (!all(is.element(c("trt", "exc", "exc.seed"), colnames(res)))) {
      print("No merge.")
      return(res)
    }
    res <-
      cbind.data.frame(
        excF = rowMeans(res[, c("exc", "exc.seed")]),
        excF.up = rowMeans(res[, c("exc.up", "exc.seed.up")]),
        excF.down = rowMeans(res[, c("exc.down", "exc.seed.down")]),
        res = rowMeans(res[, c("trt", "exc", "exc.seed")]),
        res.up = rowMeans(res[, c("trt.up", "exc.up", "exc.seed.up")]),
        res.down = rowMeans(res[, c("trt.down", "exc.down", "exc.seed.down")]),
        res
      )
    
    if (is.element("fnc", colnames(res))) {
      res <- cbind.data.frame(
        resF = res[, "res"] + res[, "fnc"],
        resF.up = res[, "res.up"] + res[, "fnc.up"],
        resF.down = res[, "res.down"] + res[, "fnc.down"],
        res
      )
    }
    
    remove.sig <- add.up.down.suffix(c("exc", "exc.seed", "fnc"))
    res <- res[, setdiff(colnames(res), remove.sig)]
    colnames(res) <- gsub("excF", "exc", colnames(res))
    return(res)
  }
  
  
  
  get.OE.bulk <- function(r,
                          gene.sign = NULL,
                          num.rounds = 1000,
                          full.flag = F) {
    set.seed(1234)
    r$genes.mean <- rowMeans(r$tpm)
    r$zscores <- sweep(r$tpm, 1, r$genes.mean, FUN = '-')
    r$genes.dist <- r$genes.mean
    r$genes.dist.q <- discretize(r$genes.dist, n.cat = 50)
    r$sig.scores <-
      matrix(
        data = 0,
        nrow = ncol(r$tpm),
        ncol = length(gene.sign)
      )
    sig.names <- names(gene.sign)
    colnames(r$sig.scores) <- sig.names
    r$sig.scores.raw <- r$sig.scores
    
    rand.flag <-
      is.null(r$rand.scores) |
      !all(is.element(names(gene.sign), colnames(r$rand.scores)))
    
    if (rand.flag) {
      print("Computing also random scores.")
      r$rand.scores <- r$sig.scores
    }
    
    for (i in sig.names) {
      b.sign <- is.element(r$genes, gene.sign[[i]])
      if (sum(b.sign) < 2) {
        next()
      }
      
      if (rand.flag) {
        rand.scores  <-
          get.semi.random.OE(r, r$genes.dist.q, b.sign, num.rounds = num.rounds)
      } else{
        rand.scores <- r$rand.scores[, i]
      }
      
      raw.scores <- colMeans(r$zscores[b.sign, ])
      final.scores <- raw.scores - rand.scores
      r$sig.scores[, i] <- final.scores
      r$sig.scores.raw[, i] <- raw.scores
      r$rand.scores[, i] <- rand.scores
    }
    
    if (full.flag) {
      return(r)
    }
    sig.scores <- r$sig.scores
    return(sig.scores)
    
  }
  
  get.semi.random.OE <- function(r,
                                 genes.dist.q,
                                 b.sign,
                                 num.rounds = 1000,
                                 full.flag = F) {
    # Previous name: get.random.sig.scores
    sign.q <- as.matrix(table(genes.dist.q[b.sign]))
    q <- rownames(sign.q)
    idx.all <- c()
    B <-
      matrix(data = F,
             nrow = length(genes.dist.q),
             ncol = num.rounds)
    Q <-
      matrix(data = 0,
             nrow = length(genes.dist.q),
             ncol = num.rounds)
    
    for (i in 1:nrow(sign.q)) {
      num.genes <- sign.q[i]
      if (num.genes > 0) {
        idx <- which(is.element(genes.dist.q, q[i]))
        for (j in 1:num.rounds) {
          idxj <- sample(idx, num.genes)
          Q[i, j] <- sum(B[idxj, j] == T)
          B[idxj, j] <- T
        }
      }
    }
    
    rand.scores <- apply(B, 2, function(x)
      colMeans(r$zscores[x, ]))
    if (full.flag) {
      return(rand.scores)
    }
    rand.scores <- rowMeans(rand.scores)
    return(rand.scores)
  }
  
  
  get.OE <- function(r,
                     sig,
                     bulk.flag = F,
                     num.rounds = 1000) {
    if (bulk.flag) {
      print("Compute bulk scores")
      scores <- get.OE.bulk(r, sig, num.rounds = num.rounds)
    } else{
      print("Compute single-cell scores")
      scores <- get.OE.sc(r, sig, num.rounds = num.rounds)
    }
    
    names(sig) <- gsub(" ", ".", names(sig))
    two.sided <-
      unique(gsub(".up", "", gsub(".down", "", names(sig))))
    b <- is.element(paste0(two.sided, ".up"), names(sig)) &
      is.element(paste0(two.sided, ".down"), names(sig))
    if (any(b)) {
      two.sided <- two.sided[b]
      scores2 <-
        as.matrix(scores[, paste0(two.sided, ".up")] - scores[, paste0(two.sided, ".down")])
      colnames(scores2) <- two.sided
      scores <- cbind(scores2, scores)
    }
    
    if (!is.null(r$cells)) {
      rownames(scores) <- r$cells
    } else{
      if (!is.null(r$samples)) {
        rownames(scores) <- r$samples
      }
    }
    return(scores)
  }
  
  discretize <- function(v, n.cat, add.prefix = "") {
    q1 <- quantile(v, seq(
      from = (1 / n.cat),
      to = 1,
      by = (1 / n.cat)
    ))
    names(v) <- paste0(add.prefix, names(v))
    return(v)
  }
  
  add.up.down.suffix <- function(v) {
    v <- lapply(v, function(x)
      x <- c(x, paste0(x, ".up"), paste0(x, ".down")))
    v <- unlist(v, use.names = F)
    return(v)
  }
  
  
  # runing the OE calculation
  
  # load signatures
  load(
    paste0(
      sig_dir,
      "ImmuneResistance-master/Results/Signatures/resistance.program.RData"
    )
  )
  load(
    paste0(
      sig_dir,
      "ImmuneResistance-master/Results/Signatures/cell.type.sig.full.RData"
    )
  )
  
  othersigs = data.table::fread(
    paste0(sig_dir, 'ImmuneResistance-master/other47sig_for_OE.txt'),
    sep = '\t',
    header = T
  )
  othersigs = as.list(othersigs)
  othersigs = sapply(othersigs, function(x)
    x[!(x %in% "")])
  
  res.sig2 = c(res.sig, othersigs)
  
  # testing if the TPM matrix came from UCEC_Xena website
  # Note: The value in the file below is: log2(norm_value+1)
  # tcgaRNApath = "/picb/bigdata/project/FengFYM/UCEC_Xena/EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena"
  # tcgaRNAmatrix = as.data.frame(data.table::fread(tcgaRNApath))
  # colnames(tcgaRNAmatrix) = substr(colnames(tcgaRNAmatrix), 1, 12)
  #
  # sample_types = read.table('data/TCGA_cancertypes_subtypes.txt',
  #                           header = T,
  #                           sep = '\t')
  # cancertypes = unique(sample_types$TCGA.Study)
  
  # for (cancer in cancertypes) {
  # tumor_samples = sample_types[sample_types$TCGA.Study == cancer_type,]
  # tcgaRNAmatrix_cancer = tcgaRNAmatrix[colnames(tcgaRNAmatrix) %in% c('sample', tumor_samples$TCGA.Participant.Barcode)]
  #
  #
  #
  # tcgaRNAmatrix_cancer =  aggregate(tcgaRNAmatrix_cancer[-which(colnames(tcgaRNAmatrix_cancer) ==
  #                                                                 "sample")],
  #                                   by = list(tcgaRNAmatrix_cancer$sample),
  #                                   median)
  
  df = list()
  df$tpm = gene_expression # tcgaRNAmatrix_cancer
  rownames(df$tpm) = rownames(df$tpm)
  # df$tpm = df$tpm[-1]
  
  if (length(which(is.na(rowMeans(df$tpm)))) > 0) {
    df$tpm = df$tpm[-which(is.na(rowMeans(df$tpm))), ]
  }
  
  df$samples = colnames(df$tpm)
  df$genes = rownames(df$tpm)
  
  df <- compute.samples.res.scores(
    r = df,
    res.sig = res.sig2,
    cell.sig = cell.sig,
    cc.sig = NULL,
    residu.flag = F,
    num.rounds = 1000
  )
  
  
  if (!file.exists(save_path)) {
    dir.create(save_path)
  }
  save_path = paste0(save_path,"sig_OE_results/")
  if(!file.exists(save_path)){dir.create(save_path)}
  
  save(df, file = paste0(save_path, 'sig_OE_', cancer,"_",sample_type, ".Rdata"))
  # return(df)
  # }
}

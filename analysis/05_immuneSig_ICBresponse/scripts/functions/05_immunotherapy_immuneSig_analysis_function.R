icb_compr <- function(result,
                      response,
                      resultdir,
                      filename) {
                        
  if(nrow(response[is.na(response$label), ])>0)response <- response[!is.na(response$label), ]

  # scale signature value (z-score)
  result_scaled <- result[-1] # as.data.frame(apply(result[-1],2, scale))
  result_scaled <- cbind(result[1], result_scaled)
  names(response) <- c('patient','label')
  res_df <- merge(response, result_scaled, by = 'patient')
  res_df2 = as.list(res_df)
 
  # number of patients ----------------------------------------------------------
  R_num <- nrow(res_df[which(res_df$label=='R'),])
  NR_num <- nrow(res_df[which(res_df$label=='NR'),])



 # signature mean value --------------------------------------------------------
  R_mean <- lapply(res_df2[names(res_df2)[-c(1,2)]], function(x){
   if(!is.na(x)){
    mean(x[which(res_df2$label=='R')])
  }else {
     NA
  } 
 })
  R_mean <- do.call(rbind, R_mean)

  NR_mean <- lapply(res_df2[names(res_df2)[-c(1,2)]], function(x){
   if(!is.na(x)){
    mean(x[which(res_df2$label=='NR')])
  } else {
     NA
  }
 })
  NR_mean <- do.call(rbind, NR_mean)

  simple_stat <- cbind(data.frame(R_num = rep(R_num, length(R_mean))), R_mean, 
                       data.frame(NR_num = rep(NR_num, length(NR_mean))),  NR_mean)

  # p-values  ------------------------------------------------------------------
  p_val <- lapply(res_df2[names(res_df2)[-c(1,2)]], function(x){
   if(!is.na(x)){
    wilcox.test(x[which(res_df2$label=='R')], x[which(res_df2$label =='NR')])$p.value
  } else {
     NA
  }
 })
  
  p_val_greater <- lapply(res_df2[names(res_df2)[-c(1,2)]], function(x){
    if(!is.na(x)){

    wilcox.test(x[which(res_df2$label=='R')], x[which(res_df2$label =='NR')], alternative = 'greater')$p.value
  }else {
     NA
  }
  })
  
  p_val_less <- lapply(res_df2[names(res_df2)[-c(1,2)]], function(x){
    if(!is.na(x)){
    wilcox.test(x[which(res_df2$label=='R')], x[which(res_df2$label =='NR')], alternative = 'less')$p.value
    }else {
     NA
  }
  })
   
  p_val <- do.call(rbind, p_val)
  p_val_greater <- do.call(rbind, p_val_greater)
  p_val_less <- do.call(rbind, p_val_less)
  p_val <- cbind(p_val, p_val_greater, p_val_less)
  colnames(p_val) = c('wilcox_P_value', 'R_greater_NR_P_value', 'R_less_NR_P_value')
  

  # signature predictive power -------------------------------------------------

  # Over simplified version
  # auc <- lapply(res_df2[names(res_df)[-c(1,2)]], function(x){
  #   tmp = data.frame(label = res_df2$label, immunesig = x)
  #   pred <- prediction( tmp$immunesig, tmp$label)
  #   aucPerf <- performance( pred, 'auc' )
  #   AUCValue<-aucPerf@y.values[[1]]
  #   return(AUCValue)
  # })
  
  # auc <- do.call(rbind, auc)
  # colnames(auc) = 'auc'
  # result_test = data.frame(cbind(p_val, auc))
  # result_test$immunesig = rownames(result_test)
  
  # update auc computing method
  # inspired by AnyuanGuo

  # res_df = res_df[order(res_df$label), ]
  # res_df$label <- factor(res_df$label,levels = c("NR","R"), order=TRUE)
  # res_df$tag <- 1
  # res_df[res_df$label == "NR",]$tag = 0
  # res_df$tag <- factor(res_df$tag,levels = c(0,1), order=TRUE)

  result_scaled <- as.data.frame(apply(result[-1],2, scale))
  result_scaled <- cbind(result[1], result_scaled)
  names(response) <- c('patient','label')
  res_df <- merge(response, result_scaled, by = 'patient')
  res_df = res_df[order(res_df$label), ]
  res_df$label <- factor(res_df$label,levels = c("NR","R"), order=TRUE)
  res_df$tag <- 1
  res_df[res_df$label == "NR",]$tag = 0
  res_df$tag <- factor(res_df$tag,levels = c(0,1), order=TRUE)

  # 1) training/testing set Stratified sampling
  library(sampling)
  library(pROC)
  set.seed(1234)

 auc_500sampling <- lapply(as.list(seq(500)), function(i){

  n = round(table(res_df$label) * 2/3)
  sub_train = strata(res_df, strataname = ("tag"), size= n, method = "srswor")
  data_train = res_df[sub_train$ID_unit, ]
  data_test = res_df[-sub_train$ID_unit, ]
  data_train_list = as.list(data_train[-c(1,2, grep("tag", names(data_train)))])
  data_test_list = as.list(data_test[-c(1,2, grep("tag", names(data_train)))])

  auc <- lapply(as.list(seq(length(data_train_list))), function(x){
    if(!is.na(c(data_train_list[[x]], data_test_list[[x]]))){
    
    sig_train = data.frame(label = data_train$tag, immunesig = data_train_list[[x]])
    sig_test = data.frame(label = data_test$tag, immunesig = data_test_list[[x]])

    # 2) glm() modeling
    mod = glm(label~., data = sig_train, family = "binomial")
    # 3) predicting
    prob <- predict(object = mod, newdata = sig_test, type = "response")
    # 3) access AUC of testing set

    pred <- prediction(prob, sig_test$label)
    aucPerf <- performance( pred, 'auc' )
    AUCValue<-aucPerf@y.values[[1]]
    return(AUCValue)
    } else {
    return( NA)
  }  
    # roc_curve <- roc(sig_test$label, prob)
  })
  names(auc) = names(data_train_list)
  auc <- do.call(rbind, auc)
  return(auc)
  })


  auc_500sampling <- do.call(cbind, auc_500sampling)
  auc_mean <- data.frame(auc = rowMeans(auc_500sampling))
  # colnames(auc_mean) = 'auc'

  # prepare estimator for meta analysis ----------------------------------------
  # REFERENCE: Meta-analysis of tumor- and T cell-intrinsic mechanisms of sensitization to checkpoint inhibition
  # https://github.com/kevlitchfield1/CPI1000_paper/blob/main/code/Figure2/Fig2_script.R
  
  ## estimation --------------------
  result_scaled <- as.data.frame(apply(result[-1],2, scale))
  result_scaled <- cbind(result[1], result_scaled)
  names(response) <- c('patient','label')
  res_df <- merge(response, result_scaled, by = 'patient')
  res_df = res_df[order(res_df$label), ]
  res_df$label <- factor(res_df$label,levels = c("NR","R"), order=TRUE)
  res_df$tag <- 1
  res_df[res_df$label == "NR",]$tag = 0
  res_df$tag <- factor(res_df$tag,levels = c(0,1), order=TRUE)
  res_df2 = as.list(res_df)

  meta_estimator <- lapply(res_df2[names(res_df2)[-c(1,2)]], function(x){

    if(!is.na(x)){
    # tmp = data.frame(label = res_df$label, immunesig = x)
    tmp = data.frame(patient = res_df$patient, label = res_df$label, immunesig = x)
    # tmp$immunesig = scale(tmp$immunesig)

    tmp$tag <- 1
    tmp[tmp$label=="NR",]$tag <- 0
    # tmp$tag <- factor(tmp$tag,levels = c(0,1), order=TRUE)

    mod = glm(tag~immunesig, data = tmp, family = "binomial")
    b <- summary(mod)
    if(nrow(b$coefficients)>1){

      TE <-b$coefficients[2,1]
      seTE <- b$coefficients[2,2]
      TE_Zval <- b$coefficients[2,3]
      TE_Pval <- b$coefficients[2,4]
    }else{      
      TE <- NA 
      seTE <- NA
      TE_Zval <- NA
      TE_Pval <- NA
    }
    res <- data.frame(TE = TE, seTE = seTE, TE_Zval = TE_Zval, TE_Pval = TE_Pval)
    }else{
    res <- data.frame(TE = NA, seTE = NA, TE_Zval = NA, TE_Pval = NA)
    }

    return(res)

  })
  meta_estimator <- do.call(rbind, meta_estimator)

  # result_test = data.frame(cbind(simple_stat, p_val, auc_mean, meta_estimator))
  simple_stat$immunesig = row.names(simple_stat)
  p_val = as.data.frame(p_val); p_val$immunesig = row.names(p_val)
  auc_mean = as.data.frame(auc_mean); auc_mean$immunesig = row.names(auc_mean)
  meta_estimator$immunesig = row.names(meta_estimator)
  # result_test = data.frame(cbind(simple_stat, p_val, auc_mean, meta_estimator))
  result_test = inner_join(simple_stat, inner_join(p_val, inner_join(auc_mean, meta_estimator)))
 
  result_test = result_test[c("immunesig", "R_num","R_mean","NR_num", "NR_mean", 
                              "wilcox_P_value", "R_greater_NR_P_value","R_less_NR_P_value","auc",
                              "TE","seTE","TE_Zval","TE_Pval")]

  # result_test$immunesig = rownames(result_test)
  

  write.csv(result_test, paste0(resultdir, filename), quote = F)
  
}


expr4cibersort <- function(expr,
                           dataset,
                           savepath){
  genenames <- data.frame(GeneSymbol = rownames(expr))
  cbind(genenames, expr) %>% write.table(
    paste0(savepath, "/", dataset, "_TPM.txt"),
    sep = "\t",
    quote = F,
    row.names = F
  )
}

expr4immuncelai <- function(expr,
                           dataset,
                           savepath){
  expr <- log2(expr + 1)
  genenames <- data.frame(GeneSymbol = rownames(expr))
  cbind(genenames, expr) %>% write.table(
    paste0(savepath, "/", dataset, "_TPM.txt"),
    sep = "\t",
    quote = F,
    row.names = F
  )

}

expr4tide <- function(expr,
                      dataset,
                      savepath){
  expr_normed <- log2(expr+1) - apply(log2(expr+1),1,mean)
  # genenames <- data.frame(GeneSymbol = rownames(expr))
  # cbind(genenames, expr_normed) 
  expr_normed %>% write.table(
    paste0(savepath, "/", dataset, "_TPM_normed.txt"),
    sep = "\t",
    quote = F,
    row.names = T
  )
}


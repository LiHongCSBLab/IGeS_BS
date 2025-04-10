glm_model <- function(result,
                      response,
                      resultdir,
                      filename) {

  rownames(result) <- result$patientID
  result <- as.data.frame(scale(result[-1]))
  result <- cbind(data.frame(patient = rownames(result)), result)
  
  names(response) <- c('patient','label')
  res_df <- merge(response, result, by = 'patient')
  res_df2 = as.list(res_df)
  
  
  
  # update auc computing method
  # inspired by AnyuanGuo
  res_df = res_df[order(res_df$label), ]
  res_df$label <- factor(res_df$label,levels = c("NR","R"), order=TRUE)
  res_df$tag <- 1
  res_df[res_df$label == "NR",]$tag = 0
  res_df$tag <- factor(res_df$tag,levels = c(0,1), order=TRUE)

  # 1) training/testing set Stratified sampling
  library(glmnet)
  library(sampling)
  library(pROC)
  set.seed(1234)

 auc_500sampling <- lapply(as.list(seq(500)), function(i){

  n = round(table(res_df$label) * 2/3)
  sub_train = strata(res_df, strataname = ("label"), size= n, method = "srswor")
  data_train = res_df[sub_train$ID_unit, ]
  data_test = res_df[-sub_train$ID_unit, ]

  data_train_mat = data_train[-c(1,2, grep("tag", names(data_train)))]
  data_test_mat = data_test[-c(1,2, grep("tag", names(data_train)))]
    
    sig_train = cbind(data_train['tag'], data_train_mat)
    sig_test = cbind(data_test['tag'], data_test_mat)

    # 2) glm() modeling
    mod = glm(tag~., data = sig_train, family = "binomial")


    # 3) predicting
    prob <- predict(object = mod, newdata = sig_test, type = "response")
    # 3) access AUC of testing set
    pred <- prediction(prob, sig_test$tag)
    aucPerf <- performance( pred, 'auc' )
    AUCValue<-aucPerf@y.values[[1]]
    return(AUCValue)

 })


  auc_500sampling <- do.call(cbind, auc_500sampling)
  auc_mean <- data.frame(auc = rowMeans(auc_500sampling))
  
}



glmnet_model <- function(result,
                      response,
                      resultdir,
                      filename) {

  # result = ITSp_gsva[-2]
  rownames(result) <- result$patient
  # result <- as.data.frame(scale(result[-1]))
  result <- cbind(data.frame(patient = rownames(result)), result)
  
  names(response) <- c('patient','label')
  res_df <- merge(response, result, by = 'patient')
  
  

  
  # update auc computing method
  # inspired by AnyuanGuo
  res_df = res_df[order(res_df$label), ]
  res_df$label <- factor(res_df$label,levels = c("NR","R"), order=TRUE)
  res_df$tag <- 1
  res_df[res_df$label == "NR",]$tag = 0
  res_df$tag <- factor(res_df$tag,levels = c(0,1), order=TRUE)

  # 1) training/testing set Stratified sampling
  library(glmnet)
  library(sampling)
  library(pROC)
  set.seed(1234)

 auc_500sampling <- lapply(as.list(seq(500)), function(i){

  n = round(table(res_df$label) * 2/3)
  sub_train = strata(res_df, strataname = ("label"), size= n, method = "srswor")
  data_train = res_df[sub_train$ID_unit, ]
  data_test = res_df[-sub_train$ID_unit, ]

  data_train_mat = data_train[-c(1,2, grep("tag", names(data_train)))]
  data_test_mat = data_test[-c(1,2, grep("tag", names(data_train)))]
    
    sig_train = cbind(data_train['tag'], data_train_mat)
    sig_test = cbind(data_test['tag'], data_test_mat)

    
    cv.fit=cv.glmnet(as.matrix(sig_train[-c(1)]),sig_train$tag,family='binomial', alpha = 0.5,type.measure="auc")

  for (i in 0:10) {
    assign(paste("fit", i, sep=""), cv.glmnet(as.matrix(sig_train[-c(1)]),sig_train$tag, type.measure="auc", 
                                              alpha=i/10,family='binomial'))
  }
  summary(fit5)
  # # Plot solution paths:
  # par(mfrow=c(3,1))
  # # For plotting options, type '?plot.glmnet' in R console
  # plot(fit10, main="LASSO")

  # plot(fit0, main="Ridge")

  # plot(fit5, main="Elastic Net")


    fit<-glmnet(as.matrix(sig_train[-c(1)]),sig_train$tag,family="binomial", alpha = 0.5, lambda = fit5$lambda)
    coefficients<-coef(fit,s=cv.fit$lambda.min)

    prob <- predict(fit5, newx = as.matrix(sig_test[-1]), alpha = 0.5, s = fit5$lambda.min)  # make predictions

    pred <- prediction(prob, sig_test$tag)
    aucPerf <- performance( pred, 'auc' )
    AUCValue<-aucPerf@y.values[[1]]
    return(AUCValue)

 })


  auc_500sampling <- do.call(cbind, auc_500sampling)
  auc_mean <- data.frame(auc = rowMeans(auc_500sampling))
  
}



svm_model <- function(result,
                      response,
                      resultdir,
                      filename) {

  rownames(result) <- result$patient
  result <- as.data.frame(scale(result[-1]))
  result <- cbind(data.frame(patient = rownames(result)), result)
  
  names(response) <- c('patient','label')
  res_df <- merge(response, result, by = 'patient')
  res_df2 = as.list(res_df)
  
  

  
  # update auc computing method
  # inspired by AnyuanGuo
  res_df = res_df[order(res_df$label), ]
  res_df$label <- factor(res_df$label,levels = c("NR","R"), order=TRUE)
  res_df$tag <- 1
  res_df[res_df$label == "NR",]$tag = 0
  res_df$tag <- factor(res_df$tag,levels = c(0,1), order=TRUE)

  # 1) training/testing set Stratified sampling
  library(e1071)
  library(sampling)
  library(pROC)
  set.seed(1234)

  auc_500sampling <- lapply(as.list(seq(500)), function(i){

  n = round(table(res_df$label) * 2/3)
  sub_train = strata(res_df, strataname = ("label"), size= n, method = "srswor")
  data_train = res_df[sub_train$ID_unit, ]
  data_test = res_df[-sub_train$ID_unit, ]

  data_train_mat = data_train[-c(1,2, grep("tag", names(data_train)))]
  data_test_mat = data_test[-c(1,2, grep("tag", names(data_train)))]
    
    sig_train = cbind(data_train['tag'], data_train_mat)
    sig_test = cbind(data_test['tag'], data_test_mat)


    
    # fit<-svm(tag ~ ., data = sig_train, kernel = "linear", cost = 10, scale = FALSE)
    fit<-svm(tag~.,data = sig_train,kernel="radial",cost=1,gamma=1/ncol(sig_train))
    summary(fit)

    prob <- predict(fit, as.matrix(sig_test[-1]))  # make predictions
    svm.table<- table(prob,sig_test$tag)
    
    library(caret)
    measurement <- confusionMatrix(svm.table)

    # pred <- prediction(prob, sig_test$tag)
    # aucPerf <- performance( pred, 'auc' )
    # AUCValue<-aucPerf@y.values[[1]]
    return(measurement$overall[1])

 })


  auc_500sampling <- do.call(cbind, auc_500sampling)
  auc_mean <- data.frame(auc = rowMeans(auc_500sampling))
  
}




rf_model <- function(result,
                      response,
                      resultdir,
                      filename) {

  rownames(result) <- result$patient
  result <- as.data.frame(scale(result[-1]))
  result <- cbind(data.frame(patient = rownames(result)), result)
  
  names(response) <- c('patient','label')
  res_df <- merge(response, result, by = 'patient')
  res_df2 = as.list(res_df)
  
  

  
  # update auc computing method
  # inspired by AnyuanGuo
  res_df = res_df[order(res_df$label), ]
  res_df$label <- factor(res_df$label,levels = c("NR","R"), order=TRUE)
  res_df$tag <- 1
  res_df[res_df$label == "NR",]$tag = 0
  res_df$tag <- factor(res_df$tag,levels = c(0,1), order=TRUE)

  # 1) training/testing set Stratified sampling
  library(randomForest)
  library(sampling)
  library(pROC)
  set.seed(1234)

  auc_500sampling <- lapply(as.list(seq(500)), function(i){

  n = round(table(res_df$label) * 2/3)
  sub_train = strata(res_df, strataname = ("label"), size= n, method = "srswor")
  data_train = res_df[sub_train$ID_unit, ]
  data_test = res_df[-sub_train$ID_unit, ]

  data_train_mat = data_train[-c(1,2, grep("tag", names(data_train)))]
  data_test_mat = data_test[-c(1,2, grep("tag", names(data_train)))]
    
    sig_train = cbind(data_train['tag'], data_train_mat)
    sig_test = cbind(data_test['tag'], data_test_mat)


    
    fit <- randomForest(tag ~ ., data = sig_train, importance = TRUE)
    # importance_fit <- data.frame(importance(fit))
    # head(importance_fit)


    prob <- predict(fit, as.matrix(sig_test[-1]))  # make predictions
    svm.table<- table(prob,sig_test$tag)
    
    library(caret)
    measurement <- confusionMatrix(svm.table)

    # pred <- prediction(prob, sig_test$tag)
    # aucPerf <- performance( pred, 'auc' )
    # AUCValue<-aucPerf@y.values[[1]]
    return(measurement$overall[1])

 })


  auc_500sampling <- do.call(cbind, auc_500sampling)
  auc_mean <- data.frame(auc = rowMeans(auc_500sampling))
  
}



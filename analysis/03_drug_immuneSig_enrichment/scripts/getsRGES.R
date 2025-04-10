#find best alpha and beta
find_alpha_beta <- function(){
  alphas <- seq(-1, 1, 0.1)
  betas <- seq(-1, 1, 0.1)
  all_values <- data.frame()
  for (alpha in alphas){
    for (beta in betas){
      lincs_drug_prediction_subset <- subset(lincs_drug_prediction, cell_id %in% c(cell_line_selected)) #HT29 MCF7
      lincs_drug_prediction_subset$RGES <- sapply(1:nrow(lincs_drug_prediction_subset), function(id){
        getsRGES(lincs_drug_prediction_subset[id,"RGES"], lincs_drug_prediction_subset[id, "pert_dose"], lincs_drug_prediction_subset[id, "pert_time"], alpha, beta)
      })
      lincs_drug_prediction_subset <- aggregate(RGES ~ pert_iname, lincs_drug_prediction_subset, mean)
      
      activity_RGES <- merge(lincs_drug_prediction_subset, lincs_drug_activity_subset, by="pert_iname")
      
      activity_RGES_summarized <- activity_RGES #aggregate(cbind(RGES, standard_value) ~ pert_iname, activity_RGES,  min)
      
      cor <- cor(activity_RGES_summarized$RGES, log(activity_RGES_summarized$standard_value, 10), method="spearman")
      all_values <- rbind(all_values, data.frame(cor, alpha, beta))
    }
  }
  return(all_values)
}

getsRGES1 <- function(RGES, cor, pert_dose, pert_time){
  sRGES <- RGES
  if (pert_time == 24){
    sRGES <- RGES + predict(lm_treatment_dose_24, data.frame(treatment_dose=round(log(as.numeric(pert_dose), 10), 1)))
  }
  if (pert_time == 6){
    sRGES <- RGES + predict(lm_treatment_dose_6, data.frame(treatment_dose=round(log(as.numeric(pert_dose), 10), 1)))
  }
  return (sRGES * cor/max_cor )
}

getsRGES2 <- function(RGES, cor, pert_dose, pert_time){
  sRGES <- RGES 
  
  #older version
  if (pert_time < 24){
    sRGES <- sRGES - 0.1
  }
  
  if (pert_dose < 10){
    sRGES <- sRGES - 0.2
  }
  return(sRGES * cor)
}

getsRGES3 <- function(RGES, cor, pert_dose, pert_time, diff, max_cor){
  
  sRGES <- RGES
  pert_time <- ifelse(pert_time < 24, "short", "long")
  pert_dose <- ifelse(pert_dose < 10, "low", "high")
  if (pert_time == "short" & pert_dose == "low"){
    sRGES <- sRGES + diff[4]
  }
  if (pert_dose ==  "low" & pert_time == "long"){
    sRGES <- sRGES + diff[2]
  }
  if (pert_dose ==  "high" & pert_time == "short"){
    sRGES <- sRGES + diff[1]
  }
  return(sRGES ) #* cor/max_cor
}

##################
####
getsRGES <- function(RGES, cor, pert_dose, pert_time, diff, max_cor){
  
  sRGES <- RGES
  pert_time <- ifelse(pert_time < 24, "short", "long")
  pert_dose <- ifelse(pert_dose < 10, "low", "high")
  if (pert_time == "short" & pert_dose == "low"){
    sRGES <- sRGES + diff[4]
  }
  if (pert_dose ==  "low" & pert_time == "long"){
    sRGES <- sRGES + diff[2]
  }
  if (pert_dose ==  "high" & pert_time == "short"){
    sRGES <- sRGES + diff[1]
  }
  
  return(sRGES * cor/max_cor) #
  # return(sRGES )
}


#summarize NES
merge_NES <- function(resfile,
                      lincs_cell_line_weight){
  # should use pert_dose > 0.01
  require(plyr)
  
  res_subset <- subset(resfile,  treatment_dose > 0 & treatment_time %in% c(6, 24))
  #pairs that share the same drug and cell id
  res_pairs <- merge(res_subset, res_subset, by=c("drug_index", "cell_id")) 
  #x is the reference
  res_pairs <- subset(res_pairs, files.x != files.y & treatment_time.x == 24 & treatment_dose.x == 10)
  
  #difference of RGES to the reference 
  res_pairs$cmap_diff <- res_pairs$NES.x - res_pairs$NES.y
  res_pairs$dose <- round(log(as.numeric(res_pairs$treatment_dose.y), 2), 1)
  
  #estimate difference
  res_pairs$dose_bin <- ifelse(res_pairs$treatment_dose.y < 10, "low", "high")
  diff <- tapply(res_pairs$cmap_diff, paste(res_pairs$dose_bin, res_pairs$treatment_time.y), mean)
  
  #ignore weighting cell lines
  pred <- merge(resfile, lincs_cell_line_weight, by.x="cell_id")
  pred$RGES <- sapply(1:nrow(pred), function(id) {
    getsRGES(
      RGES = pred$NES[id],
      cor = pred$mean_cor[id],
      pert_dose = pred$treatment_dose[id],
      pert_time = pred$treatment_time[id],
      diff = diff,
      max_cor = max(pred$mean_cor)
    )
  })
  
  # res_pairs_subset <- subset(res_pairs, treatment_time.y == 24 )
  # treatment_dose_cmap_diff_24 <- tapply(res_pairs_subset$cmap_diff, res_pairs_subset$treatment_dose.y, mean)
  # treatment_dose_cmap_diff_24 <- data.frame(treatment_dose = as.numeric(names(treatment_dose_cmap_diff_24)), cmap_diff= treatment_dose_cmap_diff_24)
  # lm_treatment_dose_24 <- lm(cmap_diff ~ treatment_dose, data = treatment_dose_cmap_diff_24)
  # coef(lm_treatment_dose_24)
  
  # res_pairs_subset <- subset(res_pairs, treatment_time.y == 6)
  # treatment_dose_cmap_diff_6 <- tapply(res_pairs_subset$cmap_diff, res_pairs_subset$treatment_dose.y, mean)
  # treatment_dose_cmap_diff_6 <- data.frame(treatment_dose = as.numeric(names(treatment_dose_cmap_diff_6)), cmap_diff= treatment_dose_cmap_diff_6)
  # lm_treatment_dose_6 <- lm(cmap_diff ~ treatment_dose, data = treatment_dose_cmap_diff_6)
  
  # pred$RGESfixtime <- sapply(1:nrow(pred), function(id) {
  #   getsRGES1(
  #     RGES = pred$NES[id],
  #     cor = pred$mean_cor[id],
  #     pert_dose = pred$treatment_dose[id],
  #     pert_time = pred$treatment_time[id]
  #   )
  # })
  # cmpd_freq <- table(pred$drug_index)
  # pred <- subset(pred, drug_index %in% names(cmpd_freq[cmpd_freq > 0]))
  
  pred_merged <- plyr::ddply(pred,  .(drug_index),  summarise,
                             mean = mean(RGES),
                             n = length(RGES),
                             median = median(RGES),
                             sd = sd(RGES))
  pred_merged$mNES <- pred_merged$mean
  pred_merged <- pred_merged[order(pred_merged$mNES), ]
  
  
  if(length(unique(resfile$treatment_time)<=1)){
    if(all(is.na(resfile$treatment_dose))|all(is.na(resfile$NES))){
      pred_merged$coef =5
    }else{
      lm_test <- lm(NES ~ as.numeric(treatment_dose), data = resfile)
      pred_merged$coef = coef(lm_test)[[2]]
    }
  }else if(length(unique(resfile$treatment_dose)<=1)){
    if(all(is.na(resfile$treatment_time))|all(is.na(resfile$NES))){
      pred_merged$coef =5
    }else{
      lm_test <- lm(NES ~ as.numeric(treatment_time), data = resfile)
      pred_merged$coef = coef(lm_test)[[2]]
    }
  }else if(length(unique(resfile$treatment_dose)>1 & length(unique(resfile$treatment_time)>1))){
    
    if(all(is.na(resfile$treatment_dose)) | all(is.na(resfile$treatment_time))|all(is.na(resfile$NES))){
      pred_merged$coef =5
    }else{
      lm_test <- lm(NES ~ as.numeric(treatment_time)*as.numeric(treatment_dose), data = resfile)
      pred_merged$coef = coef(lm_test)[[2]]
    }
  }
  
  return(pred_merged)
}



getsRGESmodified <- function(RGES, cor, pert_dose, pert_time, diff, max_cor){
  
  sRGES <- RGES
  if(length(RGES) == 1){
    if(pert_time == 24 & pert_dose == 10){
       sRGES <- sRGES
    }else if(pert_time < 24 & pert_dose < 10){
      #  sRGES <- sRGES*(pert_time/24)*(pert_dose/10)
       sRGES <- sRGES # +0.1
    }else if(pert_time > 24 & pert_dose > 10 ){
      sRGES <- sRGES # -0.1
    }else{
      sRGES <- sRGES
    }
   }else{

  pert_time <- ifelse(pert_time < 24, "short", "long")
  pert_dose <- ifelse(pert_dose < 10, "low", "high")
  
  if(length(diff)==1){
    sRGES <- sRGES + diff[1]
  }else{
    
    if (pert_time == "short" & pert_dose == "low"){
      # sRGES <- sRGES + diff[4]
      sRGES <- sRGES + diff['low 6']
    }
    if (pert_dose ==  "low" & pert_time == "long"){
      sRGES <- sRGES + diff['low 24']
    }
    
    if (pert_dose ==  "high" & pert_time == "short"){
      sRGES <- sRGES + diff['high 6']
    }

    if (pert_dose ==  "high" & pert_time == "long"){
      sRGES <- sRGES + diff['high 24']
    }
    }
  }

  return(sRGES * cor/max_cor) #
  # return(sRGES )
}

#summarize NES
# modified to deal with when there is no correlation score between cell lines 
# and tumor samples to correct the merge score
merge_NESmodified <- function(resfile,
                      cancer,
                      lincs_cell_line_weight) {
  # should use pert_dose > 0.01
  require(plyr)
  resfile$treatment_time = as.numeric(resfile$treatment_time)
  resfile$treatment_dose = as.numeric(resfile$treatment_dose)
  res_subset <- subset(resfile,  as.numeric(treatment_dose) > 0 & as.numeric(treatment_time) %in% c(6, 24))
  #pairs that share the same drug and cell id
  res_pairs <- merge(res_subset, res_subset, by=c("drug_index", "cell_id")) 
  #x is the reference
  # res_pairs <- subset(res_pairs, files.x != files.y & treatment_time.x == 24 & treatment_dose.x == 10)
  
  # #difference of RGES to the reference 
  # res_pairs$cmap_diff <- res_pairs$NES.x - res_pairs$NES.y
  # res_pairs$dose <- round(log(as.numeric(res_pairs$treatment_dose.y), 2), 1)
  
  #estimate difference
  # res_pairs$dose_bin <- ifelse(res_pairs$treatment_dose.y < 10, "low", "high")
  # diff <- tapply(res_pairs$cmap_diff, paste(res_pairs$dose_bin, res_pairs$treatment_time.y), mean)

  if(is.element(24, res_pairs$treatment_time.x) & is.element(10, res_pairs$treatment_dose.x)){
    res_pairs <- subset(res_pairs, files.x != files.y & treatment_time.x == 24 & treatment_dose.x == 10)
    #difference of RGES to the reference 
    res_pairs$cmap_diff <- res_pairs$NES.x - res_pairs$NES.y
    res_pairs$dose <- round(log(as.numeric(res_pairs$treatment_dose.y), 2), 1)
    
    #estimate difference
    res_pairs$dose_bin <- ifelse(as.numeric(res_pairs$treatment_dose.y) < 10, "low", "high")
    diff <- tapply(res_pairs$cmap_diff, paste(res_pairs$dose_bin,as.numeric( res_pairs$treatment_time.y)), mean)
    
  }else{
    res_pairs <- subset(res_pairs, files.x != files.y & treatment_time.x == max(res_pairs$treatment_time.x) & treatment_dose.x == max(res_pairs$treatment_dose.x))
    #difference of RGES to the reference 
    res_pairs$cmap_diff <- res_pairs$NES.x - res_pairs$NES.y
    res_pairs$dose <- round(log(as.numeric(res_pairs$treatment_dose.y), 2), 1)
    
    #estimate difference
    res_pairs$dose_bin <- ifelse(as.numeric(res_pairs$treatment_dose.y) < 10, "low", "high")
    diff <- tapply(res_pairs$cmap_diff, paste(res_pairs$dose_bin, as.numeric(res_pairs$treatment_time.y)), mean)
    
  }

  # modified to deal if cell line doesnt have pre-calculate correlation with Tumor samples
  if(length(intersect(resfile$cell_id, lincs_cell_line_weight$cell_id)) > 0 ){
    # original one
    pred <- merge(resfile, lincs_cell_line_weight, by.x="cell_id")
    pred$RGES <- sapply(1:nrow(pred), function(id) {
      getsRGESmodified(
        RGES = pred$NES[id],
        cor = pred$mean_cor[id],
        pert_dose = pred$treatment_dose[id],
        pert_time = pred$treatment_time[id],
        diff = diff,
        max_cor = max(pred$mean_cor)
      )
    })

   }else{
    
    pred <- resfile
    pred$mean_cor <- mean(lincs_cell_line_weight[lincs_cell_line_weight$cancertype==cancer,]$mean_cor)
    pred$max_cor <- max(lincs_cell_line_weight[lincs_cell_line_weight$cancertype==cancer,]$mean_cor)
    pred$RGES <- sapply(1:nrow(pred), function(id) {
      getsRGESmodified(
        RGES = pred$NES[id],
        cor = pred$mean_cor[id],
        pert_dose = pred$treatment_dose[id],
        pert_time = pred$treatment_time[id],
        diff = diff,
        max_cor = max(pred$mean_cor)
      )
    })
  }
  
  pred_merged <- plyr::ddply(pred,  .(drug_index),  summarise,
                             mean = mean(RGES),
                             n = length(RGES),
                             median = median(RGES),
                             sd = sd(RGES))
  pred_merged$mNES <- pred_merged$mean
  pred_merged <- pred_merged[order(pred_merged$mNES), ]
  
  
  if(length(unique(resfile$treatment_time)) <= 1){
    if(all(is.na(resfile$treatment_dose))|all(is.na(resfile$NES))){
      pred_merged$coef = 5
    }else{
      lm_test <- lm(NES ~ as.numeric(treatment_dose), data = resfile)
      pred_merged$coef = coef(lm_test)[[2]]
    }
  }else if(length(unique(resfile$treatment_dose))<=1){
    if(all(is.na(resfile$treatment_time))|all(is.na(resfile$NES))){
      pred_merged$coef =5
    }else{
      lm_test <- lm(NES ~ as.numeric(treatment_time), data = resfile)
      pred_merged$coef = coef(lm_test)[[2]]
    }
  }else if(length(unique(resfile$treatment_dose))>1 & length(unique(resfile$treatment_time)>1)){
    
    if(all(is.na(resfile$treatment_dose)) | all(is.na(resfile$treatment_time))|all(is.na(resfile$NES))){
      pred_merged$coef =5
    }else{
      lm_test <- lm(NES ~ as.numeric(treatment_time)*as.numeric(treatment_dose), data = resfile)
      pred_merged$coef = coef(lm_test)[[2]]
    }
  }
  
  return(pred_merged)
}

library(meta)
data(Fleiss93)
num_study <- nrow(Fleiss93)
# ------------------------------------------------------------------------------
# Calculation of fixed effect and random effects estimates (risk ratio, odds ratio, risk difference, arcsine difference, or diagnostic odds ratio) for meta-analyses with binary outcome data.

# method: By default, both fixed effect (also called common effect) and random effects models are considered
res1 = metabin(event.e, # Number of events in experimental group or true positives in diagnostic study.
        n.e, # Number of observations in experimental group or number of ill participants in diagnostic study.
        event.c, # Number of events in control group or false positives in diagnostic study.
        n.c, # Number of observations in control group or number of healthy participants in diagnostic study
        data = Fleiss93,
        sm = "OR") # which summary measure ("RR", "OR", "RD", "ASD", or "DOR") is to be used for pooling of studies

# P-value of heterogeneity test.
res1$pval.Q

if(res1$pval.Q < 0.05){
   # P-value of heterogeneity test < 0.05, using random effects meta-analysis  
   res2 <- update(res1, fixed = FALSE)
}else if(res1$pval.Q > 0.05){
   # P-value of heterogeneity test < 0.05, using fixed effects meta-analysis  
   res2 <- update(res1, random = FALSE)
}

res3 <- summary(res2)

# forest plot
forest(res2)

# 发表偏倚检验
funnel(res2)
# Note, the minimum number of studies is three
if(num_study <= 10){
    res_biastest_egger<-metabias(res2, method.bias="egger", k.min=3)
}else if(num_study > 10){
    res_biastest_peters<-metabias(res2, method.bias="peters", k.min=3)
}
# 该meta分析存在发表偏倚的可能性。
# 使用trim and filled或者copas模型等进行校正
tf1 <- trimfill(res2, comb.fixed=TRUE)
summary(tf1)
funnel(tf1)



# ------------------------------------------------------------------------------
# Fixed effect and random effects meta-analysis based on estimates (e.g. log hazard ratios) 
# and their standard errors. The inverse variance method is used for pooling.

# used by (2021) Meta-analysis of tumor- and T cell-intrinsic mechanisms of sensitization to checkpoint inhibition
# main function: metabin()
Fleiss93$TE <- res3$TE
Fleiss93$seTE <-res3$seTE
res_metagen1 <- metagen(TE, seTE, studlab = study,
                       sm = "OR",data = Fleiss93)

if(res_metagen1$pval.Q < 0.05){
   # P-value of heterogeneity test < 0.05, using random effects meta-analysis  
   res_metagen2 <- update(res_metagen1, fixed = FALSE)
}else if(res_metagen1$pval.Q > 0.05){
   # P-value of heterogeneity test < 0.05, using fixed effects meta-analysis  
   res_metagen2 <- update(res_metagen1, random = FALSE)
}
summary(res_metagen2)



# ------------------------------------------------------------------------------
immunesig_index <- unique(resall$immune_sig)
x = immunesig_index[4]
tmp <- resall[resall$immune_sig == x, ]
res_metagen1 <- metagen(TE, seTE, studlab = dataset, sm = "OR", data = tmp)


if(res_metagen1$pval.Q < 0.05){
  # P-value of heterogeneity test < 0.05, using random effects meta-analysis  
  res_metagen2 <- update(res_metagen1, fixed = FALSE)
  summary(res_metagen2)
  
  res <- data.frame(OR = exp(res_metagen2$TE.random),
                    lower_OR = exp(res_metagen2$lower.random),
                    upper_OR = exp(res_metagen2$upper.random),
                    Meta_Pval = res_metagen2$pval.random)       
}else if(res_metagen1$pval.Q > 0.05){
  # P-value of heterogeneity test < 0.05, using fixed effects meta-analysis  
  res_metagen2 <- update(res_metagen1, random = FALSE)
  summary(res_metagen2)
  
  res <- data.frame(OR = exp(res_metagen2$TE.fixed),
                    lower_OR = exp(res_metagen2$lower.fixed),
                    upper_OR = exp(res_metagen2$upper.fixed),
                    Meta_Pval = res_metagen2$pval.fixed)       
}


res_metacont1 <-  metacont(n.e = R_num, 
                           mean.e = R_mean, # Number of events in experimental group or true positives in diagnostic study.
                           sd.e = seTE,
                           n.c = NR_num, 
                           mean.c = NR_mean, # Number of events in control group or false positives in diagnostic study.
                           sd.c = seTE,
                           data = tmp,
                           sm = "MD") # which summary measure ("RR", "OR", "RD", "ASD", or "DOR") is to be used for pooling of studies

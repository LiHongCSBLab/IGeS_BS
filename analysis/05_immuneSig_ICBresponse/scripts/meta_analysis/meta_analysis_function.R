
ICB_combined_Pvalue <- function(resall, 
                                cancer = NULL,
                                treatment = NULL,
                                save_path){

    if(is.null(cancer)){
        immunesig_index <- unique(resall$immune_sig)
        # combined p-value -------------------------------------------------------------
        res_merged_p <- lapply(as.list(immunesig_index), function(x){
            tmp <- resall[resall$immune_sig == x, ]
            if(length(which(is.na(tmp$wilcox_P_value)))>0)tmp = tmp[-which(is.na(tmp$wilcox_P_value)),]
            if(length(which(is.na(tmp$R_greater_NR_P_value)))>0)tmp = tmp[-which(is.na(tmp$R_greater_NR_P_value)),]
            if(length(which(is.na(tmp$R_less_NR_P_value)))>0)tmp = tmp[-which(is.na(tmp$R_less_NR_P_value)),]

            print(x)
            # 'wilcox_P_value', 'R_greater_NR_P_value', 'R_less_NR_P_value'
            merged_p_naive <- nrow(tmp[tmp$wilcox_P_value < 0.05,])/nrow(tmp)
            merged_p_fisher <- merge_p_values(tmp$wilcox_P_value, method = "Fisher")
            merged_p_ACAT <- ACAT(Pvals=tmp[tmp$wilcox_P_value != 0 & tmp$wilcox_P_value != 1, ]$wilcox_P_value)

            merged_RgreaterNR_p_naive <- nrow(tmp[tmp$R_greater_NR_P_value<0.05,])/nrow(tmp)
            merged_RgreaterNR_p_fisher <- merge_p_values(tmp$R_greater_NR_P_value, method = "Fisher")
            merged_RgreaterNR_p_ACAT <- ACAT(Pvals=tmp[tmp$R_greater_NR_P_value != 0 & tmp$R_greater_NR_P_value != 1, ]$R_greater_NR_P_value)

            merged_RlessNR_p_naive <- nrow(tmp[tmp$R_less_NR_P_value < 0.05,])/nrow(tmp)
            merged_RlessNR_p_fisher <- merge_p_values(tmp$R_less_NR_P_value, method = "Fisher")
            merged_RlessNR_p_ACAT <- ACAT(Pvals=tmp[tmp$R_less_NR_P_value != 0 & tmp$R_less_NR_P_value != 1, ]$R_less_NR_P_value)

            res = data.frame(merged_p_naive = merged_p_naive,
                            merged_p_fisher = merged_p_fisher,
                            merged_p_ACAT = merged_p_ACAT,
                            merged_RgreaterNR_p_naive = merged_RgreaterNR_p_naive,
                            merged_RgreaterNR_p_fisher = merged_RgreaterNR_p_fisher,
                            merged_RgreaterNR_p_ACAT = merged_RgreaterNR_p_ACAT,
                            merged_RlessNR_p_naive = merged_RlessNR_p_naive,
                            merged_RlessNR_p_fisher = merged_RlessNR_p_fisher,
                            merged_RlessNR_p_ACAT = merged_RlessNR_p_ACAT)
            return(res)

        })
        names(res_merged_p) <- immunesig_index
        res_merged_p <- do.call(rbind, res_merged_p)
        write.csv(res_merged_p, paste0(save_path, "res_merged_p.csv"), quote = F, row.names = T)

    }else if(!is.null(cancer)) {
        res_cancer <- resall[resall$cancer == cancer, ]
        immunesig_index <- unique(res_cancer$immune_sig)
        # combined p-value -------------------------------------------------------------
        res_merged_p <- lapply(as.list(immunesig_index), function(x){
            tmp <- res_cancer[res_cancer$immune_sig == x, ]
            print(x)
            # 'wilcox_P_value', 'R_greater_NR_P_value', 'R_less_NR_P_value'
            if(length(which(is.na(tmp$wilcox_P_value)))>0)tmp = tmp[-which(is.na(tmp$wilcox_P_value)),]
            if(length(which(is.na(tmp$R_greater_NR_P_value)))>0)tmp = tmp[-which(is.na(tmp$R_greater_NR_P_value)),]
            if(length(which(is.na(tmp$R_less_NR_P_value)))>0)tmp = tmp[-which(is.na(tmp$R_less_NR_P_value)),]
            
            if(nrow(tmp) > 0){

            merged_p_naive <- nrow(tmp[tmp$wilcox_P_value < 0.05,])/nrow(tmp)
            merged_p_fisher <- merge_p_values(tmp$wilcox_P_value, method = "Fisher")
            merged_p_ACAT <- ACAT(Pvals=tmp[tmp$wilcox_P_value != 0 & tmp$wilcox_P_value != 1, ]$wilcox_P_value)

            merged_RgreaterNR_p_naive <- nrow(tmp[tmp$R_greater_NR_P_value<0.05,])/nrow(tmp)
            merged_RgreaterNR_p_fisher <- merge_p_values(tmp$R_greater_NR_P_value, method = "Fisher")
            merged_RgreaterNR_p_ACAT <- ACAT(Pvals=tmp[tmp$R_greater_NR_P_value != 0 & tmp$R_greater_NR_P_value != 1, ]$R_greater_NR_P_value)

            merged_RlessNR_p_naive <- nrow(tmp[tmp$R_less_NR_P_value < 0.05,])/nrow(tmp)
            merged_RlessNR_p_fisher <- merge_p_values(tmp$R_less_NR_P_value, method = "Fisher")
            merged_RlessNR_p_ACAT <- ACAT(Pvals=tmp[tmp$R_less_NR_P_value != 0 & tmp$R_less_NR_P_value != 1, ]$R_less_NR_P_value)

            res = data.frame(merged_p_naive = merged_p_naive,
                            merged_p_fisher = merged_p_fisher,
                            merged_p_ACAT = merged_p_ACAT,
                            merged_RgreaterNR_p_naive = merged_RgreaterNR_p_naive,
                            merged_RgreaterNR_p_fisher = merged_RgreaterNR_p_fisher,
                            merged_RgreaterNR_p_ACAT = merged_RgreaterNR_p_ACAT,
                            merged_RlessNR_p_naive = merged_RlessNR_p_naive,
                            merged_RlessNR_p_fisher = merged_RlessNR_p_fisher,
                            merged_RlessNR_p_ACAT = merged_RlessNR_p_ACAT)
            }else{

            res = data.frame(merged_p_naive = NA,
                            merged_p_fisher = NA,
                            merged_p_ACAT = NA,
                            merged_RgreaterNR_p_naive = NA,
                            merged_RgreaterNR_p_fisher = NA,
                            merged_RgreaterNR_p_ACAT = NA,
                            merged_RlessNR_p_naive = NA,
                            merged_RlessNR_p_fisher = NA,
                            merged_RlessNR_p_ACAT = NA)
            }
        return(res)                

        })
        names(res_merged_p) <- immunesig_index
        res_merged_p <- do.call(rbind, res_merged_p)
        write.csv(res_merged_p, paste0(save_path, cancer, "_res_merged_p.csv"), quote = F, row.names = T)

    }

    return(res_merged_p)
}

# meta analysis ----------------------------------------------------------------

ICB_metaAnalysis <- function(resall, 
                             cancer = NULL,
                             treatment = NULL,
                             save_path){
    
    library(meta)
    
    if(is.null(cancer)){
        immunesig_index <- unique(resall$immune_sig)
        # combined p-value -------------------------------------------------------------
        which(immunesig_index== "dc_sig160")
        # immunesig_index[-c(98, 109, 122,200,321)]
        # ERROR:  Fisher scoring algorithm did not converge. See 'help(rma)' for possible remedies
        
        res_metaRes <- lapply(as.list(immunesig_index), function(x){ # [-c(98, 109, 122,200,321)]
            # x = immunesig_index[4]
            tmp <- resall[resall$immune_sig == x, ]
            # tmp <- tmp[c(1:4,14,15)]
            print(x)
            # if(length(which(is.na(tmp$wilcox_P_value)))>0)tmp = tmp[-which(is.na(tmp$wilcox_P_value)),]
            # if(length(which(is.na(tmp$R_greater_NR_P_value)))>0)tmp = tmp[-which(is.na(tmp$R_greater_NR_P_value)),]
            # if(length(which(is.na(tmp$R_less_NR_P_value)))>0)tmp = tmp[-which(is.na(tmp$R_less_NR_P_value)),]
          
          
            # # # print(boxplot.stats(tmp$seTE)$out)
            # # print(which(tmp$seTE %in% boxplot.stats(tmp$seTE)$out))
            # if(length(boxplot.stats(tmp$seTE)$out)>0){
            #     if(boxplot.stats(tmp$seTE)$out > 1){
            #         tmp = tmp[-which(tmp$seTE %in% boxplot.stats(tmp$seTE)$out),]
            #     }
            # }
            # if(nrow(tmp[tmp$seTE > 1,])>0) tmp <- tmp[tmp$seTE <= 1,]
            res_metagen1 <- metagen(TE, seTE, studlab = dataset, sm = "OR", data = tmp,
                                    control=list(stepadj=0.5, maxiter=10000))


            # if(res_metagen1$pval.Q < 0.05){
                # P-value of heterogeneity test < 0.05, using random effects meta-analysis  
                res_metagen2 <- update(res_metagen1, fixed = FALSE)
                summary(res_metagen2)

                res <- data.frame(OR = exp(res_metagen2$TE.random),
                                lower_OR = exp(res_metagen2$lower.random),
                                upper_OR = exp(res_metagen2$upper.random),
                                Meta_Pval = res_metagen2$pval.random)       
            # }else if(res_metagen1$pval.Q > 0.05){
            #     # P-value of heterogeneity test < 0.05, using fixed effects meta-analysis  
            #     res_metagen2 <- update(res_metagen1, random = FALSE)
            #     summary(res_metagen2)

            #     res <- data.frame(OR = exp(res_metagen2$TE.fixed),
            #                     lower_OR = exp(res_metagen2$lower.fixed),
            #                     upper_OR = exp(res_metagen2$upper.fixed),
            #                     Meta_Pval = res_metagen2$pval.fixed)       
            # # }

            return(res)

        })

        names(res_metaRes) <- immunesig_index # [-c(98, 109, 122,200,321)]
        res_metaRes <- do.call(rbind, res_metaRes)
        nrow(res_metaRes[res_metaRes$Meta_Pval <0.05, ])
     
        write.csv(res_metaRes, paste0(save_path, "res_metaRes.csv"), quote = F, row.names = T)

    }else if(!is.null(cancer)) {
         res_cancer <- resall[resall$cancer == cancer, ]
         res_cancer[is.na(res_cancer)] = 0
         immunesig_index <- unique(res_cancer$immune_sig)
        # combined p-value -------------------------------------------------------------
        # which(immunesig_index== "dc_sig160")
        # immunesig_index[-c(98, 109, 122,200,321)]
        # ERROR:  Fisher scoring algorithm did not converge. See 'help(rma)' for possible remedies
        
        res_metaRes <- lapply(as.list(immunesig_index), function(x){ # [-c(98, 109, 122,200,321)]
            tmp <- res_cancer[res_cancer$immune_sig == x, ]
            
            print(x)
            if(length(which(is.na(tmp$wilcox_P_value)))>0)tmp = tmp[-which(is.na(tmp$wilcox_P_value)),]
            if(length(which(is.na(tmp$R_greater_NR_P_value)))>0)tmp = tmp[-which(is.na(tmp$R_greater_NR_P_value)),]
            if(length(which(is.na(tmp$R_less_NR_P_value)))>0)tmp = tmp[-which(is.na(tmp$R_less_NR_P_value)),]
            # print(boxplot.stats(tmp$seTE)$out)
            # print(which(tmp$seTE %in% boxplot.stats(tmp$seTE)$out))
            if(length(boxplot.stats(tmp$seTE)$out)>0){
                if(boxplot.stats(tmp$seTE)$out > 100){
                    tmp = tmp[-which(tmp$seTE %in% boxplot.stats(tmp$seTE)$out),]
                }
            }
            res_metagen1 <- metagen(TE, seTE, studlab = dataset, sm = "OR",data = tmp,
                                    control=list(stepadj=0.5, maxiter=10000))
            # if(res_metagen1$pval.Q < 0.05){
                # P-value of heterogeneity test < 0.05, using random effects meta-analysis  
                res_metagen2 <- update(res_metagen1, fixed = FALSE)
                summary(res_metagen2)

                res <- data.frame(OR = exp(res_metagen2$TE.random),
                                lower_OR = exp(res_metagen2$lower.random),
                                upper_OR = exp(res_metagen2$upper.random),
                                Meta_Pval = res_metagen2$pval.random)       
            # }else if(res_metagen1$pval.Q > 0.05){
            #     # P-value of heterogeneity test < 0.05, using fixed effects meta-analysis  
            #     res_metagen2 <- update(res_metagen1, random = FALSE)
            #     summary(res_metagen2)

            #     res <- data.frame(OR = exp(res_metagen2$TE.fixed),
            #                     lower_OR = exp(res_metagen2$lower.fixed),
            #                     upper_OR = exp(res_metagen2$upper.fixed),
            #                     Meta_Pval = res_metagen2$pval.fixed)       
            # }

            return(res)

        })

        names(res_metaRes) <- immunesig_index # [-c(98, 109, 122,200,321)]
        res_metaRes <- do.call(rbind, res_metaRes)
        nrow(res_metaRes[res_metaRes$Meta_Pval <0.05, ])
     
        write.csv(res_metaRes, paste0(save_path, cancer, "_res_metaRes.csv"), quote = F, row.names = T)
    }
    return(res_metaRes)

}



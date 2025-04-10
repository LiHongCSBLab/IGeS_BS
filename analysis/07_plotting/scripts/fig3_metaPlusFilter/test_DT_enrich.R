ITS_p_s_df <- read.table("07_plotting_v2/03_ITScharacteristic_metaPlusFilter/ITSp_allgene/ITS_p_s_df.txt", sep = '\t', header = T)
dt_s <- list(drugTarget_all=unique(drugTarget_all$genename),
             drugTarget_FDA=unique(drugTarget_FDA$genename),
             drugTarget_FDAcancer=unique(drugTarget_FDA[which(drugTarget_FDA$ifCancer_drugbank == 'yes'),]$genename),
             ITS_s_gene=unique(ITS_p_s_df$gene_name))

dt_s <- lapply(dt_s, function(x)if(length(which(is.na(x)))>0){x[-which(is.na(x))]}else{x})

length(ITS_p_s_df[ITS_p_s_df$weighted_rate > 0.01,]$gene_name)
length(intersect(ITS_p_s_df[ITS_p_s_df$weighted_rate > 0.01,]$gene_name, unique(drugTarget_FDA$genename)))
length(unique(drugTarget_FDA$genename))
length(intersect(ITS_p_s_df[ITS_p_s_df$weighted_rate > 0.01,]$gene_name, unique(drugTarget_FDA[which(drugTarget_FDA$ifCancer_drugbank == 'yes'),]$genename)))
length(unique(drugTarget_FDA[which(drugTarget_FDA$ifCancer_drugbank == 'yes'),]$genename))

num_Gene = 19960 # length(unique(unlist(dt_s)))
num_sITS = 5283+633+182
num_DT = 633+182+208+1670
num_DTcancer = 182+208
num_sITS_DT = 633+182
num_sITS_DTcancer = 182

phyper(num_sITS_DT-1, num_sITS, num_Gene-5283-633,  num_DT, lower.tail=F)
phyper(num_sITS_DTcancer-1, num_sITS, num_Gene-5283-633-182, num_DTcancer, lower.tail=F)


phyper(265-1, 1740+265, 20000-1740,  540+265,lower.tail=F)


dat <- data.frame(
  "dtc"      = c(182, 208),
  "none_dtc" = c(5283 + 633, num_Gene - (5283+633+182+208)),
  row.names = c("sITS", "Non-sITS"),
  stringsAsFactors = FALSE
)
# fisher.test(dat)
chisq.test(dat)

dat <- data.frame(
  "dtc"      = c(633+182, 208+1670),
  "none_dtc" = c(5283, num_Gene-(5283+633+182+208+1670)),
  row.names = c("sITS", "Non-sITS"),
  stringsAsFactors = FALSE
)
# fisher.test(dat)
chisq.test(dat)






gs_dt = rbind(data.frame(gs_name = "drugTarget_all", SYMBOL = unique(drugTarget_all$genename)),
              data.frame(gs_name = "drugTarget_FDA", SYMBOL = unique(drugTarget_FDA$genename)),
              data.frame(gs_name = "drugTarget_FDAcancer", SYMBOL = unique(drugTarget_FDA[which(drugTarget_FDA$ifCancer_drugbank == 'yes'),]$genename)))

gs_all = unique(rbind(gs_all, gs_dt))

tmp2 = as.vector(scale(ITS_p_s_df3$weighted_rate))
names(tmp2) = ITS_p_s_df3$gene_name

res = ITS_GS_enrich(tmp2,
                    gs=gs_dt,
                    pvalueCutoff = 1,
                    type = "custom",
                    method = 'GSEA')
data.frame(res)[c(1,6,7)]
res2=data.frame(res)

gseaplot2(res, geneSetID = c(2,3,1), pvalue_table = TRUE,  ES_geom = "dot")

library(ggplot2)
library(cowplot)

pp <- lapply(1:3, function(i) {
    anno <- res[i, c("NES", "pvalue", "p.adjust")]
    lab <- paste0(names(anno), "=",  round(anno, 3), collapse="\n")

    gsearank(res, i, res[i, 2]) + xlab(NULL) +ylab(NULL) +
        annotate("text", 10000, res[i, "enrichmentScore"] * .75, label = lab, hjust=0, vjust=0)
})
plot_grid(plotlist=pp, ncol=1)
cancer = "COAD"
cancer = "READ"
tumor_subtype = "Primary"
purity_method = "TUMERIC"
GEO_ID = "GSE160785"
drugindex = "drug989"
s1=0.4
Pval=1

setwd("/picb/bigdata/project/FengFYM/p2_1_immunotherapy_targetedtherapy/r3_individual_datasets_validation/")

savepath = paste0("result_", GEO_ID, "/drug_immunesig_", purity_method,"_GSEA_result_r",s1,"_meta_oe_original_en_padj",Pval,"_final/", cancer, "_", tumor_subtype, "/", drugindex)

df_gsea_immunesig2 = read.csv(paste0(savepath, "/CellvstreatedMouse.csv"))
dftmp = reshape2::melt(df_gsea_immunesig2, 
                         id.vars = c('immunesig', 'sign', 'x'),
                         measure.vars = c('treated_mouse', "drugindex"),
                         variable.name = c('value'))
  names(dftmp) = c('immunesig','sign','x','label','value')
  
  dftmp = dftmp[order(dftmp$x), ]
  dftmp$label = factor(dftmp$label)
  
  p = ggplot(dftmp, aes(x = x, y = value, color=label, fill=label)) + 
    geom_line() + 
    geom_hline(aes(yintercept = 0)) +
    geom_vline(aes(xintercept = 0)) +
    geom_point(size = 3)+
    theme_bw() +
    theme(axis.text.x  =  element_text(angle=-90, vjust=0.5,hjust = 0, size=12),
          axis.text.y = element_text(size = 15),
          axis.title.x = element_blank())+
    scale_x_continuous(breaks = unique(dftmp[c("x","immunesig")])$x,
                       labels = unique(dftmp[c("x","immunesig")])$immunesig)
  
  p
  ggsave(paste0(savepath, "/CellvstreatedMouse.pdf"),
         p, height = 7, width=5)
  




      df = unique(df_gsea_immunesig2[c('immune_sig','NES_diff', 'flag')])
      row.names(df) = df$immune_sig
      ds = df[grep('sensitive', df$flag), ][2]; names(ds) = "s"
      dr = df[grep('resistant', df$flag), ][2]; names(dr) = "r"
      if(nrow(dr) == 0 ){
        dr = data.frame(r = 0)
        rownames(dr) = "none"
      }
  
      ds1 = (matrix(0,nrow(dr))); rownames(ds1) = rownames(dr)
      dr1 = (matrix(0,nrow(ds))); rownames(dr1) = rownames(ds)
      ds2 = rbind(as.matrix(ds), ds1);colnames(ds2) = "s"
      dr2 = rbind(dr1, as.matrix(dr));colnames(dr2) = "r"
      
      df <- as.data.frame(cbind(ds2, dr2))

      df$max <- ceiling(max(c(max(df$s), max(df$r))))
      df$min <- floor(min(c(min(df$s), min(df$r))))

      df$mean <- 0
      df <- df[c("max","min","r","s","mean")]
      # df$V1 <-  (df$V1 - mean(df$V1))/sd(df$V1)
      df <- as.data.frame(t(df))
      
      pdf(paste0(savepath, "/", drugindex, "_treatedMouse.pdf"), width = 6, height = 6.5)
      radarchart(df,
                 axistype = 4,
                 # pfcol = c(rgb(0.2,0.2,0.2,0.9,0.9),NA,NA),
                 pfcol = c(NA,NA, NA),
                 pcol= c(rgb(0.2,0.5,0.5,0.5,0.9),
                         rgb(0.7,0.2,0.2,0.5,0.9),
                         rgb(0.22,0.22,0.22,0.5)),
                 plty = 1, 
                 plwd = 4,
                 # Customize the grid
                 cglcol = "grey70",  cglwd = 0.8,
                 # Customize the axis
                 axislabcol = "black", 
                 # Variable labels
                 cglty = 1,
                 # caxislabels = seq(floor(min(c(min(mat_s3), min(mat_r3)))),
                 #                   ceiling(max(c(max(mat_s3), max(mat_r3)))),
                 #                   3),
                 caxislabels =  c(seq(-4,4,2)),
                 vlcex = 0.6,
                 title = "Impact on Immune")
      # abline(h =0, v = 0, lty = 2)
      dev.off()




# Download data from BROAD Single cell portal, 
# https://portals.broadinstitute.org/single_cell/study/SCP11/melanoma-intra-tumor-heterogeneity

data = read.table('GSE72056_melanoma_single_cell_revised_v2.txt',header=T)

cancerx = (data[2,]==2)  # malignant
stromax = (data[2,]==1)  # malignant

plot.gene = function(g) {
  cancer= unlist(data[which(data[,1]==g),cancerx])
  stroma= unlist(data[which(data[,1]==g),stromax])
  wp = wilcox.test(cancer,stroma)$p.value
  boxplot(cancer,stroma,names=c(paste0('cancer, N=',length(cancer)),paste0('stroma, N=',length(stroma))),main=g,ylab='expression',sub=paste('p-value : ',wp))
}

par(mfrow=c(2,2))

plot.gene('ACVR2B')
plot.gene('LRP6')
plot.gene('BMPR1B')
plot.gene('SEMA6A')


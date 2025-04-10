library(gplots)
library(missMDA)
library(RColorBrewer)
library("colorspace")

# ignores columns c()
quantile_norm = function(df,ignore=c(3)){
	df_rank <- apply(df,2,rank,ties.method="min")
	df_sorted <- data.frame(apply(df, 2, sort))
	df_mean <- apply(df_sorted, 1, mean)
	if (length(ignore)>0) { df_mean <- apply(df_sorted[,-ignore], 1, mean) }
	df_final = apply(df_rank,2,function(x) df_mean[x])
	rownames(df_final) <- rownames(df)
	return(df_final)
}

#setwd('/home/umesh/b/pro/purity/tcga/ana/ana1/')

x0= read.csv('purity.csv', header=TRUE, row.names=1)
dim(x0)

x0[x0>.98]=NA
x0[x0<.10]=NA

# remove rows with all NA
x1=x0[rowSums(is.na(x0))!=4, ]
dim(x1)

#'KIRP','THYM',
#l=as.list(strsplit('BLCA BRCA', " "))
for (i in c('BLCA', 'BRCA', 'CESC', 'CRC', 'ESCA', 'GBM', 'HNSC', 'KIRC', 'KIRP', 'LGG', 'LIHC', 'LUAD', 'LUSC', 'OV', 'PAAD', 'PRAD', 'SKCM', 'STAD', 'THCA', 'UCEC') ){ 

print(i)

pdf(paste(i,'_norm.pdf',sep=''))

# remove ascat column 3 NA data in case of KIRP and THYM
if (i=='KIRP'){ x=x1[ x1$type == i, ][c(2,4,5)] }

else { x=x1[ x1$type == i, ][2:5] }

# set value threshold

# impute missing data
imp = imputePCA(x,ncp=1)

x_imp = imp$completeObs


#change colum and row names
#colnames(x_imp) = colnames(x)
#rownames(x_imp) = rownames(x)

x_cor = cor(x,use='pairwise.complete.obs')
x_cor_imp = cor(x_imp,use='pairwise.complete.obs')
mean_cor_imp = apply(x_cor_imp,1,sum)/3

print('mean_cor_imp')
mean_cor_imp
# purbayes  abscnseq  estimate 
#0.7160241 0.6568832 0.6592943 

#my_palette = colorRampPalette(c("green", "yellow", "red"))(n = 299)


heatmap.2(x_imp, 
#		cellnote = x, # same data set for cell labels
		main = "Correlation", # heat map title
#		notecol="black",
		# change font color of cell labels to black
		density.info="none", # turns off density plot inside color legend
		trace="none",
		# turns off trace lines inside the heat map
		margins =c(8,9),
		# widens margins around plot
#		col=my_palette,
		# use on color palette defined earlier
		dendrogram="col",
		# only draw a row dendrogram
		Rowv=FALSE,
		# turn off column clustering
		distfun=function(x) as.dist((1-cor(t(x),use='pairwise.complete.obs'))/2),
		hclust=function(x) hclust(x,method="complete"))



# too big, ignore
#boxplot(t(x_imp),ylab='purity',xlab='samples')


boxplot(x_imp,   ylab='purity',xlab='methods')



# column shifts one left in case of KIRP and THYM
if (i=='KIRP' || i=='THYM' ){ x_imp_qn = quantile_norm(x_imp,ignore=c(2)) }

else {x_imp_qn = quantile_norm(x_imp) }




heatmap.2(x_imp_qn,
#		cellnote = x, # same data set for cell labels
		main = "Correlation", # heat map title
		notecol="black",
		# change font color of cell labels to black
		density.info="none", # turns off density plot inside color legend
		trace="none",
		# turns off trace lines inside the heat map
		margins =c(8,9),
		# widens margins around plot
#		col=my_palette,
		# use on color palette defined earlier
		dendrogram="col",
		# only draw a row dendrogram
		Rowv=FALSE,# turn off column clustering
		distfun=function(x) as.dist((1-cor(t(x),use='pairwise.complete.obs'))/2),
		hclust=function(x) hclust(x,method="complete"))


#too big
boxplot(t(x_imp_qn),ylab='purity',xlab='samples')

boxplot(x_imp_qn,ylab='purity',xlab='methods')


plot(apply(x_imp,1,mean),apply(x_imp_qn,1,mean),xlim=c(0,1),ylim=c(0,1))


#add straight line to plot
abline(a=0,b=1)

# PCA
# imputed data
mean_pur_imp = apply(x_imp_qn,1,mean)
median_pur_imp = apply(x_imp_qn,1,median)
color=diverge_hcl(length(mean_pur_imp))[rank(mean_pur_imp)]



boxplot(mean_pur_imp,xlab='Estimated tumor purity (% cancer cells)',horizontal = TRUE)



pc = prcomp(x_imp, center = TRUE, scale = TRUE)
plot(pc$x[, 1], pc$x[, 2], main = "PCA", xlab = "PC1", ylab = "PC2", col = color,pch=20)

plot(mean_pur_imp,pc$x[, 1]) # almost perfect correlation ...

plot(median_pur_imp,pc$x[, 1]) # almost perfect correlation ...


write.table(x, paste('thres',i,sep='_'),quote = FALSE )
write.table(x_imp,paste('imp',i,sep='_'),quote = FALSE)
write.table(x_imp_qn,paste('imp_qn',i,sep='_'),quote = FALSE)
write.table(mean_pur_imp,paste('purity',i,sep='_'),quote = FALSE)

}


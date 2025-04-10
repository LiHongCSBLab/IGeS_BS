import pandas as pd
from scipy.optimize import nnls
import scipy.stats as ss

def deconv(ex,pu,cn):	
	pu['stroma']=1-pu
	# mean purity
	pum=pu.iloc[:,0].mean()

	with open('deconv.txt','w') as f:
		f.write('cancer stroma pur_exp_corr cnv2_exp_corr\n')

		# deconvolute each gene
		for g in ex.index: 
			exg=ex.ix[g]	
			flag=0
			
			if g in cn.index: 		
				# cnv data of the gene
				cng=cn.ix[g]		
					
				cng0=cng[cng!=2].index	
				# cnv!=2 samples exp	
				exg0=exg.ix[cng0]
								
				cng2=cng[cng==2].index	
				# cnv==2 samples exp
				exg2=exg.ix[cng2]		
								
				if len(cng2)>4:	
					# exp wilcoxon test for cnv0 and cnv2
					pval=ss.ranksums(exg0, exg2)[1]	
					if pval<5e-6:
						flag=1
						
			if flag==0:						
				d =nnls(pu,exg)[0]
				cor=ss.pearsonr(pu.iloc[:,0], exg)[0]
				f.write( '%s %.4f %.4f %.3f %s\n' % (g, d[0], d[1], cor, 'NA') )

			elif flag==1:
				# remove cnv2 samples to calc stroma
				d =nnls(pu.ix[cng0],exg0)[0]
				cor=ss.pearsonr(pu.iloc[:,0].loc[cng0],exg0)[0]
													
				# C: (T-S+PS)/P
				c= (exg.mean() - d[1] + pum*d[1])/pum			

				f.write( '%s %.4f %.4f %.3f %.1e\n' % (g, c, d[1], cor , pval) )


# copy number data
cn=pd.read_csv('data/all_thresholded.by_genes.txt.bz2', index_col=0, sep='\t').rename(columns=lambda x:x[:15])

# tumor purity
pu=pd.read_csv('data/GBM_purity.csv', index_col=0)

# expression log2(x+1) transformed
ex=pd.read_csv('data/GBM_logx+1.csv.bz2', index_col=0)
	
# extract overlapping samples
o=cn.columns.intersection(pu.index)	
pu=pu.ix[o]
ex=ex[o]
cn=cn[o]
	
# do deconvolution	
deconv(ex,pu,cn)
	
	


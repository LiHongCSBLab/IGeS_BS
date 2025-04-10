import pandas as pd
from scipy.optimize import nnls
from sklearn.utils import resample

# tumor purity
pu=pd.read_csv('data/GBM_purity.csv', index_col=0)
pu['stroma']=1-pu

# expression log2(x+1) transformed
ex=pd.read_csv('data/GBM_logx+1.csv.bz2', index_col=0).T

o=pd.DataFrame(index=ex.columns)

for i in range(5):
	print 'bootstrap #',i
	ex1,pu1=resample(ex,pu,random_state=i) #
	
	ca=[]
	st=[]

	for g in ex.columns:
		d=nnls(pu1, ex1[g])[0]
		ca.append(d[0])
		st.append(d[1])		
	o['C'+str(i)]=ca
	o['S'+str(i)]=st
#o.astype('float16').to_pickle('boot/'+ty+'.pkl.gz')	

out=pd.DataFrame(index=ex.columns)
ca=o.filter(regex='C')
st=o.filter(regex='S')
out['C50']=ca.quantile(.50,axis=1)
out['C95']=ca.quantile(.95,axis=1)
out['C05']=ca.quantile(.05,axis=1)
out['S50']=st.quantile(.50,axis=1)
out['S95']=st.quantile(.95,axis=1)
out['S05']=st.quantile(.05,axis=1)

out.to_csv('output/GBM_deconv_bootstrap.csv.bz2')

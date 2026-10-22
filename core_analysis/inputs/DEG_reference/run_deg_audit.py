import pandas as pd, numpy as np, gzip, zipfile, re, os
from scipy.stats import t as tdist
from statsmodels.stats.multitest import multipletests
base='/mnt/data/SONFH_DEG_Comparison'
# read normalized gene-level log2 expression
expr=pd.read_csv(base+'/expression.csv.gz')
manifest=pd.read_csv(base+'/manifest.csv')
# align manifest to expression columns
samples=[c for c in expr.columns if c.startswith('GSM')]
manifest=manifest.set_index('GSM').loc[samples].reset_index()
assert len(samples)==40 and len(manifest)==40
# Group coding: SONFH=1, exposed non-SONFH=0; sex male=1 female=0
group=(manifest['group']=='SONFH').astype(float).to_numpy()
sex=(manifest['sex']=='Male').astype(float).to_numpy()
X=np.column_stack([np.ones(len(samples)),group,sex])
Y=expr[samples].to_numpy(dtype=float) # genes x samples
# OLS estimates for all genes
XtX_inv=np.linalg.inv(X.T@X)
B=(XtX_inv@X.T@Y.T).T # genes x 3
resid=Y.T - X@B.T # samples x genes
H=np.sum((X@XtX_inv)*X,axis=1)
# HC3 sandwich: meat = X' diag(e^2/(1-h)^2) X, vectorized over genes
adj=(resid/(1-H[:,None]))**2
se=np.empty(Y.shape[0])
for j in range(Y.shape[0]):
    meat=X.T@(X*adj[:,j,None])
    cov=XtX_inv@meat@XtX_inv
    se[j]=np.sqrt(max(cov[1,1],0))
beta=B[:,1]
tstat=beta/se
pvals=2*tdist.sf(np.abs(tstat),df=len(samples)-X.shape[1])
fdr=multipletests(pvals,method='fdr_bh')[1]
res=pd.DataFrame({'gene':expr['Gene_Symbol'].astype(str),'beta_SONFH_vs_nonSONFH':beta,'HC3_SE':se,'t':tstat,'df':len(samples)-X.shape[1],'P_value':pvals,'FDR_BH_all_genes':fdr})
res=res.sort_values('P_value')
res.to_csv(base+'/DEG_all_genes_HC3_BH.csv',index=False)
# extract catabolism genes from frozen GMT
with zipfile.ZipFile('/mnt/data/SONFH_Research.zip') as z:
    gmt=z.read('01_Frozen_Master/SONFH_QC04_frozen_85_members_full.gmt').decode('utf-8')
line=[l for l in gmt.splitlines() if l.startswith('GOBP_HYDROGEN_PEROXIDE_CATABOLIC_PROCESS\t')][0]
members=line.split('\t')[2:]
matched=set(expr['Gene_Symbol'].astype(str))
cat=pd.DataFrame({'gene':members})
cat['measurable_in_final_matrix']=cat['gene'].isin(matched)
cat=cat.merge(res,on='gene',how='left')
cat=cat.sort_values('P_value',na_position='last')
cat.to_csv(base+'/DEG_hydrogen_peroxide_catabolism_members.csv',index=False)
# biosynthesis members too
line2=[l for l in gmt.splitlines() if l.startswith('GOBP_HYDROGEN_PEROXIDE_BIOSYNTHETIC_PROCESS\t')][0]
mem2=line2.split('\t')[2:]
bio=pd.DataFrame({'gene':mem2})
bio['measurable_in_final_matrix']=bio['gene'].isin(matched)
bio=bio.merge(res,on='gene',how='left').sort_values('P_value',na_position='last')
bio.to_csv(base+'/DEG_hydrogen_peroxide_biosynthesis_members.csv',index=False)
# summary
summary={
'n_samples':len(samples),'n_genes_tested':len(expr),'n_nominal_p_lt_005':int((res.P_value<.05).sum()),
'n_FDR_lt_005':int((res.FDR_BH_all_genes<.05).sum()),
'catabolism_n_listed':len(members),'catabolism_n_measurable':int(cat.measurable_in_final_matrix.sum()),
'catabolism_n_negative_beta':int((cat.beta_SONFH_vs_nonSONFH<0).sum()),
'catabolism_n_nominal_p_lt_005':int((cat.P_value<.05).sum()),
'catabolism_n_FDR_lt_005_global':int((cat.FDR_BH_all_genes<.05).sum()),
'biosynthesis_n_listed':len(mem2),'biosynthesis_n_measurable':int(bio.measurable_in_final_matrix.sum()),
'biosynthesis_n_negative_beta':int((bio.beta_SONFH_vs_nonSONFH<0).sum()),
'biosynthesis_n_nominal_p_lt_005':int((bio.P_value<.05).sum()),
'biosynthesis_n_FDR_lt_005_global':int((bio.FDR_BH_all_genes<.05).sum())}
pd.DataFrame([summary]).to_csv(base+'/DEG_summary.csv',index=False)
print(summary)
print('\nCatabolism measurable genes:')
print(cat[cat.measurable_in_final_matrix][['gene','beta_SONFH_vs_nonSONFH','P_value','FDR_BH_all_genes']].to_string(index=False))
print('\nBiosynthesis measurable genes:')
print(bio[bio.measurable_in_final_matrix][['gene','beta_SONFH_vs_nonSONFH','P_value','FDR_BH_all_genes']].to_string(index=False))

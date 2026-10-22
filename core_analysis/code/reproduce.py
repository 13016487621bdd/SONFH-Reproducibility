"""Portable independent verification from archived scores, not a GSVA reimplementation.
Run: python code/reproduce.py (paths resolve from this file).
Outputs never overwrite archived inputs. Assertions fail on numerical discrepancies.
"""
from pathlib import Path
import json, hashlib, itertools, sys, platform, warnings
import numpy as np
import pandas as pd
from scipy import stats
import scipy

ROOT = Path(__file__).resolve().parents[1]
I, O = ROOT/'inputs', ROOT/'results'
O.mkdir(exist_ok=True)
checks=[]
def check(name, value, tolerance=1e-9):
    err=float(np.max(np.abs(np.asarray(value))))
    checks.append(dict(check=name,max_absolute_error=err,tolerance=tolerance,status='PASS' if err<tolerance else 'FAIL'))
    assert err<tolerance, (name,err)
def bh(p):
    p=np.asarray(p); ix=np.argsort(p); q=np.empty(len(p))
    q[ix]=np.minimum(1,np.minimum.accumulate((p[ix]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1])
    return q
def zrows(a):
    return (a-a.mean(axis=1,keepdims=True))/a.std(axis=1,ddof=1,keepdims=True)
def regress(Y, X):
    """OLS coefficients, exact HC3 sandwich, residual-df Student t inference."""
    Y=np.atleast_2d(Y).astype(float)
    A=np.linalg.inv(X.T@X)@X.T
    b=A@Y.T; residual=Y.T-X@b
    h=np.sum(X*A.T,axis=1)
    se=np.sqrt(np.sum((A[1,:,None]*residual/(1-h[:,None]))**2,axis=0))
    df=X.shape[0]-X.shape[1]; t=b[1]/se; p=2*stats.t.sf(abs(t),df)
    delta=stats.t.ppf(.975,df)*se
    return pd.DataFrame(dict(beta=b[1],HC3_SE=se,t=t,p=p,CI_lower=b[1]-delta,CI_upper=b[1]+delta))
def save(t,name):
    t.to_csv(O/name,index=False)
def match(t,ref,name,cols=('beta','HC3_SE','t','p','CI_lower','CI_upper','FDR')):
    for col in cols: check(name+' / '+col,t[col].to_numpy()-ref[col].to_numpy())

manifest=pd.read_csv(I/'03_QC/SONFH_sample_manifest_aligned.csv').set_index('GSM')
ids=manifest.index.tolist()
assert manifest.index.is_unique and len(ids)==40
assert set(manifest.group)=={'GC_exposed_non_SONFH','SONFH'}
assert set(manifest.sex)=={'Female','Male'}
case=(manifest.group=='SONFH').to_numpy(); male=(manifest.sex=='Male').to_numpy()
assert case.sum()==30 and (~case).sum()==10
X=np.column_stack([np.ones(40),case,male]).astype(float)
assert np.linalg.matrix_rank(X)==3
expr=pd.read_csv(I/'03_QC/SONFH_QC03_gene_expression_median_log2.csv.gz',index_col=0)
assert expr.shape==(18833,40) and expr.index.is_unique
assert list(expr.columns)==ids and np.isfinite(expr.to_numpy()).all()
sets={}
for line in (I/'frozen_85.gmt').read_text().splitlines():
    a=line.split('\t'); assert a[0] not in sets; sets[a[0]]=a[2:]
assert len(sets)==85
scores={}
for method in ['GSVA','ssGSEA']:
    s=pd.read_csv(I/f'04_GeneSet_Scoring/SONFH_{method}_scores_85x40.csv',index_col=0)
    assert s.shape==(85,40) and list(s.columns)==ids and set(s.index)==set(sets)
    scores[method]=s
    r=regress(s.to_numpy(),X);r.insert(0,'gene_set',s.index);r['FDR']=bh(r.p)
    refname=f'SONFH_Layer1_{method}_85_'+('full_results.csv' if method=='GSVA' else 'results.csv')
    ref=pd.read_csv(I/'05_Statistics'/refname).set_index('gene_set').loc[s.index].reset_index()
    match(r,ref,'Layer 1 '+method);save(r,f'layer1_{method}.csv')

pairdef=pd.read_csv(I/'05_Statistics/SONFH_Layer3_11_pairs_GSVA.csv')[['Pair_ID','Set_A','Set_B']]
pairvalues={}; pairresults={}
for method,s in scores.items():
    z=pd.DataFrame(zrows(s.to_numpy()),index=s.index,columns=ids)
    values=np.array([z.loc[a].to_numpy()-z.loc[b].to_numpy() for a,b in zip(pairdef.Set_A,pairdef.Set_B)])
    r=pd.concat([pairdef,regress(values,X)],axis=1);r['FDR']=bh(r.p)
    ref=pd.read_csv(I/f'05_Statistics/SONFH_Layer3_11_pairs_{method}.csv')
    match(r,ref,'Layer 3 '+method);save(r,f'layer3_{method}.csv')
    pairvalues[method]=values;pairresults[method]=r
g=scores['GSVA']; gz=pd.DataFrame(zrows(g.to_numpy()),index=g.index,columns=ids)
cat='GOBP_HYDROGEN_PEROXIDE_CATABOLIC_PROCESS'; bio='GOBP_HYDROGEN_PEROXIDE_BIOSYNTHETIC_PROCESS'
component=regress(gz.loc[[bio,cat]].to_numpy(),X)
component.insert(0,'outcome',['Biosynthesis','Catabolism']);save(component,'components_standardized.csv')
assert abs(component.beta.iloc[0]-component.beta.iloc[1]-pairresults['GSVA'].beta.iloc[-1])<1e-12
individual=manifest[['group','sex']].copy()
individual['biosynthesis_z']=gz.loc[bio];individual['catabolism_z']=gz.loc[cat]
individual['contrast']=individual.biosynthesis_z-individual.catabolism_z
genes=['ALAS2','AHSP','GYPA','SLC4A1']
proxy=zrows(expr.loc[genes].to_numpy()).mean(axis=0)
proxy5=zrows(expr.loc[['HBB']+genes].to_numpy()).mean(axis=0)
individual['E1_four_gene_proxy']=proxy;individual['HBB_five_gene_sensitivity']=proxy5
save(individual.reset_index(),'participant_scores.csv')
old=pd.read_csv(I/'E1_reference/E1_proxy_individual_scores.csv').set_index('GSM').loc[ids]
check('Four-gene proxy vs archived E1',proxy-old.E1_four_gene_proxy)
proxyrows=[]
for outcome,y in [('Catabolism',g.loc[cat].to_numpy()),('Biosynthesis',g.loc[bio].to_numpy()),('Contrast',pairvalues['GSVA'][-1])]:
    for label,cov in [('Sex adjusted',None),('Four-gene proxy adjusted',proxy),('Five-gene proxy sensitivity',proxy5)]:
        r=regress(y,X if cov is None else np.column_stack([X,cov])).iloc[0].to_dict()
        r.update(outcome=outcome,model=label);proxyrows.append(r)
proxy_results=pd.DataFrame(proxyrows);save(proxy_results,'proxy_models.csv')
ref_summary=json.loads((I/'E1_reference/E1_consistency_closure_summary.json').read_text())
for out,key in [('Catabolism','catabolism'),('Biosynthesis','biosynthesis'),('Contrast','pair')]:
    for label,suffix in [('Sex adjusted','primary'),('Four-gene proxy adjusted','adjusted')]:
        r=proxy_results[(proxy_results.outcome==out)&(proxy_results.model==label)].iloc[0]
        for a,b in [('beta','beta'),('HC3_SE','SE_HC3'),('p','P_value_t'),('CI_lower','CI95_low'),('CI_upper','CI95_high')]:
            check('E1 '+out+' '+suffix+' '+a,r[a]-ref_summary[key+'_'+suffix][b])
panels={
 'Neutrophil-associated':['FCGR3B','CEACAM8','CSF3R','S100A8','S100A9'],
 'Monocyte-associated':['CD14','LYZ','FCGR1A','LST1','S100A12'],
 'T-cell-associated':['CD3D','CD3E','CD2','TRAC','IL7R'],
 'B-cell-associated':['MS4A1','CD79A','CD79B','CD19','BANK1'],
 'Erythroid-associated (frozen four genes)':genes,
 'Platelet-associated':['PPBP','PF4','GP9','ITGA2B','TUBB1']}
panelrows=[]
for label,requested in panels.items():
    used=[x for x in requested if x in expr.index];p=zrows(expr.loc[used].to_numpy()).mean(axis=0)
    rho,pval=stats.spearmanr(p,g.loc[cat])
    panelrows.append(dict(panel=label,genes_requested=';'.join(requested),genes_used=';'.join(used),n_genes=len(used),rho_all=rho,p_all=pval,rho_control=stats.spearmanr(p[~case],g.loc[cat].to_numpy()[~case]).statistic,rho_SONFH=stats.spearmanr(p[case],g.loc[cat].to_numpy()[case]).statistic))
save(pd.DataFrame(panelrows),'marker_panels_corrected.csv')
all_adj=regress(g.to_numpy(),np.column_stack([X,proxy]));all_adj.insert(0,'gene_set',g.index);all_adj['FDR']=bh(all_adj.p)
save(all_adj,'exploratory_all85_proxy_adjusted.csv')

removal=pd.read_csv(I/'09_KeyGene_Removal_Audit/participant_scores_by_scenario.csv').set_index('GSM').loc[ids]
rr=regress(removal.iloc[:,2:].to_numpy().T,X)
rr.insert(0,'scenario',removal.columns[2:])
ref=pd.read_csv(I/'09_KeyGene_Removal_Audit/gene_removal_HC3_results.csv').set_index('scenario').loc[rr.scenario]
for a,b in [('beta','beta'),('HC3_SE','HC3_SE'),('p','p_nominal'),('CI_lower','CI_lower'),('CI_upper','CI_upper')]:check('Removal '+a,rr[a].to_numpy()-ref[b].to_numpy())
save(rr,'gene_removal_nominal.csv')
check('Archived removal baseline matches frozen scores',removal.original-g.loc[cat])

# Refit all 40 deletions: original full-cohort Z scores remain frozen.
loob=[];looq=[]
for j in range(40):
    mask=np.arange(40)!=j
    r=regress(g.to_numpy()[:,mask],X[mask]);loob.append(r.beta.to_numpy())
    r=regress(pairvalues['GSVA'][:,mask],X[mask]);looq.append(bh(r.p))
loob=np.array(loob).T;looq=np.array(looq).T
ref=pd.read_csv(I/'05_Statistics/SONFH_Layer1_LOO_beta_85.csv').set_index('gene_set').loc[g.index]
check('Layer 1 LOO minimum',loob.min(axis=1)-ref.beta_LOO_min)
check('Layer 1 LOO maximum',loob.max(axis=1)-ref.beta_LOO_max)
oldq=pd.read_csv(I/'05_Statistics/SONFH_Layer3_LOO_FDR_11x40.csv',index_col=0)
check('All 440 Layer 3 LOO FDR values',looq-oldq.loc[pairdef.Pair_ID,ids].to_numpy())
save(pd.DataFrame(dict(removed=ids,Pair11_FDR=looq[-1])),'pair11_LOO.csv')

# Family CR2 statistics are not recalculated by this Python script.
family=pd.read_csv(I/'05_Statistics/SONFH_Layer2_21_family_omnibus_results.csv')
check('Family BH correction from archived P values',bh(family.p)-family.FDR)
check('Family P values from archived F and df',stats.f.sf(family.statistic,family.df_num,family.df_denom)-family.p)
save(family,'layer2_ARCHIVED_CR2.csv')

# Exact permutation: corrected sample-level tissue scores, not the older CSV.
tissue=pd.read_csv(I/'06_Validation/GSE345158_Pair11_SOD2_corrected_scores.csv')
tcase=tissue.Sample.str.startswith('Case-').to_numpy(); assert tcase.sum()==5
v=tissue.Pair11.to_numpy();effect=v[tcase].mean()-v[~tcase].mean()
perms=np.array([v[list(ix)].mean()-np.delete(v,ix).mean() for ix in itertools.combinations(range(9),5)])
extreme=np.sum(abs(perms)>=abs(effect)-1e-10);p=extreme/126
ref=pd.read_csv(I/'06_Validation/GSE345158_Pair11_overlap_sensitivity.csv').iloc[0]
check('Tissue corrected difference (archive discrepancy <1e-6)',effect-ref.Difference,1e-6)
check('Tissue corrected exact P',p-ref.Exact_P)
save(pd.DataFrame(dict(mean_difference=perms)),'tissue_exact_permutations.csv')
tissue['group']=np.where(tcase,'SONFH','Fracture control');save(tissue,'tissue_corrected_scores.csv')
loo=[]
for j in range(9):
    keep=np.arange(9)!=j
    loo.append(dict(removed=tissue.Sample.iloc[j],difference=v[keep&tcase].mean()-v[keep&~tcase].mean()))
save(pd.DataFrame(loo),'tissue_LOO_corrected_fixed_scores.csv')

# Genome-wide HC3 gene associations, with BH across all 18,833 genes.
gene=regress(expr.to_numpy(),X);gene.insert(0,'gene',expr.index);gene['FDR']=bh(gene.p)
ref=pd.read_csv(I/'DEG_reference/DEG_all_genes_HC3_BH.csv').set_index('gene').loc[expr.index]
for a,b in [('beta','beta_SONFH_vs_nonSONFH'),('p','P_value'),('FDR','FDR_BH_all_genes')]:check('Genome-wide gene '+a,gene[a].to_numpy()-ref[b].to_numpy())
save(gene,'genomewide_HC3_BH.csv')
save(gene[gene.gene.isin(set(sets[bio]+sets[cat]))],'peroxide_member_gene_results.csv')

# Check archived R objects against the portable inputs, without executing R.
import rdata
with warnings.catch_warnings():
    warnings.simplefilter('ignore');param=rdata.read_rds(I/'04_GeneSet_Scoring/SONFH_GSVA_parameters.rds')
check('GSVA parameter object expression matrix',param.exprData.values-expr.to_numpy())
assert list(param.exprData.coords['dim_0'].values)==expr.index.tolist()
assert list(param.exprData.coords['dim_1'].values)==ids
assert all(set(param.geneSets[k])==set(sets[k]) for k in sets)
assert param.kcdf[0]=='Gaussian' and param.tau[0]==1 and param.maxDiff[0] and not param.absRanking[0]
tobj=rdata.read_rds(I/'06_Validation/GSE345158_reproducibility_inputs.rds')
tm=tobj['gsva_matrix_sod2']; assert tm.shape==(18573,9)
assert 'SOD2' in tm.coords['dim_0'].values
pd.DataFrame(tm.values,index=tm.coords['dim_0'].values,columns=tm.coords['dim_1'].values).to_csv(O/'tissue_transformed_input_export.csv.gz')
save(pd.DataFrame(checks),'verification_checks.csv')
summary=dict(samples=40,cases=30,controls=10,genes=18833,gene_sets=85,
 layer1_significant=int((pd.read_csv(O/'layer1_GSVA.csv').FDR<.05).sum()),
 layer2_significant_archived=int((family.FDR<.05).sum()),layer3_significant=int((pairresults['GSVA'].FDR<.05).sum()),
 pair11_GSVA=pairresults['GSVA'].iloc[-1].to_dict(),pair11_ssGSEA=pairresults['ssGSEA'].iloc[-1].to_dict(),
 proxy_genes=genes,proxy_catabolism_rho=float(stats.spearmanr(proxy,g.loc[cat]).statistic),
 proxy_adjusted_85_significant=int((all_adj.FDR<.05).sum()),
 tissue_corrected_effect=float(effect),tissue_exact_extreme=int(extreme),tissue_exact_allocations=126,tissue_exact_P=float(p),
 numerical_checks_passed=len(checks),
 execution_scope='Independent Python recalculation from frozen scores and matrices; no R runtime in this environment. GSVA/ssGSEA scoring and CR2 family inference NOT rerun here.',
 python=sys.version,platform=platform.platform(),numpy=np.__version__,pandas=pd.__version__,scipy=scipy.__version__)
(O/'verification_summary.json').write_text(json.dumps(summary,indent=2))
print(json.dumps(summary,indent=2))

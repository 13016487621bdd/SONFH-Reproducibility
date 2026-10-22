import pandas as pd, numpy as np, os, json
from scipy.stats import spearmanr
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
base='/mnt/data/sonfh_e1_closure'
out='/mnt/data/SONFH_E1_Consistency_Closure'
os.makedirs(out,exist_ok=True)
manifest=pd.read_csv(base+'/03_QC/SONFH_sample_manifest_aligned.csv').set_index('GSM')
expr=pd.read_csv(base+'/03_QC/SONFH_QC03_gene_expression_median_log2.csv.gz',index_col=0)
samples=manifest.index.tolist()
X=expr[samples].T
score=pd.read_csv(base+'/04_GeneSet_Scoring/SONFH_GSVA_scores_85x40.csv').set_index('Unnamed: 0')
cat_name='GOBP_HYDROGEN_PEROXIDE_CATABOLIC_PROCESS'
bio_name='GOBP_HYDROGEN_PEROXIDE_BIOSYNTHETIC_PROCESS'
assert cat_name in score.index and bio_name in score.index
cat=score.loc[cat_name,samples].astype(float)
bio=score.loc[bio_name,samples].astype(float)
proxy_genes=['ALAS2','AHSP','GYPA','SLC4A1']
assert all(g in X.columns for g in proxy_genes)
# cohort z scores (sample SD), mean across exactly 4 genes
z=X[proxy_genes].apply(lambda s:(s-s.mean())/s.std(ddof=1))
proxy=z.mean(axis=1).reindex(samples)
# group/covariates
group=manifest.loc[samples,'group']
case=(group=='SONFH').astype(int).values
male=(manifest.loc[samples,'sex']=='Male').astype(int).values

def fit(y, include_proxy=False):
    cols=[np.ones(len(samples)),case,male]
    if include_proxy: cols.append(proxy.values)
    m=sm.OLS(np.asarray(y,dtype=float),np.column_stack(cols)).fit(cov_type='HC3',use_t=True)
    idx=1
    return {'beta':float(m.params[idx]),'SE_HC3':float(m.bse[idx]),'CI95_low':float(m.conf_int()[idx,0]),'CI95_high':float(m.conf_int()[idx,1]),'df_t':int(m.df_resid),'P_value_t':float(m.pvalues[idx])}

# E1 marker table using frozen requested/measured definition; other marker panels per current E1
panels={
'Neutrophil-associated':['FCGR3B','CEACAM8','CSF3R','S100A8','S100A9'],
'Monocyte-associated':['CD14','LYZ','FCGR1A','LST1','S100A12'],
'T-cell-associated':['CD3D','CD3E','CD2','TRAC','IL7R'],
'B-cell-associated':['MS4A1','CD79A','CD79B','CD19','BANK1'],
'Erythroid/globin-associated':['HBB','ALAS2','AHSP','GYPA','SLC4A1'],
'Platelet-associated':['PPBP','PF4','GP9','ITGA2B','TUBB1']}
rows=[]; proxies={}
for name,requested in panels.items():
    used=[g for g in requested if g in X.columns]
    if not used: continue
    pz=X[used].apply(lambda s:(s-s.mean())/s.std(ddof=1))
    p=pz.mean(axis=1).reindex(samples)
    proxies[name]=p
    rho,pval=spearmanr(p,cat)
    control=(group=='GC_exposed_non_SONFH')
    cases=(group=='SONFH')
    rho_ctrl,p_ctrl=spearmanr(p[control],cat[control])
    rho_case,p_case=spearmanr(p[cases],cat[cases])
    rows.append({'marker_panel':name,'genes_requested':';'.join(requested),'genes_measured_and_used':';'.join(used),'n_genes_used':len(used),'control_mean_score':p[control].mean(),'case_mean_score':p[cases].mean(),'case_minus_control_mean_score':p[cases].mean()-p[control].mean(),'Spearman_rho_marker_with_catabolism_all40':rho,'Spearman_P_all40':pval,'Spearman_rho_controls10':rho_ctrl,'Spearman_P_controls10':p_ctrl,'Spearman_rho_cases30':rho_case,'Spearman_P_cases30':p_case,'interpretation':'Unvalidated expression proxy; not a measured cell fraction'})
E1=pd.DataFrame(rows)
E1.to_csv(out+'/Exploratory_Table_E1_unified.csv',index=False)
# individual scores include exact proxy, 5-gene sensitivity, group/sex, and GSVA scores
p5genes=['HBB']+proxy_genes
z5=X[p5genes].apply(lambda s:(s-s.mean())/s.std(ddof=1)); proxy5=z5.mean(axis=1).reindex(samples)
ind=pd.DataFrame({'GSM':samples,'group':group.values,'sex':manifest.loc[samples,'sex'].values,'E1_four_gene_proxy':proxy.values,'E4_five_gene_proxy_including_HBB':proxy5.values,'hydrogen_peroxide_catabolism_GSVA':cat.values,'hydrogen_peroxide_biosynthesis_GSVA':bio.values})
ind['biosynthesis_minus_catabolism_standardized']=((bio-bio.mean())/bio.std(ddof=1)-(cat-cat.mean())/cat.std(ddof=1)).values
ind.to_csv(out+'/E1_proxy_individual_scores.csv',index=False)
# three core outcomes, primary vs exploratory
pair=ind['biosynthesis_minus_catabolism_standardized'].values
outcomes=[('Hydrogen_peroxide_catabolism_GSVA',cat.values),('Hydrogen_peroxide_biosynthesis_GSVA',bio.values),('Hydrogen_peroxide_biosynthesis_minus_catabolism_standardized',pair)]
modelrows=[]
for name,y in outcomes:
    rho,rhop=spearmanr(proxy,y)
    pri=fit(y,False); adj=fit(y,True)
    modelrows.append({'outcome':name,'proxy_definition':'E1 frozen four-gene proxy: ALAS2;AHSP;GYPA;SLC4A1','Spearman_rho_all40':rho,'Spearman_p_all40':rhop,'primary_beta':pri['beta'],'primary_SE_HC3':pri['SE_HC3'],'primary_CI95_low':pri['CI95_low'],'primary_CI95_high':pri['CI95_high'],'primary_P_t':pri['P_value_t'],'proxy_adjusted_beta':adj['beta'],'proxy_adjusted_SE_HC3':adj['SE_HC3'],'proxy_adjusted_CI95_low':adj['CI95_low'],'proxy_adjusted_CI95_high':adj['CI95_high'],'proxy_adjusted_P_t':adj['P_value_t'],'interpretation':'Exploratory adjustment only; proxy is not validated cell proportion and may overadjust disease-related expression'})
E2=pd.DataFrame(modelrows); E2.to_csv(out+'/Exploratory_Table_E2_E1_proxy_adjustment.csv',index=False)
# HBB inclusion sensitivity
rows=[]
for name,y in outcomes:
    rho,rhop=spearmanr(proxy5,y); adj=fit(y,True) # fit currently only uses 4-gene proxy; custom model below
    m=sm.OLS(np.asarray(y,dtype=float),np.column_stack([np.ones(40),case,male,proxy5.values])).fit(cov_type='HC3',use_t=True)
    rows.append({'outcome':name,'proxy_definition':'Sensitivity proxy: HBB;ALAS2;AHSP;GYPA;SLC4A1','n_genes':5,'Spearman_rho_with_proxy':rho,'Spearman_P':rhop,'beta_SONFH_adjusted':float(m.params[1]),'SE_HC3':float(m.bse[1]),'CI95_low':float(m.conf_int()[1,0]),'CI95_high':float(m.conf_int()[1,1]),'P_value_t':float(m.pvalues[1]),'interpretation':'Sensitivity only; HBB was not included in the frozen E1 four-gene proxy'})
pd.DataFrame(rows).to_csv(out+'/Exploratory_Table_E4_HBB_inclusion_sensitivity.csv',index=False)
# all 85 programs, primary and E1-adjusted; compute FDR separately
scores=score[samples]
rows=[]
for gs,row in scores.iterrows():
    y=row.values.astype(float)
    pri=fit(y,False); adj=fit(y,True)
    rows.append({'gene_set':gs,'beta_primary_sex_adjusted':pri['beta'],'P_primary_t':pri['P_value_t'],'beta_exploratory_E1_proxy_adjusted':adj['beta'],'P_exploratory_t':adj['P_value_t']})
E3=pd.DataFrame(rows)
E3['FDR_primary_recomputed']=multipletests(E3.P_primary_t,method='fdr_bh')[1]
E3['FDR_exploratory_proxy_adjusted']=multipletests(E3.P_exploratory_t,method='fdr_bh')[1]
E3.to_csv(out+'/Exploratory_Table_E3_all85_E1_proxy_adjustment.csv',index=False)
# all correlations group-specific for four-gene proxy vs catabolism
corrrows=[]
for label,mask in [('all40',np.ones(40,dtype=bool)),('controls10',(group=='GC_exposed_non_SONFH').values),('SONFH30',(group=='SONFH').values)]:
    rho,p=spearmanr(proxy.values[mask],cat.values[mask]); corrrows.append({'subset':label,'n':int(mask.sum()),'rho':rho,'P':p})
pd.DataFrame(corrrows).to_csv(out+'/E1_catabolism_proxy_correlations.csv',index=False)
# summary and consistency checks
summary={
 'samples':len(samples),'cases':int(case.sum()),'controls':int((case==0).sum()),
 'E1_proxy_genes':proxy_genes,'HBB_present_in_matrix':bool('HBB' in X.columns),'HBB_in_frozen_proxy':False,
 'E1_four_gene_mean_control':float(proxy[group=='GC_exposed_non_SONFH'].mean()),'E1_four_gene_mean_case':float(proxy[group=='SONFH'].mean()),
 'catabolism_proxy_rho':float(spearmanr(proxy,cat).statistic),'catabolism_proxy_p':float(spearmanr(proxy,cat).pvalue),
 'catabolism_primary':fit(cat.values,False),'catabolism_adjusted':fit(cat.values,True),
 'biosynthesis_primary':fit(bio.values,False),'biosynthesis_adjusted':fit(bio.values,True),
 'pair_primary':fit(pair,False),'pair_adjusted':fit(pair,True),
 'n_primary_fdr_lt_005':int((E3.FDR_primary_recomputed<0.05).sum()),'n_adjusted_fdr_lt_005':int((E3.FDR_exploratory_proxy_adjusted<0.05).sum()),
 'adjusted_model_label':'exploratory only; do not replace primary model'
}
with open(out+'/E1_consistency_closure_summary.json','w') as f: json.dump(summary,f,indent=2)
print(json.dumps(summary,indent=2))

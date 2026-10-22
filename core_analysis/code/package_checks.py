"""Assemble English supplementary tables, input hashes, and artifact validation."""
from pathlib import Path
import re, json, hashlib, importlib.metadata
import pandas as pd
import numpy as np
import fitz
from PIL import Image
ROOT=Path(__file__).resolve().parents[1];I=ROOT/'inputs';R=ROOT/'results';F=ROOT/'figures'
evidence_path=R/'execution_evidence.json'
evidence=json.loads(evidence_path.read_text()) if evidence_path.exists() else {}
native_ok=evidence.get('native_scoring_cr2',{}).get('status')=='PASS'
tissue_ok=evidence.get('tissue_upstream',{}).get('status')=='PASS'
labels=pd.read_csv(I/'gene_set_labels.csv')
labels=labels.rename(columns={'Biological domain':'domain','Functional family':'functional_group','Gene set':'gene_set','Biological process label (English)':'process','MSigDB gene count':'frozen_gene_count','GPL15207 mapped gene count':'platform_mapped_gene_count','Platform coverage (%)':'platform_coverage_percent'})
labels=labels[['domain','functional_group','gene_set','process','frozen_gene_count','platform_mapped_gene_count','platform_coverage_percent']]
expr=pd.read_csv(I/'03_QC/SONFH_QC03_gene_expression_median_log2.csv.gz',index_col=0)
sets={}
for line in (I/'frozen_85.gmt').read_text().splitlines():
    v=line.split('\t');sets[v[0]]=v[2:]
labels['measured_genes_in_analysis_matrix']=[sum(g in expr.index for g in sets[k]) for k in labels.gene_set]
labels['measured_gene_members']=[';'.join(g for g in sets[k] if g in expr.index) for k in labels.gene_set]
labels.to_csv(R/'Supplementary_Table_S1_gene_sets.csv',index=False)
def save(t,n):t.to_csv(R/n,index=False)
save(labels.merge(pd.read_csv(R/'layer1_GSVA.csv'),on='gene_set'), 'Supplementary_Table_S2_all85_associations.csv')
f=pd.read_csv(R/'layer2_ARCHIVED_CR2.csv')
def cleanfamily(s):
    s=s.split('.',1)[1];s=re.sub(r'^Family \d+\s*[—-]\s*','',s).split('（')[0]
    return {'胆固醇 / 甾醇':'Cholesterol / sterols','磷脂':'Phospholipids','三酰甘油 / 中性脂质':'Triglycerides / neutral lipids','神经酰胺':'Ceramides','糖脂':'Glycolipids','脂蛋白':'Lipoproteins','脂滴':'Lipid droplets','脂肪酸':'Fatty acids'}.get(s,s)
f['domain_number']=f.family.str.extract(r'Domain (\d)').astype(int);f['functional_group']=f.family.map(cleanfamily)
f=f.drop(columns='family');f['execution_status']=('CR2 refit and BH verified in user-run native R; distribution calculations independently verified in Python' if native_ok else 'Archived CR2 output; native R rerun pending; F distribution and BH verified in Python')
save(f,'Supplementary_Table_S3_family_comparisons.csv')
a=pd.read_csv(R/'layer3_GSVA.csv');b=pd.read_csv(R/'layer3_ssGSEA.csv')
save(a.merge(b,on=['Pair_ID','Set_A','Set_B'],suffixes=('_GSVA','_ssGSEA')),'Supplementary_Table_S4_opposing_contrasts.csv')
r=pd.read_csv(R/'gene_removal_nominal.csv');counts=pd.read_csv(I/'09_KeyGene_Removal_Audit/gene_removal_HC3_results.csv')
r=r.merge(counts[['scenario','n_genes_in_frozen_set','n_genes_measured']],on='scenario').rename(columns={'p':'P_nominal'})
r['execution_status']=('GSVA rescoring verified in user-run native R; regression independently verified in Python' if native_ok else 'Regression independently verified; native GSVA rescoring pending')
save(r,'Supplementary_Table_S5_membership_sensitivity.csv')
save(pd.read_csv(R/'marker_panels_corrected.csv'),'Supplementary_Table_S6_expression_proxies.csv')
p=pd.read_csv(R/'proxy_models.csv').rename(columns={'p':'P_nominal'});p['interpretation']='Exploratory expression proxy; not measured cell fractions'
save(p,'Supplementary_Table_S7_proxy_regression.csv')
save(pd.read_csv(R/'exploratory_all85_proxy_adjusted.csv'),'Supplementary_Table_S8_all85_proxy_adjusted.csv')
g=pd.read_csv(R/'peroxide_member_gene_results.csv');g['biosynthesis_member']=g.gene.isin(sets['GOBP_HYDROGEN_PEROXIDE_BIOSYNTHETIC_PROCESS']);g['catabolism_member']=g.gene.isin(sets['GOBP_HYDROGEN_PEROXIDE_CATABOLIC_PROCESS'])
save(g,'Supplementary_Table_S9_peroxide_genes.csv')
t=pd.read_csv(I/'06_Validation/GSE345158_Pair11_overlap_sensitivity.csv')
if native_ok:
    t['execution_status']=['Corrected tissue scoring and exact P verified in user-run native R and Python','Shared-gene rescoring and exact P verified in user-run native R']
else:
    t['execution_status']=['Exact permutation independently recomputed from corrected scores; GSVA rerun pending','Archived shared-gene rescoring result; native R rerun pending']
if tissue_ok:
    ts=json.loads((R/'tissue_native_R/verification_summary.json').read_text())
    t['archived_Difference']=t.Difference
    idx=t.Version=='SOD2_corrected'
    assert idx.sum()==1
    t.loc[idx,'Difference']=ts['effect']
    t.loc[idx,'Exact_P']=ts['exact_P']
    t.loc[idx,'execution_status']='Counts-to-VST, frozen annotation, GSVA and exact P verified in user-run native R; exact P also verified in Python'
save(t,'Supplementary_Table_S10_tissue_sensitivity.csv')
m=pd.read_csv(I/'03_QC/SONFH_sample_manifest_aligned.csv');rows=[]
for label,level in [('Participants',None),('Female','Female'),('Male','Male')]:
    rows.append(dict(characteristic=label,GC_exposed_without_SONFH=int(((m.group=='GC_exposed_non_SONFH')&(True if level is None else m.sex==level)).sum()),SONFH=int(((m.group=='SONFH')&(True if level is None else m.sex==level)).sum())))
save(pd.DataFrame(rows),'Table1_available_sample_characteristics.csv')
checks=[]
for n in ['Figure1','Figure2','Figure3','Figure4','Figure5','Supplementary_Figure_S1','Supplementary_Figure_S2','Supplementary_Figure_S3']:
    p=F/(n+'.pdf');d=fitz.open(p);expected=3 if n.endswith('S1') else 1
    assert len(d)==expected,(n,len(d));assert p.stat().st_size<10_000_000
    for page in d:
        assert abs(page.rect.width/72*25.4-170)<.05
        assert page.rect.height/72*25.4<200
        for block in page.get_text('dict')['blocks']:
            for line in block.get('lines',[]):
                for s in line['spans']:
                    x0,y0,x1,y1=s['bbox']
                    assert x0>=-.5 and y0>=-.5 and x1<=page.rect.width+.5 and y1<=page.rect.height+.5,(n,s['text'])
    checks.append(dict(file=p.name,pages=len(d),status='PASS'))
    if n.startswith('Figure'):
        t=F/(n+'.tiff');im=Image.open(t);im.load();assert min(im.info['dpi'])>=599
        assert t.stat().st_size<10_000_000
        assert im.info.get('compression')=='tiff_lzw'
        checks.append(dict(file=t.name,pages=1,status='PASS'))
assert len(fitz.open(ROOT/'SONFH_JTM_Figure_Review.pdf'))==10
legends=(ROOT/'documentation/Figure_legends.txt').read_text()
for section in re.split(r'\n\n(?=(?:Figure \d|Supplementary Figure S\d))',legends):
    title,body=section.split('\n',1)
    assert len(title.split())<=18,(title,'title length') # <=15 after removing Figure-number prefix
    assert len(body.split())<=300,(title,len(body.split()))
    checks.append(dict(file=title,pages=None,status=f'{len(body.split())} legend words'))
save(pd.DataFrame(checks),'artifact_validation.csv')
hashes=[]
for p in sorted(I.rglob('*')):
    if p.is_file():hashes.append(dict(path=str(p.relative_to(ROOT)),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
save(pd.DataFrame(hashes),'input_SHA256.csv')
pkgs=['numpy','pandas','scipy','matplotlib','openpyxl','rdata','PyMuPDF','Pillow']
(ROOT/'requirements-tested.txt').write_text('\n'.join(k+'=='+importlib.metadata.version(k) for k in pkgs)+'\n')
print('Validated 5 main figures, 3 supplementary figures, 10 review pages; generated 10 supplementary tables and sample counts.')

"""Regenerate figures from verified tables. Run reproduce.py first."""
from pathlib import Path
import re, textwrap, json
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
from PIL import Image

ROOT=Path(__file__).resolve().parents[1]; R=ROOT/'results'; I=ROOT/'inputs'; F=ROOT/'figures'
F.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.labelsize':8,
 'axes.titlesize':9,'xtick.labelsize':7,'ytick.labelsize':7,'pdf.fonttype':42,
 'ps.fonttype':42,'svg.fonttype':'none','axes.linewidth':.6,'lines.linewidth':1,
 'savefig.facecolor':'white','axes.spines.top':False,'axes.spines.right':False})
BLUE='#2878A5';RED='#C65848';GREY='#88939D';DARK='#243745';TEAL='#258C83';GOLD='#B78D35'
DC={1:'#41749B',2:'#B68A32',3:'#258C83',4:'#8B6BA7',5:'#A66F49',6:'#73828C'}
main=[];supp=[]
def save(fig,name,mainfigure=True):
    def atomic(ext,**kwargs):
        tmp=F/(name+'.building');fig.savefig(tmp,format=ext,**kwargs);tmp.replace(F/(name+'.'+ext))
    atomic('pdf')
    atomic('png',dpi=180)
    if mainfigure:
        atomic('tiff',dpi=600,pil_kwargs={'compression':'tiff_lzw'})
        atomic('svg')
        main.append(name)
    else:supp.append(name)
    plt.close(fig)
def panel(ax,letter,title):
    ax.set_title(title,loc='left',pad=12,fontweight='bold')
    ax.text(-.12,1.06,letter,transform=ax.transAxes,fontsize=12,fontweight='bold',va='bottom')
def clean(ax):
    ax.spines[['left','top','right']].set_visible(False)
    ax.tick_params(axis='y',length=0)
    ax.grid(axis='x',color='#E6EBEF',lw=.5,zorder=0)
def ptext(p):return f'{p:.2g}' if p<.001 else f'{p:.4f}'
lab=pd.read_csv(I/'gene_set_labels.csv').set_index('Gene set')['Biological process label (English)'].to_dict()
def label(k):
    return lab.get(k,k.replace('GOBP_','').replace('HALLMARK_','').replace('_',' ').capitalize()).replace(' (Hallmark)','')
audit=pd.read_csv(I/'05_Statistics/SONFH_Layer1_GSVA_85_audit.csv').set_index('gene_set')
dom={k:int(re.search(r'Domain (\d)',str(audit.loc[k,'Domain'])).group(1)) for k in audit.index}
l1=pd.read_csv(R/'layer1_GSVA.csv');l3=pd.read_csv(R/'layer3_GSVA.csv');fam=pd.read_csv(R/'layer2_ARCHIVED_CR2.csv')
ind=pd.read_csv(R/'participant_scores.csv');case=(ind.group=='SONFH').to_numpy()
gs=pd.read_csv(I/'04_GeneSet_Scoring/SONFH_GSVA_scores_85x40.csv',index_col=0)
cat='GOBP_HYDROGEN_PEROXIDE_CATABOLIC_PROCESS';bio='GOBP_HYDROGEN_PEROXIDE_BIOSYNTHETIC_PROCESS'
pairlabels=['Coagulation - fibrinolysis','Endothelial proliferation - apoptosis','Fatty acid synthesis - beta-oxidation','Triglyceride storage - catabolism','Cholesterol import - efflux','Phospholipid synthesis - catabolism','Ceramide synthesis - catabolism','Glycolipid synthesis - catabolism','Lipoprotein synthesis - catabolism','Lipid droplet formation - disassembly','H2O2 biosynthesis - catabolism']
def fdrkey(color):
    return [Line2D([],[],marker='o',linestyle='none',color=color,markersize=4,label='FDR < 0.05'),Line2D([],[],marker='.',linestyle='none',color=GREY,markersize=5,label='FDR >= 0.05')]

# Figure 1: analytical branches share the same source, with no sequential filter.
fig=plt.figure(figsize=(6.693,5.65));ax=fig.add_axes([0,0,1,1]);ax.set_axis_off();ax.set_xlim(0,1);ax.set_ylim(0,1)
def box(x,y,w,h,title,body,color=DARK,fs=8):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.009,rounding_size=0.012',fc='#F7F9FB',ec=color,lw=.9))
    ax.text(x+w/2,y+h*.73,title,ha='center',va='center',fontweight='bold',fontsize=fs,color=color)
    ax.text(x+w/2,y+h*.34,body,ha='center',va='center',fontsize=fs-0.4,linespacing=1.45,color=DARK)
def arrow(a,b):ax.add_patch(FancyArrowPatch(a,b,arrowstyle='-|>',mutation_scale=9,color='#70808D',lw=.8))
box(.08,.825,.84,.135,'BLOOD TRANSCRIPTOME | GSE123568','Glucocorticoid-exposed: 30 SONFH and 10 without SONFH\nPrimeView microarray | 18,833 genes')
box(.08,.62,.84,.145,'PRESPECIFIED PATHOPHYSIOLOGICAL FRAMEWORK','6 domains | 85 functional gene sets | 28 functional groups\n21 groups eligible for within-group comparison | 11 opposing-process contrasts')
arrow((.5,.825),(.5,.775))
for x,c,title,body in [(.035,BLUE,'INDIVIDUAL GENE SETS','85 sex-adjusted associations\n8 with FDR < 0.05'),(.36,TEAL,'RELATED GENE SETS','21 within-group comparisons\n0 with FDR < 0.05'),(.685,RED,'OPPOSING PROCESSES','11 prespecified contrasts\n1 with FDR < 0.05')]:
    box(x,.385,.28,.155,title,body,c,7.7);arrow((.5,.62),(x+.14,.548))
ax.text(.5,.635,'Three parallel hypothesis families; FDR correction within each',ha='center',fontsize=6.8,color='#526474')
box(.055,.12,.43,.18,'SENSITIVITY ASSESSMENT','Alternative scoring | gene membership\nIndividual-sample influence\nFrozen four-gene erythroid expression proxy',BLUE,7.8)
box(.535,.12,.41,.18,'EXPLORATORY TISSUE ASSESSMENT','GSE345158: 5 SONFH and 4 fracture controls\nExact permutation | shared-gene sensitivity\nSOD2-corrected annotation',TEAL,7.8)
for center in [.175,.5,.825]:ax.plot([center,center],[.385,.345],color='#70808D',lw=.8)
ax.plot([.175,.825],[.345,.345],color='#70808D',lw=.8)
arrow((.27,.345),(.27,.31));arrow((.74,.345),(.74,.31))
ax.text(.5,.045,'Cross-sectional associations with established SONFH',ha='center',fontsize=8,fontweight='bold',color=DARK)
save(fig,'Figure1')

# Figure 2: selected associations and observed per-participant expression patterns.
sig=l1[l1.FDR<.05].copy();sig['domain']=sig.gene_set.map(dom)
sig=sig.sort_values(['domain','beta']).reset_index(drop=True)
fig=plt.figure(figsize=(6.693,5.7));gr=fig.add_gridspec(2,1,height_ratios=[1,1.1],hspace=.53,left=.37,right=.91,top=.94,bottom=.11)
ax=fig.add_subplot(gr[0]);yy=np.arange(len(sig))
for j,row in sig.iterrows():
    ax.errorbar(row.beta,j,xerr=[[row.beta-row.CI_lower],[row.CI_upper-row.beta]],fmt='o',color=BLUE,ms=4,capsize=2,zorder=3)
    ax.text(1.06,j,ptext(row.FDR),transform=ax.get_yaxis_transform(),va='center',fontsize=7,clip_on=False)
ax.text(1.06,1.07,'FDR',transform=ax.transAxes,fontsize=7,fontweight='bold')
ax.set_yticks(yy,[label(k) for k in sig.gene_set]);ax.invert_yaxis();ax.axvline(0,color='#687A89',ls='--',lw=.8)
ax.set_xlim(-.7,.8);ax.set_xlabel('Sex-adjusted difference in GSVA score (95% CI)');clean(ax);panel(ax,'A','Eight associated gene sets')
ax=fig.add_subplot(gr[1]);data=gs.loc[sig.gene_set].to_numpy();z=(data-data.mean(axis=1,keepdims=True))/data.std(axis=1,ddof=1,keepdims=True)
im=ax.imshow(z,aspect='auto',cmap='RdBu_r',vmin=-2.5,vmax=2.5,interpolation='nearest')
ax.set_yticks(yy,[label(k) for k in sig.gene_set]);ax.set_xticks([4.5,24.5],['GC-exposed\nwithout SONFH (n = 10)','SONFH (n = 30)']);ax.tick_params(length=0)
ax.axvline(9.5,color='white',lw=1.8);panel(ax,'B','Individual expression patterns')
cbax=fig.add_axes([.94,.15,.014,.18]);fig.colorbar(im,cax=cbax,ticks=[-2,0,2]);cbax.set_title('Z',fontsize=7,pad=5)
save(fig,'Figure2')

# Figure 3: show multiplicity-adjusted evidence rather than nominal P alone.
def familylabel(k):
    s=k.split('.',1)[1];s=re.sub(r'^Family \d+\s*[—-]\s*','',s);s=s.split('（')[0].strip()
    return {'胆固醇 / 甾醇':'Cholesterol / sterols','磷脂':'Phospholipids','三酰甘油 / 中性脂质':'Triglycerides / neutral lipids','神经酰胺':'Ceramides','糖脂':'Glycolipids','脂蛋白':'Lipoproteins','脂滴':'Lipid droplets','脂肪酸':'Fatty acids'}.get(s,s)
ff=fam.sort_values('FDR').reset_index(drop=True)
fig,ax=plt.subplots(figsize=(6.693,5.2));fig.subplots_adjust(left=.46,right=.87,top=.91,bottom=.13)
for j,r in ff.iterrows():
    d=int(re.search(r'Domain (\d)',r.family).group(1));ax.scatter(-np.log10(r.FDR),j,c=GREY,s=22,zorder=3)
    ax.text(1.04,j,f'{r.FDR:.3f}',transform=ax.get_yaxis_transform(),va='center',fontsize=7)
ax.text(1.04,1.04,'FDR',transform=ax.transAxes,fontweight='bold',fontsize=7)
ax.set_yticks(range(21),[familylabel(x) for x in ff.family]);ax.invert_yaxis()
ax.axvline(-np.log10(.05),ls='--',color=RED,lw=1);ax.text(-np.log10(.05)-.03,20.9,'FDR = 0.05',ha='right',fontsize=7,color=RED)
ax.set_xlim(0,1.4);ax.set_xlabel('-log10(FDR-adjusted P)');clean(ax)
ax.text(0,1.05,'21 functional groups',transform=ax.transAxes,fontsize=8,fontweight='bold')
save(fig,'Figure3')

# Figure 4: full prespecified contrasts; components and subject-level observations.
fig=plt.figure(figsize=(6.693,6.35));gr=fig.add_gridspec(2,2,height_ratios=[1.35,1],wspace=.75,hspace=.64,left=.14,right=.93,top=.95,bottom=.1)
ax=fig.add_subplot(gr[0,:]);ax.set_position([.43,.58,.43,.36])
for j,r in l3.iterrows():
    c=RED if r.FDR<.05 else GREY
    ax.errorbar(r.beta,j,xerr=[[r.beta-r.CI_lower],[r.CI_upper-r.beta]],fmt='o',ms=4,color=c,capsize=2,zorder=3)
    ax.text(1.08,j,f'{r.FDR:.3f}',transform=ax.get_yaxis_transform(),va='center',fontsize=7)
ax.text(1.08,1.05,'FDR',transform=ax.transAxes,fontweight='bold',fontsize=7)
ax.set_yticks(range(11),pairlabels);ax.invert_yaxis();ax.axvline(0,ls='--',color=GREY,lw=.8)
ax.set_xlim(-2.1,2.8);ax.set_xlabel('Sex-adjusted contrast difference (95% CI)');clean(ax);panel(ax,'A','Eleven opposing-process contrasts')
fig.legend(handles=fdrkey(RED),frameon=False,fontsize=6.5,loc='upper center',bbox_to_anchor=(.65,.515),ncol=2,handletextpad=.3)
ax=fig.add_subplot(gr[1,0]);comp=pd.read_csv(R/'components_standardized.csv')
for j,r in comp.iterrows():ax.errorbar(r.beta,j,xerr=[[r.beta-r.CI_lower],[r.CI_upper-r.beta]],fmt='o',ms=5,color=BLUE,capsize=3)
r=l3.iloc[-1];ax.errorbar(r.beta,2,xerr=[[r.beta-r.CI_lower],[r.CI_upper-r.beta]],fmt='D',ms=5,color=RED,capsize=3)
ax.set_yticks([0,1,2],['Biosynthesis','Catabolism','Difference']);ax.invert_yaxis();ax.axvline(0,color=GREY,ls='--',lw=.8);clean(ax)
ax.set_xlabel('Standardized score coefficient\n(95% CI)');panel(ax,'B','Peroxide components')
ax=fig.add_subplot(gr[1,1]);
for j,col in enumerate(['biosynthesis_z','catabolism_z']):
    for mask,offset,c in [(~case,-.16,BLUE),(case,.16,RED)]:
        y=ind.loc[mask,col].to_numpy();x=j+offset+np.random.default_rng(30+j+int(mask.sum())).uniform(-.08,.08,len(y))
        ax.scatter(x,y,s=12,c=c,alpha=.7,edgecolors='none')
        ax.plot([j+offset-.09,j+offset+.09],[y.mean()]*2,c=c,lw=1.7)
ax.set_xticks([0,1],['Biosynthesis','Catabolism']);ax.set_ylabel('Standardized GSVA score');ax.axhline(0,color='#D8E0E5',lw=.7)
ax.scatter([],[],s=15,color=BLUE,label='Without SONFH');ax.scatter([],[],s=15,color=RED,label='SONFH')
handles,labels=ax.get_legend_handles_labels()
fig.legend(handles,labels,frameon=False,fontsize=7,loc='lower center',bbox_to_anchor=(.72,.015),ncol=2,handletextpad=.25)
panel(ax,'C','Individual scores')
save(fig,'Figure4')

# Figure 5: interpretation and cross-tissue evidence, all individual tissue samples.
fig=plt.figure(figsize=(6.693,5.8));gr=fig.add_gridspec(2,2,wspace=.63,hspace=.73,left=.12,right=.95,top=.93,bottom=.11)
ax=fig.add_subplot(gr[0,0]);
for mask,c,l in [(~case,BLUE,'Without SONFH'),(case,RED,'SONFH')]:ax.scatter(ind.loc[mask,'E1_four_gene_proxy'],gs.loc[cat].to_numpy()[mask],s=19,c=c,alpha=.8,label=l,edgecolors='none')
rho=stats.spearmanr(ind.E1_four_gene_proxy,gs.loc[cat]).statistic
ax.text(.03,.96,f'Spearman rho = {rho:.3f}',transform=ax.transAxes,va='top',fontsize=7)
ax.set_xlabel('Four-gene erythroid expression proxy');ax.set_ylabel('Catabolism GSVA score');ax.legend(frameon=False,fontsize=6.3,loc='lower right');panel(ax,'A','Erythroid-associated expression')
ax=fig.add_subplot(gr[0,1]);pm=pd.read_csv(R/'proxy_models.csv');p=pm[(pm.outcome=='Catabolism')&(pm.model!='Five-gene proxy sensitivity')]
for j,(_,r) in enumerate(p.iterrows()):
    ax.errorbar(r.beta,j,xerr=[[r.beta-r.CI_lower],[r.CI_upper-r.beta]],fmt='o',color=BLUE if j==0 else TEAL,ms=5,capsize=3)
    ax.text(.02,.92-j*.43,f'P = {ptext(r.p)}',transform=ax.transAxes,fontsize=7)
ax.set_yticks([0,1],['Sex adjusted','Proxy adjusted']);ax.set_ylim(1.6,-.6);ax.axvline(0,ls='--',color=GREY,lw=.8);clean(ax)
ax.set_xlabel('Catabolism coefficient (95% CI)');panel(ax,'B','Exploratory proxy adjustment')
ax=fig.add_subplot(gr[1,0]);t=pd.read_csv(R/'tissue_corrected_scores.csv');tc=t.group=='SONFH'
for j,mask,c in [(0,~tc,BLUE),(1,tc,RED)]:
    v=t.loc[mask,'Pair11'].to_numpy();ax.scatter(j+np.linspace(-.08,.08,len(v)),v,s=27,c=c,edgecolors='white',lw=.4,zorder=3);ax.plot([j-.14,j+.14],[v.mean()]*2,c=c,lw=1.7)
ax.set_xticks([0,1],['Fracture controls\n(n = 4)','SONFH\n(n = 5)']);ax.set_xlim(-.45,1.45);ax.set_ylabel('Biosynthesis - catabolism\n(standardized score contrast)');ax.axhline(0,c='#D8E0E5',lw=.7)
ax.text(.04,.97,'Exact P = 0.0397',transform=ax.transAxes,va='top',fontsize=7);panel(ax,'C','Femoral head tissue')
ax=fig.add_subplot(gr[1,1]);ts=pd.read_csv(I/'06_Validation/GSE345158_Pair11_overlap_sensitivity.csv')
for j,r in ts.iterrows():
    ax.scatter(r.Difference,j,c=TEAL if j==0 else GREY,s=32)
    ax.text(.04,.90-j*.45,f'Exact P = {r.Exact_P:.4f}',transform=ax.transAxes,fontsize=7)
ax.set_yticks([0,1],['Corrected\nannotation','DUOX1/2\nremoved']);ax.set_ylim(1.6,-.6);ax.set_xlim(0,1.65);clean(ax)
ax.set_xlabel('SONFH - fracture-control\nmean contrast difference');panel(ax,'D','Shared-gene sensitivity')
save(fig,'Figure5')

# Supplementary figure S1: all 85 results across three legible pages.
domainnames={1:'Vascular / microcirculation',2:'Lipid metabolism',3:'Hypoxia / oxidative stress / redox',4:'Inflammation / immunity',5:'Bone remodeling / repair',6:'Cellular injury / death'}
with PdfPages(F/'Supplementary_Figure_S1.building') as pdf:
    for part,domains in enumerate([(1,2),(3,4),(5,6)],1):
        d=l1[l1.gene_set.map(dom).isin(domains)].copy();d['domain']=d.gene_set.map(dom);d=d.sort_values(['domain','beta']).reset_index(drop=True)
        fig,ax=plt.subplots(figsize=(6.693,7.0 if len(d)>20 else 3.5));fig.subplots_adjust(left=.57,right=.94,top=.94,bottom=.12)
        for j,r in d.iterrows():ax.errorbar(r.beta,j,xerr=[[r.beta-r.CI_lower],[r.CI_upper-r.beta]],fmt='o' if r.FDR<.05 else '.',color=BLUE if r.FDR<.05 else GREY,ms=4,capsize=1.5,lw=.7)
        ax.set_yticks(range(len(d)),[textwrap.fill(label(k),43) for k in d.gene_set],fontsize=6.6);ax.invert_yaxis();ax.axvline(0,color=GREY,lw=.7,ls='--');clean(ax)
        ax.set_xlabel('Sex-adjusted GSVA coefficient (95% CI)');ax.set_title(' / '.join(domainnames[k] for k in domains),loc='right',fontsize=8)
        fig.legend(handles=fdrkey(BLUE),frameon=False,fontsize=6.5,loc='lower center',bbox_to_anchor=(.72,.018),ncol=2,handletextpad=.3)
        pdf.savefig(fig);fig.savefig(F/f'Supplementary_Figure_S1_page{part}.png',dpi=180);plt.close(fig)
(F/'Supplementary_Figure_S1.building').replace(F/'Supplementary_Figure_S1.pdf')
supp.append('Supplementary_Figure_S1')

# S2: scoring, targeted gene removal, and fixed-score leave-one-out inference.
fig=plt.figure(figsize=(6.693,6.3));gr=fig.add_gridspec(3,1,height_ratios=[.65,1.1,1],hspace=.8,left=.38,right=.93,top=.94,bottom=.09)
ax=fig.add_subplot(gr[0]);
for j,m in enumerate(['GSVA','ssGSEA']):
    r=pd.read_csv(R/f'layer3_{m}.csv').iloc[-1];ax.errorbar(r.beta,j,xerr=[[r.beta-r.CI_lower],[r.CI_upper-r.beta]],fmt='o',color=RED if j==0 else GREY,capsize=3,ms=5)
    ax.text(.98,j,f'FDR = {r.FDR:.4f}',transform=ax.get_yaxis_transform(),ha='right',va='bottom',fontsize=7)
ax.set_yticks([0,1],['GSVA','ssGSEA']);ax.set_ylim(1.6,-.6);ax.axvline(0,color=GREY,ls='--',lw=.7);clean(ax);ax.set_xlabel('Peroxide contrast coefficient (95% CI)');panel(ax,'A','Alternative scoring')
ax=fig.add_subplot(gr[1]);rr=pd.read_csv(R/'gene_removal_nominal.csv')
rlabels=['Original membership','SNCA removed','PRDX2 removed','SNCA + PRDX2 removed','MPO + EPX removed','All four genes removed']
for j,r in rr.iterrows():ax.errorbar(r.beta,j,xerr=[[r.beta-r.CI_lower],[r.CI_upper-r.beta]],fmt='o',color=TEAL,capsize=2,ms=4)
ax.set_yticks(range(6),rlabels);ax.invert_yaxis();ax.axvline(0,color=GREY,ls='--',lw=.7);clean(ax);ax.set_xlabel('Catabolism GSVA coefficient (95% CI)');panel(ax,'B','Targeted membership sensitivity')
ax=fig.add_subplot(gr[2]);lo=pd.read_csv(R/'pair11_LOO.csv')
ax.scatter(range(1,41),lo.Pair11_FDR,c=np.where(lo.Pair11_FDR<.05,RED,GREY),s=13)
ax.axhline(.05,c=RED,ls='--',lw=.8,label='FDR = 0.05');ax.axhline(l3.FDR.iloc[-1],c=BLUE,lw=.7,label='Full cohort')
ax.legend(frameon=False,fontsize=6.5,loc='upper right')
ax.set_xlabel('Removed participant (frozen sample order)');ax.set_ylabel('Pair 11 FDR');panel(ax,'C','Leave-one-out inference')
save(fig,'Supplementary_Figure_S2',False)

# S3: gene-level results with genome-wide correction, not a second selection step.
ge=pd.read_csv(R/'peroxide_member_gene_results.csv').set_index('gene');sets={}
for line in (I/'frozen_85.gmt').read_text().splitlines():
    a=line.split('\t');sets[a[0]]=a[2:]
fig,axes=plt.subplots(1,2,figsize=(6.693,5.8));fig.subplots_adjust(left=.13,right=.96,top=.91,bottom=.12,wspace=.45)
for ax,key,title,letter in zip(axes,[bio,cat],['Biosynthesis members','Catabolism members'],['A','B']):
    dd=ge.loc[[x for x in sets[key] if x in ge.index]].sort_values('beta')
    for j,(_,r) in enumerate(dd.iterrows()):ax.errorbar(r.beta,j,xerr=[[r.beta-r.CI_lower],[r.CI_upper-r.beta]],fmt='o' if r.FDR<.05 else '.',color=RED if r.FDR<.05 else GREY,ms=4,capsize=2)
    ax.set_yticks(range(len(dd)),dd.index);ax.invert_yaxis();ax.axvline(0,color=GREY,ls='--',lw=.7);clean(ax);ax.set_xlabel('Log2 expression coefficient\n(95% CI)');panel(ax,letter,title)
fig.legend(handles=fdrkey(RED),frameon=False,fontsize=7,loc='lower center',bbox_to_anchor=(.55,.015),ncol=2,handletextpad=.3)
save(fig,'Supplementary_Figure_S3',False)

# One visual review PDF; submission still uses one single-page file per main figure.
import fitz
review=fitz.open()
for n in main+supp:
    doc=fitz.open(F/(n+'.pdf'));review.insert_pdf(doc)
review.save(ROOT/'SONFH_JTM_Figure_Review.building');review.close()
(ROOT/'SONFH_JTM_Figure_Review.building').replace(ROOT/'SONFH_JTM_Figure_Review.pdf')
(F/'figure_inventory.json').write_text(json.dumps({'main':main,'supplementary':supp},indent=2))
print('Created',len(main),'main figures and',len(supp),'supplementary figures.')

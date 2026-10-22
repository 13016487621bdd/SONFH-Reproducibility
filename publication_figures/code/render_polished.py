"""Layout-only revision using verified archived values. Does not refit models.
Figure 6 preserves vector marks from the user's native-R PDF, without digitizing data.
"""
from pathlib import Path
import json,shutil
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch,FancyArrowPatch,Patch
from matplotlib.colors import LinearSegmentedColormap,ListedColormap
import fitz
from PIL import Image

ROOT=Path(__file__).resolve().parents[1];D=ROOT/'plot_data';O=ROOT/'figures';O.mkdir(exist_ok=True)
ink='#263746';muted='#65727C';orange='#B5632D';control='#368A88';case='#8368A7';light='#E0E6EA'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8.5,'axes.linewidth':.65,
 'axes.edgecolor':'#7B8992','text.color':ink,'axes.labelcolor':ink,'xtick.color':ink,
 'ytick.color':ink,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none'})
def save(fig,num):
 for ext in ['pdf','svg','png','tiff']:
  kw={'dpi':300} if ext=='png' else ({'dpi':600,'pil_kwargs':{'compression':'tiff_lzw'}} if ext=='tiff' else {})
  fig.savefig(O/f'Figure{num}.{ext}',facecolor='white',**kw)
 plt.close(fig)
def clean(ax):
 ax.spines[['top','right','left']].set_visible(False);ax.set_yticks([]);ax.tick_params(axis='x',labelsize=8,length=3)
def heading(fig,x,y,letter,title):
 fig.text(x,y,letter,fontsize=11,fontweight='bold',va='top')
 fig.text(x+.034,y,title,fontsize=9.5,fontweight='bold',va='top')
def ci(ax,r,y,color=ink,marker='o'):
 ax.errorbar(r.beta,y,xerr=[[r.beta-r.CI_lower],[r.CI_upper-r.beta]],fmt=marker,
   color=color,elinewidth=1.1,capsize=2,markersize=3.8)
def pformat(v):return f'{v:.2e}' if v<.0001 else f'{v:.4f}'

# Figure 1: concise, parallel structure; source details remain in the legend.
f,a=plt.subplots(figsize=(7.2,6.2));f.subplots_adjust(0,0,1,1);a.set(xlim=(0,7.2),ylim=(0,6.2));a.axis('off')
def text(x,y,s,size=8.5,bold=False,ha='left',color=ink):
 a.text(x,y,s,fontsize=size,fontweight='bold' if bold else 'normal',ha=ha,va='top',linespacing=1.35,color=color)
def box(x,y,w,h,title,body,fill='#EAF1F5',accent='#315D77',style='-'):
 a.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=.015,rounding_size=.04',facecolor=fill,edgecolor='#81939D',lw=.7,linestyle=style))
 text(x+.13,y+h-.13,title,8.8,True,color=accent);text(x+.13,y+h-.43,body,8.3)
def arrow(x1,y1,x2,y2):a.add_patch(FancyArrowPatch((x1,y1),(x2,y2),arrowstyle='-|>',mutation_scale=9,lw=.8,color='#9CABB4'))
text(.16,6.08,'A  Study inputs and functional definitions',10,True)
box(.16,4.55,3.31,1.13,'Blood discovery | GSE123568','30 SONFH + 10 GC-exposed non-SONFH\nDeposited RMA + GPL15207 annotation\nProbe-to-gene median → 18,833 genes × 40')
box(3.68,4.55,3.36,1.13,'Literature-guided functional panel','6 predefined domains; 85 frozen gene sets\n28 functional groups; 21 with ≥2 sets\n11 predefined opposing-process pairs')
arrow(1.82,4.53,1.82,4.33);arrow(5.36,4.53,5.36,4.33)
a.plot([1.82,5.36],[4.33,4.33],color='#9CABB4',lw=.8);arrow(3.6,4.33,3.6,4.12)
text(3.6,4.08,'GSVA scoring across all 40 blood samples',9.5,True,ha='center')
text(.16,3.73,'B  Three prespecified parallel questions',10,True)
a.plot([1.26,5.95],[3.40,3.40],color='#9CABB4',lw=.8)
for x in [1.26,3.60,5.95]:arrow(x,3.40,x,3.19)
box(.16,1.60,2.20,1.55,'Layer 1 | Association','Which functions are\nassociated with SONFH?\n\n85 sets; sex-adjusted HC3\nBH across 85 tests')
box(2.51,1.60,2.20,1.55,'Layer 2 | Effect differences','Do related components have\ndistinguishable effects?\n\n21 groups; standardized scores\nSex-adjusted CR2; BH across 21',fill='#EAF3F0',accent='#286C60',style='--')
box(4.86,1.60,2.18,1.55,'Layer 3 | Process contrasts','Are opposing-process\ndifferences associated?\n\n11 contrasts: z(A) − z(B)\nSex-adjusted HC3; BH across 11',fill='#F8F0E5',accent='#8A622E',style=':')
text(3.6,1.44,'Distinct questions; no sequential significance filter.',9,True,ha='center')
text(.16,1.09,'C  Sensitivity analyses and exploratory follow-up',10,True)
box(.16,.10,3.31,.78,'Blood-based characterization','Alternative scoring; membership; leave-one-out\nErythroid expression proxy; component analysis')
box(3.68,.10,3.36,.78,'Exploratory tissue assessment','GSE345158: 5 SONFH + 4 fracture controls\nExact permutation; shared-gene sensitivity')
save(f,1)

# Figure 2: readable physical font size, shorter label column, aligned forest/heatmap.
res=pd.read_csv(D/'layer1_all85.csv').set_index('gene_set');scores=pd.read_csv(D/'GSVA_scores.csv',index_col=0)
meta=pd.read_csv(D/'blood_manifest.csv').set_index('GSM')
names={'GOBP_HYDROGEN_PEROXIDE_CATABOLIC_PROCESS':'Hydrogen peroxide catabolism',
 'GOBP_GLUTATHIONE_METABOLIC_PROCESS':'Glutathione metabolism','GOBP_GLUTATHIONE_CATABOLIC_PROCESS':'Glutathione catabolism',
 'GOBP_SUPEROXIDE_ANION_GENERATION':'Superoxide anion generation','GOBP_T_CELL_ACTIVATION':'T-cell activation',
 'GOBP_T_CELL_RECEPTOR_SIGNALING_PATHWAY':'T-cell receptor signaling','GOBP_B_CELL_RECEPTOR_SIGNALING_PATHWAY':'B-cell receptor signaling',
 'GOBP_GLYCOLIPID_TRANSPORT':'Glycolipid transport'}
order=list(names);assert set(res.index[res.FDR<.05])==set(order)
sample_order=list(meta.index[meta.group=='GC_exposed_non_SONFH'])+list(meta.index[meta.group=='SONFH'])
s=scores.loc[order,sample_order];z=s.sub(s.mean(axis=1),axis=0).div(s.std(axis=1,ddof=1),axis=0)
f=plt.figure(figsize=(7.2,4.35))
la=f.add_axes([.018,.22,.309,.56]);ax=f.add_axes([.338,.22,.185,.56]);q=f.add_axes([.533,.22,.077,.56]);h=f.add_axes([.643,.22,.34,.56])
for aa in [la,ax,q,h]:aa.set_ylim(7.5,-.5)
for aa in [la,q]:aa.axis('off')
for i,k in enumerate(order):
 la.text(0,i,names[k],va='center',fontsize=8.5);ci(ax,res.loc[k],i);q.text(.5,i,pformat(res.loc[k,'FDR']),ha='center',va='center',fontsize=8.2)
for aa in [la,ax,q]:
 for yy in [3.5,6.5]:aa.axhline(yy,color=light,lw=.6,zorder=0)
ax.axvline(0,color='#9CA7AE',ls='--',lw=.7);ax.set_xlim(-.95,.95);ax.set_xticks([-.8,0,.8]);clean(ax)
ax.set_xlabel('GSVA score difference\nSONFH − exposed control',fontsize=8.3,labelpad=7)
cmap=LinearSegmentedColormap.from_list('blue_orange',['#2166AC','#F7F7F7','#B35806'])
im=h.imshow(z,aspect='auto',interpolation='nearest',vmin=-2,vmax=2,cmap=cmap);h.axis('off');h.axvline(9.5,color='white',lw=1)
for yy in [3.5,6.5]:h.axhline(yy,color='white',lw=.6)
bar=f.add_axes([.643,.797,.34,.023]);bar.imshow([[0]*10+[1]*30],aspect='auto',interpolation='nearest',cmap=ListedColormap([control,case]));bar.axis('off')
heading(f,.018,.975,'A','Sex-adjusted associations');heading(f,.643,.975,'B','Individual profiles')
f.text(.018,.89,'8 gene sets with FDR < 0.05 across 85 tests',fontsize=8.5,color=muted)
f.text(.346,.823,'Effect and 95% CI',fontsize=8.5);f.text(.572,.823,'FDR',ha='center',fontsize=8.5)
f.text(.643,.90,'Exposed non-SONFH (n = 10)',color=control,fontsize=8.1);f.text(.643,.852,'SONFH (n = 30)',color=case,fontsize=8.1)
cax=f.add_axes([.677,.115,.272,.025]);cb=f.colorbar(im,cax=cax,orientation='horizontal',ticks=[-2,0,2]);cb.outline.set_visible(False);cb.ax.tick_params(labelsize=8,length=2)
cb.set_label('Gene-set score z',fontsize=8.5,labelpad=3)
save(f,2);z.to_csv(D/'Figure2_display_z.csv')

# Figure 3: physical 8.2–8.5 pt labels at 7.2 inches, full prespecified families retained.
fam=pd.read_csv(D/'family_comparisons.csv').sort_values(['domain_number','functional_group']).reset_index(drop=True)
pairs=pd.read_csv(D/'layer3_GSVA.csv');assert len(fam)==21 and len(pairs)==11
short={'Angiogenesis / endothelial remodeling':'Angiogenesis / endothelium','Osteoclast differentiation / bone resorption':'Osteoclasts / bone resorption'}
plabels=['Coagulation − fibrinolysis','Endothelial proliferation\n− apoptosis','Fatty acid biosynthesis\n− beta-oxidation','Triglyceride storage\n− catabolism','Cholesterol import\n− efflux','Phospholipid biosynthesis\n− catabolism','Ceramide biosynthesis\n− catabolism','Glycolipid biosynthesis\n− catabolism','Lipoprotein biosynthesis\n− catabolism','Lipid-droplet formation\n− disassembly','Hydrogen peroxide\nbiosynthesis − catabolism']
f=plt.figure(figsize=(7.2,5.75));la=f.add_axes([.018,.15,.292,.68]);ax=f.add_axes([.321,.15,.108,.68]);qa=f.add_axes([.438,.15,.066,.68]);lb=f.add_axes([.53,.15,.218,.68]);bx=f.add_axes([.76,.15,.15,.68]);qb=f.add_axes([.922,.15,.073,.68])
for aa in [la,ax,qa]:aa.set_ylim(20.5,-.5)
for aa in [lb,bx,qb]:aa.set_ylim(10.5,-.5)
for aa in [la,qa,lb,qb]:aa.axis('off')
for i,r in fam.iterrows():
 la.text(0,i,short.get(r.functional_group,r.functional_group),fontsize=8.1,va='center');ax.plot(r.FDR,i,'o',color=ink,ms=3);qa.text(.5,i,f'{r.FDR:.4f}',fontsize=8.1,ha='center',va='center')
 if i<20 and fam.iloc[i+1].domain_number!=r.domain_number:
  for aa in [la,ax,qa]:aa.axhline(i+.5,color=light,lw=.55)
ax.set_xlim(0,1.03);ax.axvline(.05,color=orange,ls=':',lw=1.1);ax.set_xticks([0,.5,1]);clean(ax);ax.set_xlabel('FDR',fontsize=8.5,labelpad=6)
for i,r in pairs.iterrows():
 sig=r.FDR<.05;co=orange if sig else '#75838D';lb.text(0,i,plabels[i],va='center',fontsize=8.3,fontweight='bold' if sig else 'normal',color=co if sig else ink,linespacing=1.05)
 ci(bx,r,i,co,'D' if sig else 'o');qb.text(.5,i,f'{r.FDR:.4f}',ha='center',va='center',fontsize=8.1,color=co,fontweight='bold' if sig else 'normal')
bx.set_xlim(-2.05,2.9);bx.axvline(0,color='#9CA7AE',ls='--',lw=.7);bx.set_xticks([-2,0,2]);clean(bx);bx.set_xlabel('Disease-associated effect\nin z(A) − z(B)',fontsize=8.2,labelpad=6)
heading(f,.018,.975,'A','Related functional groups');heading(f,.53,.975,'B','Opposing-process contrasts')
f.text(.018,.915,'0 of 21 met FDR < 0.05',fontsize=8.5,color=muted);f.text(.53,.915,'1 of 11 met FDR < 0.05',fontsize=8.5,color=muted)
f.text(.33,.85,'FDR',fontsize=8.5);f.text(.471,.85,'FDR',fontsize=8.5,ha='center');f.text(.76,.85,'Effect and 95% CI',fontsize=8.3);f.text(.958,.85,'FDR',fontsize=8.5,ha='center')
f.text(.018,.025,'Separate prespecified correction families: 21 functional groups and 11 contrasts.',fontsize=8.5,color=muted)
save(f,3)

# Figure 4: compact, shared-scale distributions; forest rows directly aligned with numeric columns.
v=pd.read_csv(D/'Figure4_participant_values.csv',index_col=0);comp=pd.read_csv(D/'components.csv');r=pd.concat([comp,pairs.iloc[[-1]][comp.columns.intersection(pairs.columns)]],ignore_index=True)
compcols=['beta','CI_lower','CI_upper','p'];rr=pd.concat([comp[compcols],pairs.iloc[[-1]][compcols]],ignore_index=True)
f=plt.figure(figsize=(7.2,5.35));titles=['Hydrogen peroxide\nbiosynthesis','Hydrogen peroxide\ncatabolism','Biosynthesis − catabolism'];positions=[.085,.405,.725]
caseidx=(meta.loc[v.columns,'group']=='SONFH').to_numpy();yl=(np.floor(v.to_numpy().min()-.2),np.ceil(v.to_numpy().max()+.2))
for j,x in enumerate(positions):
 aa=f.add_axes([x,.47,.245,.35]);groups=[v.iloc[j,~caseidx].to_numpy(),v.iloc[j,caseidx].to_numpy()]
 # R boxplot.stats uses Tukey hinges (not NumPy percentile quartiles for n=10/30).
 def r_box(g):
  g=np.sort(g);n=len(g);n4=np.floor((n+3)/2)/2
  hinge=lambda q:(g[int(np.floor(q))-1]+g[int(np.ceil(q))-1])/2
  q1,q3=hinge(n4),hinge(n+1-n4);iqr=q3-q1
  inside=g[(g>=q1-1.5*iqr)&(g<=q3+1.5*iqr)]
  return dict(q1=q1,med=float(np.median(g)),q3=q3,whislo=float(inside.min()),whishi=float(inside.max()),fliers=[])
 b=aa.bxp([r_box(g) for g in groups],positions=[1,2],widths=.48,showfliers=False,patch_artist=True,medianprops={'color':ink,'linewidth':1},whiskerprops={'linewidth':.7},capprops={'linewidth':.7})
 for patch,col in zip(b['boxes'],[control,case]):patch.set(facecolor=col+'20',edgecolor=col,linewidth=.8)
 for k,(g,col) in enumerate(zip(groups,[control,case])):
  offs=((np.arange(len(g))%7)-3)*.035;aa.scatter(k+1+offs,g,s=9,c=col,edgecolors='white',linewidths=.25,zorder=3)
 aa.set_ylim(yl);aa.set_xlim(.55,2.45);aa.axhline(0,color=light,lw=.7,ls=':');aa.set_xticks([1,2],['Exposed\ncontrol','SONFH'],fontsize=8.3);aa.set_title(titles[j],fontsize=8.8,fontweight='bold',pad=8)
 aa.spines[['top','right','bottom']].set_visible(False);aa.tick_params(axis='x',length=0);aa.tick_params(axis='y',labelsize=8)
 if j==0:aa.set_ylabel('Component z or contrast',fontsize=8.5)
 else:aa.set_yticklabels([]);aa.spines['left'].set_visible(False);aa.tick_params(axis='y',length=0)
heading(f,.018,.975,'A','Individual score distributions')
heading(f,.018,.365,'B','Sex-adjusted associations')
la=f.add_axes([.018,.105,.335,.18]);ax=f.add_axes([.37,.105,.31,.18]);nb=f.add_axes([.70,.105,.113,.18]);pv=f.add_axes([.84,.105,.145,.18])
for aa in [la,ax,nb,pv]:aa.set_ylim(2.6,-.6)
for aa in [la,nb,pv]:aa.axis('off')
for j,name in enumerate(['Biosynthesis','Catabolism','Biosynthesis − catabolism']):
 co=orange if j==2 else ink;la.text(0,j,name,fontsize=8.8,va='center',color=co,fontweight='bold' if j==2 else 'normal');ci(ax,rr.iloc[j],j,co,'D' if j==2 else 'o');nb.text(.5,j,f'{rr.iloc[j].beta:.4f}',ha='center',va='center',fontsize=8.5,color=co);pv.text(.5,j,pformat(rr.iloc[j].p),ha='center',va='center',fontsize=8.5,color=co)
ax.set_xlim(-2.1,2.8);ax.axvline(0,color='#9CA7AE',lw=.7,ls='--');ax.set_xticks([-2,0,2]);clean(ax);ax.set_xlabel('SONFH − exposed control',fontsize=8.3,labelpad=4)
f.text(.39,.30,'Effect and 95% CI',fontsize=8.5);f.text(.755,.30,'Effect',ha='center',fontsize=8.5);f.text(.912,.30,'Nominal P',ha='center',fontsize=8.5)
f.text(.018,.018,'Contrast FDR = 0.0394 across 11 tests; component P values are exploratory.',fontsize=8.3,color=muted)
save(f,4)

# Figure 5: align the four panels, preserve each outcome's own units.
ss=pd.read_csv(D/'layer3_ssGSEA.csv');rem=pd.read_csv(D/'membership_models.csv');pm=pd.read_csv(D/'proxy_models.csv');proxy=pd.read_csv(D/'proxy_individual.csv');catpm=pm[pm.outcome=='Catabolism'].reset_index(drop=True)
f=plt.figure(figsize=(7.2,5.9))
def smallforest(labelbox,plotbox,statbox,results,labels,statkind,xlim,xticks):
 la=f.add_axes(labelbox);aa=f.add_axes(plotbox);st=f.add_axes(statbox);n=len(results)
 for q in [la,aa,st]:q.set_ylim(n-.35,-.65)
 la.axis('off');st.axis('off')
 for j,(_,row) in enumerate(results.iterrows()):
  la.text(.98,j,labels[j],ha='right',va='center',fontsize=8.3);ci(aa,row,j);st.text(.5,j,(f'{row.FDR:.4f}' if statkind=='FDR' else pformat(row.p)),ha='center',va='center',fontsize=8.1)
 aa.axvline(0,color='#9CA7AE',ls='--',lw=.7);aa.set_xlim(xlim);aa.set_xticks(xticks);clean(aa);return aa
sc=pd.concat([pairs.iloc[[-1]],ss.iloc[[-1]]],ignore_index=True)
aa=smallforest([.018,.58,.105,.26],[.135,.58,.245,.26],[.393,.58,.083,.26],sc,['GSVA','ssGSEA'],'FDR',(-.5,4.2),[0,2,4]);aa.set_xlabel('Contrast effect and 95% CI',fontsize=8.2,labelpad=5)
bb=smallforest([.518,.58,.205,.26],[.735,.58,.165,.26],[.915,.58,.08,.26],rem,['Original set','Remove SNCA','Remove PRDX2','Remove SNCA + PRDX2','Remove MPO + EPX','Remove all four'],'P',(-.58,.03),[-.5,-.25,0]);bb.set_xlabel('GSVA effect and 95% CI',fontsize=8.2,labelpad=5)
heading(f,.018,.975,'A','Alternative scoring');heading(f,.518,.975,'B','Gene membership')
f.text(.434,.858,'FDR (11)',fontsize=8.1,ha='center');f.text(.955,.858,'Nominal P',fontsize=8.1,ha='center')
cc=f.add_axes([.09,.115,.36,.265]);is_case=proxy.group=='SONFH'
for mask,col,lab in [(~is_case,control,'Exposed non-SONFH'),(is_case,case,'SONFH')]:
 cc.scatter(proxy.loc[mask,'E1_four_gene_proxy'],proxy.loc[mask,'hydrogen_peroxide_catabolism_GSVA'],s=15,c=col,edgecolors='white',linewidths=.3,label=lab)
cc.spines[['top','right']].set_visible(False);cc.tick_params(labelsize=8,length=3);cc.set_xlabel('Four-gene erythroid expression proxy',fontsize=8.2);cc.set_ylabel('Catabolism GSVA score',fontsize=8.3);cc.legend(frameon=False,fontsize=7.6,loc='lower right',handletextpad=.3,borderpad=0)
dd=smallforest([.518,.115,.205,.265],[.735,.115,.165,.265],[.915,.115,.08,.265],catpm,['Sex adjusted','Four-gene proxy adjusted','Five-gene proxy (+HBB)'],'P',(-.6,.25),[-.5,0,.2]);dd.set_xlabel('GSVA effect and 95% CI',fontsize=8.2,labelpad=5)
heading(f,.018,.47,'C','Proxy–catabolism relationship');heading(f,.518,.47,'D','Proxy adjustment')
f.text(.09,.405,'Spearman ρ = 0.769; n = 40',fontsize=8.3,color=muted);f.text(.955,.395,'Nominal P',fontsize=8.1,ha='center')
f.text(.018,.018,'Proxy adjustment is exploratory; the proxy is not a measured blood-cell fraction.',fontsize=8.3,color=muted)
save(f,5)

# Figure 6: retain native-R vector points/mean bars, relocate labels to a compact template.
# No individual values for the DUOX-removal analysis are inferred or digitized.
source_path=D/'Figure6_native_R_source.pdf'
if not source_path.exists():
 candidates=sorted((ROOT.parent/'native_R'/'Figure6'/'results').glob('*/figure/Figure6_R.pdf'),key=lambda q:q.stat().st_mtime)
 if not candidates:
  raise SystemExit('Figures 1-5 exported. For Figure 6, run native_R/Figure6/RUN_FIGURE6.R first, then rerun this optional Python layout script. The redundant archived native-source PDF is not distributed.')
 source_path=candidates[-1]
source=fitz.open(source_path);doc=fitz.open();p=doc.new_page(width=518.4,height=279)
from matplotlib.font_manager import findfont,FontProperties
font_regular=findfont(FontProperties(family='DejaVu Sans'));font_bold=findfont(FontProperties(family='DejaVu Sans',weight='bold'))
p.insert_font(fontname='Body',fontfile=font_regular);p.insert_font(fontname='Bold',fontfile=font_bold)
def ptext(x,y,txt,size=8.5,bold=False,color=(.15,.21,.27)):
 p.insert_text((x,y),txt,fontsize=size,fontname='Bold' if bold else 'Body',color=color)
def center(x,y,txt,size=8.5):
 w=fitz.Font(fontfile=font_regular).text_length(txt,fontsize=size);ptext(x-w/2,y,txt,size)
for j in range(2):
 dx=j*253;src=fitz.Rect(76.32+j*362.88,139,349.92+j*362.88,304)
 dst=fitz.Rect(50+dx,73,248+dx,211)
 p.show_pdf_page(dst,source,0,clip=src,keep_proportion=False)
 ptext(50+dx,17,('A  Original gene sets' if j==0 else 'B  Without DUOX1 / DUOX2'),9.4,True)
 ptext(50+dx,36,('Mean difference = 1.2937' if j==0 else 'Mean difference = 1.3341'))
 ptext(50+dx,52,('Exact P = 0.0397 (5/126)' if j==0 else 'Exact P = 0.0556 (7/126)'),8.5,color=(.39,.45,.49))
 for old_y,label in [(154.18,'1.0'),(181.29,'0.5'),(208.41,'0.0'),(235.53,'−0.5'),(262.65,'−1.0'),(289.76,'−1.5')]:
  yy=73+(old_y-139)/165*138;p.draw_line((47+dx,yy),(50+dx,yy),color=(.48,.54,.57),width=.6)
  ptext(26+dx,yy+3,label,8.2)
 for old_x,label,nlab in [(144.72,'Fracture controls','(n = 4)'),(281.52,'SONFH','(n = 5)')]:
  xx=50+dx+(old_x-76.32)/273.6*198;center(xx,230,label);center(xx,243,nlab)
p.insert_text((12,210),'Biosynthesis z − catabolism z',fontsize=8.5,fontname='Body',rotate=90,color=(.15,.21,.27))
ptext(50,269,'GSE345158; points, individual samples; horizontal bars, group means.',8.1,color=(.39,.45,.49))
doc.save(O/'Figure6.pdf',garbage=4,deflate=True)
p=doc[0];preview=p.get_pixmap(matrix=fitz.Matrix(300/72,300/72),alpha=False);preview.set_dpi(300,300);preview.save(O/'Figure6.png')
pix=p.get_pixmap(matrix=fitz.Matrix(600/72,600/72),alpha=False)
Image.frombytes('RGB',[pix.width,pix.height],pix.samples).save(O/'Figure6.tiff',compression='tiff_lzw',dpi=(600,600))
(O/'Figure6.svg').write_text(p.get_svg_image())
doc.close();source.close()
print('All six figures rendered. Statistical values unchanged; Figure6 vector marks retained.')

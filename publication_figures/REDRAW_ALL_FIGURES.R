# Layout-only redraw of the six checked figures. Open in RStudio and click Source.
# Statistical calculations remain in the original full-analysis packages.
# Base R/grid are used for layout. Figure 6 also needs the existing GSVA 2.4.9.
redraw_all <- function() {
  candidates<-character()
  for(fr in sys.frames())if(!is.null(fr$ofile))candidates<-c(candidates,as.character(fr$ofile))
  args<-commandArgs(FALSE);candidates<-c(candidates,sub('^--file=','',args[grepl('^--file=',args)]))
  roots<-c(if(length(candidates))dirname(normalizePath(candidates,mustWork=FALSE)),getwd())
  roots<-roots[file.exists(file.path(roots,'plot_data','GSVA_scores.csv'))]
  if(!length(roots)) {
    if(!interactive())stop('Cannot locate extracted package.')
    message('Select REDRAW_ALL_FIGURES.R in the extracted package.');roots<-dirname(file.choose())
  }
  root<-normalizePath(roots[1],winslash='/',mustWork=TRUE);d<-file.path(root,'plot_data')
  out<-file.path(root,'R_redrawn');dir.create(out,showWarnings=FALSE)
  read<-function(n,...)read.csv(file.path(d,n),check.names=FALSE,stringsAsFactors=FALSE,fileEncoding='UTF-8-BOM',...)
  meta<-read('blood_manifest.csv');gs<-as.matrix(read('GSVA_scores.csv',row.names=1))
  l1<-read('layer1_all85.csv');l3<-read('layer3_GSVA.csv');ss<-read('layer3_ssGSEA.csv')
  families<-read('family_comparisons.csv');families<-families[order(families$domain_number,families$functional_group,method='radix'),]
  comp<-read('components.csv');removal<-read('membership_models.csv');pm<-read('proxy_models.csv');proxy<-read('proxy_individual.csv')
  stopifnot(identical(colnames(gs),meta$GSM),identical(dim(gs),c(85L,40L)),nrow(l1)==85L,nrow(l3)==11L,
    nrow(families)==21L,sum(l1$FDR<.05)==8L,sum(l3$FDR<.05)==1L,sum(families$FDR<.05)==0L)
  stopifnot(max(abs(p.adjust(l1$p,'BH')-l1$FDR))<1e-12,max(abs(p.adjust(l3$p,'BH')-l3$FDR))<1e-12)
  ink<-'#263746';muted<-'#65727C';accent<-'#B5632D';control<-'#368A88';case<-'#8368A7';light<-'#E0E6EA'
  txt<-function(x,y,s,sz=8.5,bold=FALSE,col=ink,ha=0,va=1,rot=0)
    grid::grid.text(s,x=x,y=y,just=c(ha,va),rot=rot,gp=grid::gpar(fontfamily='sans',fontsize=sz,
      fontface=if(bold)'bold' else 'plain',col=col,lineheight=1.15))
  seg<-function(x1,y1,x2,y2,col=ink,lw=.7,lty=1)
    grid::grid.lines(x=c(x1,x2),y=c(y1,y2),gp=grid::gpar(col=col,lwd=lw,lty=lty))
  rect<-function(x,y,w,h,fill,border=NA,lw=.6,lty=1)
    grid::grid.rect(x=x+w/2,y=y+h/2,width=w,height=h,gp=grid::gpar(fill=fill,col=border,lwd=lw,lty=lty))
  point<-function(x,y,col=ink,r=1.8,diamond=FALSE) {
    if(diamond)grid::grid.points(x=x,y=y,pch=18,size=grid::unit(r*2.6,'pt'),gp=grid::gpar(col=col))
    else grid::grid.circle(x=x,y=y,r=grid::unit(r,'pt'),gp=grid::gpar(fill=col,col=NA))
  }
  title<-function(x,y,letter,label){txt(x,y,letter,11,TRUE);txt(x+.034,y,label,9.5,TRUE)}
  pformat<-function(p)ifelse(p<.0001,formatC(p,format='e',digits=2),sprintf('%.4f',p))
  mapx<-function(v,x,w,limits)x+(v-limits[1])/diff(limits)*w
  rowy<-function(i,n,y,h)y+h-(i-.5)/n*h
  forest_axis<-function(x,y,w,h,limits,ticks,label) {
    seg(x,y,x+w,y,col=muted)
    for(v in ticks){xx<-mapx(v,x,w,limits);seg(xx,y,xx,y-.008,col=muted);txt(xx,y-.012,format(v,trim=TRUE),8,ha=.5)}
    zero<-mapx(0,x,w,limits);if(zero>=x&&zero<=x+w)seg(zero,y,zero,y+h,col='#9CA7AE',lty=2)
    txt(x+w/2,y-.055,label,8.2,ha=.5)
  }
  ci<-function(r,i,n,x,y,w,h,limits,col=ink,diamond=FALSE) {
    yy<-rowy(i,n,y,h);lo<-mapx(r$CI_lower,x,w,limits);hi<-mapx(r$CI_upper,x,w,limits)
    seg(lo,yy,hi,yy,col,1.1);seg(lo,yy-.005,lo,yy+.005,col,.8);seg(hi,yy-.005,hi,yy+.005,col,.8)
    point(mapx(r$beta,x,w,limits),yy,col,1.8,diamond)
  }
  savefig<-function(number,width,height,fun) {
    draw<-function(){grid::grid.newpage();fun()}
    dev<-function(opener){opener();on.exit(grDevices::dev.off(),add=TRUE);draw()}
    path<-function(ext)file.path(out,paste0('Figure',number,'.',ext))
    dev(function(){if(capabilities('cairo'))grDevices::cairo_pdf(path('pdf'),width=width,height=height,family='sans') else grDevices::pdf(path('pdf'),width=width,height=height)})
    dev(function()grDevices::png(path('png'),width=width,height=height,units='in',res=300,bg='white'))
    tryCatch(dev(function()grDevices::tiff(path('tiff'),width=width,height=height,units='in',res=600,compression='lzw',bg='white')),
      error=function(e)warning('TIFF failed: ',conditionMessage(e)))
    message('Rendered Figure ',number)
  }

  # Figure 1: fixed study design; counts are checked against archived definitions.
  savefig(1,7.2,6.2,function(){
    xx<-function(x)x/7.2;yy<-function(y)y/6.2
    t<-function(x,y,s,sz=8.5,bold=FALSE,ha=0,col=ink)txt(xx(x),yy(y),s,sz,bold,col,ha)
    box<-function(x,y,w,h,heading,body,fill='#EAF1F5',color='#315D77',lty=1){
      grid::grid.roundrect(x=xx(x+w/2),y=yy(y+h/2),width=xx(w),height=yy(h),r=grid::unit(.04,'in'),gp=grid::gpar(fill=fill,col='#81939D',lwd=.7,lty=lty))
      t(x+.13,y+h-.13,heading,8.8,TRUE,col=color);t(x+.13,y+h-.43,body,8.3)
    }
    arrow<-function(x1,y1,x2,y2)grid::grid.lines(x=xx(c(x1,x2)),y=yy(c(y1,y2)),arrow=grid::arrow(length=grid::unit(.07,'in'),type='closed'),gp=grid::gpar(col='#9CABB4',fill='#9CABB4',lwd=.8))
    t(.16,6.08,'A  Study inputs and functional definitions',10,TRUE)
    box(.16,4.55,3.31,1.13,'Blood discovery | GSE123568','30 SONFH + 10 GC-exposed non-SONFH\nDeposited RMA + GPL15207 annotation\nProbe-to-gene median: 18,833 genes x 40')
    box(3.68,4.55,3.36,1.13,'Literature-guided functional panel','6 predefined domains; 85 frozen gene sets\n28 functional groups; 21 with at least 2 sets\n11 predefined opposing-process pairs')
    arrow(1.82,4.53,1.82,4.33);arrow(5.36,4.53,5.36,4.33)
    seg(xx(1.82),yy(4.33),xx(5.36),yy(4.33),col='#9CABB4');arrow(3.6,4.33,3.6,4.12)
    t(3.6,4.08,'GSVA scoring across all 40 blood samples',9.5,TRUE,.5)
    t(.16,3.73,'B  Three prespecified parallel questions',10,TRUE)
    seg(xx(1.26),yy(3.40),xx(5.95),yy(3.40),col='#9CABB4')
    for(x in c(1.26,3.6,5.95))arrow(x,3.40,x,3.19)
    box(.16,1.60,2.20,1.55,'Layer 1 | Association','Which functions are\nassociated with SONFH?\n\n85 sets; sex-adjusted HC3\nBH across 85 tests')
    box(2.51,1.60,2.20,1.55,'Layer 2 | Effect differences','Do related components have\ndistinguishable effects?\n\n21 groups; standardized scores\nSex-adjusted CR2; BH across 21','#EAF3F0','#286C60',2)
    box(4.86,1.60,2.18,1.55,'Layer 3 | Process contrasts','Are opposing-process\ndifferences associated?\n\n11 contrasts: z(A) - z(B)\nSex-adjusted HC3; BH across 11','#F8F0E5','#8A622E',3)
    t(3.6,1.44,'Distinct questions; no sequential significance filter.',9,TRUE,.5)
    t(.16,1.09,'C  Sensitivity analyses and exploratory follow-up',10,TRUE)
    box(.16,.10,3.31,.78,'Blood-based characterization','Alternative scoring; membership; leave-one-out\nErythroid expression proxy; component analysis')
    box(3.68,.10,3.36,.78,'Exploratory tissue assessment','GSE345158: 5 SONFH + 4 fracture controls\nExact permutation; shared-gene sensitivity')
  })

  ids<-c('GOBP_HYDROGEN_PEROXIDE_CATABOLIC_PROCESS','GOBP_GLUTATHIONE_METABOLIC_PROCESS',
    'GOBP_GLUTATHIONE_CATABOLIC_PROCESS','GOBP_SUPEROXIDE_ANION_GENERATION','GOBP_T_CELL_ACTIVATION',
    'GOBP_T_CELL_RECEPTOR_SIGNALING_PATHWAY','GOBP_B_CELL_RECEPTOR_SIGNALING_PATHWAY','GOBP_GLYCOLIPID_TRANSPORT')
  labels<-c('Hydrogen peroxide catabolism','Glutathione metabolism','Glutathione catabolism',
    'Superoxide anion generation','T-cell activation','T-cell receptor signaling','B-cell receptor signaling','Glycolipid transport')
  stopifnot(setequal(l1$gene_set[l1$FDR<.05],ids))
  rows<-l1[match(ids,l1$gene_set),];sample_order<-c(meta$GSM[meta$group=='GC_exposed_non_SONFH'],meta$GSM[meta$group=='SONFH'])
  z<-t(scale(t(gs[ids,sample_order,drop=FALSE])));stopifnot(max(abs(z))<=2)
  pal<-grDevices::colorRampPalette(c('#2166AC','#F7F7F7','#B35806'))(201)
  savefig(2,7.2,4.35,function(){
    title(.018,.975,'A','Sex-adjusted associations');title(.643,.975,'B','Individual profiles')
    txt(.018,.89,'8 gene sets with FDR < 0.05 across 85 tests',8.5,col=muted)
    txt(.346,.823,'Effect and 95% CI',8.5);txt(.572,.823,'FDR',8.5,ha=.5)
    txt(.643,.90,'Exposed non-SONFH (n = 10)',8.1,col=control);txt(.643,.852,'SONFH (n = 30)',8.1,col=case)
    for(i in 1:8){y<-rowy(i,8,.22,.56);txt(.018,y,labels[i],8.5,va=.5);ci(rows[i,],i,8,.338,.22,.185,.56,c(-.95,.95));txt(.572,y,pformat(rows$FDR[i]),8.2,ha=.5,va=.5)
      for(j in 1:40){col<-pal[1+round((z[i,j]+2)/4*200)];rect(.643+(j-1)*.34/40,.22+(8-i)*.56/8,.34/40,.56/8,col)}
    }
    for(y in c(.22+.56*4/8,.22+.56/8)){seg(.018,y,.610,y,col=light,lw=.55);seg(.643,y,.983,y,col='white',lw=.6)}
    seg(.643+.34/4,.22,.643+.34/4,.78,col='white',lw=.9)
    rect(.643,.797,.34/4,.023,control);rect(.643+.34/4,.797,.34*3/4,.023,case)
    forest_axis(.338,.22,.185,.56,c(-.95,.95),c(-.8,0,.8),'GSVA score difference\nSONFH - exposed control')
    for(j in 1:201)rect(.677+(j-1)*.272/201,.115,.272/201,.025,pal[j])
    for(v in c(-2,0,2))txt(.677+(v+2)/4*.272,.102,as.character(v),8,ha=.5)
    txt(.813,.055,'Gene-set score z',8.5,ha=.5)
  })

  plabels<-c('Coagulation - fibrinolysis','Endothelial proliferation\n- apoptosis','Fatty acid biosynthesis\n- beta-oxidation',
    'Triglyceride storage\n- catabolism','Cholesterol import\n- efflux','Phospholipid biosynthesis\n- catabolism',
    'Ceramide biosynthesis\n- catabolism','Glycolipid biosynthesis\n- catabolism','Lipoprotein biosynthesis\n- catabolism',
    'Lipid-droplet formation\n- disassembly','Hydrogen peroxide\nbiosynthesis - catabolism')
  flabels<-families$functional_group
  flabels[flabels=='Angiogenesis / endothelial remodeling']<-'Angiogenesis / endothelium'
  flabels[flabels=='Osteoclast differentiation / bone resorption']<-'Osteoclasts / bone resorption'
  savefig(3,7.2,5.75,function(){
    title(.018,.975,'A','Related functional groups');title(.53,.975,'B','Opposing-process contrasts')
    txt(.018,.915,'0 of 21 met FDR < 0.05',8.5,col=muted);txt(.53,.915,'1 of 11 met FDR < 0.05',8.5,col=muted)
    txt(.33,.85,'FDR',8.5);txt(.471,.85,'FDR',8.5,ha=.5);txt(.76,.85,'Effect and 95% CI',8.3);txt(.958,.85,'FDR',8.5,ha=.5)
    for(i in 1:21){y<-rowy(i,21,.15,.68);txt(.018,y,flabels[i],8.1,va=.5);point(mapx(families$FDR[i],.321,.108,c(0,1.03)),y,r=1.4);txt(.471,y,sprintf('%.4f',families$FDR[i]),8.1,ha=.5,va=.5)
      if(i<21&&families$domain_number[i]!=families$domain_number[i+1])seg(.018,y-.68/42,.504,y-.68/42,col=light,lw=.5)
    }
    forest_axis(.321,.15,.108,.68,c(0,1.03),c(0,.5,1),'FDR')
    seg(mapx(.05,.321,.108,c(0,1.03)),.15,mapx(.05,.321,.108,c(0,1.03)),.83,col=accent,lty=3)
    for(i in 1:11){r<-l3[i,];col<-if(r$FDR<.05)accent else '#75838D';y<-rowy(i,11,.15,.68)
      txt(.53,y,plabels[i],8.3,r$FDR<.05,if(r$FDR<.05)accent else ink,va=.5)
      ci(r,i,11,.76,.15,.15,.68,c(-2.05,2.9),col,r$FDR<.05);txt(.958,y,sprintf('%.4f',r$FDR),8.1,r$FDR<.05,col,ha=.5,va=.5)
    }
    forest_axis(.76,.15,.15,.68,c(-2.05,2.9),c(-2,0,2),'Disease-associated effect\nin z(A) - z(B)')
    txt(.018,.035,'Separate prespecified correction families: 21 functional groups and 11 contrasts.',8.5,col=muted)
  })

  values<-as.matrix(read('Figure4_participant_values.csv',row.names=1))
  r4<-rbind(comp[,c('beta','CI_lower','CI_upper','p')],l3[11,c('beta','CI_lower','CI_upper','p')])
  iscase<-meta$group[match(colnames(values),meta$GSM)]=='SONFH'
  yr<-c(floor(min(values)-.2),ceiling(max(values)+.2))
  savefig(4,7.2,5.35,function(){
    title(.018,.975,'A','Individual score distributions')
    positions<-c(.085,.405,.725);titles<-c('Hydrogen peroxide\nbiosynthesis','Hydrogen peroxide\ncatabolism','Biosynthesis - catabolism')
    for(j in 1:3){x<-positions[j];w<-.245;y<-.47;h<-.35;my<-function(v)y+(v-yr[1])/diff(yr)*h
      txt(x+w/2,.87,titles[j],8.8,TRUE,ha=.5)
      seg(x,my(0),x+w,my(0),col=light,lty=3)
      if(j==1){seg(x,y,x,y+h,col=muted);for(v in seq(yr[1],yr[2],by=2))txt(x-.012,my(v),as.character(v),8,ha=1,va=.5)
        txt(.028,y+h/2,'Component z or contrast',8.5,ha=.5,va=.5,rot=90)}
      for(k in 1:2){v<-values[j,if(k==1)!iscase else iscase];col<-c(control,case)[k];cx<-x+if(k==1)w*.237 else w*.763
        st<-grDevices::boxplot.stats(v)$stats;bw<-w*.25
        seg(cx,my(st[1]),cx,my(st[5]),col,.7);seg(cx-bw/4,my(st[1]),cx+bw/4,my(st[1]),col,.7);seg(cx-bw/4,my(st[5]),cx+bw/4,my(st[5]),col,.7)
        rect(cx-bw/2,my(st[2]),bw,my(st[4])-my(st[2]),grDevices::adjustcolor(col,alpha.f=.13),col,.8)
        seg(cx-bw/2,my(st[3]),cx+bw/2,my(st[3]),ink,1)
        offs<-((seq_along(v)-1)%%7-3)*.035*w/1.9
        for(i in seq_along(v))point(cx+offs[i],my(v[i]),col,1.05)
        txt(cx,.425,if(k==1)'Exposed\ncontrol' else 'SONFH',8.3,ha=.5)
      }
    }
    title(.018,.365,'B','Sex-adjusted associations')
    for(i in 1:3){col<-if(i==3)accent else ink;y<-rowy(i,3,.105,.18)
      txt(.018,y,c('Biosynthesis','Catabolism','Biosynthesis - catabolism')[i],8.8,i==3,col,va=.5)
      ci(r4[i,],i,3,.37,.105,.31,.18,c(-2.1,2.8),col,i==3)
      txt(.755,y,sprintf('%.4f',r4$beta[i]),8.5,col=col,ha=.5,va=.5);txt(.912,y,pformat(r4$p[i]),8.5,col=col,ha=.5,va=.5)
    }
    txt(.39,.31,'Effect and 95% CI',8.5);txt(.755,.31,'Effect',8.5,ha=.5);txt(.912,.31,'Nominal P',8.5,ha=.5)
    forest_axis(.37,.105,.31,.18,c(-2.1,2.8),c(-2,0,2),'SONFH - exposed control')
    txt(.018,.025,'Contrast FDR = 0.0394 across 11 tests; component P values are exploratory.',8.3,col=muted)
  })

  savefig(5,7.2,5.9,function(){
    small<-function(r,labels,lx,x,w,statx,y,h,limits,ticks,fdr=FALSE){
      n<-nrow(r)
      for(i in seq_len(n)){yy<-rowy(i,n,y,h);txt(lx,yy,labels[i],8.3,ha=1,va=.5);ci(r[i,],i,n,x,y,w,h,limits)
        txt(statx,yy,if(fdr)sprintf('%.4f',r$FDR[i]) else pformat(r$p[i]),8.1,ha=.5,va=.5)}
      forest_axis(x,y,w,h,limits,ticks,if(fdr)'Contrast effect and 95% CI' else 'GSVA effect and 95% CI')
    }
    title(.018,.975,'A','Alternative scoring');title(.518,.975,'B','Gene membership')
    txt(.434,.858,'FDR (11)',8.1,ha=.5);txt(.955,.858,'Nominal P',8.1,ha=.5)
    small(rbind(l3[11,],ss[11,]),c('GSVA','ssGSEA'),.123,.135,.245,.434,.58,.26,c(-.5,4.2),c(0,2,4),TRUE)
    small(removal,c('Original set','Remove SNCA','Remove PRDX2','Remove SNCA + PRDX2','Remove MPO + EPX','Remove all four'),.723,.735,.165,.955,.58,.26,c(-.58,.03),c(-.5,-.25,0))
    title(.018,.47,'C','Proxy-catabolism relationship');title(.518,.47,'D','Proxy adjustment')
    txt(.09,.405,'Spearman rho = 0.769; n = 40',8.3,col=muted);txt(.955,.395,'Nominal P',8.1,ha=.5)
    xx<-proxy$E1_four_gene_proxy;vv<-proxy$hydrogen_peroxide_catabolism_GSVA
    xr<-range(xx)+c(-.12,.12);yr<-range(vv)+c(-.03,.03);x<-.09;y<-.115;w<-.36;h<-.265
    seg(x,y,x+w,y,col=muted);seg(x,y,x,y+h,col=muted)
    for(v in -2:2)txt(mapx(v,x,w,xr),y-.015,as.character(v),8,ha=.5)
    for(v in seq(-.4,.4,.2))txt(x-.012,y+(v-yr[1])/diff(yr)*h,sprintf('%.1f',v),8,ha=1,va=.5)
    for(i in seq_along(xx))point(mapx(xx[i],x,w,xr),y+(vv[i]-yr[1])/diff(yr)*h,if(proxy$group[i]=='SONFH')case else control,1.7)
    txt(x+w/2,y-.055,'Four-gene erythroid expression proxy',8.2,ha=.5)
    txt(.025,y+h/2,'Catabolism GSVA score',8.3,ha=.5,va=.5,rot=90)
    point(.288,.165,control,1.5);txt(.300,.165,'Exposed non-SONFH',7.6,va=.5)
    point(.288,.140,case,1.5);txt(.300,.140,'SONFH',7.6,va=.5)
    small(pm[pm$outcome=='Catabolism',],c('Sex adjusted','Four-gene proxy adjusted','Five-gene proxy (+HBB)'),.723,.735,.165,.955,.115,.265,c(-.6,.25),c(-.5,0,.2))
    txt(.018,.025,'Proxy adjustment is exploratory; the proxy is not a measured blood-cell fraction.',8.3,col=muted)
  })

  # Figure 6 uses the same corrected matrix and frozen sets as the successful tissue audit.
  # It only rescales/scores the two sets for display; it does not repeat VST or fit blood models.
  if(!requireNamespace('GSVA',quietly=TRUE))stop('Figures 1-5 complete. GSVA is needed for Figure 6.')
  if(as.character(packageVersion('GSVA'))!='2.4.9')stop('Figure 6 requires the archived GSVA 2.4.9.')
  obj<-readRDS(file.path(d,'tissue_reproducibility_inputs.rds'))
  expected<-c(paste0('Case-',1:5),'Ctrl-1','Ctrl-2','Ctrl-4','Ctrl-5')
  bio<-'GOBP_HYDROGEN_PEROXIDE_BIOSYNTHETIC_PROCESS';catab<-'GOBP_HYDROGEN_PEROXIDE_CATABOLIC_PROCESS'
  vv<-list();summ<-read('tissue_summary.csv')
  for(j in 1:2){sets<-obj$pair11_sets;if(j==2)sets<-lapply(sets,setdiff,y=c('DUOX1','DUOX2'))
    pa<-GSVA::gsvaParam(obj$gsva_matrix_sod2,sets,kcdf='Gaussian',tau=1,maxDiff=TRUE,absRanking=FALSE,minSize=1,maxSize=Inf)
    z<-t(scale(t(GSVA::gsva(pa,verbose=FALSE))));v<-z[bio,expected]-z[catab,expected];vv[[j]]<-v
    caseidx<-startsWith(expected,'Case-');dif<-mean(v[caseidx])-mean(v[!caseidx]);ix<-combn(9,5)
    perm<-apply(ix,2,function(i)mean(v[i])-mean(v[-i]));p<-mean(abs(perm)>=abs(dif)-1e-10)
    stopifnot(abs(dif-summ$Difference[j])<2e-6,abs(p-summ$Exact_P[j])<1e-8)
  }
  old<-read('tissue_original_scores.csv');old<-old[match(expected,old$Sample),]
  stopifnot(max(abs(vv[[1]]-old$Pair11))<1e-7)
  savefig(6,7.2,3.875,function(){
    yr<-range(unlist(vv));yr<-yr+c(-1,1)*diff(yr)*.14
    for(j in 1:2){dx<-(j-1)*.488;x<-.096+dx;y<-.244;w<-.382;h<-.495
      txt(x,.955,if(j==1)'A  Original gene sets' else 'B  Without DUOX1 / DUOX2',9.4,TRUE)
      txt(x,.895,paste0('Mean difference = ',sprintf('%.4f',summ$Difference[j])),8.5)
      txt(x,.838,if(j==1)'Exact P = 0.0397 (5/126)' else 'Exact P = 0.0556 (7/126)',8.5,col=muted)
      my<-function(v)y+(v-yr[1])/diff(yr)*h
      seg(x,my(-1.5),x,my(1),col=muted);seg(x,my(0),x+w,my(0),col=light,lty=3)
      for(v in seq(-1.5,1,.5)){seg(x-.006,my(v),x,my(v),col=muted);txt(x-.013,my(v),sprintf('%.1f',v),8.2,ha=1,va=.5)}
      for(k in 1:2){g<-vv[[j]][if(k==1)6:9 else 1:5];cx<-x+w*c(.25,.75)[k]
        off<-seq(-.05,.05,length.out=length(g))*w
        for(i in seq_along(g))point(cx+off[i],my(g[i]),c(control,case)[k],2.3)
        seg(cx-w*.12,my(mean(g)),cx+w*.12,my(mean(g)),ink,1.6)
        txt(cx,.18,if(k==1)'Fracture controls\n(n = 4)' else 'SONFH\n(n = 5)',8.5,ha=.5)
      }
    }
    txt(.025,.49,'Biosynthesis z - catabolism z',8.5,ha=.5,va=.5,rot=90)
    txt(.096,.055,'GSE345158; points, individual samples; horizontal bars, group means.',8.1,col=muted)
  })
  capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'))
  writeLines('PASS_LAYOUT_REDRAW_6_FIGURES: checked plot data; no blood model refitting.',file.path(out,'REDRAW_STATUS.txt'))
  message('PASS_LAYOUT_REDRAW_6_FIGURES\nFiles saved to: ',out)
  invisible(out)
}
redraw_all()

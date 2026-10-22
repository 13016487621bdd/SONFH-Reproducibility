draw_figure5 <- function(scoring,removal,proxy_data,proxy_models,rho,out) {
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  ink<-'#263746';muted<-'#65727C';orange<-'#B5632D';control<-'#368A88';case<-'#8368A7'
  pformat<-function(v)ifelse(v<.0001,formatC(v,format='e',digits=2),sprintf('%.4f',v))
  forest<-function(r,labels,title,xlim,stat,stat_title,xlab) {
    par(mar=c(4,10,4,6),xaxs='i',yaxs='i')
    n<-nrow(r)
    plot.new();plot.window(xlim=xlim,ylim=c(n+.7,.3))
    abline(v=0,col='#A5AFB6',lty=2)
    segments(r$CI_lower,1:n,r$CI_upper,1:n,col=ink,lwd=1.3)
    segments(r$CI_lower,(1:n)-.07,r$CI_lower,(1:n)+.07,col=ink)
    segments(r$CI_upper,(1:n)-.07,r$CI_upper,(1:n)+.07,col=ink)
    points(r$beta,1:n,pch=16,col=ink,cex=.85)
    axis(1,cex.axis=.8,col=muted)
    axis(2,at=1:n,labels=labels,las=1,tick=FALSE,cex.axis=.8)
    text(xlim[2]+diff(xlim)*.1,1:n,stat,adj=0,xpd=NA,cex=.76,col=ink)
    text(xlim[2]+diff(xlim)*.1,.35,stat_title,adj=c(0,0),xpd=NA,cex=.76,col=muted)
    mtext(title,side=3,line=2.2,adj=0,font=2,cex=.96)
    mtext(xlab,side=1,line=2.5,cex=.8)
  }
  catmodels<-proxy_models[proxy_models$outcome=='Catabolism',]
  stopifnot(identical(scoring$method,c('GSVA','ssGSEA')),nrow(removal)==6L,nrow(catmodels)==3L)
  draw<-function() {
    layout(matrix(1:4,nrow=2,byrow=TRUE))
    par(oma=c(2,.1,.4,.1),family='sans')
    forest(scoring,c('GSVA','ssGSEA'),'A  Alternative scoring',c(-.8,4.2),
      sprintf('%.4f',scoring$FDR),'FDR (11)', 'Contrast coefficient and 95% CI')
    forest(removal,c('Original set','Remove SNCA','Remove PRDX2','Remove SNCA + PRDX2',
      'Remove MPO + EPX','Remove all four'), 'B  Gene membership: catabolism',c(-.58,.06),
      pformat(removal$p),'Nominal P','GSVA score coefficient and 95% CI')
    par(mar=c(4,5,4,1),xaxs='r',yaxs='r')
    colors<-ifelse(proxy_data$group=='SONFH',case,control)
    plot(proxy_data$proxy4,proxy_data$catabolism,pch=21,bg=colors,col='white',lwd=.5,
      cex=1.1,xlab='',ylab='',axes=FALSE)
    axis(1,cex.axis=.8,col=muted);axis(2,cex.axis=.8,las=1,col=muted)
    mtext('Four-gene erythroid expression proxy',side=1,line=2.5,cex=.85)
    mtext('Catabolism GSVA score',side=2,line=3,cex=.85)
    mtext('C  Proxy and catabolism scores',side=3,line=2.2,adj=0,font=2,cex=.96)
    mtext(sprintf('Spearman rho = %.3f; all 40 participants',rho),side=3,line=.7,adj=0,cex=.8,col=muted)
    legend('bottomright',legend=c('Exposed non-SONFH','SONFH'),pch=21,pt.bg=c(control,case),
           col='white',bty='n',cex=.77)
    forest(catmodels,c('Sex adjusted','Four-gene proxy adjusted','Five-gene proxy (+HBB)'),
      'D  Proxy adjustment: catabolism',c(-.6,.25),pformat(catmodels$p),'Nominal P',
      'GSVA score coefficient and 95% CI')
    mtext('Sensitivity analyses: transcriptional features; the proxy is not a measured blood-cell fraction.',
          side=1,outer=TRUE,line=.6,cex=.8,col=muted)
  }
  save_plot<-function(open_device){open_device();on.exit(dev.off(),add=TRUE);draw()}
  save_plot(function()pdf(file.path(out,'Figure5_R.pdf'),width=11.5,height=8.5))
  save_plot(function()png(file.path(out,'Figure5_R.png'),width=11.5,height=8.5,units='in',res=300,bg='white'))
  tryCatch(save_plot(function()tiff(file.path(out,'Figure5_R.tiff'),width=11.5,height=8.5,
      units='in',res=600,compression='lzw',bg='white')),error=function(e)warning('TIFF export failed: ',conditionMessage(e)))
  write.csv(scoring,file.path(out,'Figure5A_scoring.csv'),row.names=FALSE)
  write.csv(removal,file.path(out,'Figure5B_membership.csv'),row.names=FALSE)
  write.csv(proxy_data,file.path(out,'Figure5C_proxy_scatter.csv'),row.names=FALSE)
  write.csv(catmodels,file.path(out,'Figure5D_proxy_models.csv'),row.names=FALSE)
  legend<-paste('Figure 5. Sensitivity and interpretation of the hydrogen peroxide-related transcriptional features.',
    '(A) Sex-adjusted coefficients and HC3-based nominal 95% confidence intervals for the biosynthesis-minus-catabolism contrast using GSVA or ssGSEA.',
    'Constituent scores were standardized separately for each scoring method across all 40 participants; contrasts were not restandardized.',
    'FDR values were calculated separately across all 11 prespecified contrasts under each method. These coefficients are not a formal test of differences between scoring methods.',
    '(B) Catabolism-associated scores recalculated after removal of SNCA, PRDX2, both genes, MPO and EPX, or all four genes. Original GSVA score units are used; nominal P values are exploratory.',
    'All four refers to SNCA, PRDX2, MPO and EPX. The scenarios do not constitute a new primary FDR family.',
    '(C) Four-gene expression proxy, defined as the mean of expression z scores for ALAS2, AHSP, GYPA and SLC4A1, versus the catabolism-associated GSVA score; the pooled Spearman coefficient is descriptive.',
    '(D) Sex-adjusted catabolism coefficients before and after additional adjustment for the four-gene proxy or an exploratory five-gene proxy including HBB.',
    'The four genes are outside the frozen catabolism set; HBB overlaps that set, so the five-gene analysis is a separate sensitivity analysis.',
    'Intervals and P values in B and D are nominal and unadjusted for multiplicity. All models compare established SONFH with glucocorticoid-exposed participants without SONFH.',
    'Proxy-related attenuation does not establish causal confounding by cellular composition, and the proxy was not validated against measured blood-cell counts.',
    'Fixed-score leave-one-out results are provided separately; no scores or standardization were recomputed within omissions.',
    'CI, confidence interval; FDR, false discovery rate; SONFH, steroid-induced osteonecrosis of the femoral head.')
  writeLines(legend,file.path(out,'Figure5_legend.txt'))
  invisible(out)
}

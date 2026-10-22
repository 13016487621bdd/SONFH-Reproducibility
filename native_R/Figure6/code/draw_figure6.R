draw_figure6<-function(data,summary,out) {
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  versions<-c('SOD2_corrected','Remove_DUOX1_DUOX2')
  stopifnot(identical(as.character(summary$Version),versions),nrow(data)==18L)
  ink<-'#263746';muted<-'#65727C';control<-'#368A88';case<-'#8368A7'
  yr<-range(data$Contrast);pad<-diff(yr)*.14;yr<-yr+c(-pad,pad)
  draw<-function() {
    par(mfrow=c(1,2),oma=c(3,.3,4.1,.3),family='sans',xaxs='i')
    for(j in 1:2) {
      par(mar=c(4.6,5,4.8,1.2))
      d<-data[data$Version==versions[j],]
      stopifnot(sum(d$group=='Control')==4L,sum(d$group=='SONFH')==5L)
      plot.new();plot.window(xlim=c(.5,2.5),ylim=yr)
      abline(h=0,lty=3,col='#D7DEE2')
      for(k in 1:2) {
        ix<-which(d$group==c('Control','SONFH')[k])
        offset<-seq(-.10,.10,length.out=length(ix))
        points(k+offset,d$Contrast[ix],pch=21,bg=c(control,case)[k],col='white',lwd=.5,cex=1.3)
        avg<-mean(d$Contrast[ix])
        segments(k-.23,avg,k+.23,avg,col=ink,lwd=2)
      }
      axis(1,at=1:2,labels=c('Fracture controls\n(n = 4)','SONFH\n(n = 5)'),tick=FALSE,cex.axis=.9)
      axis(2,las=1,cex.axis=.85,col=muted)
      mtext('Biosynthesis z - catabolism z',side=2,line=3,cex=.9)
      mtext(c('A  SOD2-corrected original sets','B  Without shared DUOX1 / DUOX2')[j],
        side=3,line=3.1,adj=0,font=2,cex=.98)
      mtext(sprintf('Mean difference = %.4f',summary$Difference[j]),side=3,line=1.7,adj=0,cex=.87,col=ink)
      mtext(sprintf('Two-sided exact P = %.4f (%d / %d)',summary$Exact_P[j],
        summary$Extreme_assignments[j],summary$Total_assignments[j]),side=3,line=.4,adj=0,cex=.87,col=muted)
    }
    mtext('Exploratory femoral head tissue assessment | GSE345158',side=3,outer=TRUE,line=1.8,
      adj=.5,font=2,cex=1.1,col=ink)
    mtext('Points: individual samples. Horizontal bars: group means. These controls are fracture controls.',
      side=1,outer=TRUE,line=1.2,cex=.8,col=muted)
  }
  save_plot<-function(open_device){open_device();on.exit(dev.off(),add=TRUE);draw()}
  save_plot(function()pdf(file.path(out,'Figure6_R.pdf'),width=10.2,height=5.9))
  save_plot(function()png(file.path(out,'Figure6_R.png'),width=10.2,height=5.9,units='in',res=300,bg='white'))
  tryCatch(save_plot(function()tiff(file.path(out,'Figure6_R.tiff'),width=10.2,height=5.9,
     units='in',res=600,compression='lzw',bg='white')),error=function(e)warning('TIFF export failed: ',conditionMessage(e)))
  legend<-paste('Figure 6. Exploratory assessment of the hydrogen peroxide contrast in femoral head tissue.',
    'Individual biosynthesis-minus-catabolism contrasts are shown for five SONFH samples and four fracture controls from GSE345158.',
    '(A) Original frozen gene sets scored using the symbol-mapped, SOD2-corrected expression matrix.',
    '(B) Both sets rescored after removing their shared DUOX1 and DUOX2 genes.',
    'GSVA component scores were standardized separately within each version across all nine tissue samples; the differences were not restandardized.',
    'Points denote individual samples; horizontal bars denote group means. No confidence interval is implied by the bars.',
    'Mean differences compare SONFH with fracture controls. Two-sided exact P values enumerate all 126 allocations of five case labels among nine samples and count absolute mean differences at least as large as observed, including the observed allocation.',
    'P values are nominal exploratory tests and were not included in the blood primary FDR families.',
    'Removal of DUOX1 and DUOX2 did not reduce the observed mean difference, although the exact P value crossed the nominal 0.05 threshold.',
    'The small tissue sample and different control population support exploratory assessment of direction rather than matched independent validation of a blood biomarker.',
    'Tissue scores do not directly measure enzyme activity, peroxide clearance or physiological redox balance.',
    'SONFH, steroid-induced osteonecrosis of the femoral head; GSVA, gene set variation analysis.')
  writeLines(legend,file.path(out,'Figure6_legend.txt'))
  invisible(out)
}

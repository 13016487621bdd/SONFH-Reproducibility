# Panel A: unadjusted participant scores. Panel B: adjusted HC3 coefficients.
draw_figure4 <- function(values, estimates, meta, out) {
  stopifnot(identical(rownames(values), c('Biosynthesis','Catabolism','Contrast')),
            identical(colnames(values), meta$GSM), identical(estimates$outcome, rownames(values)),
            all(is.finite(values)), all(estimates$CI_lower < estimates$beta & estimates$beta < estimates$CI_upper))
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  ink <- '#263746'; muted <- '#65727C'; orange <- '#B5632D'
  control <- '#368A88'; case <- '#8368A7'
  group <- factor(meta$group, levels = c('GC_exposed_non_SONFH','SONFH'))
  titles <- c('Hydrogen peroxide biosynthesis', 'Hydrogen peroxide catabolism',
              'Biosynthesis - catabolism')
  yr <- range(values); pad <- diff(yr)*.12; yr <- yr + c(-pad,pad)
  draw <- function() {
    layout(matrix(1:6, nrow = 2, byrow = TRUE), widths = c(3.5,3.4,3.2), heights = c(3.5,2.1))
    par(oma = c(2.8,.3,4.5,.3), family = 'sans', xaxs = 'i')
    for (j in 1:3) {
      par(mar = c(3.2, if (j==1) 4 else 2, 2, 1))
      boxplot(values[j, ] ~ group, names = c('Exposed\ncontrol','SONFH'),
              col = grDevices::adjustcolor(c(control,case), alpha.f = .15), border = c(control,case),
              outline = FALSE, ylim = yr, xlim = c(.5,2.5), boxwex = .5, las = 1,
              axes = FALSE, ylab = '', xlab = '', medlwd = 1.5)
      abline(h = 0, lty = 3, col = '#D4DBDF')
      for (k in 1:2) {
        ix <- which(as.integer(group)==k)
        offset <- ((seq_along(ix)-1) %% 7 - 3)*.035
        points(k+offset, values[j, ix], pch = 21, bg = c(control,case)[k],
               col = 'white', lwd = .4, cex = .78)
      }
      axis(1, at = 1:2, labels = c('Exposed\ncontrol','SONFH'), tick = FALSE, cex.axis = .85)
      if (j==1) {
        axis(2, las = 1, cex.axis = .8, col = muted)
        mtext('Score (component z or difference)', side = 2, line = 2.7, cex = .8)
      }
      mtext(titles[j], side = 3, line = .5, cex = .88, font = 2, col = ink)
    }
    setup <- function(limits) {
      par(mar = c(3, .5, 2.5, .5))
      plot.new();plot.window(xlim = limits, ylim = c(3.7,.3))
    }
    setup(c(0,1))
    text(0,1:3,c('Biosynthesis','Catabolism','Biosynthesis - catabolism'),
         adj = c(0,.5), cex = 1, col = c(ink,ink,orange), font = c(1,1,2))
    mtext('B  Sex-adjusted associations', side = 3, line = 1.1, adj = 0, cex = 1, font = 2)
    setup(c(-2.1,2.8))
    colors <- c(muted,ink,orange)
    abline(v = 0, col = '#9AA4AB', lty = 2)
    segments(estimates$CI_lower,1:3,estimates$CI_upper,1:3,col = colors,lwd = 1.4)
    segments(estimates$CI_lower,(1:3)-.07,estimates$CI_lower,(1:3)+.07,col = colors)
    segments(estimates$CI_upper,(1:3)-.07,estimates$CI_upper,(1:3)+.07,col = colors)
    points(estimates$beta,1:3,pch = c(16,16,18),col = colors,cex = c(.9,.9,1.1))
    axis(1,at = -2:2,cex.axis = .85,col = muted)
    mtext('SONFH - exposed control',side = 1,line = 2,cex = .85)
    mtext('Coefficient and 95% CI',side = 3,line = .4,cex = .85)
    setup(c(0,1))
    ptext <- ifelse(estimates$p < .0001,formatC(estimates$p,format='e',digits=2),sprintf('%.4f',estimates$p))
    text(.28,1:3,sprintf('%.4f',estimates$beta),cex=.9,col=colors)
    text(.78,1:3,ptext,cex=.9,col=colors)
    mtext('Coefficient',side=3,line=.4,at=.28,cex=.85)
    mtext('Nominal P',side=3,line=.4,at=.78,cex=.85)
    mtext('A  Individual scores: 10 exposed controls and 30 SONFH cases',outer=TRUE,side=3,
          line=2.5,adj=0,cex=1,font=2,col=ink)
    mtext(sprintf('Contrast FDR = %.4f (BH across 11 contrasts); component P values are exploratory.',
                  estimates$primary_11_test_FDR[3]),outer=TRUE,side=1,line=1,cex=.83,col=muted)
  }
  save_plot <- function(open_device) {open_device();on.exit(dev.off(),add=TRUE);draw()}
  save_plot(function() pdf(file.path(out,'Figure4_R.pdf'),width=11,height=7.3))
  save_plot(function() png(file.path(out,'Figure4_R.png'),width=11,height=7.3,units='in',res=300,bg='white'))
  tryCatch(save_plot(function() tiff(file.path(out,'Figure4_R.tiff'),width=11,height=7.3,
     units='in',res=600,compression='lzw',bg='white')),error=function(e) warning('TIFF export failed: ',conditionMessage(e)))
  legend <- paste('Figure 4. Component-level characterization of the hydrogen peroxide contrast.',
    '(A) Individual biosynthesis-associated scores, catabolism-associated scores, and biosynthesis-minus-catabolism contrasts in 10 glucocorticoid-exposed participants without SONFH and 30 participants with established SONFH.',
    'Each component was standardized using its mean and sample standard deviation across all 40 participants; the contrast was not restandardized.',
    'Boxes show the median and interquartile range; whiskers extend to the most extreme observations within 1.5 interquartile ranges. All participant values are overlaid with deterministic horizontal offsets.',
    'These distributions are descriptive and unadjusted. Shared vertical limits permit comparison of the plotted score units.',
    '(B) Disease coefficients and HC3-based nominal 95% confidence intervals from separate sex-adjusted regressions of the two components and their contrast.',
    'Component P values are exploratory and unadjusted for multiplicity. The contrast FDR was calculated across all 11 prespecified opposing-process contrasts, not across the three displayed rows.',
    'The contrast coefficient equals the biosynthesis coefficient minus the catabolism coefficient; their standard errors and confidence intervals cannot be subtracted in this manner.',
    'The association was predominantly attributable to the lower catabolism-associated component and does not establish increased peroxide generation, impaired clearance, physiological redox imbalance, or superior diagnostic performance.',
    'CI, confidence interval; FDR, false discovery rate; SONFH, steroid-induced osteonecrosis of the femoral head.')
  writeLines(legend,file.path(out,'Figure4_legend.txt'))
  invisible(out)
}

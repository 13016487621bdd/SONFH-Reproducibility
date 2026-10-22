# Base-R Figure 2 renderer. Receives freshly fitted results from RUN_FIGURE2.R.
draw_figure2 <- function(scores, res, meta, out) {
  ids <- c('GOBP_HYDROGEN_PEROXIDE_CATABOLIC_PROCESS',
           'GOBP_GLUTATHIONE_METABOLIC_PROCESS', 'GOBP_GLUTATHIONE_CATABOLIC_PROCESS',
           'GOBP_SUPEROXIDE_ANION_GENERATION', 'GOBP_T_CELL_ACTIVATION',
           'GOBP_T_CELL_RECEPTOR_SIGNALING_PATHWAY', 'GOBP_B_CELL_RECEPTOR_SIGNALING_PATHWAY',
           'GOBP_GLYCOLIPID_TRANSPORT')
  labels <- c('Hydrogen peroxide catabolism', 'Glutathione metabolism', 'Glutathione catabolism',
              'Superoxide anion generation', 'T-cell activation', 'T-cell receptor signaling',
              'B-cell receptor signaling', 'Glycolipid transport')
  stopifnot(nrow(res) == 85L, sum(res$FDR < .05) == 8L,
            setequal(res$gene_set[res$FDR < .05], ids),
            max(abs(p.adjust(res$p, method = 'BH') - res$FDR)) < 1e-12,
            identical(colnames(scores), meta$GSM))
  sample_order <- c(meta$GSM[meta$group == 'GC_exposed_non_SONFH'], meta$GSM[meta$group == 'SONFH'])
  stopifnot(length(sample_order) == 40L, sum(meta$group == 'GC_exposed_non_SONFH') == 10L)
  estimates <- res[match(ids, res$gene_set), ]
  x <- scores[ids, sample_order, drop = FALSE]
  z <- t(scale(t(x)))
  stopifnot(all(is.finite(z)), max(abs(rowMeans(z))) < 1e-12,
            max(abs(apply(z, 1, sd) - 1)) < 1e-12, max(abs(z)) <= 2,
            all(estimates$CI_lower < estimates$beta & estimates$beta < estimates$CI_upper))
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  ink <- '#263746'; muted <- '#65727C'; control <- '#368A88'; case <- '#8368A7'
  colors <- colorRampPalette(c('#2166AC', '#F7F7F7', '#B35806'))(201)
  breaks <- seq(-2, 2, length.out = 202)
  draw <- function() {
    layout(matrix(c(1,2,3,4,0,0,0,5), nrow = 2, byrow = TRUE),
           widths = c(3.15, 2.05, .85, 3.7), heights = c(8, 1.1))
    par(oma = c(1, .5, 5, .5), family = 'sans', xaxs = 'i', yaxs = 'i')
    setup <- function(limits) {
      par(mar = c(3.5, .2, 1, .2))
      plot.new(); plot.window(xlim = limits, ylim = c(8.5, -.5))
    }
    setup(c(0, 1))
    text(0, 1:8, labels, adj = c(0, .5), cex = .9, col = ink)
    abline(h = c(4.5, 7.5), col = '#DDE2E6', lwd = .7)
    mtext('A  Sex-adjusted associations', side = 3, line = 3, adj = 0, font = 2, cex = 1)
    mtext('Eight gene sets with FDR < 0.05 across 85 tests', side = 3, line = 1.5, adj = 0, cex = .8, col = muted)
    setup(c(-.95, .95))
    abline(v = 0, col = '#9AA4AB', lty = 2)
    abline(h = c(4.5, 7.5), col = '#DDE2E6', lwd = .7)
    segments(estimates$CI_lower, 1:8, estimates$CI_upper, 1:8, col = ink, lwd = 1.3)
    segments(estimates$CI_lower, (1:8)-.09, estimates$CI_lower, (1:8)+.09, col = ink)
    segments(estimates$CI_upper, (1:8)-.09, estimates$CI_upper, (1:8)+.09, col = ink)
    points(estimates$beta, 1:8, pch = 16, col = ink, cex = .85)
    axis(1, at = c(-.8, -.4, 0, .4, .8), cex.axis = .8, col = muted)
    mtext('SONFH - exposed control', side = 1, line = 2, cex = .8)
    mtext('GSVA score difference', side = 1, line = 3, cex = .8)
    mtext('Coefficient and 95% CI', side = 3, line = .3, cex = .85)
    setup(c(0, 1))
    abline(h = c(4.5, 7.5), col = '#DDE2E6', lwd = .7)
    fdr_text <- ifelse(estimates$FDR < .0001, formatC(estimates$FDR, format = 'e', digits = 2), sprintf('%.4f', estimates$FDR))
    text(.5, 1:8, fdr_text, cex = .8, col = ink)
    mtext('FDR', side = 3, line = .3, cex = .85)
    setup(c(.5, 40.5))
    image(x = 1:40, y = 1:8, z = t(z), col = colors, breaks = breaks, add = TRUE, useRaster = TRUE)
    rect(.5, -.05, 10.5, .3, col = control, border = NA)
    rect(10.5, -.05, 40.5, .3, col = case, border = NA)
    segments(10.5, .5, 10.5, 8.5, col = 'white', lwd = 1.5)
    abline(h = c(4.5, 7.5), col = 'white', lwd = .7)
    mtext('B  Individual GSVA profiles', side = 3, line = 3, adj = 0, font = 2, cex = 1)
    mtext('GC-exposed non-SONFH (n=10)', side = 3, line = 1.7, adj = 0, col = control, cex = .8)
    mtext('SONFH (n=30)', side = 3, line = .6, adj = 0, col = case, cex = .8)
    par(mar = c(2.5, 2.3, 0, 2.3))
    plot.new(); plot.window(xlim = c(-2,2), ylim = c(0,1))
    edges <- seq(-2,2,length.out = length(colors)+1)
    rect(edges[-length(edges)], .55, edges[-1], .95, col = colors, border = NA)
    axis(1, at = c(-2,0,2), pos = .55, cex.axis = .8, tck = -.08, col = muted)
    mtext('Gene-set score z', side = 1, line = 1.3, cex = .8)
  }
  save_plot <- function(open_device) {
    open_device(); on.exit(dev.off(), add = TRUE); draw()
  }
  save_plot(function() pdf(file.path(out, 'Figure2_R.pdf'), width = 11, height = 6.4))
  save_plot(function() png(file.path(out, 'Figure2_R.png'), width = 11, height = 6.4,
                           units = 'in', res = 300, bg = 'white'))
  tryCatch(save_plot(function() tiff(file.path(out, 'Figure2_R.tiff'), width = 11, height = 6.4,
                  units = 'in', res = 600, compression = 'lzw', bg = 'white')),
           error = function(e) warning('TIFF export failed: ', conditionMessage(e)))
  estimates$process <- labels
  write.csv(estimates, file.path(out, 'Figure2_estimates_RECALCULATED.csv'), row.names = FALSE)
  write.csv(z, file.path(out, 'Figure2_display_z_RECALCULATED.csv'))
  write.csv(meta[match(sample_order, meta$GSM), c('GSM','group','sex')],
            file.path(out, 'Figure2_sample_order.csv'), row.names = FALSE)
  invisible(z)
}

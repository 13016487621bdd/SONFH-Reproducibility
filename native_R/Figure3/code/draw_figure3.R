# Figure 3: base R renderer and archived-result checks.
# Open this file in RStudio and click Source. No additional packages required.
# This script redraws archived results; it does not refit CR2/HC3 models.

figure3_main <- function(root_override = NULL) {
  # Locate this script when sourced in RStudio or run via Rscript.
  candidates <- character()
  for (fr in sys.frames()) {
    if (!is.null(fr$ofile)) candidates <- c(candidates, as.character(fr$ofile))
  }
  argv <- commandArgs(trailingOnly = FALSE)
  cli <- sub('^--file=', '', argv[grepl('^--file=', argv)])
  candidates <- c(candidates, cli)
  roots <- c(root_override, if (length(candidates)) dirname(normalizePath(candidates, mustWork = FALSE)), getwd())
  roots <- roots[file.exists(file.path(roots, 'inputs', 'Layer2_21_groups.csv'))]
  if (!length(roots)) {
    if (!interactive()) stop('Cannot find inputs. Run from the extracted package folder.')
    message('Select inputs/Layer2_21_groups.csv in the extracted package.')
    selected <- file.choose()
    roots <- dirname(dirname(selected))
  }
  root <- roots[1]
  inp <- file.path(root, 'inputs')
  out <- file.path(root, 'results_R')
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  read_input <- function(name) read.csv(file.path(inp, name), check.names = FALSE,
                                       stringsAsFactors = FALSE, fileEncoding = 'UTF-8-BOM')
  a <- read_input('Layer2_21_groups.csv')
  a <- a[order(a$domain_number, a$functional_group, method = 'radix'), ]
  b <- read_input('Layer3_11_GSVA_contrasts.csv')
  defs <- read_input('Layer3_frozen_pairs.csv')
  stopifnot(nrow(a) == 21L, all(a$n_sets >= 2), sum(a$n_sets) == 78,
            nrow(b) == 11L, nrow(defs) == 11L,
            identical(b$Pair_ID, paste0('Pair_', 1:11)),
            identical(b$Set_A, defs[['Set A']]), identical(b$Set_B, defs[['Set B']]),
            all(is.finite(a$p)), all(is.finite(b$p)),
            max(abs(p.adjust(a$p, method = 'BH') - a$FDR)) < 1e-12,
            max(abs(p.adjust(b$p, method = 'BH') - b$FDR)) < 1e-12,
            sum(a$FDR < .05) == 0L, sum(b$FDR < .05) == 1L,
            b$Pair_ID[b$FDR < .05] == 'Pair_11',
            all(b$CI_lower < b$beta & b$beta < b$CI_upper))
  labels_a <- a$functional_group
  labels_a[labels_a == 'Angiogenesis / endothelial remodeling'] <- 'Angiogenesis / endothelium'
  labels_a[labels_a == 'Osteoclast differentiation / bone resorption'] <- 'Osteoclasts / bone resorption'
  labels_b <- c('Coagulation - fibrinolysis',
                'Endothelial proliferation\n- endothelial apoptosis',
                'Fatty acid biosynthesis\n- beta-oxidation',
                'Triglyceride storage\n- catabolism',
                'Cholesterol import - efflux',
                'Phospholipid biosynthesis\n- catabolism',
                'Ceramide biosynthesis\n- catabolism',
                'Glycolipid biosynthesis\n- catabolism',
                'Lipoprotein biosynthesis\n- catabolism',
                'Lipid-droplet formation\n- disassembly',
                'Hydrogen peroxide\nbiosynthesis - catabolism')
  ink <- '#263746'; gray <- '#74818A'; orange <- '#B5632D'
  draw <- function() {
    layout(matrix(1:6, nrow = 1), widths = c(3.05, 1.25, .65, 2.55, 1.65, .65))
    par(oma = c(3, 1, 5, 1), family = 'sans', xaxs = 'i', yaxs = 'i')
    setup <- function(n, limits, bottom = FALSE) {
      par(mar = c(3.5, .15, 1.5, .15))
      plot.new(); plot.window(xlim = limits, ylim = c(n + .5, .5))
      if (bottom) axis(1, cex.axis = .8, col = gray, col.axis = ink)
    }
    separators <- which(diff(a$domain_number) != 0) + .5
    setup(21, c(0, 1))
    abline(h = separators, col = '#E2E6E9', lwd = .6)
    text(0, 1:21, labels_a, pos = 4, offset = 0, cex = .78, col = ink)
    mtext('A  Related functional groups', side = 3, adj = 0, line = 2.8, font = 2, cex = 1)
    mtext('0 of 21 met FDR < 0.05', side = 3, adj = 0, line = 1.2, cex = .8, col = gray)
    setup(21, c(0, 1.02), TRUE)
    abline(h = separators, col = '#E2E6E9', lwd = .6)
    abline(v = .05, col = orange, lty = 3, lwd = 1.5)
    points(a$FDR, 1:21, pch = 16, cex = .7, col = ink)
    mtext('FDR', side = 1, line = 2, cex = .85)
    setup(21, c(0, 1))
    abline(h = separators, col = '#E2E6E9', lwd = .6)
    text(.5, 1:21, sprintf('%.4f', a$FDR), cex = .8, col = ink)
    mtext('FDR', side = 3, line = .25, cex = .85)
    setup(11, c(0, 1))
    colors <- ifelse(b$FDR < .05, orange, gray)
    text(0, 1:11, labels_b, pos = 4, offset = 0, cex = .8,
         col = ifelse(b$FDR < .05, orange, ink), font = ifelse(b$FDR < .05, 2, 1))
    mtext('B  Opposing-process contrasts', side = 3, adj = 0, line = 2.8, font = 2, cex = 1)
    mtext('1 of 11 met FDR < 0.05', side = 3, adj = 0, line = 1.2, cex = .8, col = gray)
    setup(11, c(-2.05, 2.9), TRUE)
    abline(v = 0, lty = 2, col = '#9CA7AE')
    segments(b$CI_lower, 1:11, b$CI_upper, 1:11, col = colors, lwd = 1.3)
    segments(b$CI_lower, (1:11)-.08, b$CI_lower, (1:11)+.08, col = colors)
    segments(b$CI_upper, (1:11)-.08, b$CI_upper, (1:11)+.08, col = colors)
    points(b$beta, 1:11, pch = ifelse(b$FDR < .05, 18, 16), col = colors, cex = .9)
    mtext('Sex-adjusted difference', side = 1, line = 2, cex = .8)
    mtext('in z(A) - z(B)', side = 1, line = 3, cex = .8)
    mtext('Coefficient and 95% CI', side = 3, line = .25, cex = .8)
    setup(11, c(0, 1))
    text(.5, 1:11, sprintf('%.4f', b$FDR), col = colors, cex = .8,
         font = ifelse(b$FDR < .05, 2, 1))
    mtext('FDR', side = 3, line = .25, cex = .85)
    mtext('Separate prespecified test families: 21 groups and 11 contrasts.',
          outer = TRUE, side = 1, line = 1, cex = .8, col = gray)
  }
  save_plot <- function(open_device) {
    open_device()
    on.exit(dev.off(), add = TRUE)
    draw()
  }
  save_plot(function() pdf(file.path(out, 'Figure3_R.pdf'), width = 12, height = 7.3))
  save_plot(function() png(file.path(out, 'Figure3_R.png'), width = 12, height = 7.3,
                           units = 'in', res = 300, bg = 'white'))
  # TIFF support depends on the local R graphics build; PNG/PDF are always attempted first.
  tryCatch(save_plot(function() tiff(file.path(out, 'Figure3_R.tiff'), width = 12,
             height = 7.3, units = 'in', res = 600, compression = 'lzw', bg = 'white')),
           error = function(e) warning('TIFF export failed: ', conditionMessage(e)))
  write.csv(a, file.path(out, 'Figure3A_values.csv'), row.names = FALSE)
  b$display_label <- gsub('\n', ' ', labels_b, fixed = TRUE)
  write.csv(b, file.path(out, 'Figure3B_values.csv'), row.names = FALSE)
  report <- c('PASS: archived-result checks completed.',
              'Layer 2: 21 groups; separate BH correction; 0 significant.',
              sprintf('Lowest Layer 2 FDR: %.16g', min(a$FDR)),
              'Layer 3: 11 contrasts; separate BH correction; only Pair_11 significant.',
              sprintf('Pair_11 FDR: %.16g', b$FDR[11]),
              'Intervals are nominal 95% CIs, not multiplicity-adjusted.',
              'No native CR2/HC3 model refitting was performed by this script.')
  writeLines(report, file.path(out, 'verification_R.txt'))
  capture.output(sessionInfo(), file = file.path(out, 'sessionInfo_R.txt'))
  message(paste(report, collapse = '\n'))
  message('Figures and checks saved to: ', normalizePath(out, winslash = '/', mustWork = TRUE))
  invisible(out)
}

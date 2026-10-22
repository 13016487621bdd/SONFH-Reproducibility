# Open this file in RStudio, then click Source.
# Scope: archived gene-level expression -> GSVA -> 11 HC3 contrasts + component models -> Figure 4.
# Frozen output files are read ONLY for verification, never as model inputs.

run_figure4 <- function() {
  candidates <- character()
  for (fr in sys.frames()) if (!is.null(fr$ofile)) candidates <- c(candidates, as.character(fr$ofile))
  args <- commandArgs(FALSE)
  candidates <- c(candidates, sub('^--file=', '', args[grepl('^--file=', args)]))
  roots <- c(if (length(candidates)) dirname(normalizePath(candidates, mustWork = FALSE)), getwd())
  roots <- roots[file.exists(file.path(roots, 'inputs', 'frozen_85.gmt'))]
  if (!length(roots)) {
    if (!interactive()) stop('Cannot locate the extracted package.')
    message('Select RUN_FIGURE4.R from the extracted package.')
    roots <- dirname(file.choose())
  }
  root <- normalizePath(roots[1], winslash = '/', mustWork = TRUE)
  stamp <- paste0(format(Sys.time(), '%Y%m%d_%H%M%S'), '_', Sys.getpid())
  out <- file.path(root, 'results', stamp)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  passed <- FALSE
  audit <- data.frame(check = character(), max_absolute_error = numeric(), tolerance = numeric())
  on.exit({
    write.csv(audit, file.path(out, 'numerical_checks.csv'), row.names = FALSE)
    capture.output(sessionInfo(), file = file.path(out, 'sessionInfo.txt'))
    writeLines(if (passed) 'PASS_EXPRESSION_TO_FIGURE4' else
      'NOT_PASSED: run stopped before completion; retain the Console error.', file.path(out, 'FINAL_STATUS.txt'))
  }, add = TRUE)
  verify <- function(label, x, y, tol = 1e-8) {
    if (length(x) != length(y) || !length(x)) stop(label, ': length mismatch.')
    error <- max(abs(as.numeric(x) - as.numeric(y)))
    audit <<- rbind(audit, data.frame(check = label, max_absolute_error = error, tolerance = tol))
    if (!is.finite(error) || error >= tol) stop(label, ': numerical mismatch ', error)
    message('PASS: ', label, ' (max error ', format(error, digits = 4), ')')
  }
  required <- c(GSVA = '2.4.9')
  for (pkg in names(required)) {
    if (!requireNamespace(pkg, quietly = TRUE)) stop('Missing package: ', pkg, '. No packages were installed automatically.')
    actual <- as.character(packageVersion(pkg))
    if (actual != required[[pkg]]) stop('Archived version required: ', pkg, ' ', required[[pkg]], '; current version: ', actual)
  }
  input <- file.path(root, 'inputs')
  ref <- file.path(root, 'reference')
  read_csv <- function(path, ...) read.csv(path, check.names = FALSE, stringsAsFactors = FALSE,
                                          fileEncoding = 'UTF-8-BOM', ...)
  files <- c(list.files(input, recursive = TRUE, full.names = TRUE),
             list.files(ref, recursive = TRUE, full.names = TRUE),
             list.files(file.path(root, 'code'), full.names = TRUE), file.path(root, 'RUN_FIGURE4.R'))
  write.csv(data.frame(file = substring(files, nchar(root) + 2), MD5 = unname(tools::md5sum(files))),
            file.path(out, 'input_and_code_MD5.csv'), row.names = FALSE)

  message('Stage 1/4: expression and frozen definitions')
  meta <- read_csv(file.path(input, '03_QC/SONFH_sample_manifest_aligned.csv'))
  meta$group <- factor(meta$group, levels = c('GC_exposed_non_SONFH', 'SONFH'))
  meta$sex <- factor(meta$sex, levels = c('Female', 'Male'))
  rownames(meta) <- meta$GSM
  expr <- as.matrix(read_csv(file.path(input, '03_QC/SONFH_QC03_gene_expression_median_log2.csv.gz'), row.names = 1))
  stopifnot(identical(dim(expr), c(18833L, 40L)), nrow(meta) == 40L,
            identical(colnames(expr), meta$GSM), !anyDuplicated(meta$GSM),
            !anyDuplicated(rownames(expr)), !anyNA(meta$group), !anyNA(meta$sex), all(is.finite(expr)),
            sum(meta$group == 'SONFH') == 30L, sum(meta$group == 'GC_exposed_non_SONFH') == 10L)
  lines <- strsplit(readLines(file.path(input, 'frozen_85.gmt'), warn = FALSE), '\t')
  sets <- setNames(lapply(lines, function(x) x[-c(1, 2)]), vapply(lines, `[`, character(1), 1))
  stopifnot(length(sets) == 85L, !anyDuplicated(names(sets)))

  message('Stage 2/4: recalculate all 85 GSVA scores; please wait')
  param <- GSVA::gsvaParam(expr, sets, kcdf = 'Gaussian', tau = 1, maxDiff = TRUE,
                           absRanking = FALSE, minSize = 1, maxSize = Inf)
  scores <- GSVA::gsva(param, verbose = FALSE)
  old_scores <- as.matrix(read_csv(file.path(input, '04_GeneSet_Scoring/SONFH_GSVA_scores_85x40.csv'), row.names = 1))
  stopifnot(setequal(rownames(scores), rownames(old_scores)), setequal(colnames(scores), colnames(old_scores)))
  scores <- scores[rownames(old_scores), meta$GSM, drop = FALSE]
  verify('All 85 x 40 GSVA scores', scores, old_scores[, meta$GSM, drop = FALSE])
  write.csv(scores, file.path(out, 'GSVA_scores_RECALCULATED.csv'))
  gsva_z <- t(scale(t(scores)))
  stopifnot(all(is.finite(gsva_z)))
  write.csv(gsva_z, file.path(out, 'GSVA_z_scores_RECALCULATED.csv'))
  saveRDS(param, file.path(out, 'GSVA_parameters_RECALCULATED.rds'))

  message('Stage 3/4: refit 11 contrasts and two component models')
  hc3 <- function(y) {
    fit <- lm(y ~ group + sex, data = meta)
    x <- model.matrix(fit)
    stopifnot(fit$rank == 3L, df.residual(fit) == 37L)
    aa <- solve(crossprod(x), t(x)); h <- rowSums(x * t(aa)); e <- residuals(fit)
    se <- sqrt(sum((aa['groupSONFH', ] * e / (1 - h))^2))
    beta <- unname(coef(fit)['groupSONFH']); df <- df.residual(fit)
    c(beta = beta, HC3_SE = se, t = beta / se, p = 2 * pt(-abs(beta / se), df),
      CI_lower = beta - qt(.975, df) * se, CI_upper = beta + qt(.975, df) * se)
  }
  pairs <- read_csv(file.path(ref, 'frozen_11_pairs.csv'))
  old_pairs <- read_csv(file.path(ref, 'Layer3_archived_all11.csv'))
  stopifnot(nrow(pairs) == 11L, identical(pairs[['Set A']], old_pairs$Set_A),
            identical(pairs[['Set B']], old_pairs$Set_B))
  contrasts <- t(vapply(seq_len(11), function(j)
    gsva_z[pairs[['Set A']][j], ] - gsva_z[pairs[['Set B']][j], ], numeric(40)))
  rownames(contrasts) <- old_pairs$Pair_ID; colnames(contrasts) <- meta$GSM
  fitvals <- as.data.frame(t(apply(contrasts, 1, hc3)))
  fitvals$FDR <- p.adjust(fitvals$p, method = 'BH')
  new_pairs <- cbind(old_pairs[, c('Pair_ID', 'Set_A', 'Set_B')], fitvals)
  for (k in c('beta', 'HC3_SE', 't', 'p', 'CI_lower', 'CI_upper', 'FDR'))
    verify(paste('All11 contrasts', k), new_pairs[[k]], old_pairs[[k]])
  stopifnot(pairs[['Set A']][11] == 'GOBP_HYDROGEN_PEROXIDE_BIOSYNTHETIC_PROCESS',
            pairs[['Set B']][11] == 'GOBP_HYDROGEN_PEROXIDE_CATABOLIC_PROCESS')
  values <- rbind(Biosynthesis = gsva_z[pairs[['Set A']][11], ],
                  Catabolism = gsva_z[pairs[['Set B']][11], ],
                  Contrast = contrasts[11, ])
  components <- as.data.frame(t(apply(values[1:2, , drop = FALSE], 1, hc3)))
  components$outcome <- rownames(components)
  old_components <- read_csv(file.path(ref, 'components_standardized.csv'))
  old_components <- old_components[match(components$outcome, old_components$outcome), ]
  stopifnot(!anyNA(old_components$outcome), nrow(old_components) == 2L)
  for (k in c('beta', 'HC3_SE', 't', 'p', 'CI_lower', 'CI_upper'))
    verify(paste('Components', k), components[[k]], old_components[[k]])
  display <- rbind(components[, c('beta','HC3_SE','t','p','CI_lower','CI_upper')],
                    new_pairs[11, c('beta','HC3_SE','t','p','CI_lower','CI_upper')])
  display$outcome <- c('Biosynthesis','Catabolism','Contrast')
  display$primary_11_test_FDR <- c(NA_real_, NA_real_, new_pairs$FDR[11])
  verify('Coefficient decomposition', display$beta[3], display$beta[1] - display$beta[2], tol = 1e-12)
  write.csv(new_pairs, file.path(out, 'Layer3_all11_RECALCULATED.csv'), row.names = FALSE)
  write.csv(values, file.path(out, 'Figure4_participant_values.csv'))
  write.csv(display, file.path(out, 'Figure4_models_RECALCULATED.csv'), row.names = FALSE)

  message('Stage 4/4: draw Figure 4 from freshly recalculated values')
  plot_env <- new.env(parent = environment())
  source(file.path(root, 'code/draw_figure4.R'), local = plot_env, encoding = 'UTF-8')
  figure_folder <- file.path(out, 'figure')
  plot_env$draw_figure4(values, display, meta, figure_folder)
  stopifnot(file.exists(file.path(figure_folder, 'Figure4_R.png')),
            file.exists(file.path(figure_folder, 'Figure4_R.pdf')))
  passed <- TRUE
  message('\nPASS_EXPRESSION_TO_FIGURE4\nGSVA, all 11 contrasts, both components and coefficient decomposition passed.')
  message('Output folder: ', out)
  message('Figure folder: ', figure_folder)
  invisible(out)
}
run_figure4()


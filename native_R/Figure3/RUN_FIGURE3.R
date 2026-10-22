# Open this file in RStudio, then click Source.
# Scope: archived gene-level expression -> GSVA -> CR2/HC3 -> Figure 3.
# Frozen output files are read ONLY for verification, never as model inputs.

run_figure3 <- function() {
  candidates <- character()
  for (fr in sys.frames()) if (!is.null(fr$ofile)) candidates <- c(candidates, as.character(fr$ofile))
  args <- commandArgs(FALSE)
  candidates <- c(candidates, sub('^--file=', '', args[grepl('^--file=', args)]))
  roots <- c(if (length(candidates)) dirname(normalizePath(candidates, mustWork = FALSE)), getwd())
  roots <- roots[file.exists(file.path(roots, 'inputs', 'frozen_85.gmt'))]
  if (!length(roots)) {
    if (!interactive()) stop('Cannot locate the extracted package.')
    message('Select RUN_FIGURE3.R from the extracted package.')
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
    writeLines(if (passed) 'PASS_EXPRESSION_TO_FIGURE3' else
      'NOT_PASSED: run stopped before completion; retain the Console error.', file.path(out, 'FINAL_STATUS.txt'))
  }, add = TRUE)
  verify <- function(label, x, y, tol = 1e-8) {
    if (length(x) != length(y) || !length(x)) stop(label, ': length mismatch.')
    error <- max(abs(as.numeric(x) - as.numeric(y)))
    audit <<- rbind(audit, data.frame(check = label, max_absolute_error = error, tolerance = tol))
    if (!is.finite(error) || error >= tol) stop(label, ': numerical mismatch ', error)
    message('PASS: ', label, ' (max error ', format(error, digits = 4), ')')
  }
  required <- c(GSVA = '2.4.9', clubSandwich = '0.7.0')
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
             list.files(file.path(root, 'code'), full.names = TRUE), file.path(root, 'RUN_FIGURE3.R'))
  write.csv(data.frame(file = substring(files, nchar(root) + 2), MD5 = unname(tools::md5sum(files))),
            file.path(out, 'input_and_code_MD5.csv'), row.names = FALSE)

  message('Stage 1/5: expression and frozen definitions')
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

  message('Stage 2/5: recalculate all 85 GSVA scores; please wait')
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

  message('Stage 3/5: refit all 21 CR2 functional-group tests')
  master <- read_csv(file.path(input, '05_Statistics/SONFH_Layer1_GSVA_85_audit.csv'))
  stopifnot(nrow(master) == 85L, setequal(master$gene_set, names(sets)))
  # Only frozen identifiers/group membership are read from master; its association values are not used.
  fkey <- paste(master[['Domain']], master[['Family / structural unit（家族/结构单元）']], sep = '.')
  family_members <- split(master$gene_set, fkey)
  family_members <- family_members[lengths(family_members) >= 2]
  stopifnot(length(family_members) == 21L, sum(lengths(family_members)) == 78L)
  write.csv(data.frame(family = rep(names(family_members), lengths(family_members)),
                        gene_set = unlist(family_members, use.names = FALSE),
                        reference_set = rep(vapply(family_members, `[`, character(1), 1), lengths(family_members))),
            file.path(out, 'Layer2_members_and_references.csv'), row.names = FALSE)
  sample_aligned <- meta
  source(file.path(root, 'code/layer2_original_function.R'), local = TRUE, encoding = 'UTF-8')
  res <- do.call(rbind, lapply(names(family_members), analyze_family_v2))
  res$FDR <- p.adjust(res$p, method = 'BH')
  old <- read_csv(file.path(input, '05_Statistics/SONFH_Layer2_21_family_omnibus_results.csv'))
  stopifnot(!anyDuplicated(res$family), setequal(res$family, old$family))
  res <- res[match(old$family, res$family), ]
  stopifnot(identical(res$method, old$method), identical(as.integer(res$n_sets), as.integer(old$n_sets)))
  for (k in c('statistic', 'df_num', 'df_denom', 'p', 'FDR')) verify(paste('Layer2', k), res[[k]], old[[k]])
  write.csv(res, file.path(out, 'Layer2_CR2_RECALCULATED.csv'), row.names = FALSE)

  message('Stage 4/5: refit all 11 HC3 contrasts')
  pairs <- read_csv(file.path(ref, 'Layer3_frozen_pairs.csv'))
  old_pairs <- read_csv(file.path(ref, 'Layer3_11_GSVA_contrasts.csv'))
  stopifnot(nrow(pairs) == 11L, identical(pairs[['Set A']], old_pairs$Set_A), identical(pairs[['Set B']], old_pairs$Set_B))
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
  contrasts <- t(vapply(seq_len(11), function(j)
    gsva_z[pairs[['Set A']][j], ] - gsva_z[pairs[['Set B']][j], ], numeric(40)))
  rownames(contrasts) <- old_pairs$Pair_ID; colnames(contrasts) <- meta$GSM
  fitvals <- as.data.frame(t(apply(contrasts, 1, hc3)))
  fitvals$FDR <- p.adjust(fitvals$p, method = 'BH')
  new_pairs <- cbind(old_pairs[, c('Pair_ID', 'Set_A', 'Set_B')], fitvals)
  for (k in c('beta', 'HC3_SE', 't', 'p', 'CI_lower', 'CI_upper', 'FDR'))
    verify(paste('Layer3', k), new_pairs[[k]], old_pairs[[k]])
  write.csv(contrasts, file.path(out, 'Layer3_participant_contrasts_RECALCULATED.csv'))
  write.csv(new_pairs, file.path(out, 'Layer3_HC3_RECALCULATED.csv'), row.names = FALSE)

  message('Stage 5/5: draw Figure 3 from freshly recalculated results')
  mapping <- read_csv(file.path(input, 'family_display_map.csv'))
  stopifnot(nrow(mapping) == 21L, setequal(mapping$family, res$family))
  mapped <- mapping[match(res$family, mapping$family), ]
  plot_a <- cbind(res, mapped[, c('domain_number', 'functional_group')])
  archived_a <- read_csv(file.path(ref, 'Layer2_21_groups.csv'))
  ix <- match(paste(archived_a$domain_number, archived_a$functional_group),
              paste(plot_a$domain_number, plot_a$functional_group))
  stopifnot(!anyNA(ix), !anyDuplicated(ix))
  for (k in c('statistic', 'p', 'FDR')) verify(paste('Figure3A mapping', k), plot_a[[k]][ix], archived_a[[k]])
  stage <- file.path(out, 'figure_build')
  dir.create(file.path(stage, 'inputs'), recursive = TRUE)
  write.csv(plot_a, file.path(stage, 'inputs/Layer2_21_groups.csv'), row.names = FALSE)
  write.csv(new_pairs, file.path(stage, 'inputs/Layer3_11_GSVA_contrasts.csv'), row.names = FALSE)
  write.csv(pairs, file.path(stage, 'inputs/Layer3_frozen_pairs.csv'), row.names = FALSE)
  plot_env <- new.env(parent = environment())
  source(file.path(root, 'code/draw_figure3.R'), local = plot_env, encoding = 'UTF-8')
  figure_folder <- plot_env$figure3_main(root_override = stage)
  stopifnot(file.exists(file.path(figure_folder, 'Figure3_R.png')),
            file.exists(file.path(figure_folder, 'Figure3_R.pdf')))
  passed <- TRUE
  message('\nPASS_EXPRESSION_TO_FIGURE3\nAll statistical comparisons passed; Figure 3 uses recalculated results.')
  message('Output folder: ', out)
  message('Figure folder: ', figure_folder)
  invisible(out)
}
run_figure3()

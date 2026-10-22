# Open this file in RStudio, then click Source.
# Scope: archived gene-level expression -> GSVA -> scoring/membership/proxy/leave-one-out analyses -> Figure 5.
# Frozen output files are read ONLY for verification, never as model inputs.

run_figure5 <- function() {
  candidates <- character()
  for (fr in sys.frames()) if (!is.null(fr$ofile)) candidates <- c(candidates, as.character(fr$ofile))
  args <- commandArgs(FALSE)
  candidates <- c(candidates, sub('^--file=', '', args[grepl('^--file=', args)]))
  roots <- c(if (length(candidates)) dirname(normalizePath(candidates, mustWork = FALSE)), getwd())
  roots <- roots[file.exists(file.path(roots, 'inputs', 'frozen_85.gmt'))]
  if (!length(roots)) {
    if (!interactive()) stop('Cannot locate the extracted package.')
    message('Select RUN_FIGURE5.R from the extracted package.')
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
    writeLines(if (passed) 'PASS_EXPRESSION_TO_FIGURE5' else
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
             list.files(file.path(root, 'code'), full.names = TRUE), file.path(root, 'RUN_FIGURE5.R'))
  write.csv(data.frame(file = substring(files, nchar(root) + 2), MD5 = unname(tools::md5sum(files))),
            file.path(out, 'input_and_code_MD5.csv'), row.names = FALSE)

  message('Stage 1/7: expression and frozen definitions')
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

  message('Stage 2/7: recalculate all 85 GSVA scores; please wait')
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

  hc3 <- function(y, dat = meta, proxy = NULL) {
    if (is.null(proxy)) fit <- lm(y ~ group + sex, data = dat) else {
      dat$proxy <- proxy
      fit <- lm(y ~ group + sex + proxy, data = dat)
    }
    x <- model.matrix(fit); aa <- solve(crossprod(x),t(x)); h <- rowSums(x*t(aa)); e <- residuals(fit)
    se <- sqrt(sum((aa['groupSONFH', ]*e/(1-h))^2)); b <- unname(coef(fit)['groupSONFH']); df <- df.residual(fit)
    c(beta=b, HC3_SE=se, t=b/se, p=2*pt(-abs(b/se),df),
      CI_lower=b-qt(.975,df)*se, CI_upper=b+qt(.975,df)*se)
  }
  compare_models <- function(label, fresh, archived) {
    stopifnot(nrow(fresh)==nrow(archived))
    for (k in intersect(c('beta','HC3_SE','t','p','CI_lower','CI_upper','FDR'),names(archived)))
      verify(paste(label,k),fresh[[k]],archived[[k]])
  }
  pairs <- read_csv(file.path(ref,'frozen_11_pairs.csv'))
  ref_gs <- read_csv(file.path(ref,'Layer3_GSVA.csv'))
  stopifnot(nrow(pairs)==11L, identical(pairs[['Set A']],ref_gs$Set_A), identical(pairs[['Set B']],ref_gs$Set_B))
  bio <- 'GOBP_HYDROGEN_PEROXIDE_BIOSYNTHETIC_PROCESS'
  catab <- 'GOBP_HYDROGEN_PEROXIDE_CATABOLIC_PROCESS'
  stopifnot(pairs[['Set A']][11]==bio,pairs[['Set B']][11]==catab)
  fit_all <- function(mat) {
    r <- as.data.frame(t(apply(mat,1,hc3)));r$FDR <- p.adjust(r$p,'BH');r
  }
  make_contrasts <- function(mat) {
    zz <- t(scale(t(mat)))
    cc <- t(vapply(seq_len(11),function(j) zz[pairs[['Set A']][j],]-zz[pairs[['Set B']][j],],numeric(40)))
    rownames(cc)<-ref_gs$Pair_ID;colnames(cc)<-meta$GSM;cc
  }
  message('Stage 3/7: ssGSEA rescoring and both complete test families')
  sp <- GSVA::ssgseaParam(expr,sets,alpha=.25,normalize=FALSE,minSize=1,maxSize=Inf)
  ss <- GSVA::gsva(sp,verbose=FALSE)
  ss <- ss[rownames(scores),meta$GSM,drop=FALSE]
  old_ss <- as.matrix(read_csv(file.path(ref,'ssGSEA_scores.csv'),row.names=1))
  verify('All85 ssGSEA scores',ss,old_ss[rownames(ss),colnames(ss)])
  cg <- make_contrasts(scores);cs <- make_contrasts(ss)
  rg <- fit_all(cg);rs <- fit_all(cs)
  compare_models('GSVA all11',rg,ref_gs)
  compare_models('ssGSEA all11',rs,read_csv(file.path(ref,'Layer3_ssGSEA.csv')))
  l1g <- fit_all(scores);l1s <- fit_all(ss)
  for (method in c('GSVA','ssGSEA')) {
    rr <- if(method=='GSVA') l1g else l1s
    oo <- read_csv(file.path(ref,paste0('Layer1_',method,'.csv')))
    stopifnot(setequal(rownames(rr),oo$gene_set))
    oo <- oo[match(rownames(rr),oo$gene_set),]
    compare_models(paste(method,'all85'),rr,oo)
  }
  scoring <- rbind(rg[11,],rs[11,]);scoring$method<-c('GSVA','ssGSEA')
  write.csv(ss,file.path(out,'ssGSEA_scores_RECALCULATED.csv'))
  write.csv(rg,file.path(out,'GSVA_all11_RECALCULATED.csv'))
  write.csv(rs,file.path(out,'ssGSEA_all11_RECALCULATED.csv'))
  write.csv(l1g,file.path(out,'GSVA_all85_RECALCULATED.csv'))
  write.csv(l1s,file.path(out,'ssGSEA_all85_RECALCULATED.csv'))

  message('Stage 4/7: recalculate six gene-membership scenarios')
  scenarios <- list(original=character(), remove_SNCA='SNCA', remove_PRDX2='PRDX2',
    remove_SNCA_PRDX2=c('SNCA','PRDX2'),remove_MPO_EPX=c('MPO','EPX'),
    remove_SNCA_PRDX2_MPO_EPX=c('SNCA','PRDX2','MPO','EPX'))
  removal_scores <- matrix(NA_real_,40,6,dimnames=list(meta$GSM,names(scenarios)))
  for (j in seq_along(scenarios)) {
    pp <- GSVA::gsvaParam(expr,list(catabolism=setdiff(sets[[catab]],scenarios[[j]])),
      kcdf='Gaussian',tau=1,maxDiff=TRUE,absRanking=FALSE,minSize=1,maxSize=Inf)
    removal_scores[,j] <- GSVA::gsva(pp,verbose=FALSE)[1,meta$GSM]
  }
  old_rem <- read_csv(file.path(ref,'gene_removal_scores.csv'))
  old_rem <- old_rem[match(meta$GSM,old_rem$GSM),]
  for (n in names(scenarios)) verify(paste('Membership score',n),removal_scores[,n],old_rem[[n]])
  removal <- as.data.frame(t(apply(removal_scores,2,hc3)));removal$scenario<-names(scenarios)
  old_rem_models <- read_csv(file.path(ref,'gene_removal_models.csv'))
  old_rem_models <- old_rem_models[match(removal$scenario,old_rem_models$scenario),]
  compare_models('Membership HC3',removal,old_rem_models)
  removal$measurable_genes <- vapply(scenarios,function(v)sum(setdiff(sets[[catab]],v)%in%rownames(expr)),integer(1))
  write.csv(removal_scores,file.path(out,'Membership_scores_RECALCULATED.csv'))
  write.csv(removal,file.path(out,'Membership_models_RECALCULATED.csv'),row.names=FALSE)

  message('Stage 5/7: expression proxies and adjusted models')
  proxy_genes <- c('ALAS2','AHSP','GYPA','SLC4A1')
  stopifnot(all(c(proxy_genes,'HBB')%in%rownames(expr)),!any(proxy_genes%in%sets[[catab]]))
  gene_z <- t(scale(t(expr[c(proxy_genes,'HBB'),,drop=FALSE])))
  proxy4 <- colMeans(gene_z[proxy_genes,,drop=FALSE]);proxy5 <- colMeans(gene_z)
  old_proxy <- read_csv(file.path(ref,'proxy_individual.csv'))
  old_proxy <- old_proxy[match(meta$GSM,old_proxy$GSM),]
  verify('Four-gene proxy',proxy4,old_proxy$E1_four_gene_proxy)
  verify('Five-gene proxy',proxy5,old_proxy$E4_five_gene_proxy_including_HBB)
  rho <- cor(proxy4,scores[catab,],method='spearman')
  old_cor <- read_csv(file.path(ref,'proxy_correlations.csv'))
  verify('Proxy-catabolism Spearman rho',rho,old_cor$rho[old_cor$subset=='all40'])
  outcomes <- list(Catabolism=scores[catab,],Biosynthesis=scores[bio,],Contrast=cg[11,])
  proxy_models <- do.call(rbind,lapply(names(outcomes),function(n) {
    v<-outcomes[[n]]
    r<-as.data.frame(rbind(hc3(v),hc3(v,proxy=proxy4),hc3(v,proxy=proxy5)))
    r$outcome<-n;r$model<-c('Sex adjusted','Four-gene proxy adjusted','Five-gene proxy sensitivity');r
  }))
  old_pm <- read_csv(file.path(ref,'proxy_models.csv'))
  idx <- match(paste(proxy_models$outcome,proxy_models$model),paste(old_pm$outcome,old_pm$model))
  stopifnot(!anyNA(idx),!anyDuplicated(idx));compare_models('Proxy models',proxy_models,old_pm[idx,])
  write.csv(proxy_models,file.path(out,'Proxy_models_RECALCULATED.csv'),row.names=FALSE)
  all_adj <- as.data.frame(t(apply(scores,1,function(y)hc3(y,proxy=proxy4))))
  all_adj$FDR <- p.adjust(all_adj$p,'BH')
  old_adj <- read_csv(file.path(ref,'all85_proxy_adjusted.csv'))
  old_adj <- old_adj[match(rownames(all_adj),old_adj$gene_set),]
  stopifnot(!anyNA(old_adj$gene_set));compare_models('All85 exploratory proxy-adjusted',all_adj,old_adj)
  write.csv(all_adj,file.path(out,'All85_exploratory_proxy_adjusted_RECALCULATED.csv'))
  proxy_data <- data.frame(GSM=meta$GSM,group=meta$group,proxy4=proxy4,proxy5=proxy5,catabolism=scores[catab,])
  write.csv(proxy_data,file.path(out,'Proxy_participant_values.csv'),row.names=FALSE)

  message('Stage 6/7: fixed-score leave-one-out checks (40 omissions)')
  loo_beta <- matrix(NA_real_,85,40,dimnames=dimnames(scores))
  loo_fdr <- matrix(NA_real_,11,40,dimnames=dimnames(cg))
  for (i in seq_len(40)) {
    dm<-meta[-i,,drop=FALSE]
    r1<-t(apply(scores[,-i,drop=FALSE],1,function(y)hc3(y,dat=dm)))
    r3<-t(apply(cg[,-i,drop=FALSE],1,function(y)hc3(y,dat=dm)))
    loo_beta[,i]<-r1[,'beta'];loo_fdr[,i]<-p.adjust(r3[,'p'],'BH')
    if(i%%10==0)message('Leave-one-out: ',i,'/40')
  }
  old_l3 <- as.matrix(read_csv(file.path(ref,'Layer3_LOO_FDR.csv'),row.names=1))
  verify('All11 x 40 LOO FDR',loo_fdr,old_l3[rownames(loo_fdr),colnames(loo_fdr)])
  loomin<-apply(loo_beta,1,min);loomax<-apply(loo_beta,1,max)
  old_l1<-read_csv(file.path(ref,'Layer1_LOO_summary.csv'));old_l1<-old_l1[match(rownames(scores),old_l1$gene_set),]
  verify('All85 LOO beta minima',loomin,old_l1$beta_LOO_min)
  verify('All85 LOO beta maxima',loomax,old_l1$beta_LOO_max)
  write.csv(loo_beta,file.path(out,'LOO_fixed_score_all85_beta.csv'))
  write.csv(loo_fdr,file.path(out,'LOO_fixed_score_all11_FDR.csv'))
  writeLines(c('Scores/standardization fixed from full 40 samples; only regressions refitted.',
      paste('Pair11 significant omissions:',sum(loo_fdr[11,]<.05),'of 40')),
      file.path(out,'LOO_scope_and_summary.txt'))

  message('Stage 7/7: draw Figure 5 from recalculated results')
  pe <- new.env(parent=environment())
  source(file.path(root,'code/draw_figure5.R'),local=pe,encoding='UTF-8')
  figdir<-file.path(out,'figure')
  pe$draw_figure5(scoring,removal,proxy_data,proxy_models,rho,figdir)
  stopifnot(file.exists(file.path(figdir,'Figure5_R.png')),file.exists(file.path(figdir,'Figure5_R.pdf')))
  passed<-TRUE
  message('\nPASS_EXPRESSION_TO_FIGURE5\nScoring, membership, proxy models and fixed-score LOO checks passed.')
  message('Output folder: ',out);message('Figure folder: ',figdir)
  invisible(out)
}
run_figure5()


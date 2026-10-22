# Standalone native R verification. NOT executed in the preparation environment.
# Run from the package root: Rscript code/verify_scoring_and_CR2.R
# Required archived version: GSVA 2.4.9; clubSandwich 0.7.0.
args <- commandArgs(trailingOnly=FALSE)
script <- sub('^--file=', '', args[grepl('^--file=', args)])
root <- if(length(script)) dirname(dirname(normalizePath(script))) else normalizePath('.')
inp <- file.path(root, 'inputs'); out <- file.path(root, 'results', 'native_R')
dir.create(out, recursive=TRUE, showWarnings=FALSE)
for(pkg in c('GSVA','clubSandwich')) if(!requireNamespace(pkg,quietly=TRUE)) stop('Install required package: ',pkg)
if(as.character(packageVersion('GSVA'))!='2.4.9') warning('GSVA version differs from frozen 2.4.9; comparisons below still required.')
if(as.character(packageVersion('clubSandwich'))!='0.7.0') warning('clubSandwich version differs from frozen 0.7.0.')
audit <- data.frame(check=character(),max_absolute_error=numeric())
verify <- function(label,x,y,tol=1e-8) {
  err <- max(abs(as.numeric(x)-as.numeric(y)))
  audit <<- rbind(audit,data.frame(check=label,max_absolute_error=err))
  if(!is.finite(err) || err>=tol) stop(label, ': discrepancy ',err)
}
meta <- read.csv(file.path(inp,'03_QC/SONFH_sample_manifest_aligned.csv'))
meta$group <- factor(meta$group,levels=c('GC_exposed_non_SONFH','SONFH'))
meta$sex <- factor(meta$sex,levels=c('Female','Male'))
rownames(meta) <- meta$GSM
expr <- as.matrix(read.csv(file.path(inp,'03_QC/SONFH_QC03_gene_expression_median_log2.csv.gz'),row.names=1,check.names=FALSE))
stopifnot(identical(colnames(expr),meta$GSM))
line <- strsplit(readLines(file.path(inp,'frozen_85.gmt')),'\t')
sets <- setNames(lapply(line,function(v)v[-c(1,2)]),vapply(line,`[`,character(1),1))
stopifnot(length(sets)==85)
gs_param <- readRDS(file.path(inp,'04_GeneSet_Scoring/SONFH_GSVA_parameters.rds'))
ss_param <- readRDS(file.path(inp,'04_GeneSet_Scoring/SONFH_ssGSEA_parameters.rds'))
# Native parameter objects preserve every scoring option from the original run.
verify('GSVA parameter matrix',gs_param@exprData,expr)
stopifnot(identical(rownames(gs_param@exprData),rownames(expr)),identical(colnames(gs_param@exprData),colnames(expr)))
for(k in names(sets)) stopifnot(setequal(gs_param@geneSets[[k]],sets[[k]]))
gs <- GSVA::gsva(gs_param,verbose=FALSE)
ss <- GSVA::gsva(ss_param,verbose=FALSE)
oldgs <- as.matrix(read.csv(file.path(inp,'04_GeneSet_Scoring/SONFH_GSVA_scores_85x40.csv'),row.names=1,check.names=FALSE))
oldss <- as.matrix(read.csv(file.path(inp,'04_GeneSet_Scoring/SONFH_ssGSEA_scores_85x40.csv'),row.names=1,check.names=FALSE))
verify('All GSVA scores',gs[rownames(oldgs),colnames(oldgs)],oldgs)
verify('All ssGSEA scores',ss[rownames(oldss),colnames(oldss)],oldss)
write.csv(gs,file.path(out,'GSVA_scores.csv'));write.csv(ss,file.path(out,'ssGSEA_scores.csv'))
gsva_z <- t(scale(t(gs)))
master <- read.csv(file.path(inp,'05_Statistics/SONFH_Layer1_GSVA_85_audit.csv'),check.names=FALSE)
fkey <- paste(master[['Domain']],master[['Family / structural unit（家族/结构单元）']],sep='.')
family_members <- split(master$gene_set,fkey)
family_members <- family_members[lengths(family_members)>=2]
stopifnot(length(family_members)==21)
sample_aligned <- meta
# Supplied function now receives all its required objects explicitly above.
source(file.path(inp,'08_Scripts/SONFH_Layer2_family_analysis.R'),local=TRUE)
res <- do.call(rbind,lapply(names(family_members),analyze_family_v2))
res$FDR <- p.adjust(res$p,method='BH')
old <- read.csv(file.path(inp,'05_Statistics/SONFH_Layer2_21_family_omnibus_results.csv'))
res <- res[match(old$family,res$family),]
stopifnot(!anyNA(res$family))
for(k in c('statistic','df_num','df_denom','p','FDR')) verify(paste('CR2 family',k),res[[k]],old[[k]])
write.csv(res,file.path(out,'layer2_CR2_recalculated.csv'),row.names=FALSE)

hc3 <- function(y) {
  fit <- lm(y~group+sex,data=meta); xx <- model.matrix(fit)
  aa <- solve(crossprod(xx),t(xx)); h <- rowSums(xx*t(aa)); e <- residuals(fit)
  se <- sqrt(sum((aa['groupSONFH',]*e/(1-h))^2)); b <- coef(fit)['groupSONFH'];df <- df.residual(fit)
  c(beta=b,HC3_SE=se,p=2*pt(-abs(b/se),df),CI_lower=b-qt(.975,df)*se,CI_upper=b+qt(.975,df)*se)
}
catname <- 'GOBP_HYDROGEN_PEROXIDE_CATABOLIC_PROCESS'
scenarios <- list(original=character(),remove_SNCA='SNCA',remove_PRDX2='PRDX2',remove_SNCA_PRDX2=c('SNCA','PRDX2'),remove_MPO_EPX=c('MPO','EPX'),remove_SNCA_PRDX2_MPO_EPX=c('SNCA','PRDX2','MPO','EPX'))
oldrem <- read.csv(file.path(inp,'09_KeyGene_Removal_Audit/participant_scores_by_scenario.csv'))
oldrem <- oldrem[match(meta$GSM,oldrem$GSM),]
for(n in names(scenarios)) {
  sp <- GSVA::gsvaParam(expr,list(catabolism=setdiff(sets[[catname]],scenarios[[n]])),kcdf='Gaussian',tau=1,maxDiff=TRUE,absRanking=FALSE,minSize=1,maxSize=Inf)
  v <- GSVA::gsva(sp,verbose=FALSE)[1,meta$GSM]
  verify(paste('Gene removal score',n),v,oldrem[[n]])
}

# Tissue pipeline starts at archived transformed, symbol-mapped SOD2-corrected matrix.
# Raw-count normalization/annotation is NOT reconstructed by this entry point.
obj <- readRDS(file.path(inp,'06_Validation/GSE345158_reproducibility_inputs.rds'))
stopifnot(identical(dim(obj$gsva_matrix_sod2),c(18573L,9L)))
pairdef <- read.csv(file.path(inp,'05_Statistics/SONFH_Layer3_11_pairs_GSVA.csv'))
bioname <- pairdef$Set_A[11]; catname <- pairdef$Set_B[11]
stopifnot(setequal(names(obj$pair11_sets),c(bioname,catname)))
exact <- function(v,case) {
  d <- mean(v[case])-mean(v[!case]);ix <- combn(seq_along(v),sum(case))
  per <- apply(ix,2,function(j)mean(v[j])-mean(v[-j]))
  c(Difference=d,Exact_P=mean(abs(per)>=abs(d)-1e-10))
}
trows <- list()
for(n in c('SOD2_corrected','Remove_DUOX1_DUOX2')) {
  ts <- obj$pair11_sets
  if(n!='SOD2_corrected') ts <- lapply(ts,setdiff,y=c('DUOX1','DUOX2'))
  tp <- GSVA::gsvaParam(obj$gsva_matrix_sod2,ts,kcdf='Gaussian',tau=1,maxDiff=TRUE,absRanking=FALSE,minSize=1,maxSize=Inf)
  tv <- GSVA::gsva(tp,verbose=FALSE);tz <- t(scale(t(tv)))
  v <- tz[bioname,]-tz[catname,]
  tm <- obj$sample_meta[colnames(tv),,drop=FALSE];is_case <- tm$group=='SONFH'
  test <- exact(v,is_case);trows[[n]] <- data.frame(Version=n,Difference=test[1],Exact_P=test[2])
  write.csv(data.frame(Sample=names(v),Contrast=v),file.path(out,paste0('tissue_',n,'_scores.csv')),row.names=FALSE)
}
tr <- do.call(rbind,trows);oldtr <- read.csv(file.path(inp,'06_Validation/GSE345158_Pair11_overlap_sensitivity.csv'))
oldtr <- oldtr[match(tr$Version,oldtr$Version),]
verify('Tissue effect (archive sub-micro discrepancy)',tr$Difference,oldtr$Difference,tol=2e-6)
verify('Tissue exact P',tr$Exact_P,oldtr$Exact_P)
write.csv(tr,file.path(out,'tissue_recalculated.csv'),row.names=FALSE)
write.csv(audit,file.path(out,'native_R_verification.csv'),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'))
cat('Native R scoring, CR2 and tissue verification PASSED.\n')

# Fixed count-to-score recipe. Core operations matched in the user-run recovery audit.
# This integrated path/launcher version was syntax-checked, not rerun in native R here.
# Run from the SONFH_JTM project root.
# Frozen annotation is an input; frozen matrices/scores are validation targets only.
local({
 run <- function() {
  root <- normalizePath('.', winslash='/', mustWork=TRUE)
  inp <- file.path(root, 'inputs', '06_Validation')
  report <- file.path(root, 'results', 'TISSUE_CHECK_REPORT.txt')
  out <- file.path(root, 'results', 'tissue_reconstructed', format(Sys.time(), '%Y%m%d_%H%M%S'))
  dir.create(out, recursive=TRUE, showWarnings=FALSE)
  log <- file(report, open='wt', encoding='UTF-8')
  sink(log, split=TRUE)
  on.exit({sink();close(log)}, add=TRUE)
  cat('TISSUE UPSTREAM AUDIT\n', format(Sys.time()), '\n')
  cat('Fixed matched recipe: DESeq2::vst; design ~1; blind=TRUE; nsub=1000; parametric fit.\n')
  cat('Frozen gene annotation is an input, not queried from a current database.\n')
  cat('Only full-matrix agreement may qualify a candidate; no P-value selection.\n\n')
  passed <- FALSE
  result <- tryCatch({
    pkgs <- c('DESeq2','readxl','GSVA','SummarizedExperiment')
    missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly=TRUE)]
    if(length(missing)) stop('Missing R packages: ', paste(missing, collapse=', '),
                             '. Send TISSUE_CHECK_REPORT.txt for the next step.')
    cat('PACKAGE VERSIONS\n')
    for(p in pkgs) cat(p, as.character(utils::packageVersion(p)), '\n')
    countfile <- file.path(inp,'GSE345158_gene_count_exp.annotation.xlsx')
    mapfile <- file.path(inp,'GSE345158_Pair11_analysis_objects.rds')
    reffile <- file.path(inp,'GSE345158_reproducibility_inputs.rds')
    scorefile <- file.path(inp,'GSE345158_Pair11_SOD2_corrected_scores.csv')
    paths <- c(countfile,mapfile,reffile,scorefile)
    stopifnot(all(file.exists(paths)))
    cat('\nINPUT MD5\n'); print(tools::md5sum(paths))
    raw <- as.data.frame(readxl::read_excel(countfile), check.names=FALSE)
    expected <- c(paste0('Case-',1:5),'Ctrl-1','Ctrl-2','Ctrl-4','Ctrl-5')
    stopifnot(nrow(raw)==20345L,ncol(raw)==10L,names(raw)[1]=='gene_id',
              identical(names(raw)[-1],expected),!anyDuplicated(raw$gene_id))
    counts <- as.matrix(raw[,-1,drop=FALSE]); rownames(counts) <- raw$gene_id
    stopifnot(all(is.finite(counts)),all(counts>=0),all(counts==round(counts)))
    storage.mode(counts) <- 'integer'
    archive <- readRDS(mapfile); ref <- readRDS(reffile)
    stopifnot(identical(rownames(counts),rownames(ref$vst_matrix)),
              identical(colnames(counts),colnames(ref$vst_matrix)))
    meta <- data.frame(group=factor(c(rep('SONFH',5),rep('Control',4)),
                                  levels=c('Control','SONFH')),row.names=expected)
    stopifnot(all(as.character(ref$sample_meta[expected,'group'])==as.character(meta$group)))
    cat('\nCOUNTS: 20345 genes x 9 samples; SONFH=5, Control=4.\n')

    # Derive the symbol mapping from the complete frozen annotation table.
    # Remove invalid symbols first, then simultaneously exclude ALL rows with
    # duplicated ENSEMBL identifiers or duplicated symbols in that valid table.
    gm <- archive$gene_mapping[,c('ENSEMBL','SYMBOL'),drop=FALSE]
    valid <- !is.na(gm$SYMBOL) & nzchar(gm$SYMBOL)
    gm <- gm[valid,,drop=FALSE]
    ambiguous_id <- duplicated(gm$ENSEMBL) | duplicated(gm$ENSEMBL,fromLast=TRUE)
    ambiguous_symbol <- duplicated(gm$SYMBOL) | duplicated(gm$SYMBOL,fromLast=TRUE)
    mapping <- gm[!ambiguous_id & !ambiguous_symbol,,drop=FALSE]
    stopifnot(nrow(mapping)==18572L,
              identical(as.character(mapping$ENSEMBL),as.character(ref$unique_mapping$ENSEMBL)),
              identical(as.character(mapping$SYMBOL),as.character(ref$unique_mapping$SYMBOL)),
              identical(as.character(mapping$ENSEMBL),as.character(archive$unique_mapping$ENSEMBL)),
              identical(as.character(mapping$SYMBOL),as.character(archive$unique_mapping$SYMBOL)))
    stopifnot(!'SOD2'%in%mapping$SYMBOL,'ENSG00000112096'%in%rownames(counts))
    write.csv(mapping,file.path(out,'symbol_mapping_reconstructed.csv'),row.names=FALSE)
    cat('ANNOTATION: 18572 one-to-one mappings reproduced exactly.\n')
    cat('Frozen correction: append ENSG00000112096 as SOD2.\n')

    dds <- DESeq2::DESeqDataSetFromMatrix(countData=counts,colData=meta,design=~group)
    dds <- DESeq2::estimateSizeFactors(dds,type='ratio')
    sf <- DESeq2::sizeFactors(dds)
    oldsf <- archive$size_factors[expected]
    sferr <- max(abs(sf-as.numeric(oldsf)))
    cat('SIZE FACTORS maximum absolute error:',format(sferr,digits=17),'\n')
    stopifnot(is.finite(sferr),sferr<1e-12)
    write.csv(data.frame(Sample=expected,recalculated=sf,archived=as.numeric(oldsf)),
              file.path(out,'size_factors.csv'),row.names=FALSE)

    # Four fixed standard recipes, defined before native execution.
    # The archived VST is consistent with a parametric transform; do not fit
    # coefficients to the reference matrix or transplant an inferred trend.
    recipes <- data.frame(
      candidate=c('vst_group_blind_FALSE','vst_intercept_blind_TRUE',
                  'fullVST_group_blind_FALSE','fullVST_intercept_blind_TRUE'),
      method=c('vst','vst','fullVST','fullVST'),
      blind=c(FALSE,TRUE,FALSE,TRUE), stringsAsFactors=FALSE)
    # Fixed after the successful 2026-10-21 recovery run; no recipe search.
    recipes <- recipes[2L,,drop=FALSE]
    rows <- list(); matrices <- list()
    cat('\nFIXED RECIPE: intercept-only vst, blind=TRUE, nsub=1000, fitType=parametric.\n')
    cat('Full-matrix tolerance: 1e-8. Each trend is estimated from counts.\n')
    for(i in seq_len(nrow(recipes))) {
      label <- recipes$candidate[i]
      cat('\nRunning candidate',i,'of',nrow(recipes),':',label,'\n')
      ans <- tryCatch(withCallingHandlers({
        work <- dds
        if(recipes$blind[i]) DESeq2::design(work) <- ~1
        transformed <- if(recipes$method[i]=='vst') {
          DESeq2::vst(work,blind=recipes$blind[i],nsub=1000,fitType='parametric')
        } else {
          DESeq2::varianceStabilizingTransformation(work,blind=recipes$blind[i],fitType='parametric')
        }
        mat <- as.matrix(SummarizedExperiment::assay(transformed))
        stopifnot(identical(rownames(mat),rownames(ref$vst_matrix)),
                  identical(colnames(mat),colnames(ref$vst_matrix)),all(is.finite(mat)))
        err <- max(abs(mat-ref$vst_matrix))
        list(matrix=mat,error=err,message='')
      },warning=function(w){cat('WARNING:',conditionMessage(w),'\n');invokeRestart('muffleWarning')},
        message=function(m){cat(conditionMessage(m),'\n');invokeRestart('muffleMessage')}),
        error=function(e)list(matrix=NULL,error=NA_real_,message=conditionMessage(e)))
      matches <- is.finite(ans$error) && ans$error<1e-8
      rows[[i]] <- data.frame(candidate=label,max_absolute_error=ans$error,
                              matches_reference=matches,error_message=ans$message)
      if(matches) matrices[[label]] <- ans$matrix
      cat('Maximum absolute error:',format(ans$error,digits=17),'; matches:',matches,'\n')
      if(nzchar(ans$message)) cat('Candidate error:',ans$message,'\n')
      write.csv(do.call(rbind,rows),file.path(out,'VST_candidate_comparison.csv'),row.names=FALSE)
    }
    candidates <- do.call(rbind,rows)
    cat('\nCANDIDATE COMPARISON\n');print(candidates,row.names=FALSE)
    if(!length(matrices)) stop('No fixed VST recipe matched. Do not replace frozen data. Send this report.')
    chosen <- names(matrices)[1]
    cat('Matching recipes:',paste(names(matrices),collapse='; '),'\n')
    cat('Deterministically using the first matching recipe:',chosen,'\n')
    cat('Matrix agreement identifies an equivalent reconstruction, not unique historical settings.\n')
    vst_new <- matrices[[chosen]]
    mapped <- vst_new[match(mapping$ENSEMBL,rownames(vst_new)),,drop=FALSE]
    rownames(mapped) <- mapping$SYMBOL
    sod2 <- vst_new['ENSG00000112096',,drop=FALSE];rownames(sod2) <- 'SOD2'
    corrected <- rbind(mapped,sod2)
    stopifnot(identical(dim(corrected),c(18573L,9L)),
              identical(rownames(corrected),rownames(ref$gsva_matrix_sod2)),
              identical(colnames(corrected),colnames(ref$gsva_matrix_sod2)))
    maperr <- max(abs(corrected-ref$gsva_matrix_sod2))
    cat('SYMBOL-MAPPED MATRIX maximum error:',format(maperr,digits=17),'\n')
    stopifnot(is.finite(maperr),maperr<1e-8)

    bio <- 'GOBP_HYDROGEN_PEROXIDE_BIOSYNTHETIC_PROCESS'
    catabolism <- 'GOBP_HYDROGEN_PEROXIDE_CATABOLIC_PROCESS'
    stopifnot(setequal(names(ref$pair11_sets),c(bio,catabolism)))
    par <- GSVA::gsvaParam(corrected,ref$pair11_sets,kcdf='Gaussian',tau=1,
                           maxDiff=TRUE,absRanking=FALSE,minSize=1,maxSize=Inf)
    scores <- GSVA::gsva(par,verbose=FALSE)
    z <- t(scale(t(scores)))
    contrast <- z[bio,expected]-z[catabolism,expected]
    old <- read.csv(scorefile,check.names=FALSE)
    stopifnot(!anyDuplicated(old$Sample),setequal(old$Sample,expected))
    old <- old[match(expected,old$Sample),]
    zerr <- max(abs(c(z[bio,expected]-old$Biosynthesis_Z,
                      z[catabolism,expected]-old$Catabolism_Z,contrast-old$Pair11)))
    cat('CORRECTED Z SCORES AND PAIR11 maximum error:',format(zerr,digits=17),'\n')
    stopifnot(is.finite(zerr),zerr<1e-7)
    case <- meta$group=='SONFH'
    effect <- mean(contrast[case])-mean(contrast[!case])
    permutations <- combn(seq_along(contrast),sum(case))
    perm_effect <- apply(permutations,2,function(j)mean(contrast[j])-mean(contrast[-j]))
    exact_p <- mean(abs(perm_effect)>=abs(effect)-1e-10)
    stopifnot(abs(effect-1.293719258568)<1e-7,abs(exact_p-5/126)<1e-12)
    cat('PAIR11 EFFECT:',format(effect,digits=17),'\n')
    cat('TWO-SIDED EXACT P:',format(exact_p,digits=17),' (5/126)\n')
    saveRDS(list(vst_matrix=vst_new,gsva_matrix_sod2=corrected,scores=scores,z=z,
                 contrast=contrast,mapping=mapping,recipe=recipes[recipes$candidate==chosen,],
                 all_matching_recipes=names(matrices),effect=effect,exact_p=exact_p),
            file.path(out,'tissue_upstream_reconstructed.rds'))
    write.csv(data.frame(Sample=expected,Biosynthesis_Z=z[bio,expected],
                         Catabolism_Z=z[catabolism,expected],Pair11=contrast),
              file.path(out,'corrected_scores_recalculated.csv'),row.names=FALSE)
    passed <- TRUE
    invisible(TRUE)
  },error=function(e){cat('\nAUDIT STOPPED:',conditionMessage(e),'\n');invisible(FALSE)})
  cat('\nSESSION INFO\n');print(sessionInfo())
  cat('\nFINAL STATUS:',if(passed)'PASS_RECONSTRUCTED_COUNTS_TO_PAIR11' else 'INCOMPLETE_SEND_REPORT','\n')
  cat('Scope: supplied gene counts through VST, frozen symbol mapping and corrected Pair11.\n')
  cat('FASTQ alignment/count generation and fresh annotation database lookup are outside this scope.\n')
  cat('Report:',report,'\n')
  cat('Detailed outputs:',out,'\n')
  if(!passed) message('Please upload TISSUE_CHECK_REPORT.txt. No frozen input was changed.')
 }
 run()
})

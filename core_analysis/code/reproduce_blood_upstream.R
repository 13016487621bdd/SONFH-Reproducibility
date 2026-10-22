# Base-R reproduction from deposited probe expression and platform annotation.
# Run with working directory set to the SONFH_JTM project root:
# source("code/reproduce_blood_upstream.R", encoding="UTF-8")
# Original frozen matrix is used only for final validation and is never overwritten.
local({
  root <- normalizePath('.',winslash='/',mustWork=TRUE)
  inp <- file.path(root,'inputs','01_Blood_Source')
  out <- file.path(root,'results','blood_upstream_R')
  one <- function(pattern) {
    p <- list.files(inp,pattern=pattern,full.names=TRUE,ignore.case=TRUE)
    if(length(p)!=1L) stop('Expected exactly one source file matching ',pattern,' in ',inp)
    p
  }
  series <- one('^GSE123568_series_matrix.*\\.gz$')
  platform <- one('^GPL15207.*\\.zip$')
  reference <- file.path(root,'inputs','03_QC','SONFH_QC03_gene_expression_median_log2.csv.gz')
  stopifnot(file.exists(reference))
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  con <- gzfile(series,open='rt',encoding='UTF-8')
  found <- FALSE
  repeat {
    line <- readLines(con,n=1L,warn=FALSE)
    if(!length(line))break
    if(startsWith(line,'!series_matrix_table_begin')){found<-TRUE;break}
  }
  if(!found){close(con);stop('Series matrix start marker missing')}
  probes <- as.matrix(read.delim(con,header=TRUE,row.names=1,check.names=FALSE,comment.char='!',colClasses=c('character',rep('numeric',40))))
  close(con)
  stopifnot(identical(dim(probes),c(49293L,40L)),!anyDuplicated(rownames(probes)),all(is.finite(probes)))
  members <- unzip(platform,list=TRUE)$Name
  members <- members[grepl('(^|/)GPL15207[^/]*\\.txt$',members,ignore.case=TRUE)]
  stopifnot(length(members)==1L)
  con <- unz(platform,members[1],open='r',encoding='UTF-8')
  skip <- 0L;header <- NULL
  repeat {
    line <- readLines(con,n=1L,warn=FALSE)
    if(!length(line))break
    if(startsWith(line,'ID\t')){header<-strsplit(line,'\t',fixed=TRUE)[[1]];break}
    skip <- skip+1L
  }
  close(con)
  stopifnot(!is.null(header),all(c('ID','Gene Symbol')%in%header))
  classes <- rep('NULL',length(header));classes[header%in%c('ID','Gene Symbol','Annotation Date')]<-'character'
  con <- unz(platform,members[1],open='r',encoding='UTF-8')
  gpl <- read.delim(con,skip=skip,header=TRUE,check.names=FALSE,colClasses=classes,quote='',comment.char='',na.strings=c('NA',''),fill=TRUE)
  close(con)
  stopifnot(!anyDuplicated(gpl$ID))
  symbol <- trimws(gpl[['Gene Symbol']][match(rownames(probes),gpl$ID)])
  symbol[is.na(symbol)] <- ''
  missing <- symbol%in%c('','---')
  ambiguous <- !missing & grepl('///',symbol,fixed=TRUE)
  unambiguous <- !missing & !ambiguous
  nonstandard <- unambiguous & grepl('/',symbol,fixed=TRUE)
  keep <- unambiguous & !nonstandard
  stopifnot(sum(missing)==430,sum(ambiguous)==1994,sum(unambiguous)==46869,sum(nonstandard)==2,sum(keep)==46867)
  stopifnot(setequal(rownames(probes)[nonstandard],c('11761930_at','11762411_at')))
  x <- probes[keep,,drop=FALSE]
  gene <- symbol[keep]
  ordered_genes <- sort(unique(gene),method='radix')
  groups <- split(seq_along(gene),factor(gene,levels=ordered_genes))
  rebuilt <- t(vapply(groups,function(j) {
    if(length(j)==1L) as.numeric(x[j,]) else apply(x[j,,drop=FALSE],2,median)
  },numeric(ncol(x))))
  rownames(rebuilt)<-ordered_genes;colnames(rebuilt)<-colnames(probes)
  frozen <- as.matrix(read.csv(reference,row.names=1,check.names=FALSE))
  stopifnot(identical(dim(rebuilt),c(18833L,40L)),identical(rownames(rebuilt),rownames(frozen)),identical(colnames(rebuilt),colnames(frozen)))
  err <- abs(rebuilt-frozen);maximum_error <- max(err)
  stopifnot(all(is.finite(err)),maximum_error<1e-12)
  write.csv(data.frame(ID_REF=rownames(probes)[unambiguous],Gene_Symbol=symbol[unambiguous]),file.path(out,'QC02_probe_mapping_RECONSTRUCTED.csv'),row.names=FALSE)
  decision <- rep('RETAIN_UNIQUE_GENE',length(symbol));decision[missing]<-'EXCLUDE_UNANNOTATED';decision[ambiguous]<-'EXCLUDE_AMBIGUOUS_MULTIGENE';decision[nonstandard]<-'EXCLUDE_NONSTANDARD_SLASH_SYMBOL'
  con<-gzfile(file.path(out,'QC03_probe_decisions_RECONSTRUCTED.csv.gz'),'wt')
  write.csv(data.frame(ID_REF=rownames(probes),Gene_Symbol=symbol,decision=decision),con,row.names=FALSE);close(con)
  con<-gzfile(file.path(out,'QC03_gene_expression_REPRODUCED.csv.gz'),'wt');write.csv(rebuilt,con);close(con)
  write.csv(data.frame(sample=colnames(rebuilt),maximum_absolute_error=apply(err,2,max)),file.path(out,'samplewise_matrix_comparison.csv'),row.names=FALSE)
  paths <- c(series,platform,reference)
  write.csv(data.frame(filename=basename(paths),MD5=unname(tools::md5sum(paths))),file.path(out,'source_MD5.csv'),row.names=FALSE)
  report <- c('BLOOD UPSTREAM REPRODUCTION PASSED',
    'Scope: deposited RMA probe matrix to gene-level matrix; not CEL-to-RMA reproduction.',
    paste('Source probes:',nrow(probes)),paste('Samples:',ncol(probes)),paste('Platform rows:',nrow(gpl)),
    paste('Unannotated:',sum(missing)),paste('Ambiguous:',sum(ambiguous)),paste('Nonstandard excluded:',sum(nonstandard)),
    paste('Retained probes:',sum(keep)),paste('Final genes:',nrow(rebuilt)),paste('Expression values compared:',length(err)),
    'Gene identifiers and order: EXACT','Sample identifiers and order: EXACT',
    paste('Maximum absolute difference:',format(maximum_error,digits=17)),
    'Tolerance: 1e-12',
    'QC02/QC03 audit tables reconstructed from sources; no historical-table comparison claimed.')
  writeLines(report,file.path(out,'blood_upstream_report.txt'),useBytes=TRUE)
  capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'))
  cat(paste(report,collapse='\n'),'\n')
})

# Counts -> fixed reconstructed VST recipe -> frozen mapping -> GSVA -> Figure 6.
# Open in RStudio and click Source. No historical-recipe search or P-value selection.
run_figure6 <- function() {
  candidates<-character()
  for(fr in sys.frames())if(!is.null(fr$ofile))candidates<-c(candidates,as.character(fr$ofile))
  args<-commandArgs(FALSE);candidates<-c(candidates,sub('^--file=','',args[grepl('^--file=',args)]))
  roots<-c(if(length(candidates))dirname(normalizePath(candidates,mustWork=FALSE)),getwd())
  roots<-roots[file.exists(file.path(roots,'inputs','GSE345158_reproducibility_inputs.rds'))]
  if(!length(roots)) {
    if(!interactive())stop('Cannot locate extracted package.')
    message('Select RUN_FIGURE6.R from the extracted package.');roots<-dirname(file.choose())
  }
  root<-normalizePath(roots[1],winslash='/',mustWork=TRUE);inp<-file.path(root,'inputs')
  out<-file.path(root,'results',paste0(format(Sys.time(),'%Y%m%d_%H%M%S'),'_',Sys.getpid()))
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  passed<-FALSE;audit<-data.frame(check=character(),max_absolute_error=numeric(),tolerance=numeric())
  on.exit({
    write.csv(audit,file.path(out,'numerical_checks.csv'),row.names=FALSE)
    capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'))
    writeLines(if(passed)'PASS_COUNTS_TO_FIGURE6' else 'NOT_PASSED: see Console error.',file.path(out,'FINAL_STATUS.txt'))
  },add=TRUE)
  verify<-function(label,x,y,tol=1e-8) {
    if(length(x)!=length(y)||!length(x))stop(label,': length mismatch')
    err<-max(abs(as.numeric(x)-as.numeric(y)))
    audit<<-rbind(audit,data.frame(check=label,max_absolute_error=err,tolerance=tol))
    if(!is.finite(err)||err>=tol)stop(label,': numerical discrepancy ',err)
    message('PASS: ',label,' (max error ',format(err,digits=4),')')
  }
  for(pkg in c('DESeq2','readxl','GSVA','SummarizedExperiment'))
    if(!requireNamespace(pkg,quietly=TRUE))stop('Missing package: ',pkg,'. No packages installed automatically.')
  for(pkg in c('DESeq2','GSVA')) {
    want<-if(pkg=='DESeq2')'1.50.2' else '2.4.9'
    if(as.character(packageVersion(pkg))!=want)stop('Archived version required: ',pkg,' ',want)
  }
  files<-c(list.files(inp,full.names=TRUE),file.path(root,'RUN_FIGURE6.R'),list.files(file.path(root,'code'),full.names=TRUE))
  write.csv(data.frame(file=substring(files,nchar(root)+2),MD5=unname(tools::md5sum(files))),
    file.path(out,'input_and_code_MD5.csv'),row.names=FALSE)
  read_csv<-function(name,...)read.csv(file.path(inp,name),stringsAsFactors=FALSE,check.names=FALSE,...)
  message('Stage 1/5: raw counts and frozen mapping')
  raw<-as.data.frame(readxl::read_excel(file.path(inp,'GSE345158_gene_count_exp.annotation.xlsx')),check.names=FALSE)
  expected<-c(paste0('Case-',1:5),'Ctrl-1','Ctrl-2','Ctrl-4','Ctrl-5')
  stopifnot(nrow(raw)==20345L,ncol(raw)==10L,names(raw)[1]=='gene_id',
            identical(names(raw)[-1],expected),!anyDuplicated(raw$gene_id))
  counts<-as.matrix(raw[,-1,drop=FALSE]);rownames(counts)<-raw$gene_id
  stopifnot(all(is.finite(counts)),all(counts>=0),all(counts==round(counts)))
  storage.mode(counts)<-'integer'
  archive<-readRDS(file.path(inp,'GSE345158_Pair11_analysis_objects.rds'))
  ref<-readRDS(file.path(inp,'GSE345158_reproducibility_inputs.rds'))
  stopifnot(identical(rownames(counts),rownames(ref$vst_matrix)),identical(colnames(counts),colnames(ref$vst_matrix)))
  meta<-data.frame(group=factor(c(rep('SONFH',5),rep('Control',4)),levels=c('Control','SONFH')),row.names=expected)
  stopifnot(all(as.character(ref$sample_meta[expected,'group'])==as.character(meta$group)))
  gm<-archive$gene_mapping[,c('ENSEMBL','SYMBOL'),drop=FALSE]
  gm<-gm[!is.na(gm$SYMBOL)&nzchar(gm$SYMBOL),,drop=FALSE]
  ambiguous_id<-duplicated(gm$ENSEMBL)|duplicated(gm$ENSEMBL,fromLast=TRUE)
  ambiguous_symbol<-duplicated(gm$SYMBOL)|duplicated(gm$SYMBOL,fromLast=TRUE)
  mapping<-gm[!ambiguous_id&!ambiguous_symbol,,drop=FALSE]
  stopifnot(nrow(mapping)==18572L,
    identical(as.character(mapping$ENSEMBL),as.character(ref$unique_mapping$ENSEMBL)),
    identical(as.character(mapping$SYMBOL),as.character(ref$unique_mapping$SYMBOL)),
    !'SOD2'%in%mapping$SYMBOL,'ENSG00000112096'%in%rownames(counts))
  write.csv(mapping,file.path(out,'symbol_mapping_RECONSTRUCTED.csv'),row.names=FALSE)

  message('Stage 2/5: fixed intercept-only VST recipe and SOD2 correction')
  dds<-DESeq2::DESeqDataSetFromMatrix(countData=counts,colData=meta,design=~1)
  dds<-DESeq2::estimateSizeFactors(dds,type='ratio')
  sf<-DESeq2::sizeFactors(dds)
  verify('Size factors',sf,archive$size_factors[expected],tol=1e-12)
  vst<-DESeq2::vst(dds,blind=TRUE,nsub=1000,fitType='parametric')
  mat<-as.matrix(SummarizedExperiment::assay(vst))
  stopifnot(identical(dimnames(mat),dimnames(ref$vst_matrix)))
  verify('Full VST matrix',mat,ref$vst_matrix)
  mapped<-mat[match(mapping$ENSEMBL,rownames(mat)),,drop=FALSE];rownames(mapped)<-mapping$SYMBOL
  sod2<-mat['ENSG00000112096',,drop=FALSE];rownames(sod2)<-'SOD2'
  corrected<-rbind(mapped,sod2)
  stopifnot(identical(dim(corrected),c(18573L,9L)),identical(dimnames(corrected),dimnames(ref$gsva_matrix_sod2)))
  verify('SOD2-corrected symbol matrix',corrected,ref$gsva_matrix_sod2)
  saveRDS(list(vst_matrix=mat,gsva_matrix_sod2=corrected,mapping=mapping,size_factors=sf,
    recipe='design ~1; ratio size factors; vst blind=TRUE nsub=1000 fitType=parametric'),
    file.path(out,'upstream_RECONSTRUCTED.rds'))

  message('Stage 3/5: GSVA, standardization and all 126 label permutations')
  bio<-'GOBP_HYDROGEN_PEROXIDE_BIOSYNTHETIC_PROCESS';catab<-'GOBP_HYDROGEN_PEROXIDE_CATABOLIC_PROCESS'
  stopifnot(setequal(names(ref$pair11_sets),c(bio,catab)),
    all(c('DUOX1','DUOX2')%in%ref$pair11_sets[[bio]]),all(c('DUOX1','DUOX2')%in%ref$pair11_sets[[catab]]))
  case<-meta$group=='SONFH';assignments<-combn(seq_len(9),sum(case));stopifnot(ncol(assignments)==126L)
  versions<-c('SOD2_corrected','Remove_DUOX1_DUOX2');rows<-list();participants<-list();perms<-list();zs<-list()
  for(version in versions) {
    sets<-ref$pair11_sets
    if(version!='SOD2_corrected')sets<-lapply(sets,setdiff,y=c('DUOX1','DUOX2'))
    pp<-GSVA::gsvaParam(corrected,sets,kcdf='Gaussian',tau=1,maxDiff=TRUE,absRanking=FALSE,minSize=1,maxSize=Inf)
    score<-GSVA::gsva(pp,verbose=FALSE);z<-t(scale(t(score)));zs[[version]]<-z
    v<-z[bio,expected]-z[catab,expected]
    difference<-mean(v[case])-mean(v[!case])
    perm<-apply(assignments,2,function(j)mean(v[j])-mean(v[-j]))
    hits<-sum(abs(perm)>=abs(difference)-1e-10);p<-hits/length(perm)
    rows[[version]]<-data.frame(Version=version,Difference=difference,Exact_P=p,Extreme_assignments=hits,Total_assignments=126)
    participants[[version]]<-data.frame(Version=version,Sample=expected,group=as.character(meta$group),
      Biosynthesis_Z=z[bio,expected],Catabolism_Z=z[catab,expected],Contrast=as.numeric(v))
    perms[[version]]<-data.frame(Version=version,assignment=seq_along(perm),difference=perm)
  }
  summary<-do.call(rbind,rows);data<-do.call(rbind,participants)
  old<-read_csv('GSE345158_Pair11_SOD2_corrected_scores.csv');old<-old[match(expected,old$Sample),]
  verify('Original component z and contrast',as.matrix(data[data$Version==versions[1],c('Biosynthesis_Z','Catabolism_Z','Contrast')]),
    as.matrix(old[,c('Biosynthesis_Z','Catabolism_Z','Pair11')]),tol=1e-7)
  oldsum<-read_csv('GSE345158_Pair11_overlap_sensitivity.csv');oldsum<-oldsum[match(summary$Version,oldsum$Version),]
  # Original archived summary contains rounded values; participant-level scores above use stricter tolerance.
  verify('Both tissue effects',summary$Difference,oldsum$Difference,tol=2e-6)
  verify('Both exact P values',summary$Exact_P,oldsum$Exact_P,tol=1e-8)
  stopifnot(identical(as.integer(summary$Extreme_assignments),c(5L,7L)))
  write.csv(summary,file.path(out,'Tissue_effects_RECALCULATED.csv'),row.names=FALSE)
  write.csv(data,file.path(out,'Tissue_participant_scores_RECALCULATED.csv'),row.names=FALSE)
  write.csv(do.call(rbind,perms),file.path(out,'All126_permutations_by_version.csv'),row.names=FALSE)

  message('Stage 4/5: fixed-score tissue leave-one-out differences')
  loo<-do.call(rbind,lapply(versions,function(version){
    v<-data$Contrast[data$Version==version]
    data.frame(Version=version,removed=expected,difference=vapply(seq_len(9),function(i){
      keep<-seq_len(9)!=i;mean(v[keep&case])-mean(v[keep&!case])},numeric(1)))
  }))
  oldloo<-read_csv('tissue_LOO_reference.csv');oldloo<-oldloo[match(expected,oldloo$removed),]
  verify('Original tissue LOO effects',loo$difference[loo$Version==versions[1]],oldloo$difference)
  write.csv(loo,file.path(out,'Tissue_fixed_score_LOO.csv'),row.names=FALSE)
  message('Stage 5/5: draw Figure 6 from recalculated tissue scores')
  pe<-new.env(parent=environment());source(file.path(root,'code/draw_figure6.R'),local=pe,encoding='UTF-8')
  figdir<-file.path(out,'figure');pe$draw_figure6(data,summary,figdir)
  stopifnot(file.exists(file.path(figdir,'Figure6_R.png')),file.exists(file.path(figdir,'Figure6_R.pdf')))
  passed<-TRUE
  message('\nPASS_COUNTS_TO_FIGURE6\nFull VST matrix, symbol mapping, both tissue analyses and plotting completed.')
  message('Equivalent reconstructed workflow; not recovered unique historical source code.')
  message('Output folder: ',out);message('Figure folder: ',figdir)
  invisible(out)
}
run_figure6()

# Figure 1: study design, base R + grid only. Open in RStudio and click Source.
# Checks archived counts/definitions; does not refit statistical models.
run_figure1 <- function() {
  candidates <- character()
  for (fr in sys.frames()) if (!is.null(fr$ofile)) candidates <- c(candidates, as.character(fr$ofile))
  args <- commandArgs(FALSE)
  candidates <- c(candidates, sub('^--file=', '', args[grepl('^--file=', args)]))
  roots <- c(if (length(candidates)) dirname(normalizePath(candidates, mustWork = FALSE)), getwd())
  roots <- roots[file.exists(file.path(roots, 'inputs', 'frozen_85.gmt'))]
  if (!length(roots)) {
    if (!interactive()) stop('Cannot locate the extracted package.')
    message('Select RUN_FIGURE1.R from the extracted package.')
    roots <- dirname(file.choose())
  }
  root <- normalizePath(roots[1], winslash = '/', mustWork = TRUE)
  inp <- file.path(root, 'inputs')
  out <- file.path(root, 'results', paste0(format(Sys.time(), '%Y%m%d_%H%M%S'), '_', Sys.getpid()))
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  passed <- FALSE
  on.exit({
    capture.output(sessionInfo(), file = file.path(out, 'sessionInfo.txt'))
    writeLines(if (passed) 'PASS_FIGURE1_DEFINITIONS_AND_RENDER' else 'NOT_PASSED: see Console error.',
               file.path(out, 'FINAL_STATUS.txt'))
  }, add = TRUE)
  read_csv <- function(path, ...) read.csv(path, stringsAsFactors = FALSE, check.names = FALSE,
                                          fileEncoding = 'UTF-8-BOM', ...)
  meta <- read_csv(file.path(inp, 'sample_manifest.csv'))
  expr <- as.matrix(read_csv(file.path(inp, 'gene_expression.csv.gz'), row.names = 1))
  defs <- read_csv(file.path(inp, 'frozen_85_definitions.csv'))
  pairs <- read_csv(file.path(inp, 'frozen_11_pairs.csv'))
  tissue <- read_csv(file.path(inp, 'tissue_sample_manifest.csv'))
  gmt <- strsplit(readLines(file.path(inp, 'frozen_85.gmt'), warn = FALSE), '\t')
  set_names <- vapply(gmt, `[`, character(1), 1)
  domain <- defs[['Domain']]
  group <- paste(domain, defs[['Family / structural unit（家族/结构单元）']], sep = '.')
  group_sizes <- table(group)
  n <- nrow(meta); cases <- sum(meta$group == 'SONFH'); controls <- sum(meta$group == 'GC_exposed_non_SONFH')
  stopifnot(identical(dim(expr), c(18833L, 40L)), identical(colnames(expr), meta$GSM),
            !anyDuplicated(meta$GSM), !anyDuplicated(rownames(expr)), all(is.finite(expr)),
            n == 40L, cases == 30L, controls == 10L, length(unique(domain)) == 6L,
            nrow(defs) == 85L, !anyDuplicated(defs$gene_set), length(gmt) == 85L,
            !anyDuplicated(set_names), setequal(defs$gene_set, set_names),
            length(group_sizes) == 28L, sum(group_sizes >= 2) == 21L,
            sum(group_sizes[group_sizes >= 2]) == 78L, nrow(pairs) == 11L,
            all(pairs[['Set A']] %in% set_names), all(pairs[['Set B']] %in% set_names),
            !anyDuplicated(paste(pairs[['Set A']], pairs[['Set B']])),
            nrow(tissue) == 9L, !anyDuplicated(tissue$Sample),
            sum(tissue$Group == 'SONFH') == 5L, sum(tissue$Group == 'Control') == 4L)
  numbers <- data.frame(item = c('Blood samples','SONFH','GC-exposed non-SONFH','Genes',
     'Domains','Gene sets','Functional groups','Groups with at least 2 sets','Sets in multi-set groups',
     'Opposing-process pairs','Tissue SONFH','Tissue fracture controls'),
     value = c(n,cases,controls,nrow(expr),length(unique(domain)),nrow(defs),length(group_sizes),
               sum(group_sizes >= 2),sum(group_sizes[group_sizes >= 2]),nrow(pairs),5,4))
  write.csv(numbers, file.path(out, 'Figure1_verified_counts.csv'), row.names = FALSE)
  files <- c(list.files(inp, full.names = TRUE), file.path(root, 'RUN_FIGURE1.R'))
  write.csv(data.frame(file = substring(files, nchar(root)+2), MD5 = unname(tools::md5sum(files))),
            file.path(out, 'input_and_code_MD5.csv'), row.names = FALSE)
  ink <- '#223544'; muted <- '#526371'; line <- '#A9B6BD'
  draw <- function() {
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(xscale = c(0,9), yscale = c(0,8.1)))
    on.exit(grid::popViewport(), add = TRUE)
    label <- function(x,y,s,size=10,bold=FALSE,color=ink,center=FALSE) {
      grid::grid.text(s,x=grid::unit(x,'native'),y=grid::unit(y,'native'),
        just=c(if(center) 'centre' else 'left','top'),
        gp=grid::gpar(fontfamily='sans',fontsize=size,fontface=if(bold) 'bold' else 'plain',
                      col=color,lineheight=1.25))
    }
    box <- function(x,y,w,h,title,body,fill='#EAF1F5',accent='#315D77',fs=10,lty=1) {
      grid::grid.roundrect(x=grid::unit(x+w/2,'native'),y=grid::unit(y+h/2,'native'),
        width=grid::unit(w,'native'),height=grid::unit(h,'native'),r=grid::unit(.06,'in'),
        gp=grid::gpar(fill=fill,col='#687C89',lwd=.8,lty=lty))
      label(x+.16,y+h-.16,title,10.5,TRUE,accent)
      label(x+.16,y+h-.53,body,fs)
    }
    path <- function(x,y) grid::grid.lines(x=grid::unit(x,'native'),y=grid::unit(y,'native'),gp=grid::gpar(col=line,lwd=1.1))
    arrow <- function(x1,y1,x2,y2) grid::grid.lines(x=grid::unit(c(x1,x2),'native'),
       y=grid::unit(c(y1,y2),'native'),arrow=grid::arrow(length=grid::unit(.08,'in'),type='closed'),
       gp=grid::gpar(col=line,fill=line,lwd=1.1))
    label(.25,7.98,'A',14,TRUE)
    label(.62,7.98,'Study inputs and functional definitions',12,TRUE)
    box(.25,5.90,4.08,1.65,'Discovery cohort | GSE123568',
        sprintf('%d SONFH + %d glucocorticoid-exposed controls',cases,controls),fs=9.8)
    label(.41,6.73,'Deposited RMA matrix',9.2)
    label(.41,6.45,'GPL15207 annotation',9.2)
    path(c(2.10,2.40,2.40),c(6.64,6.64,6.50))
    path(c(2.10,2.40,2.40),c(6.36,6.36,6.50))
    arrow(2.40,6.50,2.61,6.50)
    label(2.68,6.72,'Median across\nmapped probes',9.2)
    label(.41,6.12,sprintf('Output: %s genes x %d samples',format(nrow(expr),big.mark=',',trim=TRUE),n),9.6,TRUE)
    box(4.58,6.12,4.17,1.43,'Literature-guided gene-set panel',
        paste0(length(unique(domain)),' predefined pathophysiological domains\n',nrow(defs),
               ' frozen gene sets in ',length(group_sizes),' functional groups\n',sum(group_sizes>=2),
               ' groups with >=2 sets; ',nrow(pairs),' opposing-process pairs'),fs=9.8)
    arrow(2.29,5.88,2.29,5.66);arrow(6.66,6.10,6.66,5.66)
    path(c(2.29,6.66),c(5.66,5.66));arrow(4.5,5.66,4.5,5.45)
    label(4.5,5.38,'GSVA scoring in all 40 discovery samples',11,TRUE,center=TRUE)
    label(.25,5.03,'B',14,TRUE)
    label(.62,5.03,'Three parallel questions specified before association testing',11.7,TRUE)
    path(c(1.63,7.37),c(4.68,4.68))
    for (x in c(1.63,4.5,7.37)) arrow(x,4.68,x,4.44)
    box(.25,2.49,2.73,1.91,'Layer 1 | Association',
        'Which functions are associated\nwith established SONFH?\n\n85 sets; sex-adjusted HC3\nBH across 85 tests',fs=9.5)
    box(3.14,2.49,2.73,1.91,'Layer 2 | Effect differences',
        'Do related functional components\nhave distinguishable effects?\n\n21 groups; standardized scores\nSex-adjusted CR2; BH across 21',
        fill='#EAF3F0',accent='#286C60',fs=9.3,lty=2)
    box(6.02,2.49,2.73,1.91,'Layer 3 | Process contrasts',
        'Are opposing-process differences\nassociated with SONFH?\n\n11 contrasts: z(A) - z(B)\nSex-adjusted HC3; BH across 11',
        fill='#F8F0E5',accent='#8A622E',fs=9.3,lty=3)
    label(4.5,2.30,'Three distinct questions. No sequential significance filter.',10.5,TRUE,center=TRUE)
    label(.25,1.94,'C',14,TRUE)
    label(.62,1.94,'Sensitivity analyses and exploratory candidate follow-up',11.7,TRUE)
    box(.25,.48,4.08,1.08,'Blood-based characterization',
        'Alternative scoring; gene membership; leave-one-out\nErythroid/globin expression proxy; component analysis',fs=9.2)
    box(4.58,.48,4.17,1.08,'Exploratory tissue assessment',
        'GSE345158: 5 SONFH + 4 fracture controls\nExact permutation; shared-gene sensitivity analysis',fs=9.2)
    label(4.5,.27,'Purpose: prioritize transcriptional candidates for longitudinal evaluation',10.3,TRUE,center=TRUE)
  }
  save_plot <- function(open_device) {open_device();on.exit(dev.off(),add=TRUE);draw()}
  save_plot(function() pdf(file.path(out,'Figure1_R.pdf'),width=9,height=8.1))
  save_plot(function() png(file.path(out,'Figure1_R.png'),width=9,height=8.1,units='in',res=300,bg='white'))
  tryCatch(save_plot(function() tiff(file.path(out,'Figure1_R.tiff'),width=9,height=8.1,units='in',
       res=600,compression='lzw',bg='white')),error=function(e) warning('TIFF export failed: ',conditionMessage(e)))
  if (capabilities('cairo')) save_plot(function() svg(file.path(out,'Figure1_R.svg'),width=9,height=8.1))
  legend <- paste('Figure 1. Study design and prespecified parallel analytical layers.',
    'Deposited RMA expression and archived GPL15207 annotation were combined for probe-to-gene aggregation.',
    'The blood discovery analysis used glucocorticoid-exposed participants without SONFH as controls.',
    'Individual associations, within-group differences in disease-associated effects, and opposing-process contrasts addressed distinct questions with separate Benjamini-Hochberg (BH) correction families.',
    'Constituent scores in Layers 2 and 3 were standardized across all 40 discovery samples.',
    'Follow-up analyses characterized candidate signals; the exploratory tissue comparison used fracture controls.',
    'Solid, dashed, and dotted borders distinguish the three analytical layers independently of color.',
    'RMA, robust multi-array average; GSVA, gene set variation analysis; HC3, heteroskedasticity-consistent covariance estimator; CR2, cluster-robust covariance estimator with small-sample correction.',
    'All gene-set scores represent relative transcriptional patterns, not direct measurements of physiological function, enzyme activity, or metabolic flux.')
  writeLines(legend,file.path(out,'Figure1_legend.txt'))
  passed <- TRUE
  message('PASS_FIGURE1_DEFINITIONS_AND_RENDER\nCounts and definitions checked; study-design figure generated.')
  message('No statistical models were refitted; this is a study-design figure.')
  message('Output folder: ',out)
  invisible(out)
}
run_figure1()

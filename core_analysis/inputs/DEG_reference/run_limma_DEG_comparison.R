# Exploratory conventional gene-level differential-expression comparison
# Run from the SONFH_Research project root in RStudio.
# Requires limma. This does NOT replace the frozen three-layer gene-set analysis.

suppressPackageStartupMessages({
  if (!requireNamespace('limma', quietly = TRUE)) {
    stop('Package limma is required. Install it with BiocManager::install("limma").')
  }
})

expr_file <- file.path('03_QC', 'SONFH_QC03_gene_expression_median_log2.csv.gz')
manifest_file <- file.path('03_QC', 'SONFH_sample_manifest_aligned.csv')
gmt_file <- file.path('01_Frozen_Master', 'SONFH_QC04_frozen_85_members_full.gmt')
out_dir <- file.path('09_DEG_Comparison')
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

expr_df <- read.csv(expr_file, check.names = FALSE, stringsAsFactors = FALSE)
stopifnot('Gene_Symbol' %in% names(expr_df))
rownames_expr <- make.unique(as.character(expr_df$Gene_Symbol))
expr <- as.matrix(expr_df[, setdiff(names(expr_df), 'Gene_Symbol'), drop = FALSE])
rownames(expr) <- rownames_expr
storage.mode(expr) <- 'numeric'

manifest <- read.csv(manifest_file, stringsAsFactors = FALSE)
stopifnot(all(colnames(expr) %in% manifest$GSM))
manifest <- manifest[match(colnames(expr), manifest$GSM), , drop = FALSE]
stopifnot(nrow(manifest) == ncol(expr), !anyNA(manifest$group), !anyNA(manifest$sex))

# non-SONFH is the reference; coefficient groupSONFH is SONFH minus non-SONFH.
manifest$group <- factor(ifelse(manifest$group == 'SONFH', 'SONFH', 'nonSONFH'),
                         levels = c('nonSONFH', 'SONFH'))
manifest$sex <- factor(manifest$sex)
design <- model.matrix(~ group + sex, data = manifest)
if (!('groupSONFH' %in% colnames(design))) stop('Could not identify groupSONFH coefficient.')

fit <- limma::lmFit(expr, design)
fit <- limma::eBayes(fit, trend = TRUE, robust = TRUE)
all_res <- limma::topTable(fit, coef = 'groupSONFH', number = Inf,
                           adjust.method = 'BH', sort.by = 'P')
all_res$Gene_Symbol <- rownames(all_res)
all_res <- all_res[, c('Gene_Symbol', setdiff(names(all_res), 'Gene_Symbol'))]
write.csv(all_res, file.path(out_dir, 'DEG_all_genes_limma_BH.csv'), row.names = FALSE)

# Read frozen gene sets, preserving their pre-defined membership.
gmt_lines <- readLines(gmt_file, warn = FALSE)
read_set <- function(set_name) {
  hit <- gmt_lines[startsWith(gmt_lines, paste0(set_name, '\t'))]
  if (length(hit) != 1L) stop(paste('Expected exactly one GMT row for', set_name))
  strsplit(hit, '\t', fixed = TRUE)[[1]][-(1:2)]
}
sets <- list(
  hydrogen_peroxide_catabolism = read_set('GOBP_HYDROGEN_PEROXIDE_CATABOLIC_PROCESS'),
  hydrogen_peroxide_biosynthesis = read_set('GOBP_HYDROGEN_PEROXIDE_BIOSYNTHETIC_PROCESS')
)
for (nm in names(sets)) {
  members <- data.frame(Gene_Symbol = sets[[nm]], stringsAsFactors = FALSE)
  members$Measurable_in_expression_matrix <- members$Gene_Symbol %in% rownames(expr)
  members <- merge(members, all_res, by = 'Gene_Symbol', all.x = TRUE, sort = FALSE)
  members <- members[match(sets[[nm]], members$Gene_Symbol), , drop = FALSE]
  write.csv(members, file.path(out_dir, paste0('DEG_', nm, '_members.csv')), row.names = FALSE)
}

summary <- data.frame(
  Dataset = 'GSE123568',
  N_nonSONFH = sum(manifest$group == 'nonSONFH'),
  N_SONFH = sum(manifest$group == 'SONFH'),
  N_genes_tested = nrow(expr),
  N_genes_nominal_P_lt_0_05 = sum(all_res$P.Value < 0.05, na.rm = TRUE),
  N_genes_BH_FDR_lt_0_05 = sum(all_res$adj.P.Val < 0.05, na.rm = TRUE),
  Model = 'limma linear model: expression ~ group + sex; eBayes(trend=TRUE, robust=TRUE)',
  FDR_scope = 'All genes tested in the gene-level analysis',
  stringsAsFactors = FALSE
)
write.csv(summary, file.path(out_dir, 'DEG_limma_summary.csv'), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, 'sessionInfo.txt'))
message('Finished. Review outputs in: ', out_dir)

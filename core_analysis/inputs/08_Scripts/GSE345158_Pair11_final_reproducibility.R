
library(GSVA)

# 1. Load archived analysis inputs
obj <- readRDS(
  "06_Validation/GSE345158_reproducibility_inputs.rds"
)

stopifnot(
  nrow(obj$gsva_matrix_sod2) == 18573,
  ncol(obj$gsva_matrix_sod2) == 9,
  length(obj$pair11_sets) == 2
)

# 2. Calculate GSVA scores
param <- gsvaParam(
  exprData = obj$gsva_matrix_sod2,
  geneSets = obj$pair11_sets,
  kcdf = "Gaussian",
  minSize = 1,
  maxSize = Inf
)

scores <- gsva(param, verbose = FALSE)

# 3. Standardize each program across 9 samples
z <- t(scale(t(scores)))

# 4. Calculate Pair 11
pair11 <- z[1, ] - z[2, ]

# 5. Align sample metadata
meta <- obj$sample_meta[
  colnames(z), ,
  drop = FALSE
]

stopifnot(
  identical(rownames(meta), colnames(z)),
  sum(meta$group == "SONFH") == 5,
  sum(meta$group == "Control") == 4
)

is_case <- meta$group == "SONFH"

# 6. Observed group difference
effect <- mean(pair11[is_case]) -
  mean(pair11[!is_case])

# 7. Exact permutation test
groups <- combn(seq_along(pair11), sum(is_case))

permuted <- apply(groups, 2, function(idx) {
  mean(pair11[idx]) - mean(pair11[-idx])
})

p_exact <- mean(
  abs(permuted) >= abs(effect) - 1e-10
)

# 8. Reproducibility checks
stopifnot(
  abs(effect - 1.29372) < 0.01,
  abs(p_exact - 0.03968254) < 1e-6
)

cat("Pair 11 effect:", effect, "\n")
cat("Exact permutation P:", p_exact, "\n")
cat("Reproducibility check: PASSED\n")


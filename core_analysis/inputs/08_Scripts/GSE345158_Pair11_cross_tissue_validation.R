
# =====================================================
# GSE345158 | SONFH cross-tissue validation
# Part 1: Load data and confirm sample identities
# =====================================================

library(readxl)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(DESeq2)
library(GSVA)

# 1. Read original count matrix
count_file <- paste0(
  "06_Validation/",
  "GSE345158_gene_count_exp.annotation.xlsx"
)

stopifnot(file.exists(count_file))

count_df <- read_excel(count_file)

# 2. Confirm expression matrix dimensions
stopifnot(
  nrow(count_df) == 20345,
  ncol(count_df) == 10,
  names(count_df)[1] == "gene_id"
)

# 3. Confirm exact sample names
expected_samples <- c(
  paste0("Case-", 1:5),
  "Ctrl-1", "Ctrl-2", "Ctrl-4", "Ctrl-5"
)


stopifnot(
  identical(names(count_df)[-1], expected_samples),
  !anyDuplicated(count_df$gene_id)
)

# 4. Construct sample metadata
sample_meta <- data.frame(
  sample = expected_samples,
  group = factor(
    c(rep("SONFH", 5), rep("Control", 4)),
    levels = c("Control", "SONFH")
  ),
  row.names = expected_samples
)

# 5. Confirm group sizes
print(table(sample_meta$group))

cat("Part 1 passed: data and sample checks\n")
# =====================================================
# Part 6: Save SOD2 annotation sensitivity results
# =====================================================

write.csv(
  comparison,
  "06_Validation/GSE345158_SOD2_sensitivity.csv",
  row.names = FALSE
)

write.csv(
  data.frame(
    Sample = colnames(scores_sod2),
    Group = sample_meta[colnames(scores_sod2), "group"],
    Pair11_corrected = diff_sod2
  ),
  "06_Validation/GSE345158_SOD2_corrected_scores.csv",
  row.names = FALSE
)

cat("Part 6 passed: SOD2 correction saved\n")# 1. 获取病例组和对照组的样本名称
case_samples <- rownames(sample_meta)[
  sample_meta$group == "SONFH"
]

control_samples <- rownames(sample_meta)[
  sample_meta$group == "Control"
]

# 2. 分别计算两个程序的组间均值差
component_effects <- data.frame(
  Program = c("H2O2_biosynthesis", "H2O2_catabolism"),
  SONFH_mean = rowMeans(
    z_sod2[, case_samples, drop = FALSE]
  ),
  Control_mean = rowMeans(
    z_sod2[, control_samples, drop = FALSE]
  )
)

component_effects$Difference <-
  component_effects$SONFH_mean -
  component_effects$Control_mean

# 3. 查看结果
print(round(component_effects[, -1], 3))

# 4. 检查两个程序的差值是否等于Pair 11效应
cat(
  "Pair 11组间差异：",
  round(
    component_effects$Difference[1] -
      component_effects$Difference[2],
    3
  ),
  "\n"

)

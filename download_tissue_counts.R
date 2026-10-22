# Optional download of the exact public count input used in the study.
# Open in RStudio and click Source; uses only base R.
local({
  url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE345nnn/GSE345158/suppl/GSE345158_gene_count_exp.annotation.xlsx"
  expected <- "5ae538bfa2cc1c5aee87f539e8be2768"
  dest <- file.path(getwd(), "GSE345158_gene_count_exp.annotation.xlsx")
  tmp <- tempfile(fileext = ".xlsx")
  on.exit(unlink(tmp), add = TRUE)
  download.file(url, tmp, mode = "wb", method = "libcurl", quiet = FALSE)
  actual <- unname(tools::md5sum(tmp))
  if (!identical(actual, expected)) stop("Downloaded file differs from the archived study input: ", actual)
  if (file.exists(dest)) {
    if (!identical(unname(tools::md5sum(dest)), expected)) stop("Existing local file differs; it was not overwritten.")
  } else if (!file.copy(tmp, dest)) stop("Could not save downloaded count file.")
  cat("PASS: exact public count input verified.\n", normalizePath(dest, winslash = "/"), "\n")
})

# Integrated convenience launcher; syntax-checked but not executed as one job.
# Its component computations have separate execution evidence in results/.
# From project root: source("code/run_all_R.R", encoding="UTF-8")
local({
  stopifnot(file.exists('inputs/frozen_85.gmt'))
  cat('1/3: blood deposited probes to gene expression\n')
  source('code/reproduce_blood_upstream.R',encoding='UTF-8',local=new.env(parent=globalenv()))
  cat('2/3: native gene-set scoring, CR2, removal and tissue sensitivity\n')
  source('code/verify_scoring_and_CR2.R',encoding='UTF-8',local=new.env(parent=globalenv()))
  cat('3/3: tissue integer counts to corrected Pair11, fixed matched recipe\n')
  source('code/reproduce_tissue_upstream.R',encoding='UTF-8',local=new.env(parent=globalenv()))
})

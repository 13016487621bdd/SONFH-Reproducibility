
analyze_family_v2 <- function(family_name) {
  
  sets <- family_members[[family_name]]
  x <- t(gsva_z[sets, , drop = FALSE])
  
  contrast <- sweep(
    x[, -1, drop = FALSE],
    1,
    x[, 1],
    "-"
  )
  
  # 两个基因集：单系数检验
  if (length(sets) == 2) {
    
    dat <- sample_aligned
    dat$score_diff <- as.numeric(contrast[, 1])
    
    fit <- lm(score_diff ~ group + sex, data = dat)
    
    V <- clubSandwich::vcovCR(
      fit,
      cluster = dat$GSM,
      type = "CR2"
    )
    
    test <- clubSandwich::coef_test(
      fit,
      vcov = V,
      test = "Satterthwaite"
    )
    
    d <- test[test$Coef == "groupSONFH", ]
    
    return(data.frame(
      family = family_name,
      n_sets = length(sets),
      method = "CR2_Satterthwaite",
      statistic = as.numeric(d[["tstat"]])^2,
      df_num = 1,
      df_denom = as.numeric(d[["df_Satt"]]),
      p = as.numeric(d[["p_Satt"]])
    ))
  }
  
  # 三个及以上基因集：联合检验
  fit <- lm(
    contrast ~ group + sex,
    data = sample_aligned
  )
  
  V <- clubSandwich::vcovCR(
    fit,
    cluster = sample_aligned$GSM,
    type = "CR2"
  )
  
  idx <- grep(":groupSONFH$", rownames(V))
  
  stopifnot(
    length(idx) == ncol(contrast),
    all(is.finite(V))
  )
  
  C <- diag(nrow(V))[idx, , drop = FALSE]
  
  w <- clubSandwich::Wald_test(
    fit,
    constraints = C,
    vcov = V,
    test = "HTZ"
  )
  
  data.frame(
    family = family_name,
    n_sets = length(sets),
    method = "CR2_HTZ",
    statistic = as.numeric(w$Fstat),
    df_num = as.numeric(w$df_num),
    df_denom = as.numeric(w$df_denom),
    p = as.numeric(w$p_val)
  )
}

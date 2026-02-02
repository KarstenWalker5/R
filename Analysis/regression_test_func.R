regression_check_table <- function(
    model,
    alpha = 0.05,
    corr_breaks = c(none = 0.70, low = 0.85, medium = 0.95), # >= medium => high
    vif_breaks  = c(none = 5, low = 10, medium = 20),        # >= medium => high
    influence_rules = list(
      cooks = function(n) 4 / n,
      leverage = function(n, p) 2 * p / n,
      std_resid = 3
    )
) {
  # ---- dependencies ----
  req <- c("broom", "tibble", "dplyr", "stats", "lmtest", "sandwich", "car")
  missing <- req[!vapply(req, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop(
      "Missing packages: ", paste(missing, collapse = ", "),
      ". Install them first, e.g. install.packages(c(",
      paste0('"', missing, '"', collapse = ", "),
      "))"
    )
  }
  
  # ---- helper functions ----
  risk_bucket <- function(x, breaks_named) {
    b <- breaks_named
    if (is.na(x)) return(NA_character_)
    if (x < b[["none"]]) "none"
    else if (x < b[["low"]]) "low"
    else if (x < b[["medium"]]) "medium"
    else "high"
  }
  
  mk_row <- function(check, statistic, value = NA_real_, p_value = NA_real_,
                     threshold = NA_character_, pass = NA, note = NA_character_) {
    tibble::tibble(
      check = check,
      statistic = statistic,
      value = value,
      p_value = p_value,
      threshold = threshold,
      pass = pass,
      note = note
    )
  }
  
  # ---- core objects ----
  n <- stats::nobs(model)
  p <- length(stats::coef(model))
  aug <- broom::augment(model)
  res <- stats::resid(model)
  
  out <- tibble::tibble()
  
  # ---- 1) Normality (Shapiro) ----
  sh <- stats::shapiro.test(res)
  out <- dplyr::bind_rows(
    out,
    mk_row(
      check = "Normality",
      statistic = "Shapiro-Wilk",
      value = unname(sh$statistic),
      p_value = sh$p.value,
      threshold = paste0("p > ", alpha),
      pass = sh$p.value > alpha,
      note = if (n > 5000) "Large n: Shapiro often rejects small deviations." else NA_character_
    )
  )
  
  # ---- 2) Heteroskedasticity (Breusch–Pagan) ----
  bp <- lmtest::bptest(model)
  out <- dplyr::bind_rows(
    out,
    mk_row(
      check = "Homoscedasticity",
      statistic = "Breusch-Pagan",
      value = unname(bp$statistic),
      p_value = bp$p.value,
      threshold = paste0("p > ", alpha),
      pass = bp$p.value > alpha,
      note = "Fail suggests heteroskedasticity; consider robust (HC3) SEs."
    )
  )
  
  # ---- 3) Independence (Durbin–Watson) ----
  dw <- lmtest::dwtest(model)
  out <- dplyr::bind_rows(
    out,
    mk_row(
      check = "Independence",
      statistic = "Durbin-Watson",
      value = unname(dw$statistic),
      p_value = dw$p.value,
      threshold = paste0("p > ", alpha),
      pass = dw$p.value > alpha,
      note = "Only interpretable for ordered/time-series residuals."
    )
  )
  
  # ---- 4) Outliers (car::outlierTest) ----
  ot <- tryCatch(car::outlierTest(model), error = function(e) NULL)
  
  if (is.null(ot)) {
    out <- dplyr::bind_rows(
      out,
      mk_row(
        check = "Outliers",
        statistic = "Bonferroni outlier test",
        threshold = paste0("p > ", alpha),
        pass = TRUE,
        note = "No outliers flagged by car::outlierTest()."
      )
    )
  } else {
    rstud   <- if (!is.null(ot$rstudent)) as.numeric(ot$rstudent) else NA_real_
    p_bonf  <- if (!is.null(ot$bonf.p)) as.numeric(ot$bonf.p) else NA_real_
    p_unadj <- if (!is.null(ot$unadj.p)) as.numeric(ot$unadj.p) else NA_real_
    
    pval <- if (!is.na(p_bonf)) p_bonf else p_unadj
    
    obs_label <- NA_character_
    if (!is.null(names(ot$rstudent)) && length(names(ot$rstudent)) == 1) {
      obs_label <- names(ot$rstudent)
    } else if (!is.null(ot$obs)) {
      obs_label <- as.character(ot$obs)
    }
    
    out <- dplyr::bind_rows(
      out,
      mk_row(
        check = "Outliers",
        statistic = "Bonferroni outlier test",
        value = rstud,
        p_value = pval,
        threshold = paste0("p > ", alpha),
        pass = ifelse(is.na(pval), NA, pval > alpha),
        note = ifelse(
          is.na(obs_label),
          "Most extreme obs flagged.",
          paste0("Most extreme obs: ", obs_label)
        )
      )
    )
  }
  
  # ---- 5) Influence heuristics (counts) ----
  cooks_thr <- influence_rules$cooks(n)
  lev_thr   <- influence_rules$leverage(n, p)
  sr_thr    <- influence_rules$std_resid
  
  n_cooks <- sum(aug$.cooksd > cooks_thr, na.rm = TRUE)
  n_lev   <- sum(aug$.hat > lev_thr, na.rm = TRUE)
  n_sr    <- sum(abs(aug$.std.resid) > sr_thr, na.rm = TRUE)
  
  out <- dplyr::bind_rows(
    out,
    mk_row(
      check = "Influence",
      statistic = "Count Cook's D > 4/n",
      value = n_cooks,
      threshold = paste0("0 (rule: 4/n = ", signif(cooks_thr, 3), ")"),
      pass = (n_cooks == 0),
      note = "Heuristic; investigate flagged points."
    ),
    mk_row(
      check = "Influence",
      statistic = "Count leverage > 2p/n",
      value = n_lev,
      threshold = paste0("0 (rule: 2p/n = ", signif(lev_thr, 3), ")"),
      pass = (n_lev == 0),
      note = "Heuristic; high leverage points can dominate fit."
    ),
    mk_row(
      check = "Influence",
      statistic = "Count |std resid| > 3",
      value = n_sr,
      threshold = "0 (rule: |std resid| > 3)",
      pass = (n_sr == 0),
      note = "Heuristic; large residuals may indicate outliers/misspecification."
    )
  )
  
  # ---- 6) Collinearity risk (max |corr| + max VIF/GVIF) ----
  mf <- stats::model.frame(model)
  
  x_df <- mf %>%
    dplyr::select(-1) %>%
    dplyr::select(where(is.numeric))
  
  max_abs_corr <- NA_real_
  corr_risk <- NA_character_
  
  if (ncol(x_df) >= 2) {
    cm <- stats::cor(x_df, use = "pairwise.complete.obs")
    max_abs_corr <- max(abs(cm[upper.tri(cm)]), na.rm = TRUE)
    corr_risk <- risk_bucket(max_abs_corr, corr_breaks)
  }
  
  max_vif <- NA_real_
  vif_risk <- NA_character_
  
  vif_tbl <- tryCatch(car::vif(model), error = function(e) NULL)
  if (!is.null(vif_tbl)) {
    if (is.matrix(vif_tbl) && all(c("GVIF", "Df") %in% colnames(vif_tbl))) {
      adj <- vif_tbl[, "GVIF"]^(1 / (2 * vif_tbl[, "Df"]))
      max_vif <- max(adj, na.rm = TRUE)
    } else if (is.matrix(vif_tbl)) {
      max_vif <- max(vif_tbl[, 1], na.rm = TRUE)
    } else {
      max_vif <- max(as.numeric(vif_tbl), na.rm = TRUE)
    }
    vif_risk <- risk_bucket(max_vif, vif_breaks)
  }
  
  out <- dplyr::bind_rows(
    out,
    mk_row(
      check = "Collinearity",
      statistic = "Max |correlation| among numeric predictors",
      value = max_abs_corr,
      threshold = paste0(
        "none<", corr_breaks[["none"]],
        " low<", corr_breaks[["low"]],
        " medium<", corr_breaks[["medium"]],
        " high≥", corr_breaks[["medium"]]
      ),
      pass = NA,
      note = paste0("Risk: ", corr_risk)
    ),
    mk_row(
      check = "Collinearity",
      statistic = "Max VIF (or adjusted GVIF)",
      value = max_vif,
      threshold = paste0(
        "none<", vif_breaks[["none"]],
        " low<", vif_breaks[["low"]],
        " medium<", vif_breaks[["medium"]],
        " high≥", vif_breaks[["medium"]]
      ),
      pass = NA,
      note = paste0("Risk: ", vif_risk)
    )
  )
  
  out
}

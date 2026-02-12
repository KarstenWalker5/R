pre_treatment_diagnostics <- function(data,
                                      outcome = "sales",
                                      treat = "treat",
                                      time = "date",
                                      post = "post",
                                      unit = "city",
                                      intervention_date = NULL,
                                      do_event_study = FALSE) {
  # --- tidy-eval setup ---
  outcome_sym <- rlang::sym(outcome)
  treat_sym   <- rlang::sym(treat)
  time_sym    <- rlang::sym(time)
  post_sym    <- rlang::sym(post)
  unit_sym    <- rlang::sym(unit)
  
  df <- data
  
  # If intervention_date not supplied, infer as first post==1 date
  if (is.null(intervention_date)) {
    intervention_date <- df %>%
      dplyr::filter(!!post_sym == 1) %>%
      dplyr::summarize(first_post = min(!!time_sym)) %>%
      dplyr::pull(first_post)
  }
  
  # -------------------------
  # 1) Pre-period subset
  # -------------------------
  pre_data <- df %>% dplyr::filter(!!post_sym == 0)
  
  # -------------------------
  # 2) Equal means test (t-test)
  # -------------------------
  mean_formula <- stats::reformulate(termlabels = treat, response = outcome)
  mean_test <- stats::t.test(mean_formula, data = pre_data)
  
  # -------------------------
  # 3) Variance tests (Levene + Fligner)
  #    use reformulate() to build formula safely
  # -------------------------
  group_term  <- paste0("factor(", treat, ")")
  lev_formula <- stats::reformulate(termlabels = group_term, response = outcome)
  
  # Levene (median-centered)
  lev_out <- car::leveneTest(
    lev_formula,
    data   = pre_data,
    center = median
  )
  
  # Fligner–Killeen (robust, nonparametric)
  fligner_out <- stats::fligner.test(
    lev_formula,
    data = pre_data
  )
  
  levene_p  <- lev_out[1, "Pr(>F)"]
  fligner_p <- fligner_out$p.value
  
  # -------------------------
  # 3b) KS test: distributional equality (treated vs control)
  # -------------------------
  pre_treated <- pre_data %>%
    dplyr::filter(!!treat_sym == 1) %>%
    dplyr::pull(!!outcome_sym)
  
  pre_control <- pre_data %>%
    dplyr::filter(!!treat_sym == 0) %>%
    dplyr::pull(!!outcome_sym)
  
  # Only run KS if both groups have data
  if (length(pre_treated) > 0 && length(pre_control) > 0) {
    ks_out <- stats::ks.test(pre_treated, pre_control)
    ks_p   <- ks_out$p.value
  } else {
    ks_out <- NULL
    ks_p   <- NA_real_
    warning("KS test not run: one of the groups has zero pre-period observations.")
  }
  
  # -------------------------
  # 4) Parallel trends slope test
  #    outcome ~ as.numeric(time) * factor(treat)  (pre only)
  # -------------------------
  slope_formula <- as.formula(
    paste0(
      outcome,
      " ~ as.numeric(", time, ") * factor(", treat, ")"
    )
  )
  
  slope_model   <- stats::lm(slope_formula, data = pre_data)
  slope_summary <- summary(slope_model)
  
  # interaction term is the "differential slope"
  interaction_term <- grep("as.numeric\\(", 
                           rownames(slope_summary$coefficients), 
                           value = TRUE)[1]
  
  slope_p <- slope_summary$coefficients[interaction_term, "Pr(>|t|)"]
  
  # -------------------------
  # 5) Pre-period parallel trends plot
  # -------------------------
  pre_trends_plot <- pre_data %>%
    dplyr::mutate(group = dplyr::if_else(!!treat_sym == 1, "Treated", "Control")) %>%
    dplyr::group_by(!!time_sym, group) %>%
    dplyr::summarize(avg_outcome = mean(!!outcome_sym), .groups = "drop") %>%
    ggplot2::ggplot(ggplot2::aes(x = !!time_sym, y = avg_outcome, color = group)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::labs(
      title = paste("Pre-Period Parallel Trends:", outcome),
      x     = rlang::as_label(time_sym),
      y     = paste("Average", outcome),
      color = "Group"
    ) +
    ggplot2::theme_minimal(base_size = 14)
  
  # -------------------------
  # 6) Optional: Event-study model (full sample)
  # -------------------------
  event_study_model <- NULL
  if (do_event_study) {
    df_es <- df %>%
      dplyr::mutate(rel_time = as.integer(!!time_sym - intervention_date))
    
    # Build fixest formula: outcome ~ i(rel_time, treat, ref = -1) | unit + time
    fml_str <- paste0(
      outcome,
      " ~ i(rel_time, ", treat, ", ref = -1) | ",
      unit, " + ", time
    )
    
    event_study_model <- fixest::feols(
      fml     = as.formula(fml_str),
      data    = df_es,
      cluster = as.formula(paste0("~", unit))
    )
  }
  
  # -------------------------
  # 7) Bundle everything into a list
  # -------------------------
  results <- list(
    pre_data          = pre_data,
    mean_test         = mean_test,
    levene_test       = lev_out,
    levene_p_value    = levene_p,
    fligner_test      = fligner_out,
    fligner_p_value   = fligner_p,
    ks_test           = ks_out,
    ks_p_value        = ks_p,
    slope_model       = slope_model,
    slope_p_value     = slope_p,
    pre_trends_plot   = pre_trends_plot,
    event_study_model = event_study_model
  )
  
  # Quick console summary
  cat("\n=== Pre-Treatment Diagnostics ===\n")
  cat("Equal-means t-test (", outcome, "): p-value = ",
      signif(mean_test$p.value, 3), "\n", sep = "")
  cat("Levene's test (variance equal): p-value = ",
      signif(levene_p, 3), "\n", sep = "")
  cat("Fligner–Killeen test (variance equal): p-value = ",
      signif(fligner_p, 3), "\n", sep = "")
  cat("KS test (distributional equality, treated vs control): p-value = ",
      signif(ks_p, 3), "\n", sep = "")
  cat("Parallel-trend slope test (time*treat interaction): p-value = ",
      signif(slope_p, 3), "\n", sep = "")
  if (do_event_study) {
    cat("Event study model estimated: yes (use iplot(results$event_study_model))\n")
  } else {
    cat("Event study model estimated: no (set do_event_study = TRUE to add)\n")
  }
  cat("=================================\n")
  
  invisible(results)
}

pre_treatment_diagnostics_panel <- function(data,
                                            outcome,
                                            treat,
                                            time,
                                            unit,
                                            intervention_date,
                                            transform = c("none", "log1p"),
                                            agg = mean,
                                            do_event_study = FALSE,
                                            ...) {
  transform <- match.arg(transform)
  
  stopifnot(all(c(outcome, treat, time, unit) %in% names(data)))
  stopifnot(!is.null(intervention_date))
  
  df <- data %>%
    dplyr::mutate(
      .time  = .data[[time]],
      .treat = .data[[treat]],
      .unit  = .data[[unit]],
      .y_raw = .data[[outcome]],
      .y = dplyr::if_else(transform == "log1p", log1p(.y_raw), as.numeric(.y_raw)),
      post = dplyr::if_else(.time >= intervention_date, 1L, 0L)
    )
  
  # Collapse to treated/control mean series by time
  collapsed <- df %>%
    dplyr::group_by(.time, .treat) %>%
    dplyr::summarise(
      y_mean  = agg(.y, na.rm = TRUE),          # <-- renamed from "outcome" to "y_mean"
      n_units = dplyr::n_distinct(.unit),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      post = dplyr::if_else(.time >= intervention_date, 1L, 0L),
      grp_unit = dplyr::if_else(.treat == 1, "TreatedMean", "ControlMean")
    )
  
  res <- pre_treatment_diagnostics(
    data = collapsed %>%
      dplyr::rename(!!time := .time, !!treat := .treat, !!unit := grp_unit),
    outcome = "y_mean",                          # <-- pass the new safe name
    treat   = treat,
    time    = time,
    post    = "post",
    unit    = unit,
    intervention_date = intervention_date,
    do_event_study = do_event_study,
    ...
  )
  
  list(
    intervention_date = intervention_date,
    collapsed_means = collapsed,
    results = res
  )
}

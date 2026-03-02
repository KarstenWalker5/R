pre_treatment_diagnostics <- function(data,
                                      outcome = "sales",
                                      treat = "treat",
                                      time = "date",
                                      post = "post",
                                      unit = "city",
                                      intervention_date = NULL,
                                      do_event_study = FALSE,
                                      rolling_window = 4,
                                      spline_df = 4,
                                      quantiles = c(0.1, 0.5, 0.9)) {
  
  theme_fancy <- function() {
    ggplot2::theme_minimal(base_size = 13) +
      ggplot2::theme(
        panel.grid.minor = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(face = "bold"),
        plot.subtitle = ggplot2::element_text(color = "grey35")
      )
  }
  
  stopifnot(all(c(outcome, treat, time, post, unit) %in% names(data)))
  stopifnot(!is.null(intervention_date))
  
  df <- data %>%
    dplyr::mutate(
      .time  = .data[[time]],
      .treat = .data[[treat]],
      .post  = .data[[post]],
      .unit  = .data[[unit]],
      .y     = as.numeric(.data[[outcome]])
    ) %>%
    dplyr::arrange(.unit, .time)
  
  pre_df <- df %>% dplyr::filter(.time < intervention_date)
  
  collapsed <- df %>%
    dplyr::group_by(.time, .treat) %>%
    dplyr::summarise(
      y_mean = mean(.y, na.rm = TRUE),
      n_units = dplyr::n_distinct(.unit),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      post = dplyr::if_else(.time >= intervention_date, 1L, 0L)
    )
  
  slope_df <- collapsed %>%
    dplyr::filter(post == 0) %>%
    dplyr::arrange(.treat, .time) %>%
    dplyr::group_by(.treat) %>%
    dplyr::mutate(
      dy = y_mean - dplyr::lag(y_mean),
      dt = as.numeric(.time - dplyr::lag(.time)),
      slope = dy / dt
    ) %>%
    dplyr::ungroup()
  
  p_trends <- ggplot2::ggplot(collapsed,
                              ggplot2::aes(x = .time, y = y_mean, color = factor(.treat))) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_vline(xintercept = intervention_date, linetype = "dashed") +
    ggplot2::labs(
      title = "Mean outcome over time by group",
      x = time,
      y = outcome,
      color = "treat"
    ) +
    theme_fancy()
  
  p_slopes <- ggplot2::ggplot(slope_df,
                              ggplot2::aes(x = .time, y = slope, color = factor(.treat))) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dotted") +
    ggplot2::labs(
      title = "Pre-period slopes",
      x = time,
      y = "slope",
      color = "treat"
    ) +
    theme_fancy()
  
  event_study <- NULL
  if (isTRUE(do_event_study)) {
    df_es <- df %>%
      dplyr::mutate(rel_time = as.integer(.time - intervention_date))
    
    event_study <- stats::lm(.y ~ factor(rel_time) * .treat, data = df_es)
  }
  
  list(
    collapsed = collapsed,
    slope_df = slope_df,
    plots = list(trends = p_trends, slopes = p_slopes),
    event_study = event_study
  )
}

pre_treatment_diagnostics_panel <- function(
    data,
    outcome,
    treat,
    time,
    unit,
    intervention_date,
    transform = c("log1p", "none"),
    gap_window = 14,
    do_event_study = FALSE,
    cluster_unit = NULL
) {
  # ---- packages ----
  requireNamespace("dplyr", quietly = TRUE)
  requireNamespace("tidyr", quietly = TRUE)
  requireNamespace("ggplot2", quietly = TRUE)
  requireNamespace("fixest", quietly = TRUE)
  
  has_slider <- requireNamespace("slider", quietly = TRUE)
  if (!has_slider) {
    message("Package 'slider' not found — using zoo::rollapply as fallback.")
    requireNamespace("zoo", quietly = TRUE)
  }
  
  transform <- match.arg(transform)
  if (is.null(cluster_unit)) cluster_unit <- unit
  
  # ---- checks ----
  stopifnot(is.data.frame(data))
  stopifnot(all(c(outcome, treat, time, unit) %in% names(data)))
  stopifnot(!is.null(intervention_date))
  
  # ---- prepare ----
  df <- data %>%
    dplyr::mutate(
      .time  = as.Date(.data[[time]]),
      .treat = as.integer(.data[[treat]]),
      .unit  = .data[[unit]],
      .y_raw = .data[[outcome]],
      rel_time = as.integer(as.Date(.data[[time]]) - as.Date(intervention_date)),
      post     = dplyr::if_else(as.Date(.data[[time]]) >= as.Date(intervention_date), 1L, 0L)
    )
  
  df <- df %>%
    dplyr::mutate(
      .y = dplyr::case_when(
        transform == "log1p" ~ log1p(as.numeric(.y_raw)),
        TRUE ~ as.numeric(.y_raw)
      ),
      dow = factor(strftime(.time, "%u")) # 1..7 weekday
    )
  
  # ---- collapse means (used for gap plots) ----
  collapsed_means <- df %>%
    dplyr::group_by(.time, .treat) %>%
    dplyr::summarise(
      y_mean  = mean(.y, na.rm = TRUE),
      n_units = dplyr::n_distinct(.unit),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      post = dplyr::if_else(.time >= as.Date(intervention_date), 1L, 0L),
      grp_unit = dplyr::if_else(.treat == 1, "TreatedMean", "ControlMean")
    ) %>%
    dplyr::arrange(.time)
  
  # ---- gap series (treated - control) ----
  gap_series <- collapsed_means %>%
    tidyr::pivot_wider(names_from = .treat, values_from = y_mean, names_prefix = "treat_") %>%
    dplyr::mutate(
      treat_0 = if ("treat_0" %in% names(.)) treat_0 else NA_real_,
      treat_1 = if ("treat_1" %in% names(.)) treat_1 else NA_real_
    ) %>%
    dplyr::rename(control = treat_0, treated = treat_1) %>%
    dplyr::mutate(gap = treated - control) %>%
    dplyr::select(.time, control, treated, gap) %>%
    dplyr::arrange(.time)
  
  # ---- rolling gap stats (PRE ONLY) ----
  gap_rolling_df <- gap_series %>%
    dplyr::filter(.time < as.Date(intervention_date)) %>%
    dplyr::arrange(.time)
  
  if (nrow(gap_rolling_df) > 0) {
    if (has_slider) {
      gap_rolling_df <- gap_rolling_df %>%
        dplyr::mutate(
          roll_mean = slider::slide_dbl(gap, mean, .before = gap_window - 1, .complete = TRUE),
          roll_sd   = slider::slide_dbl(gap, sd,   .before = gap_window - 1, .complete = TRUE)
        )
    } else {
      gap_rolling_df <- gap_rolling_df %>%
        dplyr::mutate(
          roll_mean = zoo::rollapply(gap, width = gap_window, FUN = mean, align = "right", fill = NA),
          roll_sd   = zoo::rollapply(gap, width = gap_window, FUN = sd,   align = "right", fill = NA)
        )
    }
  } else {
    gap_rolling_df <- dplyr::mutate(gap_rolling_df, roll_mean = NA_real_, roll_sd = NA_real_)
  }
  
  gap_rolling_plot <- ggplot2::ggplot(gap_rolling_df, ggplot2::aes(x = .time)) +
    ggplot2::geom_line(ggplot2::aes(y = roll_mean), linewidth = 0.9, na.rm = TRUE) +
    ggplot2::geom_line(ggplot2::aes(y = roll_sd), linetype = "dashed", linewidth = 0.9, na.rm = TRUE) +
    ggplot2::geom_vline(xintercept = as.Date(intervention_date), linetype = "dotted") +
    ggplot2::labs(
      title = sprintf("Gap stability (pre): %d-day rolling mean (solid) & SD (dashed)", gap_window),
      subtitle = paste0("Outcome: ", transform, "(", outcome, ")"),
      x = "Date",
      y = "Rolling statistic"
    ) +
    ggplot2::theme_minimal()
  
  # ---- helper: extract p-value from fixest coeftable ----
  get_coef_p <- function(mod, pattern) {
    if (is.null(mod)) return(NA_real_)
    ct <- tryCatch(fixest::coeftable(mod), error = function(e) NULL)
    if (is.null(ct)) return(NA_real_)
    rn <- rownames(ct)
    hit <- grep(pattern, rn)
    if (length(hit) == 0) return(NA_real_)
    # if multiple hits, return the smallest p-value (conservative)
    pv <- ct[hit, "Pr(>|t|)"]
    suppressWarnings(min(as.numeric(pv), na.rm = TRUE))
  }
  
  # ---- helper: joint wald p-value for interaction terms ----
  get_joint_p <- function(mod, keep_regex) {
    if (is.null(mod)) return(NA_real_)
    w <- tryCatch(fixest::wald(mod, keep = keep_regex), error = function(e) NULL)
    if (is.null(w)) return(NA_real_)
    # fixest wald returns an object where pvalue is often in $p.value
    if (!is.null(w$p.value)) return(as.numeric(w$p.value))
    if (!is.null(w[["p.value"]])) return(as.numeric(w[["p.value"]]))
    NA_real_
  }
  
  cluster_formula <- stats::as.formula(paste0("~", cluster_unit))
  
  # ---- PANEL-APPROPRIATE MODELS (PRE ONLY) ----
  
  # 1) Linear differential pretrend
  linear_pretrend_model <- tryCatch(
    fixest::feols(
      .y ~ rel_time * .treat | .unit + .time,
      data = df %>% dplyr::filter(rel_time < 0),
      cluster = cluster_formula
    ),
    error = function(e) NULL
  )
  
  # 2) Curvature / nonlinearity (quadratic interacted)
  curvature_model <- tryCatch(
    fixest::feols(
      .y ~ poly(rel_time, 2, raw = TRUE) * .treat | .unit + .time,
      data = df %>% dplyr::filter(rel_time < 0),
      cluster = cluster_formula
    ),
    error = function(e) NULL
  )
  
  # 3) Seasonality-adjusted differential pretrend (weekday dummies)
  seasonality_model <- tryCatch(
    fixest::feols(
      .y ~ rel_time * .treat + i(dow) | .unit + .time,
      data = df %>% dplyr::filter(rel_time < 0),
      cluster = cluster_formula
    ),
    error = function(e) NULL
  )
  
  # 4) Optional: Pre-only event study (best visual diagnostic)
  event_study_model <- NULL
  event_study_joint_test <- NULL
  if (isTRUE(do_event_study)) {
    event_study_model <- tryCatch(
      fixest::feols(
        .y ~ fixest::i(rel_time, .treat, ref = -1) | .unit + .time,
        data = df %>% dplyr::filter(rel_time < 0),
        cluster = cluster_formula
      ),
      error = function(e) NULL
    )
    event_study_joint_test <- get_joint_p(event_study_model, "rel_time::")
  }
  
  # ---- summary table (like the full function, but only appropriate checks) ----
  # p-values
  p_linear  <- get_coef_p(linear_pretrend_model, "rel_time.*:.*\\.treat|rel_time:.*treat|rel_time:.*\\.treat|rel_time.*\\.treat")
  p_curv    <- get_joint_p(curvature_model, ":\\.treat")   # joint test of quadratic interaction terms
  p_season  <- get_coef_p(seasonality_model, "rel_time.*:.*\\.treat|rel_time:.*treat|rel_time:.*\\.treat|rel_time.*\\.treat")
  p_es      <- if (isTRUE(do_event_study)) event_study_joint_test else NA_real_
  
  # rolling gap summary (descriptive)
  roll_sd_mean  <- suppressWarnings(mean(gap_rolling_df$roll_sd, na.rm = TRUE))
  roll_sd_max   <- suppressWarnings(max(gap_rolling_df$roll_sd,  na.rm = TRUE))
  roll_mean_abs_max <- suppressWarnings(max(abs(gap_rolling_df$roll_mean), na.rm = TRUE))
  
  summary_table <- data.frame(
    check = c(
      "Linear differential pretrend (rel_time × treat) p-value",
      "Curvature / nonlinearity (quadratic × treat) joint p-value",
      "Seasonality-adjusted differential pretrend p-value",
      "Event study (pre leads) joint p-value",
      sprintf("Gap rolling SD (%d-day) mean", gap_window),
      sprintf("Gap rolling SD (%d-day) max", gap_window),
      sprintf("Gap rolling mean |max| (%d-day)", gap_window)
    ),
    value = c(
      p_linear,
      p_curv,
      p_season,
      p_es,
      roll_sd_mean,
      roll_sd_max,
      roll_mean_abs_max
    ),
    stringsAsFactors = FALSE
  )
  
  # ---- print summary to console (similar to full function) ----
  cat("\n=== Pre-Treatment Diagnostics (1 treated / 1 control mode) ===\n")
  cat("Outcome:", transform, "(", outcome, ")\n", sep = "")
  cat("Linear differential pretrend p-value: ", signif(p_linear, 3), "\n", sep = "")
  cat("Curvature (joint) p-value: ", signif(p_curv, 3), "\n", sep = "")
  cat("Seasonality-adjusted pretrend p-value: ", signif(p_season, 3), "\n", sep = "")
  if (isTRUE(do_event_study)) {
    cat("Event-study pre-leads joint p-value: ", signif(p_es, 3), "\n", sep = "")
  } else {
    cat("Event-study pre-leads joint p-value: (not computed; set do_event_study = TRUE)\n")
  }
  cat(sprintf("Gap rolling SD (%d-day): mean = %s, max = %s\n",
              gap_window,
              signif(roll_sd_mean, 3),
              signif(roll_sd_max, 3)))
  cat(sprintf("Gap rolling mean (%d-day): max |mean| = %s\n",
              gap_window,
              signif(roll_mean_abs_max, 3)))
  cat("=============================================================\n")
  
  # ---- return object ----
  out <- list(
    intervention_date = as.Date(intervention_date),
    n_treated_units = length(unique(df$.unit[df$.treat == 1])),
    n_control_units = length(unique(df$.unit[df$.treat == 0])),
    collapsed_means = collapsed_means,
    gap_series = gap_series,
    gap_rolling_df = gap_rolling_df,
    gap_rolling_plot = gap_rolling_plot,
    linear_pretrend_model = linear_pretrend_model,
    curvature_model = curvature_model,
    seasonality_model = seasonality_model,
    event_study_model = event_study_model,
    event_study_joint_p_value = event_study_joint_test,
    summary_table = summary_table
  )
  
  class(out) <- c("precheck_panel", class(out))
  return(out)
}

pre_treatment_diagnostics <- function(data,
                                      outcome = "sales",
                                      treat = "treat",
                                      time = "date",
                                      post = "post",
                                      unit = "city",
                                      intervention_date = NULL,
                                      do_event_study = FALSE,
                                      # ---- new knobs for added diagnostics ----
                                      rolling_window = 4,          # window length in time points (pre only)
                                      spline_df = 4,               # degrees of freedom for ns() in curvature test
                                      quantiles = c(0.1, 0.5, 0.9) # (not required for checks below; kept for easy extension)
) {
  
  # ---- theme lives INSIDE the function so every returned plot is themed ----
  theme_fancy <- function() {
    ggplot2::theme_minimal(base_family = "Asap Condensed") +
      ggplot2::theme(
        panel.grid.minor = ggplot2::element_blank(),
        legend.position  = "bottom"
      )
  }
  
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
  
  # --- KS plots (returned) ---
  #  - ECDF comparison (treated vs control)
  #  - Density comparison (treated vs control)
  ks_ecdf_plot <- NULL
  ks_density_plot <- NULL
  
  if (!is.null(ks_out)) {
    ks_plot_df <- dplyr::bind_rows(
      dplyr::tibble(value = pre_treated, group = "Treated"),
      dplyr::tibble(value = pre_control, group = "Control")
    )
    
    ks_ecdf_plot <- ggplot2::ggplot(ks_plot_df, ggplot2::aes(x = value, color = group)) +
      ggplot2::stat_ecdf(linewidth = 1) +
      ggplot2::labs(
        title = paste0("KS Test (Pre-Period ECDF): ", outcome),
        subtitle = paste0(
          "KS D = ", signif(unname(ks_out$statistic), 3),
          " | p = ", signif(ks_p, 3)
        ),
        x = outcome,
        y = "Empirical CDF",
        color = "Group"
      ) +
      theme_fancy()
    
    ks_density_plot <- ggplot2::ggplot(ks_plot_df, ggplot2::aes(x = value, fill = group)) +
      ggplot2::geom_density(alpha = 0.35, linewidth = 0.8) +
      ggplot2::labs(
        title = paste0("Pre-Period Distribution: ", outcome),
        subtitle = paste0(
          "KS D = ", signif(unname(ks_out$statistic), 3),
          " | p = ", signif(ks_p, 3)
        ),
        x = outcome,
        y = "Density",
        fill = "Group"
      ) +
      theme_fancy()
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
  
  # --- Slope comparison plot (returned) ---
  # This visualizes treated vs control fitted slopes in the pre-period.
  slope_comparison_plot <- pre_data %>%
    dplyr::mutate(
      group = dplyr::if_else(!!treat_sym == 1, "Treated", "Control"),
      time_num = as.numeric(!!time_sym)
    ) %>%
    dplyr::group_by(!!time_sym, group) %>%
    dplyr::summarize(
      avg_outcome = mean(!!outcome_sym),
      time_num = dplyr::first(time_num),
      .groups = "drop"
    ) %>%
    ggplot2::ggplot(ggplot2::aes(x = !!time_sym, y = avg_outcome, color = group)) +
    ggplot2::geom_point(alpha = 0.8) +
    ggplot2::geom_smooth(method = "lm", se = TRUE, linewidth = 1) +
    ggplot2::labs(
      title = paste0("Pre-Period Slope Comparison: ", outcome),
      subtitle = paste0("Differential slope (time*treat interaction): p = ", signif(slope_p, 3)),
      x     = rlang::as_label(time_sym),
      y     = paste("Average", outcome),
      color = "Group"
    ) +
    theme_fancy()
  
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
      title = paste0("Pre-Period Parallel Trends: ", outcome),
      x     = rlang::as_label(time_sym),
      y     = paste("Average", outcome),
      color = "Group"
    ) +
    theme_fancy()
  
  # ============================================================
  # NEW: Additional parallel trends pre-checks (1–5 + extras)
  # ============================================================
  
  # -------------------------
  # 8) Pre-trend shape, not just slope
  #    (a) Higher-order time interactions via natural splines (pre only)
  # -------------------------
  curvature_model <- NULL
  curvature_p_value <- NA_real_
  curvature_plot <- NULL
  
  # Build a safe pre dataset with numeric time
  pre_for_shape <- pre_data %>%
    dplyr::mutate(
      group = dplyr::if_else(!!treat_sym == 1, "Treated", "Control"),
      time_num = as.numeric(!!time_sym)
    ) %>%
    dplyr::filter(!is.na(time_num))
  
  if (nrow(pre_for_shape) > 0) {
    # Model: outcome ~ ns(time_num, df=k) * factor(treat)
    curvature_formula <- as.formula(
      paste0(outcome, " ~ splines::ns(time_num, df = ", spline_df, ") * factor(", treat, ")")
    )
    curvature_model <- stats::lm(curvature_formula, data = pre_for_shape)
    
    # Joint test for all interaction spline terms (best-effort; depends on term names)
    cn <- names(stats::coef(curvature_model))
    int_terms <- grep(paste0("splines::ns\\(time_num.*\\):factor\\(", treat, "\\)"),
                      cn, value = TRUE)
    
    if (length(int_terms) > 0) {
      # car::linearHypothesis gives an F-test for joint null = 0
      lh <- car::linearHypothesis(curvature_model, int_terms)
      curvature_p_value <- lh[2, "Pr(>F)"]
    }
    
    # Plot: fitted curves by group (with smooth spline-ish fit via loess for visualization)
    curvature_plot <- pre_for_shape %>%
      dplyr::group_by(!!time_sym, group) %>%
      dplyr::summarize(avg_outcome = mean(!!outcome_sym), .groups = "drop") %>%
      ggplot2::ggplot(ggplot2::aes(x = !!time_sym, y = avg_outcome, color = group)) +
      ggplot2::geom_point(alpha = 0.7) +
      ggplot2::geom_smooth(method = "loess", se = TRUE, span = 0.8, linewidth = 1) +
      ggplot2::labs(
        title = paste0("Pre-Period Shape Check (Nonlinearity): ", outcome),
        subtitle = paste0("Spline interaction joint test p = ", signif(curvature_p_value, 3)),
        x = rlang::as_label(time_sym),
        y = paste("Average", outcome),
        color = "Group"
      ) +
      theme_fancy()
  }
  
  # -------------------------
  # 9) Seasonality / calendar alignment checks
  #    (a) Seasonality-adjusted slope test
  # -------------------------
  seasonality_slope_model <- NULL
  seasonality_slope_p_value <- NA_real_
  seasonality_slope_plot <- NULL
  
  # attempt to coerce time to Date for calendar features (safe-ish)
  pre_for_cal <- pre_data %>%
    dplyr::mutate(
      time_date = as.Date(!!time_sym),
      time_num  = as.numeric(!!time_sym),
      group     = dplyr::if_else(!!treat_sym == 1, "Treated", "Control"),
      dow       = as.factor(weekdays(as.Date(!!time_sym))),
      month     = as.factor(format(as.Date(!!time_sym), "%Y-%m"))
    )
  
  # Seasonality adjusted model only if time_date is not all NA
  if (sum(!is.na(pre_for_cal$time_date)) > 0) {
    seasonality_slope_formula <- as.formula(
      paste0(outcome, " ~ time_num * factor(", treat, ") + dow + month")
    )
    seasonality_slope_model <- stats::lm(seasonality_slope_formula, data = pre_for_cal)
    
    ss <- summary(seasonality_slope_model)
    # try to find the time_num:treat interaction
    inter2 <- grep("^time_num:factor\\(", rownames(ss$coefficients), value = TRUE)[1]
    if (!is.na(inter2) && length(inter2) == 1) {
      seasonality_slope_p_value <- ss$coefficients[inter2, "Pr(>|t|)"]
    }
    
    seasonality_slope_plot <- pre_for_cal %>%
      dplyr::group_by(!!time_sym, group) %>%
      dplyr::summarize(avg_outcome = mean(!!outcome_sym), .groups = "drop") %>%
      ggplot2::ggplot(ggplot2::aes(x = !!time_sym, y = avg_outcome, color = group)) +
      ggplot2::geom_line(linewidth = 1) +
      ggplot2::labs(
        title = paste0("Seasonality-Adjusted Parallel Trends View: ", outcome),
        subtitle = paste0("Adjusted slope interaction p = ", signif(seasonality_slope_p_value, 3)),
        x = rlang::as_label(time_sym),
        y = paste("Average", outcome),
        color = "Group"
      ) +
      theme_fancy()
  }
  
  #    (b) Balance on calendar composition (chi-square on dow x group) (pre only)
  calendar_balance_test <- NULL
  calendar_balance_p_value <- NA_real_
  calendar_balance_table <- NULL
  
  if (sum(!is.na(pre_for_cal$time_date)) > 0) {
    calendar_balance_table <- pre_for_cal %>%
      dplyr::filter(!is.na(dow)) %>%
      dplyr::count(group, dow) %>%
      tidyr::pivot_wider(names_from = dow, values_from = n, values_fill = 0)
    
    # build a contingency table for chisq
    tab <- pre_for_cal %>%
      dplyr::filter(!is.na(dow)) %>%
      dplyr::count(group, dow) %>%
      tidyr::pivot_wider(names_from = dow, values_from = n, values_fill = 0)
    
    if (nrow(tab) == 2) {
      mat <- as.matrix(tab %>% dplyr::select(-group))
      rownames(mat) <- tab$group
      calendar_balance_test <- suppressWarnings(stats::chisq.test(mat))
      calendar_balance_p_value <- calendar_balance_test$p.value
    }
  }
  
  # -------------------------
  # 10) Distributional stability over time (not just pooled pre)
  #     Rolling-window distance (KS + Wasserstein) over pre period
  # -------------------------
  rolling_distance_df <- NULL
  rolling_distance_plot <- NULL
  
  # simple 1D Wasserstein distance using quantile interpolation
  wasserstein_1d <- function(x, y, n_q = 200) {
    x <- x[is.finite(x)]
    y <- y[is.finite(y)]
    if (length(x) < 2 || length(y) < 2) return(NA_real_)
    qs <- seq(0, 1, length.out = n_q)
    qx <- stats::quantile(x, probs = qs, na.rm = TRUE, type = 7)
    qy <- stats::quantile(y, probs = qs, na.rm = TRUE, type = 7)
    mean(abs(qx - qy), na.rm = TRUE)
  }
  
  # build time-ordered pre slice with group and time key
  pre_roll <- pre_data %>%
    dplyr::mutate(
      group = dplyr::if_else(!!treat_sym == 1, "Treated", "Control")
    ) %>%
    dplyr::filter(!is.na(!!time_sym))
  
  # roll over unique time points (assumes time is sortable)
  uniq_time <- sort(unique(pre_roll %>% dplyr::pull(!!time_sym)))
  
  if (length(uniq_time) >= max(rolling_window, 2)) {
    roll_rows <- lapply(seq_len(length(uniq_time) - rolling_window + 1), function(i) {
      t_start <- uniq_time[i]
      t_end   <- uniq_time[i + rolling_window - 1]
      t_mid   <- uniq_time[i + floor((rolling_window - 1) / 2)]
      
      wdat <- pre_roll %>%
        dplyr::filter(!!time_sym >= t_start, !!time_sym <= t_end)
      
      x <- wdat %>% dplyr::filter(group == "Treated") %>% dplyr::pull(!!outcome_sym)
      y <- wdat %>% dplyr::filter(group == "Control") %>% dplyr::pull(!!outcome_sym)
      
      ks_d <- NA_real_
      ks_p <- NA_real_
      if (length(x) > 0 && length(y) > 0) {
        kso <- tryCatch(stats::ks.test(x, y), error = function(e) NULL)
        if (!is.null(kso)) {
          ks_d <- unname(kso$statistic)
          ks_p <- kso$p.value
        }
      }
      
      w1 <- wasserstein_1d(x, y)
      
      dplyr::tibble(
        window_start = t_start,
        window_end   = t_end,
        window_mid   = t_mid,
        ks_D         = ks_d,
        ks_p_value   = ks_p,
        wasserstein  = w1
      )
    })
    
    rolling_distance_df <- dplyr::bind_rows(roll_rows)
    
    rolling_distance_plot <- rolling_distance_df %>%
      tidyr::pivot_longer(cols = c("ks_D", "wasserstein"), names_to = "metric", values_to = "value") %>%
      ggplot2::ggplot(ggplot2::aes(x = window_mid, y = value, linetype = metric)) +
      ggplot2::geom_line(linewidth = 1) +
      ggplot2::labs(
        title = paste0("Rolling Pre-Period Distribution Distance: ", outcome),
        subtitle = paste0("Window length = ", rolling_window, " time points"),
        x = rlang::as_label(time_sym),
        y = "Distance",
        linetype = "Metric"
      ) +
      theme_fancy()
  }
  
  # -------------------------
  # 11) Unit-level heterogeneity & composition checks
  #     Unit-specific pre slopes distribution comparison
  # -------------------------
  unit_slope_df <- NULL
  unit_slope_t_test <- NULL
  unit_slope_ks_test <- NULL
  unit_slope_plot <- NULL
  
  pre_unit <- pre_data %>%
    dplyr::mutate(time_num = as.numeric(!!time_sym)) %>%
    dplyr::filter(!is.na(time_num))
  
  if (nrow(pre_unit) > 0) {
    # per-unit slope in pre: outcome ~ time_num
    unit_slope_df <- pre_unit %>%
      dplyr::group_by(!!unit_sym) %>%
      dplyr::summarize(
        treat_value = dplyr::first(!!treat_sym),
        n_points = dplyr::n(),
        slope = {
          if (dplyr::n() >= 2) {
            stats::coef(stats::lm(stats::reformulate("time_num", response = outcome), data = dplyr::cur_data()))[["time_num"]]
          } else {
            NA_real_
          }
        },
        .groups = "drop"
      ) %>%
      dplyr::mutate(group = dplyr::if_else(treat_value == 1, "Treated", "Control"))
    
    # compare slope distributions (treated vs control)
    sl_t <- unit_slope_df %>% dplyr::filter(group == "Treated") %>% dplyr::pull(slope)
    sl_c <- unit_slope_df %>% dplyr::filter(group == "Control") %>% dplyr::pull(slope)
    
    if (sum(is.finite(sl_t)) > 1 && sum(is.finite(sl_c)) > 1) {
      unit_slope_t_test <- stats::t.test(sl_t, sl_c)
      unit_slope_ks_test <- stats::ks.test(sl_t, sl_c)
    }
    
    unit_slope_plot <- unit_slope_df %>%
      dplyr::filter(is.finite(slope)) %>%
      ggplot2::ggplot(ggplot2::aes(x = group, y = slope, fill = group)) +
      ggplot2::geom_violin(trim = FALSE, alpha = 0.35) +
      ggplot2::geom_boxplot(width = 0.15, outlier.alpha = 0.4) +
      ggplot2::labs(
        title = "Unit-Level Pre-Period Slope Distribution",
        subtitle = paste0(
          "t-test p = ",
          if (!is.null(unit_slope_t_test)) signif(unit_slope_t_test$p.value, 3) else "NA",
          " | KS p = ",
          if (!is.null(unit_slope_ks_test)) signif(unit_slope_ks_test$p.value, 3) else "NA"
        ),
        x = "Group",
        y = "Per-unit slope (pre)"
      ) +
      theme_fancy()
  }
  
  # -------------------------
  # 12) Pre-period “outlier influence” checks
  #     (a) Robust slope test (MASS::rlm) + compare to OLS
  # -------------------------
  robust_slope_model <- NULL
  robust_interaction_est <- NA_real_
  robust_interaction_se  <- NA_real_
  robust_interaction_p   <- NA_real_
  
  if (nrow(pre_for_shape) > 0) {
    # same spec as original slope_model
    robust_slope_model <- tryCatch(
      MASS::rlm(slope_formula, data = pre_data, maxit = 100),
      error = function(e) NULL
    )
    
    if (!is.null(robust_slope_model)) {
      rs <- summary(robust_slope_model)
      rn <- rownames(rs$coefficients)
      
      # try to match the original interaction term pattern
      r_inter <- grep("as.numeric\\(", rn, value = TRUE)[1]
      
      if (!is.na(r_inter) && length(r_inter) == 1) {
        robust_interaction_est <- rs$coefficients[r_inter, "Value"]
        robust_interaction_se  <- rs$coefficients[r_inter, "Std. Error"]
        # rlm doesn't provide p-values; use normal approx as a quick diagnostic
        z <- robust_interaction_est / robust_interaction_se
        robust_interaction_p <- 2 * stats::pnorm(-abs(z))
      }
    }
  }
  
  #     (b) Influence diagnostics on slope_model (Cook's D), flag top dates/points
  influence_df <- NULL
  cooks_plot <- NULL
  top_influential <- NULL
  
  if (nrow(pre_data) > 0) {
    cd <- stats::cooks.distance(slope_model)
    
    influence_df <- pre_data %>%
      dplyr::mutate(
        .cooks_d = as.numeric(cd),
        .row_id  = dplyr::row_number(),
        .time    = !!time_sym,
        .group   = dplyr::if_else(!!treat_sym == 1, "Treated", "Control")
      )
    
    top_influential <- influence_df %>%
      dplyr::arrange(dplyr::desc(.cooks_d)) %>%
      dplyr::slice_head(n = 10) %>%
      dplyr::select(.row_id, .time, .group, .cooks_d)
    
    cooks_plot <- influence_df %>%
      ggplot2::ggplot(ggplot2::aes(x = .time, y = .cooks_d, color = .group)) +
      ggplot2::geom_point(alpha = 0.7) +
      ggplot2::geom_smooth(method = "loess", se = FALSE, span = 0.9, linewidth = 1) +
      ggplot2::labs(
        title = "Influence Check (Pre-Period): Cook's Distance for Slope Model",
        subtitle = "Large Cook’s D indicates points strongly influencing the slope interaction",
        x = rlang::as_label(time_sym),
        y = "Cook's distance",
        color = "Group"
      ) +
      theme_fancy()
  }
  
  # -------------------------
  # 13) Pre-only event-study placebo (lead coefficients + joint test)
  # -------------------------
  placebo_event_study_model <- NULL
  placebo_event_study_plot  <- NULL
  placebo_leads_joint_test  <- NULL
  placebo_leads_p_value     <- NA_real_
  
  # We estimate a placebo event study using PRE data only.
  # Note: We do NOT include a full time FE here because rel_time is a deterministic function of time.
  # We include unit FE to soak up level differences; remaining lead coefficients should be ~0 under parallel trends.
  pre_es <- pre_data %>%
    dplyr::mutate(rel_time = as.integer(!!time_sym - intervention_date))
  
  # require at least a few distinct rel_time values
  if (nrow(pre_es) > 0 && dplyr::n_distinct(pre_es$rel_time) >= 3) {
    placebo_fml_str <- paste0(
      outcome,
      " ~ i(rel_time, ", treat, ", ref = -1) | ",
      unit
    )
    
    placebo_event_study_model <- tryCatch(
      fixest::feols(
        fml     = as.formula(placebo_fml_str),
        data    = pre_es,
        cluster = as.formula(paste0("~", unit))
      ),
      error = function(e) NULL
    )
    
    if (!is.null(placebo_event_study_model)) {
      # extract coefficients and keep only *leads* (rel_time < 0) excluding ref
      b  <- stats::coef(placebo_event_study_model)
      nm <- names(b)
      
      # best-effort keep of rel_time terms
      keep <- grepl("rel_time", nm)
      est_df <- dplyr::tibble(term = nm[keep], estimate = unname(b[keep])) %>%
        dplyr::mutate(
          rel_time = suppressWarnings(as.integer(gsub(".*?(-?\\d+).*", "\\1", term)))
        ) %>%
        dplyr::filter(!is.na(rel_time)) %>%
        dplyr::filter(rel_time < 0)
      
      # joint test that all lead coefficients = 0
      # fixest::wald can test a vector of coefficients (by name)
      if (nrow(est_df) > 0) {
        coef_names <- est_df$term
        placebo_leads_joint_test <- tryCatch(
          fixest::wald(placebo_event_study_model, coef_names),
          error = function(e) NULL
        )
        if (!is.null(placebo_leads_joint_test) && "p.value" %in% names(placebo_leads_joint_test)) {
          placebo_leads_p_value <- placebo_leads_joint_test$p.value
        }
      }
      
      # build CI if possible
      ci <- tryCatch(stats::confint(placebo_event_study_model), error = function(e) NULL)
      if (!is.null(ci) && nrow(est_df) > 0) {
        ci_df <- as.data.frame(ci)
        ci_df$term <- rownames(ci_df)
        ci_df <- dplyr::as_tibble(ci_df) %>% dplyr::rename(conf_low = 1, conf_high = 2)
        est_df <- est_df %>% dplyr::left_join(ci_df, by = "term")
      } else if (nrow(est_df) > 0) {
        est_df <- est_df %>% dplyr::mutate(conf_low = NA_real_, conf_high = NA_real_)
      }
      
      if (nrow(est_df) > 0) {
        placebo_event_study_plot <- ggplot2::ggplot(est_df, ggplot2::aes(x = rel_time, y = estimate)) +
          ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
          ggplot2::geom_point(size = 2) +
          ggplot2::geom_line(linewidth = 0.8) +
          ggplot2::geom_errorbar(
            ggplot2::aes(ymin = conf_low, ymax = conf_high),
            width = 0.15,
            alpha = 0.8
          ) +
          ggplot2::labs(
            title = paste0("Pre-Only Placebo Event Study (Leads): ", outcome),
            subtitle = paste0("Joint test of lead coefficients p = ", signif(placebo_leads_p_value, 3)),
            x = "Relative time (lead periods, < 0)",
            y = "Estimated effect"
          ) +
          theme_fancy()
      }
    }
  }
  
  # -------------------------
  # 6) Optional: Event-study model (full sample)
  # -------------------------
  event_study_model <- NULL
  event_study_plot  <- NULL
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
    
    # --- Event study plot (returned) ---
    # Robust manual extraction + ggplot so we can theme and title reliably.
    coefs <- stats::coef(event_study_model)
    nm <- names(coefs)
    
    # Keep only i(rel_time, treat, ...) terms (best-effort across naming variants)
    keep <- grepl("rel_time", nm)
    
    est_df <- dplyr::tibble(term = nm[keep], estimate = unname(coefs[keep])) %>%
      dplyr::mutate(
        rel_time = suppressWarnings(as.integer(gsub(".*?(-?\\d+).*", "\\1", term)))
      ) %>%
      dplyr::filter(!is.na(rel_time))
    
    ci <- tryCatch(stats::confint(event_study_model), error = function(e) NULL)
    if (!is.null(ci)) {
      ci_df <- as.data.frame(ci)
      ci_df$term <- rownames(ci_df)
      ci_df <- dplyr::as_tibble(ci_df) %>%
        dplyr::rename(conf_low = 1, conf_high = 2)
      
      est_df <- est_df %>% dplyr::left_join(ci_df, by = "term")
    } else {
      est_df <- est_df %>% dplyr::mutate(conf_low = NA_real_, conf_high = NA_real_)
    }
    
    event_study_plot <- ggplot2::ggplot(est_df, ggplot2::aes(x = rel_time, y = estimate)) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
      ggplot2::geom_point(size = 2) +
      ggplot2::geom_line(linewidth = 0.8) +
      ggplot2::geom_errorbar(
        ggplot2::aes(ymin = conf_low, ymax = conf_high),
        width = 0.15,
        alpha = 0.8
      ) +
      ggplot2::labs(
        title = paste0("Event Study: ", outcome),
        subtitle = "Coefficients from i(rel_time, treat, ref = -1)",
        x = "Relative time (periods from intervention)",
        y = "Estimated effect"
      ) +
      theme_fancy()
  }
  
  # -------------------------
  # 7) Bundle everything into a list
  # -------------------------
  results <- list(
    pre_data                   = pre_data,
    mean_test                  = mean_test,
    levene_test                = lev_out,
    levene_p_value             = levene_p,
    fligner_test               = fligner_out,
    fligner_p_value            = fligner_p,
    ks_test                    = ks_out,
    ks_p_value                 = ks_p,
    slope_model                = slope_model,
    slope_p_value              = slope_p,
    pre_trends_plot            = pre_trends_plot,
    # returned plots requested (original + earlier additions)
    ks_ecdf_plot               = ks_ecdf_plot,
    ks_density_plot            = ks_density_plot,
    slope_comparison_plot      = slope_comparison_plot,
    event_study_model          = event_study_model,
    event_study_plot           = event_study_plot,
    # ---- new additions (1–5 + extras) ----
    curvature_model            = curvature_model,
    curvature_p_value          = curvature_p_value,
    curvature_plot             = curvature_plot,
    seasonality_slope_model    = seasonality_slope_model,
    seasonality_slope_p_value  = seasonality_slope_p_value,
    seasonality_slope_plot     = seasonality_slope_plot,
    calendar_balance_test      = calendar_balance_test,
    calendar_balance_p_value   = calendar_balance_p_value,
    calendar_balance_table     = calendar_balance_table,
    rolling_distance_df        = rolling_distance_df,
    rolling_distance_plot      = rolling_distance_plot,
    unit_slope_df              = unit_slope_df,
    unit_slope_t_test          = unit_slope_t_test,
    unit_slope_ks_test         = unit_slope_ks_test,
    unit_slope_plot            = unit_slope_plot,
    robust_slope_model         = robust_slope_model,
    robust_interaction_est     = robust_interaction_est,
    robust_interaction_se      = robust_interaction_se,
    robust_interaction_p       = robust_interaction_p,
    influence_df               = influence_df,
    cooks_plot                 = cooks_plot,
    top_influential            = top_influential,
    placebo_event_study_model  = placebo_event_study_model,
    placebo_leads_joint_test   = placebo_leads_joint_test,
    placebo_leads_p_value      = placebo_leads_p_value,
    placebo_event_study_plot   = placebo_event_study_plot
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
  cat("Nonlinear shape check (spline interaction): p-value = ",
      signif(curvature_p_value, 3), "\n", sep = "")
  cat("Seasonality-adjusted slope test: p-value = ",
      signif(seasonality_slope_p_value, 3), "\n", sep = "")
  cat("Calendar composition balance (chi-square over DOW): p-value = ",
      signif(calendar_balance_p_value, 3), "\n", sep = "")
  cat("Robust slope interaction (rlm, normal approx): p-value = ",
      signif(robust_interaction_p, 3), "\n", sep = "")
  cat("Pre-only placebo event study (lead joint test): p-value = ",
      signif(placebo_leads_p_value, 3), "\n", sep = "")
  if (do_event_study) {
    cat("Event study model estimated: yes (use results$event_study_plot)\n")
  } else {
    cat("Event study model estimated: no (set do_event_study = TRUE to add)\n")
  }
  cat("=================================\n")
  
  invisible(results)
}

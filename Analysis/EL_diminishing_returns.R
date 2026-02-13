library(dplyr)
library(tidyr)
library(ggplot2)
library(mgcv)
library(purrr)
library(scales)

###### How to choose option and interpret results ######
# If using method = "pct_of_max" (default), the line marks the smallest total where predicted retention ≥ 90% of that curve’s maximum.
# Ex: Maximum predicted 7d retention might be ~48%, 90% of that is ~43%. At 23 questions answered, predicted retention first reaches 43%.
# After ~23 questions, you’ve captured most of the achievable retention lift.
# If using method = "marginal" the line marks the first point where the marginal retention gain per additional question falls below your threshold.
# So if it says 31 (<0.1pp/action) that means that at ~31 questions each additional question increases retention by less than 0.1 percentage points.
# Interpretation: beyond ~31 questions, effort per question is barely moving retention.
# The label is positioned at x = threshold_total, y = predicted retention at that point.

# Definition: Smallest x where predicted retention reaches pct × the curve’s max (e.g., 95%).
# Interpretation: “How much activity gets you most of the achievable retention level on this curve?”
# Pros
#   Easy to explain (“we’re at 95% of plateau”).
#   Stable when there’s a real plateau.
#   Good for describing the upper bound region.
# Cons
#   If the curve doesn’t truly plateau (keeps creeping upward), the “max” is mostly determined by the far-right tail → threshold can land very high.
#   Sensitive to the chosen max range (max_total) and sparse tail behavior.
#   Often yields a “power-user” threshold that’s not a practical KPI.
# Use it when you have a clear saturation/plateau and want a “near-plateau” benchmark, you’re describing “how far to go for almost all the gain”
# especially for mature funnels or power-user programs, or you’re okay with thresholds that might be high if the curve keeps rising.
# Ex: “Users who answer ~25 flashcards reach ~45% 7d retention.”

plot_diminishing_returns_multi <- function(df,
                                           total_col = flashcards_questions_answered,
                                           retention_cols = c("7d_retained", "28d_retained", "90d_retained"),
                                           n_bins = 50,
                                           min_bin_n = 20,
                                           max_total = NULL,
                                           # thresholding:
                                           method = c("pct_of_max", "marginal"),
                                           pct = 0.90,                 # for pct_of_max
                                           marginal_pp = 0.001,        # for marginal: retention per +1 action (0.001 = 0.1pp)
                                           grid_n = 400,
                                           show_points = TRUE) {
  
  method <- match.arg(method)
  
  dat <- df %>%
    transmute(
      total = as.numeric({{ total_col }}),
      across(all_of(retention_cols), ~ as.integer(as.logical(.x)))
    ) %>%
    filter(!is.na(total), total >= 0) %>%
    drop_na(all_of(retention_cols))
  
  if (!is.null(max_total)) dat <- dat %>% filter(total <= max_total)
  
  probs <- unique(quantile(dat$total, probs = seq(0, 1, length.out = n_bins + 1), na.rm = TRUE))
  if (length(probs) < 6) stop("Not enough unique totals to bin reliably. Try fewer bins (e.g., n_bins = 15).")
  
  binned <- dat %>%
    mutate(bin = cut(total, breaks = probs, include.lowest = TRUE, right = TRUE)) %>%
    group_by(bin) %>%
    summarise(
      total_mid = median(total),
      n = n(),
      across(all_of(retention_cols), ~ mean(.x), .names = "{.col}"),
      .groups = "drop"
    ) %>%
    filter(n >= min_bin_n) %>%
    pivot_longer(
      cols = all_of(retention_cols),
      names_to = "window",
      values_to = "retention"
    ) %>%
    mutate(window = factor(window, levels = retention_cols))
  
  # Fit per-window smooth + prediction grid + compute threshold
  model_data <- binned %>%
    group_by(window) %>%
    nest()
  
  results <- model_data %>%
    mutate(
      fit = map(data, ~{
        k_val <- min(12, nrow(.x) - 1)
        mgcv::gam(
          retention ~ s(log1p(total_mid), k = k_val),
          data = .x,
          family = quasibinomial()
        )
      }),
      curve = map2(data, fit, ~{
        grid <- tibble(
          total_mid = seq(min(.x$total_mid), max(.x$total_mid), length.out = grid_n)
        )
        grid %>%
          mutate(pred = predict(.y, newdata = ., type = "response"))
      }),
      # numeric slope in *raw total* units: d(pred)/d(total)
      curve_slope = map(curve, ~{
        .x %>%
          arrange(total_mid) %>%
          mutate(slope = c(NA_real_, diff(pred) / diff(total_mid)))
      }),
      threshold_total = pmap_dbl(list(curve, curve_slope, window), function(curve, curve_slope, window){
        if (method == "pct_of_max") {
          target <- pct * max(curve$pred, na.rm = TRUE)
          # first total where pred crosses target
          idx <- which(curve$pred >= target)[1]
          if (is.na(idx)) return(NA_real_)
          return(curve$total_mid[idx])
        } else {
          # first total where marginal gain per +1 action drops below marginal_pp
          idx <- which(curve_slope$slope <= marginal_pp)[1]
          if (is.na(idx)) return(NA_real_)
          return(curve_slope$total_mid[idx])
        }
      }),
      threshold_pred = map2_dbl(curve, threshold_total, ~{
        if (is.na(.y)) return(NA_real_)
        # interpolate predicted retention at threshold
        approx(x = .x$total_mid, y = .x$pred, xout = .y, rule = 2)$y
      })
    )
  
  curve_df <- results %>% select(window, curve) %>% unnest(curve)
  
  ann_df <- results %>%
    transmute(
      window,
      threshold_total,
      threshold_pred,
      label = case_when(
        is.na(threshold_total) ~ "no threshold",
        method == "pct_of_max" ~ paste0(round(threshold_total), " (", round(100 * pct), "% max)"),
        TRUE ~ paste0(round(threshold_total), " (<", round(marginal_pp * 100, 2), "pp/action)")
      )
    ) %>%
    filter(!is.na(threshold_total), !is.na(threshold_pred))
  
  # Prettier legend labels
  window_labels <- setNames(
    nm = retention_cols,
    object = gsub("_retained$", " retained", retention_cols)
  )
  
  p <- ggplot() +
    { if (show_points)
      geom_point(
        data = binned,
        aes(x = .data$total_mid, y = .data$retention, color = .data$window),
        alpha = 0.30
      ) } +
    geom_line(
      data = curve_df,
      aes(x = .data$total_mid, y = .data$pred, color = .data$window),
      linewidth = 1.1
    ) +
    # vertical lines at thresholds (per window)
    geom_vline(
      data = ann_df,
      aes(xintercept = .data$threshold_total, color = .data$window),
      linetype = "dashed",
      alpha = 0.8
    ) +
    # labels near the curve at the threshold point
    geom_text(
      data = ann_df,
      aes(x = .data$threshold_total, y = .data$threshold_pred, label = .data$label, color = .data$window),
      vjust = -0.8,
      size = 3.3,
      show.legend = FALSE
    ) +
    scale_color_discrete(name = "Retention", labels = window_labels) +
    scale_x_continuous(
      trans = "log1p",
      breaks = c(0, 1, 2, 5, 10, 20, 50, 100, 200, 500, 800),
      labels = label_number()
    ) +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    labs(
      x = "Flashcards questions answered (log scale: log1p)",
      y = "Retention rate",
      title = "Diminishing returns of flashcards questions answered on retention",
      subtitle = if (method == "pct_of_max") {
        paste0("Threshold = smallest total reaching ", round(100 * pct), "% of max predicted retention (per curve)")
      } else {
        paste0("Threshold = smallest total where marginal gain < ", marginal_pp * 100, " percentage points per action (per curve)")
      }
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom")
  
  return(p)
}

# Example:
plot_diminishing_returns_multi(el_data%>%filter(flashcards_questions_answered>0))                           # default: 90% of max
plot_diminishing_returns_multi(el_data, pct = 0.80)               # stricter plateau
plot_diminishing_returns_multi(el_data, method="marginal",marginal_pp = 0.0005)             # 0.05pp per action

# Plot diminshing returns inflection to find 

plot_diminishing_returns_inflection <- function(df,
                                                total_col = flashcards_questions_answered,
                                                retention_cols = c("7d_retained", "28d_retained", "90d_retained"),
                                                n_bins = 50,
                                                min_bin_n = 20,
                                                max_total = NULL,
                                                slope_fraction = 0.25,
                                                grid_n = 500) {
  
  dat <- df %>%
    transmute(
      total = as.numeric({{ total_col }}),
      across(all_of(retention_cols), ~ as.integer(as.logical(.x)))
    ) %>%
    filter(!is.na(total), total >= 0) %>%
    drop_na(all_of(retention_cols))
  
  if (!is.null(max_total)) dat <- dat %>% filter(total <= max_total)
  
  probs <- unique(quantile(dat$total, probs = seq(0, 1, length.out = n_bins + 1), na.rm = TRUE))
  
  binned <- dat %>%
    mutate(bin = cut(total, breaks = probs, include.lowest = TRUE)) %>%
    group_by(bin) %>%
    summarise(
      total_mid = median(total),
      n = n(),
      across(all_of(retention_cols), ~ mean(.x), .names = "{.col}"),
      .groups = "drop"
    ) %>%
    filter(n >= min_bin_n) %>%
    pivot_longer(
      cols = all_of(retention_cols),
      names_to = "window",
      values_to = "retention"
    ) %>%
    mutate(window = factor(window, levels = retention_cols))
  
  model_data <- binned %>%
    group_by(window) %>%
    nest()
  
  results <- model_data %>%
    mutate(
      fit = map(data, ~{
        mgcv::gam(
          retention ~ s(log1p(total_mid), k = min(12, nrow(.x)-1)),
          data = .x,
          family = quasibinomial()
        )
      }),
      curve = map2(data, fit, ~{
        grid <- tibble(
          total_mid = seq(min(.x$total_mid), max(.x$total_mid), length.out = grid_n)
        )
        grid %>%
          mutate(pred = predict(.y, newdata = ., type = "response"))
      }),
      curve = map(curve, ~{
        .x %>%
          arrange(total_mid) %>%
          mutate(
            slope = c(NA, diff(pred)/diff(total_mid))
          )
      }),
      threshold = map_dbl(curve, ~{
        initial_slope <- max(.x$slope, na.rm = TRUE)
        cutoff <- slope_fraction * initial_slope
        idx <- which(.x$slope <= cutoff)[1]
        if (is.na(idx)) return(NA_real_)
        .x$total_mid[idx]
      })
    )
  
  curve_df <- results %>%
    select(window, curve) %>%
    unnest(curve)
  
  threshold_df <- results %>%
    transmute(
      window,
      threshold
    )
  
  window_labels <- setNames(
    gsub("_retained$", " retained", retention_cols),
    retention_cols
  )
  
  ggplot() +
    geom_line(
      data = curve_df,
      aes(x = total_mid, y = pred, color = window),
      linewidth = 1.2
    ) +
    geom_vline(
      data = threshold_df,
      aes(xintercept = threshold, color = window),
      linetype = "dashed"
    ) +
    geom_text(
      data = threshold_df,
      aes(x = threshold, y = max(curve_df$pred),
          label = paste0("~", round(threshold))),
      angle = 90,
      vjust = -0.4,
      show.legend = FALSE
    ) +
    scale_color_discrete(name = "Retention", labels = window_labels) +
    scale_x_continuous(
      trans = "log1p",
      labels = label_number()
    ) +
    scale_y_continuous(labels = percent_format()) +
    labs(
      x = "Flashcards questions answered (log scale)",
      y = "Retention rate",
      title = "Inflection / Early-Saturation Threshold",
      subtitle = paste0("Threshold where marginal gain falls below ",
                        slope_fraction*100, "% of initial slope")
    ) +
    theme_minimal(base_size = 13) +
    theme(legend.position = "bottom")
}

plot_diminishing_returns_inflection(el_data%>%filter(flashcards_questions_answered>0), slope_fraction = 0.20)

plot_diminishing_returns_inflection(el_data%>%filter(flashcards_questions_answered>0), slope_fraction = 0.30)

# Interpretation: “Where do additional actions stop being as efficient as they were early on?”
# Pros
#   Focuses on the steep part of the curve (activation zone), not the extreme tail.
#   Better aligned with “diminishing returns” as a rate concept.
#   More practical for KPIs because it tends to land in a reachable range.
# Cons
#   Requires choices: what counts as “early” (init_quantile) and what decay ratio (slope_fraction).
#   Can be noisy if the curve is wiggly;
#   Less intuitive to explain than “95% of max” unless you translate it (“returns drop by ~65%”).
# Use it when you want an activation KPI (“good enough engagement”) rather than a near-maximum benchmark.
plot_diminishing_returns_inflection <- function(df,
                                                total_col = flashcards_questions_answered,
                                                retention_cols = c("7d_retained", "28d_retained", "90d_retained"),
                                                n_bins = 50,
                                                min_bin_n = 20,
                                                max_total = NULL,
                                                slope_fraction = 0.25,
                                                # how to define "initial" slope:
                                                init_quantile = 0.10,   # use first 10% of x-range (excluding 0)
                                                min_x_for_threshold = 1,# don't allow threshold at 0
                                                consec = 5,             # require slope below cutoff for consec points
                                                grid_n = 600,
                                                show_points = TRUE,
                                                shade_after_threshold = TRUE) {
  
  dat <- df %>%
    transmute(
      total = as.numeric({{ total_col }}),
      across(all_of(retention_cols), ~ as.integer(as.logical(.x)))
    ) %>%
    filter(!is.na(total), total >= 0) %>%
    drop_na(all_of(retention_cols))
  
  if (!is.null(max_total)) dat <- dat %>% filter(total <= max_total)
  
  probs <- unique(quantile(dat$total, probs = seq(0, 1, length.out = n_bins + 1), na.rm = TRUE))
  if (length(probs) < 6) stop("Not enough unique totals to bin reliably; try fewer bins (e.g., 15).")
  
  binned <- dat %>%
    mutate(bin = cut(total, breaks = probs, include.lowest = TRUE, right = TRUE)) %>%
    group_by(bin) %>%
    summarise(
      total_mid = median(total),
      n = n(),
      across(all_of(retention_cols), ~ mean(.x), .names = "{.col}"),
      .groups = "drop"
    ) %>%
    filter(n >= min_bin_n) %>%
    pivot_longer(cols = all_of(retention_cols), names_to = "window", values_to = "retention") %>%
    mutate(window = factor(window, levels = retention_cols))
  
  model_data <- binned %>% group_by(window) %>% nest()
  
  results <- model_data %>%
    mutate(
      fit = map(data, ~{
        k_val <- min(12, nrow(.x) - 1)
        mgcv::gam(retention ~ s(log1p(total_mid), k = k_val),
                  data = .x, family = quasibinomial())
      }),
      curve = map2(data, fit, ~{
        grid <- tibble(total_mid = seq(min(.x$total_mid), max(.x$total_mid), length.out = grid_n))
        grid %>% mutate(pred = predict(.y, newdata = grid, type = "response"))
      }),
      curve = map(curve, ~{
        .x %>%
          arrange(total_mid) %>%
          mutate(
            slope = c(NA_real_, diff(pred) / diff(total_mid)),
            slope_pos = pmax(slope, 0)   # we care about diminishing *positive* gains
          )
      }),
      init_slope = map_dbl(curve, ~{
        # define initial slope as median positive slope in early x-range excluding very small x
        x0 <- max(min_x_for_threshold, quantile(.x$total_mid, probs = 0.02, na.rm = TRUE))
        x1 <- quantile(.x$total_mid, probs = init_quantile, na.rm = TRUE)
        early <- .x %>% filter(total_mid >= x0, total_mid <= x1, is.finite(slope_pos))
        if (nrow(early) == 0) return(NA_real_)
        median(early$slope_pos, na.rm = TRUE)
      }),
      cutoff = init_slope * slope_fraction,
      threshold_total = pmap_dbl(list(curve, cutoff), function(curve, cutoff) {
        if (!is.finite(cutoff) || cutoff <= 0) return(NA_real_)
        
        dfc <- curve %>%
          filter(total_mid >= min_x_for_threshold, is.finite(slope_pos)) %>%
          mutate(below = slope_pos <= cutoff)
        
        if (nrow(dfc) < consec) return(NA_real_)
        
        # first index where we have `consec` consecutive points below cutoff
        run_ok <- sapply(seq_len(nrow(dfc) - consec + 1),
                         function(i) all(dfc$below[i:(i + consec - 1)]))
        idx <- which(run_ok)[1]
        if (is.na(idx)) return(NA_real_)
        dfc$total_mid[idx]
      }),
      threshold_pred = map2_dbl(curve, threshold_total, ~{
        if (is.na(.y)) return(NA_real_)
        approx(x = .x$total_mid, y = .x$pred, xout = .y, rule = 2)$y
      })
    )
  
  curve_df <- results %>% select(window, curve) %>% unnest(curve)
  
  ann_df <- results %>%
    transmute(
      window,
      threshold_total,
      threshold_pred,
      label = ifelse(is.na(threshold_total), "no threshold",
                     paste0(round(threshold_total), " (", round(100 * slope_fraction), "% early slope)"))
    ) %>%
    filter(!is.na(threshold_total), !is.na(threshold_pred))
  
  # shading regions AFTER threshold (per window)
  shade_df <- ann_df %>%
    transmute(window, xmin = threshold_total, xmax = Inf)
  
  window_labels <- setNames(gsub("_retained$", " retained", retention_cols), retention_cols)
  
  gg <- ggplot() +
    { if (shade_after_threshold && nrow(shade_df) > 0)
      geom_rect(
        data = shade_df,
        aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = window),
        alpha = 0.08,
        inherit.aes = FALSE
      ) } +
    { if (show_points)
      geom_point(
        data = binned,
        aes(x = total_mid, y = retention, color = window),
        alpha = 0.30
      ) } +
    geom_line(
      data = curve_df,
      aes(x = total_mid, y = pred, color = window),
      linewidth = 1.2
    ) +
    geom_vline(
      data = ann_df,
      aes(xintercept = threshold_total, color = window),
      linetype = "dashed",
      alpha = 0.9
    ) +
    geom_text(
      data = ann_df,
      aes(x = threshold_total, y = threshold_pred, label = label, color = window),
      vjust = -0.8,
      size = 3.3,
      show.legend = FALSE
    ) +
    scale_color_discrete(name = "Retention", labels = window_labels) +
    scale_fill_discrete(name = "Retention", labels = window_labels, guide = "none") +
    scale_x_continuous(
      trans = "log1p",
      breaks = c(0, 1, 2, 5, 10, 20, 50, 100, 200, 500, 800),
      labels = label_number()
    ) +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    labs(
      x = "Flashcards questions answered (log scale: log1p)",
      y = "Retention rate",
      title = "Early-saturation threshold from marginal effect decay",
      subtitle = paste0(
        "Threshold = first x where marginal gain stays ≤ ",
        round(100 * slope_fraction), "% of early marginal gain (early = first ",
        round(100 * init_quantile), "% of x-range)"
      )
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom")
  
  gg
}

# Example:
plot_diminishing_returns_inflection(el_data, 
                                    slope_fraction = 0.25, 
                                    init_quantile = 0.10, 
                                    min_x_for_threshold = 1)

# Each marked point represents the smallest engagement threshold (T) at which the marginal improvement in retention meaningfully decays relative 
# to the early part of the curve. In other words, beyond this threshold, raising the bar further yields disproportionately smaller gains in 
# retention among qualifying users.

plot_threshold_tradeoff_decision <- function(df,
                                             total_col = flashcards_questions_answered,
                                             retention_cols = c("retained7d", "retained28d", "retained90d"),
                                             n_thresholds = 200,
                                             max_total = NULL,
                                             min_above_n = 200,
                                             x_max = 50,
                                             x_log1p = FALSE,
                                             # ---- diminishing-returns marker params ----
                                             dr_fraction = 0.25,     # "diminishing" = marginal gain falls below 25% of early gain
                                             early_quantile = 0.15,  # early region = first 15% of thresholds (after filtering)
                                             consec = 5,             # require condition holds for N consecutive steps
                                             show_dr_vline = TRUE) {
  
  dat <- df %>%
    transmute(
      total = as.numeric({{ total_col }}),
      across(all_of(retention_cols), ~ as.integer(as.logical(.x)))
    ) %>%
    filter(!is.na(total), total >= 0) %>%
    drop_na(all_of(retention_cols))
  
  if (!is.null(max_total)) dat <- dat %>% filter(total <= max_total)
  
  n_total <- nrow(dat)
  
  thresholds <- unique(round(
    quantile(dat$total, probs = seq(0, 0.99, length.out = n_thresholds), na.rm = TRUE)
  )) |> sort()
  
  thresholds <- thresholds[thresholds >= 1]
  if (!is.null(x_max)) thresholds <- thresholds[thresholds <= x_max]
  
  # Retention curves (conditional on being above threshold)
  retention_df <- lapply(retention_cols, function(window) {
    lapply(thresholds, function(T) {
      above <- dat$total >= T
      n_above <- sum(above)
      if (n_above < min_above_n) return(NULL)
      
      tibble(
        window = window,
        threshold = T,
        retained_rate = mean(dat[[window]][above], na.rm = TRUE),
        n_above = n_above
      )
    }) |> bind_rows()
  }) |> bind_rows()
  
  if (nrow(retention_df) == 0) {
    stop("No thresholds meet min_above_n. Lower min_above_n or increase x_max / data size.")
  }
  
  retention_df <- retention_df %>%
    mutate(
      window = factor(window, levels = retention_cols),
      qualifying_rate = n_above / n_total
    )
  
  # Qualifying curve (same for all windows; use unique thresholds)
  qualify_df <- retention_df %>%
    group_by(threshold) %>%
    summarise(qualifying_rate = max(qualifying_rate), .groups = "drop") %>%
    arrange(threshold)
  
  # ----- Compute diminishing-returns threshold per window -----
  dr_df <- retention_df %>%
    group_by(window) %>%
    arrange(threshold) %>%
    mutate(
      d_ret = c(NA_real_, diff(retained_rate)),   # marginal gain per step in threshold grid
      idx = row_number()
    ) %>%
    filter(is.finite(d_ret)) %>%
    group_modify(~{
      dfw <- .x
      
      # Define "early" marginal gain as median of early portion of the curve
      early_end <- max(2, floor(nrow(dfw) * early_quantile))
      early_gain <- median(dfw$d_ret[1:early_end], na.rm = TRUE)
      
      # If early_gain is tiny/negative, the curve doesn't have a meaningful "increasing" region
      if (!is.finite(early_gain) || early_gain <= 0) {
        return(tibble(dr_threshold = NA_real_, dr_retained = NA_real_, early_gain = early_gain, cutoff = NA_real_))
      }
      
      cutoff <- dr_fraction * early_gain
      
      ok <- dfw$d_ret <= cutoff
      
      # Require `consec` consecutive ok's
      if (length(ok) < consec) {
        return(tibble(dr_threshold = NA_real_, dr_retained = NA_real_, early_gain = early_gain, cutoff = cutoff))
      }
      
      run_ok <- sapply(seq_len(length(ok) - consec + 1),
                       function(i) all(ok[i:(i + consec - 1)]))
      j <- which(run_ok)[1]
      
      if (is.na(j)) {
        tibble(dr_threshold = NA_real_, dr_retained = NA_real_, early_gain = early_gain, cutoff = cutoff)
      } else {
        tibble(
          dr_threshold = dfw$threshold[j],
          dr_retained  = dfw$retained_rate[j],
          early_gain = early_gain,
          cutoff = cutoff
        )
      }
    }) %>%
    ungroup()
  
  # ----- Secondary axis scaling for qualifying rate -----
  y1_min <- min(retention_df$retained_rate, na.rm = TRUE)
  y1_max <- max(retention_df$retained_rate, na.rm = TRUE)
  if (y1_max - y1_min < 1e-6) { y1_min <- 0; y1_max <- 1 }
  
  q_to_y <- function(q) y1_min + q * (y1_max - y1_min)
  y_to_q <- function(y) (y - y1_min) / (y1_max - y1_min)
  
  window_labels <- setNames(gsub("_retained$", " retained", retention_cols), retention_cols)
  
  p <- ggplot() +
    # qualifying (secondary axis)
    geom_line(
      data = qualify_df,
      aes(x = threshold, y = q_to_y(qualifying_rate)),
      linewidth = 1.0,
      alpha = 0.6
    ) +
    # retention curves
    geom_line(
      data = retention_df,
      aes(x = threshold, y = retained_rate, color = window),
      linewidth = 1.2
    ) +
    geom_point(
      data = retention_df,
      aes(x = threshold, y = retained_rate, color = window),
      alpha = 0.25
    ) +
    # diminishing-returns markers (points)
    geom_point(
      data = dr_df %>% filter(!is.na(dr_threshold)),
      aes(x = dr_threshold, y = dr_retained, color = window),
      size = 3.2
    ) +
    geom_text(
      data = dr_df %>% filter(!is.na(dr_threshold)),
      aes(x = dr_threshold, y = dr_retained,
          label = paste0("DR @ ", dr_threshold),
          color = window),
      vjust = -0.9,
      size = 3.2,
      show.legend = FALSE
    ) +
    { if (show_dr_vline)
      geom_vline(
        data = dr_df %>% filter(!is.na(dr_threshold)),
        aes(xintercept = dr_threshold, color = window),
        linetype = "dashed",
        alpha = 0.7
      )
    } +
    scale_color_discrete(name = "Retention", labels = window_labels) +
    scale_y_continuous(
      labels = percent_format(accuracy = 1),
      name = "% retained among users with total ≥ T",
      sec.axis = sec_axis(~ y_to_q(.),
                          name = "% of users qualifying (total ≥ T)",
                          labels = percent_format(accuracy = 1))
    ) +
    labs(
      x = "Engagement Threshold, Flashcards questions answered",
      title = "Engagement Threshold Trade-off",
      subtitle = paste0(
        "Marginal retention gain drops to ≤ ",
        round(100 * dr_fraction), "% of early marginal gain (early = first ",
        round(100 * early_quantile), "% of thresholds, sustained for ", consec, " steps)."
      )
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom")
  
  if (x_log1p) {
    p <- p + scale_x_continuous(trans = "log1p", labels = label_number())
  } else {
    p <- p + scale_x_continuous(labels = label_number())
  }
  
  p
}

# Example:
plot_threshold_tradeoff_decision(
  el_data%>%
    filter(flashcards_questions_answered>0),
  x_max = 100, # maximum limit
  dr_fraction = 0.25, # earliest threshold where the marginal gain falls to ≤ 25% of its early marginal gain. Lower value marker happens later
  early_quantile = 0.15, # use the first 15% of evaluated thresholds to compute the baseline marginal gain.
  consec = 3 # diminishing condition must hold for 3 consecutive thresholds.
)

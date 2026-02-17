# ============================================================
# DMA/City feature builder + exclusions + ranking (GENERIC)
# - No dma_code anywhere
# - Metric-generic outputs: <metric>_mean, <metric>_sd, <metric>_cv, <metric>_trend
# - College proxy (Change #2): regression-based "back-to-school lift"
#     y ~ fall_indicator + linear_time_trend
#   college_idx = beta_fall / mean(y)   (normalized lift)
# ============================================================

library(dplyr)
library(lubridate)
library(stringr)
library(rlang)
library(purrr)

# -----------------------------
# Aggregate to weekly (generic)
# -----------------------------
to_weekly <- function(df,
                      week_start = 1,
                      date_col = "date",
                      group_cols = c("dma_name"),
                      value_cols = NULL) {
  
  stopifnot(date_col %in% names(df))
  stopifnot(all(group_cols %in% names(df)))
  
  # Default: sum all numeric columns not in date/group cols
  if (is.null(value_cols)) {
    value_cols <- df %>% select(where(is.numeric)) %>% names()
    value_cols <- setdiff(value_cols, c(group_cols, date_col))
  }
  stopifnot(length(value_cols) > 0, all(value_cols %in% names(df)))
  
  df %>%
    mutate(
      .date = as.Date(.data[[date_col]]),
      week  = floor_date(.date, unit = "week", week_start = week_start)
    ) %>%
    group_by(across(all_of(group_cols)), week) %>%
    summarise(
      across(all_of(value_cols), ~ sum(.x, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    rename(date = week)
}

# -----------------------------
# Build baseline features + regression-based college index
# (Change #2)
# -----------------------------
build_dma_features <- function(df,
                               pre_start,
                               pre_end,
                               dma_col = "dma_name",
                               date_col = "date",
                               value_cols,
                               weekly = TRUE,
                               week_start = 1,
                               min_pre_periods = 6,
                               # College index settings
                               college_metric = NULL,            # which metric to use for the index
                               fall_months = c(8, 9, 10),         # Aug–Oct tends to work better than Sep–Oct
                               min_non_na_obs = 10,               # guardrail for regression stability
                               add_metric_specific_college_col = TRUE) {
  
  required <- c(dma_col, date_col, value_cols)
  missing  <- setdiff(required, names(df))
  if (length(missing) > 0) {
    stop(
      "Missing columns in df: ", paste0(missing, collapse = ", "),
      "\nAvailable columns: ", paste0(names(df), collapse = ", ")
    )
  }
  
  if (is.null(college_metric)) college_metric <- value_cols[1]
  stopifnot(college_metric %in% value_cols)
  
  # Pre-period subset + optional weekly aggregation
  d <- df %>%
    mutate(.date = as.Date(.data[[date_col]])) %>%
    filter(.date >= as.Date(pre_start), .date <= as.Date(pre_end)) %>%
    { if (weekly)
      to_weekly(
        .,
        week_start = week_start,
        date_col   = date_col,
        group_cols = dma_col,
        value_cols = value_cols
      )
      else .
    } %>%
    mutate(
      dma_name = .data[[dma_col]],
      date  = as.Date(date),
      month = month(date)
    )
  
  # Eligibility / periods
  base_summary <- d %>%
    group_by(dma_name) %>%
    summarise(
      n_periods = n(),
      eligible  = n_periods >= min_pre_periods,
      .groups = "drop"
    )
  
  # Per-metric summary tables, then join by dma_name (avoids dma_name...1 issues)
  metric_tables <- purrr::map(value_cols, function(metric) {
    d %>%
      group_by(dma_name) %>%
      summarise(
        !!paste0(metric, "_mean") := mean(.data[[metric]], na.rm = TRUE),
        !!paste0(metric, "_sd")   := sd(.data[[metric]], na.rm = TRUE),
        !!paste0(metric, "_cv")   := sd(.data[[metric]], na.rm = TRUE) /
          pmax(mean(.data[[metric]], na.rm = TRUE), 1e-9),
        !!paste0(metric, "_trend") := {
          tt <- as.numeric(date - min(date))
          y  <- .data[[metric]]
          if (sum(!is.na(y)) < 3 || length(unique(tt)) < 2) NA_real_
          else coef(lm(y ~ tt))[2]
        },
        .groups = "drop"
      )
  })
  metric_summaries <- purrr::reduce(metric_tables, left_join, by = "dma_name")
  
  # -----------------------------
  # College proxy (Change #2): regression-based fall lift
  # y ~ fall + tt ; normalized by mean(y)
  # -----------------------------
  college_df <- d %>%
    group_by(dma_name) %>%
    summarise(
      college_idx = {
        dd <- cur_data()
        y  <- dd[[college_metric]]
        ok <- !is.na(y)
        
        if (sum(ok) < min_non_na_obs) NA_real_ else {
          fall <- dd$month %in% fall_months
          tt   <- as.numeric(dd$date - min(dd$date))
          
          if (length(unique(tt[ok])) < 2) NA_real_ else {
            fit <- lm(y ~ fall + tt)
            sm  <- summary(fit)
            
            # standardized fall lift (t-stat)
            t_fall <- unname(sm$coefficients["fallTRUE", "t value"])
            t_fall
          }
        }
      },
      .groups = "drop"
    )
  
  # Optionally add metric-specific college index name (e.g., active_users_college_idx)
  if (isTRUE(add_metric_specific_college_col)) {
    college_df <- college_df %>%
      mutate(!!paste0(college_metric, "_college_idx") := college_idx)
  }
  
  base_summary %>%
    left_join(metric_summaries, by = "dma_name") %>%
    left_join(college_df, by = "dma_name")
}

# -----------------------------
# Identify exclusions
# - NYC/LA by name pattern
# - Top-K by size metric
# - College-heavy by quantile / threshold
# -----------------------------
identify_exclusions <- function(features,
                                size_metric,
                                top_k_exclude = 10,
                                exclude_name_patterns = c("NEW YORK", "LOS ANGELES"),
                                college_col = "college_idx",
                                college_exclude_quantile = 0.90,
                                college_idx_min = NULL) {
  
  stopifnot(size_metric %in% names(features))
  stopifnot(college_col %in% names(features))
  
  f <- features %>% filter(eligible)
  
  by_name <- f %>%
    filter(str_detect(str_to_upper(dma_name),
                      str_c(exclude_name_patterns, collapse = "|"))) %>%
    transmute(dma_name, reason = "excluded_by_name")
  
  by_topk <- f %>%
    arrange(desc(.data[[size_metric]])) %>%
    slice_head(n = top_k_exclude) %>%
    transmute(dma_name, reason = paste0("top_", top_k_exclude, "_by_", size_metric))
  
  q_cut <- quantile(f[[college_col]], probs = college_exclude_quantile, na.rm = TRUE)
  by_college <- f %>%
    filter(.data[[college_col]] >= q_cut) %>%
    transmute(dma_name, reason = paste0("college_top_", round((1 - college_exclude_quantile) * 100), "pct"))
  
  if (!is.null(college_idx_min)) {
    by_college2 <- f %>%
      filter(.data[[college_col]] >= college_idx_min) %>%
      transmute(dma_name, reason = paste0("college_ge_", college_idx_min))
    by_college <- bind_rows(by_college, by_college2)
  }
  
  bind_rows(by_name, by_topk, by_college) %>%
    distinct(dma_name, reason)
}

# -----------------------------
# Rank remaining DMAs as candidates
# - Generic: pick a base metric, and we derive columns:
#     <metric>_mean, <metric>_cv, <metric>_trend
# -----------------------------
rank_candidates <- function(features,
                            exclusions,
                            base_metric,
                            college_col = "college_idx",
                            min_mean = 0) {
  
  mean_col  <- paste0(base_metric, "_mean")
  cv_col    <- paste0(base_metric, "_cv")
  trend_col <- paste0(base_metric, "_trend")
  
  stopifnot(all(c(mean_col, cv_col, trend_col, college_col) %in% names(features)))
  
  excl_names <- unique(exclusions$dma_name)
  
  features %>%
    filter(eligible) %>%
    filter(!(dma_name %in% excl_names)) %>%
    filter(.data[[mean_col]] >= min_mean) %>%
    mutate(
      score =
        log1p(.data[[mean_col]]) -
        2.0 * .data[[cv_col]] -
        0.0001 * abs(.data[[trend_col]]) -
        0.5 * pmax(.data[[college_col]], 0)  # penalize positive fall lift
    ) %>%
    arrange(desc(score))
}

# -----------------------------
# Suggest treated set: diversify across size tiers
# -----------------------------
suggest_treated_set <- function(ranked_candidates,
                                treated_n = 8,
                                size_metric,
                                tiers = 5,
                                seed = 1) {
  set.seed(seed)
  stopifnot(size_metric %in% names(ranked_candidates))
  
  d <- ranked_candidates %>%
    mutate(size_tier = ntile(.data[[size_metric]], tiers))
  
  per_tier <- max(1, floor(treated_n / tiers))
  
  picked <- d %>%
    group_by(size_tier) %>%
    slice_head(n = per_tier) %>%
    ungroup() %>%
    slice_head(n = treated_n)
  
  if (nrow(picked) < treated_n) {
    need <- treated_n - nrow(picked)
    picked <- bind_rows(
      picked,
      d %>% filter(!(dma_name %in% picked$dma_name)) %>% slice_head(n = need)
    )
  }
  
  picked %>% slice_head(n = treated_n)
}

# ============================================================
# Example usage (your case)
# ============================================================
features <- build_dma_features(
  df,
  pre_start  = "2024-01-01",
  pre_end    = "2025-08-01",
  dma_col    = "city",
  date_col   = "date",
  value_cols = c("active_users", "signups"),
  college_metric = "active_users",   # or "signups" if you insist
  fall_months = c(8, 9, 10)          # tweak if needed
)

exclusions <- identify_exclusions(
  features,
  size_metric = "signups_mean",      # or "active_users_mean"
  top_k_exclude = 10,
  college_col = "college_idx",
  college_exclude_quantile = 0.90
)

ranked <- rank_candidates(
  features, exclusions,
  base_metric = "signups",           # uses signups_mean/cv/trend
  college_col = "college_idx",
  min_mean = 0
)

treated <- suggest_treated_set(
  ranked,
  treated_n = 10,
  size_metric = "signups_mean",
  tiers = 5
)

################################################################################
# Threshold sweeps for engagement intensity metrics vs retention outcomes
# Outputs: consolidated table with recommended threshold k, knee point, and
# stability flags (by cohort week) + definition quality indicators.
#
# HOW TO USE (high level):
# 1) Prepare df with:
#    - user_id
#    - cohort_date (Date or YYYY-MM-DD string)
#    - retention outcomes: retained_7_14d, retained_28d (0/1)
#    - baseline covariates: sessions_7d, days_active_7d, platform, tenure_days, etc.
#    - engagement intensity metrics to sweep: counts like qchat_messages_7d, etc.
# 2) Set engagement_metrics, outcomes, covariates
# 3) Run consolidated <- run_all_sweeps(...)
################################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(splines)
  library(rlang)
})

################################################################################
# 0) Helper functions
################################################################################

# Knee finder:
# Identify diminishing returns point along a lift curve as the first k where the
# marginal gain falls below (frac * max marginal gain). This gives a "knee" that
# is often a good, interpretable threshold.
find_knee <- function(lift_by_k, ks, frac = 0.2) {
  # lift_by_k must be aligned with ks and ordered by ks ascending
  d <- diff(lift_by_k)
  if (length(d) == 0 || all(!is.finite(d))) return(NA_integer_)
  max_d <- max(d, na.rm = TRUE)
  if (!is.finite(max_d) || max_d <= 0) return(NA_integer_)
  idx <- which(d <= frac * max_d)[1]
  if (is.na(idx)) return(NA_integer_)
  ks[idx + 1]  # +1 because diff shifts indices
}

# Stability label:
# Uses (1) variability of lift across cohort weeks and (2) sign flips
# to flag brittle definitions.
stability_flag <- function(lift_sd, sign_flip_rate,
                           sd_thresh = 0.01, flip_thresh = 0.25) {
  if (!is.finite(lift_sd)) return("unknown")
  if (is.finite(sign_flip_rate) && sign_flip_rate >= flip_thresh) return("unstable_sign")
  if (lift_sd >= sd_thresh) return("unstable_magnitude")
  "stable"
}

# Definition quality label:
# A "good" definition should:
#  - have meaningful lift (not just noise)
#  - not be vanishingly small coverage
#  - be stable across cohort weeks
definition_quality <- function(lift, coverage, stability,
                               min_lift = 0.003,      # 0.3pp absolute lift default
                               min_coverage = 0.02) { # at least 2% of users qualify
  if (!is.finite(lift) || !is.finite(coverage)) return("unknown")
  if (coverage < min_coverage) return("low_coverage")
  if (lift < min_lift) return("low_lift")
  if (stability != "stable") return("brittle")
  "good"
}

# Choose a sensible k grid for each metric automatically:
# Uses the 95th percentile (capped) to avoid long-tail thresholds with tiny sample.
default_k_grid <- function(x, k_max_cap = 30) {
  x <- x[is.finite(x)]
  if (length(x) < 100) return(1:10)
  k_max <- min(k_max_cap, as.integer(stats::quantile(x, 0.95, na.rm = TRUE)))
  k_max <- max(k_max, 5)
  1:k_max
}

################################################################################
# 1) Core threshold sweep (one metric × one outcome)
################################################################################
# For each threshold k:
#  - Define engaged_k = 1(metric >= k)
#  - Fit logistic regression for outcome with covariates
#  - Estimate adjusted lift using marginal standardization:
#       lift(k) = E[p(outcome|engaged_k=1)] - E[p(outcome|engaged_k=0)]
#  - Compute coverage(k) = share of users engaged at k
#  - Compute score(k) = lift(k) * coverage(k)  (a simple lift-coverage tradeoff)

run_threshold_sweep <- function(df,
                                metric_col,
                                outcome_col,
                                covariates,
                                ks = NULL) {
  
  metric <- sym(metric_col)
  outcome <- sym(outcome_col)
  
  # Basic cleaning and canonical columns
  d0 <- df %>%
    filter(!is.na(!!metric), is.finite(!!metric), !is.na(!!outcome)) %>%
    mutate(
      .metric_val = as.numeric(!!metric),
      .outcome = as.integer(!!outcome)
    )
  
  # If too small, return empty
  if (nrow(d0) < 200) {
    return(tibble(
      metric = metric_col, outcome = outcome_col,
      k = integer(), lift = numeric(), coverage = numeric(), score = numeric(),
      n = integer()
    ))
  }
  
  # Default k grid is based on metric distribution
  if (is.null(ks)) ks <- default_k_grid(d0$.metric_val)
  
  # Build model formula once
  cov_str <- paste(covariates, collapse = " + ")
  f_str <- paste0(".outcome ~ engaged_k",
                  ifelse(nchar(cov_str) > 0, paste0(" + ", cov_str), ""))
  f <- as.formula(f_str)
  
  # Iterate over thresholds
  map_dfr(ks, function(k) {
    d <- d0 %>% mutate(engaged_k = as.integer(.metric_val >= k))
    covg <- mean(d$engaged_k == 1, na.rm = TRUE)
    
    # Skip degenerate thresholds where almost everyone or almost no one qualifies
    if (!is.finite(covg) || covg < 0.01 || covg > 0.99) {
      return(tibble(k = k, lift = NA_real_, coverage = covg, score = NA_real_, n = nrow(d)))
    }
    
    # Fit logistic regression, protect against errors (perfect separation, etc.)
    fit <- tryCatch(glm(f, data = d, family = binomial()),
                    error = function(e) NULL)
    if (is.null(fit)) {
      return(tibble(k = k, lift = NA_real_, coverage = covg, score = NA_real_, n = nrow(d)))
    }
    
    # Marginal standardization:
    # "What would predicted retention be if everyone were engaged vs not engaged?"
    p1 <- mean(predict(fit, newdata = mutate(d, engaged_k = 1), type = "response"), na.rm = TRUE)
    p0 <- mean(predict(fit, newdata = mutate(d, engaged_k = 0), type = "response"), na.rm = TRUE)
    lift <- p1 - p0
    
    tibble(k = k, lift = lift, coverage = covg, score = lift * covg, n = nrow(d))
  }) %>%
    mutate(metric = metric_col, outcome = outcome_col)
}

################################################################################
# 2) Stability sweeps by cohort week
################################################################################
# We compute the same sweep overall + per cohort_week:
#   cohort_week = week-start date (Monday) for cohort_date
#
# Then we can estimate how "brittle" each threshold is:
#  - lift_sd(k): sd of lift across cohort weeks
#  - sign_flip_rate(k): fraction of cohort weeks with lift(k) < 0

run_sweep_with_stability <- function(df,
                                     metric_col,
                                     outcome_col,
                                     covariates,
                                     ks = NULL,
                                     cohort_col = "cohort_date") {
  
  # Robust week start (Monday). Works without lubridate.
  # Convert to Date first.
  d <- df %>%
    mutate(.cohort_date = as.Date(.data[[cohort_col]])) %>%
    filter(!is.na(.cohort_date)) %>%
    mutate(
      # Monday-based week start:
      cohort_week = .cohort_date - ((as.integer(format(.cohort_date, "%u")) - 1))
    )
  
  overall <- run_threshold_sweep(d, metric_col, outcome_col, covariates, ks)
  
  by_week <- d %>%
    group_by(cohort_week) %>%
    group_modify(~ run_threshold_sweep(.x, metric_col, outcome_col, covariates, ks)) %>%
    ungroup()
  
  list(overall = overall, by_week = by_week)
}

################################################################################
# 3) Summarize one sweep into a single "definition row"
################################################################################
# Produces:
#  - recommended_k: best threshold under an objective (score, penalized by instability)
#  - knee_k: diminishing returns point on lift curve
#  - stability: stable / unstable_sign / unstable_magnitude / unknown
#  - definition_quality: good / low_lift / low_coverage / brittle / unknown

summarize_sweep <- function(sweep_overall,
                            sweep_by_week,
                            knee_frac = 0.2,
                            stability_sd_thresh = 0.01,
                            stability_flip_thresh = 0.25,
                            stability_penalty_lambda = 0.0,
                            quality_min_lift = 0.003,
                            quality_min_coverage = 0.02) {
  
  # Keep only valid rows
  o <- sweep_overall %>%
    filter(is.finite(lift), is.finite(score), is.finite(coverage)) %>%
    arrange(k)
  
  metric_name <- unique(sweep_overall$metric)
  outcome_name <- unique(sweep_overall$outcome)
  
  if (nrow(o) == 0) {
    return(tibble(
      metric = metric_name,
      outcome = outcome_name,
      recommended_k = NA_integer_,
      knee_k = NA_integer_,
      best_score = NA_real_,
      lift_at_reco = NA_real_,
      coverage_at_reco = NA_real_,
      lift_sd_at_reco = NA_real_,
      sign_flip_rate_at_reco = NA_real_,
      stability = "unknown",
      definition_quality = "unknown"
    ))
  }
  
  # Knee point on the lift curve (diminishing returns)
  knee_k <- find_knee(o$lift, o$k, frac = knee_frac)
  
  # Stability per k across cohort weeks
  stab <- sweep_by_week %>%
    filter(is.finite(lift)) %>%
    group_by(k) %>%
    summarise(
      lift_sd = sd(lift, na.rm = TRUE),
      sign_flip_rate = mean(lift < 0, na.rm = TRUE),
      .groups = "drop"
    )
  
  o2 <- o %>% left_join(stab, by = "k")
  
  # Penalized objective (optional): score - lambda * lift_sd
  o2 <- o2 %>% mutate(score_adj = score - stability_penalty_lambda * coalesce(lift_sd, 0))
  
  reco_row <- o2 %>% slice_max(order_by = score_adj, n = 1, with_ties = FALSE)
  
  stab_label <- stability_flag(
    lift_sd = reco_row$lift_sd[[1]],
    sign_flip_rate = reco_row$sign_flip_rate[[1]],
    sd_thresh = stability_sd_thresh,
    flip_thresh = stability_flip_thresh
  )
  
  qual_label <- definition_quality(
    lift = reco_row$lift[[1]],
    coverage = reco_row$coverage[[1]],
    stability = stab_label,
    min_lift = quality_min_lift,
    min_coverage = quality_min_coverage
  )
  
  tibble(
    metric = metric_name,
    outcome = outcome_name,
    recommended_k = reco_row$k[[1]],
    knee_k = knee_k,
    best_score = reco_row$score[[1]],
    lift_at_reco = reco_row$lift[[1]],
    coverage_at_reco = reco_row$coverage[[1]],
    lift_sd_at_reco = reco_row$lift_sd[[1]],
    sign_flip_rate_at_reco = reco_row$sign_flip_rate[[1]],
    stability = stab_label,
    definition_quality = qual_label
  )
}

################################################################################
# 4) Run sweeps across all metrics × outcomes and consolidate
################################################################################

run_all_sweeps <- function(df,
                           engagement_metrics,
                           outcomes,
                           covariates,
                           ks_map = NULL,                 # optional: named list metric->vector ks
                           cohort_col = "cohort_date",
                           stability_penalty_lambda = 0.0,
                           knee_frac = 0.2,
                           stability_sd_thresh = 0.01,
                           stability_flip_thresh = 0.25,
                           quality_min_lift = 0.003,
                           quality_min_coverage = 0.02) {
  
  crossing(metric = engagement_metrics, outcome = outcomes) %>%
    mutate(res = pmap(list(metric, outcome), function(metric, outcome) {
      
      # Choose a k grid for this metric, if provided
      ks <- if (!is.null(ks_map) && metric %in% names(ks_map)) ks_map[[metric]] else NULL
      
      # Compute overall + cohort-week sweeps
      sweeps <- run_sweep_with_stability(
        df = df,
        metric_col = metric,
        outcome_col = outcome,
        covariates = covariates,
        ks = ks,
        cohort_col = cohort_col
      )
      
      # Summarize to one row with reco k, knee, stability, quality
      summarize_sweep(
        sweep_overall = sweeps$overall,
        sweep_by_week = sweeps$by_week,
        knee_frac = knee_frac,
        stability_sd_thresh = stability_sd_thresh,
        stability_flip_thresh = stability_flip_thresh,
        stability_penalty_lambda = stability_penalty_lambda,
        quality_min_lift = quality_min_lift,
        quality_min_coverage = quality_min_coverage
      )
    })) %>%
    select(res) %>%
    unnest(res) %>%
    arrange(outcome, desc(best_score))
}

################################################################################
# 5) Example configuration (EDIT THIS TO MATCH YOUR DATA)
################################################################################

# 5a) Pick engagement intensity metrics (counts / volumes), NOT pre-thresholded flags.
# Option 1 (recommended): explicit list
# engagement_metrics <- c("qchat_messages_7d", "learn_items_7d", "test_questions_7d")

# Option 2: regex-driven auto-selection (adjust to your naming conventions)
# engagement_metrics <- names(df) %>%
#   keep(~ grepl("(messages|questions|items|rounds|_7d$)", .x)) %>%
#   setdiff(c("sessions_7d", "days_active_7d"))  # don't sweep covariates by accident

# 5b) Outcomes: you said focus on 7–14d and 28d retention.
# outcomes <- c("retained_7_14d", "retained_28d")

# 5c) Covariates: keep these lightweight and widely available.
# Use ns(tenure_days, df=4) only if tenure_days exists and is numeric.
# covariates <- c(
#   "sessions_7d",
#   "days_active_7d",
#   "platform",
#   "ns(tenure_days, df = 4)",
#   "entry_surface",
#   "geo"
# )

# 5d) Optional: per-metric k grids if some metrics have very different scales
# ks_map <- list(
#   qchat_messages_7d = 1:20,
#   learn_items_7d = c(1:20, 25, 30, 40, 50)
# )

################################################################################
# 6) Run + output consolidated table
################################################################################

# consolidated <- run_all_sweeps(
#   df = df,
#   engagement_metrics = engagement_metrics,
#   outcomes = outcomes,
#   covariates = covariates,
#   cohort_col = "cohort_date",
#   stability_penalty_lambda = 0.0, # raise (e.g., 0.5-2.0) to prefer stable thresholds
#   knee_frac = 0.2,
#   stability_sd_thresh = 0.01,     # tune based on typical lift scale and sample sizes
#   stability_flip_thresh = 0.25,
#   quality_min_lift = 0.003,       # tune: 0.3pp absolute lift
#   quality_min_coverage = 0.02     # tune: >=2% of users qualify
# )
#
# print(consolidated, n = 200)

################################################################################
# 7) Optional: quick sanity filters (recommended)
################################################################################

# consolidated_filtered <- consolidated %>%
#   filter(definition_quality %in% c("good", "low_lift", "brittle")) %>%
#   arrange(outcome, desc(best_score))
#
# print(consolidated_filtered, n = 200)
################################################################################

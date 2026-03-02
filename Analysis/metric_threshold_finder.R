library(dplyr)
library(purrr)
library(tibble)

find_engagement_thresholds_all <- function(df,
                                           retained_col,
                                           metric_cols = NULL,
                                           exclude_cols = c("user_id"),
                                           exclude_regex = "(^user_|_id$|^id$)",
                                           n_thresholds = 200,
                                           min_above_n = 200,
                                           dr_fraction = 0.25,
                                           early_quantile = 0.15,
                                           consec = 3,
                                           min_unique = 20,
                                           max_total = NULL) {
  
  if (!retained_col %in% names(df)) stop("retained_col not found in df")
  
  # Pick candidate metrics
  if (is.null(metric_cols)) {
    metric_cols <- names(df)[sapply(df, is.numeric)]
    metric_cols <- setdiff(metric_cols, retained_col)
  }
  
  metric_cols <- setdiff(metric_cols, exclude_cols)
  metric_cols <- metric_cols[!grepl(exclude_regex, metric_cols, ignore.case = TRUE)]
  
  retained <- as.integer(as.logical(df[[retained_col]]))
  
  map_dfr(metric_cols, function(metric_name) {
    
    metric <- df[[metric_name]]
    
    if (!is.numeric(metric)) {
      return(tibble(
        metric = metric_name,
        optimal_threshold = NA_real_,
        status = "skipped",
        reason = "non_numeric",
        n_non_na = sum(!is.na(metric)),
        n_unique = NA_integer_,
        early_gain = NA_real_,
        cutoff = NA_real_,
        max_threshold_scanned = NA_real_
      ))
    }
    
    dat <- tibble(metric = as.numeric(metric), retained = retained) %>%
      filter(!is.na(metric), !is.na(retained))
    
    n_non_na <- nrow(dat)
    n_unique <- dplyr::n_distinct(dat$metric)
    
    if (n_non_na < 2 * min_above_n) {
      return(tibble(
        metric = metric_name,
        optimal_threshold = NA_real_,
        status = "no_result",
        reason = "too_few_rows_after_na_drop",
        n_non_na = n_non_na,
        n_unique = n_unique,
        early_gain = NA_real_,
        cutoff = NA_real_,
        max_threshold_scanned = NA_real_
      ))
    }
    
    if (n_unique < min_unique) {
      return(tibble(
        metric = metric_name,
        optimal_threshold = NA_real_,
        status = "no_result",
        reason = "too_few_unique_values",
        n_non_na = n_non_na,
        n_unique = n_unique,
        early_gain = NA_real_,
        cutoff = NA_real_,
        max_threshold_scanned = NA_real_
      ))
    }
    
    if (!is.null(max_total)) {
      dat <- dat %>% filter(metric <= max_total)
      n_non_na <- nrow(dat)
      n_unique <- dplyr::n_distinct(dat$metric)
      if (n_non_na < 2 * min_above_n) {
        return(tibble(
          metric = metric_name,
          optimal_threshold = NA_real_,
          status = "no_result",
          reason = "too_few_rows_after_max_total_cap",
          n_non_na = n_non_na,
          n_unique = n_unique,
          early_gain = NA_real_,
          cutoff = NA_real_,
          max_threshold_scanned = NA_real_
        ))
      }
    }
    
    # Candidate thresholds on a quantile grid (robust to long tails)
    thresholds <- unique(dat$metric |> quantile(probs = seq(0, 0.99, length.out = n_thresholds), na.rm = TRUE))
    thresholds <- sort(unique(thresholds))
    thresholds <- thresholds[thresholds > min(dat$metric, na.rm = TRUE)]
    
    # Build conditional retention curve r(T) = P(retained | metric >= T)
    curve <- map_dfr(thresholds, function(T) {
      above <- dat$metric >= T
      n_above <- sum(above)
      if (n_above < min_above_n) return(NULL)
      tibble(
        threshold = T,
        retained_rate = mean(dat$retained[above], na.rm = TRUE),
        n_above = n_above
      )
    })
    
    if (nrow(curve) < 10) {
      return(tibble(
        metric = metric_name,
        optimal_threshold = NA_real_,
        status = "no_result",
        reason = "too_few_threshold_points_after_min_above_n",
        n_non_na = n_non_na,
        n_unique = n_unique,
        early_gain = NA_real_,
        cutoff = NA_real_,
        max_threshold_scanned = ifelse(nrow(curve) == 0, NA_real_, max(curve$threshold, na.rm = TRUE))
      ))
    }
    
    curve <- curve %>%
      arrange(threshold) %>%
      mutate(d_ret = c(NA_real_, diff(retained_rate))) %>%
      filter(is.finite(d_ret))
    
    if (nrow(curve) < 10) {
      return(tibble(
        metric = metric_name,
        optimal_threshold = NA_real_,
        status = "no_result",
        reason = "too_few_points_after_derivative",
        n_non_na = n_non_na,
        n_unique = n_unique,
        early_gain = NA_real_,
        cutoff = NA_real_,
        max_threshold_scanned = max(curve$threshold, na.rm = TRUE)
      ))
    }
    
    early_end <- max(2, floor(nrow(curve) * early_quantile))
    early_gain <- median(curve$d_ret[1:early_end], na.rm = TRUE)
    
    if (!is.finite(early_gain) || early_gain <= 0) {
      return(tibble(
        metric = metric_name,
        optimal_threshold = NA_real_,
        status = "no_result",
        reason = "early_gain_non_positive",
        n_non_na = n_non_na,
        n_unique = n_unique,
        early_gain = early_gain,
        cutoff = NA_real_,
        max_threshold_scanned = max(curve$threshold, na.rm = TRUE)
      ))
    }
    
    cutoff <- dr_fraction * early_gain
    ok <- curve$d_ret <= cutoff
    
    if (length(ok) < consec) {
      return(tibble(
        metric = metric_name,
        optimal_threshold = NA_real_,
        status = "no_result",
        reason = "not_enough_points_for_consec",
        n_non_na = n_non_na,
        n_unique = n_unique,
        early_gain = early_gain,
        cutoff = cutoff,
        max_threshold_scanned = max(curve$threshold, na.rm = TRUE)
      ))
    }
    
    run_ok <- sapply(seq_len(length(ok) - consec + 1),
                     function(i) all(ok[i:(i + consec - 1)]))
    idx <- which(run_ok)[1]
    
    if (is.na(idx)) {
      return(tibble(
        metric = metric_name,
        optimal_threshold = NA_real_,
        status = "no_result",
        reason = "no_diminishing_point_found",
        n_non_na = n_non_na,
        n_unique = n_unique,
        early_gain = early_gain,
        cutoff = cutoff,
        max_threshold_scanned = max(curve$threshold, na.rm = TRUE)
      ))
    }
    
    tibble(
      metric = metric_name,
      optimal_threshold = curve$threshold[idx],
      status = "ok",
      reason = NA_character_,
      n_non_na = n_non_na,
      n_unique = n_unique,
      early_gain = early_gain,
      cutoff = cutoff,
      max_threshold_scanned = max(curve$threshold, na.rm = TRUE)
    )
  })
}

th_7d  <- find_engagement_thresholds_all(el_data, retained_col = "retained7d")


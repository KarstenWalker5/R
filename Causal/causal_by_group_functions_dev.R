# These functions allow one to run a CausalImpact analysis by group, plot a nice summary table of results,
# and save a ggPlot formatted PDF of the model plots.
# Additionally there is a function that finds the top 10 best control time series for each test series and then
# creates a synthetic control metric total. This allows us to run CausalImpact at the group level.

# WARNING, DO NOT CLICK ON THE run_causal_impact RESULT DATAFRAME IN R STUDIO. THE CAUSAL IMPACT MODELS ARE VERY LARGE
# AND WILL BE MOVED FROM THE TEMPORARY sink() MEMORY INTO RAM AND FREEZE YOUR SESSION. PLEASE ADD THE FOLLOWING LINE OF CODE
# WHEN YOU CALL THE FUNCTION

# Run CausalImpact by group
# Requires a nested data frame that has an intervention_date, response variable, optional control variable, and end date.
run_causal_impact <- function(data,
                              intervention_date_col,
                              response_var,
                              control_vars,
                              max_post_period_date = NULL) {
  library(dplyr)
  library(zoo)
  library(CausalImpact)
  library(lubridate)
  
  data <- data %>% arrange(date)
  
  # Convert column name to character and extract values
  intervention_date_col <- as.character(intervention_date_col)
  intervention_date <- unique(data[[intervention_date_col]])
  if (!inherits(intervention_date, "Date")) {
    intervention_date <- ymd(intervention_date)
  }
  
  # Determine post period end
  if (!is.null(max_post_period_date)) {
    post_end <- if (is.character(max_post_period_date) && max_post_period_date %in% names(data)) {
      val <- unique(data[[max_post_period_date]])
      if (!inherits(val, "Date")) lubridate::ymd(val) else val
    } else {
      max_post_period_date
    }
  } else {
    post_end <- max(data$date)
  }
  
  # Set pre and post periods
  pre_period <- c(min(data$date), intervention_date - 1)
  post_period <- c(intervention_date, post_end)
  
  # Safely build model vars
  model_vars <- if (is.null(control_vars) || length(control_vars) == 0) {
    response_var
  } else {
    c(response_var, control_vars)
  }
  
  # Build zoo object
  ts_data <- tryCatch({
    zoo::zoo(data %>% select(all_of(model_vars)), order.by = data$date)
  }, error = function(e) return(NULL))
  
  if (is.null(ts_data)) {
    return(list(list(
      model = NULL,
      summary_text = "Zoo creation failed",
      printed_summary = "Zoo creation failed",
      p_value = NA_real_,
      rel_effect = NA_real_,
      rel_effect_lower = NA_real_,
      rel_effect_upper = NA_real_
    )))
  }
  
  temp_file <- tempfile()
  sink(temp_file)
  
  result <- tryCatch({
    impact <- CausalImpact(ts_data, pre_period, post_period)
    s <- summary(impact)$summary
    
    sink()
    
    printed_summary <- paste(capture.output(print(impact)), collapse = "\n")
    summary_text <- paste(capture.output(summary(impact, "report")), collapse = "\n")
    
    list(
      model = impact,
      summary_text = summary_text,
      printed_summary = printed_summary,
      p_value = s["p", "Average"],
      rel_effect = s["Relative.effect", "Average"],
      rel_effect_lower = s["Relative.effect", "Lower"],
      rel_effect_upper = s["Relative.effect", "Upper"]
    )
  }, error = function(e) {
    sink()
    list(
      model = NULL,
      summary_text = paste("CausalImpact failed:", e$message),
      printed_summary = paste("CausalImpact failed:", e$message),
      p_value = NA_real_,
      rel_effect = NA_real_,
      rel_effect_lower = NA_real_,
      rel_effect_upper = NA_real_
    )
  }, finally = {
    if (sink.number() > 0) sink()
    unlink(temp_file)
  })
  
  return(result) # Ensure output is length-1 list for mutate(pmap(...))
}

# summary function
summarize_causal_results <- function(df, group_col = "group") {
  df %>%
    mutate(
      summary_text = map_chr(causal_impact, function(x) {
        if (!is.null(x$summary_text)) x$summary_text else "No summary"
      }),
      printed_summary = map_chr(causal_impact, function(x) {
        if (!is.null(x$printed_summary)) x$printed_summary else "No printed summary"
      }),
      p_value = map_dbl(causal_impact, function(x) {
        if (!is.null(x$model$summary$p)) x$model$summary$p[1] else NA_real_
      }),
      rel_effect = map_dbl(causal_impact, function(x) {
        if (!is.null(x$model$summary$RelEffect)) x$model$summary$RelEffect[1] else NA_real_
      }),
      rel_effect_lower = map_dbl(causal_impact, function(x) {
        if (!is.null(x$model$summary$RelEffect.lower)) x$model$summary$RelEffect.lower[1] else NA_real_
      }),
      rel_effect_upper = map_dbl(causal_impact, function(x) {
        if (!is.null(x$model$summary$RelEffect.upper)) x$model$summary$RelEffect.upper[1] else NA_real_
      }),
      abs_effect = map_dbl(causal_impact, function(x) {
        if (!is.null(x$model$summary$AbsEffect)) x$model$summary$AbsEffect[1] else NA_real_
      }),
      abs_effect_lower = map_dbl(causal_impact, function(x) {
        if (!is.null(x$model$summary$AbsEffect.lower)) x$model$summary$AbsEffect.lower[1] else NA_real_
      }),
      abs_effect_upper = map_dbl(causal_impact, function(x) {
        if (!is.null(x$model$summary$AbsEffect.upper)) x$model$summary$AbsEffect.upper[1] else NA_real_
      }),
      rel_effect_percent = rel_effect * 100,
      abs_effect_percent = abs_effect * 100,
      rel_effect_ci = sprintf("%.1f%% [%.1f%%, %.1f%%]",
                              rel_effect * 100,
                              rel_effect_lower * 100,
                              rel_effect_upper * 100),
      abs_effect_ci = sprintf("%.1f%% [%.1f%%, %.1f%%]",
                              abs_effect * 100,
                              abs_effect_lower * 100,
                              abs_effect_upper * 100),
      uplift_vs_baseline = map_dbl(causal_impact, function(x) {
        predicted <- tryCatch({
          if (!is.null(x$model$summary$Pred) && length(x$model$summary$Pred) > 0) {
            x$model$summary$Pred[1]
          } else {
            NA_real_
          }
        }, error = function(e) NA_real_)
        
        effect <- tryCatch({
          if (!is.null(x$model$summary$AbsEffect) && length(x$model$summary$AbsEffect) > 0) {
            x$model$summary$AbsEffect[1]
          } else {
            NA_real_
          }
        }, error = function(e) NA_real_)
        
        if (!is.na(predicted) && predicted != 0 && !is.na(effect)) {
          effect / predicted
        } else {
          NA_real_
        }
      }),
      uplift_vs_baseline_percent = uplift_vs_baseline * 100,
      significant = !is.na(p_value) & p_value < 0.05
    )%>%
    select(-causal_impact)
}

# Plotting function
save_causal_impact_ggplots <- function(df, output_file = "causal_impact_plots.pdf") {
  library(ggplot2)
  library(patchwork)
  
  pdf(output_file, width = 12, height = 10)
  
  for (i in seq_len(nrow(df))) {
    group_label <- df$group[i]
    ci_obj <- df$causal_impact[[i]]
    
    if (!is.null(ci_obj) && !is.null(ci_obj$model$series)) {
      df_plot <- as.data.frame(ci_obj$model$series) %>%
        tibble::rownames_to_column("date") %>%
        dplyr::mutate(date = as.Date(date))
      
      p1 <- ggplot(df_plot, aes(x = date)) +
        geom_line(aes(y = response)) +
        geom_line(aes(y = point.pred), linetype = "dashed", color = "blue") +
        geom_ribbon(aes(ymin = point.pred.lower, ymax = point.pred.upper), alpha = 0.2) +
        labs(title = "Observed vs. Predicted") +
        theme_causal()
      
      p2 <- ggplot(df_plot, aes(x = date)) +
        geom_line(aes(y = point.effect), color = "darkgreen") +
        geom_ribbon(aes(ymin = point.effect.lower, ymax = point.effect.upper), alpha = 0.2) +
        labs(title = "Pointwise Causal Effect") +
        theme_causal()
      
      p3 <- ggplot(df_plot, aes(x = date)) +
        geom_line(aes(y = cum.effect), color = "darkred") +
        geom_ribbon(aes(ymin = cum.effect.lower, ymax = cum.effect.upper), alpha = 0.2) +
        labs(title = "Cumulative Effect") +
        theme_causal()
      
      full_plot <- (p1 / p2 / p3) + patchwork::plot_annotation(title = group_label)
      
      print(full_plot)  # Force patchwork rendering
    }
  }
  
  dev.off()
}

# This function replaces manually looping through each data frame after matching. It finds the best control group
# for each test group, maps each group to the identified control group, and then loops through each group and calculates control metric.
build_synthetic_control <- function(data,
                                    matching_metric = "sales",
                                    date_col        = "date",
                                    city_col        = "city",
                                    treat_col       = "treat",
                                    match_window    = c("2024-06-01", "2024-09-15"),
                                    n_matches       = 3) {
  suppressWarnings({
    library(MarketMatching)
    library(dplyr)
    library(lubridate)
    library(purrr)
    library(tibble)
    library(glue)
    
    # 0) Ensure types: date as Date, city as character ----------------------
    data <- data %>%
      mutate(
        !!date_col := ymd(.data[[date_col]]),
        !!city_col := as.character(.data[[city_col]])
      )
    
    # 1) Identify treated cities (we want synthetic controls for these) -----
    treated_cities <- data %>%
      filter(.data[[treat_col]] == 1) %>%
      pull(.data[[city_col]]) %>%
      unique()
    
    # candidates for donors (controls)
    control_cities <- data %>%
      filter(.data[[treat_col]] == 0) %>%
      pull(.data[[city_col]]) %>%
      unique()
    
    if (length(treated_cities) == 0 || length(control_cities) == 0) {
      stop("Need at least one treated and one control city.")
    }
    
    # 2) Data for best_matches: all cities, but only the needed columns -----
    matching_data <- data %>%
      select(
        !!city_col,
        !!date_col,
        !!matching_metric
      )
    
    # 3) Run best_matches for all treated cities at once --------------------
    mm <- best_matches(
      data                  = matching_data,
      markets_to_be_matched = treated_cities,
      id_variable           = city_col,
      date_variable         = date_col,
      matching_variable     = matching_metric,
      parallel              = FALSE,
      warping_limit         = 2,
      dtw_emphasis          = 0.5,
      start_match_period    = as.Date(match_window[1]),
      end_match_period      = as.Date(match_window[2]),
      matches               = 10        # ask for plenty, we'll trim to n_matches
    )
    
    best_tab <- as_tibble(mm$BestMatches)
    
    # The column with treated IDs has the same name as city_col (e.g. "city")
    id_col_name <- mm$MarketID
    
    names(best_tab)[names(best_tab) == id_col_name] <- "treated_city"
    
    # 4) Build treated → control mapping, enforcing donors are controls -----
    match_mapping <- best_tab %>%
      mutate(
        treated_city = as.character(treated_city),
        BestControl  = as.character(BestControl)
      ) %>%
      filter(
        treated_city %in% treated_cities,
        BestControl  %in% control_cities
      ) %>%
      group_by(treated_city) %>%
      # choose top n_matches by whatever ranking column is there
      # Often there is a 'RelativeDistance' or 'Correlation' column.
      arrange(RelativeDistance, .by_group = TRUE) %>%
      slice_head(n = n_matches) %>%
      ungroup() %>%
      rename(control_city = BestControl)
    
    if (nrow(match_mapping) == 0) {
      message("No valid treated→control matches after filtering to control cities.")
      return(tibble())
    }
    
    # 5) Build the panel: one row per date × treated_city with treat + synth totals
    synth_panel <- map_dfr(unique(match_mapping$treated_city), function(cty) {
      matched_controls <- match_mapping %>%
        filter(treated_city == cty) %>%
        pull(control_city)
      
      # treated series
      df_treat <- data %>%
        filter(.data[[city_col]] == cty) %>%
        group_by(.data[[date_col]]) %>%
        summarise(
          treat_sales = sum(.data[[matching_metric]], na.rm = TRUE),
          post        = max(.data[["post"]]),   # assuming 0/1 per date
          .groups     = "drop"
        ) %>%
        rename(date = !!date_col)
      
      # synthetic control = sum across matched control cities (treat == 0)
      df_ctrl <- data %>%
        filter(.data[[city_col]] %in% matched_controls,
               .data[[treat_col]] == 0) %>%
        group_by(.data[[date_col]]) %>%
        summarise(
          synth_sales = median(.data[[matching_metric]], na.rm = TRUE),
          .groups     = "drop"
        ) %>%
        rename(date = !!date_col)
      
      df_treat %>%
        left_join(df_ctrl, by = "date") %>%
        mutate(treated_city = cty)
    })
    
    synth_panel
  })
}

# Build Synthetic Control using Ridge Regression
build_synth_controls_ridge <- function(data,
                                 top_k_donors = 5,   # how many best-matching donors per treated unit
                                 alpha = 0           # 0 = ridge, 1 = lasso, in between = elastic net
) {
  # Packages
  stopifnot(all(c("date", "city", "treat", "post", "sales") %in% names(data)))
  suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(glmnet)
  })
  
  # Ensure date is ordered
  data <- data %>%
    mutate(date = as.Date(date)) %>%
    arrange(date, city)
  
  # Get treated vs donor cities
  treat_status <- data %>%
    distinct(city, treat)
  
  treated_cities <- treat_status %>%
    filter(treat == 1) %>%
    pull(city)
  
  donor_cities <- treat_status %>%
    filter(treat == 0) %>%
    pull(city)
  
  if (length(treated_cities) == 0) stop("No treated cities (treat == 1) found.")
  if (length(donor_cities)  == 0) stop("No donor cities (treat == 0) found.")
  
  # Get post flag by date (assumed common across cities)
  post_by_date <- data %>%
    distinct(date, post)
  
  # Wide format: one row per date, one column per city, plus post flag
  wide <- data %>%
    select(date, city, sales) %>%
    tidyr::pivot_wider(names_from = city, values_from = sales) %>%
    left_join(post_by_date, by = "date") %>%
    arrange(date)
  
  date_vec <- wide$date
  post_vec <- wide$post
  
  pre_idx  <- which(post_vec == 0)
  post_idx <- which(post_vec == 1)
  
  if (length(pre_idx) < 5) {
    stop("Very few pre-period observations; synthetic control will be unstable.")
  }
  
  # Helper to build synthetic control for a single treated city
  build_for_one_city <- function(city_name) {
    # Outcome for treated city
    if (!city_name %in% names(wide)) {
      warning("City ", city_name, " not found as a column in wide data; skipping.")
      return(NULL)
    }
    
    y_treated <- wide[[city_name]]
    
    # Matrix of all donors
    donor_cols_present <- intersect(donor_cities, names(wide))
    if (length(donor_cols_present) == 0) {
      warning("No donor columns present in wide data for treated city ", city_name, "; skipping.")
      return(NULL)
    }
    X_all_donors <- as.matrix(wide[donor_cols_present])
    
    # Pre-period slices
    y_pre <- y_treated[pre_idx]
    X_pre <- X_all_donors[pre_idx, , drop = FALSE]
    
    # Handle any missing data in pre-period
    cc <- complete.cases(cbind(y_pre, X_pre))
    y_pre_cc <- y_pre[cc]
    X_pre_cc <- X_pre[cc, , drop = FALSE]
    
    if (length(y_pre_cc) < 5) {
      warning("Too few complete pre-period observations for city ", city_name, "; skipping.")
      return(NULL)
    }
    
    # 1. Compute correlations with each donor in pre-period
    cor_vec <- sapply(seq_len(ncol(X_pre_cc)), function(j) {
      tryCatch(
        cor(y_pre_cc, X_pre_cc[, j], use = "complete.obs"),
        error = function(e) NA_real_
      )
    })
    names(cor_vec) <- colnames(X_pre_cc)
    
    # Order donors by absolute correlation (best matches first)
    ord <- order(abs(cor_vec), decreasing = TRUE, na.last = TRUE)
    donor_ranked <- names(cor_vec)[ord]
    
    # Choose top_k_donors (or all if fewer)
    k <- min(top_k_donors, length(donor_ranked))
    donors_used <- donor_ranked[seq_len(k)]
    
    # ---- SAFER SUBSETTING: use positions, drop missing ----
    donor_idx <- match(donors_used, colnames(X_pre_cc))
    donor_idx <- donor_idx[!is.na(donor_idx)]
    
    if (length(donor_idx) == 0) {
      warning("No matching donor columns found for city ", city_name, "; skipping.")
      return(NULL)
    }
    
    X_pre_k <- X_pre_cc[, donor_idx, drop = FALSE]
    X_all_k <- X_all_donors[, donor_idx, drop = FALSE]
    donors_used <- colnames(X_pre_k)  # aligned names
    
    # 2. Fit regularized regression in pre-period
    cv_fit <- cv.glmnet(
      x = X_pre_k,
      y = y_pre_cc,
      alpha = alpha,
      intercept = FALSE,
      standardize = TRUE
    )
    
    best_lambda <- cv_fit$lambda.min
    
    fit <- glmnet(
      x = X_pre_k,
      y = y_pre_cc,
      alpha = alpha,
      lambda = best_lambda,
      intercept = FALSE,
      standardize = TRUE
    )
    
    coef_vec <- as.numeric(coef(fit))   # includes intercept (0) + donors
    w <- coef_vec[-1]                   # drop intercept
    names(w) <- donors_used
    
    # 3. Build synthetic control series for all dates
    control_series <- as.numeric(X_all_k %*% w)
    
    tibble(
      date          = date_vec,
      city          = city_name,
      sales         = y_treated,
      control_sales = control_series
    )
  }
  
  # Apply to each treated city and bind together
  synth_list <- lapply(treated_cities, build_for_one_city)
  synth_df <- bind_rows(Filter(Negate(is.null), synth_list)) %>%
    arrange(city, date)
  
  synth_df
}

# Using tidysynth
build_synth_controls_tidysynth <- function(data, ...) {

  # ---- Packages & checks ----
  if (!requireNamespace("tidysynth", quietly = TRUE)) {
    stop("Package 'tidysynth' is required. Please install it with install.packages('tidysynth').")
  }
  
  stopifnot(all(c("date", "city", "treat", "post", "sales") %in% names(data)))
  
  library(dplyr)
  library(purrr)
  library(tidysynth)
  
  # Normalize / order
  data <- data %>%
    mutate(date = as.Date(date)) %>%
    arrange(date, city)
  
  treated_cities <- data %>%
    filter(treat == 1) %>%
    distinct(city) %>%
    pull(city)
  
  if (length(treated_cities) == 0) {
    stop("No treated cities (treat == 1) found.")
  }
  
  # ---- Inner function: run tidysynth for ONE treated city ----
  run_for_city <- function(city_name) {
    # This treated city + all donors
    d_city <- data %>%
      filter(city == city_name | treat == 0)
    
    # Treatment date: first post-period date for this treated city
    treat_date <- d_city %>%
      filter(city == city_name, post == 1) %>%
      summarise(treat_date = min(date)) %>%
      pull(treat_date)
    
    if (is.na(treat_date)) {
      warning("No post-period for treated city ", city_name, "; skipping.")
      return(NULL)
    }
    
    # Pre-period time window
    pre_window <- d_city %>%
      filter(date <= treat_date) %>%
      pull(date) %>%
      unique() %>%
      sort()
    
    if (length(pre_window) < 5) {
      warning("Very few pre-period observations for city ", city_name, "; skipping.")
      return(NULL)
    }
    
    # tidysynth pipeline
    synth_obj <- d_city %>%
      synthetic_control(
        outcome = sales,
        unit    = city,
        time    = date,
        i_unit  = city_name,
        i_time  = treat_date,
        generate_placebos = FALSE
      ) %>%
      # simple example predictor; add more if you like
      generate_predictor(
        time_window = pre_window,
        mean_sales  = mean(sales, na.rm = TRUE)
      ) %>%
      generate_weights(
        optimization_window = pre_window
      ) %>%
      generate_control()
    
    sc_df <- grab_synthetic_control(synth_obj)
    
    sc_df %>%
      transmute(
        date          = as.Date(time_unit),
        city          = city_name,
        sales         = real_y,
        control_sales = synth_y
      )
  }
  
  # ---- Map over all treated cities, bind results ----
  synth_list <- map(treated_cities, run_for_city)
  
  bind_rows(compact(synth_list)) %>%
    arrange(city, date)
}

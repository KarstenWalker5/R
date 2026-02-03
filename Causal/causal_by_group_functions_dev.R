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
  
  theme_causal<- function() {
    ggthemes::theme_fivethirtyeight() %+replace% 
      theme(
        plot.title = element_text(hjust=0.5, size=22,face="bold", color='#ebba34'),
        plot.subtitle = element_text(hjust=0.5, size=12, face="italic"),
        legend.position = "bottom",
        legend.text=element_text(size=16, face="bold"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,size=12,face="bold"),
        axis.text.y = element_text(size=10,face="bold"),
        axis.title = element_text(size=11,face="bold"), 
        axis.title.x = element_blank(),
        plot.background=element_rect(fill="white"),
        panel.background=element_rect(fill="white"),
        legend.background=element_rect(fill="white")
      )
  }
  
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

# Plots to the RStudio viewer using ggplot
plot_causal_impact_ggplots <- function(df, pause = TRUE) {
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tibble)
  
  theme_causal<- function() {
    ggthemes::theme_fivethirtyeight() %+replace% 
      theme(
        plot.title = element_text(hjust=0.5, size=22,face="bold", color='#ebba34'),
        plot.subtitle = element_text(hjust=0.5, size=12, face="italic"),
        legend.position = "bottom",
        legend.text=element_text(size=16, face="bold"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,size=12,face="bold"),
        axis.text.y = element_text(size=10,face="bold"),
        axis.title = element_text(size=11,face="bold"), 
        axis.title.x = element_blank(),
        plot.background=element_rect(fill="white"),
        panel.background=element_rect(fill="white"),
        legend.background=element_rect(fill="white")
      )
  }
  
  for (i in seq_len(nrow(df))) {
    group_label <- df$group[i]
    ci_obj <- df$causal_impact[[i]]
    
    if (!is.null(ci_obj) && !is.null(ci_obj$model$series)) {
      
      df_plot <- as.data.frame(ci_obj$model$series) %>%
        rownames_to_column("date") %>%
        mutate(date = as.Date(date))
      
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
      
      full_plot <- (p1 / p2 / p3) +
        patchwork::plot_annotation(title = group_label)
      
      print(full_plot)
      
      if (pause && i < nrow(df)) {
        message("Press <Enter> to view next plot...")
        invisible(readline())
      }
    }
  }
}


# This function replaces manually looping through each data frame after matching. It finds the best control group
# for each test group, maps each group to the identified control group, and then loops through each group and calculates control metric.
build_synthetic_control <- function(data,
                                    matching_metric = "sales",
                                    date_col        = "date",
                                    city_col        = "city",
                                    treat_col       = "treat",
                                    n_matches       = 3,
                                    match_window    = NULL,
                                    intervention_date_col = NULL,
                                    pre_period_days = NULL,
                                    per_city_window = FALSE,
                                    warping_limit   = 2,
                                    dtw_emphasis    = 0.5,
                                    parallel        = FALSE) {
  suppressWarnings({
    library(MarketMatching)
    library(dplyr)
    library(lubridate)
    library(purrr)
    library(tibble)
    
    treat_name <- paste0("treatment_", matching_metric)
    synth_name <- paste0("synthetic_", matching_metric)
    
    eps <- 1e-9
    
    data <- data %>%
      mutate(
        !!date_col := ymd(.data[[date_col]]),
        !!city_col := as.character(.data[[city_col]])
      )
    
    treated_cities <- data %>%
      filter(.data[[treat_col]] == 1) %>%
      pull(.data[[city_col]]) %>%
      unique()
    
    control_cities <- data %>%
      filter(.data[[treat_col]] == 0) %>%
      pull(.data[[city_col]]) %>%
      unique()
    
    if (length(treated_cities) == 0 || length(control_cities) == 0) {
      stop("Need at least one treated and one control city.")
    }
    
    matching_data <- data %>%
      select(!!city_col, !!date_col, !!matching_metric)
    
    compute_start_end <- function(df, end_date) {
      end_date <- as.Date(end_date)
      start_date <- min(df[[date_col]], na.rm = TRUE)
      if (!is.null(pre_period_days)) {
        start_date <- max(start_date, end_date - as.integer(pre_period_days) + 1)
      }
      list(start = start_date, end = end_date)
    }
    
    use_match_window <- !is.null(match_window) && length(match_window) == 2
    use_intervention <- !is.null(intervention_date_col)
    
    if (!use_match_window && !use_intervention) {
      stop("Provide either match_window (length 2) OR intervention_date_col.")
    }
    
    # helper to compute match mapping from a best_matches result table
    build_mapping <- function(best_tab) {
      best_tab %>%
        mutate(
          treated_city = as.character(treated_city),
          BestControl  = as.character(BestControl)
        ) %>%
        filter(
          treated_city %in% treated_cities,
          BestControl  %in% control_cities
        ) %>%
        group_by(treated_city) %>%
        arrange(RelativeDistance, .by_group = TRUE) %>%
        slice_head(n = n_matches) %>%
        ungroup() %>%
        rename(control_city = BestControl) %>%
        group_by(treated_city) %>%
        mutate(
          raw_w  = 1 / (RelativeDistance + eps),
          weight = raw_w / sum(raw_w)
        ) %>%
        ungroup() %>%
        select(treated_city, control_city, RelativeDistance, weight)
    }
    
    # -------------------------
    # Build match_mapping
    # -------------------------
    if (!per_city_window) {
      
      if (use_match_window) {
        start_match_period <- as.Date(match_window[1])
        end_match_period   <- as.Date(match_window[2])
      } else {
        # shared end date uses earliest intervention among treated cities
        end_dates <- data %>%
          filter(.data[[city_col]] %in% treated_cities) %>%
          pull(.data[[intervention_date_col]]) %>%
          unique()
        
        end_dates <- ymd(end_dates)
        end_match_period <- min(end_dates, na.rm = TRUE) - 1
        
        se <- compute_start_end(data, end_match_period)
        start_match_period <- se$start
        end_match_period   <- se$end
      }
      
      mm <- best_matches(
        data                  = matching_data,
        markets_to_be_matched = treated_cities,
        id_variable           = city_col,
        date_variable         = date_col,
        matching_variable     = matching_metric,
        parallel              = parallel,
        warping_limit         = warping_limit,
        dtw_emphasis          = dtw_emphasis,
        start_match_period    = start_match_period,
        end_match_period      = end_match_period,
        matches               = 10
      )
      
      best_tab <- as_tibble(mm$BestMatches)
      id_col_name <- mm$MarketID
      names(best_tab)[names(best_tab) == id_col_name] <- "treated_city"
      
      match_mapping <- build_mapping(best_tab)
      
    } else {
      
      if (!use_intervention) {
        stop("per_city_window = TRUE requires intervention_date_col.")
      }
      
      # lookup intervention date per treated city
      city_iv <- data %>%
        filter(.data[[city_col]] %in% treated_cities) %>%
        select(!!city_col, !!intervention_date_col) %>%
        distinct() %>%
        mutate(intervention_date = ymd(.data[[intervention_date_col]])) %>%
        select(treated_city = !!city_col, intervention_date)
      
      match_mapping <- map_dfr(treated_cities, function(cty) {
        end_date <- city_iv %>%
          filter(treated_city == cty) %>%
          pull(intervention_date)
        
        if (length(end_date) == 0 || is.na(end_date[1])) return(tibble())
        
        end_match_period <- as.Date(end_date[1]) - 1
        se <- compute_start_end(data, end_match_period)
        
        mm <- best_matches(
          data                  = matching_data,
          markets_to_be_matched = cty,
          id_variable           = city_col,
          date_variable         = date_col,
          matching_variable     = matching_metric,
          parallel              = parallel,
          warping_limit         = warping_limit,
          dtw_emphasis          = dtw_emphasis,
          start_match_period    = se$start,
          end_match_period      = se$end,
          matches               = 10
        )
        
        best_tab <- as_tibble(mm$BestMatches)
        id_col_name <- mm$MarketID
        names(best_tab)[names(best_tab) == id_col_name] <- "treated_city"
        
        build_mapping(best_tab) %>%
          filter(treated_city == cty)
      })
    }
    
    if (nrow(match_mapping) == 0) {
      message("No valid treated→control matches after filtering to control cities.")
      return(tibble())
    }
    
    # -------------------------
    # Build synthetic panel
    # -------------------------
    synth_panel <- map_dfr(unique(match_mapping$treated_city), function(cty) {
      
      w_tbl <- match_mapping %>%
        filter(treated_city == cty) %>%
        arrange(desc(weight)) %>%
        transmute(control_city, weight)
      
      donor_cities  <- w_tbl$control_city
      donor_weights <- w_tbl$weight
      donor_weights_named <- setNames(donor_weights, donor_cities)
      
      df_treat <- data %>%
        filter(.data[[city_col]] == cty) %>%
        group_by(.data[[date_col]]) %>%
        summarise(
          !!treat_name := sum(.data[[matching_metric]], na.rm = TRUE),
          post         = max(.data[["post"]]),
          .groups      = "drop"
        ) %>%
        rename(date = !!date_col)
      
      df_ctrl <- data %>%
        filter(.data[[treat_col]] == 0) %>%
        mutate(control_city = as.character(.data[[city_col]])) %>%
        inner_join(w_tbl, by = "control_city") %>%
        group_by(.data[[date_col]]) %>%
        summarise(
          !!synth_name := weighted.mean(.data[[matching_metric]], w = weight, na.rm = TRUE),
          .groups      = "drop"
        ) %>%
        rename(date = !!date_col)
      
      df_treat %>%
        left_join(df_ctrl, by = "date") %>%
        mutate(
          treated_city = cty,
          synthetic_cities        = list(donor_cities),
          synthetic_weights       = list(donor_weights),
          synthetic_weights_named = list(donor_weights_named)
        )
    })
    
    synth_panel
  })
}


# Return Top N matches per city, no weighting, no synthetic control
best_matches_all_cities <- function(data,
                                    matching_metric      = "sales",
                                    date_col             = "date",
                                    city_col             = "city",
                                    treat_col            = "treat",
                                    matches_per_city     = 10,
                                    donors_controls_only = TRUE,
                                    compute_weights      = FALSE,
                                    warping_limit        = 2,
                                    dtw_emphasis         = 0.5,
                                    parallel             = FALSE,
                                    intervention_date_col = NULL,
                                    pre_period_days      = NULL,
                                    per_city_window      = FALSE,
                                    match_window         = NULL) {
  suppressWarnings({
    library(MarketMatching)
    library(dplyr)
    library(lubridate)
    library(tibble)
    library(purrr)
    
    eps <- 1e-9
    
    data <- data %>%
      mutate(
        !!date_col := ymd(.data[[date_col]]),
        !!city_col := as.character(.data[[city_col]])
      )
    
    all_cities <- data %>% pull(.data[[city_col]]) %>% unique()
    if (length(all_cities) < 2) stop("Need at least 2 cities to compute matches.")
    
    control_cities <- data %>%
      filter(.data[[treat_col]] == 0) %>%
      pull(.data[[city_col]]) %>%
      unique()
    
    matching_data <- data %>% select(!!city_col, !!date_col, !!matching_metric)
    
    # Helper: compute match start/end given a df and a chosen end date
    compute_start_end <- function(df, end_date) {
      end_date <- as.Date(end_date)
      start_date <- min(df[[date_col]], na.rm = TRUE)
      if (!is.null(pre_period_days)) {
        start_date <- max(start_date, end_date - as.integer(pre_period_days) + 1)
      }
      list(start = start_date, end = end_date)
    }
    
    # Decide matching windows
    use_match_window <- !is.null(match_window) && length(match_window) == 2
    use_intervention <- !is.null(intervention_date_col)
    
    if (!use_match_window && !use_intervention) {
      stop("Provide either match_window (length 2) OR intervention_date_col.")
    }
    
    # Approach A: one shared window for all cities
    if (!per_city_window) {
      
      if (use_match_window) {
        start_match_period <- as.Date(match_window[1])
        end_match_period <- as.Date(match_window[2])
      } else {
        # Shared end date = earliest intervention date among treated cities
        treated_cities <- data %>%
          filter(.data[[treat_col]] == 1) %>%
          pull(.data[[city_col]]) %>%
          unique()
        
        if (length(treated_cities) == 0) stop("No treated cities found for intervention-based matching.")
        
        end_date <- data %>%
          filter(.data[[city_col]] %in% treated_cities) %>%
          pull(.data[[intervention_date_col]]) %>%
          unique()
        
        end_date <- ymd(end_date)
        end_date <- min(end_date, na.rm = TRUE) - 1  # end of pre period
        
        se <- compute_start_end(data, end_date)
        start_match_period <- se$start
        end_match_period <- se$end
      }
      
      mm <- best_matches(
        data                  = matching_data,
        markets_to_be_matched = all_cities,
        id_variable           = city_col,
        date_variable         = date_col,
        matching_variable     = matching_metric,
        parallel              = parallel,
        warping_limit         = warping_limit,
        dtw_emphasis          = dtw_emphasis,
        start_match_period    = start_match_period,
        end_match_period      = end_match_period,
        matches               = matches_per_city
      )
      
      best_tab <- as_tibble(mm$BestMatches)
      id_col_name <- mm$MarketID
      names(best_tab)[names(best_tab) == id_col_name] <- "city"
      
      out <- best_tab %>%
        mutate(city = as.character(city), BestControl = as.character(BestControl)) %>%
        filter(city != BestControl)
      
      if (donors_controls_only) out <- out %>% filter(BestControl %in% control_cities)
      
      out <- out %>%
        group_by(city) %>%
        arrange(RelativeDistance, .by_group = TRUE) %>%
        slice_head(n = matches_per_city) %>%
        mutate(rank = row_number()) %>%
        ungroup() %>%
        rename(match_city = BestControl)
      
      if (compute_weights) {
        out <- out %>%
          group_by(city) %>%
          mutate(raw_w = 1 / (RelativeDistance + eps), weight = raw_w / sum(raw_w)) %>%
          ungroup()
      }
      
      return(out %>% select(city, match_city, RelativeDistance, rank, any_of("weight")))
    }
    
    # Approach B: per-city window (most correct)
    # We run best_matches one city at a time so each city uses its own pre period end date
    if (use_match_window) {
      stop("per_city_window=TRUE is intended for intervention_date_col based matching. Use per_city_window=FALSE for match_window.")
    }
    
    # Build a lookup of intervention dates per city
    city_interventions <- data %>%
      select(!!city_col, !!intervention_date_col) %>%
      distinct() %>%
      mutate(
        intervention_date = ymd(.data[[intervention_date_col]])
      ) %>%
      select(city = !!city_col, intervention_date)
    
    out <- map_dfr(all_cities, function(cty) {
      end_date <- city_interventions %>%
        filter(city == cty) %>%
        pull(intervention_date)
      
      if (length(end_date) == 0 || is.na(end_date[1])) {
        return(tibble())  # skip if no intervention date for this city
      }
      
      end_date <- as.Date(end_date[1]) - 1
      se <- compute_start_end(data, end_date)
      
      mm <- best_matches(
        data                  = matching_data,
        markets_to_be_matched = cty,
        id_variable           = city_col,
        date_variable         = date_col,
        matching_variable     = matching_metric,
        parallel              = parallel,
        warping_limit         = warping_limit,
        dtw_emphasis          = dtw_emphasis,
        start_match_period    = se$start,
        end_match_period      = se$end,
        matches               = matches_per_city
      )
      
      best_tab <- as_tibble(mm$BestMatches)
      id_col_name <- mm$MarketID
      names(best_tab)[names(best_tab) == id_col_name] <- "city"
      
      best_tab %>%
        mutate(city = as.character(city), BestControl = as.character(BestControl)) %>%
        filter(city == cty, city != BestControl) %>%
        { if (donors_controls_only) dplyr::filter(., BestControl %in% control_cities) else . } %>%
        arrange(RelativeDistance) %>%
        slice_head(n = matches_per_city) %>%
        mutate(rank = row_number()) %>%
        rename(match_city = BestControl) %>%
        { if (compute_weights)
          dplyr::mutate(., raw_w = 1 / (RelativeDistance + eps)) %>%
            dplyr::mutate(weight = raw_w / sum(raw_w))
          else .
        } %>%
        select(city, match_city, RelativeDistance, rank, any_of("weight"))
    })
    
    out
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

# Uses DiD to run incrementality tests by city
incrementality_by_city <- function(df_city, treat_col, synth_col) {
  df_city <- df_city %>%
    arrange(date) %>%
    mutate(
      gap = .data[[treat_col]] - .data[[synth_col]]
    ) %>%
    filter(!is.na(gap), !is.na(post))
  
  pre  <- df_city %>% filter(post == 0)
  post <- df_city %>% filter(post == 1)
  
  # Guardrails
  if (nrow(pre) < 2 || nrow(post) < 2) {
    return(tibble(
      n_pre = nrow(pre), n_post = nrow(post),
      pre_mean_gap = NA_real_, post_mean_gap = NA_real_,
      did = NA_real_, did_p = NA_real_,
      post_mean = NA_real_, post_p = NA_real_
    ))
  }
  
  # 1) DID on time series gaps: (mean post gap) - (mean pre gap)
  did <- mean(post$gap) - mean(pre$gap)
  
  # Simple OLS on gap ~ post (robust SE would be nicer, but keeping deps minimal)
  fit <- tryCatch(stats::lm(gap ~ post, data = df_city), error = function(e) NULL)
  did_p <- if (!is.null(fit)) summary(fit)$coefficients["post", "Pr(>|t|)"] else NA_real_
  
  # 2) Paired t-test on post gaps vs 0 (are post gaps positive?)
  tt <- tryCatch(stats::t.test(post$gap, mu = 0), error = function(e) NULL)
  
  tibble(
    n_pre = nrow(pre),
    n_post = nrow(post),
    pre_mean_gap = mean(pre$gap),
    post_mean_gap = mean(post$gap),
    did = did,
    did_p = did_p,
    post_mean = if (!is.null(tt)) unname(tt$estimate) else NA_real_,
    post_p = if (!is.null(tt)) tt$p.value else NA_real_,
    post_ci_low = if (!is.null(tt)) tt$conf.int[1] else NA_real_,
    post_ci_high = if (!is.null(tt)) tt$conf.int[2] else NA_real_
  )
}

# Placebo Test Single City
placebo_test_one_city <- function(df_city,
                                  treat_col,
                                  synth_col,
                                  post_col = "post",
                                  n_placebos = 500,
                                  seed = 1,
                                  one_sided = TRUE,
                                  keep_dist = FALSE) {
  set.seed(seed)
  
  df_city <- df_city %>%
    dplyr::arrange(date) %>%
    dplyr::mutate(gap = .data[[treat_col]] - .data[[synth_col]]) %>%
    dplyr::filter(!is.na(date), !is.na(gap), !is.na(.data[[post_col]]))
  
  post_idx <- which(df_city[[post_col]] == 1)
  n_post <- length(post_idx)
  
  if (n_post < 5 || nrow(df_city) < (n_post + 10)) {
    return(tibble::tibble(
      n_post = n_post,
      observed_cum_lift = NA_real_,
      observed_avg_lift = NA_real_,
      placebo_p = NA_real_
    ))
  }
  
  observed <- sum(df_city$gap[post_idx], na.rm = TRUE)
  observed_avg <- mean(df_city$gap[post_idx], na.rm = TRUE)
  
  n <- nrow(df_city)
  max_start <- n - n_post + 1
  starts <- sample.int(max_start, size = n_placebos, replace = TRUE)
  
  placebo_lifts <- vapply(
    starts,
    function(s) sum(df_city$gap[s:(s + n_post - 1)], na.rm = TRUE),
    numeric(1)
  )
  
  placebo_p <- if (one_sided) {
    mean(placebo_lifts >= observed)
  } else {
    mean(abs(placebo_lifts) >= abs(observed))
  }
  
  tibble::tibble(
    n_post = n_post,
    observed_cum_lift = observed,
    observed_avg_lift = observed_avg,
    
    placebo_p = placebo_p,
    z_score = (observed - mean(placebo_lifts)) / sd(placebo_lifts),
    snr = observed / sd(placebo_lifts),
    
    ci_90_low = quantile(placebo_lifts, 0.05),
    ci_90_high = quantile(placebo_lifts, 0.95),
    ci_95_low = quantile(placebo_lifts, 0.025),
    ci_95_high = quantile(placebo_lifts, 0.975),
    
    rank_percentile = mean(placebo_lifts < observed),
    sign_stability = mean(sign(placebo_lifts) == sign(observed)),
    
    placebo_lifts = if (keep_dist) list(placebo_lifts) else NULL
  )
}

# Placebo test, multiple cities. Requires single city version to be run first.

run_placebo_tests_by_city <- function(df,
                                      group_col = "treated_city",
                                      treat_col,
                                      synth_col,
                                      post_col = "post",
                                      n_placebos = 500,
                                      seed = 1,
                                      one_sided = TRUE,
                                      keep_dist = FALSE) {
  
  df %>%
    dplyr::group_by(.data[[group_col]]) %>%
    dplyr::group_modify(~{
      placebo_test_one_city(
        df_city = .x,
        treat_col = treat_col,
        synth_col = synth_col,
        post_col = post_col,
        n_placebos = n_placebos,
        seed = seed,
        one_sided = one_sided,
        keep_dist = keep_dist
      )
    }) %>%
    dplyr::ungroup() %>%
    dplyr::rename(!!group_col := .data[[group_col]])
}

# MDE/Power test single city
# Estimates how big an incremental effect you could reliably detect for a single treated city, given historical noise and post-period length
mde_one_city <- function(df_city,
                         treat_col,
                         synth_col,
                         post_col = "post",
                         alpha = 0.05,
                         target_power = 0.80,
                         two_sided = TRUE) {
  
  df_city <- df_city %>%
    dplyr::arrange(date) %>%
    dplyr::mutate(gap = .data[[treat_col]] - .data[[synth_col]]) %>%
    dplyr::filter(!is.na(date), !is.na(gap), !is.na(.data[[post_col]]))
  
  pre  <- df_city %>% dplyr::filter(.data[[post_col]] == 0)
  post <- df_city %>% dplyr::filter(.data[[post_col]] == 1)
  
  n_pre  <- nrow(pre)
  n_post <- nrow(post)
  pre_sd <- stats::sd(pre$gap)
  
  if (is.na(pre_sd) || n_post <= 1) {
    return(tibble::tibble(
      n_pre = n_pre,
      n_post = n_post,
      pre_sd_gap = pre_sd,
      mde_avg = NA_real_,
      mde_cum = NA_real_,
      alpha = alpha,
      target_power = target_power,
      two_sided = two_sided
    ))
  }
  
  z_alpha <- if (two_sided) stats::qnorm(1 - alpha / 2) else stats::qnorm(1 - alpha)
  z_beta  <- stats::qnorm(target_power)
  
  mde_avg <- (z_alpha + z_beta) * pre_sd / sqrt(n_post)
  mde_cum <- mde_avg * n_post
  
  tibble::tibble(
    n_pre = n_pre,
    n_post = n_post,
    pre_sd_gap = pre_sd,
    mde_avg = mde_avg,
    mde_cum = mde_cum,
    alpha = alpha,
    target_power = target_power,
    two_sided = two_sided
  )
}

# MDE by group
run_mde_by_city <- function(df,
                            group_col = "treated_city",
                            treat_col,
                            synth_col,
                            post_col = "post",
                            alpha = 0.05,
                            target_power = 0.80,
                            two_sided = TRUE) {
  
  df %>%
    dplyr::group_by(.data[[group_col]]) %>%
    dplyr::group_modify(~{
      mde_one_city(
        df_city = .x,
        treat_col = treat_col,
        synth_col = synth_col,
        post_col = post_col,
        alpha = alpha,
        target_power = target_power,
        two_sided = two_sided
      )
    }) %>%
    dplyr::ungroup() %>%
    dplyr::rename(!!group_col := .data[[group_col]])
}

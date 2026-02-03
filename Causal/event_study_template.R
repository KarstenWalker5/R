# Event-study on the GAP (treatment - synthetic)

library(dplyr)
library(fixest)
library(broom)
library(ggplot2)

# -Compute gap and per-unit event_time ----
# Replace these names if your dataset uses different columns:
#   unit_col = unit id, time_col = date or integer time, y_col = observed outcome,
#   y_synth_col = synthetic counterfactual, g_col = first-treated time (NA if never treated)
unit_col <- "unit"
time_col <- "time"
y_col <- "y"
y_synth_col <- "y_synth"   # synthetic counterfactual series (from your SC step)
g_col <- "g"               # first treated time per unit (NA for never-treated)

# Build dataframe with gap and event time
df_es <- panel %>%
  arrange(.data[[unit_col]], .data[[time_col]]) %>%
  mutate(
    # gap = observed - synthetic (incremental effect should be gap>0 post)
    gap = .data[[y_col]] - .data[[y_synth_col]],
    # If g is not present compute from D/post: first time where D==1 (if you have D)
    !!g_col := if (!g_col %in% names(.) && "D" %in% names(.)) {
      ave(ifelse(D==1, as.integer(.data[[time_col]]), NA_integer_), .data[[unit_col]],
          FUN = function(x) { v <- x[!is.na(x)]; if (length(v)) min(v) else NA_integer_ })
    } else .data[[g_col]]
  ) %>%
  group_by(.data[[unit_col]]) %>%
  mutate(
    # Compute event time for treated units; NA for never-treated
    event_time = ifelse(!is.na(.data[[g_col]]), as.integer(.data[[time_col]] - .data[[g_col]]), NA_integer_)
  ) %>%
  ungroup()

# Remove units without a valid g (optional), or keep them but mark treated_flag
df_es <- df_es %>%
  mutate(treated_flag = ifelse(!is.na(.data[[g_col]]), 1L, 0L))

# Bin tails (choose K leads/lags to display; cap extremes) ----
K <- 8     # adjust to how many leads / lags you want visible
df_es <- df_es %>%
  mutate(event_time_binned = ifelse(is.na(event_time), NA_integer_,
                                    pmin(pmax(event_time, -K), K)))

# Check counts by event_time 
counts_ev <- df_es %>% group_by(event_time_binned) %>% summarise(n = n())
print(counts_ev)

# Sun & Abraham on the GAP 
# Use unit + time FE, cluster at unit level
es_fit_gap <- feols(
  gap ~ sunab(event_time_binned, treated_flag, ref = -1) | 
    paste0(unit_col) + paste0(time_col),
  data = df_es,
  cluster = ~get(unit_col)
)

# Quick summary and aggregated ATT
print(summary(es_fit_gap, agg = "att"))

# Coefficients for plotting 
es_tidy <- broom::tidy(es_fit_gap, conf.int = TRUE) %>%
  # term names for sunab might look like "sunab::event_time_binned::k" or "event_time_binned::k"
  filter(grepl("event_time_binned", term) | grepl("sunab", term)) %>%
  mutate(
    # extract numeric event_time from term
    evt = as.integer(gsub(".*(?:::|\\()([^:()]+)\\).*", "\\1", term))
  ) %>%
  arrange(evt)

ggplot(es_tidy, aes(x = evt, y = estimate)) +
  geom_point() +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) +
  geom_vline(xintercept = -0.5, linetype = "dashed") +
  labs(x = "Event time", y = "Estimate (gap)", title = "Event study (Sun & Abraham on gap)") +
  theme_minimal()

# Pre-trend joint (Wald) test
# Identify pre-event terms (evt < 0)
pre_terms <- es_tidy %>% filter(evt < 0)
if (nrow(pre_terms) >= 2) {
  coef_names <- pre_terms$term
  coefs <- coef(es_fit_gap)[coef_names]
  V <- vcov(es_fit_gap, cluster = ~get(unit_col))[coef_names, coef_names, drop = FALSE]
  wald_stat <- as.numeric(t(coefs) %*% solve(V, coefs))
  df_w <- length(coefs)
  pval_pre <- 1 - pchisq(wald_stat, df_w)
  cat("Pre-trend joint Wald p-value:", signif(pval_pre,3), "\n")
} else {
  message("Not enough pre-coefs to run joint pre-trend test.")
}

# ---- If few clusters: wild cluster bootstrap (fixest::boottest) ----
# Example: test for joint effect at post window or aggregated ATT
# install.packages("fixest") # if needed
# boottest(es_fit_gap, cnames = "sunab::event_time_binned::0", cluster = df_es[[unit_col]])
# or run boottest for aggregate ATT term if available



# Difference-in-Differences template
# - Synthetic panel dataset generator
# - Parallel trends plots
# - Baseline DiD (TWFE)
# - Event study / dynamic DiD
# - Basic robustness / diagnostics

# Load Packages, set theme
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(fixest)  
  library(stringr)
})

set.seed(123)

theme_fancy <- function() {
  theme_minimal(base_family = "Asap Condensed") +
    theme(panel.grid.minor = element_blank(),
          legend.position = "bottom")
}

# Synthetic data
# Design:
# - Units: cohorts or schools or geo buckets (e.g., "school_id")
# - Daily observations
# - Metric: engagement proxy (e.g., sessions per active learner)
# - Treatment: product change rolled out on a specific date
# - Effect ramps up over ~2 weeks then stabilizes

n_units <- 150
n_days  <- 180   

units <- paste0("school_", str_pad(1:n_units, width = 3, pad = "0"))
dates <- seq.Date(
  from = as.Date("2024-01-01"),
  by   = "day",
  length.out = n_days
)

# Intervention date
treat_date <- as.Date("2024-04-01")

# 40% of units treated
treated_units <- sample(units, size = round(0.4 * n_units), replace = FALSE)

df <- tidyr::crossing(unit = units, date = dates) %>%
  mutate(
    treated = as.integer(unit %in% treated_units),
    post    = as.integer(date > treat_date),
    
    # Event time in days relative to intervention
    event_time = as.integer(date - treat_date)
  )


# 2. Fixed Effects
# Unit fixed effects 
unit_fe <- tibble(
  unit = units,
  alpha_i = rnorm(n_units, mean = 0, sd = 0.8))

# Daily time effects:
# - global trend
# - weekly seasonality (weekday/weekend)
# - mild noise
time_fe <- tibble(
  date = dates) %>%
  mutate(
    t = as.integer(date - min(date)),
    # Weekly pattern (higher usage mid-week, lower weekends)
    weekly = 0.6 * sin(2 * pi * t / 7),
    # Slow product growth / seasonality
    trend = 0.01 * t,
    gamma_t = trend + weekly + rnorm(n(), 0, 0.15)
  ) %>%
  select(date, gamma_t)

# Dynamic treatment effect (ramps over 14 days)
dynamic_effect <- function(k) {
  ifelse(
    k <= 0, 0,
    ifelse(k <= 14, 0.05 * k, 0.05 * 14)
  )
}

df <- df %>%
  left_join(unit_fe, by = "unit") %>%
  left_join(time_fe, by = "date") %>%
  mutate(
    tau = dynamic_effect(event_time),
    # Example covariate:
    # share of mobile users or new learners that day
    x_mobile = plogis(rnorm(n(), -0.3, 0.8)),
    # Outcome: daily engagement metric
    # (e.g., sessions per active learner)
    y = 2.5 +
      alpha_i +
      gamma_t +
      0.4 * x_mobile +
      (treated * post) * tau +
      rnorm(n(), 0, 0.4)
  )

# 3) Aggregate daily trends treated vs control
avg_daily <- df %>%
  group_by(date, treated) %>%
  summarise(
    y_mean = mean(y),
    y_se   = sd(y) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(group = ifelse(treated == 1, "Treated", "Control"))

ggplot(avg_daily, aes(x = date, y = y_mean, color = group)) +
  geom_line(linewidth = 1) +
  geom_ribbon(
    aes(ymin = y_mean - 1.96 * y_se, ymax = y_mean + 1.96 * y_se, fill = group),
    alpha = 0.15,
    color = NA ) +
  geom_vline(xintercept = treat_date, linetype = "dashed") +
  labs(
    title = "Daily engagement metric over time",
    subtitle = "Treated vs Control (Quizlet-style online metric)",
    x = "Date",
    y = "Mean daily engagement",
    color = NULL,
    fill = NULL) +
  theme_fancy()

# 4. Baseline DiD (TWFE)
m_did <- feols(
  y ~ treated * post | unit + date,
  data = df,
  cluster = ~ unit
)

summary(m_did)

# 5. Event study (daily dynamic effects)
# Reference period: event_time = -1 (day before launch)

m_es <- feols(
  y ~ i(event_time, treated, ref = -1) | unit + date,
  data = df,
  cluster = ~ unit
)

summary(m_es)

# Helper to extract coefficients
extract_fixest <- function(model) {
  ct <- fixest::coeftable(model)
  tibble(
    term = rownames(ct),
    estimate = ct[, "Estimate"],
    std_error = ct[, "Std. Error"]
  )
}

es_tbl <- extract_fixest(m_es) %>%
  filter(str_detect(term, "^event_time::")) %>%
  mutate(
    event_time = as.integer(str_extract(term, "-?\\d+")),
    ci_low = estimate - 1.96 * std_error,
    ci_high = estimate + 1.96 * std_error
  ) %>%
  arrange(event_time)

ggplot(es_tbl, aes(x = event_time, y = estimate)) +
  geom_hline(yintercept = 0, linewidth = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.3) +
  labs(
    title = "Daily event study (dynamic DiD)",
    subtitle = "Effect relative to day −1; dashed line = launch day",
    x = "Days since launch",
    y = "Estimated treatment effect"
  ) +
  theme_fancy()

# 6. DiD with covariate (precision gain)
m_did_x <- feols(
  y ~ treated * post + x_mobile | unit + date,
  data = df,
  cluster = ~ unit
)

summary(m_did_x)

# 7. Placebo test (fake earlier launch)

placebo_date <- as.Date("2024-03-01")

df_placebo <- df %>%
  mutate(
    post_p = as.integer(date > placebo_date),
    event_time_p = as.integer(date - placebo_date)
  )

m_placebo <- feols(
  y ~ treated * post_p | unit + date,
  data = df_placebo,
  cluster = ~ unit
)

summary(m_placebo)

# Placebo event study
m_es_placebo <- feols(
  y ~ i(event_time_p, treated, ref = -1) | unit + date,
  data = df_placebo,
  cluster = ~ unit
)

es_placebo_tbl <- extract_fixest(m_es_placebo) %>%
  filter(str_detect(term, "^event_time_p::")) %>%
  mutate(
    event_time = as.integer(str_extract(term, "-?\\d+")),
    ci_low = estimate - 1.96 * std_error,
    ci_high = estimate + 1.96 * std_error
  ) %>%
  arrange(event_time)

ggplot(es_placebo_tbl, aes(x = event_time, y = estimate)) +
  geom_hline(yintercept = 0, linewidth = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.3) +
  labs(
    title = "Placebo daily event study",
    subtitle = "Fake launch date — effects should be ~0 pre and post",
    x = "Days since placebo date",
    y = "Estimated treatment effect") +
  theme_fancy()

# 8. Pooled DiD regression (NO fixed effects)

m_did_pooled <- feols(
  y ~ treated * post,
  data = df,
  cluster = ~ unit
)

summary(m_did_pooled)

# The DiD estimate is the interaction term:
# treated:post  (this is the ATT under parallel trends)

#  DiD means (Control/Treated x Pre/Post)
did_means <- df %>%
  mutate(
    group = ifelse(treated == 1, "Treated", "Control"),
    period = ifelse(post == 1, "Post", "Pre")
  ) %>%
  group_by(group, period) %>%
  summarise(
    mean_y = mean(y),
    se_y = sd(y) / sqrt(n()),
    .groups = "drop"
  )

ggplot(did_means, aes(x = period, y = mean_y, group = group, color = group)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean_y - 1.96 * se_y, ymax = mean_y + 1.96 * se_y),
                width = 0.15) +
  labs(
    title = "Difference-in-Differences (no fixed effects)",
    x = NULL,
    y = "Mean outcome (y)",
    color = NULL
  ) +
  theme_fancy()

# 10 Add covariates, still no FE. Use only covariates not affected by treatment.

m_did_pooled_x <- feols(
  y ~ treated * post + x_mobile,
  data = df,
  cluster = ~ unit
)

summary(m_did_pooled_x)
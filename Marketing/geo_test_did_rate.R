# Geo-test analysis w/rate-based variance reduction
#
# Assumes:
# * Randomization was at DMA level (cluster randomized)
# * KPI is upgrades; exposure proxy is active_users
# * We use upgrade rate = upgrades / active_users for variance reduction
#
# Report:
# * Incremental upgrade rate per active user per day
# * % lift in rate
# * Convert rate effect -> incremental upgrades (optional but useful)
# * 95% CI + p-value
#
# Inputs expected in environment:
# assignment  # data frame with dma + group ("TEST","CONTROL","BAU")
# df2         # daily data with columns: dma, date, upgrades, active_users
# pre_start_date, pre_end_date, test_start_date, test_end_date (Date)

# Libraries and seed
library(lubridate)
library(dplyr)
library(sandwich)
library(lmtest)
library(ggplot2)

# Source plotting theme, you can also use theme_minimal. Theme link here: https://github.com/KarstenWalker5/R/blob/main/Helpers/themes.R
source("/Users/karstenwalker/Documents/GitHub/R/Helpers/themes.R")

set.seed(123)

# Sample BigQuery data load if you're not readng in a flat .csv
# Store your Google Cloud project ID
bq_auth()

project_id <- "compact-sylph-785"

# Define your SQL query (example using a public dataset)
sql <- "SELECT * 
        FROM `YOUR TABLE`"

# Run the query
tb <- bq_project_query(project_id, sql)

geo_test_data<- bq_table_download(tb)

# 1) Build DMA-level pre/post aggregates

analysis_df <- df2 %>%
  filter(date >= pre_start_date & date <= test_end_date) %>%
  mutate(period = ifelse(date <= pre_end_date, "pre", "post")) %>%
  group_by(dma, period) %>%
  summarise(
    sum_up = sum(upgrades, na.rm = TRUE),
    sum_au = sum(active_users, na.rm = TRUE),
    n_days = n(),
    .groups = "drop"
  ) %>%
  mutate(rate = ifelse(sum_au > 0, sum_up / sum_au, NA_real_)) %>%
  select(dma, period, sum_up, sum_au, n_days, rate) %>%
  pivot_wider(
    names_from = period,
    values_from = c(sum_up, sum_au, n_days, rate)
  ) %>%
  left_join(assignment, by = "dma") %>%
  mutate(treatment = ifelse(group == "TEST", 1, 0)) %>%
  filter(group %in% c("TEST", "CONTROL")) %>%
  filter(is.finite(rate_pre), is.finite(rate_post))

# 2) Balance check

t.test(rate_pre ~ treatment, data = analysis_df)

# Test 2: Pre-period slope parallel trends test
# Tests whether Test and Control DMAs have different linear time trends during the pre-period.
# A statistically significant treatment:time interaction suggests potential violation of parallel trends, which can bias DiD estimates.

pre_daily_df <- df2 %>%
  filter(date >= pre_start_date & date <= pre_end_date) %>%
  left_join(assignment, by = "dma") %>%
  filter(group %in% c("TEST", "CONTROL")) %>%
  mutate(
    treatment = ifelse(group == "TEST", 1, 0),
    time_index = as.numeric(as.Date(date) - min(as.Date(date))),
    rate = ifelse(active_users > 0, upgrades / active_users, NA_real_)
  ) %>%
  filter(is.finite(rate))

pre_trend_model <- lm(rate ~ treatment * time_index, data = pre_daily_df)
summary(pre_trend_model)

# Test 3: Cluster-robust inference for pre-period slope test
# Re-computes inference for the pre-period slope test using standard errors clustered at the DMA level.
# Clustering matches the unit of randomization and helps avoid overstating precision when daily observations within DMA are correlated.
library(sandwich)
library(lmtest)

coeftest(pre_trend_model, vcov = vcovCL(pre_trend_model, cluster = ~ dma))

# Test 4: Pre-period trend visualization
# Plots average daily rate over the pre-period for Test vs Control DMAs.
# Visual inspection complements the slope test and can reveal non-linear deviations from parallel trends.
library(ggplot2)

pre_avg <- pre_daily_df %>%
  group_by(date, treatment) %>%
  summarise(mean_rate = mean(rate, na.rm = TRUE), .groups = "drop")

ggplot(pre_avg, aes(x = date, y = mean_rate, color = factor(treatment))) +
  geom_line(size = 1.2) +
  labs(
    title = "Pre-Period Trends: Test vs Control (Rate Metric)",
    x = "Date",
    y = "Mean Daily Upgrade Rate",
    color = "Treatment (1=Test)"
  )

# Test 5: Denominator stability checks for rate metric
# Validates that active_users (the denominator) is balanced and evolves similarly in Test and Control during pre-period.
# If the denominator differs or trends differently, rate changes can reflect composition shifts rather than true numerator effects.
denom_pre_df <- df2 %>%
  filter(date >= pre_start_date & date <= pre_end_date) %>%
  left_join(assignment, by = "dma") %>%
  filter(group %in% c("TEST", "CONTROL")) %>%
  group_by(dma) %>%
  summarise(pre_active_users = sum(active_users, na.rm = TRUE), .groups = "drop") %>%
  left_join(assignment, by = "dma") %>%
  mutate(treatment = ifelse(group == "TEST", 1, 0))

t.test(pre_active_users ~ treatment, data = denom_pre_df)

denom_pre_daily <- df2 %>%
  filter(date >= pre_start_date & date <= pre_end_date) %>%
  left_join(assignment, by = "dma") %>%
  filter(group %in% c("TEST", "CONTROL")) %>%
  mutate(
    treatment = ifelse(group == "TEST", 1, 0),
    time_index = as.numeric(as.Date(date) - min(as.Date(date)))
  )

denom_slope_model <- lm(active_users ~ treatment * time_index, data = denom_pre_daily)
coeftest(denom_slope_model, vcov = vcovCL(denom_slope_model, cluster = ~ dma))

# Test 6: Pre/post day coverage check (composition check)
# Checks whether DMAs have consistent day coverage across the analysis window using the aggregated day counts.
# Unequal day coverage can distort pre/post aggregates and create artificial differences unrelated to treatment.
coverage_check <- analysis_df %>%
  transmute(
    dma,
    group,
    treatment,
    n_days_pre = n_days_pre,
    n_days_post = n_days_post
  )

coverage_check


# 3) Main model 

model_rate <- lm(rate_post ~ treatment + rate_pre, data = analysis_df)

summary(model_rate)

# Test 7: Heteroskedasticity-robust inference for main DiD regression (rate metric)
# Computes robust (HC1) standard errors for the ANCOVA-style rate regression (rate_post ~ treatment + rate_pre).
# This reduces sensitivity to heteroskedasticity across DMAs, which is common when markets differ in scale and volatility.
coeftest(model_rate, vcov = vcovHC(model_rate, type = "HC1"))

# Test 8: Cluster-robust inference for main DiD regression (rate metric)
# Computes standard errors clustered at the DMA level to align inference with DMA-level randomization.
# This guards against overstating precision when DMAs have correlated outcomes across days.
coeftest(model_rate, vcov = vcovCL(model_rate, cluster = ~ dma))


# Effect in rate units (Δ upgrades per active user per day)
beta_rate <- coef(model_rate)["treatment"]

# % lift in rate vs control post rate
avg_control_post_rate <- analysis_df %>%
  filter(treatment == 0) %>%
  summarise(mean(rate_post)) %>%
  pull()

pct_lift_rate <- beta_rate / avg_control_post_rate
pct_lift_rate

# 95% CI on rate effect + CI on % lift
ci_beta_rate <- confint(model_rate)["treatment", ]
ci_pct_lift_rate <- ci_beta_rate / avg_control_post_rate

ci_beta_rate
ci_pct_lift_rate

# Test 9: Placebo DiD test using only pre-period data (rate metric)
# Splits the pre-period into two halves and tests for a spurious treatment effect on a fake pre2-pre1 rate change.
# A significant placebo effect suggests underlying differences or artifacts unrelated to the intervention.
mid_pre <- as.Date(pre_start_date) + floor((as.Date(pre_end_date) - as.Date(pre_start_date)) / 2)

placebo_df <- df2 %>%
  filter(date >= pre_start_date & date <= pre_end_date) %>%
  mutate(period = ifelse(date <= mid_pre, "pre1", "pre2")) %>%
  group_by(dma, period) %>%
  summarise(
    sum_up = sum(upgrades, na.rm = TRUE),
    sum_au = sum(active_users, na.rm = TRUE),
    rate = ifelse(sum_au > 0, sum_up / sum_au, NA_real_),
    .groups = "drop"
  ) %>%
  select(dma, period, rate) %>%
  pivot_wider(names_from = period, values_from = rate) %>%
  left_join(assignment, by = "dma") %>%
  mutate(treatment = ifelse(group == "TEST", 1, 0)) %>%
  filter(group %in% c("TEST", "CONTROL")) %>%
  mutate(diff_placebo = pre2 - pre1) %>%
  filter(is.finite(diff_placebo))

t.test(diff_placebo ~ treatment, data = placebo_df)


# 4) Optional: Convert rate effect -> incremental upgrades over post period (helpful for business readout)

# Interpretation:
#   beta_rate is (Δ upgrades) / (active_user) over the post window (since we used sums over post window in rate_post)
# Because rate_post = sum_up_post / sum_au_post (over the post period),
# beta_rate * sum_au_post gives incremental upgrades over the post period (per DMA),
# and summing across TEST DMAs gives total incremental upgrades.

inc_upgrades_total <- analysis_df %>%
  filter(treatment == 1) %>%
  summarise(total_inc_upgrades = beta_rate * sum(sum_au_post, na.rm = TRUE)) %>%
  pull()

inc_upgrades_total

# CI version
inc_upgrades_ci <- analysis_df %>%
  filter(treatment == 1) %>%
  summarise(
    inc_low  = ci_beta_rate[1] * sum(sum_au_post, na.rm = TRUE),
    inc_high = ci_beta_rate[2] * sum(sum_au_post, na.rm = TRUE)
  )

inc_upgrades_ci


# 5) Equivalent DiD t-test on rate differences (sanity check)

analysis_df <- analysis_df %>%
  mutate(rate_diff = rate_post - rate_pre)

t.test(rate_diff ~ treatment, data = analysis_df)

# Test 10: Regression on (rate_post - rate_pre) with robust and clustered standard errors
# Estimates the DiD effect using the differenced rate outcome (rate_diff ~ treatment), which should align with the main regression.
# Robust and clustered SEs provide an additional check that inference is not sensitive to model specification.
model_diff <- lm(rate_diff ~ treatment, data = analysis_df)
summary(model_diff)
coeftest(model_diff, vcov = vcovHC(model_diff, type = "HC1"))
coeftest(model_diff, vcov = vcovCL(model_diff, cluster = ~ dma))

# Test 11: Randomization inference (permutation test) on the DiD rate_diff
# Computes a p-value by repeatedly re-assigning treatment labels and comparing permuted effects to the observed effect.
# This is especially appropriate for cluster-randomized experiments because it leverages the randomization directly.
set.seed(123)
B <- 2000

observed_effect <- coef(model_diff)["treatment"]

perm_effects <- replicate(B, {
  fake_treatment <- sample(analysis_df$treatment)
  coef(lm(rate_diff ~ fake_treatment, data = analysis_df))[2]
})

perm_p_value <- mean(abs(perm_effects) >= abs(observed_effect))
perm_p_value

# Test 12: Influence / outlier diagnostics at the DMA level
# Checks whether a small number of DMAs have outsized influence on the estimated treatment effect.
# If a few DMAs dominate, consider reporting sensitivity analyses that exclude them or winsorize extreme outcomes.
influence_measures(model_rate)
plot(cooks.distance(model_rate), type = "h", main = "Cook's Distance for DMA-Level DiD Regression (Rate Metric)")


# ROAS estimation
# Make sure spend_df includes only that channel’s spend and that your “spend removed” is what would have been spent absent the test. 
# If the platform reallocates budget automatically, you need to capture the counterfactual spend (often “planned budget”) rather than observed.
# Estimate ROAS from the DiD output by turning the estimated incremental impact into incremental value, then dividing by incremental spend.
# Because our test is a spend shutoff (treatment = $0), the most natural “baseline ROAS” we’ll report is actually:
#  * ROAS of the channel = (incremental value lost when shutting it off) / (spend removed), which is equivalent to incremental value per $1 of spend.
# 
# Requires:
# - beta_rate, ci_beta_rate from model_rate
# - analysis_df contains sum_au_post and treatment==1 for test DMAs
# - spend_df: dma, date, spend  (channel spend)
# - assignment: dma, group
# - test_start_date, test_end_date
# - value_per_upgrade (set below)

# 1) Total post active_users in test DMAs
test_post_active_users <- analysis_df %>%
  filter(treatment == 1) %>%
  summarise(total_au_post = sum(sum_au_post, na.rm = TRUE)) %>%
  pull(total_au_post)

# 2) Incremental upgrades over post period (point estimate + CI)
inc_upgrades <- beta_rate * test_post_active_users
inc_upgrades_ci <- c(ci_beta_rate[1], ci_beta_rate[2]) * test_post_active_users

# 3) Spend removed in test DMAs over post period
spend_test_post <- spend_df %>%
  filter(date >= test_start_date & date <= test_end_date) %>%
  inner_join(assignment %>% filter(group == "TEST"), by = "dma") %>%
  summarise(total_spend = sum(spend, na.rm = TRUE)) %>%
  pull(total_spend)

# 4) Value per upgrade, you have to define this
# Examples:
# - expected LTV per upgrade
# - gross profit LTV per upgrade
value_per_upgrade <- 100  # <-- CHANGE

# 5) Incremental value (point estimate + CI)
inc_value <- inc_upgrades * value_per_upgrade
inc_value_ci <- inc_upgrades_ci * value_per_upgrade

# 6) ROAS
# For shutoff tests, inc_value is often negative (lost value).
# "Baseline ROAS of channel" is typically reported as positive: -inc_value / spend_removed
roas_point <- -inc_value / spend_test_post
roas_ci <- -inc_value_ci / spend_test_post

tibble(
  spend_removed_test_post = spend_test_post,
  inc_upgrades_post = inc_upgrades,
  inc_upgrades_post_ci_low = inc_upgrades_ci[1],
  inc_upgrades_post_ci_high = inc_upgrades_ci[2],
  value_per_upgrade = value_per_upgrade,
  inc_value_post = inc_value,
  inc_value_post_ci_low = inc_value_ci[1],
  inc_value_post_ci_high = inc_value_ci[2],
  roas_point = roas_point,
  roas_ci_low = roas_ci[1],
  roas_ci_high = roas_ci[2]
)

# If roas_point = 2.5: every $1 of channel spend drives ~$2.50 of incremental value.
# If the estimate is noisy, the CI might cross 0 → interpret as “not statistically distinguishable”.

# Test 13: Optional spillover / contamination robustness check (requires adjacency or proximity definition)
# Tests sensitivity to possible spillovers where exposure in BAU or nearby treated DMAs could affect control outcomes.
# Provide a DMA adjacency mapping (neighbors) and optionally exclude controls adjacent to treated DMAs, then re-run the primary model.
# Example expected input: neighbors_df with columns (dma, neighbor_dma)
# If you do not have adjacency data, skip this section.

# neighbors_df <- read.csv("dma_neighbors.csv")  # <-- provide your mapping
# treated_dmas <- assignment %>% filter(group == "TEST") %>% pull(dma)
# adjacent_to_treated <- neighbors_df %>% filter(neighbor_dma %in% treated_dmas) %>% pull(dma) %>% unique()
#
# analysis_df_no_adjacent <- analysis_df %>%
#   filter(!(dma %in% adjacent_to_treated))
#
# model_no_adjacent <- lm(rate_post ~ treatment + rate_pre, data = analysis_df_no_adjacent)
# summary(model_no_adjacent)
# coeftest(model_no_adjacent, vcov = vcovCL(model_no_adjacent, cluster = ~ dma))


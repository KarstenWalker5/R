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

library(dplyr)
library(tidyr)

# 1) Build DMA-level pre/post aggregates (CHANGED: keep sums + compute rate)

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
  # --- CHANGE: compute daily rate using sums (more stable than mean(up/au)) ---
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

# 2) Balance check (CHANGED: check pre RATE balance)

t.test(rate_pre ~ treatment, data = analysis_df)

# 3) Main model (CHANGED: regression on rate with pre-rate adjustment)

model_rate <- lm(rate_post ~ treatment + rate_pre, data = analysis_df)

summary(model_rate)

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
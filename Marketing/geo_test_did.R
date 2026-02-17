
# Assumes
# * Randomization was at DMA level
# * Outcome is aggregated at DMA level
# * Observations are independent at DMA level
# * Pre-period included as covariate improves precision

# Report:
# * Incremental upgrades per DMA per day
# * Total incremental upgrades over test
# * % lift
# * 95% CI
# * p-value
# * Example executive sentence: Shutting off the channel reduced upgrades by 8.4% (95% CI: 4.1%–12.7%, p=0.01).

# df structure must contain
# * assignment  # dma + group
# * df2         # daily upgrades
# * test_start_date
# * test_end_date
# * pre_start_date
# * pre_end_date

# Why we use DiD and not something "fancier"
# * This is a classic randomized cluster experiment.
# * Randomized DMA assignment
# * 25–68 treated DMAs
# * Clean pre-period
# * Balanced arms
# * No staggered timing

# Build DMA level dataset
library(dplyr)

analysis_df <- df2 %>%
  filter(date >= pre_start_date & date <= test_end_date) %>%
  mutate(period = ifelse(date <= pre_end_date, "pre", "post")) %>%
  group_by(dma, period) %>%
  summarise(mean_up = mean(upgrades, na.rm=TRUE), .groups="drop") %>%
  tidyr::pivot_wider(names_from = period, values_from = mean_up) %>%
  left_join(assignment, by="dma") %>%
  mutate(treatment = ifelse(group=="TEST", 1, 0)) %>%
  filter(group %in% c("TEST","CONTROL"))

# Run before test,should not be statistically different. If it is, your assignment wasn’t well balanced.
t.test(pre ~ treatment, data = analysis_df)

# run regression
model <- lm(post ~ treatment + pre, data = analysis_df)

summary(model)

# convert to % lift
beta <- coef(model)["treatment"]

avg_control_post <- analysis_df %>%
  filter(treatment==0) %>%
  summarise(mean(post)) %>%
  pull()

pct_lift <- beta / avg_control_post
pct_lift

confint(model)["treatment", ]

# Gives same result
analysis_df <- analysis_df %>%
  mutate(diff = post - pre)

t.test(diff ~ treatment, data = analysis_df)

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
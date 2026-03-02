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
# * signups         # daily upgrades
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

# Build DMA level dataset

analysis_df <- df2 %>%
  filter(date >= pre_start_date & date <= test_end_date) %>%
  mutate(period = ifelse(date <= pre_end_date, "pre", "post")) %>%
  group_by(dma, period) %>%
  summarise(mean_up = mean(upgrades, na.rm=TRUE), .groups="drop") %>%
  tidyr::pivot_wider(names_from = period, values_from = mean_up) %>%
  left_join(assignment, by="dma") %>%
  mutate(treatment = ifelse(group=="TEST", 1, 0)) %>%
  filter(group %in% c("TEST","CONTROL"))

# Run before test, should not be stat sig. If it is, your assignment wasn’t well balanced.
t.test(pre ~ treatment, data = analysis_df)

# Test 2: Pre-period slope parallel trends test
# Tests whether Test and Control DMAs have different linear time trends during the pre-period.
# A statistically significant treatment:time interaction suggests potential violation of parallel trends, which can bias DiD estimates.

pre_daily_df <- df2 %>%
  filter(date >= pre_start_date & date <= pre_end_date) %>%
  left_join(assignment, by = "dma") %>%
  filter(group %in% c("TEST", "CONTROL")) %>%
  mutate(
    treatment = ifelse(group == "TEST", 1, 0),
    time_index = as.numeric(as.Date(date) - min(as.Date(date)))
  )

pre_trend_model <- lm(upgrades ~ treatment * time_index, data = pre_daily_df)
summary(pre_trend_model)

# Test 3: Cluster-robust inference for pre-period slope test
# This re-computes inference for the pre-period slope test using standard errors clustered at the DMA level.
# Clustering matches the unit of randomization and helps avoid overstating precision when outcomes are correlated within DMA.

coeftest(pre_trend_model, vcov = vcovCL(pre_trend_model, cluster = ~ dma))

# Test 4: Pre-period trend visualization
# This plots average daily upgrades over the pre-period for Test vs Control DMAs.
# Visual inspection complements the slope test and can reveal non-linear or localized deviations from parallel trends.

pre_avg <- pre_daily_df %>%
  group_by(date, treatment) %>%
  summarise(mean_up = mean(upgrades, na.rm = TRUE), .groups = "drop")

ggplot(pre_avg, aes(x = date, y = mean_up, color = factor(treatment))) +
  geom_line(size = 1.2) +
  labs(
    title = "Pre-Period Trends: Test vs Control",
    x = "Date",
    y = "Mean Daily Upgrades",
    color = "Treatment (1=Test)"
  )+
  theme_fancy()

# Test 5: Pre/post day coverage check 
# This checks whether DMAs have consistent day coverage across the analysis window.
# Missing days can distort pre/post means and create artificial differences unrelated to treatment.
coverage_check <- df2 %>%
  filter(date >= pre_start_date & date <= test_end_date) %>%
  mutate(period = ifelse(date <= pre_end_date, "pre", "post")) %>%
  group_by(dma, period) %>%
  summarise(n_days = n(), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = period, values_from = n_days) %>%
  left_join(assignment, by = "dma") %>%
  filter(group %in% c("TEST", "CONTROL"))

coverage_check

# run regression
model <- lm(post ~ treatment + pre, data = analysis_df)

summary(model)

# Test 6: Heteroskedasticity-robust inference for main DiD regression
# This computes robust HC1 standard errors for the main ANCOVA-style DiD regression (post ~ treatment + pre).
# It reduces sensitivity to heteroskedasticity across DMAs, which is common when markets differ in scale and volatility.
coeftest(model, vcov = vcovHC(model, type = "HC1"))

# Test 7: Cluster-robust inference for main DiD regression
# This computes standard errors clustered at the DMA level for the main regression.
# Clustering aligns inference with the DMA-level randomization and guards against within-DMA correlation.
coeftest(model, vcov = vcovCL(model, cluster = ~ dma))

# convert to % lift
beta <- coef(model)["treatment"]

avg_control_post <- analysis_df %>%
  filter(treatment==0) %>%
  summarise(mean(post)) %>%
  pull()

pct_lift <- beta / avg_control_post

pct_lift

confint(model)["treatment", ]

# Test 8: Placebo DiD test using only pre-period data
# This splits the pre-period into two halves and runs a "fake" DiD (pre2 - pre1) between Test and Control.
# A significant placebo effect suggests underlying differences unrelated to the intervention (e.g., non-parallel pre trends or data artifacts).
mid_pre <- as.Date(pre_start_date) + floor((as.Date(pre_end_date) - as.Date(pre_start_date)) / 2)

placebo_df <- df2 %>%
  filter(date >= pre_start_date & date <= pre_end_date) %>%
  mutate(period = ifelse(date <= mid_pre, "pre1", "pre2")) %>%
  group_by(dma, period) %>%
  summarise(mean_up = mean(upgrades, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = period, values_from = mean_up) %>%
  left_join(assignment, by = "dma") %>%
  mutate(treatment = ifelse(group == "TEST", 1, 0)) %>%
  filter(group %in% c("TEST", "CONTROL")) %>%
  mutate(diff_placebo = pre2 - pre1)

t.test(diff_placebo ~ treatment, data = placebo_df)

# Gives same result
analysis_df <- analysis_df %>%
  mutate(diff = post - pre)

t.test(diff ~ treatment, data = analysis_df)

# Test 9: Regression on post-pre with robust and clustered standard errors
# This estimates the same DiD effect using the differenced outcome (diff ~ treatment), which should align closely with the main regression.
# Reporting robust and clustered SEs provides an additional check that inference is not sensitive to model specification.

model_diff <- lm(diff ~ treatment, data = analysis_df)

summary(model_diff)

coeftest(model_diff, vcov = vcovHC(model_diff, type = "HC1"))

coeftest(model_diff, vcov = vcovCL(model_diff, cluster = ~ dma))

# Test 10: Permutation test on the DiD diff
# This computes a p-value by repeatedly re-assigning treatment labels and comparing the permuted effect to the observed effect.
# It is especially appropriate for cluster-randomized experiments because it relies directly on the randomization scheme rather than asymptotic approximations.
B <- 2000

observed_effect <- coef(model_diff)["treatment"]

perm_effects <- replicate(B, {
  fake_treatment <- sample(analysis_df$treatment)
  coef(lm(diff ~ fake_treatment, data = analysis_df))[2]
})

perm_p_value <- mean(abs(perm_effects) >= abs(observed_effect))

perm_p_value

# Test 11: Outlier diagnostics at the DMA level
# Checks whether a small number of DMAs have outsized influence on the estimated treatment effect.
# If a few DMAs dominate, consider reporting sensitivity analyses that exclude them or winsorize extreme outcomes.
influence_measures(model)

plot(cooks.distance(model), type = "h", main = "Cook's Distance for DMA-Level DiD Regression")

# ROAS estimation
# Make sure spend_df includes only that channel’s spend and that “spend removed” is what would have been spent absent the test. 
# If the platform reallocates budget automatically, you need to capture the counterfactual spend (often “planned budget”) rather than observed.
# Estimate ROAS from the DiD output by turning the estimated incremental impact into incremental value, then dividing by incremental spend.
# Because our test is a spend shutoff w/ treatment of $0, the most natural “baseline ROAS” we’ll report is actually
#  * ROAS of the channel = (incremental value lost when shutting it off) / (spend removed), which is equivalent to incremental value per $1 of spend.
# 
# Requires:
# * beta_rate, ci_beta_rate from model_rate
# * analysis_df contains sum_au_post and treatment==1 for test DMAs
# * spend_df: dma, date, spend  (channel spend)
# * assignment: dma, group
# * test_start_date, test_end_date
# *value_per_upgrade (set below)

# 1. Total post active_users in test DMAs
test_post_active_users <- analysis_df %>%
  filter(treatment == 1) %>%
  summarise(total_au_post = sum(sum_au_post, na.rm = TRUE)) %>%
  pull(total_au_post)

# 2. Incremental upgrades over post period (point estimate + CI)
inc_upgrades <- beta_rate * test_post_active_users
inc_upgrades_ci <- c(ci_beta_rate[1], ci_beta_rate[2]) * test_post_active_users

# 3. Spend removed in test DMAs over post period
spend_test_post <- spend_df %>%
  filter(date >= test_start_date & date <= test_end_date) %>%
  inner_join(assignment %>% filter(group == "TEST"), by = "dma") %>%
  summarise(total_spend = sum(spend, na.rm = TRUE)) %>%
  pull(total_spend)

# 4.  Value per upgrade, you have to define this
# Examples:
# * expected LTV per upgrade
# * gross profit LTV per upgrade

value_per_upgrade <- 100  # <-- CHANGE

# 5. Incremental value (point estimate + CI)
inc_value <- inc_upgrades * value_per_upgrade

inc_value_ci <- inc_upgrades_ci * value_per_upgrade

# 6. ROAS
# For spend-down tests, inc_value is often negative (lost value).
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

# Test 12: Optional spillover / contamination robustness check (requires adjacency or proximity definition)
# This tests sensitivity to possible spillovers where spend in BAU or nearby treated DMAs could affect control outcomes.
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
# model_no_adjacent <- lm(post ~ treatment + pre, data = analysis_df_no_adjacent)
# summary(model_no_adjacent)
# coeftest(model_no_adjacent, vcov = vcovCL(model_no_adjacent, cluster = ~ dma))
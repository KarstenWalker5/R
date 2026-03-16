library(tidyverse)
library(lubridate)
library(CausalImpact)
library(broom)

# Params
rollout_date <- as.Date("2026-01-05")

pre_period_end <- rollout_date - 1

#Load data
teacher_data<-read.csv(file="/Users/karstenwalker/Downloads/new_teacher_data.csv")%>%
  rename(date=dt)%>%
  mutate(date = as.Date(date)) %>%
  arrange(date)

###### Cannibalization ######
rollout_date <- as.Date("2026-01-05")

teacher_data <- teacher_data %>%
  mutate(date = as.Date(date)) %>%
  arrange(date) %>%
  mutate(
    dow = factor(lubridate::wday(date, week_start = 1)),
    woy = factor(lubridate::isoweek(date)),
    log_paid = log(paid)
  )

pre_period  <- c(min(teacher_data$date), rollout_date - 1)

post_period <- c(rollout_date, max(teacher_data$date))

X <- model.matrix(~ dow + woy, data = teacher_data)

# Drop intercept + force safe names for BSTS
X <- X[, -1, drop = FALSE]

colnames(X) <- make.names(colnames(X), unique = TRUE)

stopifnot(!any(colnames(X) == ""), !any(is.na(colnames(X))))

impact_mat <- zoo(cbind(log_paid = teacher_data$log_paid, X), order.by = teacher_data$date)

impact_paid <- CausalImpact(
  impact_mat,
  pre_period,
  post_period,
  model.args = list(nseasons = 7, season.duration = 1))

summary(impact_paid)

summary(impact_paid, "report")

# Extract results to simple plot
paid_results <- as.data.frame(impact_paid$series) %>%
  rownames_to_column("date") %>%
  mutate(date = as.Date(date))%>%
  filter(date>="2025-08-01")

ggplot(paid_results, aes(x = date)) +
  geom_line(aes(y = response), size = 1) +
  geom_line(aes(y = point.pred), linetype = "dashed") +
  geom_ribbon(aes(ymin = point.pred.lower,
                  ymax = point.pred.upper),
              alpha = 0.2) +
  geom_vline(xintercept = rollout_date,
             linetype = "dotted") +
  labs(title = "Cannibalizaton of Verification on Paid Teachers",
       y = "Paid Teachers",
       x = NULL) +
  theme_minimal()

# Using base plot
cutoff_date <- as.Date("2025-10-01")

# Pull the CausalImpact series and filter dates
ci_df <- as.data.frame(impact_paid$series) %>%
  rownames_to_column("date") %>%
  mutate(date = as.Date(date)) %>%
  filter(date >= cutoff_date)

# Theme
causal_theme <- ggthemes::theme_fivethirtyeight() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold", color = "#ebba34"),
    plot.subtitle = element_text(hjust = 0.5, size = 16, face = "italic"),
    legend.position = "bottom",
    legend.text = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 12, face = "bold"),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_blank(),
    plot.background = element_rect(fill = "white"),
    panel.background = element_rect(fill = "white"),
    legend.background = element_rect(fill = "white")
  )

# Original (observed vs predicted)
p1 <- ggplot(ci_df, aes(x = date)) +
  geom_ribbon(aes(ymin = point.pred.lower, ymax = point.pred.upper), alpha = 0.25) +
  geom_line(aes(y = point.pred, linetype = "Predicted"), linewidth = 0.8) +
  geom_line(aes(y = response, linetype = "Observed"), linewidth = 0.9) +
  geom_vline(xintercept = rollout_date, linetype = "dotted") +
  scale_linetype_manual(values = c("Observed" = "solid", "Predicted" = "dashed")) +
  labs(
    x = "Date",
    y = "Weighted Avg Ranking",
    title = "CausalImpact Estimate of \nTeacher Verification Cannibalization") +
  causal_theme

# Pointwise effect 
p2 <- ggplot(ci_df, aes(x = date)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_ribbon(aes(ymin = point.effect.lower, ymax = point.effect.upper), alpha = 0.25) +
  geom_line(aes(y = point.effect), linewidth = 0.9) +
  geom_vline(xintercept = rollout_date, linetype = "dotted") +
  labs(x = NULL, y = "Pointwise Effect", title = NULL) +
  causal_theme+
  theme(legend.position = "none")

# Cumulative effect
p3 <- ggplot(ci_df, aes(x = date)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_ribbon(aes(ymin = cum.effect.lower, ymax = cum.effect.upper), alpha = 0.25) +
  geom_line(aes(y = cum.effect), linewidth = 0.9) +
  geom_vline(xintercept = rollout_date, linetype = "dotted") +
  labs(x = NULL, y = "Cumulative Effect", title = NULL) +
  causal_theme +
  theme(legend.position = "none")

# Stack panels
(p1 / p2 / p3) + plot_layout(heights = c(1.2, 1, 1))

###### Upgrade Rate ######
teacher_data <- teacher_data %>%
  mutate(upgrade_rate = paid / (paid + free))

pre_period  <- c(min(teacher_data$date), rollout_date - 1)

post_period <- c(rollout_date, max(teacher_data$date))

# Calendar covariates
X <- model.matrix(~ dow + woy, data = teacher_data)
X <- X[, -1, drop = FALSE]
colnames(X) <- make.names(colnames(X), unique = TRUE)

impact_mat <- zoo(
  cbind(upgrade_rate = teacher_data$upgrade_rate, X),
  order.by = teacher_data$date
)

impact_upgrade <- CausalImpact(
  impact_mat,
  pre_period,
  post_period,
  model.args = list(
    nseasons = 7,
    season.duration = 1))

summary(impact_upgrade)

# Paid ↓, Verified ↑, Total flat	Strong substitution
# Paid flat, Verified ↑, Total up, Mild substition
# Paid ↓, Verified ↑, Total ↑	Partial substitution
# Paid flat but on trend, Verified ↑, Total ↑	Mostly incremental
# Paid ↑ slower than trend	Soft cannibalization

# Extract results to simple plot
paid_results <- as.data.frame(impact_upgrade$series) %>%
  rownames_to_column("date") %>%
  mutate(date = as.Date(date))%>%
  filter(date>="2025-10-01")

ggplot(paid_results, aes(x = date)) +
  geom_line(aes(y = response), size = 1) +
  geom_line(aes(y = point.pred), linetype = "dashed") +
  geom_ribbon(aes(ymin = point.pred.lower,
                  ymax = point.pred.upper),
              alpha = 0.2) +
  geom_vline(xintercept = rollout_date,
             linetype = "dotted") +
  labs(title = "Cannibalizaton of Verification on Paid Teachers",
       y = "Paid Teachers",
       x = NULL) +
  theme_minimal()

# Using base plot
cutoff_date <- as.Date("2025-10-01")

# Pull the CausalImpact series and filter dates
ci_df <- as.data.frame(impact_upgrade$series) %>%
  rownames_to_column("date") %>%
  mutate(date = as.Date(date)) %>%
  filter(date >= cutoff_date)

# Original (observed vs predicted)
p1 <- ggplot(ci_df, aes(x = date)) +
  geom_ribbon(aes(ymin = point.pred.lower, ymax = point.pred.upper), alpha = 0.25) +
  geom_line(aes(y = point.pred, linetype = "Predicted"), linewidth = 0.8) +
  geom_line(aes(y = response, linetype = "Observed"), linewidth = 0.9) +
  geom_vline(xintercept = rollout_date, linetype = "dotted") +
  scale_linetype_manual(values = c("Observed" = "solid", "Predicted" = "dashed")) +
  labs(
    x = "Date",
    y = "Weighted Avg Ranking",
    title = "CausalImpact Estimate of \nTeacher Verification on Upgrade Rate") +
  causal_theme

# Pointwise effect 
p2 <- ggplot(ci_df, aes(x = date)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_ribbon(aes(ymin = point.effect.lower, ymax = point.effect.upper), alpha = 0.25) +
  geom_line(aes(y = point.effect), linewidth = 0.9) +
  geom_vline(xintercept = rollout_date, linetype = "dotted") +
  labs(x = NULL, y = "Pointwise Effect", title = NULL) +
  causal_theme+
  theme(legend.position = "none")

# Cumulative effect
p3 <- ggplot(ci_df, aes(x = date)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_ribbon(aes(ymin = cum.effect.lower, ymax = cum.effect.upper), alpha = 0.25) +
  geom_line(aes(y = cum.effect), linewidth = 0.9) +
  geom_vline(xintercept = rollout_date, linetype = "dotted") +
  labs(x = NULL, y = "Cumulative Effect", title = NULL) +
  causal_theme +
  theme(legend.position = "none")

# Stack panels
(p1 / p2 / p3) + plot_layout(heights = c(1.2, 1, 1))


###### Premium Access ######
teacher_data <- teacher_data %>%
  mutate(premium_access = paid + verified)

pre_period  <- c(min(teacher_data$date), rollout_date - 1)

post_period <- c(rollout_date, max(teacher_data$date))

# Calendar covariates
X <- model.matrix(~ dow + woy, data = teacher_data)
X <- X[, -1, drop = FALSE]
colnames(X) <- make.names(colnames(X), unique = TRUE)

impact_mat <- zoo(
  cbind(premium_access = teacher_data$premium_access, X),
  order.by = teacher_data$date
)

impact_premium <- CausalImpact(
  impact_mat,
  pre_period,
  post_period,
  model.args = list(
    nseasons = 7,
    season.duration = 1))

summary(impact_premium)

###### Substitution Effect ######

# prep data
df <- teacher_data %>%
  mutate(date = as.Date(date)) %>%
  arrange(date) %>%
  mutate(total_teachers = free + paid + verified)

pre_period  <- c(min(df$date), rollout_date - 1)

post_period <- c(rollout_date, max(df$date))

# paid impact using log scale and calendar covariates
X_paid <- model.matrix(~ dow + woy, data = df)

X_paid <- X_paid[, -1, drop = FALSE]

colnames(X_paid) <- make.names(colnames(X_paid), unique = TRUE)

paid_impact_mat <- zoo(
  cbind(log_paid = df$log_paid, X_paid),
  order.by = df$date)

impact_paid <- CausalImpact(
  paid_impact_mat,
  pre_period,
  post_period,
  model.args = list(nseasons = 7, season.duration = 1))

print(summary(impact_paid))

# extract paid effects (log scale)
paid_shortfall_cum_log <- as.numeric(impact_paid$summary["Cumulative", "AbsEffect"])

paid_shortfall_pct_log <- as.numeric(impact_paid$summary["Cumulative", "RelEffect"]) * 100
# Paid grew by 8.4%

paid_shortfall_avg_log <- as.numeric(impact_paid$summary["Average", "AbsEffect"])

paid_shortfall_avg_pct_log <- as.numeric(impact_paid$summary["Average", "RelEffect"]) * 100

# substitution accounting using level growth comparison
post_df <- df %>% 
  filter(date >= rollout_date)

verified_growth <- with(post_df, last(verified) - first(verified))

paid_growth_obs <- with(post_df, last(paid) - first(paid))

# substitution ratio is only meaningful if paid is below counterfactual aka negative.
# Our paid is positive so this should return 0
substitution_ratio <- if (paid_shortfall_pct_log < 0 && verified_growth > 0) {
  baseline_paid <- df %>% 
    filter(date == rollout_date) %>% 
    pull(paid) %>% 
    .[1]
  
  approx_paid_shortfall_units <- abs(
    (paid_shortfall_pct_log / 100) * baseline_paid * nrow(post_df)
  )
  
  approx_paid_shortfall_units / verified_growth
} else {
  0
}

# descriptive paid share post rollout
post_paid_share <- post_df %>%
  mutate(paid_share = paid / (paid + verified)) %>%
  summarise(
    paid_share_start = first(paid_share),
    paid_share_end   = last(paid_share),
    paid_share_delta = paid_share_end - paid_share_start
  )

# incrementality check using total teachers
X_total <- model.matrix(~ dow + woy, data = df)

X_total <- X_total[, -1, drop = FALSE]

colnames(X_total) <- make.names(colnames(X_total), unique = TRUE)

total_impact_mat <- zoo(
  cbind(total_teachers = df$total_teachers, X_total),
  order.by = df$date
)

impact_total <- CausalImpact(
  total_impact_mat,
  pre_period,
  post_period,
  model.args = list(nseasons = 7, season.duration = 1)
)

summary(impact_total)

total_effect_cum <- as.numeric(impact_total$summary["Cumulative", "AbsEffect"])

total_effect_pct <- as.numeric(impact_total$summary["Cumulative", "RelEffect"]) * 100

# summary table
results <- tibble(
  rollout_date = rollout_date,
  pre_start = pre_period[1],
  pre_end   = pre_period[2],
  post_start = post_period[1],
  post_end   = post_period[2],
  paid_effect_cumulative_log = paid_shortfall_cum_log,
  paid_effect_cumulative_pct = paid_shortfall_pct_log,
  paid_effect_average_log = paid_shortfall_avg_log,
  paid_effect_average_pct = paid_shortfall_avg_pct_log,
  verified_growth_post = verified_growth,
  observed_paid_growth_post = paid_growth_obs,
  substitution_ratio = substitution_ratio,
  paid_share_start_post = post_paid_share$paid_share_start,
  paid_share_end_post   = post_paid_share$paid_share_end,
  paid_share_change_post = post_paid_share$paid_share_delta,
  total_teachers_effect_cumulative = total_effect_cum,
  total_teachers_effect_relative_pct = total_effect_pct
)

print(results)

###### Summary ######

# We find no evidence of paid cannibalization following Verified launch. 
# Paid subscriptions are estimated to be ~8% above counterfactual trajectory post-launch. 
# Verified growth does not correspond to a decline in paid.

# Boilerplate Interrupted Time Series Template. Same pre-checks as DiD and TWFE, but is used when
# we care about the characteristics of the entire dataset pre/post not just the aggregate conclusion.

###### Steps ######
# 1. Basic visual inspection of pre/post campaign to examine parallel trends assumption.
# 2. Analysis of variance in pre/post sales data.
# 3. Optional: Event study comparison to determine the difference between between pre/post, which helps us
#    determine if the difference is a function of surrounding time periods OR the intervention.
# 4. We use linear regression to compare pre-intervention trends between groups.
# 5. Finally, we use TWFE to estimate if the campaign had an impact on sales. We can do this with or without
#    spend as an input.
# 6. Our last step is a light weight intro to Hetergeneous Treatment Effects (HTE), weighting, and additional
#    methods for improving the power of our estimand.

library(dplyr)
library(tidyr)
library(lubridate)
library(lmtest)
library(sandwich)
library(fixest)
library(car)
library(feasts)
library(infer)
library(corrr)
library(CausalImpact)
library(zoo)

set.seed(7)

# Theme
theme_fancy <- function() {
  theme_minimal(base_family = "Asap Condensed") +
    theme(panel.grid.minor = element_blank(),
          legend.position = "bottom")
}

# Load causal pre-check function
source("/Users/karstenwalker/Documents/R/Causal/causal_pre_check.R")

# Load optional group by functions
source("/Users/karstenwalker/Documents/R/Causal/causal_by_group_functions_dev.R")

## Generate synthetic data for example
# Params
n_days = 200
n_cities = 10
start_date = as.Date("2024-01-01")
intervention_date = start_date + 100

treated_cities = paste0("City_", 1:5)
control_cities = paste0("City_", 6:10)

beta_spend = 2
spend_jump = 8
epsilon_sd = 3
beta_treat_direct = 4            # deterministic immediate/direct lift
beta_city_shock = 5              # idiosyncratic city×date lift (survives FE)

sales_data =
  expand_grid(
    date = seq.Date(start_date, by = "day", length.out = n_days),
    city = paste0("City_", 1:n_cities)
  ) %>%
  mutate(
    city = factor(city),
    treat = if_else(city %in% treated_cities, 1L, 0L),
    post = if_else(date >= intervention_date, 1L, 0L),
    day_index = as.integer(date - start_date),
    city_effect = rnorm(n_cities, sd = 4)[as.integer(city)],
    common_trend = 0.05 * day_index,
    spend = rnorm(n(), mean = 20, sd = 2),
    epsilon = rnorm(n(), sd = epsilon_sd)
  ) %>%
  group_by(post, treat) %>%
  mutate(spend_grp_mean = mean(spend)) %>%
  ungroup() %>%
  group_by(post) %>%
  mutate(spend_overall_pre = if_else(post == 0, mean(spend[post == 0]), NA_real_)) %>%
  ungroup() %>%
  mutate(
    spend_balanced = if_else(post == 0,
                             spend - (spend_grp_mean - spend_overall_pre),
                             spend
    ),
    spend = spend_balanced + spend_jump * treat * post,
    treat_city_shock = beta_city_shock * rnorm(n(), mean = 1, sd = 0.2) * treat * post,
    sales =
      50 +
      city_effect +
      common_trend +
      beta_spend * spend +
      beta_treat_direct * treat * post +
      treat_city_shock +
      epsilon,
    rel_time = as.integer(date - intervention_date)
  ) %>%
  select(date, city, treat, post, rel_time, spend, sales)%>%
  mutate(treat=as.factor(treat))

# Inspect
glimpse(sales_data)

head(sales_data)

report::report(sales_data%>%
                 select(-date))

# ALWAYS check for nulls
sales_data %>%
  group_by(treat) %>%
  summarise(pct_na = across(everything(), ~ sum(is.na(.x)) / n()),
            num_na = across(everything(), ~ sum(is.na(.x))))

## Summary stats

# Pre/post
sales_data %>%
  group_by(treat) %>%
  summarize(
    mean_sales = mean(sales, na.rm = TRUE),
    med_sales = median(sales, na.rm = TRUE),
    max_sales= max(sales, na.rm = TRUE),
    sd_sales = sd(sales),
    min_date = min(date),
    max_date = max(date)
  )

# Pre/post by city
sales_data %>%
  group_by(treat, city) %>%
  summarize(
    mean_sales = mean(sales, na.rm = TRUE),
    med_sales = median(sales, na.rm = TRUE),
    min_sales= min(sales, na.rm = TRUE),
    max_sales= max(sales, na.rm = TRUE),
    sd_sales = sd(sales),
    min_date = min(date),
    max_date = max(date)
  )%>%
  arrange(city, treat)

# Aggregate pre/post
sales_data%>%
  group_by(date)%>%
  summarize(sales=sum(sales))%>%
  ggplot(aes(date, sales)) +
  geom_line() +
  geom_vline(xintercept = intervention_date, linetype = "dashed", color='red')+
  labs(title = "Sales Time Series")+
  theme_fancy()

## Parallel Trends Check, Event Study
# Basic viz
sales_data %>%
  mutate(group = if_else(treat == 1, "Treated", "Control")) %>%
  group_by(date, group) %>%
  summarize(avg_sales = mean(sales), .groups = "drop") %>%
  ggplot(aes(x = date, y = avg_sales, color = group)) +
  geom_line(linewidth = 1) +
  geom_vline(xintercept = intervention_date, linetype = "dashed") +
  labs(
    title = "Parallel Trends Check: Treated vs Control",
    x = "Date",
    y = "Average Daily Sales",
    color = "Group"
  ) +
  theme_fancy()

# Add smoothing
sales_data %>%
  mutate(group = if_else(treat == 1, "Treated", "Control")) %>%
  group_by(date, group) %>%
  summarize(avg_sales = mean(sales), .groups = "drop") %>%
  ggplot(aes(x = date, y = avg_sales, color = group)) +
  geom_line(linewidth = 1) +
  geom_vline(xintercept = intervention_date, linetype = "dashed") +
  geom_smooth()+
  labs(
    title = "Parallel Trends Check: Treated vs Control",
    x = "Date",
    y = "Average Daily Sales",
    color = "Group"
  ) +
  theme_fancy()

# By city
sales_data %>%
  mutate(group = if_else(treat == 1, "Treated", "Control")) %>%
  ggplot(aes(x = date, y = sales, color = group)) +
  geom_line(alpha = 0.6) +
  geom_vline(xintercept = intervention_date, linetype = "dashed") +
  facet_wrap(~ city, ncol = 5) +
  labs(
    title = "City-Level Sales Paths by Treatment Group",
    x = "Date",
    y = "Daily Sales",
    color = "Group"
  ) +
  theme_fancy()

# Spend
sales_data %>%
  mutate(group = if_else(treat==1,"Treated","Control")) %>%
  group_by(date, group) %>%
  summarize(avg_spend = mean(spend), .groups="drop") %>%
  ggplot(aes(date, avg_spend, color=group)) +
  geom_line(size=1.1) +
  geom_vline(xintercept = intervention_date, linetype="dashed") +
  labs(title="Spend Jump at Intervention", y="Avg Daily Spend") +
  theme_fancy()

# Sales Response 
sales_data %>%
  mutate(group = if_else(treat==1,"Treated","Control")) %>%
  group_by(date, group) %>%
  summarize(avg_sales = mean(sales), .groups="drop") %>%
  ggplot(aes(date, avg_sales, color=group)) +
  geom_line(size=1.1) +
  geom_vline(xintercept = intervention_date, linetype="dashed") +
  labs(title="Sales Response to Spend Increase", y="Avg Daily Sales") +
  theme_fancy()

# Distribution of spend

ggplot(sales_data%>%
         filter(post==0), aes(x = as.factor(treat), y = sales, fill=factor(treat)))+
  geom_violin()+
  theme(legend.title=element_blank())+
  theme_fancy()

# Pre-period variance
leveneTest(sales ~ factor(treat), data = sales_data%>%
             filter(post==0), center = median)

leveneTest(spend ~ factor(treat), data = sales_data%>%
             filter(post==0), center = median)

# Non-normal data
fligner.test(sales ~ factor(treat), data = sales_data%>%
               filter(post==0))

ggplot(sales_data%>%
         filter(post==0), aes(x = factor(treat), y = sales, fill=factor(treat))) +
  geom_boxplot() +
  labs( title = "Pre-Period Sales Variance: Treated vs Control",
        x = "Group (0 = Control, 1 = Treated)",
        y = "Sales") +
  theme(legend.title=element_blank())+
  theme_fancy()

# Using helper function
diag_results <- pre_treatment_diagnostics(
  data = sales_data,
  outcome = "sales",
  treat = "treat",
  time = "date",
  post = "post",
  intervention_date = intervention_date,  
  do_event_study = TRUE                  
)

## Seasonality checks

# Day-of-week
sales_data %>%
  mutate(group = if_else(treat==1,"Treated","Control"),
         wday = weekdays(date)) %>%
  group_by(wday, group) %>%
  summarize(mean_sales = mean(sales)) %>%
  ggplot(aes(wday, mean_sales, fill = group)) +
  geom_col(position = position_dodge())+
  labs( title = "Day of Week Sales",
        x = "Group (0 = Control, 1 = Treated)",
        y = "Sales") +
  theme_fancy()

# Month
sales_data %>%
  mutate(group = if_else(treat==1,"Treated","Control"),
         month = format(date, "%m")) %>%
  group_by(month, group) %>%
  summarize(mean_sales = mean(sales)) %>%
  ggplot(aes(month, mean_sales, group=group, color=group)) + 
  geom_line(size=1.1)+
  geom_vline(xintercept = 3, linetype="dashed", color="red") +
  labs( title = "Monthly Sales",
        x = "Group (0 = Control, 1 = Treated)",
        y = "Sales") +
  theme_fancy()

## Autocorrelation
# Aggregate
acf(sales_data$sales)

pacf(sales_data$sales)

# By group
library(feasts) 

sales_ts <- sales_data %>%
  mutate(group = if_else(treat==1,"Treated","Control"))%>%
  group_by(date, group)%>%
  summarize(sales=sum(sales))%>%
  ungroup()%>%
  as_tsibble(key = group, index = date)   

# ACF
acf_df <- sales_ts %>%
  ACF(sales)

# PACF
pacf_df <- sales_ts %>%
  PACF(sales)

# Plot
acf_fig <- acf_df %>%
  autoplot() +
  labs(title = "ACF by Treatment Group")+
  theme_fancy()

pacf_fig <- pacf_df %>%
  autoplot()+
  labs(title = "PACF by Treatment Group")+
  theme_fancy()

acf_fig

pacf_fig

# Correlation Matrix
sales_data %>%
  select(sales, spend) %>%
  corrr::correlate(use = "pairwise.complete.obs")

# Prettier visualization of correlations
corr_mat <- sales_data %>%
  corrr::correlate() %>%     
  rearrange()    

# Same result
corr_mat

rplot(corr_mat,
      colours = c("red", "white", "blue"),
      print_cor = TRUE)+
  labs(title="Correlation Matrix")+
  theme_fancy()

ggplot(corr_mat %>% 
         stretch(na.rm = TRUE), aes(x = x, y = y, fill = r)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(
    low = "red",
    mid = "white",
    high = "blue",
    midpoint = 0
  ) +
  labs(title = "Correlation Matrix",
       x = "", y = "", fill = "Correlation") +
  theme_fancy() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank())

# Show pre-period parallel trends plot
diag_results$pre_trends_plot

# Inspect slope model
summary(diag_results$slope_model)

# Inspect Levene’s test object
diag_results$levene_test

## If the data is IID and satisfies SUTVA we can use difference-in-means and a T-test to measure pre/post diff

# Calculate the difference in means

diff_means <- sales_data %>%
  group_by(date, treat)%>%
  summarize(sales=sum(sales))%>%
  ungroup()%>%
  specify(sales ~ treat) %>%
  calculate("diff in means")

diff_means

# Generate a bootstrapped distribution of the difference in means based on our sample and calculate 
# the confidence interval:

boot_means <- sales_data %>%
  group_by(date, treat)%>%
  summarize(sales=sum(sales))%>%
  ungroup()%>%
  specify(sales ~ treat)%>%
  infer::generate(reps = 1000, type = "bootstrap") %>% 
  calculate("diff in means")

boostrapped_confint <- boot_means%>%
  get_confidence_interval()

boot_means %>% 
  visualize() + 
  shade_confidence_interval(boostrapped_confint,
                            color = "#8bc5ed", fill = "#85d9d2") +
  geom_vline(xintercept = diff_means$stat, size = 1, color = "#77002c") +
  labs(title = "Bootstrapped distribution of differences in means",
       x = "Treat - Untreated", y = "Count",
       subtitle = "Red line shows observed difference; shaded area shows 95% confidence interval") +
  theme_fancy()

# Step 2: Invent a world where δ is null
diffs_null <- sales_data %>%
  group_by(date, treat)%>%
  summarize(sales=sum(sales))%>%
  ungroup()%>%
  specify(sales ~ treat) %>% 
  infer::hypothesize(null = "independence") %>% 
  infer::generate(reps = 5000, type = "permute") %>% 
  calculate("diff in means")

# Step 3: Put actual observed δ in the null world and see if it fits
diffs_null %>% 
  visualize() + 
  geom_vline(xintercept = diff_means$stat, size = 1, color = "#77002c") +
  scale_y_continuous(labels = comma) +
  labs(x = "Simulated difference in sales (Treated - Untreated)", y = "Count",
       title = "Simulation-based null distribution of differences in means",
       subtitle = "Red line shows observed difference") +
  theme_fancy()

# Step 4: Calculate probability that observed δ could exist in null world
diffs_null %>% 
  infer::get_p_value(obs_stat = diff_means, direction = "both")

## Plotting KS stat
# A KS test is a comparison of the CDF of both datasets. We use this test because we want to know about
# more than just the mean, but how the tails of the dataset differ. This tells us if the datasets are
# structurally different across the entire range of observations. Additionally the KS test makes no
# assumptions about normality, equal variances, or linear relationships. It shows WHERE the groups differ.

# Grid of sales values
grid <- sales_data %>% 
  distinct(sales) %>% 
  arrange(sales)

# Define ECDF functions for each group
ecdf_treated_fun <- ecdf(sales_data$sales[sales_data$treat == 1])
ecdf_control_fun <- ecdf(sales_data$sales[sales_data$treat == 0])

# Wide ECDF data (one row per sales value, columns for each group's ECDF)
ecdf_wide <- grid %>%
  mutate(ecdf_treated = ecdf_treated_fun(sales),
         ecdf_control = ecdf_control_fun(sales),
         diff = abs(ecdf_treated - ecdf_control)
  )

# KS point (where the ECDF difference is maximal)
ks_point <- ecdf_wide %>% 
  slice_max(diff, n = 1)

# Long-format ECDF data frame used in ggplot
ecdf_wide %>%
  pivot_longer(cols  = c(ecdf_treated, ecdf_control),
               names_to  = "treat",
               values_to = "ecdf",
               names_prefix = "ecdf_" ) %>%
  mutate(
    treat = dplyr::case_when(
      treat == "treated" ~ "Treated",
      treat == "control" ~ "Control",
      TRUE ~ treat))%>%
  ggplot() +
  geom_step(
    aes(x = sales, y = ecdf, color = treat),
    linewidth = 1) +
  geom_segment(
    data = ks_point,
    aes(x = sales, xend = sales,
        y = ecdf_treated, yend = ecdf_control),
    color = "black", linetype = "dashed", linewidth=1) +
  annotate(
    "text",
    x = ks_point$sales,
    y = max(ks_point$ecdf_treated, ks_point$ecdf_control),
    label = paste0("KS = ", round(ks_point$diff, 3)),
    vjust = -0.5, hjust = -0.2,
    size = 5) +
  labs(
    title = "Kolmogorov–Smirnov Test Visualization",
    subtitle = "ECDF Comparison: Treated vs Control",
    x = "Sales",
    y = "ECDF",
    color = "Group") +
  theme_fancy() +
  scale_color_manual(
    values = c("Control" = "#1f78b4",
               "Treated" = "#e31a1c"))

## Event Study
# Estimate the treatment effect in treated vs control for every relative time period around the intervention.
# Pre coefs should avg ~0

# Create relative time (days to/from intervention)
sales_data <- sales_data %>%
  mutate(rel_time = as.integer(date - intervention_date))

es_nospend <- feols(
  sales ~ i(rel_time, treat, ref = -1) | city + date,
  data = sales_data, cluster = "city"
)

es_spend <- feols(
  sales ~ spend + i(rel_time, treat, ref = -1) | city + date,
  data = sales_data, cluster = "city"
)

event_df_spend <- tidy(es_spend, conf.int = TRUE) %>%
  filter(str_detect(term, "rel_time")) %>%
  mutate(rel_time = as.integer(str_extract(term, "-?\\d+"))) %>%
  arrange(rel_time)

ggplot(event_df_spend, aes(x = rel_time, y = estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point() +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) +
  labs(
    title = "Event-Study: Pre/Post Treatment Effects (controlling for spend)",
    x = "Days Relative to Campaign",
    y = "Estimated Effect on Sales"
  ) +
  theme_fancy()

# Regression, no spend
pre_data <- sales_data %>% filter(post == 0)

pre_trend_test <- lm(sales ~ as.numeric(date) * treat,
                     data = pre_data)

summary(pre_trend_test)

report::report(pre_trend_test)

# Summarize daily data
sales_control<-sales_data%>%
  filter(treat==0)%>%
  group_by(date)%>%
  summarize(control_sales=sum(sales),
            control_spend=sum(spend))%>%
  ungroup()

sales_treat<-sales_data%>%
  filter(treat==1)%>%
  group_by(date)%>%
  summarize(treat_sales=sum(sales),
            treat_spend=sum(spend))%>%
  ungroup()

# Define pre- and post-intervention periods 
dates <- seq.Date(as.Date(min(sales_treat$date)),
                  by=1, 
                  length.out = length(sales_treat$date))

pre_period <- as.Date(c(min(sales_treat$date), intervention_date-1))

post_period <- as.Date(c(intervention_date, max(sales_treat$date)))

sales_zoo<-zoo(cbind(sales_treat$treat_sales, 
                     sales_control$control_sales), order.by=dates)

# Causal for non smoothed data
sales_impact <- CausalImpact(sales_zoo, pre_period, post_period)

# Results summary
summary(sales_impact)

summary(sales_impact, "report")

# Plot
plot(sales_impact)+
  theme_fancy()+
  labs(x="Date",
       y="Effect",
       title="CausalImpact Estimate of Increase in Marketing Spend",
       subtitle=paste(round(sales_impact$summary$RelEffect,3)*100,
                      "% lift",
                      ",", 
                      "p=",
                      round(sales_impact$summary$p,3)))

## By group
control_cities<-sales_data%>%
  filter(treat==0)

treat_cities<-sales_data%>%
  filter(treat==1)

sales_data <- sales_data %>%
  mutate(date = ymd(date))  # ensure Date class

treat_df   <- sales_data %>% filter(treat == 1)
control_df <- sales_data %>% filter(treat == 0)

treated_cities <- unique(treat_df$city)
control_cities <- unique(control_df$city)

library(MarketMatching)

# Pre-period window for matching (adjust if you want)
pre_start <- "2024-01-01"
pre_end   <- "2024-05-31"

# Data that goes into best_matches = treated + control
library(MarketMatching)
library(dplyr)
library(lubridate)

sales_data <- sales_data %>%
  mutate(date = ymd(date))

treat_df   <- sales_data %>% filter(treat == 1)
control_df <- sales_data %>% filter(treat == 0)

synthetic_panel <- build_synthetic_control(
  data           = sales_data,
  matching_metric = "sales",
  date_col        = "date",
  city_col        = "city",
  treat_col       = "treat",
  match_window    = c("2024-01-01", "2024-05-31"),   # your pre-period
  n_matches       = 3
)

dplyr::glimpse(synthetic_panel)
# Expect: date, treat_sales, post, synth_sales, treated_city

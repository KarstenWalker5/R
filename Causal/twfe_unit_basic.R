
###### READ ME ######
# This script is a boilerplate Two-Way Fixed Effects (TWFE) implementation of Diff-in-Diff that should be used when
# comparing before/after effects of a marketing campaign. A normal T-test is not appropriate since the
# intervention is not randomized. A T-test assumes independent draws from populations with similar variance, 
# but we cannot assume that post-campaign variance is IDD. A T-test also ignores auto-correlation.

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
source("/Users/karstenwalker/Documents/R/Causal/causal_by_group_functions.R")

## Generate synthetic data for example
# Params
n_days <- 200
n_cities <- 10
start_date <- as.Date("2024-01-01")
intervention_date <- start_date + 100

treated_cities <- paste0("City_", 1:5)
control_cities <- paste0("City_", 6:10)

beta_spend <- 2        # spend → sales effect
spend_jump <- 8        # increase in spend for treated cities post-intervention
epsilon_sd <- 3        # noise level
beta_treat_direct <- 0 # optional direct treatment effect

sales_data <- 
  expand_grid(
    date = seq.Date(start_date, by = "day", length.out = n_days),
    city = paste0("City_", 1:n_cities)) %>%
  mutate(
    city  = factor(city),
    treat = if_else(city %in% treated_cities, 1L, 0L),
    post  = if_else(date >= intervention_date, 1L, 0L),
    day_index = as.integer(date - start_date),
    # City effects + common trend = parallel pre-trends
    city_effect  = rnorm(n_cities, sd = 4)[as.integer(city)],
    common_trend = 0.05 * day_index,
    # Raw spend + noise
    spend   = rnorm(n(), mean = 20, sd = 2),
    epsilon = rnorm(n(), sd = epsilon_sd)) %>%
  # balance spend in PRE-period across treated/control
  group_by(post, treat) %>%
  mutate(spend_grp_mean = mean(spend)) %>%
  ungroup() %>%
  group_by(post) %>%
  mutate(
    # overall mean spend in pre-period (post == 0)
    spend_overall_pre = if_else(post == 0, mean(spend[post == 0]), NA_real_)) %>%
  ungroup() %>%
  mutate(
    # adjust pre-period so treated & control have same mean spend
    spend_balanced = if_else(
      post == 0,
      spend - (spend_grp_mean - spend_overall_pre),
      spend ),
    # apply spend jump for treated cities post-intervention
    spend = spend_balanced + spend_jump * treat * post,
    # define sales as a function of spend
    sales = 50 +
      city_effect +
      common_trend +
      beta_spend * spend +
      beta_treat_direct * treat * post +
      epsilon) %>%
  select(date, city, treat, post, spend, sales)%>%
  mutate(treat=as.factor(treat))

# Inspect
glimpse(sales_data)

head(sales_data)

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
  cor(use = "pairwise.complete.obs")

# Show pre-period parallel trends plot
diag_results$pre_trends_plot

# Inspect slope model
summary(diag_results$slope_model)

# Inspect Levene’s test object
diag_results$levene_test

## Event Study
# Estimate the treatment effect in treated vs control for every relative time period around the intervention.
# Pre coefs should avg ~0

sales_data <- sales_data %>%
  mutate(
    rel_time = as.integer(date - intervention_date)
  )

event_study <- feols(
  sales ~ i(rel_time, treat, ref = -1) | city + date,
  data = sales_data,
  cluster = "city"
)

summary(event_study)

ggfixest::ggiplot(event_study,
      main = "Event-Study: Pre/Post Treatment Effects",
      xlab = "Days Relative to Campaign",
      ylab = "Estimated Effect on Sales")+
  theme_fancy()


# Regression, no spend
pre_data <- sales_data %>% filter(post == 0)

pre_trend_test <- lm(sales ~ as.numeric(date) * treat,
  data = pre_data)

summary(pre_trend_test)

# If spend is present
pre_data <- sales_data %>% filter(post == 0)

pre_spend_model<- lm(sales ~ treat * date, 
                     data = pre_data)

summary(pre_spend_model)

# TWFE
# All reported SEs and p-values are already cluster-robust at the city level.
did_twfe <- feols(
  sales ~ treat*post+spend | city + date,
  data = sales_data,
  cluster = "city"  # cluster SEs at unit level
)

summary(did_twfe)

# Cluster robust coef table
coeftest(did_twfe, vcov. = vcov_cl_city)

# treat:post coefficient = causal estimate of the campaign
# City FE remove static differences
# Date FE remove shocks common to all cities
# Clustering stabilizes inference

## City specific linear trends
# Compares treated vs. control relative to each city’s own linear drift
# Lowers statistical power
twfe_city_trend <- feols(
  sales ~ treat * post | city[trend] + date,
  data = sales_data,
  cluster = "city"
)
summary(twfe_city_trend)
twfe_city_trend <- feols(
  sales ~ treat * post | city[trend] + date,
  data = sales_data,
  cluster = "city"
)
summary(twfe_city_trend)

## Add DOW and WOM trends
sales_data <- sales_data %>%
  mutate(
    dow = weekdays(date),
    month = format(date, "%Y-%m")
  )

twfe_season <- feols(
  sales ~ treat * post | city + date + dow,
  data = sales_data,
  cluster = "city"
)

## Add weights
twfe_weighted <- feols(
  sales ~ treat * post + spend | city + date,
  data = sales_data,
  weights = ~ some_weight,
  cluster = "city"
)

## Introduce Hetergeneous Treatment Effects
city_spend_pre <- sales_data %>%
  filter(post == 0) %>%
  group_by(city) %>%
  summarize(mean_spend_pre = mean(spend), .groups = "drop")

sales_data <- sales_data %>%
  left_join(city_spend_pre, by = "city") %>%
  mutate(high_spend_city = mean_spend_pre > median(mean_spend_pre))

twfe_het <- feols(
  sales ~ treat * post * high_spend_city | city + date,
  data = sales_data,
  cluster = "city"
)
summary(twfe_het)

###### Interrupted Time Series ######
library(CausalImpact)
library(zoo)

# Insert logic comments

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

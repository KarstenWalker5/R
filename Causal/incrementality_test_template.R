
library(MarketMatching)
library(dplyr)
library(lubridate)
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
source("/Users/karstenwalker/Documents/GitHub/R/Causal/causal_pre_check.R")

# Load optional group by functions
source("/Users/karstenwalker/Documents/Github/R/Causal/causal_by_group_functions_dev.R")

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

# Visualize pre/post for each city
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

# Power Analysis
# Estimates how big an incremental effect you could reliably detect for a single treated city, given historical noise and post-period length.
# If your observed lift is smaller than mde_cum, a null result is not surprising.
mde_results <- run_mde_by_city(
  df = synthetic_panel,
  group_col = "treated_city",
  treat_col = treat_col,
  synth_col = synth_col,
  post_col = "post",
  alpha = 0.05,
  target_power = 0.80,
  two_sided = TRUE
)

# Create grouped DF
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

# Pre-period window for matching (adjust if you want)
pre_start <- "2024-01-01"
pre_end   <- "2024-05-31"

# Data that goes into best_matches = treated + control
sales_data <- sales_data %>%
  mutate(date = ymd(date))

treat_df   <- sales_data %>% filter(treat == 1)

control_df <- sales_data %>% filter(treat == 0)

# If you want to specify dates
synthetic_panel <- build_synthetic_control(
  data           = sales_data2,
  matching_metric = "sales",
  date_col        = "date",
  city_col        = "city",
  treat_col       = "treat",
  match_window    = c("2024-01-01", "2024-05-31"),   # your pre-period
  n_matches       = 3
)

# If you use an indicotor
synthetic_panel <- build_synthetic_control(
  data           = sales_data2,
  matching_metric = "sales",
  date_col        = "date",
  city_col        = "city",
  treat_col       = "treat",
  intervention_date_col     = 'intervention_date',   # your pre-period
  n_matches       = 3
)


dplyr::glimpse(synthetic_panel)

# Add intervention_date
synthetic_panel2<-synthetic_panel%>%
  group_by(treated_city)%>%
  mutate(intervention_date = min(date[post == 1]),
         end_date=max(date[post == 1]))%>%
  ungroup()

# If all we want are the top 10 controls
sales_data2<-sales_data%>%
  mutate(intervention_date = min(date[post == 1]),
         end_date=max(date[post == 1]))

top10 <- best_matches_all_cities(
  data = sales_data2,
  matching_metric = "sales",
  matches_per_city = 10,
  intervention_date_col = "intervention_date",
  pre_period_days = 90,       
  per_city_window = FALSE,    
  dtw_emphasis = 0.5
)


# Run CausalImpact by group, absolutely do not open the result dataframe
ci_nested <- synthetic_panel %>%
  group_by(treated_city) %>%
  arrange(date, .by_group = TRUE) %>%
  mutate(intervention_date = min(date[post == 1]),
         end_date=max(date[post == 1])) %>%
  nest() %>%
  mutate(
    causal_impact = map(data, ~ run_causal_impact(
      data = .x,
      intervention_date_col = "intervention_date",
      response_var = "treatment_sales",
      control_vars = c("synthetic_sales")
    ))
  ) %>%
  ungroup()

# Extract Results
ci_summary <- summarize_causal_results(ci_nested, 
                                       group_col = "treated_city")

# For safety you could do

ci_summary <- summarize_causal_results(synthetic_panel %>%
                                         group_by(treated_city) %>%
                                         arrange(date, .by_group = TRUE) %>%
                                         mutate(intervention_date = min(date[post == 1]),
                                                end_date=max(date[post == 1])) %>%
                                         nest() %>%
                                         mutate(
                                           causal_impact = map(data, ~ run_causal_impact(
                                             data = .x,
                                             intervention_date_col = "intervention_date",
                                             response_var = "treatment_sales",
                                             control_vars = c("synthetic_sales")
                                           ))
                                         ) %>%
                                         ungroup(), 
                                       group_col = "treated_city")

# View plots
plot_causal_impact_ggplots(ci_nested, pause = TRUE)

# One at a time
plot_causal_impact_ggplots(ci_nested%>% 
                             filter(group == "City_1"), pause = FALSE)

# Incrementality by city using DiD
inc_results <- synthetic_panel %>%
  group_by(treated_city) %>%
  group_modify(~ incrementality_by_city(
    .x,
    treat_col = "treatment_sales",
    synth_col = "synthetic_sales"
  )) %>%
  ungroup()

# Placebo Test
# This computes a distribution of placebo cumulative lifts by picking random contiguous windows of length n_post, 
# then returns a one-sided p-value (how often placebo ≥ observed).
# Tests whether the observed post-period lift for a city is unusual compared to fake “no-treatment” windows.

metric <- "sales"
treat_col <- paste0("treatment_", metric)
synth_col <- paste0("synthetic_", metric)

placebo_results <- run_placebo_tests_by_city(
  df = synthetic_panel,
  group_col = "treated_city",
  treat_col = treat_col,
  synth_col = synth_col,
  post_col = "post",
  n_placebos = 1000,
  seed = 1,
  one_sided = TRUE,
  keep_dist = FALSE   # set TRUE if you want distributions
)

## Final Reports

# Incrementality test
final_report <- lift_results %>%
  left_join(placebo_results, by = "treated_city") %>%
  left_join(mde_results, by = "treated_city")

# Causal Impact
ci_report <- ci_summary %>%
  select(
    treated_city,
    p_value_ci = p_value,
    rel_effect_ci = rel_effect,
    rel_effect_lower_ci = rel_effect_lower,
    rel_effect_upper_ci = rel_effect_upper)%>%
  left_join(
    placebo_results %>%
      select(
        treated_city,
        observed_cum_lift,
        placebo_p,
        ci_95_low,
        ci_95_high,
        z_score,
        snr), by = "treated_city")%>%
  left_join(
    mde_results %>%
      select(
        treated_city,
        n_post,
        mde_avg,
        mde_cum,
        pre_sd_gap),by = "treated_city")%>%
  mutate( sig_ci = !is.na(p_value_ci) & p_value_ci < 0.05,
    sig_placebo = !is.na(placebo_p) & placebo_p < 0.05,
    powered = !is.na(mde_cum) & abs(observed_cum_lift) >= mde_cum,
    strong_result = sig_ci & sig_placebo & powered,
    fragile_result = sig_ci & !sig_placebo,
    underpowered_null = !sig_ci & !powered,
    lift_vs_mde = observed_cum_lift / mde_cum )

# Prettier view
ci_report_summary<-ci_report %>%
  select(
    treated_city,
    observed_cum_lift,
    rel_effect_ci,
    p_value_ci,
    placebo_p,
    mde_cum,
    lift_vs_mde,
    sig_ci,
    sig_placebo,
    powered,
    strong_result,
    fragile_result,
    underpowered_null
  ) %>%
  arrange(placebo_p) 

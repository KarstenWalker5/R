
# This script uses the framework from this post and paper for diangnosing TWFE models: https://www.andrewheiss.com/blog/2021/08/25/twfe-diagnostics/
# Accompanying paper: https://arxiv.org/abs/2103.13229

library(tidyverse)     # For ggplot2, dplyr, and friends
library(haven)         # For reading Stata files
library(fixest)        # One way of doing fixed effects regression
library(estimatr)      # Another way of doing fixed effects regression
library(broom)         # For converting model objects to data frames
library(kableExtra)    # For pretty tables
library(modelsummary)  # For pretty regression tables
library(patchwork)     # For combining plots
library(scales)        # For nice formatting functions

# Custom theme 
# Get the font at https://fonts.google.com/specimen/Barlow+Semi+Condensed
theme_clean <- function() {
  theme_minimal(base_family = "Barlow Semi Condensed") +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold"),
          axis.title = element_text(family = "Barlow Semi Condensed Medium"),
          strip.text = element_text(family = "Barlow Semi Condensed",
                                    face = "bold", size = rel(1), hjust = 0),
          strip.background = element_rect(fill = "grey80", color = NA),
          plot.caption = element_text(hjust = 0))
}

# Simulate data, outcome represents a 0-100 range i.e. test scores
set.seed(1234)  

n_rows <- 500

fake_data <- tibble(treatment = rbeta(n_rows, shape1 = 3, shape2 = 7)) %>%
  mutate(treatment = round(treatment * 100)) %>%
  # Build the outcome based on some baseline level of outcome + a boost in
  # outcome that happens because of the treatment + some noise
  mutate(outcome_baseline = rnorm(n_rows, mean = 800, sd = 200),
         treatment_boost = 10 * treatment,
         outcome = outcome_baseline + treatment_boost + rnorm(n_rows, 200, 100)) %>%
  dplyr::select(outcome, treatment)

head(fake_data)

# Plot the data
hist_out <- ggplot(fake_data, aes(x = outcome)) +
  geom_histogram(binwidth = 100, color = "white",
                 boundary = 0, fill = "#FFA615") +
  scale_x_continuous(labels = dollar_format()) +
  labs(title = "Distribution of outcome",
       x = "Outcome", y = "Count") +
  coord_cartesian(ylim = c(0, 80)) +
  theme_clean() +
  theme(panel.grid.major.x = element_blank())

hist_trt <- ggplot(fake_data, aes(x = treatment)) +
  geom_histogram(binwidth = 5, color = "white",
                 boundary = 0, fill = "#A4BD0A") +
  labs(title = "Distribution of treatment",
       x = "Treatment", y = "Count") +
  coord_cartesian(ylim = c(0, 80)) +
  theme_clean() +
  theme(panel.grid.major.x = element_blank())

plot_trt_out <- ggplot(fake_data, aes(x = treatment, y = outcome)) +
  geom_point(size = 0.75) +
  geom_smooth(method = "lm", color = "#0599B0") +
  scale_y_continuous(labels = dollar_format()) +
  labs(title = "Effect of treatment on outcome",
       x = "Treatment", y = "Outcome") +
  theme_clean()

(hist_out + hist_trt) / plot_spacer() / plot_trt_out +
  plot_layout(heights = c(0.28, 0.02, 0.7))

# Fit a univariate regression to find the ATE:
model_effect <- lm(outcome ~ treatment, data = fake_data)

tidy(model_effect)

# Based on this model, a 1-unit increase in treatment causes a $9.55 increase in the outcome
# Rewrite the OLS estimate of based on the residuals of a model that predicts treatment,
# Since this is a univariate model, the residuals here are really just the deviations from the average value of the treatment
# We can confirm this by running an intercept-only model and comparing it 

# Intercept only model
trt_resid <- lm(treatment ~ 1, data = fake_data)

fake_data_resid <- fake_data %>%
  mutate(treatment_resid = residuals(trt_resid)) %>%
  mutate(mean_treat = mean(treatment),
         diff = treatment - mean_treat)

head(fake_data_resid)

# The treatment_resid and diff columns here are identical. Now that we have tghe residuqls, we can calculate beta

fake_data_resid %>%
  summarize(beta = sum(outcome * (treatment_resid / sum(treatment_resid^2))))

# beta is our original ATE

# Looking at residuals like this creates inherent weights for each observation. Essentially, betais a weighted sum of 
# the outcome variable, with the weights calculated with (treatment_resid / sum(treatment_resid^2))
# Calculate the weights for each observation:

fake_data_with_weights <- fake_data_resid %>%
  mutate(treatment_weight = treatment_resid / sum(treatment_resid^2))

head(fake_data_with_weights)

# given that this is overly perfect simulated data, the weights are really small. 
# In theory, the weights should sum up to 0. Let’s check that really quick:
  
fake_data_with_weights %>%
  summarize(total_weights = sum(treatment_weight))

# Visualize the distribution of weights too:
ggplot(fake_data_with_weights, aes(x = treatment_weight)) +
  geom_histogram(binwidth = 0.00005, color = "white",
                 boundary = 0, fill = "#F6C848") +
  geom_vline(xintercept = 0, color = "#FF2E00", linewidth = 1) +
  labs(x = "Residualized treatment weight", y = "Count") +
  scale_x_continuous(labels = comma_format()) +
  theme_clean() +
  theme(panel.grid.major.x = element_blank())

# In any univariate OLS regression of an outcome on a continuous measure of treatment intensity, 
# observations with below mean treatment intensity receive negative weight,and may be thought of as part of the comparison group.
# Since treatment here is continuous, there’s no clear division between treatment and control groups, so mathematically, 
# that threshold becomes the average value of treatment: observations where the treatment value is less than the average residual
# receive negative weight and can conceptually be considered part of the control/comparison group, while observations above
# average are in the treatment group. That’s apparent in the data too. Look at the first few rows of fake_data_with_weights above:
# observations where treatment_resid is negative have negative weights.

# To demonstrate issues with TWFE models, weights, and differential timing, Jakiela looks at the effect of eliminating 
# primary school fees on school enrollment in 15 African countries, based on data from the World Bank.

fpe_raw <- read_dta("/Users/karstenwalker/Downloads/WDI-FPE-data.dta")

# Remove rows where primary school enrollment is missing
fpe_primary <- fpe_raw %>%
  filter(!is.na(primary))

# Remove rows where secondary school enrollment is missing
fpe_secondary <- fpe_raw %>%
  filter(!is.na(secondary))

# Fixed Effects
model_lm <- lm(primary ~ treatment + country + factor(year),
               data = fpe_primary)

tidy(model_lm) %>%
  filter(!str_detect(term, "country"), !str_detect(term, "year"))

# Fixed Effects with cluster robust SE
# apply coeftest to our base model to obtain cluster robust SE
df_primary <- fpe_primary %>% 
  distinct(country) %>% 
  nrow()

model_lm_clustered <- lmtest::coeftest(model_lm,
                                       vcov = sandwich::vcovCL,
                                       cluster = ~country,
                                       df = df_primary - 1,
                                       # Keep original model so modelsummary shows R2
                                       save = TRUE)

tidy(model_lm_clustered) %>%
  filter(!str_detect(term, "country"), !str_detect(term, "year"))

# Use the estimatr to define special fixed effects that are automatically omitted from results tables 
# and returning robust and optionally clustered standard errors:
model_lm_robust <- lm_robust(primary ~ treatment,
                             fixed_effects = ~ country + year,
                             data = fpe_primary,
                             clusters = country, se_type = "stata")

tidy(model_lm_robust)

glance(model_lm_robust)

# Use fixest package to do the same
model_feols <- feols(
  primary ~ treatment | country + year,
  data    = fpe_primary,
  cluster = ~ country
)

tidy(model_feols)

glance(model_feols)

# We can show these all in a side-by-side table using the modelsummary package. 
# Conveniently, modelsummary() also lets you adjust standard errors on the fly with the vcov argument, 
# so we could theoretically handle all the clustering here instead of inside lmtest::coeftest(), lm_robust(), or feols(). 
# But since we already specified the clusters above, we’ll just use those.
# Look at secondary schools too

model_lm_sec <- lm(secondary ~ treatment + country + factor(year),
                   data = fpe_secondary)

df_secondary <- fpe_secondary %>% distinct(country) %>% nrow()

model_lm_clustered_sec <- lmtest::coeftest(model_lm_sec,
                                           vcov = sandwich::vcovCL,
                                           cluster = ~country,
                                           df = df_secondary - 1,
                                           save = TRUE)

model_lm_robust_sec <- lm_robust(secondary ~ treatment,
                                 fixed_effects = ~ country + year,
                                 data = fpe_secondary,
                                 clusters = country, se_type = "stata")

model_feols_sec <- feols(
  primary ~ treatment | country + year,
  data    = fpe_secondary,
  cluster = ~ country
)

# Define the goodness-of-fit stats to include
gof_stuff <- tribble(
  ~raw, ~clean, ~fmt,
  "nobs", "N", 0,
  "r.squared", "R²", 3
)

# Define extra rows at the end of the table
extra_rows <- tribble(
  ~term, ~a, ~b, ~c, ~d, ~e, ~f,
  "Country fixed effects", "•", "•", "•", "•", "•", "•",
  "Year fixed effects", "•", "•", "•", "•", "•", "•"
)

modelsummary(list("<code>lm()</code>" = model_lm_clustered,
                  "<code>lm_robust()</code>" = model_lm_robust,
                  "<code>feols()</code>" = model_feols,
                  "<code>lm()</code>" = model_lm_clustered_sec,
                  "<code>lm_robust()</code>" = model_lm_robust_sec,
                  "<code>feols()</code>" = model_feols_sec),
             coef_rename = c("treatment" = "Treatment"),
             estimate = "{estimate}",
             statistic = c("s.e. = {std.error}", "p = {p.value}"),
             coef_omit = "^country|^factor|Intercept",
             gof_map = gof_stuff,
             add_rows = extra_rows,
             escape = FALSE, output = "kableExtra") %>%
  add_header_above(c(" " = 1, "Primary enrollment" = 3, "Secondary enrollment" = 3)) %>%
  kable_styling(htmltable_class = "table table-sm")

# All these models show the same result: eliminating primary school fees caused primary school 
# enrollment to increase by 20.4 percentage points.

# TWFE with residualized treatment weights
# Countries have staggered start times

plot_fpe_start <- fpe_raw %>%
  filter(year == fpe_year) %>%
  arrange(fpe_year, country) %>%
  mutate(country = fct_inorder(country))

ggplot(plot_fpe_start, aes(y = fct_rev(country))) +
  geom_segment(aes(x = fpe_year, xend = 2015, yend = country),
               size = 3, color = "#D6EB52") +
  geom_text(aes(x = fpe_year + 0.2), label = "▶", family = "Arial Unicode MS",
            size = 8, color = "#A4BD0A") +
  labs(x = "Year free primary education implemented", y = NULL) +
  theme_clean()

# In the simple model earlier, the treatment residuals were just the residuals from lm(treatment ~ 1). 
# In TWFE situations, though, treatment residuals need to account for both country and year fixed effects, 
# or lm(treatment ~ country + year)

trt_resid_primary <- lm(treatment ~ country + factor(year), data = fpe_primary)

trt_resid_secondary <- lm(treatment ~ country + factor(year), data = fpe_secondary)

fpe_primary_weights <- fpe_primary %>%
  mutate(treatment_resid = residuals(trt_resid_primary)) %>%
  mutate(treatment_weight = treatment_resid / sum(treatment_resid^2))

fpe_secondary_weights <- fpe_secondary %>%
  mutate(treatment_resid = residuals(trt_resid_secondary)) %>%
  mutate(treatment_weight = treatment_resid / sum(treatment_resid^2))

fpe_primary_weights %>%
  summarize(twfe_beta_primary = sum(primary * treatment_weight))

fpe_secondary_weights %>%
  summarize(twfe_beta_secondary = sum(secondary * treatment_weight))

# With TWFE, some observations’ weights switch directions. There are systematic reasons for this. A
# Negative weights in treated observations are more likely in (1) early adopter countries, since the country-level 
# treatment mean is high, and (2) later years, since the year-level treatment mean is higher.

# 2 diagnostic questions
#   1. Do any treated units get negative weight when calculating beta of TWFE? Check this by looking at the weights
#   2. Can we reject the hypothesis that the treatment effects are homogenous? Check this by looking at the relationship between 
#       Y and the residusals. The slope shouldn’t be different.

# Investigate patterns in the residualized weights. Some of the treated country-year observations have negative weight:
# Total treated in primary data
n_treated_primary <- fpe_primary_weights %>%
  filter(treatment == 1) %>%
  nrow()

# Total treated in secondary data
n_treated_secondary <- fpe_secondary_weights %>%
  filter(treatment == 1) %>%
  nrow()

# Negatively weighted treated observations in the primary data
n_treated_negative_primary <- fpe_primary_weights %>%
  filter(treatment_weight < 0 & treatment == 1) %>%
  nrow()

n_treated_negative_primary

n_treated_negative_primary / n_treated_primary

# Negatively weighted treated observations in the secondary data
n_treated_negative_secondary <- fpe_secondary_weights %>%
  filter(treatment_weight < 0 & treatment == 1) %>%
  nrow()

n_treated_negative_secondary

# Roughly a quarter of the treated observations in both the primary and secondary enrollment 
# data are negatively weighted. That might be bad.

# Visualize distribution of weights
plot_weights_primary <- fpe_primary_weights %>%
  mutate(treatment_fct = factor(treatment, labels = c("Untreated", "Treated"), ordered = TRUE)) %>%
  mutate(oh_no = (treatment == 1 & treatment_weight < 0) | (treatment == 0 & treatment_weight > 0))

plot_weights_secondary <- fpe_secondary_weights %>%
  mutate(treatment_fct = factor(treatment, labels = c("Untreated", "Treated"), ordered = TRUE)) %>%
  mutate(oh_no = (treatment == 1 & treatment_weight < 0) | (treatment == 0 & treatment_weight > 0))

hist_primary <- ggplot(plot_weights_primary,
                       aes(x = treatment_weight, fill = treatment_fct, alpha = oh_no)) +
  geom_histogram(binwidth = 0.002, color = "white", boundary = 0,
                 position = position_identity()) +
  geom_vline(xintercept = 0, color = "#FF4136", size = 0.5) +
  scale_alpha_manual(values = c(0.6, 1)) +
  scale_fill_viridis_d(option = "rocket", begin = 0.2, end = 0.8) +
  labs(x = "Residualized treatment", y = "Count",
       title = "Primary school enrollment") +
  guides(fill = "none", alpha = "none") +
  facet_wrap(vars(treatment_fct), ncol = 1) +
  theme_clean()

hist_secondary <- ggplot(plot_weights_secondary,
                         aes(x = treatment_weight, fill = treatment_fct, alpha = oh_no)) +
  geom_histogram(binwidth = 0.002, color = "white", boundary = 0,
                 position = position_identity()) +
  geom_vline(xintercept = 0, color = "#FF4136", size = 0.5) +
  scale_alpha_manual(values = c(0.6, 1)) +
  scale_fill_viridis_d(option = "rocket", begin = 0.2, end = 0.8) +
  labs(x = "Residualized treatment", y = "Count",
       title = "Secondary school enrollment") +
  guides(fill = "none", alpha = "none") +
  facet_wrap(vars(treatment_fct), ncol = 1) +
  theme_clean()

hist_primary + hist_secondary

n_treated_negative_secondary / n_treated_secondary

# Which treated country/years are getting negative weights?
plot_waffle_primary <- fpe_primary_weights %>%
  mutate(weight_fill = case_when(
    treatment_resid < 0 & treatment ~ "Treatment observations, negative weight",
    treatment_resid > 0 & treatment ~ "Treatment observations, positive weight",
    !treatment ~ "Comparison observations"
  )) %>%
  arrange(desc(fpe_year), desc(country)) %>%
  mutate(country = fct_inorder(country)) %>%
  mutate(weight_fill = fct_inorder(weight_fill))

plot_waffle_secondary <- fpe_secondary_weights %>%
  mutate(weight_fill = case_when(
    treatment_resid < 0 & treatment ~ "Treatment observations, negative weight",
    treatment_resid > 0 & treatment ~ "Treatment observations, positive weight",
    !treatment ~ "Comparison observations"
  )) %>%
  arrange(desc(fpe_year), desc(country)) %>%
  mutate(country = fct_inorder(country)) %>%
  mutate(weight_fill = fct_inorder(weight_fill))

waffle_primary <- ggplot(plot_waffle_primary,
                         aes(x = year, y = country, fill = weight_fill)) +
  geom_tile(color = "white", size = 0.5) +
  scale_fill_manual(values = c("grey50", "#0074D9", "#FF4136"),
                    guide = guide_legend(reverse = TRUE),
                    name = NULL) +
  scale_x_continuous(expand = expansion(add = 0.5),
                     breaks = seq(1980, 2015, 5)) +
  labs(x = NULL, y = NULL, title = "Primary school enrollment") +
  coord_equal() +
  theme_clean() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        legend.key.size = unit(0.8, "lines"))

waffle_secondary <- ggplot(plot_waffle_secondary,
                           aes(x = year, y = country, fill = weight_fill)) +
  geom_tile(color = "white", size = 0.5) +
  scale_fill_manual(values = c("grey50", "#0074D9", "#FF4136"),
                    guide = guide_legend(reverse = TRUE),
                    name = NULL) +
  scale_x_continuous(expand = expansion(add = 0.5),
                     breaks = seq(1980, 2015, 5)) +
  labs(x = NULL, y = NULL, title = "Secondary school enrollment") +
  coord_equal() +
  theme_clean() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        legend.key.size = unit(0.8, "lines"))

waffle_primary / waffle_secondary

# We thus have issues with negative weighting here for early adopting countries and later country-years. 
# However, as long as the treatment effects are homogenous—or that the elimination of school fees has the 
# same effect across time and country—this negative weighting isn’t an issue. If treatment effects are 
# heterogenous—especially if the effects change over time within treated countries—the TWFE estimate will be severely biased.

# Second diagnostic tests this assumption of homogeneity based on the mathematical relationship between the residuals 
# of the outcome variableand the residuals of the treatment variable . Essentially, if there’s no difference in slopes 
# across treated and untreated observations in a regression there’s no evidence for heterogeneity and all is well.

# Build models for the residualized outcomes
out_resid_primary <- lm(primary ~ country + factor(year), data = fpe_primary)

out_resid_secondary <- lm(secondary ~ country + factor(year), data = fpe_secondary)

# Add residuals to data with weights
fpe_primary_weights <- fpe_primary_weights %>%
  mutate(out_resid = residuals(out_resid_primary))

fpe_secondary_weights <- fpe_secondary_weights %>%
  mutate(out_resid = residuals(out_resid_secondary))

# Plot the residuals
plot_out_trt_primary <- ggplot(fpe_primary_weights,
                               aes(x = treatment_resid, y = out_resid, color = factor(treatment))) +
  geom_point(size = 0.75, alpha = 0.5) +
  geom_smooth(aes(linetype = "Loess"), method = "loess", size = 1, se = FALSE, alpha = 0.5) +
  geom_smooth(aes(linetype = "OLS"), method = "lm", se = FALSE) +
  scale_color_viridis_d(option = "rocket", begin = 0.2, end = 0.8,
                        labels = c("Untreated", "Treated")) +
  scale_linetype_manual(values = c("OLS" = "solid", "Loess" = "21"),
                        guide = guide_legend(override.aes = list(color = "grey30"))) +
  labs(x = "Residualized treatment", y = "Residualized outcome",
       title = "Primary school enrollment", color = NULL, linetype = NULL) +
  theme_clean() +
  theme(legend.position = "bottom")

plot_out_trt_secondary <- ggplot(fpe_secondary_weights,
                                 aes(x = treatment_resid, y = out_resid, color = factor(treatment))) +
  geom_point(size = 0.75, alpha = 0.5) +
  geom_smooth(aes(linetype = "Loess"), method = "loess", size = 1, se = FALSE, alpha = 0.5) +
  geom_smooth(aes(linetype = "OLS"), method = "lm", se = FALSE) +
  scale_color_viridis_d(option = "rocket", begin = 0.2, end = 0.8,
                        labels = c("Untreated", "Treated")) +
  scale_linetype_manual(values = c("OLS" = "solid", "Loess" = "21"),
                        guide = guide_legend(override.aes = list(color = "grey30"))) +
  labs(x = "Residualized treatment", y = "Residualized outcome",
       title = "Secondary school enrollment", color = NULL, linetype = NULL) +
  theme_clean() +
  theme(legend.position = "bottom")

plot_out_trt_primary | plot_out_trt_secondary

# For primary enrollment, the lines for treated and untreated observations look like they have similar slopes. 
# When using Loess lines, the two groups don’t align perfectly, and the relationship between residualized treatment 
# and residualized outcome doesn’t look perfectly linear. When secondary enrollment is the outcome, the lines for 
# the treated and untreated groups seem to switch directions—the slope is negative for untreated observations but 
# positive for treated. This is a bad sign for the homogeneity assumption.

# We can statistically test if the slopes in two groups are the same or not by using an interaction term in a regression model:

check_slopes_primary <- lm(out_resid ~ treatment_resid * factor(treatment),
                           data = fpe_primary_weights)

check_slopes_secondary <- lm(out_resid ~ treatment_resid * factor(treatment),
                             data = fpe_secondary_weights)

modelsummary(list("Primary enrollment" = check_slopes_primary,
                  "Secondary enrollment" = check_slopes_secondary),
             coef_rename = c("treatment_resid" = "Residualized treatment",
                             "factor(treatment)1" = "Treatment group",
                             "treatment_resid:factor(treatment)1" = "Treatment group × residualized treatment",
                             "(Intercept)" = "Intercept"),
             estimate = "{estimate}",
             statistic = c("s.e. = {std.error}", "p = {p.value}"),
             gof_map = tribble(
               ~raw, ~clean, ~fmt,
               "nobs", "N", 0,
               "r.squared", "R²", 3
             ),
             escape = FALSE, output = "kableExtra") %>%
  kable_styling(htmltable_class = "table table-sm")

# The coefficient for “Treatment group × residualized treatment” shows the change in slope between the two groups. 
# For primary enrollment, being in the treatment group reduces the slope by 7.8 (so from 23.761 to 15.961), 
# but the standard errors around that change are huge and the p-value is not significant (p = 0.199). For secondary enrollment, 
# though, being in the treatment group increases the slope by 5.3 (from -2.9 to positive 2.4), and that change is statistically 
# significant (p = 0.009).

# So for primary enrollment, there’s not enough evidence to reject the hypothesis that the slopes are the same 
# (or the hypothesis that the effect is homogenous), so we’ll treat the effect as homogenous. 
# That means we don’t need to worry so much about the negatively-weighted treated country-years. 
# For secondary enrollment, though, there’s pretty strong evidence that the slopes aren’t the same, which means the treatment 
# effect is likely heterogenous, which also means that the treated country-years with negative weights are biasing the results 
# substantially and we should thus be worried.

###### Robustness Checks ######
# As long as we assume that treatment effects are homogenous, we can safely drop some observations from the data and still 
# find the same result. We can also strategically drop observations to check if negative treatment weights influence the results. 
# Three different robustness checks that all involve dropping different categories of observations:
#   1. Exclude later years
#   2. Limit how many post-treatment years are kept
#   3. Exclude individual countries

## Exclude end years
different_max_years_primary <- tibble(end_year = 2000:2015) %>%
  # Nest a filtered dataset in a cell for each year
  mutate(data = map(end_year, ~filter(fpe_primary, year <= .))) %>%
  # Calculate the TWFE estimate for each dataset
  mutate(model = map(data, ~lm_robust(primary ~ treatment,
                                      fixed_effects = ~ country + year,
                                      data = .,
                                      clusters = country, se_type = "stata"))) %>%
  # Extract a data frame of the results
  mutate(tidied = map(model, ~tidy(., conf.int = TRUE))) %>%
  # Calculate residuals and treatment weights
  mutate(model_resid = map(data, ~lm(treatment ~ country + factor(year), data = .)),
         treatment_resid = map(model_resid, ~residuals(.)),
         treatment_weight = map(treatment_resid, ~ . / sum(.^2)))

# Calculate how many treated observations have negative weights
prop_negative_primary <- different_max_years_primary %>%
  unnest(c(data, treatment_resid, treatment_weight)) %>%
  group_by(end_year) %>%
  summarize(n_treated_negative_weight = sum(treatment_weight < 0 & treatment == 1),
            n_treated = sum(treatment == 1),
            prop_treated_negative_weight = n_treated_negative_weight / n_treated) %>%
  mutate(prop_nice = percent(prop_treated_negative_weight, accuracy = 1))

# Extract tidied results for plotting
coefs_to_plot_primary <- different_max_years_primary %>%
  unnest(tidied)

different_max_years_secondary <- tibble(end_year = 2000:2015) %>%
  mutate(data = map(end_year, ~filter(fpe_secondary, year <= .))) %>%
  mutate(model = map(data, ~lm_robust(secondary ~ treatment,
                                      fixed_effects = ~ country + year,
                                      data = .,
                                      clusters = country, se_type = "stata"))) %>%
  mutate(tidied = map(model, ~tidy(., conf.int = TRUE))) %>%
  mutate(model_resid = map(data, ~lm(treatment ~ country + factor(year), data = .)),
         treatment_resid = map(model_resid, ~residuals(.)),
         treatment_weight = map(treatment_resid, ~ . / sum(.^2)))

prop_negative_secondary <- different_max_years_secondary %>%
  unnest(c(data, treatment_resid, treatment_weight)) %>%
  group_by(end_year) %>%
  summarize(n_treated_negative_weight = sum(treatment_weight < 0 & treatment == 1),
            n_treated = sum(treatment == 1),
            prop_treated_negative_weight = n_treated_negative_weight / n_treated) %>%
  mutate(prop_nice = percent(prop_treated_negative_weight, accuracy = 1))

coefs_to_plot_secondary <- different_max_years_secondary %>%
  unnest(tidied)

# Coefficient plot
ate_primary <- ggplot(coefs_to_plot_primary, aes(x = end_year, y = estimate)) +
  geom_hline(yintercept = 0, color = "#FF2E00", size = 1) +
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high),
                  color = "#0599B0", size = 1, fatten = 2) +
  labs(x = NULL, y = "TWFE-based treatment effect",
       title = "Primary school enrollment") +
  theme_clean() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank())

ate_secondary <- ggplot(coefs_to_plot_secondary, aes(x = end_year, y = estimate)) +
  geom_hline(yintercept = 0, color = "#FF2E00", size = 1) +
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high),
                  color = "#0599B0", size = 1, fatten = 2) +
  labs(x = NULL, y = "TWFE-based treatment effect",
       title = "Secondary school enrollment") +
  theme_clean() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank())

# Bar plot
prop_neg_primary <- ggplot(filter(prop_negative_primary,
                                  prop_treated_negative_weight > 0),
                           aes(x = end_year, y = prop_treated_negative_weight)) +
  geom_col(fill = "#FF4136") +
  geom_text(aes(label = prop_nice),
            nudge_y = 0.02, size = 2.5,
            family = "Barlow Semi Condensed Bold", color = "#FF4136") +
  coord_cartesian(xlim = c(2000, 2015)) +
  labs(x = "Last year included in data", y = NULL,
       subtitle = "% of treated observations with negative weight") +
  theme_clean() +
  theme(panel.grid.major.y = element_blank(),
        axis.text.y = element_blank())

prop_neg_secondary <- ggplot(filter(prop_negative_secondary,
                                    prop_treated_negative_weight > 0),
                             aes(x = end_year, y = prop_treated_negative_weight)) +
  geom_col(fill = "#FF4136") +
  geom_text(aes(label = prop_nice),
            nudge_y = 0.02, size = 2.5,
            family = "Barlow Semi Condensed Bold", color = "#FF4136") +
  coord_cartesian(xlim = c(2000, 2015)) +
  labs(x = "Last year included in data", y = NULL,
       subtitle = "% of treated observations with negative weight") +
  theme_clean() +
  theme(panel.grid.major.y = element_blank(),
        axis.text.y = element_blank())

((ate_primary / prop_neg_primary) + plot_layout(heights = c(0.8, 0.2))) |
  ((ate_secondary / prop_neg_secondary) + plot_layout(heights = c(0.8, 0.2)))

## Drop post-treatment years
# Dropping rows based on end years like this is neat, but it doesn’t take into account the differential treatment timing.

fpe_primary_years_since <- fpe_primary %>%
  mutate(years_after = year - fpe_year)

fpe_secondary_years_since <- fpe_secondary %>%
  mutate(years_after = year - fpe_year)

different_years_after_primary <- tibble(post_trt_years = 2:22) %>%
  # Nest a filtered dataset in a cell for each year
  mutate(data = map(post_trt_years, ~filter(fpe_primary_years_since, years_after <= .))) %>%
  # Calculate the TWFE estimate for each dataset
  mutate(model = map(data, ~lm_robust(primary ~ treatment,
                                      fixed_effects = ~ country + year,
                                      data = .,
                                      clusters = country, se_type = "stata"))) %>%
  # Extract a data frame of the results
  mutate(tidied = map(model, ~tidy(., conf.int = TRUE))) %>%
  # Calculate residuals and treatment weights
  mutate(model_resid = map(data, ~lm(treatment ~ country + factor(year), data = .)),
         treatment_resid = map(model_resid, ~residuals(.)),
         treatment_weight = map(treatment_resid, ~ . / sum(.^2)))

# Calculate how many treated observations have negative weights
prop_negative_primary <- different_years_after_primary %>%
  unnest(c(data, treatment_resid, treatment_weight)) %>%
  group_by(post_trt_years) %>%
  summarize(n_treated_negative_weight = sum(treatment_weight < 0 & treatment == 1),
            n_treated = sum(treatment == 1),
            prop_treated_negative_weight = n_treated_negative_weight / n_treated) %>%
  mutate(prop_nice = percent(prop_treated_negative_weight, accuracy = 1))

# Extract tidied results for plotting
coefs_to_plot_primary <- different_years_after_primary %>%
  unnest(tidied)

different_years_after_secondary <- tibble(post_trt_years = 2:22) %>%
  mutate(data = map(post_trt_years, ~filter(fpe_secondary_years_since, years_after <= .))) %>%
  mutate(model = map(data, ~lm_robust(secondary ~ treatment,
                                      fixed_effects = ~ country + year,
                                      data = .,
                                      clusters = country, se_type = "stata"))) %>%
  mutate(tidied = map(model, ~tidy(., conf.int = TRUE))) %>%
  mutate(model_resid = map(data, ~lm(treatment ~ country + factor(year), data = .)),
         treatment_resid = map(model_resid, ~residuals(.)),
         treatment_weight = map(treatment_resid, ~ . / sum(.^2)))

prop_negative_secondary <- different_years_after_secondary %>%
  unnest(c(data, treatment_resid, treatment_weight)) %>%
  group_by(post_trt_years) %>%
  summarize(n_treated_negative_weight = sum(treatment_weight < 0 & treatment == 1),
            n_treated = sum(treatment == 1),
            prop_treated_negative_weight = n_treated_negative_weight / n_treated) %>%
  mutate(prop_nice = percent(prop_treated_negative_weight, accuracy = 1))

coefs_to_plot_secondary <- different_years_after_secondary %>%
  unnest(tidied)

# Coefficient plot
ate_primary <- ggplot(coefs_to_plot_primary, aes(x = post_trt_years, y = estimate)) +
  geom_hline(yintercept = 0, color = "#FF2E00", size = 1) +
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high),
                  color = "#0599B0", size = 1, fatten = 2) +
  labs(x = NULL, y = "TWFE-based treatment effect",
       title = "Primary school enrollment") +
  theme_clean() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank())

ate_secondary <- ggplot(coefs_to_plot_secondary, aes(x = post_trt_years, y = estimate)) +
  geom_hline(yintercept = 0, color = "#FF2E00", size = 1) +
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high),
                  color = "#0599B0", size = 1, fatten = 2) +
  labs(x = NULL, y = "TWFE-based treatment effect",
       title = "Secondary school enrollment") +
  theme_clean() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank())

# Bar plot
prop_neg_primary <- ggplot(filter(prop_negative_primary,
                                  prop_treated_negative_weight > 0),
                           aes(x = post_trt_years, y = prop_treated_negative_weight)) +
  geom_col(fill = "#FF4136") +
  geom_text(aes(label = prop_nice),
            nudge_y = 0.02, size = 2.5,
            family = "Barlow Semi Condensed Bold", color = "#FF4136") +
  coord_cartesian(xlim = c(2, 22)) +
  labs(x = "Post-treatment years included", y = NULL,
       subtitle = "% of treated observations with negative weight") +
  theme_clean() +
  theme(panel.grid.major.y = element_blank(),
        axis.text.y = element_blank())

prop_neg_secondary <- ggplot(filter(prop_negative_secondary,
                                    prop_treated_negative_weight > 0),
                             aes(x = post_trt_years, y = prop_treated_negative_weight)) +
  geom_col(fill = "#FF4136") +
  geom_text(aes(label = prop_nice),
            nudge_y = 0.02, size = 2.5,
            family = "Barlow Semi Condensed Bold", color = "#FF4136") +
  coord_cartesian(xlim = c(2, 22)) +
  labs(x = "Last year included in data", y = NULL,
       subtitle = "% of treated observations with negative weight") +
  theme_clean() +
  theme(panel.grid.major.y = element_blank(),
        axis.text.y = element_blank())

((ate_primary / prop_neg_primary) + plot_layout(heights = c(0.8, 0.2))) |
  ((ate_secondary / prop_neg_secondary) + plot_layout(heights = c(0.8, 0.2)))

## Exclude individual countries
fpe_omit_countries_primary <- fpe_primary %>%
  arrange(fpe_year, country) %>%
  mutate(country_start = paste0(country, " (", fpe_year, ")"),
         country_start = fct_inorder(country_start)) %>%
  distinct(country, country_start) %>%
  mutate(data = map(country, ~filter(fpe_primary, country != .))) %>%
  mutate(model = map(data, ~lm_robust(primary ~ treatment,
                                      fixed_effects = ~ country + year,
                                      data = .,
                                      clusters = country, se_type = "stata"))) %>%
  mutate(tidied = map(model, ~tidy(., conf.int = TRUE)))

coefs_to_plot_primary <- fpe_omit_countries_primary %>%
  unnest(tidied)

fpe_omit_countries_secondary <- fpe_secondary %>%
  arrange(fpe_year, country) %>%
  mutate(country_start = paste0(country, " (", fpe_year, ")"),
         country_start = fct_inorder(country_start)) %>%
  distinct(country, country_start)  %>%
  mutate(data = map(country, ~filter(fpe_secondary, country != .))) %>%
  mutate(model = map(data, ~lm_robust(secondary ~ treatment,
                                      fixed_effects = ~ country + year,
                                      data = .,
                                      clusters = country, se_type = "stata"))) %>%
  mutate(tidied = map(model, ~tidy(., conf.int = TRUE)))

coefs_to_plot_secondary <- fpe_omit_countries_secondary %>%
  unnest(tidied)

ate_primary <- ggplot(coefs_to_plot_primary, aes(x = estimate, y = fct_rev(country_start))) +
  geom_vline(xintercept = 0, color = "#FF2E00", size = 1) +
  geom_pointrange(aes(xmin = conf.low, xmax = conf.high),
                  color = "#0599B0", size = 1, fatten = 2) +
  labs(x = "TWFE-based treatment effect", y = NULL,
       title = "Primary school enrollment",
       subtitle = "Each point represents the ATE when omitting the specified country",
       caption = "Year that fees were eliminated is shown in parentheses") +
  coord_cartesian(xlim = c(-10, 50)) +
  theme_clean() +
  theme(panel.grid.major.y = element_blank())

ate_secondary <- ggplot(coefs_to_plot_secondary, aes(x = estimate, y = fct_rev(country_start))) +
  geom_vline(xintercept = 0, color = "#FF2E00", size = 1) +
  geom_pointrange(aes(xmin = conf.low, xmax = conf.high),
                  color = "#0599B0", size = 1, fatten = 2) +
  labs(x = "TWFE-based treatment effect", y = NULL,
       title = "Secondary school enrollment") +
  coord_cartesian(xlim = c(-10, 50)) +
  theme_clean() +
  theme(panel.grid.major.y = element_blank())

ate_primary | ate_secondary

###### Summary ######
# 1. We can think of OLS as a weighted sum of outcome values - these weights are typically negative for untreated 
# observations and positive for treated observations
# 2. In TWFE situations, these weights take country and year into account
# 3.Because of mathy reasons, treated observations can actually get negative weights—especially countries that are early adopters, 
# and country-year observations in later years
# 4. We can see if these negative weights on treated observations are an issue by using a couple simple diagnostics: 
# see how many and which treated observations have negative weights, and see if the treatment effect is homogenous across treated countries
# 5. We can look at which treated observations have negative weights by making some plots and exploring the data
# 6. We can check for the homogeniety of the treatment effect by running a regression that uses the residualized 
# treatment and treatment status to explain the residualized outcome. If the slopes for treated and comparison 
# observations are indistinguishable, we can safely assume treatment homogeneity.
# 7. Finally, we can check how robust the TWFE estimate is to negative treatment weights and the assumption of homogeneity 
# by dropping specific types of observations and checking the ATE across these modified datasets
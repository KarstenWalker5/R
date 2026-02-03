# Based on this post: https://www.andrewheiss.com/blog/2024/03/21/demystifying-ate-att-atu/#:~:text=Contents,like%20the%20average%20treatment%20effect.

library(tidyverse)
library(ggtext)
library(ggdag)
library(dagitty)
library(gt)
library(broom)
library(marginaleffects)
library(WeightIt)

# Define a nice color palette from {MoMAColors}
# https://github.com/BlakeRMills/MoMAColors
clrs <- MoMAColors::moma.colors("ustwo")

# Download Mulish from https://fonts.google.com/specimen/Mulish
theme_nice <- function() {
  theme_minimal(base_family = "Mulish") +
    theme(
      panel.grid.minor = element_blank(),
      plot.background = element_rect(fill = "white", color = NA),
      plot.title = element_text(face = "bold"),
      axis.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold"),
      strip.background = element_rect(fill = "grey80", color = NA),
      legend.title = element_text(face = "bold")
    )
}

theme_set(theme_nice())

update_geom_defaults("text", list(family = "Mulish", fontface = "plain"))
update_geom_defaults("label", list(family = "Mulish", fontface = "plain"))
update_geom_defaults(ggdag:::GeomDagText, list(family = "Mulish", fontface = "plain"))
update_geom_defaults(ggtext::GeomRichText, list(family = "Mulish", fontface = "plain"))
# Generate synth data
withr::with_seed(1234, {
  n <- 1500
  
  nets <- tibble(
    # Generate exogenous variables
    id = 1:n,
    income = rnorm(n, 500, 100)
  ) |> 
    # Generate health, which depends on income
    mutate(
      health = (0.1 * income) + rnorm(n, 0, 10),
      health = pmin(pmax(health, 0), 100)  # Force values to 0-100 range
    ) |> 
    # Generate net usage (treatment), which depends on health and income
    mutate(
      # Make people with high health and high income a lot more likely to
      # self-select into treatment
      income_high = income > quantile(income, 0.8),
      health_high = health > quantile(health, 0.8),
      
      # Create latent propensity scores based on logit formula
      net_prob = plogis(
        -3 + (income / 300) + (health / 40) + 
          (0.9 * income_high) + (0.5 * health_high)),
      
      # Assign to treatment
      net = rbinom(n, 1, 
                   # Fancy little trick---increase the latent propensity score by 5
                   # percentage points for people already more likely (i.e. p > 50%) to
                   # receive treatment and decrease the propensity score by 8 percentage
                   # points for those already less likely to receive treatment
                   prob = ifelse(
                     net_prob > 0.5, 
                     pmin(net_prob + 0.05, 0.98), 
                     pmax(net_prob - 0.08, 0.02))
      )
    ) |> 
    # Generate potential and actual outcomes
    mutate(
      # Individual outcomes in a world where individuals are not treated
      malaria_risk_0 = 90 + (-0.2 * health) + (-0.05 * income) + rnorm(n, 0, 6),
      malaria_risk_0 = pmin(pmax(malaria_risk_0, 0), 100),  # Force values to 0-100
      
      # Individual outcomes in a world where individuals are treated
      # Basically risk_0 + a baseline treatment effect + extra treatment effect
      # boosts because of income and health
      malaria_risk_1 = malaria_risk_0 - 
        rnorm(n, 3, 1) - (0.015 * income) - (0.09 * health),
      malaria_risk_1 = pmin(pmax(malaria_risk_1, 0), 100),  # Force values to 0-100
      
      # Round stuff
      malaria_risk_0 = round(malaria_risk_0, 1),
      malaria_risk_1 = round(malaria_risk_1, 1),
      
      ice = malaria_risk_1 - malaria_risk_0,
      
      # Actual realized outcome
      malaria_risk = ifelse(net == 1, malaria_risk_1, malaria_risk_0)
    ) |> 
    # Round more stuff
    mutate(across(c(income, health), ~round(., 0))) |> 
    # Only keep some columns
    select(
      id, income, income_high, health, health_high, net,
      starts_with("malaria"), ice
    )
})

write_csv(nets, "mosquito_nets_v2.csv")

# Overview of data
nets |> glimpse()

###### Potential Outcomes Table ######
excerpt <- nets |> 
  slice(c(3, 14, 15, 29, 4, 35, 37, 73)) |> 
  arrange(id)

excerpt |> 
  add_row(id = NA) |> 
  select(
    id, income, health, net, 
    malaria_risk_1, malaria_risk_0, ice, malaria_risk
  ) |> 
  gt() |> 
  sub_missing(missing_text = "…") |>
  fmt_number(
    columns = c(starts_with("malaria_risk"), ice),
    decimals = 1
  ) |> 
  fmt_number(columns = c(income, health), decimals = 0) |> 
  
  # Column labels
  cols_label(
    id = "ID",
    income = md("$Z_{1_i}$"),
    health = md("$Z_{2_i}$"),
    net = md("$X_i$"),
    malaria_risk_0 = md("$Y^0_i$"),
    malaria_risk_1 = md("$Y^1_i$"),
    malaria_risk = md("$Y_i$"),
    ice = md("$Y^1_i - Y^0_i$")
  ) |>
  
  # Level 1 spanner labels
  tab_spanner(
    label = "Income", columns = income, 
    level = 1, id = "level1_a"
  ) |> 
  tab_spanner(
    label = "Health", columns = health, 
    level = 1, id = "level1_b"
  ) |> 
  tab_spanner(
    label = "Net use", columns = net, 
    level = 1, id = "level1_c"
  ) |> 
  tab_spanner(
    label = "Potential outcomes",
    columns = c(malaria_risk_1, malaria_risk_0),
    level = 1, id = "level1_d"
  ) |> 
  tab_spanner(
    label = "ICE or \\(\\delta_i\\)", columns = ice, 
    level = 1, id = "level1_e"
  ) |> 
  tab_spanner(
    label = "Outcome", columns = malaria_risk, 
    level = 1, id = "level1_f"
  ) |> 
  
  # Level 2 spanner labels
  tab_spanner(
    label = "Confounders",
    columns = c(income, health),
    level = 2, id = "level2_a"
  ) |> 
  tab_spanner(
    label = "Treatment", columns = net, 
    level = 2, id = "level2_b"
  ) |> 
  tab_spanner(
    label = "Unobservable",
    columns = c(malaria_risk_1, malaria_risk_0, ice), 
    level = 2, id = "level2_c"
  ) |> 
  tab_spanner(
    label = "Realized", columns = malaria_risk, 
    level = 2, id = "level2_d") |> 
  
  # Style stuff
  tab_style(
    style = cell_text(align = "center"),
    locations = cells_column_labels()
  ) |> 
  tab_style(
    style = cell_text(align = "center"),
    locations = cells_body()
  ) |> 
  tab_style(
    style = cell_text(weight = "bold"),
    locations = list(
      cells_column_spanners(spanners = starts_with("level1")),
      cells_column_labels(columns = "id")
    )
  ) |> 
  tab_style(
    style = cell_text(style = "italic"),
    locations = cells_column_spanners(spanners = starts_with("level2"))
  ) |> 
  
  tab_style(
    style = list(
      cell_fill(color = clrs[4], alpha = 0.5)
    ),
    locations = cells_body(rows = net == 1)
  ) |> 
  tab_footnote(
    footnote = "ICE = individual causal effect",
    locations = cells_column_spanners(spanners = "level1_e")
  ) |> 
  opt_footnote_marks(marks = "standard") |> 
  opt_horizontal_padding(scale = 3) |> 
  opt_table_font(font = "Jost")

###### Calculate Effect Types ######

# ATT and ATU
effect_types <- nets |> 
  group_by(net) |> 
  summarize(
    effect = mean(ice),
    n = n()
  ) |> 
  mutate(
    prop = n / sum(n),
    weighted_effect = effect * prop
  ) |> 
  mutate(estimand = case_match(net, 0 ~ "ATU", 1 ~ "ATT"), .before = 1)

# ATE
effect_types |> 
  summarize(ATE = sum(weighted_effect))

###### Weighting ######
# Goal: Find ATE using ITPW

# Option 1: Use Logistic Regression
# Logistic regression model to predict net usage
model_treatment <- glm(
  net ~ income + health, 
  data = nets, 
  family = binomial(link = "logit")
)

# Plug original data into the model to generate predicted probabilities
# (Here I'm using broom::augment_columns(), but you can also use 
#  marginaleffects::predictions() or even base R's predict())

# Create a new column for the weights
nets_with_weights <- augment_columns(
  model_treatment,
  data = nets,
  type.predict = "response"
) |> 
  rename(propensity = .fitted) |> 
  mutate(wts_ate = (net / propensity) + ((1 - net) / (1 - propensity)))

# See if it worked
nets_with_weights |> 
  select(id, income, health, net, malaria_risk, propensity, wts_ate)

# Option 2: Use WeightIt
# Get ATE-focused weights for the treatment model
weights_ate <- weightit(
  net ~ income + health, 
  data = nets, 
  method = "glm", 
  estimand = "ATE"
)

# Add the weights as a column in the data
nets_with_weights$wts_ate_automatic <- weights_ate$weights

# They're the same!
nets_with_weights |> 
  select(id, wts_ate, wts_ate_automatic) |> 
  head()

# ATE plots
plot_data_weights_ate <- tibble(
  propensity = weights_ate$ps,
  weight = weights_ate$weights,
  treatment = weights_ate$treat
)

ggplot() + 
  geom_histogram(data = filter(plot_data_weights_ate, treatment == 1), 
                 bins = 50, aes(x = propensity, fill = "Treated people")) + 
  geom_histogram(data = filter(plot_data_weights_ate, treatment == 0), 
                 bins = 50, aes(x = propensity, y = -after_stat(count), fill = "Untreated people")) +
  geom_hline(yintercept = 0, color = "white", linewidth = 0.25) +
  scale_x_continuous(labels = scales::label_percent()) +
  scale_y_continuous(label = abs) +
  scale_fill_manual(values = c(clrs[1], clrs[6]), guide = guide_legend(reverse = TRUE)) +
  labs(x = "Propensity", y = "Count", fill = NULL) +
  theme(
    legend.position = "top",
    legend.key.size = unit(0.65, "lines")
  )

# Show pre-weighting highlighted
ggplot() + 
  geom_histogram(data = filter(plot_data_weights_ate, treatment == 1), 
                 bins = 50, aes(x = propensity, fill = "Treated people")) + 
  geom_histogram(data = filter(plot_data_weights_ate, treatment == 0), 
                 bins = 50, aes(x = propensity, y = -after_stat(count), fill = "Untreated people")) +
  annotate(
    geom = "rect", xmin = 0, xmax = 0.55, ymin = -1, ymax = 23,
    fill = alpha(clrs[4], 0.1), color = clrs[4], linewidth = 1
  ) +
  annotate(
    geom = "text", x = 0.01, y = 21.5, 
    label = "Treated people who were\nunlikely to be treated.\nWeird!",
    color = clrs[3], fontface = "bold", hjust = 0, vjust = 1, lineheight = 1
  ) +
  annotate(
    geom = "rect", xmin = 0.58, xmax = 1, ymin = 1, ymax = -27.5,
    fill = alpha(clrs[4], 0.1), color = clrs[4], linewidth = 1
  ) +
  annotate(
    geom = "text", x = 0.99, y = -26, 
    label = "Weird!\nUntreated people who\nwere likely to be treated.",
    color = clrs[3], fontface = "bold", hjust = 1, vjust = 0, lineheight = 1
  ) +
  geom_hline(yintercept = 0, color = "white", linewidth = 0.25) +
  scale_x_continuous(labels = scales::label_percent()) +
  scale_y_continuous(label = abs) +
  scale_fill_manual(values = c(clrs[1], clrs[6]), guide = guide_legend(reverse = TRUE)) +
  labs(x = "Propensity", y = "Count", fill = NULL) +
  theme(
    legend.position = "top",
    legend.key.size = unit(0.65, "lines")
  )

# Show post-weighting
ggplot() + 
  geom_histogram(data = filter(plot_data_weights_ate, treatment == 1), 
                 bins = 50, aes(x = propensity, weight = weight, fill = "Treated pseudo-population")) + 
  geom_histogram(data = filter(plot_data_weights_ate, treatment == 0), 
                 bins = 50, aes(x = propensity, weight = weight, y = -after_stat(count), fill = "Untreated psuedo-population")) +
  geom_histogram(data = filter(plot_data_weights_ate, treatment == 1), 
                 bins = 50, aes(x = propensity, fill = "Treated people")) + 
  geom_histogram(data = filter(plot_data_weights_ate, treatment == 0), 
                 bins = 50, aes(x = propensity, y = -after_stat(count), fill = "Untreated people")) +
  annotate(
    geom = "text", x = 0.5, y = 60, label = "More weight here", 
    hjust = 0.5, fontface = "bold", color = clrs[3]
  ) +
  annotate(
    geom = "segment", x = 0.48, xend = 0.17, y = 55, yend = 25, color = "white", 
    arrow = arrow(angle = 15, length = unit(0.5, "lines")), linewidth = 2.5
  ) +
  annotate(
    geom = "segment", x = 0.48, xend = 0.17, y = 55, yend = 25, color = clrs[3], 
    arrow = arrow(angle = 15, length = unit(0.5, "lines")), linewidth = 0.5
  ) +
  annotate(
    geom = "segment", x = 0.52, xend = 0.8, y = 55, yend = -20, color = "white", 
    arrow = arrow(angle = 15, length = unit(0.5, "lines")), linewidth = 2.5
  ) +
  annotate(
    geom = "segment", x = 0.52, xend = 0.8, y = 55, yend = -20, color = clrs[3], 
    arrow = arrow(angle = 15, length = unit(0.5, "lines")), linewidth = 0.5
  ) +
  geom_hline(yintercept = 0, color = "white", linewidth = 0.25) +
  scale_x_continuous(labels = scales::label_percent()) +
  scale_y_continuous(label = abs) +
  scale_fill_manual(
    values = c(clrs[1], colorspace::lighten(clrs[1], 0.5), clrs[6], colorspace::lighten(clrs[6], 0.65)), 
    guide = guide_legend(reverse = FALSE, nrow = 2)) +
  labs(x = "Propensity", y = "Count", fill = NULL) +
  theme(
    legend.position = "top",
    legend.key.size = unit(0.65, "lines")
  )

# To find the ATE using these weights, we can fit a regression model for the outcome. 
# We’ll then use an approach called g-computation to find the marginal causal effect 
# (i.e. the causal effect of flipping net from “no” to “yes”). G-computation is a neat approach 
# that feels more like the original potential outcomes framework we looked at earlier. It involves a few steps:
# 1. Fit a regression model predicting the outcome, like lm(malaria_risk ~ net, weights = w_ate, ...)
# 2. Use the model to generate a set of predictions where we pretend that every person used a net (i.e. set net = 1)
# 3. Use the model to generate a set of predictions where we pretend that every person did not use a net (i.e. set net = 0)
# 4. Calculate the difference in the two sets of predictions to create a sort of predicted individual causal effect, 
# then find the average of that. This is the ATE.
# {marignaleffects} makes this incredibly easy 
# All we have to do is use avg_comparisons(model_outcome, variables = list(net = 0:1)), or even more simply a
# vg_comparisons(model_outcome, variables = "net"). This will do all the work of setting everyone’s treatment 
# to 0, to 1, and finding the difference. We can do fancier things with it too, like finding robust and/or 
# clustered standard errors, or bootstrapping the standard errors.

# G-computation
# Fit an outcome model using the inverse probability weights
model_outcome_ate <- lm(
  malaria_risk ~ net, 
  data = nets_with_weights, 
  weights = wts_ate_automatic
)

# Automatic g-computation with {marginaleffects}! This sets the "net" column to
# each possible value of net (0 and 1) for the whole dataset, then calculates
# the difference between the two sets of predictions

avg_comparisons(model_outcome_ate, variables = "net")

# Based on our inverse probability weighting approach, the ATE from observational data 
# (i.e. not using the unknowable individual potential outcomes) is −14.7, which is super close to the true ATE of −15.0

# Finding ATT or ATU
# ATE shows the effect of the net program for everyone, even people who have no need for a net. 
# If we’re interested in making this a universal program, the ATE is useful. If we’re interested 
# in what the program is currently doing for people using it, or what would happen if we expanded 
# it to just those not using it, we’d need to find the ATT or ATU instead.
# When reweighting the data for estimating the ATT, we should not do anything that messes with the weights 
# of the treated group. Since we’re interested in the effect of the program on just the treated people, 
# we need to create a pseudo-control group that looks like the treated group.

# Get ATT-focused weights for the treatment model
weights_att <- weightit(
  net ~ income + health, 
  data = nets, 
  method = "glm", 
  estimand = "ATT"
)

# Add the weights as a column in the data
nets_with_weights$wts_att <- weights_att$weights

# Mirrored histogram comparing the distribution of treated people with the reweighted pseudo-population 
# of untreated people; because this is for the ATT, treated people were not reweighted

plot_data_weights_att <- tibble(
  propensity = weights_att$ps,
  weight = weights_att$weights,
  treatment = weights_att$treat
)

ggplot() + 
  geom_histogram(data = filter(plot_data_weights_att, treatment == 0), 
                 bins = 50, aes(x = propensity, weight = weight, y = -after_stat(count), fill = "Untreated psuedo-population")) +
  geom_histogram(data = filter(plot_data_weights_att, treatment == 1), 
                 bins = 50, aes(x = propensity, fill = "Treated people")) + 
  geom_histogram(data = filter(plot_data_weights_att, treatment == 0), 
                 bins = 50, alpha = 0.5, 
                 aes(x = propensity, y = -after_stat(count), fill = "Untreated people")) +
  geom_hline(yintercept = 0, color = "white", linewidth = 0.25) +
  scale_x_continuous(labels = scales::label_percent()) +
  scale_y_continuous(label = abs) +
  scale_fill_manual(
    values = c(clrs[1], "grey50", colorspace::lighten(clrs[6], 0.65)), 
    guide = guide_legend(reverse = FALSE, nrow = 2)) +
  labs(x = "Propensity", y = "Count", fill = NULL) +
  theme(
    legend.position = "top",
    legend.key.size = unit(0.65, "lines")
  )

# The main difference is that when we do the g-computation, we only look at treated people. 
# Using only this part of the population, we’ll set all their values to 1, then set all their values to 0, 
# and then find the difference.

# Do g-computation *only* on treated observations
avg_comparisons(
  model_outcome_att, 
  variables = "net", 
  newdata = filter(nets, net == 1)
)

# ATU
# If we’re interested in what would happen if we expanded the program to people not using it, 
# we can find the ATU. The process is the same as the ATT, just in reverse. 

# Get ATU-focused weights for the treatment model
weights_atu <- weightit(
  net ~ income + health, 
  data = nets, 
  method = "glm", 
  estimand = "ATC"  # WeightIt calls this ATC instead of ATU
)

# Add the weights as a column in the data
nets_with_weights$wts_atu <- weights_atu$weights

# Add to our histogram
plot_data_weights_atu <- tibble(
  propensity = weights_atu$ps,
  weight = weights_atu$weights,
  treatment = weights_atu$treat
)

ggplot() + 
  geom_histogram(data = filter(plot_data_weights_atu, treatment == 1), 
                 bins = 50, aes(x = propensity, weight = weight, fill = "Treated pseudo-population")) + 
  geom_histogram(data = filter(plot_data_weights_atu, treatment == 1), 
                 bins = 50, alpha = 0.5, aes(x = propensity, fill = "Treated people")) + 
  geom_histogram(data = filter(plot_data_weights_atu, treatment == 0), 
                 bins = 50, aes(x = propensity, y = -after_stat(count), fill = "Untreated people")) +
  geom_hline(yintercept = 0, color = "white", linewidth = 0.25) +
  scale_x_continuous(labels = scales::label_percent()) +
  scale_y_continuous(label = abs) +
  scale_fill_manual(
    values = c("grey50", colorspace::lighten(clrs[1], 0.5), clrs[6], colorspace::lighten(clrs[6], 0.65)), 
    guide = guide_legend(reverse = FALSE, nrow = 2)) +
  labs(x = "Propensity", y = "Count", fill = NULL) +
  theme(
    legend.position = "top",
    legend.key.size = unit(0.65, "lines")
  )

# The outcome model process is the same as before—we just use the ATU-specific weights. 
# When doing the g-computation, we only use untreated people.

# Fit an outcome model using the ATU weights
model_outcome_atu <- lm(
  malaria_risk ~ net, 
  data = nets_with_weights, 
  weights = wts_atu
)

# Do g-computation *only* on untreated observations
avg_comparisons(
  model_outcome_atu, 
  variables = "net", 
  newdata = filter(nets, net == 0)
)


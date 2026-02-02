# General strategy for counterfactual estimation works in 4 steps:
# 1. Use the untreated observations to fit a model of the outcome .
# 2. Use this model to predict the counterfactual potential outcomes for treated observations.
# 3. Estimate the individual treatment effects in the subset of actually treated observations:
# 4. Average the individual estimates to obtain an ATT of interest.
# In Step 1, analysts can use the model of their choice, and each model will carry its own identifying 
# assumptions. In their article, Liu, Wang, and Xu (2020) consider three possibilities: 
#     (1) Two-way fixed effects
#     (2) Interactive fixed effects
#     (3) Matrix completion.
# Step 4 allows considerable flexibility in interpretation: We can average individual estimates in a given 
# year to estimate the year-specific ATT, or within cohorts to obtain a cohort-specific estimate. 
# This will allow us to draw nice plots.
# We will compute counter factual estimates using all three of the approaches

library('data.table')
library('fixest')
library('did')
library('fect')

set.seed(1000)


theme_fancy <- function() {
  theme_minimal(base_family = "Asap Condensed") +
    theme(panel.grid.minor = element_blank(),
          legend.position = "bottom")
}

# Simulate data and preview first 5 rows
source("/Users/karstenwalker/Documents/R/Causal/hte_simulated_data.R")

dat = simulate()

dat[1:5, .(unit, time, Y, D, X1, X2, δ)]

# Compute a simple two-way fixed effects model using the fixest::feols function,

mod = feols(Y ~ D + X1 + X2 | unit + time, dat)

coef(mod)

# Save this result and the capital-T “Truth” in a data.table so we can plot them later:

att = dat[, .(Truth = mean(δ)), by = `Time after treatment`]

att[`Time after treatment` >= 0, `2-way FE` := coef(mod)["D"]][
  `Time after treatment` < 0, Truth := 0] 

# Compute counterfactual outcomes under no treatment for units which have actually received the treatment, 
# and use those counterfactual outcomes to estimate individual-level treatment effects which can then be aggregated.
# Calculates the gap between  and the counterfactual for every unit, including the untreated one. In principle, 
# the estimator should produce estimates of 0 in the pre-treatment period. This idea allows us to conduct 
#diagnostic tests.

###### Two-Way Fixed Effects Counterfactual Estimator (FEct) ######
# Step 1: Use the untreated observations to fit a model of the outcome

mod = feols(Y ~ X1 + X2 | unit + time, dat[D == 0])

# Step 2: Use this model to predict the counterfactual potential outcomes for treated observations.

dat[, Y0hat := predict(mod, newdata = dat)]

# Step 3: Estimate the individual treatment effects in the subset of actually treated observations:
  
dat[, treatment_effect := Y - Y0hat]

# Step 4: Average the individual estimates to obtain an ATT. Here, we calculate a specific ATT for 
# each time period relative to the first period when a unit is treated:
  
tmp = dat[, .(`LWX: Fixed Effects` = mean(treatment_effect)), by = `Time after treatment`]

att = merge(att, tmp, all = TRUE)

###### Interactive fixed-effect counterfactual estimator ######
mod_ife = fect(
  Y ~ D + X1 + X2, 
  data = dat, 
  index = c("unit", "time"),
  method = "ife")

# save results for plotting later
tmp = data.table(
  `Time after treatment` = mod_ife$time,
  `LWX: Interacted Fixed Effects` = mod_ife$att)

att = tmp[att, on = .(`Time after treatment`)]

###### Matrix completion counterfactual estimator ######
mod_mc = fect(
  Y ~ D + X1 + X2, 
  data = dat, 
  index = c("unit", "time"),
  method = "mc")

# save results for plotting later
tmp = data.table(
  `Time after treatment` = mod_mc$time,
  `LWX: Matrix completion` = mod_mc$att)
att = tmp[att, on = .(`Time after treatment`)]

###### Group-Time Average Treatment Effect ######
# To estimate the group-time average treatment effect, Callaway and Sant’Anna (2020) impose a few assumptions: 
#   (1) irreversibility of treatment
#   (2) random sampling
#   (3) limited treatment anticipation
#   (4) overlap
#   (5) conditional parallel trends based on a “never-treated” or a “not-yet-treated” group.
# If these assumptions hold, the authors show that the  can be recovered by carefully selecting 
# the groups that are leveraged for comparison. The authors propose three alternative approaches 
# to doing this: outcome regression, inverse probability weighting, or doubly-robust. 
# In a first step, the analyst estimates an outcome (and/or treatment assignment) model. 
# Then, they plug-in the fitted estimates from that model into the sample analogue of the  of interest.

dat[, treat := D][
  , first.treat := fifelse(T_0i == 35, 0, T_0i)]

mod = att_gt(
  yname = "Y",
  gname = "first.treat",
  idname = "unit",
  tname = "time",
  xformla = ~1 + X1 + X2,
  data = dat)

# One of the main benefits of the Group-Time Average Treatment Effect is that it can be aggregated 
# to provide meaningful inference at different levels. For instance, in the previous section we calculated 
# ATTs at the time-period level. We can do the same here with the aggte function:
  
es <- aggte(mod, type = "dynamic")

# save results for plotting later
tmp = data.table(
  `Time after treatment` = es$egt, 
  `Callaway & Sant'Anna` = es$att.egt)
att = tmp[att, on = .(`Time after treatment`)]

###### Plot Comparison of Methods ######

tmp = melt(att, id.vars = "Time after treatment", 
           variable.name = "Model", value.name = "ATT")

palette = c("#CC6677", "#332288", "#DDCC77", "#117733", "#88CCEE", "#882255", "#44AA99")

ggplot(tmp, aes(`Time after treatment`, ATT, color = Model)) +
  geom_line() +
  scale_x_continuous(breaks = seq(-30, 15, by = 5)) +
  scale_color_manual(values = palette)  +
  theme_fancy() +
  theme(panel.grid.minor = element_blank())
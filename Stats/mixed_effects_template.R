
# Background
# We run an A/B test on a new learning feature. Outcome is whether a user completes a study session
# successfully (completed, binary). Users have repeated sessions across time; sessions are nested 
# within users, and users belong to cohorts/platforms. We use random intercepts (and optionally random slopes)
# to account for user heterogeneity.

# Extensions
# For count outcomes (e.g., #cards studied), use glmer(..., family=poisson or nbinom2 via glmmTMB).
# Consider adding time trends (date as numeric) or splines.
# For robust inference / Bayesian: brms or rstanarm.
# For more formal marginal effects: marginaleffects package (works well with many model types).


# Load Packages and set theme
library(tidyverse)
library(lubridate)
library(lme4)
library(broom.mixed)
library(modelsummary)
library(performance)

set.seed(42)

theme_fancy <- function() {
  theme_minimal(base_family = "Asap Condensed") +
    theme(panel.grid.minor = element_blank(),
          legend.position = "bottom")
}

# Synthetic data 
# A/B test for a new "Smart Round" feature on Quizlet.
# Outcome: completed (0/1) for each study session.
# Repeated sessions per user over multiple days.

n_users   <- 800
n_days    <- 28
platforms <- c("ios", "android", "web")
cohorts   <- c("new", "returning")
subjects  <- c("Spanish", "Biology", "History", "Math")

# User-level table
users <- tibble(
  user_id      = sprintf("u%04d", 1:n_users),
  platform     = sample(platforms, n_users, replace = TRUE, prob = c(0.35, 0.35, 0.30)),
  cohort       = sample(cohorts, n_users, replace = TRUE, prob = c(0.55, 0.45)),
  age_days     = pmax(1L, rpois(n_users, lambda = 30)),
  baseline_eng = rnorm(n_users, mean = 0, sd = 1)   # latent engagement propensity
) %>%
  mutate(
    # randomize at user-level
    treatment = rbinom(n_users, 1, 0.5),
    treatment = factor(treatment, levels = c(0, 1), labels = c("control", "smart_round"))
  )

# Session-level table
sessions <- users %>%
  # each user has a Poisson-ish number of sessions across the window
  mutate(n_sessions = pmax(1L, rpois(n_users, lambda = 10))) %>%
  select(user_id, treatment, platform, cohort, age_days, baseline_eng, n_sessions) %>%
  uncount(n_sessions, .id = "session_in_user") %>%
  mutate(
    date = as_date("2025-10-01") + sample.int(n_days, n(), replace = TRUE) - 1L,
    dow  = wday(date, label = TRUE, abbr = TRUE),
    subject = sample(subjects, n(), replace = TRUE, prob = c(0.35, 0.25, 0.20, 0.20)),
    session_length_min = pmax(1, rlnorm(n(), meanlog = log(6), sdlog = 0.5)),
    prior_sessions_7d  = rpois(n(), lambda = 2)
  )

# Create user random intercepts and (optional) random slopes for treatment
user_re <- users %>%
  transmute(
    user_id,
    b0_user = rnorm(n(), 0, 0.9),   # random intercept
    b1_trt  = rnorm(n(), 0, 0.25)   # random slope for treatment heterogeneity
  )

# Data generating process (logit scale)
# Fixed effects: treatment + covariates + some interactions
# Random effects: (1 + treatment | user_id)
dat <- sessions %>%
  left_join(user_re, by = "user_id") %>%
  mutate(
    # encode some fixed effects in a transparent way
    trt = if_else(treatment == "smart_round", 1, 0),
    is_new = if_else(cohort == "new", 1, 0),
    is_mobile = if_else(platform %in% c("ios", "android"), 1, 0),
    
    # Subject difficulty (latent)
    subj_effect = case_when(
      subject == "Spanish" ~  0.10,
      subject == "Biology" ~ -0.05,
      subject == "History" ~ -0.10,
      subject == "Math"    ~ -0.15,
      TRUE ~ 0
    ),
    
    # True linear predictor
    eta =
      -0.4 +
      0.18 * trt +                         # average treatment effect
      0.22 * baseline_eng +
      0.03 * log1p(prior_sessions_7d) +
      0.10 * log1p(session_length_min) +
      (-0.12) * is_new +                   # new users slightly less likely to complete
      0.08 * is_mobile +
      0.10 * trt * is_new +                # treatment helps new users a bit more
      subj_effect +
      b0_user +                            # random intercept
      b1_trt * trt                         # random slope component
  ) %>%
  mutate(
    p_completed = plogis(eta),
    completed   = rbinom(n(), 1, p_completed)
  ) %>%
  select(
    user_id, date, dow, platform, cohort, subject, treatment,
    session_length_min, prior_sessions_7d,
    completed
  ) %>%
  mutate(
    across(c(platform, cohort, subject, dow), as.factor)
  )

glimpse(dat)

# EDA 
# Completion rate by treatment
dat %>%
  group_by(treatment) %>%
  summarise(
    n = n(),
    completion_rate = mean(completed),
    .groups = "drop"
  ) %>%
  ggplot(aes(x = treatment, y = completion_rate)) +
  geom_col() +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  labs(
    title = "Completion rate by treatment",
    x = NULL, y = "Completion rate"
  ) +
  theme_fancy()

# Completion rate over time (smoothed)
dat %>%
  group_by(date, treatment) %>%
  summarise(completion_rate = mean(completed), n = n(), .groups = "drop") %>%
  ggplot(aes(x = date, y = completion_rate, color = treatment)) +
  geom_line(alpha = 0.6) +
  geom_smooth(se = FALSE, method = "loess", formula = y ~ x) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  labs(
    title = "Completion rate over time",
    x = NULL, y = "Completion rate", color = NULL
  ) +
  theme_fancy()

# Heterogeneity view: by cohort
dat %>%
  group_by(cohort, treatment) %>%
  summarise(completion_rate = mean(completed), n = n(), .groups = "drop") %>%
  ggplot(aes(x = cohort, y = completion_rate, fill = treatment)) +
  geom_col(position = position_dodge(width = 0.8)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  labs(
    title = "Completion rate by cohort and treatment",
    x = NULL, y = "Completion rate", fill = NULL
  ) +
  theme_fancy()

# 3 Mixed effects models 
# Outcome is binary -> use glmer with binomial(link="logit")
# Baseline model: random intercept by user

# Most Quizlet experiment analyses need:
# - A credible estimate of treatment effect (ATE)
# - Optional: evidence for heterogeneity by cohort (or another segment)
# - Correct handling of repeated sessions per user

# Random intercept (1 | user_id) is non-negotiable
# Random slope is optional (only if supported)
# treatment * cohort is optional (only if it improves fit / is needed)

# Choose random effects structure first. Random effects determine how much pooling/shrinkage
# you assume and whether you’re estimating user-level treatment heterogeneity.

# rescale data
dat_scaled <- dat %>%
  mutate(
    session_length_z = as.numeric(scale(log1p(session_length_min))),
    prior_sessions_z = as.numeric(scale(log1p(prior_sessions_7d)))
  )

# Random intercept, no interaction, unscaled covariates, no dow, baseline ATE
# Redundat to m1_int, should be dropped due to scaling
m1 <- glmer(
  completed ~ treatment + cohort + platform + subject +
    log1p(session_length_min) + log1p(prior_sessions_7d) +
    (1 | user_id),
  data = dat_scaled,
  family = binomial()
)

# Same as model 1 with uncorrelated covariates, no dow, still ATE only
m1_int<- glmer(
  completed ~ treatment + cohort + platform + subject +
    session_length_z + prior_sessions_z +
    (1 | user_id),
  data = dat_scaled,
  family = binomial()
)

# No random slope, has dow, baseline ATE model
# m2_int without log1p covariates, redundant to m2 and should be dropped due to scaling
m1_ran <- glmer(
  completed ~ treatment * cohort + platform + subject + dow +
    session_length_z + prior_sessions_z +
    (1 | user_id),
  data = dat_scaled,
  family = binomial(),
  control = glmerControl(optimizer = "bobyqa")
)

# Add key interaction + day-of-week + random slope for treatment by user
# Full correlated random slope, interaction + dowm, unscaled covariates
# Highest complexity andrisk of identifiability issues

m2 <- glmer(
  completed ~ treatment * cohort + platform + subject + dow +
    log1p(session_length_min) + log1p(prior_sessions_7d) +
    (1 + treatment | user_id),
  data = dat_scaled,
  family = binomial(),
  control = glmerControl(optimizer = "bobyqa")
)

# With uncorrelated random slopes, scaled covariates, interaction + dow
m2_ran <- glmer(
  completed ~ treatment * cohort + platform + subject + dow +
    session_length_z + prior_sessions_z +
    (1 | user_id) + (0 + treatment | user_id),
  data = dat_scaled,
  family = binomial(),
  control = glmerControl(optimizer = "bobyqa")
)

# Random intercept only, interaction + dow, scaled covariates, clean, stable heterogeneity model
# Numerically stable, correctly accounts for repeated sessions, allows decision-relevant heterogeneity
# avoids overfitting user-level slope variance, easy to explain and defend
m2_int <- glmer(
  completed ~ treatment * cohort + platform + subject + dow +
    session_length_z + prior_sessions_z +
    (1 | user_id),
  data = dat_scaled,
  family = binomial()
)

# 4. Model summaries 
# For logistic mixed models, coefficients are on log-odds scale by default.
# You can also exponentiate to get odds ratios.

# Compare random-effects structures using AIC + stability, preferably keeping fixed effects the same.
# Rule of thumb:
# - If random slope improves AIC by < ~2–4, it’s not worth the complexity.
# - If it improves AIC by > ~10 and the fit is non-singular, it may be worth it.
# - If you used correlated slopes (1 + treatment | user_id) and it’s unstable, prefer uncorrelated slopes:
#  (1 | user_id) + (0 + treatment | user_id).
# If your org wants a single number (“did it work?”), choose the model that yields the cleanest ATE 
# without overfitting. If your org is going to act differently by cohort (e.g., ship only for new users),
# you need treatment * cohort (i.e., an m2 variant), even if global fit gain is modest.

re_compare <- AIC(m1, m1_int, m1_ran, m2, m2_int, m2_ran) %>% 
  as_tibble()

re_compare

# Ranking
candidates <- list(
  m1 = m1, m1_ran = m1_ran, m1_int = m1_int,
  m2 = m2, m2_ran = m2_ran, m2_int = m2_int
)

ranked <- imap_dfr(candidates, ~ tibble(
  model = .y,
  AIC = AIC(.x),
  singular = isTRUE(performance::check_singularity(.x))
)) %>%
  filter(!singular) %>%
  arrange(AIC)

ranked

# Likelihood ratio test, chisq is how much better the model fits, P>.05 no sig effect
# P value shows there is effect but it is not significant
anova(m1_int, m2_int, test = "Chisq")

# Check correlation between intercept and slope, values approaching 1 are bad
vc_table <- list(
  m1     = m1,
  m1_ran = m1_ran,
  m1_int = m1_int,
  m2     = m2,
  m2_ran = m2_ran,
  m2_int = m2_int) %>%
  imap_dfr(
    ~ tidy(.x, effects = "ran_pars") %>%
      mutate(model = .y),
    .id = NULL) %>%
  select(
    model,
    group,
    term,
    estimate) %>%
  arrange(model, group, term)

vc_table

models <- list(
  "Random intercept" = m1,
  "Scaled Random intercept" = m1_ran,
  "Intercept only"= m1_int,
  "Interaction + random intercept+slope"  = m2,
  "Scaled Interaction + random intercept+slope"  = m2_ran,
  "Intercept and interaction"  = m2_int
)

modelsummary(
  models,
  stars = TRUE,
  statistic = "({std.error})",
  output = "markdown"
)

# Odds ratio table (exponentiated coefficients)
modelsummary(
  models,
  exponentiate = TRUE,
  stars = TRUE,
  statistic = "({std.error})",
  output = "markdown"
)

# non modelsummary, all models
models <- list(
  m1 = m1_intercept_only,
  m2 = m2_intercept_only
)

or_multi <- imap_dfr(
  models,
  ~ broom.mixed::tidy(
    .x,
    effects = "fixed",
    conf.int = TRUE ) %>%
    filter(term != "(Intercept)") %>%
    mutate(
      odds_ratio = exp(estimate),
      conf.low   = exp(conf.low),
      conf.high  = exp(conf.high),
      model      = .y
    ),
  .id = NULL) %>%
  select(
    model,
    term,
    odds_ratio,
    conf.low,
    conf.high,
    p.value) %>%
  arrange(term, model)

or_multi

# Wide format for printing
or_multi_wide <- or_multi %>%
  mutate(
    or_ci = sprintf(
      "%.2f [%.2f, %.2f]",
      odds_ratio, conf.low, conf.high
    )
  ) %>%
  select(model, term, or_ci) %>%
  pivot_wider(
    names_from = model,
    values_from = or_ci
  )

or_multi_wide

# 5. Tidy fixed effects & plot estimates 
tidy_m2 <- broom.mixed::tidy(m2, effects = "fixed", conf.int = TRUE) %>%
  filter(term != "(Intercept)")

tidy_m2 %>%
  ggplot(aes(x = estimate, y = fct_reorder(term, estimate))) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.6) +
  geom_point() +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
  labs(
    title = "Fixed effect estimates (log-odds) with 95% CI",
    x = "Estimate (log-odds)", y = NULL
  ) +
  theme_fancy()

# If you prefer odds ratios:
tidy_m2_or <- tidy_m2 %>%
  mutate(
    estimate = exp(estimate),
    conf.low = exp(conf.low),
    conf.high = exp(conf.high)
  )

tidy_m2_or %>%
  ggplot(aes(x = estimate, y = fct_reorder(term, estimate))) +
  geom_vline(xintercept = 1, linetype = "dashed", alpha = 0.6) +
  geom_point() +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
  scale_x_log10() +
  labs(
    title = "Fixed effect estimates (odds ratios) with 95% CI",
    x = "Odds ratio (log scale)", y = NULL
  ) +
  theme_fancy()

# 6. Random effects diagnostics 
# Random intercepts by user
ranef_m2 <- ranef(m2)$user_id %>%
  as_tibble(rownames = "user_id")

ranef_m2 %>%
  ggplot(aes(x = `(Intercept)`)) +
  geom_histogram(bins = 40) +
  labs(
    title = "Distribution of user random intercepts",
    x = "Random intercept (user)", y = "Count"
  ) +
  theme_fancy()

# Random slopes for treatment (if present)
if ("treatmentsmart_round" %in% names(ranef_m2)) {
  ranef_m2 %>%
    ggplot(aes(x = treatmentsmart_round)) +
    geom_histogram(bins = 40) +
    labs(
      title = "Distribution of user random slopes for treatment",
      x = "Random slope (treatment)", y = "Count"
    ) +
    theme_minimal()
}

# 7. Predicted probabilities (marginal effects-style) 
# Create a small prediction grid, holding other vars at typical values
pred_grid <- dat %>%
  summarise(
    session_length_min = median(session_length_min),
    prior_sessions_7d  = median(prior_sessions_7d)
  ) %>%
  crossing(
    treatment = factor(c("control", "smart_round"), levels = levels(dat$treatment)),
    cohort    = levels(dat$cohort),
    platform  = factor("web", levels = levels(dat$platform)),
    subject   = factor("Spanish", levels = levels(dat$subject)),
    dow       = factor("Mon", levels = levels(dat$dow)),
    user_id   = dat$user_id[1]  # needed for RE structure; we will set re.form=NA below
  )

pred_grid <- pred_grid %>%
  mutate(
    p = predict(m2, newdata = ., type = "response", re.form = NA)  # population-level
  )

pred_grid %>%
  ggplot(aes(x = cohort, y = p, color = treatment, group = treatment)) +
  geom_point(position = position_dodge(width = 0.2), size = 2) +
  geom_line(position = position_dodge(width = 0.2)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  labs(
    title = "Predicted completion probability by cohort (population-level)",
    x = NULL, y = "Predicted P(completed)", color = NULL
  ) +
  theme_fancy()

###### Reporting outcome ####
# Statistically
#   1. Strong evidence for user-level baseline heterogeneity
#   2. Strong evidence for a single, stable average treatment effect
# No evidence supporting:
#   1. cohort-specific treatment effects
#   2. day-of-week effects (conditional on other covariates)
#   3. user-level treatment heterogeneity
# Substantively / product-wise
# The treatment effect is consistent across cohorts and users. There is no justification for segmenting 
# rollout or interpretation by cohort. This is a positive, actionable finding, not a null result.
# One-paragraph methods justification (
# We evaluated six candidate mixed-effects logistic models varying in fixed-effect complexity and 
# random-effects structure. Models with random slopes increased complexity without improving fit. 
# We then tested whether allowing cohort-specific treatment effects and day-of-week controls improved model 
# fit using a likelihood ratio test. The interaction model did not improve fit 
# (χ² = 4.94, df = 7, p = 0.67) and had a higher AIC. We therefore selected a parsimonious 
# random-intercept model with scaled covariates, which provided the best balance of stability, fit, 
# and interpretability.
# Senior-level takeaway (what you should say out loud)
# “We explicitly tested for heterogeneity and temporal effects. The data doesn’t support them, so we report a single, stable average treatment effect.”
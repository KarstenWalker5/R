# IV Boilerplate (Synthetic Data + DAG + 2SLS) 
# 1. Generates a synthetic dataset with endogeneity:
#    - D is endogenous because an unobserved U affects both D and Y.
#    - Z is a valid instrument: it affects D (relevance) and affects Y only through D (exclusion).
# 2. Visualizes the assumed causal structure with a DAG.
# 3. Estimates OLS vs IV/2SLS, then runs common IV diagnostics:
#    - First-stage strength (weak instrument check)
#    - Wu–Hausman endogeneity test (whether OLS is inconsistent / D is endogenous)
#    - Sargan overidentification test (only if you have >1 instrument)
## Notes for adapting to real work:
# - Define D precisely (eligibility, exposure, usage intensity).
# - Choose Z with a rock-solid story (rollout/randomization/threshold/encouragement).
# - Add clustering/FE as needed (classroom, teacher, school, cohort, time).
# - Pre-specify estimand (LATE vs ATE), compliance definition, exclusions, guardrails.
# - Validate: instrument relevance, balance checks, sensitivity analyses.

# Load packages
suppressPackageStartupMessages({
  library(tidyverse)
  library(broom)
  library(AER)        # ivreg
  library(sandwich)   # robust vcov
  library(lmtest)     # coeftest
  library(dagitty)    # DAG specification
  library(ggdag)      # DAG visualization (ggplot2)
})

# Theme
theme_fancy <- function() {
  theme_minimal(base_family = "Asap Condensed") +
    theme(panel.grid.minor = element_blank(),
          legend.position = "bottom")
}

set.seed(42)

# 1. Synthetic Dataset Generator
# Story:
#  Treatment D: "Personalized Practice" exposure (e.g., feature eligibility/uptake)
#  Outcome Y: Weekly Active Learner indicator / minutes / completion score (continuous for demo)
#  Confounder U: latent motivation/ability (unobserved)
#  Instrument Z: exogenous "rollout/assignment pressure" (e.g., server-side eligibility flag),
#                  affects D but not Y except through D.
#  Observables X: prior activity, device type, country, etc.
# We know the TRUE causal effect of D on Y (set below to 0.55), so we can verify that OLS is biased
# when D is endogenous (due to U) and IV/2SLS recovers the causal effect under IV assumptions.

n <- 8000

df <- tibble(
  user_id = 1:n,
  prior_sessions = rpois(n, lambda = 3),
  device_mobile = rbinom(n, 1, 0.65),
  country_tier = sample(c("tier1", "tier2", "tier3"), n, replace = TRUE, prob = c(0.5, 0.3, 0.2)),
  # latent confounder (unobserved in real life; we keep it to validate in sim)
  U = rnorm(n),
  # Instrument: e.g., randomized-ish rollout intensity at the edge (exogenous)
  Z = rbinom(n, 1, 0.5)
) %>%
  mutate(
    # Map country tiers to numeric effect
    country_tier_num = case_when(
      country_tier == "tier1" ~ 0.20,
      country_tier == "tier2" ~ 0.00,
      country_tier == "tier3" ~ -0.15
    ),
    
    # Treatment assignment / uptake equation (first stage):
    # D depends on Z strongly (relevance), also depends on U and X (confounding)
    D = 0.0 +
      0.70 * Z +                        # strong relevance
      0.18 * log1p(prior_sessions) +
      0.10 * device_mobile +
      0.25 * U +                        # confounding: motivated learners more likely to take treatment
      0.20 * country_tier_num +
      rnorm(n, sd = 0.9),
    
    # Outcome equation:
    # Y depends on D causally (true effect = 0.55), plus U and X
    Y = 0.0 +
      0.55 * D +                        # TRUE causal effect we hope IV recovers
      0.30 * log1p(prior_sessions) +
      0.12 * device_mobile +
      0.35 * U +                        # confounding: affects outcome too
      0.15 * country_tier_num +
      rnorm(n, sd = 1.0)
  ) %>%
  # Drop U to mimic real-world unobserved confounding
  select(-U, -country_tier_num)

glimpse(df)

# 2. Quick EDA Plots help to see whether D has reasonable support and whether Z shifts D.
ggplot(df, aes(x = D)) +
  geom_histogram(bins = 50) +
  labs(title = "Treatment (D) distribution", 
       x = "D", 
       y = "Count")+
  theme_fancy()

ggplot(df, aes(x = factor(Z), y = D)) +
  geom_boxplot() +
  labs(title = "First-stage relevance check: D by instrument Z", 
       x = "Z", 
       y = "D")+
  theme_fancy()

ggplot(df, aes(x = D, y = Y)) +
  geom_point(alpha = 0.12) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(title = "Naive relationship: Y vs D (confounded)", 
       x = "D", 
       y = "Y")+
  theme_fancy()

# 3. DAG 
# Z instruments D (Z -> D)
# D causally affects Y (D -> Y)
# U confounds D and Y (U -> D and U -> Y)
# X are observed covariates that affect D and Y (X -> D, X -> Y)

# IV assumptions:
# Exclusion restriction: no direct arrow Z -> Y (Z affects Y only through D).
# Independence: Z is as-good-as random w.r.t. U (no arrow U -> Z, no common causes).
# In practice, you’ll justify Z carefully (rollout, threshold, randomness, etc.)

dag <- dagitty("dag {
  Z -> D
  D -> Y
  X -> D
  X -> Y
  U -> D
  U -> Y
}")

# Add X as a single node standing in for observed covariates (prior_sessions, device, country)
coordinates(dag) <- list(
  x = c(Z = 0, X = 0, U = 1, D = 1, Y = 2),
  y = c(Z = 1, X = 0, U = 1, D = 0.5, Y = 0.5)
)

ggdag(dag, text = FALSE, use_labels = "name") +
  theme_dag() +
  ggtitle("IV DAG: Z instruments D; U confounds D and Y; X observed covariates")

# 4. Baseline OLS 
# Why include OLS at all?
# Default: regress outcome on treatment + covariates.
# If D is endogenous (correlated with omitted causes of Y), OLS is biased/inconsistent.
# Comparing OLS vs IV highlights the magnitude/direction of endogeneity bias.

# Use log1p(prior_sessions) instead of raw prior_sessions b/c prior_sessions is a count 
# feature often right-skewed with diminishing marginal impact. log1p(x) = log(1 + x) is 
# defined at x=0 (unlike log(x)) and compresses large values, reducing leverage/outlier 
#influence and making a linear model a better approximation.

ols <- lm(Y ~ D + log1p(prior_sessions) + device_mobile + country_tier, data = df)

# Robust SE (HC2 is a decent default; cluster at user_id not meaningful here, 
# but in real data you may cluster)/ Product data often violate homoskedasticity 
# (constant variance of errors). If heteroskedasticity exists, OLS coefficients 
# are still unbiased (under exogeneity), but conventional SEs/p-values can be wrong.

# HC0/HC1/HC2/HC3 are heteroskedasticity-consistent covariance estimators. HC2 adjusts 
# for leverage (influence of high-leverage points) by inflating variance more for points 
# with high hat-values. Often a good default in moderate samples.

# vcovHC() computes a robust variance-covariance matrix, coeftest() combines the model’s 
# coefficient estimates with that vcov to produce robust SEs, t-stats, and p-values. 
# Point estimates are unchanged—only inference changes.

ols_robust <- lmtest::coeftest(ols, vcov. = sandwich::vcovHC(ols, type = "HC2"))

# summary() isn't typically meaningful bc ols_robust is a matrix
print(ols_robust)

# 5. IV / 2SLS with AER::ivreg
# 2SLS (two-stage least squares
# Stage 1: predict D using Z and covariates (D_hat = E[D | Z, X]).
# Stage 2: regress Y on D_hat and covariates, using only the Z-induced variation in D.
# Under valid IV assumptions, the coefficient on D is causal (often a LATE in imperfect compliance).
# D is endogenous; Z is the instrument

iv <- AER::ivreg(
  Y ~ D + log1p(prior_sessions) + device_mobile + country_tier |
    Z + log1p(prior_sessions) + device_mobile + country_tier,
  data = df
)

# Robust SE for IV
iv_robust <- lmtest::coeftest(iv, 
                              vcov. = sandwich::vcovHC(iv, 
                                                       type = "HC2"))

print(iv_robust)

summary(iv_robust)

# 6. Common IV Diagnostics
# (A) Weak instruments test (first-stage strength; often an F-statistic / “weak IV” diagnostic)
#  - Checks whether Z meaningfully predicts D conditional on X.
#  - Weak instruments can make IV biased (toward OLS) and very noisy.
#  - Rule-of-thumb is a “first-stage F” comfortably above ~10 (context-dependent; stronger is better).
#  - Interpretation: if weak, consider a better instrument, redesign, or alternative estimand.
# (B) Wu–Hausman endogeneity test (a.k.a. Durbin–Wu–Hausman in some outputs)
#  - Tests whether OLS and IV estimates differ more than expected by chance.
#  - Null hypothesis: D is exogenous (OLS consistent); IV not needed for consistency.
#  - Small p-value (e.g., < 0.05) suggests endogeneity => prefer IV.
#  - If instrument is weak, this test can have low power / be unreliable.
# (C) Sargan (overidentification) test — ONLY when overidentified (>1 instrument)
#  - Tests whether instruments are jointly valid (uncorrelated with error term),
#    leveraging the fact you have “extra” instruments beyond what’s needed.
#  - Null hypothesis: instruments are valid (satisfy exclusion restrictions).
#  - Large p-value is “good” (fail to reject validity); small p-value
#    suggests at least one instrument violates exclusion.
#  - With only 1 instrument (exactly identified), Sargan is not available.

summary(iv, diagnostics = TRUE)

# Extract first stage explicitly (useful for reporting)
# In this simple setup, D is continuous; for binary D you’d consider LATE interpretation & alternatives.
first_stage <- lm(D ~ Z + log1p(prior_sessions) + device_mobile + country_tier, data = df)

fs_robust <- lmtest::coeftest(first_stage, vcov. = sandwich::vcovHC(first_stage, type = "HC2"))

print(fs_robust)

summary(fs_robust)

# First-stage partial F for Z (roughly): compare model with vs without Z
# Restricted model: D ~ X (no Z)
# Full model: D ~ Z + X
# F-stat tests whether adding Z significantly improves prediction of D (conditional on X).
# A large F-statistic and small p-value implies Z is relevant (good), if F is small (weak), 
# IV estimates can be unstable and biased; consider redesign or alternatives.
first_stage_restricted <- lm(D ~ log1p(prior_sessions) + device_mobile + country_tier, data = df)

fs_anova <- anova(first_stage_restricted, first_stage)

print(fs_anova)

# 7. Coefficient Comparison Table
# We filter to term == "D" to focus on the treatment effect and compare OLS vs IV side-by-side.
est_tbl <- bind_rows(
  tidy(ols) %>% mutate(model = "OLS"),
  tidy(iv)  %>% mutate(model = "IV (2SLS)")) %>%
  filter(term == "D") %>%
  select(model, term, estimate, std.error, statistic, p.value)

print(est_tbl)

# 8. Plot OLS vs IV estimate for D 
# (We’ll use robust SEs we computed above for a fair comparison.)

get_coef_ci <- function(ct, term = "D") {
  # ct is a matrix from coeftest
  est <- ct[term, "Estimate"]
  se  <- ct[term, "Std. Error"]
  tibble(
    term = term,
    estimate = est,
    conf.low = est - 1.96 * se,
    conf.high = est + 1.96 * se
  )
}

coef_plot_df <- bind_rows(
  get_coef_ci(ols_robust) %>% mutate(model = "OLS (robust)"),
  get_coef_ci(iv_robust)  %>% mutate(model = "IV (2SLS, robust)")
)

ggplot(coef_plot_df, aes(x = model, y = estimate)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.12) +
  coord_flip() +
  labs(title = "Estimated causal effect of D on Y",
    x = NULL,
    y = "Estimate (95% CI)")+
  theme_fancy()

# Interpretation guide 
# - If OLS and IV differ meaningfully and Wu–Hausman p < 0.05, evidence D is endogenous; 
#   IV is usually preferred (if instrument is strong/credible).
# - If first-stage F is weak (rule-of-thumb ~<10), IV may be unreliable; redesign instrument 
#.  or consider alternative approaches.
# - If overidentified and Sargan p is small at least one instrument likely violates exclusion; revisit instrument set.



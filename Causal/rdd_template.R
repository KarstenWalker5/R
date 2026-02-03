# Regression Discontinuity (RD) Boilerplate in R
# - Synthetic "Quizlet-like" panel data generation
# - Estimates the RD effect using rdrobust (local polynomials with robust bias-corrected inference).
# - Sharp RD (threshold-based feature unlock)
# - Runs core RD diagnostics:
#     (1) McCrary-style density test (rddensity) for manipulation/bunching at cutoff
#     (2) Covariate balance checks (placebo outcomes) around the cutoff
# - Demonstrates bandwidth sensitivity and a panel outcome extension.
# - Optional: "propensity score threshold" RD use case
# If units cannot precisely manipulate the running variable around the cutoff, then observations
#   just above vs just below the cutoff are “as-if randomized,” and any jump in the outcome at
#   the cutoff can be interpreted causally.
# Based on this article: https://arelbundock.com/posts/panel_heterogeneous_effects/index.html

# 1. Load packages and set theme
library(tidyverse)
library(fixest)
library(rdrobust)
library(rddensity)
library(broom)

set.seed(123)

theme_fancy <- function() {
  theme_minimal(base_family = "Asap Condensed") +
    theme(panel.grid.minor = element_blank(),
          legend.position = "bottom")
}

# 2. Synthetic data generation
# Core idea:
# Users have a baseline "study intent" score near a threshold (running variable)
# If they cross a threshold (e.g., "Onboarding mastery score >= 70"),
# they unlock a feature (e.g., "Smart Review"), causing a jump in engagement.

# We'll generate panel data: user i observed weekly t = 1..T
# Outcome: weekly_study_minutes (or sessions). Treatment: feature_unlocked.
# Running variable is measured at/near onboarding, with small noise.
# Also generate a second outcome: cumulative LTV-ish over horizon with a jump at cutoff.

simulate_quizlet_rd_panel <- function(
    n_users = 8000,
    n_weeks = 12,
    cutoff = 70,
    running_mean = 70,
    running_sd = 12,
    # treatment effect at cutoff (jump)
    tau_minutes = 12,
    # slope (smooth) of running variable effect
    beta_x = 0.6,
    # time trend (learning habit)
    beta_t = 2.0,
    # user heterogeneity
    user_sd = 18,
    # idiosyncratic noise
    eps_sd = 25
) {
  
  users <- tibble(
    user_id = 1:n_users,
    # "Onboarding mastery score" is the running variable (forcing variable)
    x_raw = rnorm(n_users, mean = running_mean, sd = running_sd),
    # some covariates for balance checks
    prior_experience = rbinom(n_users, 1, 0.35),
    device_mobile = rbinom(n_users, 1, 0.65),
    # user random effect
    u_i = rnorm(n_users, 0, user_sd)
  ) %>%
    mutate(
      # Treatment assignment: Sharp RD at cutoff
      treat = as.integer(x_raw >= cutoff),
      # Centered running variable (for nicer modeling/plots)
      x_c = x_raw - cutoff
    )
  
  panel <- tidyr::expand_grid(
    user_id = users$user_id,
    week = 1:n_weeks
  ) %>%
    left_join(users, by = "user_id") %>%
    mutate(
      # smooth seasonal / habit component
      season = 6 * sin(2 * pi * week / n_weeks),
      # baseline smooth outcome: depends on running variable smoothly + time
      mu = 50 +
        beta_x * x_c +
        beta_t * week +
        8 * prior_experience +
        4 * device_mobile +
        season +
        u_i,
      # add discontinuity: feature unlock increases weekly minutes
      weekly_study_minutes = pmax(
        0,
        mu + tau_minutes * treat + rnorm(n(), 0, eps_sd)
      ),
      # A binary “active” engagement measure (another common outcome)
      active_week = rbinom(
        n(),
        size = 1,
        prob = plogis(-0.6 + 0.015 * weekly_study_minutes)
      )
    )
  
  # user-level LTV-ish outcome over horizon
  # Here: revenue relates to engagement, plus a discrete jump for treated.
  user_ltv <- panel %>%
    group_by(user_id) %>%
    summarise(
      x_raw = first(x_raw),
      x_c = first(x_c),
      treat = first(treat),
      prior_experience = first(prior_experience),
      device_mobile = first(device_mobile),
      total_minutes = sum(weekly_study_minutes),
      weeks_active = sum(active_week),
      .groups = "drop"
    ) %>%
    mutate(
      # Map engagement to $ with noise, plus a discontinuity at cutoff
      ltv_12w = pmax(
        0,
        8 + 0.18 * total_minutes + 3.5 * weeks_active +
          10 * prior_experience +
          5 * device_mobile +
          25 * treat +                  # discontinuity in LTV at cutoff
          rnorm(n(), 0, 35)
      )
    )
  
  list(panel = panel, user = user_ltv, cutoff = cutoff)
}

sim <- simulate_quizlet_rd_panel()

panel_df <- sim$panel

user_df  <- sim$user

cutoff   <- sim$cutoff

# Helper: clean RD binned plot
rd_binned_plot <- function(df, x, y, cutoff, 
                           binwidth = 2, 
                           xlim = c(-40, 40), 
                           ylab = NULL, title = NULL) {
  # x should be centered at cutoff if you pass x = x_c; else pass x_raw and cutoff will be subtracted.
  df_plot <- df %>%
    mutate(xc = {{ x }}) %>%
    filter(xc >= xlim[1], xc <= xlim[2]) %>%
    mutate(bin = floor(xc / binwidth) * binwidth)
  
  binned <- df_plot %>%
    group_by(bin) %>%
    summarise(
      x_bin = mean(xc),
      y_bar = mean({{ y }}, na.rm = TRUE),
      n = n(),
      .groups = "drop"
    )
  
  ggplot() +
    geom_point(
      data = binned,
      aes(x = x_bin, y = y_bar, size = n),
      alpha = 0.75
    ) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    scale_size_continuous(range = c(1.5, 6)) +
    labs(
      x = "Running variable (centered at cutoff)",
      y = ylab %||% rlang::as_label(rlang::enquo(y)),
      title = title
    ) +
    theme_minimal()
}

# Can also load helper_function
source("/Users/karstenwalker/Documents/R/Causal/tidy_rd.R")

# 3. Visualize discontinuity: user-level LTV (cross-sectional RD)
# RD is about a local jump at the cutoff, so you always want to see:
# 1. Whether the conditional mean looks smooth away from the cutoff
# 2. Whether there’s an apparent jump at x=0.

p_ltv <- rd_binned_plot(
  df = user_df,
  x = x_c,
  y = ltv_12w,
  cutoff = CUTOFF,
  binwidth = 2,
  xlim = c(-40, 40),
  ylab = "12-week LTV (synthetic $)",
  title = "RD intuition plot: LTV vs. onboarding score threshold"
)

# Add smooth fits on each side (local polynomial via geom_smooth)
p_ltv_smooth <- user_df %>%
  filter(x_c >= -40, x_c <= 40) %>%
  ggplot(aes(x = x_c, y = ltv_12w)) +
  geom_point(alpha = 0.12) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_smooth(
    data = ~ dplyr::filter(.x, x_c < 0),
    method = "lm",
    formula = y ~ poly(x, 2, raw = TRUE),
    se = FALSE
  ) +
  geom_smooth(
    data = ~ dplyr::filter(.x, x_c >= 0),
    method = "lm",
    formula = y ~ poly(x, 2, raw = TRUE),
    se = FALSE
  ) +
  labs(
    x = "Running variable (centered at cutoff)",
    y = "12-week LTV (synthetic $)",
    title = "RD plot with separate polynomial fits (visual only)"
  ) +
  theme_fancy()

p_lt

p_ltv_smooth

# 4. Main RD estimation
# Outcome: ltv_12w
# Running variable: x_raw (or x_c; rdrobust uses the raw forcing variable; cutoff specified)
# Estimates local linear RD, default MSE-optimal bandwidth, robust bias-corrected inference
# rdrobust implements modern RD best practices: data-driven bandwidth selection (MSE-optimal / bias-aware)
# and robust bias-corrected inference (more reliable coverage near the boundary)

# Sharp RD vs Fuzzy RD: This script is SHARP RD (treat = 1[x >= c]) with deterministic assignment.
# In real product settings, you sometimes have FUZZY RD (probability of treatment jumps),
# which uses the cutoff as an instrument for actual uptake.

rd_ltv <- rdrobust(
  y = user_df$ltv_12w,
  x = user_df$x_raw,
  c = CUTOFF,
  p = 1 # local linear
)

print(rd_ltv)

# Nice printable summary table
summary(rd_ltv)

# Helper: safely pull a value from a matrix/data.frame by column name, else NA
as_tidy_rd <- function(fit) {
  stopifnot(inherits(fit, "rdrobust"))
  
  # rdrobust always returns:
  # - Estimate: conventional + robust
  # - se:       same ordering
  # - ci:       lower/upper bounds
  # - bws:      left/right bandwidths
  
  est  <- fit$Estimate
  se   <- fit$se
  ci   <- fit$ci
  bws  <- fit$bws
  
  # Standardize positions (this is stable across rdrobust versions)
  tibble(
    estimand = c("conventional", "robust"),
    tau      = as.numeric(est[1, ]),
    se       = as.numeric(se[1, ]),
    ci_lo    = as.numeric(ci[1, c(1, 3)]),
    ci_hi    = as.numeric(ci[1, c(2, 4)]),
    h_left   = bws[1, 1],
    h_right  = bws[1, 2]
  )
}

# Can also use helper_function
rd_ltv_tidy <- tidy_rd(rd_ltv)

rd_ltv_tidy

# 5. Bandwidth sensitivity: grid of symmetric bandwidths and plot estimates.
# RD is local: results can depend on the bandwidth (how close to the cutoff you look).
# A good RD result is reasonably stable across a plausible range of bandwidths.
# If τ changes sign or swings wildly as h varies, the estimate may be fragile.
# If τ is stable but SE grows as h shrinks, that’s normal (less data -> more variance).

bw_grid <- c(5, 7.5, 10, 12.5, 15, 20, 25)

sens <- map_dfr(bw_grid, function(hh) {
  fit <- rdrobust(
    y = user_df$ltv_12w,
    x = user_df$x_raw,
    c = CUTOFF,
    p = 1,
    h = c(hh, hh)
  )
  
# Robust/RBC is generally the inferential target you report.
  tidy_rd(fit) %>%
    # keep a single row per fit (your printed output has "RD Effect (robust)")
    filter(grepl("RD Effect", estimand, ignore.case = TRUE) |
             grepl("robust|rb", estimand, ignore.case = TRUE)) %>%
    slice(1) %>%               # in case multiple rows match
    transmute(
      h = hh,
      tau_rb = tau,
      se_rb  = se,
      ci_lo  = ci_lo,
      ci_hi  = ci_hi
    )
})

sens

sens %>%
  ggplot(aes(x = h, y = tau_rb)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    x = "Symmetric bandwidth (h)",
    y = "RD estimate (robust bias-corrected)",
    title = "Bandwidth sensitivity: RD effect on LTV"
  ) +
  theme_fancy()

# 6. Diagnostics
# RD is credible when units cannot precisely manipulate x around the cutoff and potential outcomes
# are smooth in x at the cutoff absent treatment.
# This section implements two standard checks:
# 1. Density (manipulation) test
# 2. Covariate balance “placebo outcomes”

# McCrary-style density test around cutoff (manipulation check) tests whether the density of the 
# running variable jumps at the cutoff. A discontinuity in density suggests sorting/manipulation 
# near the threshold, which breaks the “as-if randomized near cutoff” intuition.
# Null hypothesis: density is continuous at cutoff (NO manipulation).
# If the test returns a small p-value (e.g., < 0.05), that’s evidence of bunching/manipulation.

dens <- rddensity(X = user_df$x_raw, c = CUTOFF)

print(summary(dens))

# Plot density 
user_df %>%
  ggplot(aes(x = x_raw)) +
  geom_histogram(bins = 60) +
  geom_vline(xintercept = CUTOFF, linetype = "dashed") +
  labs(
    x = "Running variable (onboarding score)",
    y = "Count",
    title = "Density around cutoff (visual manipulation check)"
  ) +
  theme_fancy()

# Covariate balance around cutoff (placebo outcomes)
# If RD is valid, predetermined covariates should not jump at cutoff.
# If covariates jump, it suggests non-comparability around cutoff (sorting or measurement artifacts).
# Look for estimated discontinuity ~ 0, with CIs covering 0, for each predetermined covariate.
balance_covariates <- function(df, covar, cutoff, bw = 15) {
  sub <- df %>%
    dplyr::filter(abs(x_raw - cutoff) <= bw)
  
  fit <- rdrobust(
    y = sub[[covar]],
    x = sub$x_raw,
    c = cutoff,
    p = 1
  )
  
  # Use tidy_rd() instead of indexing fit$Estimate / fit$se / fit$ci
  tr <- tidy_rd(fit) %>%
    dplyr::filter(
      grepl("RD Effect", estimand, ignore.case = TRUE) |
        grepl("robust|rb", estimand, ignore.case = TRUE)
    ) %>%
    dplyr::slice(1)
  
  tibble::tibble(
    covariate = covar,
    bw = bw,
    tau_rb = tr$tau,
    se_rb  = tr$se,
    ci_lo  = tr$ci_lo,
    ci_hi  = tr$ci_hi
  )
}

# Can also user helper_function
bal <- dplyr::bind_rows(
  balance_covariates(user_df, "prior_experience", CUTOFF, bw = 15),
  balance_covariates(user_df, "device_mobile",   CUTOFF, bw = 15)
)

bal

bal %>%
  ggplot(aes(x = covariate, y = tau_rb)) +
  geom_point() +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    x = NULL, y = "Estimated discontinuity (should be ~0)",
    title = "Covariate balance around cutoff (RD placebo checks)"
  ) +
  theme_fancy()

# 7. Panel outcome example
# RD is about the discontinuity at assignment threshold.
# With panel outcomes, two clean approaches are:
# 1. Collapse to a user-level summary over a window (e.g., weeks 1–4) and run RD once.
# 2. Estimate RD separately by time period to show dynamics (effect over time).
# Why not just run a big panel regression with an interaction?
# - You can run a big panel regression with interacton, but rdrobust is designed for 
# cross-sectional RD; collapsing or week-by-week RD is often simpler and aligns with standard RD inference.

user_weekly <- panel_df %>%
  group_by(user_id, week) %>%
  summarise(
    x_raw = first(x_raw),
    x_c = first(x_c),
    treat = first(treat),
    weekly_study_minutes = mean(weekly_study_minutes),
    .groups = "drop"
  )

# 7a. Plot discontinuity in week 2 (or any week)
user_weekly %>%
  filter(week == 2) %>%
  rd_binned_plot(
    x = x_c, y = weekly_study_minutes,
    cutoff = CUTOFF, binwidth = 2, xlim = c(-40, 40),
    ylab = paste0("Weekly study minutes (week ", week_k, ")"),
    title = paste0("RD plot: engagement jump in week ", week_k)
  )

# 7b. Aggregate weeks 1..4 and run RD
user_w14 <- panel_df %>%
  filter(week %in% 1:4) %>%
  group_by(user_id) %>%
  summarise(
    x_raw = first(x_raw),
    x_c = first(x_c),
    treat = first(treat),
    minutes_w14 = mean(weekly_study_minutes),
    .groups = "drop"
  )

rd_binned_plot(
  df = user_w14,
  x = x_c,
  y = minutes_w14,
  cutoff = CUTOFF,
  binwidth = 2,
  xlim = c(-40, 40),
  ylab = "Mean weekly study minutes (weeks 1–4)",
  title = "RD intuition plot: early engagement (weeks 1–4)"
)

rdrobust(
  y = user_w14$minutes_w14,
  x = user_w14$x_raw,
  c = CUTOFF,
  p = 1
)

# 7c. Effect over time (estimate RD separately by week, then plot)
# Look for whether the effect is immediate vs delayed, and whether it decays.
# If effects are only present in some weeks, that’s still causal—just time-varying.
rd_by_week <- map_dfr(sort(unique(user_weekly$week)), function(w) {
  sub <- user_weekly %>% filter(week == w)
  
  tr <- tidy_rd(rdrobust(y = sub$weekly_study_minutes, x = sub$x_raw, c = CUTOFF, p = 1)) %>%
    filter(grepl("RD Effect", estimand, ignore.case = TRUE) |
             grepl("robust|rb", estimand, ignore.case = TRUE)) %>%
    slice(1)
  
  tibble(
    week = w,
    tau_rb = tr$tau,
    se_rb  = tr$se,
    ci_lo  = tr$ci_lo,
    ci_hi  = tr$ci_hi
  )
})

rd_by_week %>%
  ggplot(aes(x = week, y = tau_rb)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    x = "Week",
    y = "RD estimate (weekly minutes jump at cutoff)",
    title = "RD effect over time (estimated separately each week)"
  ) +
  theme_fancy()

# Optional use case: “propensity score threshold” RD (template)
# Caution: In real work, RD on an estimated propensity score is tricky:
#  1. Estimation error in the running variable can bias inference
#  2. Often better framed as a policy threshold on a *score* used operationally
# This template is if you *do* have a deterministic score used for gating.

# Caution 
# RD on an estimated propensity score is generally problematic because the running variable
# includes estimation error (generated regressor), complicating inference.
# RD is best justified when the *operational score* used for gating is deterministic and
# policy-relevant (e.g., “if score >= 0.65, unlock X”), not when you fit a score ad-hoc.

simulate_propensity_threshold <- function(n = 10000, cutoff = 0.65) {
  df <- tibble(
    user_id = 1:n,
    age = pmax(13, round(rnorm(n, 22, 6))),
    prior_sessions = rpois(n, 5),
    marketing_touch = rbinom(n, 1, 0.4)
  ) %>%
    mutate(
      # "Operational score" (think: deterministic scoring model output used for eligibility)
      score = plogis(-2.0 + 0.06 * (age - 18) + 0.12 * prior_sessions + 0.5 * marketing_touch + rnorm(n, 0, 0.2)),
      treat = as.integer(score >= cutoff),
      score_c = score - cutoff,
      # LTV jumps at threshold + smooth relationship with score
      ltv = pmax(0, 40 + 30 * score_c + 55 * treat + rnorm(n, 0, 35))
    )
  df
}

ps_df <- simulate_propensity_threshold()

rd_binned_plot(
  df = ps_df,
  x = score_c,
  y = ltv,
  cutoff = 0.65,
  binwidth = 0.02,
  xlim = c(-0.25, 0.25),
  ylab = "LTV ($)",
  title = "RD intuition plot: LTV vs score threshold (template)"
)

rd_ps <- rdrobust(y = ps_df$ltv,
                  x = ps_df$score,
                  c = 0.65, 
                  p = 1)

rd_ps

######
# RD Simulation: "Pass a Mastery Exam" with repeated attempts, optional manipulation/bunching near cutoff,
# and shows how rddensity (McCrary-style) "lights up"

# Assumptions:
# - Each user has latent ability a_i
# - They can attempt the mastery exam multiple times (up to K)
# - Each attempt yields a noisy score
# - The recorded running variable is the FIRST passing score (or last score if never pass)
# - Treatment is "feature unlocked" if recorded_score >= cutoff (sharp RD on recorded_score)

# Manipulation mechanism:
# - Some users with scores just below cutoff are more likely to re-take
# - This creates bunching just above cutoff in the RECORDED score distribution
#   (even if each attempt score is not manipulable directly)
######

############################################################
# RD Simulation: "Pass a Mastery Exam" with repeated attempts
# + optional manipulation/bunching near cutoff
# + show how rddensity (McCrary-style) "lights up"
#
# What this sim does:
# - Each user has latent ability a_i
# - They can attempt the mastery exam multiple times (up to K)
# - Each attempt yields a noisy score
# - The recorded running variable is the FIRST passing score (or last score if never pass)
# - Treatment is "feature unlocked" if recorded_score >= cutoff (sharp RD on recorded_score)
#
# Manipulation mechanism (optional):
# - Some users with scores just below cutoff are more likely to re-take
# - This creates bunching just above cutoff in the RECORDED score distribution
#   (even if each attempt score is not manipulable directly)
############################################################

simulate_mastery_exam_rd <- function(
    n_users = 12000,
    cutoff = 70,
    max_attempts = 5,
    # ability distribution
    ability_mean = 0,
    ability_sd = 1,
    # score model per attempt
    base_score = 62,
    ability_to_score = 10,
    attempt_learning = 1.5,   # each additional attempt slightly increases expected score
    score_sd = 7,
    # outcome model after unlock
    tau_minutes = 12,
    beta_score_smooth = 0.5,
    eps_sd = 25,
    # manipulation controls
    manipulation = FALSE,
    # if manipulation=TRUE: probability to retake increases when near cutoff and below
    retake_base = 0.25,         # baseline chance to retake after a failed attempt
    retake_near_boost = 0.55,   # extra boost when close to cutoff
    near_window = 4,            # "close to cutoff" window (points below cutoff)
    # optional: discourage retakes if far below cutoff
    retake_far_shrink = 0.10
) {
  users <- tibble(
    user_id = 1:n_users,
    ability = rnorm(n_users, ability_mean, ability_sd),
    prior_experience = rbinom(n_users, 1, 0.35),
    device_mobile = rbinom(n_users, 1, 0.65)
  )
  
  # Simulate attempt-by-attempt scores and stopping/retake behavior
  attempt_tbl <- map_dfr(users$user_id, function(i) {
    a <- users$ability[users$user_id == i]
    
    scores <- numeric(max_attempts)
    took   <- logical(max_attempts)
    stop_after <- max_attempts
    
    for (k in 1:max_attempts) {
      took[k] <- TRUE
      mu_k <- base_score + ability_to_score * a + attempt_learning * (k - 1)
      scores[k] <- rnorm(1, mean = mu_k, sd = score_sd)
      
      # if pass, stop (feature unlock)
      if (scores[k] >= cutoff) {
        stop_after <- k
        break
      }
      
      # decide whether to retake (only if not last attempt)
      if (k < max_attempts) {
        if (!manipulation) {
          retake_p <- retake_base
        } else {
          # manipulation: more likely to retake when narrowly below cutoff
          gap <- cutoff - scores[k]  # positive if below cutoff
          near_below <- (gap > 0) && (gap <= near_window)
          far_below  <- (gap > near_window)
          
          retake_p <- retake_base +
            (if (near_below) retake_near_boost else 0) -
            (if (far_below)  retake_far_shrink else 0)
          
          retake_p <- pmin(pmax(retake_p, 0.01), 0.99)
        }
        
        if (runif(1) > retake_p) {
          stop_after <- k
          break
        }
      }
    }
    
    tibble(
      user_id = i,
      attempt = 1:stop_after,
      score = scores[1:stop_after],
      passed = score >= cutoff
    )
  })
  
  # Recorded running variable:
  # - first passing score if ever passed
  # - else last observed score
  recorded <- attempt_tbl %>%
    group_by(user_id) %>%
    summarise(
      attempts = n(),
      pass_any = any(passed),
      score_recorded = if (any(passed)) score[which(passed)[1]] else score[n()],
      .groups = "drop"
    ) %>%
    left_join(users, by = "user_id") %>%
    mutate(
      x_raw = score_recorded,
      x_c = x_raw - cutoff,
      treat = as.integer(x_raw >= cutoff)
    )
  
  # Outcome: weekly study minutes in the first week after exam
  # (you can extend this to a panel easily)
  outcomes <- recorded %>%
    mutate(
      mu = 50 +
        beta_score_smooth * x_c +
        8 * prior_experience +
        4 * device_mobile,
      weekly_study_minutes = pmax(0, mu + tau_minutes * treat + rnorm(n(), 0, eps_sd))
    )
  
  list(
    users = users,
    attempts = attempt_tbl,
    df = outcomes,
    cutoff = cutoff
  )
}

# 1. Simulate: no manipulation vs manipulation
sim_nom <- simulate_mastery_exam_rd(manipulation = FALSE)

sim_man <- simulate_mastery_exam_rd(manipulation = TRUE)

df_nom <- sim_nom$df

df_man <- sim_man$df

cutoff <- sim_nom$cutoff

# 2. Visualize density around cutoff 

plot_density_around_cutoff <- function(df, cutoff, window = 25, bins = 60, title = NULL) {
  df %>%
    filter(x_raw >= cutoff - window, x_raw <= cutoff + window) %>%
    ggplot(aes(x = x_raw)) +
    geom_histogram(bins = bins) +
    geom_vline(xintercept = cutoff, linetype = "dashed") +
    labs(
      x = "Recorded mastery exam score (running variable)",
      y = "Count",
      title = title
    ) +
    theme_minimal()
}

p_dens_nom <- plot_density_around_cutoff(
  df_nom, cutoff,
  title = "No manipulation: recorded score density near cutoff"
)

p_dens_man <- plot_density_around_cutoff(
  df_man, cutoff,
  title = "Manipulation via retakes: bunching in recorded scores near cutoff"
)

p1<-p_dens_nom

p2<-p_dens_man

library(patchwork)

p1/p2

# 3 Density test: rddensity (McCrary-style)

dens_nom <- rddensity(X = df_nom$x_raw, c = cutoff)

dens_man <- rddensity(X = df_man$x_raw, c = cutoff)

cat("\n--- rddensity summary: no manipulation ---\n")
print(summary(dens_nom))

cat("\n--- rddensity summary: manipulation via retakes ---\n")
print(summary(dens_man))

# 4. RD estimation (to see the effect still estimates, but diagnostics differ)

rd_nom <- rdrobust(y = df_nom$weekly_study_minutes, 
                   x = df_nom$x_raw, 
                   c = cutoff, 
                   p = 1)

rd_man <- rdrobust(y = df_man$weekly_study_minutes,
                   x = df_man$x_raw,
                   c = cutoff, p = 1)

cat("\n--- RD summary: no manipulation ---\n")
print(summary(rd_nom))

cat("\n--- RD summary: manipulation via retakes ---\n")
print(summary(rd_man))

# 5. RD binned plot for outcome 

p3<-rd_binned_plot(
  df = df_nom%>%
    mutate(x_c = x_raw - cutoff),
  x = x_c,
  y = weekly_study_minutes,
  cutoff = cutoff,           
  binwidth = 1,
  xlim = c(-25, 25),
  ylab = "Weekly study minutes",
  title = "No manipulation: RD binned outcome plot"
  )

p4<-rd_binned_plot(
  df = df_man%>%
    mutate(x_c = x_raw - cutoff),
  x = x_c,
  y = weekly_study_minutes,
  cutoff = cutoff,
  binwidth = 1,
  xlim = c(-25, 25),
  ylab = "Weekly study minutes",
  title = "Manipulation via retakes: RD binned outcome plot"
  )

p3/p4

# 6. Show attempt counts vs score: how retakes drive bunching
# ---------------------------

plot_attempts_vs_score <- function(df, cutoff, window = 25, binwidth = 1, title = NULL) {
  df %>%
    filter(x_raw >= cutoff - window, x_raw <= cutoff + window) %>%
    mutate(x_c = x_raw - cutoff,
           bin = floor(x_c / binwidth) * binwidth) %>%
    group_by(bin) %>%
    summarise(
      x_bin = mean(x_c),
      avg_attempts = mean(attempts),
      n = n(),
      .groups = "drop"
    ) %>%
    ggplot(aes(x = x_bin, y = avg_attempts)) +
    geom_line() +
    geom_point() +
    geom_vline(xintercept = 0, linetype = "dashed") +
    labs(
      x = "Score (centered at cutoff)",
      y = "Average # of attempts",
      title = title
    ) +
    theme_fancy()
}

p5<-plot_attempts_vs_score(df_nom, 
                           cutoff, 
                           title = "No manipulation: attempts vs score near cutoff")

p6<-plot_attempts_vs_score(df_man, 
                           cutoff, 
                           title = "Manipulation: attempts spike just below cutoff")

p5/p6

###### Purpose ######
# 1. Simulates a confounded dataset
# 2. Fits one propensity model (A ~ x1+x2+x3) and reuses it for both methods
# 3. Fits one propensity model (A ~ x1+x2+x3) and reuses it for both methods
# 4. IPTW (ATE) path
# 5. PSM (ATT) path
# 6. Results reporting including plots and tables

# Use when:
# 1. you need a population-level effect vs treated-group effect
# The code estimates:
#  1.IPTW (ATE) — “What would happen if everyone received the treatment vs no one?”
#  2. PSM (ATT) — “What happened to those who actually got treated compared to if they hadn’t?”
# Example:
# 1. If leadership wants the overall lift from rolling out a feature to all users → ATE.
# 2. If you want to understand how much benefit current adopters got → ATT.
#
# 2. When confounding is expected and must be adjusted
# 3. If treatment assignment is related to user characteristics (e.g., engaged users are more likely to adopt a product)
# This pipeline lets you estimate propensity scores based on pre-treatment covariates, check covariate balance with SMDs,
# and correct for differences with weighting (IPTW) or matching (PSM).

# ---- Packages ----
library(tidyverse)
library(MatchIt)      # PSM
library(cobalt)       # balance diagnostics (bal.tab)
library(survey)       # weighted estimation for IPTW
library(broom)        # tidy model outputs
library(knitr)
library(kableExtra)
library(ggplot2)
library(forcats)

set.seed(42)

# ---- 1) Simulate data with confounding ----
n  <- 4000
x1 <- rnorm(n)                      # continuous
x2 <- rbinom(n, 1, 0.4)             # binary
x3 <- rnorm(n, mean = x2, sd = 0.8) # correlated with x2

# Treatment assignment depends on covariates (confounding)
lin_ps   <- -0.3 + 0.8*x1 + 0.9*x2 - 0.5*x3
p_treat  <- plogis(lin_ps)
A        <- rbinom(n, 1, p_treat)

# Outcome depends on treatment + covariates + noise
tau <- 2.0                      # true treatment effect
Y   <- tau*A + 1.2*x1 + 0.6*x2 - 0.8*x3 + rnorm(n, sd = 1.5)

dat <- tibble(Y, A, x1, x2, x3)

# ---- 2) Propensity model (shared) ----
ps_mod <- glm(A ~ x1 + x2 + x3, data = dat, family = binomial())

dat    <- dat %>% mutate(ps = predict(ps_mod, type = "response"))

# ========================== IPTW (ATE) ==========================
# ---- 3) Stabilized IPTW weights ----
pA <- mean(dat$A)

dat <- dat %>%
  mutate(w_iptw = if_else(A == 1, pA/ps, (1-pA)/(1-ps)))

# ---- 4) Balance diagnostics (request numeric SMDs explicitly) ----
bal_iptw <- cobalt::bal.tab(
  A ~ x1 + x2 + x3,
  data      = dat,
  weights   = dat$w_iptw,
  method    = "weighting",
  estimand  = "ATE",
  un        = TRUE,
  stats     = "m",          # standardized mean differences (Diff.*)
  s.d.denom = "pooled"
)

# ---- 5) Effect estimation (ATE) with survey-weighted regression ----
des       <- svydesign(ids = ~1, data = dat, weights = ~w_iptw)

fit_iptw  <- svyglm(Y ~ A, design = des)

ate_tbl   <- tidy(fit_iptw, conf.int = TRUE) %>%
  filter(term == "A") %>%
  transmute(Method = "IPTW (ATE)",
            Estimate = estimate, SE = std.error,
            `95% CI Low` = conf.low, `95% CI High` = conf.high,
            `p-value` = p.value)

# Optional: weight diagnostics
w_summ <- dat %>% summarize(
  Method = "IPTW (stabilized)",
  Mean = mean(w_iptw), SD = sd(w_iptw),
  P1 = quantile(w_iptw, 0.01), P5 = quantile(w_iptw, 0.05),
  P50 = quantile(w_iptw, 0.50), P95 = quantile(w_iptw, 0.95),
  P99 = quantile(w_iptw, 0.99), Max = max(w_iptw)
)

# ========================== PSM (ATT) ==========================
# ---- 6) Nearest-neighbor matching on the same ps ----
m.out <- matchit(A ~ x1 + x2 + x3,
                 data = dat, distance = dat$ps,
                 method = "nearest", ratio = 1, estimand = "ATT")

# ---- 7) Balance diagnostics for PSM (numeric SMDs) ----
bal_psm <- cobalt::bal.tab(
  m.out,
  un        = TRUE,
  stats     = "m",
  s.d.denom = "pooled"
)

# ---- 8) Effect estimation (ATT) on matched sample ----
matched <- match.data(m.out)

att_fit <- lm(Y ~ A, data = matched, weights = weights)

att_tbl <- tidy(att_fit, conf.int = TRUE) %>%
  filter(term == "A") %>%
  transmute(Method = "PSM (ATT)",
            Estimate = estimate, SE = std.error,
            `95% CI Low` = conf.low, `95% CI High` = conf.high,
            `p-value` = p.value)

# =================== Minimal, robust ggplot Love Plot ===================
# Tidy bal.tab into a ggplot-friendly DF (SMDs only)
bal_to_df <- function(bt, method_label) {
  b <- bt$Balance |> as.data.frame() |> tibble::rownames_to_column("covariate")
  tibble(
    covariate  = b$covariate,
    Unadjusted = abs(as.numeric(b$Diff.Un)),
    Adjusted   = abs(as.numeric(b$Diff.Adj))
  ) |>
    tidyr::pivot_longer(c(Unadjusted, Adjusted),
                        names_to = "stage", values_to = "smd") |>
    dplyr::filter(!is.na(smd)) |>
    dplyr::mutate(method = method_label)
}

smd_iptw_df <- bal_to_df(bal_iptw, "IPTW (ATE)")

smd_psm_df  <- bal_to_df(bal_psm,  "PSM (ATT)")

love_df <- bind_rows(smd_iptw_df, smd_psm_df) %>%
  group_by(method, covariate) %>%
  mutate(order_key = max(ifelse(stage == "Unadjusted", smd, NA_real_), na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(covariate = fct_reorder(covariate, order_key, .desc = TRUE))

love_plot <- ggplot(love_df, aes(x = smd, y = covariate, color = stage, shape = stage)) +
  geom_point(size = 3) +
  geom_vline(xintercept = c(0.10, 0.20), linetype = "dashed") +
  facet_wrap(~ method, ncol = 2, scales = "free_y") +
  labs(title = "Covariate Balance (Absolute Standardized Mean Differences)",
       x = "Absolute SMD", y = "Covariate", color = "", shape = "") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

print(love_plot)

# =================== Printable Tables (Rmd-friendly) ===================
# SMD tables per method
smd_iptw_tbl <- smd_iptw_df %>%
  pivot_wider(names_from = stage, values_from = smd) %>%
  arrange(desc(Unadjusted))

smd_psm_tbl <- smd_psm_df %>%
  pivot_wider(names_from = stage, values_from = smd) %>%
  arrange(desc(Unadjusted))

# Effect estimates table
estimates_tbl <- bind_rows(ate_tbl, att_tbl)

# Matched sample summary
# ---- Replace your match_tbl block with this ----
m_summ <- summary(m.out)

nn <- m_summ$nn

stopifnot(is.matrix(nn))

# helper safely pulls by row/col name (returns 0 if the row doesn't exist)
pull_nn <- function(nn, r, c) {
  if (!(r %in% rownames(nn))) return(0)
  if (!(c %in% colnames(nn))) return(0)
  nn[r, c]
}

treated_matched <- pull_nn(nn, "Matched",   "Treated")
control_matched <- pull_nn(nn, "Matched",   "Control")
discard_ctrl    <- pull_nn(nn, "Unmatched", "Control")  # may be 0 if none

match_tbl <- tibble::tibble(
  Metric = c("Treated (matched)", "Control (matched)", "Discarded controls"),
  Value  = c(treated_matched,      control_matched,      discard_ctrl)
)


# Render with kableExtra
kbl(smd_iptw_tbl, digits = 3, caption = "IPTW (ATE): Absolute SMDs Before/After Weighting") %>%
  kable_minimal(full_width = FALSE)

kbl(smd_psm_tbl, digits = 3, caption = "PSM (ATT): Absolute SMDs Before/After Matching") %>%
  kable_minimal(full_width = FALSE)

kbl(estimates_tbl, digits = 3, caption = "Treatment Effect Estimates (ATE vs ATT)") %>%
  kable_minimal(full_width = FALSE)

kbl(w_summ, digits = 3, caption = "IPTW Weight Summary (Stabilized)") %>%
  kable_minimal(full_width = FALSE)

kbl(match_tbl, digits = 0, caption = "PSM Match Counts") %>%
  kable_minimal(full_width = FALSE)

# ----------------- Quick comparison summary (console) -----------------
cat("\n----- Quick Takeaways -----\n")
cat("IPTW (ATE) estimate:", round(coef(fit_iptw)['A'], 3), "\n")
cat("PSM (ATT) estimate :", round(coef(att_fit)['A'], 3),  "\n")

###### Extending for your own use ######
# 1. Swap the simulation with your data frame (say df), and rename: Outcome: Y → your metric (e.g., spend, MDU, CTR),
# Treatment: A → your binary exposure (owned device, saw feature). Covariates: add all pre-treatment confounders you trust.

ps_mod <- glm(treat ~ pre_var1 + pre_var2 + ... , data=df, family=binomial())
df$ps  <- predict(ps_mod, type="response")

# 2. Improve propensity estimation
# Use nonlinear learners for PS if needed (often improves balance): gbm, xgboost, or twang::ps(); 
# or WeightIt for flexible weighting.
# Always check overlap: plot PS histograms by treatment; poor overlap → trim or restrict common support.
# Stabilize / trim IPTW. If weights are heavy-tailed, trim at 1st–99th (or 5th–95th) percentiles:

q <- quantile(df$w_iptw, c(.01,.99))

df$w_iptw <- pmax(pmin(df$w_iptw, q[2]), q[1])

# Re-run bal.tab and re-estimate ATE.

# 3. Refine matching
# Add a caliper to reduce poor matches:

cal <- 0.2*sd(qlogis(df$ps))

m.out <- matchit(treat ~ covars..., data=df, distance=df$ps,
                 method="nearest", ratio=1, estimand="ATT", caliper=cal)

# Consider ratio > 1 (e.g., 1:2) for more power if good controls exist.

# Optional: Doubly-robust estimation (recommended). Add an outcome model on top of IPTW (AIPW) or use TMLE:
# Packages: AIPW, drtmle, tlverse/tmlenet, causalToolbox.
# Benefit: consistent if either PS or outcome model is correct.

# 5. Heterogeneity & reporting
# Explore CATEs (treatment effect by segment):

svyglm(Y ~ A*segment + ... , design=des)          # IPTW with interaction

lm(Y ~ A*segment, data=match.data(m.out), weights=weights)  # PSM

# Or use causal forests (grf) for nonparametric heterogeneity.

# 6. Robust inference & clustering
# If units repeat (e.g., user-month), cluster SEs by unit/geo:

sandwich::vcovCL(fit_iptw, cluster = df$user_id) |> lmtest::coeftest(fit_iptw, vcov=_)

# 7.Time-varying treatments
# If treatment changes over time, move to Marginal Structural Models with time-varying IPTW 
# (stabilized product of inverse probabilities).

# 8. RMarkdown outputs
# You already have kable tables. For executive decks, consider gt, flextable, or modelsummary to 
# standardize formatting across pages.

# 9. Guardrails & sensitivity
# Balance target: aim for |SMD| < 0.10 (often < 0.05 post-adjustment).
# Rosenbaum sensitivity (unobserved confounding): rbounds::psens (for matched data).

###### Quick checklist for your real analysis ######
# 1. Define estimand (ATE vs ATT) tied to the business question.
# 2. Specify pre-treatment covariates; confirm no leakage.
# 3. Fit PS; inspect overlap; adjust (trim/restrict) if needed.
# 4. Achieve balance (SMD < 0.1) via weighting or matching settings.
# 5. Estimate effect with appropriate SEs (survey design or matched weights).
# 6. Report ATE and/or ATT with CIs, plus love plot and diagnostics.
# 7. (Optional) Add AIPW/TMLE for robustness and CATEs for insights.
# 8. Placebos: test outcomes that shouldn’t move; test pre-period “effects” to catch residual confounding.
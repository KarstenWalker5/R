
# SYNTHETIC CONTROLS + TWFE : CODE EXAMPLES IN THIS SCRIPT

# 1) Synthetic Difference-in-Differences (SDiD)
#    Purpose:
#      - Formal hybrid of Synthetic Controls + TWFE
#      - Relaxes parallel trends via unit & time weighting
#    Use when:
#      - Treatment timing is common (or data can be subset)
#      - You want a clean ATE-style estimand
#      - Reviewers are familiar with SDiD
# 2) Synthetic Controls as preprocessing → DiD / TWFE-on-gap
#    Purpose:
#      - Classic SC intuition with simple inference
#      - Construct a synthetic counterfactual, then estimate post–pre gap
#      - Synthetic series extracted directly from dp$Y0plot %*% w
#    Use when:
#      - Single treated unit
#      - High narrative value of a “synthetic twin”
#      - Full transparency of counterfactual construction is important
# 3) SC-style preprocessing → TWFE (restricted-sample regression)
#    Purpose:
#      - Enforce overlap / approximate parallel trends, then run TWFE
#    Use when:
#      - Stakeholders want regression output
#      - You need a fast, explainable approximation
#      - Full synthetic weighting is unnecessary
# 4) Event study with synthetic-style reweighting (QP-based)
#    Purpose:
#      - Modern applied-econ workhorse
#      - Balance pre-period outcomes, then estimate dynamic effects
#    Use when:
#      - Multiple treated units
#      - Dynamic treatment effects are required
#      - Parallel trends are clearly violated
# 5) Generalized Synthetic Control (gsynth)
#    Purpose:
#      - Latent-factor / interactive fixed effects alternative
#    Use when:
#      - You want the lowest-friction SC-like workflow
#      - Many units and moderate-to-long panels

set.seed(123)

library(dplyr)
library(tidyr)
library(ggplot2)
library(gsynth)
library(Synth)
library(fixest)
library(quadprog)
library(synthdid)

##### Generate Panel Data ######
# Panel dimensions
N  <- 50   # units
Tt <- 40   # time periods
t0 <- 26   # treatment starts at time >= t0 for treated units
n_treated <- 6

treated_units <- sample(1:N, n_treated)

# Create a balanced panel
panel <- expand_grid(unit = 1:N, time = 1:Tt) %>%
  mutate(
    treated_unit = unit %in% treated_units,
    post = time >= t0,
    D = as.integer(treated_unit & post)
  )

# Data generating process:
# y_it = unit FE + time FE + latent factor structure + noise + treatment effect
unit_fe <- rnorm(N, 0, 1.0)
time_fe <- sin((1:Tt) / 6) + rnorm(Tt, 0, 0.15)

# Latent factors (interactive FE / non-parallel trends)
K <- 2
f_t <- cbind(
  scale((1:Tt))[,1],              # trending factor
  cos((1:Tt) / 5)                 # cyclical factor
)
lambda_i <- matrix(rnorm(N*K, 0, 0.8), nrow = N, ncol = K)

# Make treated units slightly different in factor loadings (harder parallel trends)
lambda_i[treated_units, 1] <- lambda_i[treated_units, 1] + 0.8

tau <- 1.2 # true treatment effect

panel <- panel %>%
  mutate(
    u = unit_fe[unit],
    g = time_fe[time],
    lf = rowSums(lambda_i[unit, , drop = FALSE] * f_t[time, , drop = FALSE]),
    eps = rnorm(n(), 0, 0.6),
    y = 2 + u + g + lf + tau * D + eps
  ) %>%
  select(unit, time, y, treated_unit, post, D)

# Quick sanity plot: average outcome by treated vs control over time
panel %>%
  group_by(time, treated_unit) %>%
  summarise(y_bar = mean(y), .groups = "drop") %>%
  ggplot(aes(time, y_bar, linetype = treated_unit)) +
  geom_line() +
  geom_vline(xintercept = t0, linetype = "dashed") +
  labs(linetype = "Treated unit", title = "Synthetic panel: treated vs control mean outcome",
       subtitle = "Dashed line = treatment start") +
  theme_minimal()

###### Synth DiD ######

# Ensure balanced, sorted panel 
panel2 <- panel %>%
  arrange(unit, time)

# Sanity checks: balanced panel and no missing outcomes
stopifnot(!anyNA(panel2$y))
stopifnot(n_distinct(panel2$unit) * n_distinct(panel2$time) == nrow(panel2))

# Define unit/time order explicitly
units <- sort(unique(panel2$unit))
times <- sort(unique(panel2$time))

N  <- length(units)
Tt <- length(times)

# Build numeric outcome matrix Y with rows=units (in 'units' order), cols=time (in 'times' order)
# Because panel2 is arranged(unit,time) and balanced, this is safe:
Y <- matrix(panel2$y, nrow = N, ncol = Tt, byrow = TRUE)

# Identify treated vs control units
treated_units <- sort(unique(panel2$unit[panel2$treated_unit]))
control_units <- setdiff(units, treated_units)

# N0 = number of control units
N0 <- length(control_units)
stopifnot(N0 > 0, N0 < N)

# T0 = number of pre-treatment periods (last pre period index)
# This assumes "post" is defined correctly in your panel (as in the synthetic generator).
T0 <- max(panel2$time[panel2$post == FALSE])
stopifnot(T0 > 0, T0 < Tt)

# Reorder rows so controls are first, treated last (this is what synthdid expects)
row_order <- c(match(control_units, units), match(treated_units, units))
stopifnot(!anyNA(row_order))

Y_ord <- Y[row_order, , drop = FALSE]

# More safety checks before estimation
stopifnot(is.numeric(Y_ord))
stopifnot(nrow(Y_ord) == N, ncol(Y_ord) == Tt)

# Estimate SDiD
tau_hat <- synthdid_estimate(Y_ord, N0 = N0, T0 = T0)
tau_hat

# SEs via bootstrap (often a good default)
se_hat <- sqrt(vcov(tau_hat, method = "bootstrap"))
c(estimate = as.numeric(tau_hat), se = se_hat)

# Plot (diagnostic visualization)
synthdid_plot(tau_hat, overlay = 1) +
  ggtitle("Synthetic DiD: observed vs synthetic counterfactual")

##### Synthetic controls as preprocessing → TWFE ######
# Build synthetic control weights (typically per treated unit), then run a TWFE regression on a 
# reweighted donor pool / restricted donors. One treated unit (pick the first treated unit), 
# build SC using Synth, then create a synthetic control series and estimate a DiD/TWFE-style effect 
# vs that synthetic.

# Data structure assumptions
# You have a long panel called `panel` with columns: unit, time, y, treated_unit (TRUE/FALSE), 
# D (0/1)  (or at least treated_unit + post)
# If your columns differ, rename them once up front.

# Standardize column names if needed 
panel_std <- panel %>%
  rename_with(~"unit",  .cols = any_of(c("unit_id", "unit", "id", "geo", "group"))) %>%
  rename_with(~"time",  .cols = any_of(c("time_id", "time", "t", "date_index"))) %>%
  rename_with(~"y",     .cols = any_of(c("y", "outcome", "Y"))) %>%
  mutate(
    unit = as.integer(unit),
    time = as.integer(time),
    y = as.numeric(y)
  ) %>%
  arrange(unit, time)

Tt <- max(panel_std$time)

pre_period <- sort(unique(panel_std$time[panel_std$time < t0]))

# Prepare Synth input using base data.frame + numeric IDs
dat_sc <- panel_std %>%
  transmute(
    unit_id = as.integer(unit),
    time_id = as.integer(time),
    y = as.numeric(y)
  ) %>%
  as.data.frame()

all_units <- sort(unique(dat_sc$unit_id))

control_units <- setdiff(all_units, tr_unit)

# Simplest stable configuration: use lagged y in pre-period as special predictors
dp <- dataprep(
  foo = dat_sc,
  dependent = "y",
  unit.variable = "unit_id",
  time.variable = "time_id",
  treatment.identifier = tr_unit,
  controls.identifier = control_units,
  predictors = "y",
  predictors.op = "mean",
  time.predictors.prior = pre_period,
  special.predictors = lapply(pre_period, function(tt) list("y", tt, "mean")),
  time.optimize.ssr = pre_period,
  time.plot = sort(unique(dat_sc$time_id))
)

# Fit synthetic control 
sc_out <- synth(dp)

# Simplest way to extract weights: use the control outcome matrix used by Synth (Y0plot) 
# and multiply by weights.

w_full <- as.numeric(sc_out$solution.w) 

Y0plot <- dp$Y0plot                      

# Synthetic outcome series for all plotted times:
y_synth <- tibble(
  time = as.integer(rownames(Y0plot)),
  y_synth = as.numeric(Y0plot %*% w_full))

# Treated series:
Y1plot <- dp$Y1plot   

y_treated <- tibble(
  time = as.integer(rownames(Y1plot)),
  y_treated = as.numeric(Y1plot))

# Gap series:
sc_series <- y_treated %>%
  left_join(y_synth, by = "time") %>%
  mutate(
    post = time >= t0,
    gap = y_treated - y_synth)

# Simple DiD-on-gap estimate (post shift in the treated-synth gap)
did_sc <- lm(gap ~ post, data = sc_series)

summary(did_sc)

# Donor-weight table (robust extraction)
# Donor IDs are best read from colnames of Y0plot when available.
donor_ids <- colnames(Y0plot)

# Plot treated vs synthetic

plot_df <- sc_series %>%
  select(time, y_treated, y_synth) %>%
  tidyr::pivot_longer(cols = c(y_treated, y_synth), names_to = "series", values_to = "value")

ggplot(plot_df, aes(time, value, linetype = series)) +
  geom_line() +
  geom_vline(xintercept = t0, linetype = "dashed") +
  theme_minimal() +
  labs(
    title = "Synthetic Control Preprocessing: Treated vs Synthetic",
    subtitle = paste("Treated unit =", tr_unit, "| Treatment start t0 =", t0),
    linetype = "")+
  theme_fancy()

# Plot gap
ggplot(sc_series, aes(time, gap)) +
  geom_line() +
  geom_vline(xintercept = t0, linetype = "dashed") +
  theme_minimal() +
  labs(
    title = "Gap (Treated - Synthetic)",
    subtitle = "If fit is good pre-treatment, post shift estimates treatment effect")+
  theme_fancy()

# 3. SC-style preprocessing → TWFE (regression-based)
# Avoids weight extraction headache

# Keep only units with good pre-period overlap (example heuristic)
pre_means <- panel %>%
  filter(post == FALSE) %>%
  group_by(unit) %>%
  summarise(pre_mean = mean(y), .groups = "drop")

treated_pre <- pre_means %>%
  filter(unit %in% unique(panel$unit[panel$treated_unit])) %>%
  summarise(mu = mean(pre_mean)) %>%
  pull(mu)

# Keep controls within a reasonable pre-period band
panel_restricted <- panel %>%
  left_join(pre_means, by = "unit") %>%
  filter(
    treated_unit |
      abs(pre_mean - treated_pre) <= sd(pre_means$pre_mean)
  )

# TWFE on the restricted sample
twfe_sc_like <- feols(
  y ~ D | unit + time,
  data = panel_restricted,
  vcov = ~unit
)

summary(twfe_sc_like)

# 4. Event study with synthetic-style reweighting: balance pre-period → run weighted event study.

# Pre-period
# Pre-period data
pre <- panel %>% filter(post == FALSE)

treated_units <- sort(unique(panel$unit[panel$treated_unit]))

control_units <- setdiff(sort(unique(panel$unit)), treated_units)

# Time index for pre-period
pre_times <- sort(unique(pre$time))

Tpre <- length(pre_times)

# Average treated pre-period path: length Tpre
treated_pre <- pre %>%
  filter(unit %in% treated_units) %>%
  group_by(time) %>%
  summarise(y = mean(y), .groups = "drop") %>%
  arrange(time) %>%
  pull(y)

# Build control outcome matrix as TIME x CONTROLS (Tpre x J)
# 1) controls x time wide
Yc_wide <- pre %>%
  filter(unit %in% control_units) %>%
  select(unit, time, y) %>%
  tidyr::pivot_wider(names_from = time, values_from = y) %>%
  arrange(unit)

# Ensure pre_times columns exist and are in correct order
# pivot_wider will name columns like "1", "2", ... if time is integer
Yc_mat_controls_by_time <- Yc_wide %>%
  select(-unit) %>%
  as.matrix()

# Now transpose to TIME x CONTROLS
Yc_pre <- t(Yc_mat_controls_by_time)

# Sanity checks
J <- length(control_units)
stopifnot(nrow(Yc_pre) == Tpre)
stopifnot(ncol(Yc_pre) == J)

# Solve for weights w:
# minimize || treated_pre - Yc_pre %*% w ||^2
# with w >= 0, sum(w) = 1
Dmat <- 2 * (t(Yc_pre) %*% Yc_pre)         

dvec <- 2 * (t(Yc_pre) %*% treated_pre)    

Amat <- cbind(diag(J), rep(1, J), rep(-1, J))

bvec <- c(rep(0, J), 1, -1)

qp <- solve.QP(Dmat = Dmat + diag(1e-8, J), dvec = dvec, Amat = Amat, bvec = bvec)
w <- pmax(qp$solution, 0)
w <- w / sum(w)

weights_df <- tibble(unit = control_units, w = w)

# Attach weights: treated weight = 1; controls weight = w
panel_w <- panel %>%
  left_join(weights_df, by = "unit") %>%
  mutate(
    w_final = ifelse(treated_unit, 1, w),
    w_final = ifelse(is.na(w_final), 0, w_final)
  )

# Define event time using first treated post period (global t0)
t0 <- min(panel$time[panel$D == 1])
panel_w <- panel_w %>%
  mutate(event_time = time - t0)

# Weighted event study with time FE (reference = -1)
es_fit <- feols(
  y ~ i(event_time, treated_unit, ref = -1) | time,
  data = panel_w,
  weights = ~ w_final,
  vcov = "HC1"
)

summary(es_fit)
iplot(es_fit, main = "Weighted event-study (synthetic-style reweighting)", 
      xlab = "Event time (t - t0)")

# 6. gsynth: 
# A. single common adoption
# B. staggered adoption
# Requires a long panel `panel` with at least: unit (unit id), time (time id), y (outcome),
# and either:
#   (A) a treatment indicator D (0/1) with a common adoption date across treated units, OR
#   (B) a treatment indicator D (0/1) with staggered adoption (unit-specific first treated time)
# This script supports BOTH by:
#   - Running model (A) directly from D (common adoption)
#   - Building `g` = first treated time per unit and running model (B) with gname="g" (staggered)

panel_gs <- panel %>%
  transmute(
    unit = as.integer(unit),
    time = as.integer(time),
    y    = as.numeric(y),
    D    = as.integer(D)
  ) %>%
  arrange(unit, time)


# A. Single common adoption date if design truly has one common adoption date, D is sufficient.

fit_single <- gsynth(
  y ~ D,
  data = panel_gs,
  index = c("unit", "time"),
  force = "two-way",
  r = c(0, 5),
  CV = TRUE,
  se = TRUE,
  inference = "parametric"   # or "bootstrap"
)

# the F-statistic it tries to report during inference needs multiple treated units 
# (so it can form a cross-sectional variance estimate for that test). If you have only 1 
# treated unit (or very few), gsynth can still estimate ATT/counterfactuals, but it can’t 
# compute that particular F test.

summary(fit_single)

plot(fit_single, type = "gap") # ATT over time

plot(fit_single, type = "counterfactual") # observed vs counterfactual

# Boostrapped
fit_single_boot <- gsynth(
  y ~ D,
  data = panel_gs,
  index = c("unit", "time"),
  force = "two-way",
  r = c(0, 5),
  CV = TRUE,
  se = TRUE,
  inference = "bootstrap",  # <— instead of parametric
  nboots = 500              # increase if you want smoother bands
)

summary(fit_single_boot)

plot(fit_single_boot, type = "att")

# ----------------------------
# B) STAGGERED ADOPTION DATE
# ----------------------------
# Build g = first treated time per unit (unit-specific adoption date).
# Convention:
#   - For treated units: g = first time where D==1
#   - For never-treated units: g = 0 (or NA; gsynth typically expects 0 for never-treated)

g_map <- panel_gs %>%
  group_by(unit) %>%
  summarise(
    g = ifelse(any(D == 1), min(time[D == 1]), 0L),
    .groups = "drop"
  )

panel_stag <- panel_gs %>%
  left_join(g_map, by = "unit")

# Sanity checks
stopifnot(!anyNA(panel_stag$g))
stopifnot(all(panel_stag$g %in% c(0L, sort(unique(panel_stag$time)))))

# Run gsynth with gname (staggered adoption)
# Notes:
# - With gname, gsynth uses g to define treated periods unit-by-unit.
# - You can still keep D in the formula or set y ~ 1; using y ~ 1 is common with gname.
#   We'll use y ~ 1 to let gsynth infer treatment timing from g.

fit_staggered <- gsynth(
  y ~ 1,
  data = panel_stag,
  index = c("unit", "time"),
  force = "two-way",
  r = c(0, 5),
  CV = TRUE,
  se = TRUE,
  inference = "parametric",  # or "bootstrap"
  gname = "g"
)

summary(fit_staggered)
plot(fit_staggered, type = "att") # ATT over time (aligned by calendar time)
plot(fit_staggered, type = "counterfactual") # observed vs counterfactual

# Optional: inspect fit objects to locate stored estimates/series for programmatic use
# str(fit_single, max.level = 2)
# str(fit_staggered, max.level = 2)


# Event-time alignment for STAGGERED adoption (using gsynth fit if available)
# Preconditions:
# - fit_staggered created with gsynth(..., gname="g")
# - panel_stag has columns: unit, time, y, D, g (first treated time; 0 for never-treated)

stopifnot(exists("fit_staggered"))
stopifnot(all(c("unit","time","y","D","g") %in% names(panel_stag)))

# Helper: build event time (relative time) for treated units only
panel_evt <- panel_stag %>%
  mutate(
    event_time = ifelse(g > 0, time - g, NA_integer_),
    treated_ever = g > 0
  )

# Try to extract a unit-time GAP from fit_staggered
# Different versions expose different slots. 
nms <- names(fit_staggered)

# Candidates (common patterns in fect/gsynth objects)
candidates <- c("gap", 
                "att", 
                "effect", 
                "eff", 
                "tau", 
                "DID", 
                "Y.ct", 
                "Y0", 
                "Y0.hat", 
                "Yhat", 
                "Y.counterfactual")

present <- candidates[candidates %in% nms]
message("Slots found among common candidates: ", 
        paste(present, collapse = ", "))

# The most direct scenario: fit_staggered$gap is a matrix/data.frame with unit/time
# If you know your object has a specific slot, replace this logic with that.
gap_obj <- NULL
if ("gap" %in% nms) gap_obj <- fit_staggered$gap

if (!is.null(gap_obj)) {
  # Try to coerce to long format
  # Commonly gap_obj is a matrix with rows=units and cols=times OR vice-versa.
  if (is.matrix(gap_obj)) {
    # Attempt to interpret dimnames
    dn <- dimnames(gap_obj)
    # If both dimnames exist and look numeric, use them
    if (!is.null(dn[[1]]) && !is.null(dn[[2]])) {
      gap_long <- as.data.frame(as.table(gap_obj)) %>%
        rename(unit = Var1, time = Var2, gap = Freq) %>%
        mutate(
          unit = suppressWarnings(as.integer(as.character(unit))),
          time = suppressWarnings(as.integer(as.character(time)))
        )
    } else {
      stop("fit_staggered$gap matrix lacks dimnames; cannot reliably map to unit/time.")
    }
  } else if (is.data.frame(gap_obj)) {
    gap_long <- gap_obj
  } else {
    stop("fit_staggered$gap exists but is not a matrix/data.frame.")
  }
  
  # Join g to compute event time and aggregate
  gap_evt <- gap_long %>%
    left_join(panel_stag %>% distinct(unit, g), by = "unit") %>%
    filter(g > 0) %>%
    mutate(event_time = time - g)
  
  # Average effect by event time
  att_event <- gap_evt %>%
    group_by(event_time) %>%
    summarise(att = mean(gap, na.rm = TRUE), n = n(), .groups = "drop")
  
  # Plot
  ggplot(att_event, aes(event_time, att)) +
    geom_line() +
    geom_vline(xintercept = -1, linetype = "dashed") +
    theme_minimal() +
    labs(
      title = "Staggered adoption: ATT aligned by event time (from gsynth gap)",
      x = "Event time (t - g_i)",
      y = "Average gap (treated - counterfactual)"
    )
  
} else {
  message(
    "Could not find fit_staggered$gap. ",
    "Use the fallback event-study alignment below (Section B)."
  )
}




# =============================================================================
# IMPROVED EVENT STUDY BLOCK
# This block runs an event study on the gap between observed and synthetic outcomes
# It uses Sun and Abraham style interactions via fixest::sunab when possible
# It includes unit and time fixed effects and clustered standard errors
# It performs a joint pre trend test and shows alternative estimators and checks
# =============================================================================

library(dplyr)
library(fixest)
library(broom)
library(ggplot2)

# Replace these names if your dataset uses different columns:
unit_col <- "treated_city"   # unit id column
time_col <- "date"          # time column, Date or integer
y_col <- "treatment_sales"  # observed outcome column
y_synth_col <- "synthetic_sales"  # synthetic counterfactual column
g_col <- "intervention_date"  # first treated time per unit, NA for never treated

# Prepare data frame and compute gap and event time per unit
df_es <- synthetic_panel %>%
  arrange(.data[[unit_col]], .data[[time_col]]) %>%
  # compute gap as observed minus synthetic
  mutate(gap = .data[[y_col]] - .data[[y_synth_col]]) %>%
  # ensure intervention date exists per unit, compute if missing and post exists
  group_by(.data[[unit_col]]) %>%
  mutate(
    # if intervention date column is not present or is NA compute as first date where post == 1
    !!g_col := ifelse(is.na(.data[[g_col]]) & "post" %in% names(.),
                      ifelse(any(post == 1), min(.data[[time_col]][post == 1]), as.Date(NA)),
                      .data[[g_col]])
  ) %>%
  ungroup() %>%
  mutate(
    # compute event time for treated units; NA for never treated
    event_time = ifelse(!is.na(.data[[g_col]]), as.integer(.data[[time_col]] - .data[[g_col]]), NA_integer_),
    treated_flag = ifelse(!is.na(.data[[g_col]]), 1L, 0L)
  )

# Bin extreme leads and lags to stabilize coefficients
K <- 8  # adjust number of leads and lags to display
df_es <- df_es %>%
  mutate(event_time_binned = ifelse(is.na(event_time), NA_integer_, pmin(pmax(event_time, -K), K)))

# Quick counts per event time to check sample sizes
counts_ev <- df_es %>% group_by(event_time_binned) %>% summarise(n = n(), .groups = "drop")
print(counts_ev)

# Compute pre period rmspe per unit to flag poor synthetic fits
rmspe_by_unit <- df_es %>%
  filter(!is.na(.data[[g_col]])) %>%
  group_by(.data[[unit_col]]) %>%
  filter(event_time < 0) %>%
  summarise(rmspe = sqrt(mean((.data[[y_col]] - .data[[y_synth_col]])^2, na.rm = TRUE)),
            n_pre = n(), .groups = "drop") %>%
  arrange(desc(rmspe))
print(rmspe_by_unit)

# Remove units without a valid intervention date for the event study
df_evt <- df_es %>% filter(!is.na(.data[[g_col]]))

# Use Sun and Abraham style event study on the gap
es_fit_gap <- feols(
  gap ~ sunab(event_time_binned, treated_flag, ref = -1) |
    get(unit_col) + get(time_col),
  data = df_evt,
  cluster = ~get(unit_col)
)

# Print summary and aggregated ATT
print(summary(es_fit_gap, agg = "att"))

# Tidy coefficients for plotting
es_tidy <- broom::tidy(es_fit_gap, conf.int = TRUE)

# Extract event time terms. Term naming can vary by fixest version.
# Here we attempt to parse terms that encode event_time_binned
evt_terms <- es_tidy %>% filter(grepl("event_time_binned", term) | grepl("sunab", term))
# Attempt to extract numeric event time from term names
extract_evt <- function(term) {
  # remove non digits and decimal and minus sign but keep numbers like -3
  m <- regmatches(term, regexpr("(-?\\d+)$", term))
  if (length(m) && nchar(m) > 0) return(as.integer(m))
  # fallback NA
  return(NA_integer_)
}
evt_vals <- vapply(evt_terms$term, extract_evt, integer(1))
evt_plot <- evt_terms %>% mutate(evt = evt_vals) %>% arrange(evt)

# Plot event study coefficients with confidence intervals
ggplot(evt_plot, aes(x = evt, y = estimate)) +
  geom_point() +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) +
  geom_vline(xintercept = -0.5, linetype = "dashed") +
  labs(x = "Event time", y = "Estimate on gap", title = "Event study on gap using Sun Abraham") +
  theme_minimal()

# Joint pre trend test using a Wald statistic
pre_terms <- evt_plot %>% filter(evt < 0)
if (nrow(pre_terms) >= 2) {
  coef_names <- pre_terms$term
  coefs <- coef(es_fit_gap)[coef_names]
  V <- vcov(es_fit_gap, cluster = ~get(unit_col))[coef_names, coef_names, drop = FALSE]
  wald_stat <- as.numeric(t(coefs) %*% solve(V, coefs))
  df_w <- length(coefs)
  pval_pre <- 1 - pchisq(wald_stat, df_w)
  cat("Pre trend joint Wald p value:", signif(pval_pre, 3), "\n")
} else {
  message("Not enough pre event coefficients to run joint pre trend test")
}

# If number of clusters is small consider wild cluster bootstrap
# Example usage comment only
# fixest::boottest(es_fit_gap, clustid = df_evt[[unit_col]], param = "sunab::event_time_binned::0")

# =============================================================================
# Alternatives and extra checks
# =============================================================================

# Alternative 1: Event study using Callaway Sant Anna via did package
if (requireNamespace("did", quietly = TRUE)) {
  library(did)
  df_did <- df_evt %>%
    mutate(unit = as.character(.data[[unit_col]]),
           time_num = as.integer(as.Date(.data[[time_col]]))) %>%
    filter(!is.na(.data[[g_col]]))
  # gname should be numeric time of first treatment
  att_gt_res <- att_gt(
    yname = "gap",
    tname = "time_num",
    idname = "unit",
    gname = as.integer(as.Date(df_did[[g_col]])),
    data = df_did,
    xformla = ~1,
    clustervars = "unit"
  )
  agg_res <- aggte(att_gt_res, type = "dynamic")
  print(summary(agg_res))
  plot(agg_res)
} else {
  message("Package did not load. Install did to run Callaway Sant Anna estimator")
}

# Alternative 2: Run event study on raw outcome as a robustness check
# Note: Running on raw y may double use pre period data used to construct synthetic
# This is only a robustness check not the primary estimator
es_fit_raw <- feols(
  .data[[y_col]] ~ sunab(event_time_binned, treated_flag, ref = -1) |
    get(unit_col) + get(time_col),
  data = df_evt,
  cluster = ~get(unit_col)
)
print(summary(es_fit_raw, agg = "att"))

# Extra check 1: placebo event study using earlier fake interventions
# This selects a range of placebo starts in the pre period and computes placebo ATTs
placebo_event_starts <- df_evt %>%
  group_by(.data[[unit_col]]) %>%
  filter(event_time < 0) %>%
  summarise(min_evt = min(event_time), max_evt = max(event_time), .groups = "drop") %>%
  pull(min_evt)

# Run a simple placebo loop for a few candidate placebo shifts
placebo_results <- list()
for (shift in seq(-4, -1)) {
  df_placebo <- df_evt %>%
    mutate(event_time_placebo = ifelse(!is.na(event_time), event_time - shift, NA_integer_),
           event_time_placebo = pmin(pmax(event_time_placebo, -K), K))
  fit_pl <- feols(
    gap ~ sunab(event_time_placebo, treated_flag, ref = -1) |
      get(unit_col) + get(time_col),
    data = df_placebo,
    cluster = ~get(unit_col)
  )
  placebo_results[[as.character(shift)]] <- summary(fit_pl, agg = "att")
}

# Print placebo summaries
print(placebo_results)

# Extra check 2: compare gap based event study to CausalImpact and placebo test outputs
# If you have existing objects named ci_summary or placebo_results_by_city join them here
# Example join code comment only
# final_report <- ci_summary %>% left_join(placebo_results_by_city, by = "treated_city") %>% left_join(mde_results, by = "treated_city")

# End of improved event study block
# =============================================================================

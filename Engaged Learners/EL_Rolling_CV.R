# ============================================================
# Train-only Rolling-Origin CV + Final Holdout Test (used once at the end)
# Data: el_data (one row per user per day)
# Target: retained7d (0/1; NA rows dropped)
# ============================================================

library(tidyverse)
library(tidymodels)
library(doParallel)

# ---------------------------
# 0) Inputs
# ---------------------------
outcome  <- "retained7d"
time_col <- "date"     # <-- change if needed
id_cols  <- c("user_id")        # <-- optional; can be c()

# Final holdout size (days)
holdout_days <- 28              # last 28 days as final test set

# Rolling CV parameters (in days) on TRAIN ONLY
lookback_days <- 120
assess_days   <- 28
step_days     <- 28
cumulative    <- TRUE

# ---------------------------
# 1) Helper: day-based rolling splits (calendar days, not rows)
# ---------------------------
build_day_rolling_splits <- function(data, date_col,
                                     lookback_days, assess_days, step_days, cumulative) {
  dates <- data %>%
    distinct(.data[[date_col]]) %>%
    arrange(.data[[date_col]]) %>%
    pull(.data[[date_col]])
  
  splits <- list()
  ids <- character(0)
  
  min_date <- min(dates, na.rm = TRUE)
  max_date <- max(dates, na.rm = TRUE)
  
  first_assess_start <- min_date + lookback_days
  assess_start_dates <- seq.Date(from = first_assess_start, to = max_date, by = step_days)
  
  fold <- 0L
  for (assess_start in assess_start_dates) {
    assess_end <- assess_start + (assess_days - 1)
    if (assess_end > max_date) break
    
    train_start <- if (cumulative) min_date else (assess_start - lookback_days)
    train_end   <- assess_start - 1
    
    analysis_idx <- which(data[[date_col]] >= train_start & data[[date_col]] <= train_end)
    assess_idx   <- which(data[[date_col]] >= assess_start & data[[date_col]] <= assess_end)
    
    if (length(analysis_idx) == 0 || length(assess_idx) == 0) next
    
    fold <- fold + 1L
    splits[[fold]] <- rsample::make_splits(
      list(analysis = analysis_idx, assessment = assess_idx),
      data
    )
    ids <- c(ids, paste0("Fold_", fold, "_", assess_start, "_to_", assess_end))
  }
  
  rsample::manual_rset(splits = splits, ids = ids)
}

# ---------------------------
# 2) Prep el_data (make outcome a factor BEFORE any splitting)
# ---------------------------
stopifnot(exists("el_data"))
stopifnot(outcome %in% names(el_data))
stopifnot(time_col %in% names(el_data))

el_data<- el_data %>%
  ungroup() %>%
  mutate(
    !!time_col := as.Date(.data[[time_col]]),
    !!outcome := case_when(
      .data[[outcome]] %in% c(1, "1", TRUE, "TRUE", "T", "true")   ~ "1",
      .data[[outcome]] %in% c(0, "0", FALSE, "FALSE", "F", "false") ~ "0",
      TRUE ~ NA_character_
    ),
    !!outcome := factor(.data[[outcome]], levels = c("0", "1"))
  ) %>%
  filter(!is.na(.data[[outcome]])) %>%
  arrange(.data[[time_col]])

stopifnot(nrow(df) > 0)

# ---------------------------
# 3) Create final holdout test window (UNTOUCHED until the end)
# ---------------------------
all_dates <- df %>%
  distinct(.data[[time_col]]) %>%
  arrange(.data[[time_col]]) %>%
  pull(.data[[time_col]])

max_date <- max(all_dates, na.rm = TRUE)
holdout_start <- max_date - (holdout_days - 1)

train_df <- df %>% filter(.data[[time_col]] < holdout_start)
test_df  <- df %>% filter(.data[[time_col]] >= holdout_start)

# Sanity checks
range(train_df[[time_col]])
range(test_df[[time_col]])

# ---------------------------
# 4) Rolling-origin CV on TRAIN ONLY
# ---------------------------
roll_cv <- build_day_rolling_splits(
  data = train_df,
  date_col = time_col,
  lookback_days = lookback_days,
  assess_days = assess_days,
  step_days = step_days,
  cumulative = cumulative
)

# ---------------------------
# 5) Recipe (sensible for boosted trees; rare-level bucketing reduces memory)
# ---------------------------
clf_recipe <-
  recipe(formula(paste(outcome, "~ .")), data = train_df) %>%
  update_role(all_of(c(id_cols, time_col)), new_role = "id") %>%
  step_other(all_nominal_predictors(), threshold = 0.005) %>%  # helps memory
  step_unknown(all_nominal_predictors()) %>%
  step_novel(all_nominal_predictors()) %>%
  step_impute_median(all_numeric_predictors()) %>%
  step_impute_mode(all_nominal_predictors()) %>%
  step_dummy(all_nominal_predictors(), one_hot = TRUE) %>%
  step_zv(all_predictors())

# ---------------------------
# 6) XGBoost spec (CV) - fixed params
# ---------------------------
xgb_spec_cv <-
  boost_tree(
    trees = 400,           # start smaller for speed; scale up later
    learn_rate = 0.1,
    tree_depth = 3,
    min_n = 20,
    loss_reduction = 0.0,
    sample_size = 0.8,
    mtry = 0.8
  ) %>%
  set_mode("classification") %>%
  set_engine(
    "xgboost",
    objective = "binary:logistic",
    eval_metric = "auc",
    nthread = 1,       # parallelize over folds, not within each fit
    counts = FALSE
  )

clf_wf <- workflow() %>%
  add_recipe(clf_recipe) %>%
  add_model(xgb_spec_cv)

# ---------------------------
# 7) Parallelized rolling CV (TRAIN ONLY)
# ---------------------------
n_workers <- max(1, parallel::detectCores(logical = FALSE) - 1)
cl <- parallel::makePSOCKcluster(n_workers)
doParallel::registerDoParallel(cl)

cv_res <- fit_resamples(
  clf_wf,
  resamples = roll_cv,
  metrics = metric_set(roc_auc, pr_auc, accuracy, mn_log_loss),
  control = control_resamples(
    save_pred = FALSE,          # memory-safe
    extract = NULL,             # memory-safe
    parallel_over = "resamples",
    verbose = TRUE
  )
)

parallel::stopCluster(cl)
doParallel::stopImplicitCluster()

cv_metrics <- collect_metrics(cv_res)
cv_metrics

# ---------------------------
# 8) Final model fit on ALL TRAIN data (still not using holdout)
#    Allow XGBoost to use CPU threads since it's a single fit
# ---------------------------
xgb_spec_final <- xgb_spec_cv %>%
  set_engine(
    "xgboost",
    objective = "binary:logistic",
    eval_metric = "auc",
    nthread = max(1, parallel::detectCores(logical = FALSE) - 1),
    counts = FALSE
  )

final_wf <- workflow() %>%
  add_recipe(clf_recipe) %>%
  add_model(xgb_spec_final)

final_fit <- fit(final_wf, data = train_df)

# ---------------------------
# 9) Final evaluation on HOLDOUT TEST (used once at the end)
# ---------------------------
test_probs <- predict(final_fit, new_data = test_df, type = "prob")
test_class <- predict(final_fit, new_data = test_df, type = "class")

test_scored <- bind_cols(
  test_df %>% select(all_of(c(id_cols, time_col, outcome))),
  test_class,
  test_probs
)

prob_col <- ".pred_1"   # because levels are c("0","1")
event_level <- "second"

test_metrics <- tibble(
  roc_auc  = roc_auc(test_scored, truth = all_of(outcome), estimate = .data[[prob_col]], event_level = event_level)$.estimate,
  pr_auc   = pr_auc(test_scored,  truth = all_of(outcome), estimate = .data[[prob_col]], event_level = event_level)$.estimate,
  accuracy = accuracy(test_scored, truth = all_of(outcome), estimate = .pred_class)$.estimate,
  log_loss = mn_log_loss(test_scored, truth = all_of(outcome), estimate = .data[[prob_col]])$.estimate
)

test_metrics
conf_mat(test_scored, truth = all_of(outcome), estimate = .pred_class)
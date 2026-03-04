library(tidymodels)
library(rsample)
library(dplyr)

set.seed(7)

# ---------------------------
# 0) Ensure correct types
# ---------------------------
el_data <- el_data %>%
  mutate(
    date = as.Date(date),
    !!outcome := factor(.data[[outcome]], levels = c("0", "1"))
  )

event_level <- "second"
prob_col <- ".pred_1"

# ---------------------------
# 1) 3-fold grouped CV by user_id
# ---------------------------
folds <- group_vfold_cv(
  data  = el_data,
  group = user_id,
  v     = 3
)

# ---------------------------
# 2) Recipe
# ---------------------------
clf_recipe <- recipe(formula(paste(outcome, "~ .")), data = el_data) %>%
  update_role(all_of(c(id_cols, time_col)), new_role = "id") %>%
  step_other(all_nominal_predictors(), threshold = 0.01, other = "other_level") %>%
  step_unknown(all_nominal_predictors(), new_level = "unknown") %>%
  step_novel(all_nominal_predictors()) %>%
  step_dummy(all_nominal_predictors(), one_hot = TRUE) %>%
  step_zv(all_predictors())

# ---------------------------
# 3) XGBoost model spec
# ---------------------------
cores_to_use <- max(1, parallel::detectCores(logical = FALSE) - 1)

xgb_spec <- boost_tree(
  trees = 500,
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
    objective   = "binary:logistic",
    eval_metric = "auc",
    counts      = FALSE,
    nthread     = cores_to_use,
    verbose     = 0
  )

# ---------------------------
# 4) Workflow
# ---------------------------
clf_wf <- workflow() %>%
  add_recipe(clf_recipe) %>%
  add_model(xgb_spec)

# ---------------------------
# 5) Fit resamples
# ---------------------------
prob_metrics <- metric_set(roc_auc, pr_auc, mn_log_loss)

cv_res <- fit_resamples(
  object    = clf_wf,
  resamples = folds,
  metrics   = prob_metrics,
  control   = control_resamples(
    save_pred = TRUE,
    verbose = TRUE,
    parallel_over = "resamples"
  )
)

# ---------------------------
# 6) Collect fold-level metrics (force event_level = "second")
# ---------------------------
preds <- collect_predictions(cv_res)

auc_by_fold <- preds %>%
  group_by(id) %>%
  roc_auc(truth = !!sym(outcome), .pred_1, event_level = event_level) %>%
  ungroup()

pr_by_fold <- preds %>%
  group_by(id) %>%
  pr_auc(truth = !!sym(outcome), .pred_1, event_level = event_level) %>%
  ungroup()

logloss_by_fold <- preds %>%
  group_by(id) %>%
  mn_log_loss(truth = !!sym(outcome), .pred_1, event_level = event_level) %>%
  ungroup()

acc_by_fold <- preds %>%
  group_by(id) %>%
  accuracy(truth = !!sym(outcome), estimate = .pred_class) %>%
  ungroup()

# ---------------------------
# 7) Summary table
# ---------------------------
summary_tbl <- bind_rows(
  auc_by_fold %>% mutate(metric = "roc_auc"),
  pr_by_fold  %>% mutate(metric = "pr_auc"),
  logloss_by_fold %>% mutate(metric = "mn_log_loss"),
  acc_by_fold %>% mutate(metric = "accuracy")
) %>%
  group_by(metric) %>%
  summarise(
    mean = mean(.estimate),
    sd   = sd(.estimate),
    min  = min(.estimate),
    max  = max(.estimate),
    .groups = "drop"
  )

summary_tbl
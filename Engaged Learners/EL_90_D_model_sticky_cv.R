
library(tidymodels)
library(dplyr)

#### 0) Data prep (unchanged) ####
el_data <- el_data %>%
  ungroup() %>%
  mutate(
    !!time_col := as.Date(.data[[time_col]]),
    !!outcome  := as.character(.data[[outcome]])
  ) %>%
  mutate(
    !!outcome := case_when(
      .data[[outcome]] %in% c("1", "TRUE", "T", "true", 1, TRUE)  ~ "1",
      .data[[outcome]] %in% c("0", "FALSE", "F", "false", 0, FALSE) ~ "0",
      TRUE ~ NA_character_
    ),
    !!outcome := factor(.data[[outcome]], levels = c("0", "1"))
  ) %>%
  filter(!is.na(.data[[outcome]])) %>%   # drop censored days
  arrange(.data[[time_col]]) %>%
  select(-retained7d, -retained90d, -retained28d, -retained28d_sticky, -week_start)

#### USER LEVEL SPLIT (unchanged) ####
user_split_1 <- initial_split(el_data %>% distinct(user_id), prop = 0.80)
train_users <- training(user_split_1)$user_id
tmp_users   <- testing(user_split_1)$user_id
user_split_2 <- initial_split(tibble(user_id = tmp_users), prop = 0.50)
validation_users <- training(user_split_2)$user_id
test_users       <- testing(user_split_2)$user_id

train_data <- el_data %>% filter(user_id %in% train_users)
validation_data <- el_data %>% filter(user_id %in% validation_users)
test_data <- el_data %>% filter(user_id %in% test_users)

stopifnot(length(intersect(unique(train_data$user_id), unique(validation_data$user_id))) == 0)
stopifnot(length(intersect(unique(train_data$user_id), unique(test_data$user_id))) == 0)
stopifnot(length(intersect(unique(validation_data$user_id), unique(test_data$user_id))) == 0)

rm(el_data)

# Class weight (kept for engine param)
spw <- sum(train_data[[outcome]] == "1") / sum(train_data[[outcome]] == "0")

#### 3) Recipe (kept your original; you can simplify if no categorical predictors) ####
clf_recipe <- recipe(formula(paste(outcome, "~ .")), data = train_data) %>%
  update_role(all_of(c(id_cols, time_col)), new_role = "id") %>%
  step_unknown(all_nominal_predictors(), new_level = "unknown") %>%
  step_novel(all_nominal_predictors()) %>%
  step_other(all_nominal_predictors(), threshold = 0.01, other = "other_level") %>%
  step_dummy(all_nominal_predictors(), one_hot = TRUE) %>%
  step_zv(all_predictors())

#### 4) Fast model spec for CV (reduce trees for speed) ####
cores_to_use <- max(1, parallel::detectCores(logical = FALSE) - 1)

# Use fewer trees during CV to speed things up
xgb_spec_cv <- boost_tree(
  trees = 500,          # smaller for fast CV
  learn_rate = 0.05,    # slightly lower learning rate (you can match baseline if desired)
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
    eval_metric = c("auc", "aucpr"),
    counts = FALSE,
    scale_pos_weight = spw,
    nthread = cores_to_use,
    verbose = 0,
    tree_method = "hist"
  )

# Workflow for CV
clf_wf_cv <- workflow() %>%
  add_recipe(clf_recipe) %>%
  add_model(xgb_spec_cv)

#### 5) 3-fold CV (stratified by outcome) ####
# This evaluates the model quickly across 3 folds to reduce variance vs a single split.
folds <- vfold_cv(train_data, v = 5, strata = !!sym(outcome))

metrics <- metric_set(roc_auc, pr_auc, mn_log_loss)

# control: save predictions to inspect later; parallel_over = "everything" uses future backend if set
ctrl <- control_resamples(save_pred = TRUE, parallel_over = "everything", verbose = TRUE)

# Run CV (this will be ~3x a single fit with trees=300)
res_cv <- fit_resamples(
  clf_wf_cv,
  resamples = folds,
  metrics = metrics,
  control = ctrl
)

# Quick summary of CV results
cv_metrics <- collect_metrics(res_cv)

cv_metrics

# Optional: inspect CV prediction-level performance (slower)
# cv_preds <- collect_predictions(res_cv)

#### 6) Quick final fit on TRAIN + VALIDATION using the same fast spec ####
train_val_data <- bind_rows(train_data, validation_data)

# No finalize_workflow() needed because xgb_spec_cv has no tune() parameters
final_wf_fast <- clf_wf_cv

final_fit_fast <- fit(final_wf_fast, data = train_val_data)

#### 7) Evaluate fast final model on TEST ####
prob_col <- ".pred_1"
event_level <- "second"

test_probs <- predict(final_fit_fast, new_data = test_data, type = "prob")
test_class <- predict(final_fit_fast, new_data = test_data, type = "class")

test_scored <- bind_cols(
  test_data %>% select(all_of(c(id_cols, time_col, outcome))),
  test_class,
  test_probs
) %>%
  mutate(!!outcome := factor(.data[[outcome]], levels = c("0", "1")))

metrics_class <- metric_set(accuracy)(
  test_scored,
  truth = !!sym(outcome),
  estimate = .pred_class
)

metrics_prob <- metric_set(roc_auc, pr_auc, mn_log_loss)(
  test_scored,
  truth = !!sym(outcome),
  .data[[prob_col]],
  event_level = event_level
)

test_metrics <- bind_rows(metrics_class, metrics_prob)

test_metrics

# Retention Deciles
test_scored %>%
  arrange(desc(.pred_1)) %>%
  mutate(rank = row_number(),
         pct = rank / n()) %>%
  mutate(decile = ntile(.pred_1, 10)) %>%
  group_by(decile) %>%
  summarise(
    retention = mean(.data[[outcome]] == "1"),
    n = n()
  )

# Lift curve
test_scored %>%
  arrange(desc(.pred_1)) %>%
  mutate(pct = row_number() / n()) %>%
  ggplot(aes(pct, .pred_1)) +
  geom_line() +
  labs(
    x = "Percent of Users Ranked by Retention Probability",
    y = "Predicted Retention Probability",
    title = "Retention Model Ranking Curve"
  ) +
  theme_minimal()
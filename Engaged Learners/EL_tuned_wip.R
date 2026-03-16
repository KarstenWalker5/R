library(tidymodels)
library(readr)
library(dplyr)
library(finetune)   # for tune_race_anova()
library(dials)
library(scales)

#### 0) Load Data + Set Columns ####
el_data <- read_csv("/Users/karstenwalker/Documents/Datasets/el_training_3_11.csv", show_col_types = FALSE) %>%
  select(-1)

# Set these
outcome  <- "retained90d_sticky"
time_col <- "date"
id_cols  <- c("user_id")

#### 0a) Data prep (outcome + time) ####
el_data <- el_data %>%
  ungroup() %>%
  mutate(
    !!time_col := as.Date(.data[[time_col]]),
    !!outcome  := as.character(.data[[outcome]])
  ) %>%
  mutate(
    !!outcome := case_when(
      .data[[outcome]] %in% c("1", "TRUE", "T", "true", 1, TRUE)     ~ "1",
      .data[[outcome]] %in% c("0", "FALSE", "F", "false", 0, FALSE)  ~ "0",
      TRUE ~ NA_character_
    ),
    !!outcome := factor(.data[[outcome]], levels = c("0", "1"))
  ) %>%
  filter(!is.na(.data[[outcome]])) %>%   # drop censored days
  arrange(.data[[time_col]]) %>%
  select(
    -retained7d, -retained90d, -retained28d, -retained28d_sticky,
    -week_start, -contains("engaged")
  )

#### 1) Setup ####
event_level <- "second"
prob_col <- ".pred_1"
cores_to_use <- max(1, parallel::detectCores(logical = FALSE) - 1)

#### 2) USER LEVEL SPLIT ####
# Sample users into train/val/test (80/10/10)

# 2a) Split USERS into train / temp (80/20)
user_split_1 <- initial_split(el_data %>% distinct(user_id), prop = 0.80)
train_users <- training(user_split_1)$user_id
tmp_users   <- testing(user_split_1)$user_id

# 2b) Split temp USERS into validation / test (50/50 => 10/10 overall)
user_split_2 <- initial_split(tibble(user_id = tmp_users), prop = 0.50)
validation_users <- training(user_split_2)$user_id
test_users       <- testing(user_split_2)$user_id

# 2c) Filter full rows by user sets
train_data <- el_data %>% filter(user_id %in% train_users)
validation_data <- el_data %>% filter(user_id %in% validation_users)
test_data <- el_data %>% filter(user_id %in% test_users)

# 2d) Sanity check: no user overlap
stopifnot(length(intersect(unique(train_data$user_id), unique(validation_data$user_id))) == 0)
stopifnot(length(intersect(unique(train_data$user_id), unique(test_data$user_id))) == 0)
stopifnot(length(intersect(unique(validation_data$user_id), unique(test_data$user_id))) == 0)

# Drop original data to save memory
rm(el_data)

# Positive weightw
spw <- sum(train_data[[outcome]] == "1") / sum(train_data[[outcome]] == "0")

#### 4) Recipe ####
clf_recipe <- recipe(formula(paste(outcome, "~ .")), data = train_data) %>%
  update_role(all_of(c(id_cols, time_col)), new_role = "id") %>%
  step_zv(all_predictors())

#### 5) Targeted XGBoost Spec (Tunable) ####
xgb_spec_tune <-
  boost_tree(
    trees = 300,            # IMPORTANT: keep tuning fast
    learn_rate = 0.05,
    tree_depth = tune(),
    min_n = tune(),
    loss_reduction = tune(),
    sample_size = tune(),   # integer count of rows
    mtry = tune()           # integer count of predictors
  ) %>%
  set_mode("classification") %>%
  set_engine(
    "xgboost",
    objective = "binary:logistic",
    eval_metric = c("auc", "aucpr"),
    tree_method = "hist",
    nthread = cores_to_use,
    verbose = 1
  )

#### 6) Workflow ####
xgb_wf <- workflow() %>%
  add_recipe(clf_recipe) %>%
  add_model(xgb_spec_tune)

#### 7) Tuning Parameter Ranges (Working) ####
xgb_params <- parameters(
  tree_depth(range = c(3L, 5L)),
  min_n(range = c(30L, 120L)),
  loss_reduction(range = c(0, 3)),
  
  # boost_tree(sample_size=) expects a proportion, so use sample_prop()
  sample_size = sample_prop(range = c(0.6, 1.0)),
  
  # With ~24 predictors, this is fine as an integer count
  mtry(range = c(8L, 24L))
)

xgb_params <- finalize(xgb_params, training = train_data)

#### 8) Tuning Grid ####
grid <- grid_latin_hypercube(xgb_params, size = 12)

#### 9) Create a single validation resample (no CV) ####
train_val_data <- bind_rows(train_data, validation_data)

val_rs <- validation_split(
  train_val_data,
  prop = nrow(train_data) / nrow(train_val_data),
  strata = !!sym(outcome)
)

#### 10) Tune on the single validation resample ####
metrics <- metric_set(roc_auc, pr_auc)

ctrl <- control_grid(
  save_pred = TRUE,
  parallel_over = "resamples"
)

res <- tune_grid(
  xgb_wf,
  resamples = val_rs,   # single validation resample
  grid = grid,
  metrics = metrics,
  control = ctrl
)

#### 11) Select Best Params (PR AUC) ####
show_best(res, metric = "pr_auc", n = 10)

best_params <- select_best(res, metric = "pr_auc")

final_wf <- finalize_workflow(xgb_wf, best_params)

#### 12) Fit Final Model on TRAIN + VALIDATION ####
final_fit <- fit(final_wf, data = train_val_data)

#### 13) Evaluate on TEST ####
test_probs <- predict(final_fit, new_data = test_data, type = "prob")

test_class <- predict(final_fit, new_data = test_data, type = "class")

test_scored <- bind_cols(
  test_data %>% select(all_of(c(id_cols, time_col, outcome))),
  test_class,
  test_probs
)

test_scored <- test_scored %>%
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
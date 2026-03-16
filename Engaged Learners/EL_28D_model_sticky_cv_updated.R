library(tidyverse)
library(tidymodels)
library(doParallel)
library(DALEX)
library(bigrquery)
library(ingredients)
library(data.table)
library(vip)

# Store your Google Cloud project ID
bq_auth()

project_id <- "compact-sylph-785"

# Define your SQL query (example using a public dataset)
sql <- "SELECT * 
        FROM `compact-sylph-785.Karsten.EL_1yr_1pct_LI_complete_3_13`  
WHERE MOD(ABS(FARM_FINGERPRINT(CAST(user_id AS STRING))), 100) < 75"

# Run the query
tb <- bq_project_query(project_id, sql)

# Load helpers
source("/Users/karstenwalker/Documents/GitHub/R/Helpers/themes.R")
source("/Users/karstenwalker/Documents/GitHub/R/Helpers/df_transformations.R")

# Set Seed
set.seed(7)

# Load data
user_clusters_3_4 <- read_csv(
  "/Users/karstenwalker/Documents/Modeling/Artifacts/user_clusters_3_4.csv",
  show_col_types = FALSE
) %>%
  select(-1) %>%
  mutate(
    week_start_date = as.Date(week_start_date),
    week_start = floor_date(week_start_date, "week", week_start = 1)
  ) %>%
  select(user_id, week_start, cluster)

el_data <- bq_table_download(tb)

el_data <- el_data %>%
  select(-uid, -reported_user_type, -country_code, -age, -platform, -account_type,-cumulative_lifetime_sessions, 
    -is_teacher_with_students, -session_count, -first_session_channel,
    -last_session_channel, -lifetime_sets_created, -never_created,
    -session_count_7d_avg, -set_pageviews_count_7d_avg, -had_first_session,
    -unique_sets_viewed, -set_pageviews_count, -course_count
  ) %>%
  filter(flashcards_questions_answered < quantile(flashcards_questions_answered, 0.995))


# Ensure Date types
el_data <- el_data %>%
  mutate(week_start = floor_date(date, "week", week_start = 1))

# Join clusters, no imputation
el_data <- el_data %>%
  left_join(user_clusters_3_4, by = c("user_id", "week_start"))

# data.table fill
DT <- as.data.table(el_data)
setorder(DT, user_id, date)
DT[, cluster := nafill(cluster, type = "locf"), by = user_id]
DT[, cluster := nafill(cluster, type = "nocb"), by = user_id]
DT[is.na(cluster), cluster := -1L]
DT[, cluster := as.integer(cluster)]

el_data <- as_tibble(DT)

cat(
  "Users with NO cluster coverage at all:",
  nrow(
    DT %>%
      group_by(user_id) %>%
      summarise(all_missing = all(is.na(cluster))) %>%
      filter(all_missing)
  ),
  "\n"
)

cat(
  "Users with -1 cluster assignment:",
  nrow(
    DT %>%
      filter(cluster == -1) %>%
      summarise(n_missing = n_distinct(user_id))
  ),
  "\n"
)

rm(DT)

# Optional: read in saved file
el_data <- read.csv(file = "/Users/karstenwalker/Documents/Datasets/el_training_3_10_sticky.csv") %>%
  select(-week_start)

# Set these
outcome  <- "retained28d_sticky"
time_col <- "date"
id_cols  <- c("user_id")

# 1) Data prep
el_data <- el_data %>%
  ungroup() %>%
  mutate(
    !!time_col := as.Date(.data[[time_col]]),
    !!outcome  := as.character(.data[[outcome]])
  ) %>%
  mutate(
    !!outcome := case_when(
      .data[[outcome]] %in% c("1", "TRUE", "T", "true", 1, TRUE) ~ "1",
      .data[[outcome]] %in% c("0", "FALSE", "F", "false", 0, FALSE) ~ "0",
      TRUE ~ NA_character_
    ),
    !!outcome := factor(.data[[outcome]], levels = c("0", "1"))
  ) %>%
  filter(!is.na(.data[[outcome]])) %>%
  arrange(.data[[time_col]]) %>%
  select(-week_start)

#### USER LEVEL SPLIT ####
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

write_csv(el_data, file = "/Users/karstenwalker/Documents/Datasets/el_training_week_3_10_sticky.csv")
rm(el_data)

# Class weight
spw <- sum(train_data[[outcome]] == "1") / sum(train_data[[outcome]] == "0")

#### 3) Recipe ####
clf_recipe <- recipe(formula(paste(outcome, "~ .")), data = train_data) %>%
  update_role(all_of(c(id_cols, time_col)), new_role = "id") %>%
  step_unknown(all_nominal_predictors(), new_level = "unknown") %>%
  step_novel(all_nominal_predictors()) %>%
  step_other(all_nominal_predictors(), threshold = 0.01, other = "other_level") %>%
  step_dummy(all_nominal_predictors(), one_hot = TRUE) %>%
  step_zv(all_predictors())

#### 4) Fast model spec for CV ####
cores_to_use <- max(1, parallel::detectCores(logical = FALSE) - 1)

xgb_spec_cv <- boost_tree(
  trees = 500,
  learn_rate = 0.05,
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

clf_wf_cv <- workflow() %>%
  add_recipe(clf_recipe) %>%
  add_model(xgb_spec_cv)

#### 5) 5-fold CV ####
folds <- vfold_cv(train_data, v = 5, strata = !!sym(outcome))
metrics <- metric_set(roc_auc, pr_auc, mn_log_loss)
ctrl <- control_resamples(save_pred = TRUE, parallel_over = "everything", verbose = TRUE)

res_cv <- fit_resamples(
  clf_wf_cv,
  resamples = folds,
  metrics = metrics,
  control = ctrl
)

cv_metrics <- collect_metrics(res_cv)
print(cv_metrics)

cv_preds <- collect_predictions(res_cv)
cv_pr_pooled <- pr_auc(cv_preds, truth = !!sym(outcome), .pred_1)
cv_roc_pooled <- roc_auc(cv_preds, truth = !!sym(outcome), .pred_1)
print(cv_pr_pooled)
print(cv_roc_pooled)

#### 6) Final fit on TRAIN + VALIDATION ####
train_val_data <- bind_rows(train_data, validation_data)
final_wf_fast <- clf_wf_cv
final_fit_fast <- fit(final_wf_fast, data = train_val_data)

#### 7) Evaluate on TEST ####
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
print(test_metrics)

#### Diagnostics ####
engine_fit <- workflows::extract_fit_engine(final_fit_fast)

vip::vip(engine_fit, num_features = 50) +
  labs(title = "28D Retention Model Feature Importance") +
  theme_minimal()

conf_mat(
  test_scored,
  truth = !!sym(outcome),
  estimate = .pred_class
)

roc_df <- roc_curve(
  test_scored,
  truth = !!sym(outcome),
  .data[[prob_col]],
  event_level = event_level
)

pr_df <- pr_curve(
  test_scored,
  truth = !!sym(outcome),
  .data[[prob_col]],
  event_level = event_level
)

ggplot(roc_df, aes(x = 1 - specificity, y = sensitivity)) +
  geom_path(linewidth = 1) +
  geom_abline(linetype = 2) +
  coord_equal() +
  labs(title = "ROC Curve (Test)", x = "1 - Specificity", y = "Sensitivity") +
  theme_minimal()

ggplot(pr_df, aes(x = recall, y = precision)) +
  geom_path(linewidth = 1) +
  labs(title = "Precision-Recall Curve (Test)", x = "Recall", y = "Precision") +
  theme_minimal()

cal_df <- test_scored %>%
  mutate(
    p = .data[[prob_col]],
    y = as.integer(.data[[outcome]] == "1"),
    bin = ntile(p, 10)
  ) %>%
  group_by(bin) %>%
  summarise(
    p_mean = mean(p, na.rm = TRUE),
    y_rate = mean(y, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )

ggplot(cal_df, aes(x = p_mean, y = y_rate)) +
  geom_point(aes(size = n)) +
  geom_abline(linetype = 2) +
  labs(
    title = "Calibration Plot (Test)",
    x = "Mean Predicted Probability (by decile)",
    y = "Observed Event Rate"
  ) +
  theme_minimal()

test_scored %>%
  arrange(desc(.pred_1)) %>%
  mutate(rank = row_number(), pct = rank / n(), decile = ntile(.pred_1, 10)) %>%
  group_by(decile) %>%
  summarise(
    retention = mean(.data[[outcome]] == "1"),
    n = n(),
    .groups = "drop"
  ) %>%
  print()

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

#### Save artifacts ####
artifact_path <- "/Users/karstenwalker/Documents/Modeling/Artifacts"
dir.create(artifact_path, showWarnings = FALSE, recursive = TRUE)

saveRDS(final_fit_fast, file = file.path(artifact_path, "el_fit28_cv.rds"))
saveRDS(res_cv, file = file.path(artifact_path, "el_resamples28_cv.rds"))
write_csv(cv_metrics, file = file.path(artifact_path, "el28_cv_metrics.csv"))
write_csv(test_metrics, file = file.path(artifact_path, "el28_test_metrics_cv.csv"))
write_csv(test_scored, file = file.path(artifact_path, "el28_test_scored_cv.csv"))
write_csv(train_data, file = file.path(artifact_path, "el28_train_cv.csv"))
write_csv(validation_data, file = file.path(artifact_path, "el28_val_cv.csv"))
write_csv(train_val_data, file = file.path(artifact_path, "el28_train_val_cv.csv"))

# Define predictor columns for raw workflow inputs (exclude outcome + id + time)
x_cols <- setdiff(names(train_val_data), c(outcome, id_cols, time_col))

model_meta <- list(
  outcome     = outcome,
  prob_col    = prob_col,
  event_level = event_level,
  id_cols     = id_cols,
  time_col    = time_col,
  x_cols      = x_cols
)

saveRDS(model_meta, file = file.path(artifact_path, "model_meta28_cv.rds"))

# Save SHAP background + observation samples from train+validation
source_df <- train_val_data %>%
  mutate(
    !!time_col := as.Date(.data[[time_col]]),
    !!outcome := factor(.data[[outcome]], levels = c("0", "1"))
  )

shap_bg <- source_df %>%
  select(all_of(x_cols)) %>%
  slice_sample(n = min(3000, n()))

shap_obs <- source_df %>%
  select(all_of(x_cols)) %>%
  slice_sample(n = min(2000, n()))

write_csv(shap_bg, file = file.path(artifact_path, "shap_bg28_cv.csv"))
write_csv(shap_obs, file = file.path(artifact_path, "shap_obs28_cv.csv"))
write_csv(test_data, file = file.path(artifact_path, "el28_test_raw_cv.csv"))

#### Optional reload block for SHAP / PDP ####
final_fit_fast <- readRDS(file.path(artifact_path, "el_fit28_cv.rds"))
model_meta <- readRDS(file.path(artifact_path, "model_meta28_cv.rds"))
train_val_data <- read_csv(file.path(artifact_path, "el28_train_val_cv.csv"), show_col_types = FALSE)
test_data <- read_csv(file.path(artifact_path, "el28_test_raw_cv.csv"), show_col_types = FALSE)
shap_bg <- read_csv(file.path(artifact_path, "shap_bg28_cv.csv"), show_col_types = FALSE)
shap_obs <- read_csv(file.path(artifact_path, "shap_obs28_cv.csv"), show_col_types = FALSE)

outcome     <- model_meta$outcome
prob_col    <- model_meta$prob_col
event_level <- model_meta$event_level
id_cols     <- model_meta$id_cols
time_col    <- model_meta$time_col
x_cols      <- model_meta$x_cols

train_val_data <- train_val_data %>%
  mutate(
    !!time_col := as.Date(.data[[time_col]]),
    !!outcome := factor(.data[[outcome]], levels = c("0", "1"))
  )

test_data <- test_data %>%
  mutate(
    !!time_col := as.Date(.data[[time_col]]),
    !!outcome := factor(.data[[outcome]], levels = c("0", "1"))
  )

stopifnot(identical(colnames(shap_bg), colnames(shap_obs)))

#### SHAP + DALEX using final CV-trained model ####
train_x <- train_val_data %>% select(all_of(x_cols))
train_y <- train_val_data %>% pull(!!sym(outcome))

pred_fun_prob1 <- function(model, newdata) {
  p <- predict(model, new_data = newdata, type = "prob")
  as.numeric(p[[prob_col]])
}

expl <- DALEX::explain(
  model = final_fit_fast,
  data = train_x,
  y = as.numeric(train_y == "1"),
  predict_function = pred_fun_prob1,
  label = "retention_model_cv",
  verbose = FALSE
)

shap_global <- predict_parts(
  explainer = expl,
  new_observation = shap_obs,
  type = "shap",
  B = 25
)

plot(shap_global, show_boxplots = TRUE, max_vars = 20) +
  ggplot2::labs(title = "Global SHAP (sampled observations)") +
  theme_minimal()

example_id <- test_data[[id_cols[1]]][1]

one_obs <- test_data %>%
  filter(.data[[id_cols[1]]] == example_id) %>%
  slice(1) %>%
  select(all_of(x_cols))

shap_local <- predict_parts(
  explainer = expl,
  new_observation = one_obs,
  type = "shap",
  B = 50
)

plot(shap_local, max_vars = 30) +
  ggplot2::labs(title = paste("Local SHAP for", example_id))

#### PDP ####
vars <- c(
  "learn_mode_answer_rate_7d_avg",
  "sets_studied_rate_7d_avg",
  "study_modes_used_7d_avg"
)

vars <- vars[vars %in% x_cols]

pdp <- model_profile(
  explainer = expl,
  variables = vars,
  type = "partial"
)

plot(pdp) +
  ggtitle("Partial Dependence Profiles") +
  theme_minimal()

if ("learn_mode_answer_rate_7d_avg" %in% names(train_val_data)) {
  ggplot(train_val_data, aes(x = learn_mode_answer_rate_7d_avg)) +
    geom_histogram(bins = 30) +
    theme_minimal()
}

#### Permutation importance stability ####
B_perm <- 10
R_reps <- 10

imp_reps <- map(1:R_reps, function(r) {
  set.seed(100 + r)
  idx <- sample.int(nrow(shap_bg), size = min(2000, nrow(shap_bg)))
  expl_r <- expl
  expl_r$data <- shap_bg[idx, , drop = FALSE]
  DALEX::model_parts(expl_r, B = B_perm, type = "difference")
})

parts_df <- bind_rows(imp_reps, .id = "rep") %>%
  filter(!variable %in% c("_baseline_", "_full_model_"))

stability <- parts_df %>%
  group_by(variable) %>%
  summarise(
    mean_drop = mean(dropout_loss, na.rm = TRUE),
    sd_drop = sd(dropout_loss, na.rm = TRUE),
    cv = sd_drop / pmax(abs(mean_drop), 1e-9),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_drop))

print(stability %>% slice_head(n = 30))

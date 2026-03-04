library(tidyverse)
library(tidymodels)
library(doParallel)
library(DALEX)
library(dplyr)
library(arrow)

# ---------------------------
# 0) Inputs
# ---------------------------
# Store your Google Cloud project ID
bq_auth()

project_id <- "compact-sylph-785"

# Define your SQL query (example using a public dataset)
sql <- "SELECT * 
        FROM `compact-sylph-785.Karsten.EL_1pct_1yr_LI_sum_creator`
        WHERE MOD(ABS(FARM_FINGERPRINT(CAST(user_id AS STRING))), 100) < 75"

# Run the query
tb <- bq_project_query(project_id, sql)

# Load helpers

source("/Users/karstenwalker/Documents/GitHub/R/Helpers/themes.R")

source("/Users/karstenwalker/Documents/GitHub/R/Helpers/df_transformations.R")

# Set Seed
set.seed(7)

# Download the results into an R data frame
el_data <- bq_table_download(tb)%>%
  select(-uid, -reported_user_type, -country_code, -retained7d, -retained90d, 
         -days_until_next_session, -cumulative_lifetime_sessions, -age, -platform, -account_type,
         -is_teacher_with_students, -session_count, -first_session_channel, -last_session_channel,
         -lifetime_sets_created,-never_created,-session_count_7d_avg, -set_pageviews_count_7d_avg, 
         -had_first_session, -unique_sets_viewed, -set_pageviews_count, -course_count)

# Optional: read in saved file
el_data<- read.csv(file="/Users/karstenwalker/Documents/Datasets/el_training_no_avg_3_3.csv")

# Note creater_of_content is defined as "User viewed their own content that day"
# Remove course count since it is not as of that day

# Set these
outcome  <- "retained7d"

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
      .data[[outcome]] %in% c("1", "TRUE", "T", "true", 1, TRUE)  ~ "1",
      .data[[outcome]] %in% c("0", "FALSE", "F", "false", 0, FALSE) ~ "0",
      TRUE ~ NA_character_
    ),
    !!outcome := factor(.data[[outcome]], levels = c("0", "1"))
  ) %>%
  filter(!is.na(.data[[outcome]])) %>%   # drop censored days
  arrange(.data[[time_col]])

# 2) Time-aware 80/10/10 split 
split_1 <- rsample::initial_time_split(el_data, prop = 0.80, lag = 0)

train_data <- training(split_1)

tmp_20<- testing(split_1)

split_2 <- rsample::initial_time_split(tmp_20, prop = 0.50, lag = 0)

validation_data <- training(split_2)

test_data <- testing(split_2)

train_data<-train_data%>%
  select(-creator_of_content)

test_data<-test_data%>%
  select(-creator_of_content)

rm(tmp_20)

#### USER LEVEL SPLIT ####
# 1) Sample users into train/val/test (80/10/10)
user_split <- initial_split(
  el_data %>% distinct(user_id),
  prop = 0.80
)

train_users <- training(user_split)$user_id
temp_users  <- testing(user_split)$user_id

set.seed(7)
val_test_split <- initial_split(
  tibble(user_id = temp_users),
  prop = 0.50
)

val_users  <- training(val_test_split)$user_id
test_users <- testing(val_test_split)$user_id

# 2) Filter full data by user sets
train_data <- el_data %>% filter(user_id %in% train_users)
validation_data <- el_data %>% filter(user_id %in% val_users)
test_data  <- el_data %>% filter(user_id %in% test_users)

# 3) Drop columns evenly
train_data <- train_data %>% select(-contains("_7d"), -creator_of_content)
validation_data <- validation_data %>% select(-contains("_7d"), -creator_of_content)
test_data <- test_data %>% select(-contains("_7d"), -creator_of_content)

# Sanity check
range(train_data[[time_col]], na.rm = TRUE)

range(validation_data[[time_col]], na.rm = TRUE)

range(test_data[[time_col]], na.rm = TRUE)

# Drop original data to save memory
write.csv(el_data, file="/Users/karstenwalker/Documents/Datasets/el_training_no_avg_3_3.csv", row.names = FALSE)

rm(el_data)

# 3) Recipe (boosted-tree friendly; avoid unnecessary transforms)
# ---------------------------
clf_recipe <- recipe(formula(paste(outcome, "~ .")), data = train_data) %>%
  update_role(all_of(c(id_cols, time_col)), new_role = "id") %>%
  step_other(all_nominal_predictors(), threshold = 0.01, other = "other_level") %>%
  step_unknown(all_nominal_predictors(), new_level = "unknown") %>%
  step_novel(all_nominal_predictors()) %>%
  step_dummy(all_nominal_predictors(), one_hot = TRUE) %>%
  step_zv(all_predictors())

# ---------------------------
# 4) XGBoost model (fixed params, no tuning)
# ---------------------------
cores_to_use <- max(1, parallel::detectCores(logical = FALSE) - 1)

xgb_spec <- boost_tree(
  trees = 400,          
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
    counts = FALSE,     
    nthread = cores_to_use, 
    verbose = 1
  )

clf_wf <- workflow() %>%
  add_recipe(clf_recipe) %>%
  add_model(xgb_spec)

fit_train <- fit(clf_wf, data = train_data)

# 5) Train on TRAIN, evaluate on test
event_level <- "second"
prob_col <- ".pred_1"

# Score test set
test_probs <- predict(fit_train, new_data = test_data, type = "prob")

test_class <- predict(fit_train, new_data = test_data, type = "class")

test_scored <- bind_cols(
  test_data %>% select(all_of(c(id_cols, time_col, outcome))),
  test_class,   
  test_probs    
)

# Ensure outcome is a factor with levels c("0","1")
test_scored <- test_scored %>%
  mutate(!!outcome := factor(.data[[outcome]], levels = c("0", "1")))

# Class-based metrics 
metrics_class <- metric_set(accuracy)(
  test_scored,
  truth = !!sym(outcome),
  estimate = .pred_class
)

# Prob-based metrics
metrics_prob <- metric_set(roc_auc, pr_auc, mn_log_loss)(
  test_scored,
  truth = !!sym(outcome),
  .data[[prob_col]],
  event_level = event_level
)

test_metrics <- bind_rows(metrics_class, metrics_prob)

test_metrics

# 1st run
# 1 accuracy    binary         0.702
# 2 roc_auc     binary         0.699
# 3 pr_auc      binary         0.828
# 4 mn_log_loss binary         0.570

# 2nd run with more creator flags
# 1 accuracy    binary         0.703
# 2 roc_auc     binary         0.700
# 3 pr_auc      binary         0.829
# 4 mn_log_loss binary         0.570

# 3rd run no averages
# 1 accuracy    binary         0.691
# 2 roc_auc     binary         0.686
# 3 pr_auc      binary         0.823
# 4 mn_log_loss binary         0.577

# 3rd run remove 7d avgs

# 3) Confusion matrix
conf_mat(
  test_scored,
  truth = !!sym(outcome),
  estimate = .pred_class
)

# 4) ROC + PR curves
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
  labs(title = "ROC Curve (Test)", x = "1 - Specificity", y = "Sensitivity")+
  theme_minimal()

ggplot(pr_df, aes(x = recall, y = precision)) +
  geom_path(linewidth = 1) +
  labs(title = "Precision–Recall Curve (Test)", x = "Recall", y = "Precision")+
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
  labs(title = "Calibration Plot (Test)", x = "Mean Predicted Probability (by decile)", y = "Observed Event Rate")+
  theme_minimal()

# Feature importance
engine_fit <- workflows::extract_fit_engine(fit_train)

vip::vip(engine_fit, num_features = 50) +
  labs(title = "Feature Importance (XGBoost - Gain)")+
  theme_minimal()

# Save objects to remove and preserve memory

artifact_path <- "/Users/karstenwalker/Documents/Modeling/Artifacts"

dir.create(artifact_path, showWarnings = FALSE, recursive = TRUE)

dir.create("artifacts", showWarnings = FALSE)

# Save fitted workflow / model
saveRDS(
  fit_train,
  file = file.path(artifact_path, "fit_train.rds")
)

# Save metadata you'll need later
model_meta <- list(
  outcome     = outcome,
  prob_col    = prob_col,
  event_level = event_level,
  id_cols     = id_cols,
  time_col    = time_col
)

saveRDS(
  model_meta,
  file = file.path(artifact_path, "model_meta.rds")
)

write_parquet(
  test_scored,
  sink = file.path(artifact_path, "test_scored.parquet")
)

# RDS version
saveRDS(
  test_scored,
  file = file.path(artifact_path, "test_scored.rds")
)

# Remove things to save memory
rm(
  el_data,
  train_data,
  test_data,
  test_scored,
  test_probs,
  test_class
)

gc()

# Optional: Pre-save only what SHAP needs

x_cols <- setdiff(
  names(test_scored),
  c(model_meta$outcome, ".pred_class", ".pred_0", ".pred_1")
)

shap_df <- test_scored %>%
  select(all_of(x_cols)) %>%
  slice_sample(n = min(3000, n()))

write_parquet(
  shap_df,
  file.path(artifact_path, "shap_df.parquet")
)

# Only data SHAP needs
fit_train  <- readRDS(file.path(artifact_path, "fit_train.rds"))

model_meta <- readRDS(file.path(artifact_path, "model_meta.rds"))

test_scored <- read_parquet(
  file.path(artifact_path, "test_scored.parquet")
)

# Load only SHAP shit
fit_train <- readRDS(file.path(artifact_path, "fit_train.rds"))
shap_df   <- read_parquet(file.path(artifact_path, "shap_df.parquet"))

# SHAP + DALEX

event_level <- "second" 
prob_col <- ".pred_1"
outcome <- outcome  

# Pick the feature columns used by the model
x_cols <- setdiff(names(train_data), outcome)

# Ensure outcome is factor with levels c("0","1")
train_x <- train_data %>%
  mutate(!!outcome := factor(.data[[outcome]], levels = c("0","1"))) %>%
  select(all_of(x_cols))

train_y <- train_data %>%
  mutate(!!outcome := factor(.data[[outcome]], levels = c("0","1"))) %>%
  pull(!!sym(outcome))

# Define a predict function that returns numeric probability of class "1"
pred_fun_prob1 <- function(model, newdata) {
  p <- predict(model, new_data = newdata, type = "prob")
  as.numeric(p[[prob_col]])
}

# DALEX is happiest with numeric y for prob models
expl <- DALEX::explain(
  model = fit_train,
  data  = train_x,
  y     = as.numeric(train_y == "1"),   
  predict_function = pred_fun_prob1,
  label = "retention_model",
  verbose = FALSE
)

# Global SHAP sample on subset to save memory
# B=5 is fine for a quick look, but I’d use B=25–50 on a sample.
shap_n <- 2000

shap_rows <- sample.int(nrow(train_x), size = min(shap_n, nrow(train_x)))

shap_global <- predict_parts(
  explainer = expl,
  new_observation = train_x[shap_rows, , drop = FALSE],
  type = "shap",
  B = 25)

plot(shap_global, show_boxplots = TRUE, max_vars = 20) +
  ggplot2::labs(title = "Global SHAP (sampled observations)")+
  theme_minimal()

# Local SHAP
example_id <- test_data$user_id[1]  #replace with any individual ID

x_test <- test_data %>%
  mutate(!!outcome := factor(.data[[outcome]], levels = c("0","1"))) %>%
  select(all_of(x_cols))

one_obs <- x_test %>% filter(user_id == example_id) %>% slice(1)

shap_local <- predict_parts(
  explainer = expl,
  new_observation = one_obs,
  type = "shap",
  B = 50
)

plot(shap_local, max_vars = 30) +
  ggplot2::labs(title = paste("Local SHAP for", example_id))

# PDP
# pick 1–3 important numeric features
vars <- c("flashcards_answers_7d", "answer_view_ratio_7d")  # example names

pdp <- model_profile(
  explainer = expl,
  variables = vars,
  type = "partial"   # PDP
)

plot(pdp) + ggplot2::labs(title = "Partial Dependence Profiles")

# Feature importance stability (Permutation Importance)

cv_parts <- purrr::map(cv_models, ~{
  expl_cv <- DALEX::explain(
    model = .x,
    data  = train_x,
    y     = as.numeric(train_y == "1"),
    predict_function = pred_fun_prob1,
    label = "cv",
    verbose = FALSE
  )
  DALEX::model_parts(expl_cv, B = 10)   # permutation importance
})

# Combine + summarize stability
parts_df <- bind_rows(cv_parts, .id = "fold") %>%
  filter(variable != "_baseline_" & variable != "_full_model_")

stability <- parts_df %>%
  group_by(variable) %>%
  summarise(
    mean_drop = mean(dropout_loss, na.rm = TRUE),
    sd_drop   = sd(dropout_loss, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_drop))

stability %>% slice_head(n = 30)
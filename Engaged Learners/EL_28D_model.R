library(dplyr)
library(lubridate)
library(readr)

# Load data
el_data_sub <- read_csv("el_data_sub.csv", show_col_types = FALSE)
user_clusters_3_4 <- read_csv("user_clusters_3_4.csv", show_col_types = FALSE)

# Ensure Date types
el_data_sub <- el_data_sub %>%
  mutate(
    date = as.Date(date),
    week_start = floor_date(date, "week", week_start = 1)  # Monday start
  )

user_clusters_3_4 <- df_clustered_uw %>%
  mutate(
    week_start_date = as.Date(week_start_date),
    week_start = floor_date(week_start_date, "week", week_start = 1)
  ) %>%
  select(user_id, week_start, cluster)

# Check for duplicate weekly clusters per user
dup_check <- user_clusters_3_4 %>%
  count(user_id, week_start) %>%
  filter(n > 1)

cat("Duplicate cluster rows per user-week:", nrow(dup_check), "\n")

# Perform LEFT JOIN (no imputation)
joined <- el_data_sub %>%
  left_join(user_clusters_3_4,
            by = c("user_id", "week_start"))

# ---- Failure Summary ----

total_rows <- nrow(joined)
failed_rows <- sum(is.na(joined$cluster))
success_rows <- total_rows - failed_rows

cat("Total rows: ", total_rows, "\n")
cat("Successful joins: ", success_rows, "\n")
cat("Failed joins (NA cluster): ", failed_rows, "\n")
cat("Failure rate: ", round(100 * failed_rows / total_rows, 2), "%\n")

# ---- Breakdown by user ----

user_failures <- joined %>%
  group_by(user_id) %>%
  summarise(
    total_days = n(),
    failed_days = sum(is.na(cluster)),
    failure_rate = failed_days / total_days
  ) %>%
  arrange(desc(failed_days))

# Top 10 users with most failed joins
head(user_failures, 10)

# ---- Breakdown by week ----

week_failures <- joined %>%
  group_by(week_start) %>%
  summarise(
    total_days = n(),
    failed_days = sum(is.na(cluster)),
    failure_rate = failed_days / total_days
  ) %>%
  arrange(desc(failed_days))

head(week_failures, 10)

# ---- Sanity check: users with zero cluster coverage ----

users_no_cluster <- joined %>%
  group_by(user_id) %>%
  summarise(all_missing = all(is.na(cluster))) %>%
  filter(all_missing)

cat("Users with NO cluster coverage at all:", nrow(users_no_cluster), "\n")

el_final <- joined %>%
  arrange(user_id, date) %>%
  group_by(user_id) %>%
  mutate(cluster= tidyr::fill(tibble(cluster), cluster, .direction = "down")$cluster,
         cluster = tidyr::fill(tibble(cluster), cluster, .direction = "up")$cluster
  ) %>%
  ungroup() %>%
  mutate(
    cluster = if_else(is.na(cluster), -1L, as.integer(cluster))
  )


total_rows <- nrow(el_final)
failed_rows <- sum(is.na(el_final$cluster))
success_rows <- total_rows - failed_rows

cat("Total rows: ", total_rows, "\n")
cat("Successful joins: ", success_rows, "\n")
cat("Failed joins (NA cluster): ", failed_rows, "\n")
cat("Failure rate: ", round(100 * failed_rows / total_rows, 2), "%\n")

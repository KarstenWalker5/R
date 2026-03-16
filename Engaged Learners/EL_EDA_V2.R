library(dplyr)
library(ggplot2)
library(mgcv)
library(readr)
library(scales)
library(bigrquery)
library(tidyverse)
library(irlba)
library(ClusterR)
library(openxlsx)
library(stringr)
library(purrr)
library(broom)
library(corrr)

# Store your Google Cloud project ID
bq_auth()

project_id <- "compact-sylph-785"

# Define your SQL query (example using a public dataset)
sql <- "SELECT * 
        FROM `compact-sylph-785.Karsten.EL_1pct_1yr_LI_dev2`"

# Run the query
tb <- bq_project_query(project_id, sql)

# Load helpers

source("/Users/karstenwalker/Documents/GitHub/R/Helpers/themes.R")

source("/Users/karstenwalker/Documents/GitHub/R/Helpers/glm_summaries.R")

source("/Users/karstenwalker/Documents/GitHub/R/Helpers/df_transformations.R")

source("/Users/karstenwalker/Documents/GitHub/R/Helpers/df_transformations.R")

# Set Seed
set.seed(7)

# Download the results into an R data frame
el_data <- bq_table_download(tb)%>%
  select(-uid, -reported_user_type, -country_code,
         -cumulative_lifetime_sessions,-is_teacher_with_students,
         -lifetime_sets_created, -contains("engaged"))%>%
  group_by(account_type)%>%
  mutate(z_score = scale(flashcards_questions_answered)) %>%
  filter(abs(z_score)<3)%>%
  ungroup()

# %>%
#   filter(account_type!="Admin")%>%
#   select(-social_science_courses, -science_courses, -arts_and_humanities_courses,
#          -math_courses, -language_courses, -not_found_courses, -other_courses,
#          -primary_subject_courses, -new_social_science_courses_in_week, -new_science_courses_in_week,
#          -new_arts_and_humanities_courses_in_week, -new_math_courses_in_week, -new_languages_courses_in_week,
#          -new_not_found_courses_in_week, -new_other_courses_in_week)%>%
#   group_by(account_type) %>%
#   mutate(z_score = scale(flashcards_questions_answered)) %>%
#   filter(abs(z_score)<3,flashcards_questions_answered>0)%>%
#   ungroup()
# 
# el_data <- el_data%>%
#   select(-social_science_courses, -science_courses, -arts_and_humanities_courses,
#          -math_courses, -languages_courses, -not_found_courses, -other_courses,
#          -primary_subject_category, -new_social_science_courses_in_week, -new_science_courses_in_week,
#          -new_arts_and_humanities_courses_in_week, -new_math_courses_in_week, -new_languages_courses_in_week,
#          -new_not_found_courses_in_week, -new_other_courses_in_week)

# Count NA by column
na_count <- el_data %>% 
  summarise(across(everything(), ~ sum(is.na(.))))

# How many users only have 1 sesh?
nrow(el_data%>%
       group_by(user_id)%>%
       summarize(total_sessions=sum(session_count))%>%
       ungroup()%>%
       filter(total_sessions>1))

# Distinct user_id should match
n_distinct(el_data$user_id)

el_data%>%
  filter(country_code=="US")%>%
  group_by(date, account_type)%>%
  summarize(avg_sessions=mean(session_count),
            avg_answers=mean(flashcards_questions_answered))%>%
  ungroup()%>%
  ggplot(aes(x=date, y=avg_answers, group=account_type, color=account_type))+
  geom_line()+
  geom_smooth()+
  theme_fancy()

el_data%>%
  filter(country_code=="US")%>%
  group_by(account_type) %>%
  mutate(z_score = scale(flashcards_questions_answered)) %>%
  filter(abs(z_score)<3) %>%
  ungroup()%>%
  ggplot(aes(flashcards_questions_answered,color=account_type))+
  geom_density()+
  theme_fancy()

# Retention buckets
el_data%>%
  group_by(account_type, retained7d) %>%
  summarize(users_7d=n_distinct(user_id))%>%
  rename(retained=retained7d)%>%
  left_join(el_data%>%
              group_by(account_type, retained28d) %>%
              summarize(users_28d=n_distinct(user_id))%>%
              rename(retained=retained28d), by=c("account_type", "retained"))%>%
  left_join(el_data%>%
              group_by(account_type, retained90d) %>%
              summarize(users_90d=n_distinct(user_id))%>%
              rename(retained=retained90d), by=c("account_type", "retained"))%>%
  left_join(el_data%>%
              group_by(account_type, retained28d_sticky) %>%
              summarize(users_28d_sticky=n_distinct(user_id))%>%
              rename(retained=retained28d_sticky), by=c("account_type", "retained"))%>%
  left_join(el_data%>%
              group_by(account_type, retained90d_sticky) %>%
              summarize(users_90d_sticky=n_distinct(user_id))%>%
              rename(retained=retained90d_sticky), by=c("account_type", "retained"))

#averages
el_data%>%
  filter(flashcards_questions_answered>0)%>%
  group_by(account_type, retained7d) %>%
  summarize(avg=mean(flashcards_questions_answered),
            med=median(flashcards_questions_answered),
            dev=sd(flashcards_questions_answered))

### Transformations
numeric_vars<-colnames(el_data%>%
                         select(contains("count"), contains("submitted"), contains("created"),
                                contains("answered"), contains("minutes_active"), contains("viewed"),
                                contains("visits"), contains("reveals"), contains("clicks")))

el_ecdf <- add_ecdf_columns(el_data%>%
                              ungroup(), numeric_vars)%>%
  select(-one_of(numeric_vars))

el_logged <- log_transform_skewed(el_data%>%
                                    ungroup()%>%
                                    select(-z_score), skew_threshold = 1)

### Correlations
cor_matrix <- el_data %>%
  filter(flashcards_questions_answered>0)%>%
  select(-user_id, -uid, -country_code, -contains("date"))%>%
  correlate(use = "complete.obs", method = "pearson")

cor_matrix_spear <- el_data %>%
  filter(flashcards_questions_answered>0)%>%
  select(-user_id, -uid, -country_code, -contains("date"))%>%
  correlate(use = "complete.obs", method = "spearman")

cor_log<-el_logged%>%
  filter(flashcards_questions_answered_log>0)%>%
  select(-contains("week"))%>%
  correlate(use = "complete.obs", method = "pearson")

cor_ecdf<-el_ecdf%>%
  select(-user_id, -uid, -contains("date"))%>%
  correlate(use = "complete.obs", method = "pearson")

# convert to long format and filter strong correlations
cor_long <- cor_matrix %>%
  stretch() %>%
  filter(abs(r) > 0.5, x != y)

cor_long

cor_retention<-cor_matrix %>%
  stretch()%>%
  filter(str_detect(x, "retained"))

cor_retention_spear<-cor_matrix %>%
  stretch()%>%
  filter(str_detect(x, "spearman"))

cor_retention_log<-cor_log %>%
  stretch()%>%
  filter(str_detect(x, "retained"))

cor_retention_ecdf<-cor_ecdf %>%
  stretch()%>%
  filter(str_detect(x, "retained"))

###### PCA + Clustering ######
# 1 row per user/wk, probably redundant but a good check
ids <- el_logged%>%
  group_by(user_id, week)%>%
  filter(row_number()==1)%>%
  ungroup()%>%
  select(user_id, week)

# Replace NA
X <- el_logged%>%
  group_by(user_id, week)%>%
  filter(row_number()==1)%>%
  ungroup()%>%
  ungroup() %>%
  select(-contains("retained"), -contains("week"), -year,
         -age) %>%
  select(-user_id) %>%
  select(where(is.numeric)) %>%
  mutate(across(everything(), ~ replace_na(.x, 0)))

# Drop zero-variance (and near-constant) feature columns

keep_sd <- sapply(X, function(col) sd(col) > 0)

X <- X[, keep_sd, drop = FALSE]

# Scale features
X_scaled <- scale(as.matrix(X))

stopifnot(!anyNA(X_scaled))

### PCA
# Fit PCA on a sample for speed, then project full dataset
n_for_pca <- min(nrow(X_scaled), 200000)

idx_pca <- sample(seq_len(nrow(X_scaled)), n_for_pca)

X_pca_fit <- X_scaled[idx_pca, , drop = FALSE]

max_pcs <- min(ncol(X_scaled), 50)

pca_fit <- irlba::prcomp_irlba(X_pca_fit, n = max_pcs, center = FALSE, scale. = FALSE)

var_explained <- (pca_fit$sdev^2) / sum(pca_fit$sdev^2)

cum_var <- cumsum(var_explained)

# Moving this number up delivers more PCs
target_var <- 0.95

n_pcs <- which(cum_var >= target_var)[1]

message(sprintf("Using %d PCs (%.1f%% cumulative variance on PCA sample).",
                n_pcs, 100 * cum_var[n_pcs]))

X_pcs_full <- X_scaled %*% pca_fit$rotation[, 1:n_pcs, drop = FALSE]

stopifnot(!anyNA(X_pcs_full))

# Elbow selection
n_for_elbow <- min(nrow(X_pcs_full), 200000)

idx_elbow <- sample(seq_len(nrow(X_pcs_full)), n_for_elbow)

X_elbow <- X_pcs_full[idx_elbow, , drop = FALSE]

k_grid <- 2:20

# Compute WSS for each k using MiniBatchKmeans on the elbow sample
elbow_df <- purrr::map_dfr(k_grid, function(k) {
  set.seed(1 + k)
  mbk_tmp <- ClusterR::MiniBatchKmeans(
    data = X_elbow,
    clusters = k,
    batch_size = 20000,
    num_init = 5,
    max_iters = 100,
    init_fraction = 1.0,
    initializer = "kmeans++",
    verbose = TRUE
  )
  tibble(k = k, wss = sum(mbk_tmp$WCSS_per_cluster))
})

# Plot elbow curve
ggplot(elbow_df, aes(x = k, y = wss)) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  labs(
    title = "Elbow plot (MiniBatchKmeans on PCA-reduced data)",
    subtitle = paste0("WSS vs k on n=", n_for_elbow, " sampled rows; PCs=", n_pcs,
                      " (", round(100 * cum_var[n_pcs], 1), "% var)"),
    x = "Number of clusters (k)",
    y = "Total within-cluster sum of squares (WSS)"
  )

# Step changes to select K
elbow_df <- elbow_df %>%
  arrange(k) %>%
  mutate(
    delta = wss - lead(wss),
    pct_improvement = delta / wss
  )

# After 6 marginal gains stabilizes at <3%, noise coming to the far right, 5-6 is prob the sweet spot

### Clustering on PCs
k_final <- 6

mbk <- ClusterR::MiniBatchKmeans(
  data = X_pcs_full,
  clusters = k_final,
  batch_size = 30000,   # 10k–50k typical; increase for stability
  num_init = 5,         # increase for stability, decrease for speed
  max_iters = 100,
  init_fraction = 1.0,
  initializer = "kmeans++",
  verbose = TRUE
)

pred <- ClusterR::predict_MBatchKMeans(
  X_pcs_full,
  mbk$centroids,
  fuzzy = FALSE
)

clusters <- if (is.list(pred)) pred$clusters else pred

# Remove stuff no longer needed
rm(X_scaled)
rm(X_pcs_full)
rm(X)

# Build mapping table
df_clustered_uw <- ids %>%
  mutate(cluster = clusters) %>%
  distinct(user_id, week, .keep_all = TRUE)   

# Join back to full dataset
df_full_with_clusters <- el_logged%>%
  ungroup() %>%
  inner_join(df_clustered_uw, by = c("user_id", "week"))

# Save clusters and IDs to join to modeling df
write.csv(df_full_with_clusters%>%
            select(user_id, cluster, week_start_date, week_end_date), file="/Users/karstenwalker/Documents/Modeling/Artifacts/user_clusters_3_10.csv")

write.csv(df_full_with_clusters, file="/Users/karstenwalker/Documents/Modeling/Artifacts/cluster_full_df_3_5.csv")

# Join clusters back to a frame that includes retention
retention_summary <- el_data %>%
  summarise(
    retention7d_pct = mean(retained7d, na.rm = TRUE) * 100,
    retention28d_pct = mean(retained28d, na.rm = TRUE) * 100,
    retention90d_pct = mean(retained90d, na.rm = TRUE) * 100,
    retention28d_sticky_pct = mean(retained28d_sticky, na.rm = TRUE) * 100,
    retention90d_sticky_pct = mean(retained90d_sticky, na.rm = TRUE) * 100
  )

retention_by_cluster <- df_full_with_clusters %>%
  group_by(cluster) %>%
  summarise(
    n = n(),
    retention_rate7d = mean(retained7d, na.rm = TRUE),
    retention_rate28d = mean(retained28d, na.rm = TRUE),
    retention_rate90d = mean(retained90d, na.rm = TRUE),
    retention_rate28d_sticky = mean(retained28d_sticky, na.rm = TRUE),
    retention_rate90d_sticky = mean(retained90d_sticky, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(cluster))

df_full_with_clusters %>%
  summarise(
    n = n(),
    retention_rate7d = mean(retained7d, na.rm = TRUE),
    retention_rate28d = mean(retained28d, na.rm = TRUE),
    retention_rate90d = mean(retained90d, na.rm = TRUE),
    retention_rate28d_sticky = mean(retained28d_sticky, na.rm = TRUE),
    retention_rate90d_sticky = mean(retained90d_sticky, na.rm = TRUE),
    .groups = "drop"
  )

retention_by_cluster

# Lift by cluster
# Calculates the % vs average for the cluster for retention
overall <- df_full_with_clusters %>%
  summarise(
    r7  = mean(retained7d,  na.rm = TRUE),
    r28 = mean(retained28d, na.rm = TRUE),
    r90 = mean(retained90d, na.rm = TRUE),
    r28_sticky = mean(retained28d_sticky, na.rm = TRUE),
    r90_sticky = mean(retained90d_sticky, na.rm = TRUE)
  )

lift_by_cluster<-df_full_with_clusters %>%
  group_by(cluster) %>%
  summarise(
    n = n(),
    r7  = mean(retained7d,  na.rm = TRUE),
    r28 = mean(retained28d, na.rm = TRUE),
    r90 = mean(retained90d, na.rm = TRUE),
    r28_sticky = mean(retained28d_sticky, na.rm = TRUE),
    r90_sticky = mean(retained90d_sticky, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    lift7  = r7  / overall$r7,
    lift28 = r28 / overall$r28,
    lift90 = r90 / overall$r90,
    r28_sticky / overall$r28_sticky,
    r90_sticky / overall$r90_sticky
  ) %>%
  arrange(desc(lift7))

# Plotting
# using existing PCA projection (if available)
plot_x <- X_pcs_full

plot_df <- data.frame(
  PC1 = plot_x[, 1],
  PC2 = plot_x[, 2],
  cluster = factor(clusters)
)

# Sample for plotting 
plot_n <- min(nrow(plot_df), 200000)

plot_df_s <- plot_df%>%
  slice_sample(n = plot_n)

ggplot(plot_df_s, aes(PC1, PC2, color = cluster)) +
  geom_point(alpha = 0.35, size = 0.7) +
  theme_fancy() +
  labs(
    title = "Clusters, PC1 vs. PC2",
    subtitle = paste0("Sample n=", plot_n),
    color = "Cluster"
  )

# Adds elipses
ggplot(plot_df_s, aes(PC1, PC2, color = cluster)) +
  geom_point(alpha = 0.35, size = 0.7) +
  stat_ellipse(aes(fill = cluster),
               geom = "polygon",
               type = "norm",
               level = 0.90,      # ~covers 90% if roughly elliptical/normal
               alpha = 0.15,
               linewidth = 0.6,
               show.legend = FALSE) +
  theme_fancy() +
  labs(
    title = "Clusters, PC1 vs. PC2",
    subtitle = paste0("Sample n=", plot_n),
    color = "Cluster"
  )

### Try without PCA 
mbk2 <- ClusterR::MiniBatchKmeans(
  data = X_scaled,
  clusters = k_final,
  batch_size = 30000,
  num_init = 5,
  max_iters = 100,
  init_fraction = 1.0,
  initializer = "kmeans++",
  verbose = TRUE
)

# Predict hard cluster assignments
pred2 <- ClusterR::predict_MBatchKMeans(
  X_scaled,
  mbk2$centroids,
  fuzzy = FALSE
)

clusters2 <- if (is.list(pred2) && !is.null(pred2$clusters)) pred2$clusters else pred2

df_clustered_uw2 <- ids %>%
  mutate(cluster = clusters2) %>%
  select(user_id, week, cluster)

# Join to subset
df_full_with_clusters2 <- el_sub %>%
  inner_join(df_clustered_uw2, by = c("user_id", "week"))

# Evaluate against retention

retention_by_cluster2 <- df_full_with_clusters2 %>%
  group_by(cluster) %>%
  summarise(
    n = n(),
    retention_rate7d  = mean(retained7d,  na.rm = TRUE),
    retention_rate28d = mean(retained28d, na.rm = TRUE),
    retention_rate90d = mean(retained90d, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(cluster))

retention_by_cluster2

overall <- df_full_with_clusters %>%
  summarise(
    r7  = mean(retained7d,  na.rm = TRUE),
    r28 = mean(retained28d, na.rm = TRUE),
    r90 = mean(retained90d, na.rm = TRUE)
  )

df_full_with_clusters2 %>%
  group_by(cluster) %>%
  summarise(
    n = n(),
    r7  = mean(retained7d,  na.rm = TRUE),
    r28 = mean(retained28d, na.rm = TRUE),
    r90 = mean(retained90d, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    lift7  = r7  / overall$r7,
    lift28 = r28 / overall$r28,
    lift90 = r90 / overall$r90
  ) %>%
  arrange(desc(lift7))

### Compute cluster stats

df_full_with_clusters %>%
  group_by(cluster) %>%
  summarise(
    pct_rows_zero = mean(flashcards_questions_answered_log == 0, na.rm = TRUE) * 100,
    n = n()
  )

# Build a "global mean" named vector (after back-transforming *_log)
global_means <- df_full_with_clusters %>%
  summarise(
    across(
      where(is.numeric) & !contains("retained"),
      ~ if (grepl("_log$", cur_column())) {
        mean(exp(.x) - 1, na.rm = TRUE)
      } else {
        mean(.x, na.rm = TRUE)
      },
      .names = "mean_{sub('_log$', '', .col)}"
    )
  ) %>%
  unlist(use.names = TRUE)

# Cluster means (same transform logic)
cluster_means <- df_full_with_clusters %>%
  group_by(cluster) %>%
  summarise(
    across(
      where(is.numeric) & !contains("retained"),
      ~ if (grepl("_log$", cur_column())) {
        mean(exp(.x) - 1, na.rm = TRUE)
      } else {
        mean(.x, na.rm = TRUE)
      },
      .names = "mean_{sub('_log$', '', .col)}"
    ),
    n = n(),
    .groups = "drop"
  )

# Tidy long table: cluster mean, global mean, % lift vs global
cluster_summary_tidy <- cluster_means %>%
  pivot_longer(
    cols = -c(cluster, n),
    names_to = "metric",
    values_to = "cluster_mean"
  ) %>%
  mutate(
    global_mean = global_means[metric],
    pct_lift_vs_global = dplyr::if_else(
      is.na(global_mean) | global_mean == 0,
      NA_real_,
      100 * (cluster_mean / global_mean - 1)
    )
  ) %>%
  arrange(cluster, desc(abs(pct_lift_vs_global)))

cluster_summary_tidy

# By Teacher/Student
# Deprecated, excluding teachers from clustering
cluster_profile <- df_full_with_clusters %>%
  group_by(cluster) %>%
  summarise(
    n = n(),
    pct_teacher = mean(reported_user_type == "Teacher", na.rm = TRUE),
    pct_paid_any = mean(account_type %in% c("Plus", "PaidTeacher"), na.rm = TRUE),
    pct_free = mean(account_type == "FreeLoggedIn", na.rm = TRUE),
    pct_paid_teacher = mean(account_type == "PaidTeacher", na.rm = TRUE),
    pct_plus = mean(account_type == "Plus", na.rm = TRUE),
    retention_7d = mean(retained7d, na.rm = TRUE),
    retention_28d = mean(retained28d, na.rm = TRUE),
    retention_90d = mean(retained90d, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(retention_7d))

cluster_profile

### Cluster stability insights
# If clusters are "types" we can use as a feature for modeling
# If clusters are states we can as for leading indicators

# How much variance do the clusters contain?
total_var <- sum(apply(X_pcs_full, 2, var))

within_var <- sum(mbk$WCSS_per_cluster) / nrow(X_pcs_full)

between_ratio <- 1 - within_var / total_var

between_ratio
#.99

# Retention variance explained by cluster
model_cluster_only <- glm(retained7d ~ factor(cluster),
                          data = df_full_with_clusters,
                          family = binomial)

summary(model_cluster_only)

# Test cluster stability
cluster_runs <- replicate(5, {
  set.seed(sample(1:10000, 1))
  mbk_tmp <- ClusterR::MiniBatchKmeans(
    data = X_pcs_full,
    clusters = 6,
    batch_size = 20000,
    num_init = 5,
    max_iters = 100,
    init_fraction = 1.0,
    initializer = "kmeans++",
    verbose = FALSE
  )
  pred <- ClusterR::predict_MBatchKMeans(
    X_pcs_full,
    mbk_tmp$centroids,
    fuzzy = FALSE
  )
  if (is.list(pred)) pred$clusters else pred
}, simplify = FALSE)

mclust::adjustedRandIndex(cluster_runs[[1]], cluster_runs[[2]])
# Val of 1 means extremely stable clusters across iteratons using different seeds

# How often do users move between clusters
df_full_with_clusters %>%
  arrange(user_id, week) %>%
  group_by(user_id) %>%
  mutate(next_cluster = lead(cluster)) %>%
  count(cluster, next_cluster)

transitions <- df_full_with_clusters %>%
  arrange(user_id, week) %>%
  group_by(user_id) %>%
  mutate(next_cluster = lead(cluster)) %>%
  ungroup() %>%
  filter(!is.na(next_cluster)) %>%
  count(cluster, next_cluster) %>%
  group_by(cluster) %>%
  mutate(
    total = sum(n),
    transition_prob = n / total
  ) %>%
  ungroup()

transition_matrix <- transitions %>%
  select(cluster, next_cluster, transition_prob) %>%
  pivot_wider(
    names_from = next_cluster,
    values_from = transition_prob,
    values_fill = 0
  )

transition_matrix

# Transition Matrix Summary (P(cluster_t+1 | cluster_t))

# Overall Structure
# * Cluster 5 is the most stable state (~67% stay probability).
# * Clusters 6 and 3 show moderate stability (~50% stay).
# * Cluster 4 is moderately persistent (~45% stay).
# * Cluster 2 has lower stability (~37% stay).
# * Cluster 1 is the most transitional (~29% stay).

# Behavioral Interpretation
# * Cluster 5 represents the dominant low engagement absorbing state.
# * Cluster 6 represents the high engagement anchor state.
# * Clusters 2, 3, and 4 represent mid level engagement tiers.
# * Cluster 1 functions as a volatile redistribution state.

# Major Flow Patterns
# * Clusters 2, 3, 4, and 6 all feed meaningfully into Cluster 5.
# * Cluster 4 transitions to Cluster 5 at ~38%.
# * Cluster 6 transitions to Cluster 5 at ~29%.
# * Cluster 1 redistributes broadly and does not act as an absorbing sink.
# * Very little flow moves back into Cluster 1 from higher states.

# Structural Assessment
# * Clear absorbing low engagement state in Cluster 5.
# * Clear high engagement anchor in Cluster 6.
# * Mid tier states transition predictably toward Cluster 5.
# * System resembles a coherent behavioral lifecycle funnel.

# Strategic Implications
# * Retention is weakest in Cluster 5 and strongest in Cluster 6.
# * Preventing movement into Cluster 5 should improve retention.
# * Cluster 1 users require monitoring due to high volatility.
# * Movement patterns are directionally consistent with retention gradients.

# Conclusion
# * k = 6 is structurally clean and interpretable.
# * The cluster system shows ordered engagement tiers and meaningful lifecycle dynamics.

stability <- transitions %>%
  filter(cluster == next_cluster) %>%
  select(cluster, transition_prob) %>%
  rename(stay_prob = transition_prob)

stability

# * Cluster 5 is the most stable state (~67% persistence).
# * Clusters 6 and 3 show moderate stability (~50% persistence).
# * Cluster 4 has moderate persistence (~45%).
# * Cluster 2 is less stable (~37%).
# * Cluster 1 is the most volatile (~29% persistence).

# * Cluster 5 represents a persistent low engagement state.
# * Cluster 6 represents a stable high engagement state.
# * Clusters 3 and 4 are mid tier states with moderate stickiness.
# * Cluster 2 is semi transitional.
# * Cluster 1 is highly transitional and redistributive.

# How many unique clusters each user occupies across all weeks
user_cluster_span <- df_full_with_clusters %>%
  group_by(user_id) %>%
  summarise(
    n_clusters_visited = n_distinct(cluster),
    .groups = "drop"
  )

# table
table(user_cluster_span$n_clusters_visited)

# * Most users occupy persistent engagement tiers.
# * A minority of users exhibit meaningful behavioral mobility.
# * Cluster transitions likely reflect real engagement change rather than instability.

# At some point most users had not activated
pct_visited<-user_cluster_span %>%
  count(n_clusters_visited) %>%
  mutate(
    pct_users = 100 * n / sum(n)
  ) %>%
  arrange(n_clusters_visited)

pct_visited

# User Cluster Span Summary
# Stability Interpretation
# * The majority of users remain in one or two clusters over time.
# * Only a small fraction visit 4 or more clusters.
# * Very few users traverse the entire state space.

# Structural Assessment
# * Clusters represent relatively stable behavioral regimes.
# * Movement exists but is not excessive.
# * Low incidence of full span traversal suggests limited noise fragmentation.

# Lifecycle Implication
# * Most users occupy persistent engagement tiers.
# * A minority of users exhibit meaningful behavioral mobility.
# * Cluster transitions likely reflect real engagement change rather than instability.

# Conclusion
# * The 6-cluster system shows reasonable longitudinal stability.
# * Behavioral states appear structured rather than randomly fluctuating.

# Movement rate
# < 20% → mostly stable types
# 20–40% → moderate state behavior
# 50% → highly dynamic states
# 28.7%

movement_rate <- transitions %>%
  summarise(
    total_transitions = sum(n),
    moves = sum(n[cluster != next_cluster])
  ) %>%
  mutate(move_rate = moves / total_transitions)

movement_rate

# Heatmap
ggplot(transitions, aes(x = factor(cluster), 
                        y = factor(next_cluster), 
                        fill = transition_prob)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red") +
  labs(
    x = "Cluster at week t",
    y = "Cluster at week t+1",
    fill = "Transition Prob",
    title = "Cluster Transition Matrix"
  ) +
  theme_minimal()

# 37% probability a user leaves a cluster
initial_cluster <- df_full_with_clusters %>%
  arrange(user_id, week) %>%
  group_by(user_id) %>%
  summarise(
    first_cluster = first(cluster),
    ever_changed = any(cluster != first(cluster)),
    .groups = "drop"
  )

mean(initial_cluster$ever_changed)

# Add “prev_cluster” and “transition” features
df_states <- df_full_with_clusters %>%
  arrange(user_id, week) %>%
  group_by(user_id) %>%
  mutate(
    prev_cluster = lag(cluster),
    changed_cluster = as.integer(!is.na(prev_cluster) & cluster != prev_cluster),
    transition = if_else(is.na(prev_cluster), NA_character_,
                         paste0(prev_cluster, "->", cluster))
  ) %>%
  ungroup()

# State vs. retention
state_vs_ret<-df_states %>%
  filter(!is.na(transition)) %>%
  group_by(transition) %>%
  summarise(
    n = n(),
    r7 = mean(retained7d, na.rm = TRUE)
  ) %>%
  arrange(desc(n))

print(state_vs_ret%>%
        arrange(desc(r7)), n=50)

# Transition vs 7D Retention Summary

# What This Code Does
# - Constructs user week-to-week states:
#   - prev_cluster: prior week's cluster
#   - changed_cluster: indicator for cluster_t != cluster_t-1 (excluding first observation)
#   - transition: "prev->current" label for each observed move
# - Aggregates by transition to compute:
#   - n  = number of observed transitions
#   - r7 = mean 7-day retention for that transition

# Destination cluster dominates retention outcomes.
# * All transitions ending in Cluster 6 have the highest retention (0.90 to 0.93).
# * 6->6 is the strongest state with ~93% 7d retention.

# Cluster 2 represents a strong upper mid state.
# * 2->2 shows ~82% retention.
# * Transitions into 2 generally perform well (0.75 to 0.79 range).

# Cluster 4 and 3 represent mid tier states.
# * 4->4 ~70% retention.
# * 3->3 ~65% retention.
# * Transitions among 3 and 4 are moderate but clearly below Cluster 6.

# Cluster 5 is the clear low retention sink.
# * 5->5 ~47% retention with very large volume.
# * All transitions into 5 range from ~39% to ~54%.
# * Cluster 5 materially drives overall churn risk.

# Directional pattern is consistent.
# * Moves into 6 sharply increase retention.
# * Moves into 5 sharply decrease retention.
# * Downward transitions toward 5 show the lowest outcomes.

# Structural Assessment
# * Retention tiers are cleanly ordered by destination cluster.
# * High volume transitions reinforce the lifecycle funnel into Cluster 5.
# * Cluster 6 functions as the high engagement anchor state.

# Simple regression of current cluster vs prev
cluster_move_glm1<-glm(retained7d ~ factor(cluster) + factor(prev_cluster) + changed_cluster,
                       data = df_states,
                       family = binomial)

summary(cluster_move_glm1)

# Baseline = users in cluster 1, previously in cluster 1, who did NOT change clusters.
# Implied baseline retention ≈ plogis(0.72) ≈ 0.67.

# Effect of Current Cluster (dominant driver)
# * Cluster 6: strongly positive (β ≈ +1.46) → dramatically higher retention than cluster 1 (odds ~4.3x).
# * Cluster 2: positive (β ≈ +0.52) → meaningfully higher retention than cluster 1 (odds ~1.7x).
# * Cluster 3: slightly negative (β ≈ -0.13) → modestly lower retention than cluster 1.
# * Cluster 4: slightly negative and marginally significant.
# * Cluster 5: strongly negative (β ≈ -0.99) → substantially lower retention (odds ~0.37x).

# Interpretation:
# * Destination cluster is the strongest predictor of retention.
# * Cluster 6 is the clear high-retention state.
# * Cluster 5 is the clear low-retention state.

# Effect of Previous Cluster
# * Previous cluster 6 has a strong positive effect (β ≈ +0.56).
# * Previous clusters 2, 4, and 5 have moderate positive effects.
# * Previous cluster 3 is not significant.
# * Past state matters, but less than current cluster.

# Effect of Changing Clusters
# * changed_cluster β ≈ -0.16 (highly significant).
# * Odds ratio ≈ exp(-0.16) ≈ 0.85.
# * Changing clusters reduces retention odds by ~15%, holding current and previous clusters constant.
# * Behavioral instability independently lowers retention.

# Overall Conclusion
# * Current cluster is the primary determinant of retention.
# * Cluster 6 drives strong retention, Cluster 5 drives churn risk.
# * Movement between clusters carries a meaningful retention penalty.
# * Retention is largely determined by present state, with prior state and stability adding secondary signal.

# Replace changed_cluster with directional movement
df_states <- df_states %>%
  mutate(
    direction = case_when(
      prev_cluster < cluster ~ "up",
      prev_cluster > cluster ~ "down",
      prev_cluster == cluster ~ "same",
      TRUE ~ NA_character_
    )
  )

cluster_direction_glm<-glm(retained7d ~ factor(cluster) + factor(prev_cluster) + factor(direction),
                           data = df_states,
                           family = binomial)

summary(cluster_direction_glm)
# Baseline = users in cluster 1, previously in cluster 1, who moved down.
# Implied baseline retention ≈ plogis(0.60) ≈ 0.65.

# Current cluster strongly predicts 7 day retention.
# * Cluster 6 has the highest retention (β ≈ +1.54, odds ~4.7x vs cluster 1).
# * Cluster 2 meaningfully increases retention (β ≈ +0.52, odds ~1.7x).
# * Cluster 3 is slightly lower than cluster 1.
# * Cluster 4 is not significantly different from cluster 1.
# * Cluster 5 is substantially lower retention (β ≈ -0.94, odds ~0.39x).

# Previous cluster has smaller but meaningful effects.
# * Previous cluster 6 has a strong positive effect.
# * Previous clusters 2 and 4 increase retention odds.
# * Previous cluster 5 has a small positive effect.
# * Previous cluster 3 is not significant.
# * Current state remains much more important than prior state.

# Compared to users who move down:
# * Staying in the same cluster increases retention odds by ~13% (exp(0.124) ≈ 1.13).
# * Moving up decreases retention odds by ~7% (exp(-0.071) ≈ 0.93) after controlling for destination cluster.

# Direction adds signal, but current destination cluster captures most of the engagement effect.
# Retention is primarily determined by where the user is now, with movement direction providing modest incremental information.

# Interaction model
glm(retained7d ~ factor(prev_cluster) * factor(cluster),
    family = binomial,
    data = df_states)

# Retention is driven by the full transition, not just current cluster or direction.

# Cluster 6 is the highest retention destination overall (strong positive main effect), while Cluster 5 is the lowest (large negative main effect).

# Several transitions materially outperform their main effects.
# * Moves from cluster 2 into cluster 2 show strong positive interaction effects.
# * Moves from cluster 5 into cluster 6 and from cluster 2 into cluster 6 materially exceed what the main cluster effects alone would predict.
# * Some moves into cluster 2 from clusters 3, 4, and 5 also outperform baseline expectations.

# Some transitions underperform relative to their main effects.
# * Certain moves from cluster 6 into clusters 3 and 4 show negative interaction terms.
# * Some transitions into cluster 5 underperform even after accounting for its already low baseline retention.
# * Remaining in cluster 6 from cluster 6 shows a negative interaction relative to its strong main effect, indicating diminishing incremental lift.

# The interaction model shows that not all upward or downward moves behave the same.
# * Specific origin to destination pairs meaningfully alter retention odds.
# * The effect of moving into a cluster depends strongly on where the user came from.

# Large reduction in deviance versus null confirms strong explanatory power.
# Lower AIC relative to simpler models suggests the full transition specification fits better.

# Overall, retention depends on the exact cluster_t to cluster_t+1 transition.
# Modeling the full interaction captures transition specific effects that are masked in simpler direction or change indicators.

# Markov transition matrix
P<- prop.table(table(df_states$prev_cluster,
                     df_states$cluster), 1)

transition_df <- as.data.frame(as.table(P)) %>%
  rename(
    prev_cluster = Var1,
    cluster = Var2,
    prob = Freq
  ) %>%
  mutate(
    prev_cluster = factor(prev_cluster),
    cluster = factor(cluster)
  )

ggplot(transition_df,
       aes(x = cluster,
           y = prev_cluster,
           fill = prob)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", prob)),
            size = 3) +
  scale_fill_gradient(
    low = "#f7fbff",
    high = "red",
    name = "Transition\nProbability"
  ) +
  labs(
    title = "Cluster Transition Matrix",
    subtitle = "P(cluster_t+1 | cluster_t)",
    x = "Current Cluster (t+1)",
    y = "Previous Cluster (t)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(color = "black")
  )

# Clusters 5 and 6 are the most stable states with indicating strong behavioral stickiness.
# Clusters 3 and 4 show moderate stability around 50%, suggesting semi-persistent engagement states.
# Clusters 1 and 2 are the most transitional, with only ~29% to 37% staying, indicating higher volatility.
# Cluster 1 redistributes broadly, especially into clusters 4 and 5, and does not act as a churn sink.
# Cluster 4 frequently transitions into cluster 5, suggesting sequential upward engagement movement.
# Clusters 5 and 6 show limited outward leakage, reinforcing their role as stable high-engagement states.
# Overall, the system reflects a structured behavioral lifecycle with clear stable states and identifiable transitional pathways.

# Retention by current state and by transition
state_ret7 <- df_states %>%
  group_by(cluster) %>%
  summarise(n = n(),
            r7 = mean(retained7d, na.rm = TRUE),
            .groups = "drop")

transition_ret7 <- df_states %>%
  filter(!is.na(transition)) %>%
  group_by(transition) %>%
  summarise(n = n(),
            r7 = mean(retained7d, na.rm = TRUE),
            .groups = "drop") %>%
  arrange(desc(n))

state_ret28 <- df_states %>%
  group_by(cluster) %>%
  summarise(n = n(),
            r28 = mean(retained28d, na.rm = TRUE),
            .groups = "drop")

transition_ret28 <- df_states %>%
  filter(!is.na(transition)) %>%
  group_by(transition) %>%
  summarise(n = n(),
            r28 = mean(retained28d, na.rm = TRUE),
            .groups = "drop") %>%
  arrange(desc(n))

state_ret28_sticky <- df_states %>%
  group_by(cluster) %>%
  summarise(n = n(),
            r28 = mean(retained28d_sticky, na.rm = TRUE),
            .groups = "drop")

transition_ret28_sticky <- df_states %>%
  filter(!is.na(transition)) %>%
  group_by(transition) %>%
  summarise(n = n(),
            r28 = mean(retained28d_sticky, na.rm = TRUE),
            .groups = "drop") %>%
  arrange(desc(n))

# Number of weeks in same cluster
df_runs <- df_states %>%
  arrange(user_id, week) %>%
  group_by(user_id) %>%
  mutate(
    new_run = is.na(prev_cluster) | cluster != prev_cluster,
    run_id = cumsum(new_run)
  ) %>%
  group_by(user_id, run_id) %>%
  summarise(
    cluster = first(cluster),
    run_length_weeks = n(),
    retained7d_mean = mean(retained7d, na.rm = TRUE),
    .groups = "drop"
  )

### Cluster drivers
# Compare
# * 6 vs 5 (retention high/low respectively)
# * 6 vs 2 (high vs. upper mid)
# * 2 vs 5 (upper mid vs. low)
# * 3 vs 4 (middle clusters)
# * 1 vs 5 (lower volatile vs. low sink)

cluster_drivers <- cluster_summary_tidy %>%
  filter(abs(pct_lift_vs_global) > 20) %>%   # only meaningful differences
  group_by(cluster) %>%
  slice_max(abs(pct_lift_vs_global), n = 10) %>%
  arrange(cluster, desc(abs(pct_lift_vs_global)))

# Why is cluster 6 so much higher?
comparison <- cluster_summary_tidy %>%
  filter(cluster %in% c(6, 2)) %>%
  select(cluster, metric, pct_lift_vs_global) %>%
  pivot_wider(names_from = cluster, values_from = pct_lift_vs_global) %>%
  mutate(
    diff_6_vs_2 = `6` - `2`
  )

comparison %>%
  arrange(desc(abs(diff_6_vs_2))) %>%
  slice_head(n = 20)

# function for contrasting clusters so I don't copy/paste this shit
contrast_clusters <- function(df, c1, c2, top_n = 15) {
  df %>%
    filter(cluster %in% c(c1, c2)) %>%
    select(cluster, metric, pct_lift_vs_global) %>%
    pivot_wider(names_from = cluster, values_from = pct_lift_vs_global) %>%
    mutate(diff = .[[as.character(c1)]] - .[[as.character(c2)]]) %>%
    arrange(desc(abs(diff))) %>%
    slice_head(n = top_n)
}

# Check contrasts
contrast_6_2<- contrast_clusters(cluster_summary_tidy, 6, 2)

contrast_6_5<- contrast_clusters(cluster_summary_tidy, 6, 5)

contrast_6_2<-contrast_clusters(cluster_summary_tidy, 6, 2)

contrast_2_5<-contrast_clusters(cluster_summary_tidy, 2, 5)

contrast_3_4<-contrast_clusters(cluster_summary_tidy, 3, 4)

contrast_1_5<- contrast_clusters(cluster_summary_tidy, 1, 5)

# Compare highest retention clusters
high <- c(6, 2, 1)

wide_high <- cluster_summary_tidy %>%
  filter(cluster %in% high) %>%
  select(cluster, metric, pct_lift_vs_global) %>%
  pivot_wider(names_from = cluster, values_from = pct_lift_vs_global)

# Top 3 retention cluster diffs
# Top 3 retention cluster diffs
high_archetype_1 <- wide_high %>%
  mutate(d62 = `6` - `2`, d21 = `2` - `1`) %>%
  filter(d62 > 0, d21 > 0) %>%
  arrange(desc(d62 + d21)) %>%
  slice_head(n = 20)

# Now 4 to 1 contrast
high_archetype_2 <- wide_high %>%
  mutate(d14 = `1` - `4`, d43 = `4` - `3`) %>%
  filter(d14 > 0, d43 > 0) %>%
  arrange(desc(d14 + d43)) %>%
  slice_head(n = 20)

high_archetype_1

high_archetype_2

# Composite score of metric "styles"
# Need to add more items to the list
axis_defs <- list(
  creation = c("sets_created", "manual_sets_created", "class_creation", "folders_created", "items_added_to_folder"),
  consumption = c("questions_viewed", "flashcards_questions_answered", "minutes_active"),
  expert_solutions = c("expert_solutions"),
  organization = c("folders_", "items_added_to_folder", "organize|organization")
)

# Build axis scores using available columns (works with *_log too)
add_axis_scores <- function(df, axis_defs) {
  feature_cols <- get_feature_cols(df)
  for (ax in names(axis_defs)) {
    pats <- axis_defs[[ax]]
    cols <- feature_cols[reduce(pats, ~ .x | str_detect(feature_cols, .y), .init = FALSE)]
    if (length(cols) == 0) next
    df <- df %>%
      mutate(!!paste0("axis_", ax) := rowSums(across(all_of(cols), ~ replace_na(.x, 0)), na.rm = TRUE))
  }
  df
}

df_with_axes <- add_axis_scores(df_full_with_clusters, axis_defs)

axis_profile <- df_with_axes %>%
  group_by(cluster) %>%
  summarise(across(starts_with("axis_"), ~ mean(.x, na.rm = TRUE)),
            n = n(),
            r7 = mean(retained7d, na.rm = TRUE),
            .groups = "drop")

axis_profile

# Identify retention-neutral clusters w/ similar retention, diff behavior
# Find similar r7, show top diffs
ret_by_cluster <- df_full_with_clusters %>%
  group_by(cluster) %>%
  summarise(r7 = mean(retained7d, na.rm = TRUE), n = n(), .groups = "drop")

# Pick pairs within +/- 0.01 retention
tol <- 0.01

pairs <- tidyr::crossing(c1 = ret_by_cluster$cluster, c2 = ret_by_cluster$cluster) %>%
  filter(c1 < c2) %>%
  left_join(ret_by_cluster, by = c("c1" = "cluster")) %>%
  rename(r7_1 = r7, n1 = n) %>%
  left_join(ret_by_cluster, by = c("c2" = "cluster")) %>%
  rename(r7_2 = r7, n2 = n) %>%
  mutate(diff_r7 = abs(r7_1 - r7_2)) %>%
  filter(diff_r7 <= tol) %>%
  arrange(diff_r7)

pairs %>% slice_head(n = 10)

# Show contrasts for the best pair (if any)
best_pair_contrasts<-if (nrow(pairs) > 0) {
  best <- pairs[1, ]
  contrast_clusters(cluster_summary_tidy, best$c1, best$c2, top_n = 20)
}

# Behavioral diversity: feature breadth + entropy
# Breadth = # of features used/on-zero in a week.

feature_cols <- get_feature_cols(df_full_with_clusters)

num_feature_cols <- feature_cols[sapply(df_full_with_clusters[feature_cols], is.numeric)]

df_div <- df_full_with_clusters %>%
  mutate(
    feature_breadth = rowSums(
      across(all_of(num_feature_cols), ~ replace_na(.x, 0) > 0)
    )
  )

breadth_by_cluster <- df_div %>%
  group_by(cluster) %>%
  summarise(
    n = n(),
    mean_breadth = mean(feature_breadth, na.rm = TRUE),
    r7 = mean(retained7d, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_breadth))

breadth_by_cluster

# Entropy
# identify feature columns automatically
feature_cols <- df_full_with_clusters %>%
  select(-user_id, -week, -cluster) %>%
  select(where(is.numeric)) %>%
  select(-starts_with("retained")) %>%
  colnames()

# select count vars only
count_like <- feature_cols[
  str_detect(feature_cols, 
             "count|minutes|questions|views|clicks|created|added|reveals")
]

# eEntropy func
entropy_row <- function(x) {
  x <- pmax(x, 0)
  s <- sum(x)
  if (s == 0) return(NA_real_)
  p <- x / s
  -sum(p[p > 0] * log(p[p > 0]))
}

# compute entropy per user-week
df_entropy <- df_full_with_clusters %>%
  mutate(
    entropy_usage = apply(
      select(., all_of(count_like)) %>%
        mutate(across(everything(), ~ replace_na(.x, 0))) %>%
        as.matrix(),
      1,
      entropy_row
    )
  )

# compare by cluster
entropy_compare<-df_entropy %>%
  group_by(cluster) %>%
  summarise(
    n = n(),
    mean_entropy = mean(entropy_usage, na.rm = TRUE),
    r7 = mean(retained7d, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_entropy))

entropy_compare

# Cluster output table
cluster_fact_sheet <- df_full_with_clusters %>%
  group_by(cluster) %>%
  summarise(
    n = n(),
    r7 = mean(retained7d, na.rm = TRUE),
    r28 = mean(retained28d, na.rm = TRUE),
    r90 = mean(retained90d, na.rm = TRUE),
    pct_paid_any = mean(account_type %in% c("Plus", "PaidTeacher"), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(r7))

# Attach top drivers (by abs lift vs global) as a readable string per cluster
drivers_table <- cluster_summary_tidy %>%
  filter(!is.na(pct_lift_vs_global)) %>%
  group_by(cluster) %>%
  slice_max(abs(pct_lift_vs_global), n = 6) %>%
  summarise(
    top_drivers = paste0(
      str_replace(metric, "^mean_", ""),
      " (", sprintf("%+.0f%%", pct_lift_vs_global), ")",
      collapse = "; "
    ),
    .groups = "drop"
  )

summary_table <- cluster_fact_sheet %>%
  left_join(drivers_table, by = "cluster")

summary_table

# Write workbook of all tables
# cluster_fact_sheet, drivers_table, summary_table, cluster_summary_tidy, retention_by_cluster, lift_by_cluster, global_means, cluster_means, cluster_profile, axis_profile,
# transition_matrix, stability, pct_visited, state_ret, transition_ret, cluster_drivers, comparison, best_pair_contrasts,
# contrast_4_6, contrast_4_3, contrast_2_5, contrast_2_6, wide_high, high_archetype_1, high_archetype_2, breadth_by_cluster,
# entropy_compare

# Helper to make unique Excel-safe sheet names
make_unique_sheet_names <- function(names_in) {
  out <- character(length(names_in))
  used <- character(0)
  
  for (i in seq_along(names_in)) {
    base <- str_trim(names_in[i])
    base <- str_replace_all(base, "[:\\\\/?*\\[\\]]", " ")
    base <- str_replace_all(base, "\\s+", " ")
    base <- substr(base, 1, 31)
    
    candidate <- base
    j <- 1
    
    while (tolower(candidate) %in% tolower(used) || candidate == "") {
      suffix <- paste0(" ", j)
      candidate <- substr(base, 1, 31 - nchar(suffix))
      candidate <- paste0(candidate, suffix)
      j <- j + 1
    }
    
    out[i] <- candidate
    used <- c(used, candidate)
  }
  
  out
}

wb <- createWorkbook()

table_names <- c(
  "cluster_fact_sheet", "drivers_table", "summary_table", "cluster_summary_tidy",
  "retention_by_cluster", "lift_by_cluster", "global_means", "cluster_means",
  "cluster_profile", "axis_profile", "transition_matrix", "stability", "pct_visited",
  "state_ret", "transition_ret", "cluster_drivers", "comparison", "best_pair_contrasts",
  "contrast_4_6", "contrast_4_3", "contrast_2_5", "contrast_2_6", "wide_high",
  "high_archetype_1", "high_archetype_2", "breadth_by_cluster", "entropy_compare"
)

models <- list(
  "GLM Cluster Only" = model_cluster_only,
  "GLM Cluster + Prev + Changed" = cluster_move_glm1,
  "GLM Cluster + Prev + Direction" = cluster_direction_glm
)

model_tables <- purrr::imap(models, ~ glm_export_bundle_fast(.x))

tbl_sheet_names <- table_names %>%
  str_replace_all("_", " ") %>%
  str_to_title()

tbl_sheet_names <- make_unique_sheet_names(tbl_sheet_names)

for (i in seq_along(table_names)) {
  tbl_name <- table_names[i]
  sheet_name <- tbl_sheet_names[i]
  
  if (exists(tbl_name)) {
    tbl <- get(tbl_name)
    
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet = sheet_name, x = tbl)
    freezePane(wb, sheet = sheet_name, firstRow = TRUE)
  } else {
    message(paste("Skipping:", tbl_name, "- object not found"))
  }
}

model_sheets <- list()

for (model_name in names(model_tables)) {
  bundle <- model_tables[[model_name]]
  
  for (part in names(bundle)) {
    df_out <- bundle[[part]]
    
    proposed <- paste(model_name, "-", str_to_title(part)) %>%
      str_replace_all("_", " ") %>%
      str_replace_all("\\s+", " ") %>%
      str_trim()
    
    model_sheets[[proposed]] <- df_out
  }
}

model_sheet_names_unique <- make_unique_sheet_names(names(model_sheets))

for (i in seq_along(model_sheets)) {
  sheet_name <- model_sheet_names_unique[i]
  df_out <- model_sheets[[i]]
  
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet = sheet_name, x = df_out)
  setColWidths(wb, sheet = sheet_name, cols = 1:max(1, ncol(as.data.frame(df_out))), widths = "auto")
  freezePane(wb, sheet = sheet_name, firstRow = TRUE)
}

# save workbook
saveWorkbook(wb, "/Users/karstenwalker/Documents/Cluster_Analysis_3_5.xlsx", overwrite = TRUE)

### Next Steps
# Remove pure volume metrics temporarily
# Exclude features like sessions, questions_answered, minutes_active, and lifetime_sets_created
# look at contrasts again.
# If clusters still have meaningful sep then we might have style-based segments.
# If differences disappear then segs are primary volume/intensity based

###### Clustering using ECDF ######
# Removes magnitude, not as useful for prediction, but useful for archetypes

to_ecdf <- function(x) {
  if (all(x %in% c(0, 1, NA))) return(x)           # keep binary flags
  r <- rank(x, na.last = "keep", ties.method = "average")
  (r - 1) / (sum(!is.na(x)) - 1)                   # 0..1
}

X_ecdf <- df_entity %>%
  select(-user_id) %>%
  select(where(is.numeric)) %>%
  mutate(across(everything(), ~ replace_na(.x, 0))) %>%
  mutate(across(everything(), to_ecdf)) %>%
  as.matrix()

X_ecdf_scaled <- scale(X_ecdf)
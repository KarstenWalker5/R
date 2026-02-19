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

# Store your Google Cloud project ID
bq_auth()

project_id <- "compact-sylph-785"

# Define your SQL query (example using a public dataset)
sql <- "SELECT * 
        FROM `compact-sylph-785.Karsten.EL_1pct_1yr_LI_wk`"

# Run the query
tb <- bq_project_query(project_id, sql)

# Load helpers

source("/Users/karstenwalker/Documents/GitHub/R/Helpers/themes.R")

source("/Users/karstenwalker/Documents/GitHub/R/Helpers/df_transformations.R")

# Set Seed
set.seed(7)

# Download the results into an R data frame
el_data <- bq_table_download(tb)%>%
  filter(account_type!="Admin")%>%
  group_by(account_type) %>%
  mutate(z_score = scale(flashcards_questions_answered)) %>%
  filter(abs(z_score)<3,flashcards_questions_answered>0)%>%
  ungroup()

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
  filter(country_code=="US")%>%
  group_by(account_type, retained7d) %>%
  summarize(users_7d=n_distinct(user_id))%>%
  rename(retained=retained7d)%>%
  left_join(el_data%>%
              filter(country_code=="US")%>%
              group_by(account_type, retained28d) %>%
              summarize(users_28d=n_distinct(user_id))%>%
                          rename(retained=retained28d), by=c("account_type", "retained"))%>%
  left_join(el_data%>%
              filter(country_code=="US")%>%
              group_by(account_type, retained90d) %>%
              summarize(users_90d=n_distinct(user_id))%>%
              rename(retained=retained90d), by=c("account_type", "retained"))

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
  filter(flashcards_questions_answered_log>0)%>%
  select(-user_id, -uid, -country_code, -contains("date"))%>%
  correlate(use = "complete.obs", method = "pearson")

cor_matrix_spear <- el_data %>%
  filter(flashcards_questions_answered>0)%>%
  select(-user_id, -uid, -country_code, -contains("date"))%>%
  correlate(use = "complete.obs", method = "spearman")

cor_log<-el_logged%>%
  select(-user_id, -uid, -country_code, -contains("date"))%>%
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
el_sub<-el_logged%>%
  group_by(user_id, week)%>%
  filter(row_number()==1, reported_user_type=="Student")%>%
  ungroup()

# Drop columns not used for clustering
df_entity <- el_sub%>%
  ungroup() %>%
  select( -uid, -country_code, -platform, -account_type,
          -contains("retained"), -contains("week"), -year,
          -age, -reported_user_type
  )

# Create an ID lookup for post-clustering
ids <- el_logged%>%
  group_by(user_id, week)%>%
  filter(row_number()==1, reported_user_type=="Student")%>%
  ungroup()%>%
  select(user_id, week)

# Replace NA
X <- df_entity %>%
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
set.seed(1)

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

# After 6 marginal gains stabilizes at <3%, noise coming to the far right, 6-8 is prob the sweet spot

### Clustering on PCs
k_final <- 6

mbk <- ClusterR::MiniBatchKmeans(
  data = X_pcs_full,
  clusters = k_final,
  batch_size = 20000,   # 10k–50k typical; increase for stability
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

# Build mapping table
df_clustered_uw <- ids %>%
  mutate(cluster = clusters) %>%
  distinct(user_id, week, .keep_all = TRUE)   

# Join back to full dataset
df_full_with_clusters <- el_sub%>%
  ungroup() %>%
  inner_join(df_clustered_uw, by = c("user_id", "week"))

# Join clusters back to a frame that includes retention
retention_by_cluster <- df_full_with_clusters %>%
  group_by(cluster) %>%
  summarise(
    n = n(),
    retention_rate7d = mean(retained7d, na.rm = TRUE),
    retention_rate28d = mean(retained28d, na.rm = TRUE),
    retention_rate90d = mean(retained90d, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(cluster))

df_full_with_clusters %>%
  summarise(
    n = n(),
    retention_rate7d = mean(retained7d, na.rm = TRUE),
    retention_rate28d = mean(retained28d, na.rm = TRUE),
    retention_rate90d = mean(retained90d, na.rm = TRUE),
    .groups = "drop"
  )

retention_by_cluster

# Lift by cluster
overall <- df_full_with_clusters %>%
  summarise(
    r7  = mean(retained7d,  na.rm = TRUE),
    r28 = mean(retained28d, na.rm = TRUE),
    r90 = mean(retained90d, na.rm = TRUE)
  )

lift_by_cluster<-df_full_with_clusters %>%
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

### Try without PCA 
mbk2 <- ClusterR::MiniBatchKmeans(
  data = X_scaled,
  clusters = k_final,
  batch_size = 20000,
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

### Computer cluster stats

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

# Big flows that look like a ladder:
# 1 > 6 is huge (0.382)
# 6 > 1 is also sizable (0.205)
# 2 <> 5 <> 6 seems meaningful
# Might imply that this is state/engagement based, not a unique entity.

# Stability score per cluster
# Stay prob > 0.7 → cluster behaves like a type
# Stay prob 0.4–0.7 → semi-stable
# Stay prob < 0.4 → cluster is a state

stability <- transitions %>%
  filter(cluster == next_cluster) %>%
  select(cluster, transition_prob) %>%
  rename(stay_prob = transition_prob)

stability

# How many unique clusters each user occupies across all weeks
user_cluster_span <- df_full_with_clusters %>%
  group_by(user_id) %>%
  summarise(
    n_clusters_visited = n_distinct(cluster),
    .groups = "drop"
  )

summary(user_cluster_span$n_clusters_visited)

# table
table(user_cluster_span$n_clusters_visited)

# 58.6% in 1 cluster, 26.7% in 2
pct_visited<-user_cluster_span %>%
  count(n_clusters_visited) %>%
  mutate(
    pct_users = 100 * n / sum(n)
  ) %>%
  arrange(n_clusters_visited)

pct_visited

# Movement rate
# < 20% → mostly stable types
# 20–40% → moderate state behavior
# 50% → highly dynamic states
# 46.3%

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

# Probability a user leaves a cluster
# 43%
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
# Cluster 6 is danger zone, 4 is peak retained
# Moving from 6 to 1 nearly doubles retention
# 6 to 4 puts you into top tier retention.
# 1 to 6 cuts retention dramatically.
# 4 or 5 to 6 is catastrophic.
# Users who move 6 to 1 retain at ~87%, stay in 6 retain at ~49%, fall 4 → 6 retain at ~56%.

state_vs_ret<-df_states %>%
  filter(!is.na(transition)) %>%
  group_by(transition) %>%
  summarise(
    n = n(),
    r7 = mean(retained7d, na.rm = TRUE)
  ) %>%
  arrange(desc(n))

# Simple regression of current cluster vs prev
cluster_move_glm1<-glm(retained7d ~ factor(cluster) + factor(prev_cluster) + changed_cluster,
    data = df_states,
    family = binomial)

# 15% of deviance explained by just state + previous state + change indicator.
# baseline is cluster 1 since ommitted
# Even after controlling for current cluster, previous cluster still matters.
# Coefs
#   2	-0.395	lower retention
#   3	+0.270	higher
#   4	+0.510	much higher
#   5	-0.226	slightly lower
#   6	-1.958	dramatically lower

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

# plogis(2.0492) ≈ 0.886
# Users who moved down into cluster 1 from cluster 1’s baseline scenario retain ~88.6%
# odds ratio: exp(-2.1451) = 0.117
# So being in cluster 6 reduces odds of retention by ~88% relative to cluster 1.
# Cluster 4: +0.463
# odds ratio:exp(0.463)= 1.59, cluster 4 increases odds by ~59%.
# prev_cluster 3: +0.418 even after controlling for current cluster and direction.
# Staying in the same cluster increases odds by 23%, moving up increases odds by 37%.

# Retention by current state and by transition
state_ret <- df_states %>%
  group_by(cluster) %>%
  summarise(n = n(),
            r7 = mean(retained7d, na.rm = TRUE),
            .groups = "drop")

transition_ret <- df_states %>%
  filter(!is.na(transition)) %>%
  group_by(transition) %>%
  summarise(n = n(),
            r7 = mean(retained7d, na.rm = TRUE),
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

# Summary
# * There is a dominant behavioral pattern per user.
# * Some users move up/down activation levels.
# * Movement is limited (mostly 1–3 clusters).
# * Clusters likely represent engagement tiers or phases.
# * Should not use for modeling since they represent states, could use prior cluster as lagging indicator
# Move rate = 0.463 → ~46% of week-to-week transitions change cluster.
# Ever changed = 0.414 → ~41% of users change cluster at least once.
# Clusters visited: 58.6% stay in 1 cluster; 26.7% visit 2; 11% visit 3.
# Stay probabilities
#   * c6: 0.625 (stickiest)
#   * c2: 0.574
#   * c4: 0.539
#   * c1: 0.467
#   * c3: 0.340
#   * c5: 0.342
# This is not user personality types, which would show much higher diagonals and much lower move rate.
# Looks like engagement states where people move between levels with some stickiness for the big states (especially cluster 6).

###### Initial insights ######
# primarily separating:
# * Low activation free users, casual/low-commitment users who never fully activate (cluster 6)
# * Moderate free users (cluster 2, 7)
# * High activation users, fully activated habitual learners (cluster 4, 3)
# * Paid-heavy engaged users (cluster 8, 5)
# The retention structure aligns mostly with:
# * Engagement intensity
# * Some monetization influence
# * Not strongly with role
#
# Key Takeaways
# * Clusters are not personality types. This is a state system, not identity.
#   * 46% week-to-week movement
#   * 41% of users ever change cluster
#   * Stay probabilities mostly 0.34–0.62
# *  Clusters represent engagement regimes
#   * Likely something like Cluster 6 > Low activation, Cluster 2/5 > Moderate activation, 
#     Cluster 4 > High creator activation, Cluster 3 > High consumer activation
# * Retention is state-dependent.
# * Your cluster-only logistic regression  explains ~15% deviance, large for behavioral data.

### Cluster drivers
# Compare
# * 4 vs 6 (retention high/low respectively)
# * 4 vs 3 (both high retention)
# * 3 vs 6 (similar retention, different paid)
# * 2 vs 4 (middle clusters)
# * 1 vs 2 (lower-mid retention)

cluster_drivers <- cluster_summary_tidy %>%
  filter(abs(pct_lift_vs_global) > 20) %>%   # only meaningful differences
  group_by(cluster) %>%
  slice_max(abs(pct_lift_vs_global), n = 10) %>%
  arrange(cluster, desc(abs(pct_lift_vs_global)))

# Why is cluster 4 so much higher?
comparison <- cluster_summary_tidy %>%
  filter(cluster %in% c(4, 6)) %>%
  select(cluster, metric, pct_lift_vs_global) %>%
  pivot_wider(names_from = cluster, values_from = pct_lift_vs_global) %>%
  mutate(
    diff_4_vs_6 = `4` - `6`
  )

comparison %>%
  arrange(desc(abs(diff_4_vs_6))) %>%
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
contrast_4_6<- contrast_clusters(cluster_summary_tidy, 4, 6)

contrast_4_3<-contrast_clusters(cluster_summary_tidy, 4, 3)

contrast_2_5<-contrast_clusters(cluster_summary_tidy, 2, 5)

contrast_2_6<-contrast_clusters(cluster_summary_tidy, 2, 6)

# Compare highest retention clusters
high <- c(4, 3, 1)

wide_high <- cluster_summary_tidy %>%
  filter(cluster %in% high) %>%
  select(cluster, metric, pct_lift_vs_global) %>%
  pivot_wider(names_from = cluster, values_from = pct_lift_vs_global)

# Top 3 retention cluster diffs
high_archetype_1 <- wide_high %>%
  mutate(d43 = `4` - `3`, d31 = `3` - `1`) %>%
  filter(d43 > 0, d31 > 0) %>%
  arrange(desc(d43 + d31)) %>%
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

df_div <- df_full_with_clusters %>%
  mutate(
    feature_breadth = rowSums(across(all_of(feature_cols), ~ as.integer(replace_na(.x, 0) > 0)))
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
saveWorkbook(wb, "/Users/karstenwalker/Documents/Cluster_Analysis_V1.xlsx", overwrite = TRUE)

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
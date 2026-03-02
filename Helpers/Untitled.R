

# 1 row per user/wk, probably redundant but a good check
el_sub<-el_logged%>%
  group_by(user_id, week)%>%
  filter(row_number()==1)%>%
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
  filter(row_number()==1)%>%
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

### Clustering on PCs
k_final <- 8  # choose via elbow or stability tests

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

retention_by_cluster

# Lift by cluster
overall <- df_full_with_clusters %>%
  summarise(
    r7  = mean(retained7d,  na.rm = TRUE),
    r28 = mean(retained28d, na.rm = TRUE),
    r90 = mean(retained90d, na.rm = TRUE)
  )

df_full_with_clusters %>%
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
      .names = "mean_{.col}"
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
      .names = "mean_{.col}"
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

### Cluster drivers
# Compare
# * 4 vs 8 (retention high/low respectively)
# * 4 vs 3 (both high retention)
# * 3 vs 8 (similar retention, different paid)
# * 2 vs 7 (middle clusters)
# * 1 vs 2 (lower-mid retention)

cluster_drivers <- cluster_summary_tidy %>%
  filter(abs(pct_lift_vs_global) > 20) %>%   # only meaningful differences
  group_by(cluster) %>%
  slice_max(abs(pct_lift_vs_global), n = 8) %>%
  arrange(cluster, desc(abs(pct_lift_vs_global)))

# Why is cluster 4 so much higher?
comparison <- cluster_summary_tidy %>%
  filter(cluster %in% c(4, 8)) %>%
  select(cluster, metric, pct_lift_vs_global) %>%
  pivot_wider(names_from = cluster, values_from = pct_lift_vs_global) %>%
  mutate(
    diff_4_vs_8 = `4` - `8`
  )

comparison %>%
  arrange(desc(abs(diff_4_vs_8))) %>%
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
contrast_4_8<- contrast_clusters(cluster_summary_tidy, 4, 8)

contrast_4_3<- contrast_clusters(cluster_summary_tidy, 4, 3)

contrast_3_8<-contrast_clusters(cluster_summary_tidy, 3, 8)

contrast_2_7<-contrast_clusters(cluster_summary_tidy, 2, 7)

# Remove pure volume metrics temporarily
# Exclude features like sessions, questions_answered, minutes_active, and lifetime_sets_created
# look at contrasts again.
# If clusters still have meaningful sep then we might have style-based segments.
# If differences disappear then segs are primary volume/intensity based

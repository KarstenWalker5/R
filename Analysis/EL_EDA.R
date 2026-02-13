library(dplyr)
library(ggplot2)
library(mgcv)
library(readr)
library(scales)
library(bigrquery)
library(tidyverse)

# Store your Google Cloud project ID
bq_auth()

project_id <- "compact-sylph-785"

# Define your SQL query (example using a public dataset)
sql <- "SELECT * 
        FROM `compact-sylph-785.Karsten.EL_1pct_1yr_LI`"

# Run the query
tb <- bq_project_query(project_id, sql)

# Download the results into an R data frame
el_data <- bq_table_download(tb)%>%
  filter(account_type!="Admin")%>%
  group_by(account_type) %>%
  mutate(z_score = scale(flashcards_questions_answered)) %>%
  filter(abs(z_score)<3)

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
  
### Correlations
cor_matrix <- el_data %>%
  filter(flashcards_questions_answered>0)%>%
  select(-user_id, -date, -uid, -country_code, -z_score)%>%
  correlate(use = "complete.obs", method = "pearson")

el_transformed<-el_data%>%
  ungroup()%>%
  filter(flashcards_questions_answered>0)%>%
  select(-user_id, -date, -uid, -country_code, -z_score,-platform, -account_type, -contains('retained'))%>%
  mutate(across(is.numeric, log1p))%>%
  bind_cols(el_data%>%
              ungroup()%>%
              filter(flashcards_questions_answered>0)%>%
              select(contains('retained')))

cor_transformed<-el_transformed%>%
  correlate(use = "complete.obs", method = "pearson")

# convert to long format and filter strong correlations
cor_long <- cor_matrix %>%
  stretch() %>%
  filter(abs(r) > 0.5, x != y)

cor_long

cor_retention<-cor_matrix %>%
  stretch()%>%
  filter(str_detect(x, "retained"))

cor_retention_trans<-cor_matrix %>%
  stretch()%>%
  filter(str_detect(x, "retained"))

### Clustering
# install.packages(c("dbscan", "dplyr", "ggplot2"))
library(dbscan)

# 1) Prepare clustering frame
df_cluster <- el_data %>%
  ungroup() %>%
  filter(flashcards_questions_answered > 0, retained7d == 1) %>%
  select(-date, -uid, -country_code, -z_score, -platform, 
         -account_type,-contains("retained"),-age, 
         -days_until_next_session, -reported_user_type)

ids <- df_cluster %>% select(user_id)

X <- df_cluster %>%
  select(-user_id) %>%
  mutate(across(where(is.integer), as.numeric)) %>%
  select(where(is.numeric)) %>%
  mutate(across(everything(), ~ replace_na(.x, 0)))  # counts + flags

# Drop zero-variance columns (prevents NAs from scale())
X <- X[, sapply(X, sd) > 0, drop = FALSE]

# 2) Scale features
X_scaled <- scale(as.matrix(X))

stopifnot(!anyNA(X_scaled))

# 3) (Optional) Choose eps via kNN distance plot
kNNdistplot(X_scaled, k = 5)

abline(h = 2, col = "red")  # move this line to the elbow you see

# 4) Fit DBSCAN (adjust eps based on the elbow)
fit <- dbscan(X_scaled, eps = 2, minPts = 5)

# 5) Create joinable cluster labels + join back
cluster_labels <- ids %>%
  mutate(cluster_id = fit$cluster)  # 0 = noise

el_with_clusters <- el_data %>%
  left_join(cluster_labels, by = "user_id")

# 6) Plot clusters with ggplot (PCA to 2D for visualization)
pca <- prcomp(X_scaled, center = FALSE, scale. = FALSE)

plot_df <- cluster_labels %>%
  bind_cols(as.data.frame(pca$x[, 1:2])) %>%
  mutate(
    cluster_id = factor(cluster_id),
    is_noise = cluster_id == "0"
  )

ggplot(plot_df, aes(x = PC1, y = PC2, color = cluster_id)) +
  geom_point(aes(alpha = is_noise), size = 2) +
  scale_alpha_manual(values = c(`TRUE` = 0.35, `FALSE` = 0.9), guide = "none") +
  theme_minimal() +
  labs(
    title = "DBSCAN clusters (visualized with PCA)",
    subtitle = "Cluster 0 = noise",
    color = "Cluster"
  )
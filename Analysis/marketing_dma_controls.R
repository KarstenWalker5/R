
# Install bigrquery if you haven't already: install.packages("bigrquery")
library(bigrquery)
library(dplyr)
library(lubridate)
library(MarketMatching) 

# Replace with your actual GCP Project ID
project_id <- "compact-sylph-785" 

# Assuming your BigQuery script was saved as a table called `sample_dma_summary`
sql <- "SELECT * FROM `celdredge.sample_dma_summary`" 

bq_data <- bq_project_query(project_id, sql)

df <- bq_table_download(bq_data)

str(df)

df_sub<-df%>%
  filter(date>="2025-08-01")

city_summary<-df_sub%>%
  group_by(city)%>%
  mutate_if(is.integer, as.numeric)%>%
  summarize(users=sum(active_users),
            sub_ratio=(sum(upgrades)/sum(active_users))*100,
            total_subs=sum(upgrades),
            min_subs=min(upgrades),
            max_subs=max(upgrades),
            sd_subs=sd(upgrades),
            google_spend=sum(google_spend),
            meta_spend=sum(meta_spend))%>%
  arrange(desc(users))%>%
  mutate(rank=row_number())

pre_spend<-df%>%
  filter(date<="2025-08-01")%>%
  group_by(city)%>%
  mutate_if(is.integer, as.numeric)%>%
  summarize(users=sum(active_users),
            sub_ratio=(sum(upgrades)/sum(active_users))*100,
            total_subs=sum(upgrades),
            min_subs=min(upgrades),
            max_subs=max(upgrades),
            sd_subs=sd(upgrades),
            google_spend=sum(google_spend),
            meta_spend=sum(meta_spend))%>%
  arrange(desc(users))%>%
  mutate(rank=row_number())

treatment_cities<-as.data.frame(c('charlotte',
                                  'denver',
                                  'detroit',
                                  'indianapolis',
                                  'orlando-daytona beach-melbourne',
                                  'pittsburgh',
                                  'sacramento-stockton-modesto',
                                  'seattle-tacoma'))%>%
  rename(treat=1)

test<-df%>%
  filter(city %in% treatment_cities$treat)

matching_df<-df%>%
  mutate(treat=ifelse(city %in% treatment_cities$treat,1,0))

synthetic_panel <- best_matches_all_cities(
  data           = matching_df,
  matching_metric = "upgrades",
  date_col        = "date",
  city_col        = "city",
  treat_col       = "treat",
  match_window    = c("2024-01-01", "2025-08-01"),   # your pre-period
  matches_per_city     = 10,
  parallel=TRUE,
  match_treated_only = TRUE
)

dplyr::glimpse(synthetic_panel)

synthetic_panel <- build_synthetic_control(
  data           = matching_df,
  matching_metric = "upgrades",
  date_col        = "date",
  city_col        = "city",
  treat_col       = "treat",
  match_window    = c("2024-01-01", "2024-05-31"),
  n_matches       = 10,
  post_col        = NULL
)

#Print all treated cities and their donors
synthetic_panel %>%
  dplyr::distinct(treated_city, synthetic_cities) %>%
  dplyr::mutate(donors = purrr::map_chr(synthetic_cities, ~ paste(.x, collapse = ", "))) %>%
  dplyr::select(treated_city, donors)

control_list<-synthetic_panel %>%
  dplyr::distinct(treated_city, synthetic_cities) %>%
  tidyr::unnest_longer(synthetic_cities)
# Geo test power table rate-base variance reduction
# Assumes df2 has columns: dma, date, upgrades, active_users
# Produces ONE table with 9 permutations:
#   top10_included ∈ {0,5,10} × designs ∈ {25/25,50/50,68/68}
# Eligibility rule:
#   - Rank DMAs by total upgrades (descending).
#   - If top10_included = 0: force all top10 to BAU.
#   - If top10_included = 5: force top 5 (largest of top10) to BAU; include ranks 6–10.
#   - If top10_included = 10: include all top10.
#   - After exclusion, take NEXT K DMAs by total upgrades, where K = n_test + n_ctrl.
#     (These K are eligible for random assignment; everyone else BAU.)
# Power model:
#   - DMA-level DiD metric is now upgrade RATE:
#       rate = sum(upgrades) / sum(active_users)
#     diff = post_rate - pre_rate
#   - SD estimated as median SD across windows
#   - Two-sample t-test power with unequal n supported (pwr.t2n.test)
# Outputs:
#   - detectable_drop_rate_per_day (MDE in rate units)
#   - detectable_pct_drop_rate = MDE / avg_daily_rate_in_eligible_pool

library(dplyr)
library(tidyr)
library(purrr)
library(pwr)
library(bigrquery)
library(dplyr)
library(lubridate)

###### Load data ######

# Replace with your actual GCP Project ID
project_id <- "compact-sylph-785" 

# Assuming your BigQuery script was saved as a table called `sample_dma_summary`
sql <- "SELECT * FROM `celdredge.sample_dma_summary`" 

bq_data <- bq_project_query(project_id, sql)

df <- bq_table_download(bq_data)

set.seed(123)

###### Power Analysis ######

# Settings
durations_weeks <- c(2, 3, 4, 6, 8)
alpha <- 0.05
target_power <- 0.80

designs <- tibble::tribble(
  ~design, ~n_test, ~n_ctrl,
  "25/25", 25L,     25L,
  "50/50", 50L,     50L,
  "68/68", 68L,     68L
)

top10_included_vals <- c(0L, 5L, 10L)

# 1) Rank DMAs by total upgrades once (unchanged)

ranked <- df2 %>%
  group_by(dma) %>%
  summarise(total_upgrades = sum(upgrades, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(total_upgrades)) %>%
  mutate(rank = row_number())

top10 <- ranked %>% filter(rank <= 10) %>% pull(dma)

# Deterministic: exclude the largest 5 of the top 10 when top10_included=5
top10_exclude_5 <- top10[1:5]

# 2) Build eligible pool (unchanged)

get_eligible_nextK <- function(df, ranked_tbl, top10, top10_exclude_5,
                               top10_included, eligible_k_after_exclusion) {
  
  excluded <- if (top10_included == 10L) {
    character(0)
  } else if (top10_included == 0L) {
    top10
  } else if (top10_included == 5L) {
    top10_exclude_5
  } else {
    stop("top10_included must be 0, 5, or 10")
  }
  
  ranked_post <- ranked_tbl %>% filter(!dma %in% excluded)
  
  keep_dmas <- ranked_post %>%
    slice(1:eligible_k_after_exclusion) %>%
    pull(dma)
  
  df %>% filter(dma %in% keep_dmas)
}

# 3) Estimate SD of DMA-level DiD metric (CHANGED: use RATE)
#    rate = sum(upgrades)/sum(active_users)

estimate_sd_did <- function(df, T_weeks) {
  T_days <- T_weeks * 7
  df <- df %>% arrange(dma, date)
  
  did <- df %>%
    group_by(dma) %>%
    mutate(idx = row_number(),
           window = floor((idx - 1) / (2 * T_days))) %>%
    group_by(dma, window) %>%
    filter(n() == 2 * T_days) %>%
    mutate(period = ifelse(row_number() <= T_days, "pre", "post")) %>%
    group_by(dma, window, period) %>%
    summarise(
      sum_up = sum(upgrades, na.rm = TRUE),
      sum_au = sum(active_users, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      # --- CHANGE: rate outcome for variance reduction ---
      rate = ifelse(sum_au > 0, sum_up / sum_au, NA_real_)
    ) %>%
    select(dma, window, period, rate) %>%
    pivot_wider(names_from = period, values_from = rate) %>%
    filter(is.finite(pre), is.finite(post)) %>%
    mutate(diff = post - pre)
  
  sd_by_window <- did %>%
    group_by(window) %>%
    summarise(sd_did = sd(diff, na.rm = TRUE), .groups = "drop") %>%
    filter(is.finite(sd_did))
  
  median(sd_by_window$sd_did, na.rm = TRUE)
}

# 4) Power & MDE helpers (unchanged)

power_from_sd <- function(sd_did, mde, n_test, n_ctrl, alpha = 0.05) {
  d <- mde / sd_did
  pwr.t2n.test(n1 = n_test, n2 = n_ctrl, d = d, sig.level = alpha)$power
}

mde_from_sd <- function(sd_did, n_test, n_ctrl, target_power = 0.8, alpha = 0.05) {
  f <- function(mde) power_from_sd(sd_did, mde, n_test, n_ctrl, alpha) - target_power
  uniroot(f, lower = 1e-6, upper = 50 * sd_did)$root
}

# 5) Run one power table (CHANGED: baseline uses avg DAILY RATE)

run_power_table <- function(df_eligible, durations_weeks, n_test, n_ctrl,
                            target_power = 0.8, alpha = 0.05) {
  
  n_eligible <- n_distinct(df_eligible$dma)
  
  # --- CHANGE: baseline average daily upgrade RATE in eligible pool ---
  avg_daily_rate <- df_eligible %>%
    group_by(dma) %>%
    summarise(
      sum_up = sum(upgrades, na.rm = TRUE),
      sum_au = sum(active_users, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(rate = ifelse(sum_au > 0, sum_up / sum_au, NA_real_)) %>%
    summarise(avg_rate = mean(rate, na.rm = TRUE)) %>%
    pull(avg_rate)
  
  map_dfr(durations_weeks, function(T) {
    sd_did <- estimate_sd_did(df_eligible, T)
    mde    <- mde_from_sd(sd_did, n_test, n_ctrl, target_power, alpha)
    
    tibble(
      weeks = T,
      sd_did_rate = sd_did,
      detectable_drop_rate_per_day = mde,
      detectable_pct_drop_rate = mde / avg_daily_rate,
      avg_daily_rate = avg_daily_rate,
      eligible_dmas = n_eligible
    )
  })
}

# 6) Build the ONE combined table (9 permutations)

results_all_rate <- tidyr::crossing(
  top10_included = top10_included_vals,
  designs
) %>%
  mutate(K = n_test + n_ctrl) %>%
  mutate(
    eligible_df = pmap(
      list(top10_included, K),
      \(ti, k) get_eligible_nextK(df2, ranked, top10, top10_exclude_5, ti, k)
    )
  ) %>%
  mutate(
    power_tbl = pmap(
      list(eligible_df, n_test, n_ctrl),
      \(ed, nt, nc) run_power_table(ed, durations_weeks, nt, nc, target_power, alpha)
    )
  ) %>%
  select(top10_included, design, n_test, n_ctrl, K, power_tbl) %>%
  unnest(power_tbl) %>%
  mutate(
    n_bau = 278 - K
  ) %>%
  select(
    top10_included, design, K, n_test, n_ctrl, n_bau,
    weeks, sd_did_rate, avg_daily_rate,
    detectable_drop_rate_per_day, detectable_pct_drop_rate,
    eligible_dmas
  ) %>%
  arrange(top10_included, design, weeks)

results_all


# Sanity check which DMAs are in the eligible pool for a scenario
# Example: top10_included=0, design=50/50 => K=100 => should correspond to ranks 11–110
elig_example <- get_eligible_nextK(df2, 
                                   ranked, 
                                   top10, 
                                   top10_exclude_5, 
                                   top10_included = 0L, 
                                   eligible_k_after_exclusion = 100L)

ranked %>% 
  filter(dma %in% unique(elig_example$dma)) %>% 
  arrange(rank) %>% 
  slice(1:50)

###### Now identify which go into test/control ######
assign_test_control <- function(df, ranked_tbl, top10, top10_exclude_5,
                                top10_included, n_test, n_ctrl) {
  K <- n_test + n_ctrl
  
  elig_df <- get_eligible_nextK(df, ranked_tbl, top10, top10_exclude_5,
                                top10_included, eligible_k_after_exclusion = K)
  
  elig_dmas <- sort(unique(elig_df$dma))
  
  test_dmas <- sample(elig_dmas, n_test)
  remaining <- setdiff(elig_dmas, test_dmas)
  ctrl_dmas <- sample(remaining, n_ctrl)
  
  assignment <- tibble(dma = sort(unique(df$dma))) %>%
    mutate(group = case_when(
      dma %in% test_dmas ~ "TEST",
      dma %in% ctrl_dmas ~ "CONTROL",
      TRUE ~ "BAU"
    ))
  
  list(
    assignment = assignment,
    test_dmas = test_dmas,
    control_dmas = ctrl_dmas,
    eligible_dmas = elig_dmas
  )
}

# Example: top10_included=0 and 50/50 design
out <- assign_test_control(df2, ranked, top10, top10_exclude_5,
                           top10_included = 0L, n_test = 50L, n_ctrl = 50L)

# Final assignment table
assignment_table<-out$assignment%>%
  left_join(df2%>%
              select(dma, active_users,signups, google_spend, meta_spend)%>%
              group_by(dma)%>%
              summarise(across(everything(), list(sum))), by="dma")


# 7) Balance report
# Requires: df2 (dma,date,upgrades,active_users), assignment (dma,group)

#  pre-perod window
pre_start_date <- as.Date("2025-01-01")  # <-- change

pre_end_date   <- as.Date("2025-01-28")  # <-- change

#  Build DMA-level pre-period summary for treatment/control 
balance_df <- df2 %>%
  filter(date >= pre_start_date & date <= pre_end_date) %>%
  group_by(dma) %>%
  summarise(
    pre_upgrades = sum(upgrades, na.rm = TRUE),
    pre_active_users = sum(active_users, na.rm = TRUE),
    pre_days = n(),
    pre_rate = ifelse(pre_active_users > 0, pre_upgrades / pre_active_users, NA_real_),
    .groups = "drop"
  ) %>%
  left_join(assignment_table, by = "dma") %>%
  filter(group %in% c("TEST", "CONTROL")) %>%
  mutate(treatment = ifelse(group == "TEST", 1, 0)) %>%
  filter(is.finite(pre_rate))

#  Summary table (means/medians by group) 
balance_summary <- balance_df %>%
  group_by(group) %>%
  summarise(
    n_dmas = n(),
    mean_pre_rate = mean(pre_rate, na.rm = TRUE),
    median_pre_rate = median(pre_rate, na.rm = TRUE),
    mean_pre_upgrades = mean(pre_upgrades, na.rm = TRUE),
    median_pre_upgrades = median(pre_upgrades, na.rm = TRUE),
    mean_pre_active_users = mean(pre_active_users, na.rm = TRUE),
    median_pre_active_users = median(pre_active_users, na.rm = TRUE),
    .groups = "drop"
  )

balance_summary

# Treatment Vs. Control t-tests (pre-period balance checks) 
tt_pre_rate <- t.test(pre_rate ~ treatment, data = balance_df)

tt_pre_upgrades <- t.test(pre_upgrades ~ treatment, data = balance_df)

tt_pre_active_users <- t.test(pre_active_users ~ treatment, data = balance_df)

tt_pre_rate

tt_pre_upgrades

tt_pre_active_users

# Optional: standardized mean differences (SMD) for quick “how imbalanced?” 
# SMD < 0.10 is excellent balance; 0.10–0.20 is usually fine; >0.20 consider re-randomizing.
# Look at your primary metric

smd <- function(x, g) {
  m1 <- mean(x[g == 1], na.rm = TRUE); m0 <- mean(x[g == 0], na.rm = TRUE)
  s1 <- sd(x[g == 1], na.rm = TRUE);   s0 <- sd(x[g == 0], na.rm = TRUE)
  (m1 - m0) / sqrt((s1^2 + s0^2) / 2)
}

balance_smd <- tibble(
  metric = c("pre_rate", "pre_upgrades", "pre_active_users"),
  smd = c(
    smd(balance_df$pre_rate, balance_df$treatment),
    smd(balance_df$pre_upgrades, balance_df$treatment),
    smd(balance_df$pre_active_users, balance_df$treatment)
  )
)

balance_smd


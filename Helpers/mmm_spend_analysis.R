# This script analyzes our Marketing spend data and identifies issues that could complicate building a reliable MMM

# Load packages
  library(tidyverse)
  library(lubridate)
  library(scales)
  library(janitor)
  library(glue)


# Load data
file_path <- "/Users/karstenwalker/Downloads/mmm_spend.csv"

# Load themes
source(file="/Users/karstenwalker/Documents/GitHub/R/Helpers/themes.R")

# Set date range
expected_start <- as.Date("2022-02-03")
expected_end   <- as.Date("2026-02-03")

# Helper functions
trim_names_safe <- function(df) {
  # trims whitespace and standardizes to snake_case
  df %>% janitor::clean_names()
}

parse_currency_safe <- function(x) {
  # robust parsing of "$1,234.56" or "1,234.56" or "" to numeric
  # returns NA for truly missing strings
  readr::parse_number(as.character(x), na = c("", "NA", "N/A", "null", "NULL"))
}

iqr_outlier_flags <- function(x, k = 3) {
  # returns logical vector: TRUE if beyond k*IQR from Q1/Q3
  q1 <- quantile(x, 0.25, na.rm = TRUE)
  q3 <- quantile(x, 0.75, na.rm = TRUE)
  i  <- q3 - q1
  if (is.na(i) || i == 0) return(rep(FALSE, length(x)))
  (x < (q1 - k * i)) | (x > (q3 + k * i))
}

# Read data
raw <- readr::read_csv(
  file_path,
  col_types = readr::cols(.default = readr::col_character())
)

# Keep a copy of original names to detect trailing/odd whitespace issues
original_names <- names(raw)

# Detect column name issues (leading/trailing whitespace)
name_whitespace_issues <- tibble(
  original_name = original_names,
  has_leading_ws = str_detect(original_name, "^\\s+"),
  has_trailing_ws = str_detect(original_name, "\\s+$"),
  has_internal_double_space = str_detect(original_name, "\\s{2,}")) %>%
  filter(has_leading_ws | has_trailing_ws | has_internal_double_space)

# Clean names for downstream work
df0 <- raw %>% trim_names_safe()

if (!"day" %in% names(df0)) {
  stop("No 'day' column found after clean_names(). Please verify the date column name.")
}

# Detect & remove summary/non-date rows
# Mark rows where "day" cannot be parsed as a Date
df1 <- df0 %>%
  mutate(
    day_raw = day,
    day_parsed = suppressWarnings(lubridate::ymd(day_raw, quiet = TRUE)),
    is_non_date_row = is.na(day_parsed)
  )

non_date_rows <- df1 %>%
  filter(is_non_date_row) %>%
  select(day_raw, everything())

# Keep only valid daily rows
df <- df1 %>%
  filter(!is_non_date_row) %>%
  mutate(day = day_parsed) %>%
  select(-day_raw, -day_parsed, -is_non_date_row)

# Spend columns and parse currency/commas to numeric
spend_cols <- setdiff(names(df), "day")

# Parse numerics
df_num <- df %>%
  mutate(across(all_of(spend_cols), parse_currency_safe))

# Flag parsing failures: values that were non-empty strings but became NA
parsing_failures <- df %>%
  mutate(across(all_of(spend_cols), ~ as.character(.x))) %>%
  pivot_longer(cols = all_of(spend_cols), names_to = "channel", values_to = "raw_value") %>%
  filter(!is.na(raw_value), str_trim(raw_value) != "") %>%
  mutate(parsed_value = parse_currency_safe(raw_value)) %>%
  filter(is.na(parsed_value)) %>%
  group_by(channel) %>%
  summarise(
    n_unparseable = n(),
    examples = paste(head(unique(raw_value), 5), collapse = " | "),
    .groups = "drop"
  ) %>%
  arrange(desc(n_unparseable))

# Basic date range diagnostics

date_summary <- df_num %>%
  summarise(
    min_day = min(day, na.rm = TRUE),
    max_day = max(day, na.rm = TRUE),
    n_days = n_distinct(day),
    expected_start = expected_start,
    expected_end = expected_end,
    covers_expected_start = min_day <= expected_start,
    covers_expected_end = max_day >= expected_end
  )

# Missing days check (gaps in daily sequence)
all_days <- tibble(day = seq.Date(min(df_num$day), max(df_num$day), by = "day"))

missing_days <- anti_join(all_days, df_num %>% distinct(day), by = "day")

# Missing vs zero ambiguity & coverage diagnostics

coverage <- df_num %>%
  select(day, all_of(spend_cols)) %>%
  pivot_longer(cols = all_of(spend_cols), names_to = "channel", values_to = "spend") %>%
  group_by(channel) %>%
  summarise(
    n_days = n(),
    n_missing = sum(is.na(spend)),
    pct_missing = n_missing / n_days,
    n_zero = sum(!is.na(spend) & spend == 0),
    pct_zero = n_zero / n_days,
    n_positive = sum(!is.na(spend) & spend > 0),
    pct_positive = n_positive / n_days,
    first_non_missing = suppressWarnings(min(day[!is.na(spend)], na.rm = TRUE)),
    last_non_missing  = suppressWarnings(max(day[!is.na(spend)], na.rm = TRUE)),
    first_positive = suppressWarnings(min(day[!is.na(spend) & spend > 0], na.rm = TRUE)),
    last_positive  = suppressWarnings(max(day[!is.na(spend) & spend > 0], na.rm = TRUE)),
    .groups = "drop")%>%
  arrange(desc(pct_missing), desc(pct_zero))

# Identify channels that appear/disappear (missing at start or end)
coverage_flags <- coverage %>%
  mutate(
    missing_somewhere = n_missing > 0,
    appears_late = !is.infinite(as.numeric(first_non_missing)) & first_non_missing > min(df_num$day),
    ends_early   = !is.infinite(as.numeric(last_non_missing)) & last_non_missing < max(df_num$day),
    extremely_sparse = pct_zero >= 0.80 | n_positive < 60  # tweak threshold if desired
  )

# Outliers & regime shifts

# Identify top daily spend outliers per channel (IQR method + top-N)
spend_long <- df_num %>%
  select(day, all_of(spend_cols)) %>%
  pivot_longer(cols = all_of(spend_cols), names_to = "channel", values_to = "spend")

outlier_table <- spend_long %>%
  group_by(channel) %>%
  mutate(
    is_outlier = iqr_outlier_flags(spend, k = 3),
    rank_desc = dense_rank(desc(spend))) %>%
  filter(!is.na(spend)) %>%
  summarise(
    max_spend = max(spend, na.rm = TRUE),
    p99_spend = quantile(spend, 0.99, na.rm = TRUE),
    n_outliers_iqr3 = sum(is_outlier, na.rm = TRUE),
    worst_days = paste(day[rank_desc <= 3] %>% as.character(), collapse = ", "),
    .groups = "drop") %>%
  arrange(desc(max_spend))

# 28-day rolling mean; flag big jumps vs prior period
roll_28 <- spend_long %>%
  arrange(channel, day) %>%
  group_by(channel) %>%
  mutate(
    roll_mean_28 = slider::slide_dbl(spend, ~ mean(.x, na.rm = TRUE), .before = 27, .complete = TRUE),
    roll_mean_28_lag = lag(roll_mean_28, 28),
    roll_change = roll_mean_28 - roll_mean_28_lag,
    roll_change_pct = if_else(!is.na(roll_mean_28_lag) & roll_mean_28_lag != 0,
                              roll_change / roll_mean_28_lag,
                              NA_real_)) %>%
  ungroup()

regime_shift_flags <- roll_28 %>%
  filter(!is.na(roll_change_pct)) %>%
  group_by(channel) %>%
  summarise(
    biggest_jump_pct = max(roll_change_pct, na.rm = TRUE),
    biggest_drop_pct = min(roll_change_pct, na.rm = TRUE),
    day_biggest_jump = day[which.max(roll_change_pct)],
    day_biggest_drop = day[which.min(roll_change_pct)],
    .groups = "drop") %>%
  arrange(desc(abs(biggest_jump_pct)))

# Collinearity diagnostics (correlation matrix)

min_non_na <- 30

spend_wide <- df_num %>% select(all_of(spend_cols))

col_screen <- spend_wide %>%
  summarise(across(
    everything(),
    list(
      n_non_na = ~ sum(!is.na(.x)),
      sd_non_na = ~ sd(.x, na.rm = TRUE)
    ),
    .names = "{.col}__{.fn}"
  )) %>%
  pivot_longer(cols = everything(), names_to = "name", values_to = "value") %>%
  separate(name, into = c("channel", "metric"), sep = "__", extra = "merge", fill = "right") %>%
  pivot_wider(names_from = metric, values_from = value) %>%
  mutate(
    n_non_na = as.integer(n_non_na),
    sd_non_na = replace_na(sd_non_na, 0)
  )

keep_cols <- col_screen %>%
  filter(n_non_na >= min_non_na, sd_non_na > 0) %>%
  pull(channel)

spend_pairwise_ok <- spend_wide %>% select(all_of(keep_cols))

corr_pairwise <- suppressWarnings(cor(spend_pairwise_ok, use = "pairwise.complete.obs"))

if (is.null(rownames(corr_pairwise)) || is.null(colnames(corr_pairwise))) {
  rownames(corr_pairwise) <- colnames(spend_pairwise_ok)
  colnames(corr_pairwise) <- colnames(spend_pairwise_ok)
}

corr_long_pairwise <- as.data.frame(corr_pairwise, stringsAsFactors = FALSE) %>%
  rownames_to_column("channel1") %>%
  pivot_longer(-channel1, names_to = "channel2", values_to = "corr") %>%
  filter(channel1 < channel2) %>%
  filter(!is.na(corr), is.finite(corr))

ggplot(corr_long_pairwise, 
       aes(x = channel1, y = channel2, fill = corr)) +
  geom_tile() +
  scale_fill_gradient2(limits = c(-1, 1)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Channel Spend Correlation",
    x = "",
    y = "",
    fill = "Correlation")+
  theme_fancy()

# Colinearity, do channels still move together once we strip out calendar/trend effects?
spend_long <- df_num %>%
  select(day, all_of(spend_cols)) %>%
  pivot_longer(-day, names_to = "channel", values_to = "spend") %>%
  mutate(
    spend0 = replace_na(spend, 0),
    wday = wday(day, label = TRUE, week_start = 1),
    t = as.numeric(day - min(day))
  )

# Fit a separate “weekday + trend” model per channel and keep residuals
# For each channel we fit a linear model, wday captures weekly seasonality ns(t, df=6) is a smooth nonlinear time trend 
# (natural spline with 6 degrees of freedom), and outputs the residuals: the part of spend not explained by weekday + smooth trend.

resid_long <- spend_long %>%
  group_by(channel) %>%
  do({
    d <- .
    m <- lm(spend0 ~ wday + splines::ns(t, df = 6), data = d)
    tibble(day = d$day, channel = d$channel, resid = resid(m))
  }) %>%
  ungroup()

resid_wide <- resid_long %>%
  pivot_wider(names_from = channel, values_from = resid)

# compute correlation between channels using residuals
corr_resid <- suppressWarnings(cor(resid_wide 
                                   %>% select(-day), use = "pairwise.complete.obs"))

corr_long_resid <- as.data.frame(corr_resid) %>%
  rownames_to_column("channel1") %>%
  pivot_longer(-channel1, names_to = "channel2", values_to = "corr") %>%
  filter(channel1 < channel2) %>%
  filter(is.finite(corr))

corr_long_resid %>% 
  arrange(desc(abs(corr))) %>% 
  slice_head(n = 15)

corr_long_resid %>%
  rename(`Channel 1`=channel1,
         `Channel 2`=channel2,
         Correlation=corr)%>%
  kbl(
    caption = "Top Residual Channel Correlations After Removing Weekday and Trend Effects",
    align = "llrr",
    booktabs = TRUE
  ) %>%
  kable_styling(
    full_width = FALSE,
    bootstrap_options = c("striped", "hover", "condensed")
  ) 

# VIF
library(car)

wide0 <- df_num %>%
  select(day, all_of(spend_cols)) %>%
  mutate(across(all_of(spend_cols), ~ replace_na(.x, 0)))

wide0 <- wide0 %>% 
  mutate(y = rowSums(across(all_of(spend_cols)), na.rm = TRUE))

fit <- lm(y ~ . - day, data = wide0)

vif_tbl <- tibble(channel = names(car::vif(fit)), vif = as.numeric(car::vif(fit))) %>%
  arrange(desc(vif))

kappa_val <- kappa(model.matrix(fit))

vif_tbl

kappa_val

# Spend spikes
spike_tbl <- spend_long %>%
  group_by(channel) %>%
  arrange(day) %>%
  mutate(
    med_28 = rollapply(spend0, 28, median, fill = NA, align = "right"),
    mad_28 = rollapply(spend0, 28, mad, fill = NA, align = "right"),
    robust_z = (spend0 - med_28) / (mad_28 + 1e-6),
    is_spike = robust_z >= 8 ) %>%
  ungroup() %>%
  filter(is_spike) %>%
  arrange(desc(robust_z))

spike_tbl %>% 
  select(channel, day, spend0, med_28, mad_28, robust_z)%>%
  slice_head(n = 50)

# Weekday/weekend pattern diagnostics
weekday_summary <- df_num %>%
  mutate( wday = wday(day, label = TRUE, week_start = 1),
    is_weekend = wday(day) %in% c(6, 7)) %>%
  select(day, wday, is_weekend, all_of(spend_cols)) %>%
  pivot_longer(cols = all_of(spend_cols), names_to = "channel", values_to = "spend") %>%
  group_by(channel, wday, is_weekend) %>%
  summarise(
    avg_spend = mean(spend, na.rm = TRUE),
    med_spend = median(spend, na.rm = TRUE),
    pct_zero = mean(!is.na(spend) & spend == 0),
    .groups = "drop")

weekend_vs_weekday <- df_num %>%
  mutate(is_weekend = wday(day, week_start = 1) %in% c(6, 7)) %>%
  select(is_weekend, all_of(spend_cols)) %>%
  pivot_longer(cols = all_of(spend_cols), names_to = "channel", values_to = "spend") %>%
  group_by(channel, is_weekend) %>%
  summarise(avg_spend = mean(spend, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = is_weekend, values_from = avg_spend, names_prefix = "is_weekend_") %>%
  mutate(weekend_ratio = is_weekend_TRUE / is_weekend_FALSE) %>%
  arrange(weekend_ratio)

# Plots (ggplot2)
# spend over time by channel (facet)
spend_long %>%
  ggplot(aes(x = day, y = spend)) +
  geom_line(na.rm = TRUE) +
  facet_wrap(~ channel, scales = "free_y", ncol = 2) +
  scale_y_continuous(labels = dollar_format()) +
  labs(title = "Daily Spend by Channel",
    x = "Date",
    y = "Spend") +
  theme_fancy()

# missingness heatmap 
df_num %>%
  select(day, all_of(spend_cols)) %>%
  pivot_longer(cols = all_of(spend_cols), names_to = "channel", values_to = "spend") %>%
  mutate(present = !is.na(spend)) %>%
  ggplot(aes(x = day, y = fct_rev(channel), fill = present)) +
  geom_tile() +
  scale_fill_manual(values = c(`TRUE` = "grey20", `FALSE` = "grey85")) +
  labs(title = "Data Presence Heatmap (Missing vs Present)",
    x = "Date",
    y = "Channel",
    fill = "Present") +
  theme_fancy()

# zero-spend rate by channel
coverage %>%
  ggplot(aes(x = fct_reorder(channel, pct_zero), y = pct_zero)) +
  geom_col() +
  coord_flip() +
  scale_y_continuous(labels = percent_format()) +
  labs(title = "Zero-Spend Share by Channel",
    x = "Channel",
    y = "Share of days with spend = 0") +
  theme_awesome()

# Outliers: top 30 spend days overall
spend_long %>%
  filter(!is.na(spend)) %>%
  arrange(desc(spend)) %>%
  slice_head(n = 30) %>%
  ggplot(aes(x = reorder(paste(channel, day), spend), y = spend)) +
  geom_col() +
  coord_flip() +
  scale_y_continuous(labels = dollar_format()) +
  labs(title = "Top 30 Largest Daily Spends (Potential Outliers / Billing Dumps)",
    x = "Channel/Day",
    y = "Spend") +
  theme_fancy()

# correlation
library(tidyverse)

min_non_na <- 30

spend_wide <- df_num %>% select(all_of(spend_cols))

spend_zeroimp <- spend_wide %>%
  mutate(across(everything(), ~ replace_na(.x, 0)))

col_screen_0 <- spend_zeroimp %>%
  summarise(across(
    everything(),
    list(
      n_non_na = ~ sum(!is.na(.x)),
      sd = ~ sd(.x)
    ),
    .names = "{.col}__{.fn}"
  )) %>%
  pivot_longer(cols = everything(), names_to = "name", values_to = "value") %>%
  separate(name, into = c("channel", "metric"), sep = "__", extra = "merge", fill = "right") %>%
  pivot_wider(names_from = metric, values_from = value) %>%
  mutate(sd = replace_na(sd, 0))

keep_cols_0 <- col_screen_0 %>%
  filter(sd > 0) %>%
  pull(channel)

spend_zeroimp_ok <- spend_zeroimp %>%
  select(all_of(keep_cols_0))

corr_zeroimp <- suppressWarnings(
  cor(spend_zeroimp_ok, use = "pairwise.complete.obs")
)

if (is.null(rownames(corr_zeroimp)) || is.null(colnames(corr_zeroimp))) {
  rownames(corr_zeroimp) <- colnames(spend_zeroimp_ok)
  colnames(corr_zeroimp) <- colnames(spend_zeroimp_ok)
}

corr_long_zeroimp <- as.data.frame(corr_zeroimp, stringsAsFactors = FALSE) %>%
  rownames_to_column("channel1") %>%
  pivot_longer(-channel1, names_to = "channel2", values_to = "corr") %>%
  filter(channel1 < channel2) %>%
  filter(!is.na(corr), is.finite(corr))

ggplot(corr_long_zeroimp, aes(x = channel1, y = channel2, fill = corr)) +
  geom_tile() +
  scale_fill_gradient2(limits = c(-1, 1)) +
  labs(
    title = "Channel Spend Correlation",
    x = "Channel",
    y = "Channel",
    fill = "Correlation" ) +
  theme_fancy()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key = element_rect(color = NA, fill = NA))

# Average spend by weekday
weekday_summary %>%
  ggplot(aes(x = wday, y = avg_spend, group = channel)) +
  geom_line(na.rm = TRUE) +
  geom_point(na.rm = TRUE) +
  facet_wrap(~ channel, scales = "free_y", ncol = 2) +
  scale_y_continuous(labels = dollar_format()) +
  labs(
    title = "Average Spend by Weekday",
    x = "Weekday",
    y = "Average spend" ) +
  theme_awesome()

# print issue report
cat("\n====================\nMMM Spend Data QA Report\n====================\n")

cat("\n1) Column name hygiene issues (whitespace, etc.)\n")
if (nrow(name_whitespace_issues) == 0) {
  cat("  None detected.\n")
} else {
  print(name_whitespace_issues)
}

cat("\n2) Non-date rows detected in 'day' (likely summary rows)\n")
cat(glue("  Count: {nrow(non_date_rows)}\n"))
if (nrow(non_date_rows) > 0) {
  cat("  Examples:\n")
  print(non_date_rows %>% select(day_raw) %>% distinct() %>% slice_head(n = 10))
}

cat("\n3) Parsing failures after numeric conversion\n")
if (nrow(parsing_failures) == 0) {
  cat("  None detected.\n")
} else {
  print(parsing_failures)
}

cat("\n4) Date range coverage\n")
print(date_summary)

cat("\n5) Missing days (gaps in daily sequence)\n")
cat(glue("  Missing day count: {nrow(missing_days)}\n"))
if (nrow(missing_days) > 0) {
  cat("  First 10 missing days:\n")
  print(missing_days %>% slice_head(n = 10))
}

cat("\n6) Channel coverage / missing vs zero\n")
print(coverage_flags %>% select(
  channel, n_days, n_missing, pct_missing, n_zero, pct_zero, n_positive, pct_positive,
  first_non_missing, last_non_missing, first_positive, last_positive,
  appears_late, ends_early, extremely_sparse
))

cat("\n7) Outliers & extreme values summary\n")
print(outlier_table)

cat("\n8) Top correlated channel pairs (pairwise complete obs)\n")
print(corr_long_pairwise %>% slice_head(n = 10))

cat("\n9) Top correlated channel pairs (NA imputed as 0)\n")
print(corr_long_zeroimp %>% slice_head(n = 10))

cat("\n10) Weekend vs weekday spend ratio (avg weekend / avg weekday)\n")
print(weekend_vs_weekday)


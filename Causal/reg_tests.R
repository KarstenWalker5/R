
# This script implements basic tests that should be performed on all regression models

# Packages 
library(tidyverse)
library(broom)
library(nortest)     # lillie.test (Lilliefors)
library(stats)       # qqnorm / qqline (base) if you prefer
library(lmtest)      # bptest, dwtest
library(car)         # vif, outlierTest
library(performance) # optional: check_model, check_collinearity

# Load dataset 

if (!requireNamespace("rattle", quietly = TRUE)) {
  install.packages("rattle")
}
data(wine, package = "rattle")

df <- as_tibble(wine)

# rattle::wine uses "Type" as class label; convert to 0/1/2
df <- df %>%
  rename(target = Type) %>%
  mutate(target = as.integer(target) - 1)

# The python notebook renamed 'od280/od315_of_diluted_wines' -> 'test_diluted_wines'.
# In rattle::wine the column name differs; map if present.
# If you already have a column with that meaning, rename it here:
# df <- df %>% rename(test_diluted_wines = <your_column_name>)

df %>% dim()
df %>% glimpse()

# N/A check
df %>%
  summarise(across(everything(), ~ sum(is.na(.x)))) %>%
  pivot_longer(everything(), names_to = "column", values_to = "n_missing") %>%
  arrange(desc(n_missing))

# Regression model 
model <- lm(target ~ ., data = df)

summary(model)

broom::tidy(model, conf.int = TRUE)

broom::glance(model)

# Extract residuals 
model_resid <- resid(model)

# Normality tests

# Lilliefors / KS-normality

lillie <- nortest::lillie.test(model_resid)

lillie

cat(ifelse(lillie$p.value < 0.05, "Not normal", "Normal"),
    "| p-value:", lillie$p.value, "\n")

# Anderson–Darling 
ad <- nortest::ad.test(model_resid)
ad

# Residual histogram 
tibble(resid = model_resid) %>%
  ggplot(aes(resid)) +
  geom_histogram(bins = 30) +
  theme_minimal() +
  labs(title = "Residual histogram", x = "Residuals", y = "Count")

# QQ plot
augment(model) %>%
  ggplot(aes(sample = .std.resid)) +
  stat_qq() +
  stat_qq_line() +
  theme_minimal() +
  labs(title = "QQ Plot (standardized residuals)")

# Homoscedasticity (scatter + Breusch–Pagan)
# Residuals vs predicted scatter
augment(model) %>%
  ggplot(aes(x = .fitted, y = .resid)) +
  geom_hline(yintercept = 0) +
  geom_point(alpha = 0.7) +
  theme_minimal() +
  labs(
    title = "Residuals vs Predicted (fitted) values",
    x = "Predicted values",
    y = "Residuals"
  )

# Breusch–Pagan test (sms.het_breuschpagan)
bp <- lmtest::bptest(model)
bp
broom::tidy(bp)

# Outliers 
car::outlierTest(model)

# tidy table of influential points:
infl_tbl <- augment(model) %>%
  mutate(row = row_number()) %>%
  dplyr::select(row, .fitted, .resid, .std.resid, .hat, .cooksd) %>%
  arrange(desc(.cooksd))

infl_tbl %>% slice_head(n = 10)

# Independence Durbin–Watson
dw <- lmtest::dwtest(model)

dw

broom::tidy(dw)

# Multicollinearity (corr + VIF) 
# Correlation matrix 
variables <- df %>% 
  dplyr::select(-target)

# matrix form 
cor_mat <- cor(variables, use = "pairwise.complete.obs")
cor_mat

# tidy long form 
cor_long <- cor_mat %>%
  as.data.frame() %>%
  rownames_to_column("var1") %>%
  pivot_longer(-var1, names_to = "var2", values_to = "correlation")

cor_long %>% slice_head(n = 20)

# 11b) VIF (variance_inflation_factor)
vif_tbl <- car::vif(model)

vif_tbl

enframe(vif_tbl, name = "term", value = "vif") %>%
  arrange(desc(vif))

# Optional richer collinearity report
performance::check_collinearity(model) %>% as_tibble()

# Model summary + predictions column -
summary(model)

df <- df %>%
  mutate(predictions = fitted(model))

df %>% slice_tail(n = 20)

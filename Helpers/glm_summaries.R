library(broom)
library(dplyr)
library(stringr)
library(purrr)

#  create a tidy coefficient table with ORs + CIs 
glm_coef_table <- function(model, conf_level = 0.95) {
  broom::tidy(model, conf.int = TRUE, conf.level = conf_level) %>%
    mutate(
      odds_ratio = exp(estimate),
      or_low = exp(conf.low),
      or_high = exp(conf.high)
    ) %>%
    select(term, estimate, std.error, statistic, p.value, conf.low, conf.high, odds_ratio, or_low, or_high) %>%
    arrange(desc(abs(estimate)))
}

# Fast GLM export using Wald confidence intervals
glm_coef_table_fast <- function(model, conf_level = 0.95) {
  z <- qnorm(1 - (1 - conf_level) / 2)
  
  broom::tidy(model) %>%
    mutate(
      conf.low  = estimate - z * std.error,
      conf.high = estimate + z * std.error,
      odds_ratio = exp(estimate),
      or_low = exp(conf.low),
      or_high = exp(conf.high)
    ) %>%
    select(term, estimate, std.error, statistic, p.value,
           conf.low, conf.high, odds_ratio, or_low, or_high) %>%
    arrange(desc(abs(estimate)))
}

glm_fit_table <- function(model) {
  tibble::tibble(
    n = stats::nobs(model),
    df_null = model$df.null,
    df_residual = model$df.residual,
    null_deviance = model$null.deviance,
    residual_deviance = model$deviance,
    pseudo_r2 = 1 - (model$deviance / model$null.deviance),
    AIC = stats::AIC(model),
    BIC = stats::BIC(model)
  )
}

glm_intercept_prob <- function(model) {
  b0 <- unname(coef(model)["(Intercept)"])
  tibble::tibble(
    intercept = b0,
    baseline_prob = plogis(b0)
  )
}

glm_export_bundle_fast <- function(model) {
  list(
    coefficients = glm_coef_table_fast(model),
    fit = glm_fit_table(model),
    baseline = glm_intercept_prob(model)
  )
}

# reate a one-row model fit summary 
glm_fit_table <- function(model) {
  tibble::tibble(
    n = stats::nobs(model),
    df_null = model$df.null,
    df_residual = model$df.residual,
    null_deviance = model$null.deviance,
    residual_deviance = model$deviance,
    pseudo_r2 = 1 - (model$deviance / model$null.deviance),
    AIC = stats::AIC(model),
    BIC = stats::BIC(model)
  )
}

# create an "interpretation" table: baseline prob from intercept 
glm_intercept_prob <- function(model) {
  coefs <- coef(model)
  if (!("(Intercept)" %in% names(coefs))) return(tibble::tibble())
  tibble::tibble(
    intercept = unname(coefs["(Intercept)"]),
    baseline_prob = stats::plogis(unname(coefs["(Intercept)"]))
  )
}

#  wrapper: returns a named list of 2-3 dataframes to export 
glm_export_bundle <- function(model) {
  list(
    coefficients = glm_coef_table(model),
    fit = glm_fit_table(model),
    baseline = glm_intercept_prob(model)
  )
}

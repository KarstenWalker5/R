library(dplyr)
library(e1071)

###### Log1p transformation ######
log_transform_skewed <- function(df, skew_threshold = 1) {
  
  numeric_cols <- df %>%
    select(where(is.numeric)) %>%
    colnames()
  
  if (length(numeric_cols) == 0) return(df)
  
  # Only consider columns with:
  # - no NAs
  # - at least 3 distinct values
  # - finite values
  numeric_ok <- numeric_cols[
    sapply(df[numeric_cols], function(x) {
      !anyNA(x) &&
        all(is.finite(x)) &&
        length(unique(x)) >= 3
    })
  ]
  
  if (length(numeric_ok) == 0) {
    message("No eligible numeric columns (non-NA, finite, >=3 unique values). Returning df unchanged.")
    return(df)
  }
  
  # Compute skewness safely
  skew_vals <- sapply(df[numeric_ok], function(x) {
    e1071::skewness(x, type = 2)  # na.rm not needed because we've excluded NAs
  })
  
  # Drop NA/NaN/Inf skewness results (can happen for near-constant columns)
  skew_vals <- skew_vals[is.finite(skew_vals)]
  
  if (length(skew_vals) == 0) {
    message("Skewness could not be computed for any eligible columns. Returning df unchanged.")
    return(df)
  }
  
  skewed_vars <- names(skew_vals[skew_vals > skew_threshold])
  
  # Defensive: remove any accidental NAs and ensure they exist in df
  skewed_vars <- skewed_vars[!is.na(skewed_vars)]
  skewed_vars <- intersect(skewed_vars, colnames(df))
  
  if (length(skewed_vars) == 0) {
    message("No variables exceed skew_threshold. Returning df unchanged.")
    return(df)
  }
  
  df %>%
    mutate(across(all_of(skewed_vars), ~ log1p(.x), .names = "{.col}_log")) %>%
    select(-all_of(skewed_vars))
}


###### ECDF Transformation ######

add_ecdf_columns <- function(df, cols) {
  cols <- intersect(cols, names(df))
  if (length(cols) == 0) return(df)
  
  df %>%
    mutate(across(
      all_of(cols),
      ~ {
        x <- .x
        out <- rep(NA_real_, length(x))
        ok <- !is.na(x)
        # ECDF-like percentile rank; ties get average rank
        out[ok] <- rank(x[ok], ties.method = "average") / sum(ok)
        out
      },
      .names = "{.col}_ecdf"
    ))
}

###### Identify feature columns ######
get_feature_cols <- function(df,
                             id_cols = c("user_id", "week"),
                             cluster_col = "cluster",
                             retention_regex = "^retained") {
  setdiff(
    names(df),
    c(id_cols, cluster_col, names(df)[str_detect(names(df), retention_regex)])
  )
}

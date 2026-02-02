
# This script uses the declareDesign package for sample size calculations and fabricatr
# to create synthetic data
# Based on this post: https://arelbundock.com/posts/money_and_power/index.html

library(MASS)
library(tidyverse)
library(fabricatr)
library(patchwork)
library(DeclareDesign)
library(marginaleffects)

set.seed(123)

theme_fancy <- function() {
  theme_minimal(base_family = "Asap Condensed") +
    theme(panel.grid.minor = element_blank(),
          legend.position = "bottom")
}

# Simulate data
# In DeclareDesign, the parameters that we intend to vary—like the sample_size and treatment_effect variables—must be defined 
# before and outside the declare_model() call.
# declare_model creates a function that we can repeatedly use to simulate new datasets from the DGP.

sample_size = 100

treatment_effect = 0.2

dgp = declare_model(
  N = sample_size,
  # standard normal noise
  e = rnorm(N, mean = 0, sd = 1),
  # random dummy treatment with equal probability of treatment and control
  D = rbinom(N, size = 1, prob = 0.5),
  # intercept
  b0 = 0,
  # outcome equation
  Y = b0 + treatment_effect * D + e
)

# The fit() function defined below accepts a data frame, fits a linear regression model, 
# and returns a 1-row data frame as results.

fit = function(data) {
  model = lm(Y ~ D, data = data)
  out = data.frame(
    "term" = "OLS",
    "estimate" = coef(model)["D"]
  )
  return(out)
}

# simulate
d = dgp()

# estimate
fit(d)

# We can use the avg_comparisons() function from the marginaleffects package. In this simple model, 
# the estimate produced by avg_comparisons() will identical to the coefficient estimated by OLS. 
# In more complex models, with interactions or non-linear components, avg_comparisons() will be a 
# convenient way to estimate the average treatment effect (ATE), and to return results in the “tidy” 
# format required by DeclareDesign.

fit = function(data) {
  model = lm(Y ~ D, data = data)
  out = avg_comparisons(model, variables = "D")
  return(out)
}

fit(d)

# Now, we use DeclareDesign to create an R object that represents a complete “research design”. 
# This object is created by “adding” or “chaining” together calls to different DeclareDesign functions, 
# using the + operator. In particular, we will call three functions:
#     1. declare_model() (Data generating process)
#     2. declare_inquiry() (True value of the quantity of interest)
#     3. declare_estimator() (Estimation procedure)
# Finally, once we have defined the research design, we can call the diagnose_design() function to assess 
# its properties. This function will run simulations, estimate the treatment effect, and compute power, bias, etc.

sample_size = 100
treatment_effect = 0.2

dgp = declare_model(
  N = sample_size,
  e = rnorm(N, mean = 0, sd = 1),
  D = rbinom(N, size = 1, prob = 0.5),
  Y = treatment_effect * D + e
)

estimator = declare_estimator(handler = 
                                function(data) {
                                  model = lm(Y ~ D, data = data)
                                  out = avg_comparisons(model, variables = "D")
                                  return(out)
                                }
)

truth = declare_inquiry(ATE = treatment_effect)

design = dgp + estimator + truth

diagnose_design(design)

# To increase the power, we can increase the sample size, the treatment effect, or both.
# We use the redesign() function from DeclareDesign to create a list of new research designs, 
# with different combinations of effect and sample sizes. We then call the diagnose_design() 
# function to compare the power of the different research designs.

design_list = redesign(design,
                       sample_size = c(100, 500),
                       treatment_effect = c(0.2, 0.5)
)

diagnose_design(design_list)

###### Comparing estimators ######
sample_size = 100
treatment_effect = 0.2

dgp = declare_model(
  N = sample_size,
  e = rnorm(N, mean = 0, sd = 1),
  D = rbinom(N, size = 1, prob = 0.5),
  
  # random normal covariates
  fabricatr::draw_multivariate(c(X1, X2) ~ MASS::mvrnorm(
    n = N,
    mu = c(0, 0),
    Sigma = matrix(c(1, 0.5, 0.5, 1), 2, 2)
  )),
  
  Y = 1 + treatment_effect * D + 1 * X1 + 1 * X2 + e
)

fit = function(data) {
  # Linear model with regression adjustement for covariates
  adjusted =  lm(Y ~ D * (X1 + X2), data = data)
  adjusted = avg_comparisons(adjusted, variables = "D", vcov = "HC3")
  
  # Linear model with a single binary predictor
  unadjusted =  lm(Y ~ D, data = data)
  unadjusted = avg_comparisons(unadjusted, variables = "D", vcov = "HC3")
  
  # Combine the results
  out = rbind(adjusted, unadjusted)
  
  # Label the results
  out$term = c("Adjusted", "Unadjusted")
  return(out)
}

# We can fit this model directly to the output of dgp():
  
dgp() %>% fit()

# Diagnose the research design as before:
  
truth = declare_inquiry(ATE = treatment_effect)

estimators = declare_estimator(handler = fit)

design = dgp + truth + estimators

diagnose_design(design)
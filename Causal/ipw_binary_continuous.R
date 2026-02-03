
# Generating inverse probability weights for both binary and continuous treatments
# Source: https://www.andrewheiss.com/blog/2020/12/01/ipw-binary-continuous/
# NOTE: This script contains ONLY code from the article, assembled into one file.

library(tidyverse)
library(broom)
library(scales)
library(ggdag)
library(dagitty)
library(truncnorm)
library(ipw)
library(WeightIt)


mosquito_dag <- dagify(
  mal ~ net + inc + hlth,
  net ~ inc + hlth,
  hlth ~ inc,
  coords = list(
    x = c(mal = 4, net = 1, inc = 2, hlth = 3),
    y = c(mal = 1, net = 1, inc = 2, hlth = 2)
  ),
  exposure = "net",
  outcome = "mal"
)

ggdag_status(mosquito_dag) +
  guides(color = "none") +
  theme_dag()

# Make this randomness consistent
set.seed(1234)

# Simulate 1138 people (just for fun)
n_people <- 1138
net_data <- tibble(
  # Make an ID column (not necessary, but nice to have)
  id = 1:n_people,
  # Generate income variable: normal, 500 ± 300
  income = rnorm(n_people, mean = 500, sd = 75)
) %>%
  # Generate health variable: beta, centered around 70ish
  mutate(
    health_base = rbeta(n_people, shape1 = 7, shape2 = 4) * 100,
    # Health increases by 0.02 for every dollar in income
    health_income_effect = income * 0.02,
    # Make the final health score and add some noise
    health = health_base + health_income_effect + rnorm(n_people, mean = 0, sd = 3),
    # Rescale so it doesn't go above 100
    health = rescale(health, to = c(min(health), 100))
  ) %>%
  # Generate net variable based on income, health, and random noise
  mutate(
    net_score = (0.5 * income) + (1.5 * health) + rnorm(n_people, mean = 0, sd = 15),
    # Scale net score down to 0.05 to 0.95 to create a probability of using a net
    net_probability = rescale(net_score, to = c(0.05, 0.95)),
    # Randomly generate a 0/1 variable using that probability
    net = rbinom(n_people, 1, net_probability)
  ) %>%
  # Finally generate a malaria risk variable based on income, health, net use,
  # and random noise
  mutate(
    malaria_risk_base = rbeta(n_people, shape1 = 4, shape2 = 5) * 100,
    # Risk goes down by 10 when using a net. Because we rescale things,
    # though, we have to make the effect a lot bigger here so it scales
    # down to -10. Risk also decreases as health and income go up. I played
    # with these numbers until they created reasonable coefficients.
    malaria_effect = (-30 * net) + (-1.9 * health) + (-0.1 * income),
    # Make the final malaria risk score and add some noise
    malaria_risk = malaria_risk_base + malaria_effect + rnorm(n_people, 0, sd = 3),
    # Rescale so it doesn't go below 0,
    malaria_risk = rescale(malaria_risk, to = c(5, 70))
  ) %>%
  select(-c(
    health_base, health_income_effect, net_score, net_probability,
    malaria_risk_base, malaria_effect
  ))

head(net_data)

# Wrong correlation-is-not-causation effect
model_net_naive <- lm(malaria_risk ~ net, data = net_data)
tidy(model_net_naive)

adjustmentSets(mosquito_dag)

# Logit model to predict net use
model_predict_net <- glm(
  net ~ income + health,
  family = binomial(link = "logit"),
  data = net_data
)

# Generate propensity scores and IPWs
net_data_ipw <- augment_columns(model_predict_net, net_data, type.predict = "response") %>%
  rename(propensity = .fitted) %>%
  mutate(ipw = (net / propensity) + ((1 - net) / (1 - propensity)))

net_data_ipw %>%
  select(id, income, health, net, malaria_risk, propensity, ipw) %>%
  head()

model_net_ipw <- lm(malaria_risk ~ net, data = net_data_ipw, weights = ipw)
tidy(model_net_ipw)

# IPW with the ipw package, binary treatment
# ipwpoint() can't handle tibbles! Force net_data to be a data.frame
weights_ipwpoint <- ipwpoint(
  exposure = net,
  family = "binomial",  # The treatment is binary
  link = "logit",
  denominator = ~ income + health,
  data = as.data.frame(net_data)
)

# They're the same!
head(weights_ipwpoint$ipw.weights)
head(net_data_ipw$ipw)

net_data_ipwpoint <- net_data %>%
  mutate(ipw = weights_ipwpoint$ipw.weights)

model_net_ipwpoint <- lm(malaria_risk ~ net,
                         data = net_data_ipwpoint, weights = ipw)
tidy(model_net_ipwpoint)

# IPW with the WeightIt package, binary treatment
weights_weightit <- weightit(net ~ income + health,  # Model net use with confounders
                             data = net_data,
                             estimand = "ATE",  # Find the ATE
                             method = "ps")  # Build weights with propensity scores
weights_weightit

# See even more details here
# summary(weights_weightit)

# Same as the other methods!
head(weights_weightit$weights)

net_data_weightit <- net_data %>%
  mutate(ipw = weights_weightit$weights)

model_net_weightit <- lm(malaria_risk ~ net,
                         data = net_data_weightit, weights = ipw)
tidy(model_net_weightit)

# ----------------------------------------------------------------------------
# Continuous treatments
# ----------------------------------------------------------------------------

grant_dag <- dagify(
  mal ~ grant + inc + hlth,
  grant ~ inc + hlth,
  hlth ~ inc,
  coords = list(
    x = c(mal = 4, grant = 1, inc = 2, hlth = 3),
    y = c(mal = 1, grant = 1, inc = 2, hlth = 2)
  ),
  exposure = "grant",
  outcome = "mal"
)

ggdag_status(grant_dag) +
  guides(color = "none") +
  theme_dag()

# Make this randomness consistent
set.seed(1234)

# Simulate 1504 people
n_people <- 1504
grant_data <- tibble(
  # Make an ID column (not necessary, but nice to have)
  id = 1:n_people,
  # Generate income variable: normal, 500 ± 300
  income = rnorm(n_people, mean = 500, sd = 75)
) %>%
  # Generate health variable: beta, centered around 70ish
  mutate(
    health_base = rbeta(n_people, shape1 = 7, shape2 = 4) * 100,
    # Health increases by 0.02 for every dollar in income
    health_income_effect = income * 0.02,
    # Make the final health score and add some noise
    health = health_base + health_income_effect + rnorm(n_people, mean = 0, sd = 3),
    # Rescale so it doesn't go above 100
    health = rescale(health, to = c(min(health), 100))
  ) %>%
  # Generate grant variable
  mutate(
    grant_base = rtruncnorm(n_people, mean = 18, sd = 10, a = 5, b = 40),
    # Grants are higher for people with lower incomes; higher for people with lower health
    grant_effect = (income * -0.25) + (health * -0.5),
    # Make the final grant amount + noise + rescale it back down
    grant = grant_base + grant_effect + rnorm(n_people, mean = 0, sd = 8),
    grant = round(rescale(grant, to = c(5, 40)), 0)
  ) %>%
  # Finally generate a malaria risk variable based on income, health, grant amount,
  # and random noise
  mutate(
    malaria_risk_base = rbeta(n_people, shape1 = 4, shape2 = 5) * 100,
    # Risk goes down as grant money goes up. I played with these numbers
    # until they created reasonable coefficients.
    malaria_effect = (-40 * grant) + (-25 * health) + (-0.05 * income),
    # Make the final malaria risk score and add some noise
    malaria_risk = malaria_risk_base + malaria_effect + rnorm(n_people, 0, sd = 3),
    # Rescale so it doesn't go below 0,
    malaria_risk = rescale(malaria_risk, to = c(5, 70))
  ) %>%
  select(-c(
    health_base, health_income_effect, grant_base, grant_effect,
    malaria_risk_base, malaria_effect
  ))

head(grant_data)

# Wrong correlation-is-not-causation effect
model_grant_naive <- lm(malaria_risk ~ grant, data = grant_data)
tidy(model_grant_naive)

adjustmentSets(grant_dag)

# The numerator is the probability distribution of just the treatment variable.
# We'll use a normal distribution for it (hence dnorm()). We need to feed
# dnorm() the grant amount for each person, the predicted value from a simple
# grant ~ 1 model, and the sd of the residuals from that model
model_num <- lm(grant ~ 1, data = grant_data)
num <- dnorm(grant_data$grant,
             predict(model_num),
             sd(model_num$residuals))

# The denominator is the probability distribution of the treatment variable
# explained by the confounders. We'll again use a normal distribution for it.
# We'll feed dnorm() the grant amount, the predicted value from a model that
# includes the confounders, and the sd of the residuals from that model
model_den <- lm(grant ~ health + income, data = grant_data)
den <- dnorm(grant_data$grant,
             predict(model_den),
             sd(model_den$residuals))

# Finally, we make actual IPW weights by building the fraction
grant_data_ipw <- grant_data %>%
  mutate(ipw = num / den)

head(grant_data_ipw)

model_grant_ipw <- lm(malaria_risk ~ grant, data = grant_data_ipw, weights = ipw)
tidy(model_grant_ipw)

# IPW with the ipw package, continuous treatment
weights_continuous_ipwpoint <- ipwpoint(
  exposure = grant,
  family = "gaussian",
  numerator = ~ 1,
  denominator = ~ health + income,
  data = as.data.frame(grant_data)
)

# Same values!
head(grant_data_ipw$ipw)
head(weights_continuous_ipwpoint$ipw.weights)

grant_data_ipwpoint <- grant_data %>%
  mutate(ipw = weights_continuous_ipwpoint$ipw.weights)

model_grant_ipwpoint <- lm(malaria_risk ~ grant,
                           data = grant_data_ipwpoint, weights = ipw)
tidy(model_grant_ipwpoint)

# IPW with the WeightIt package, continuous treatment
weights_weightit <- weightit(grant ~ income + health,  # Model grant amount with confounders
                             data = grant_data,
                             stabilize = TRUE)
weights_weightit

# See even more details here
# summary(weights_weightit)

# Not the same as the other methods :(
head(weights_weightit$weights)

# Manual weights
head(grant_data_ipw$ipw)

# weightit() weights / 2
head(weights_weightit$weights) / 2

grant_data_weightit <- grant_data %>%
  mutate(ipw = weights_weightit$weights)

model_grant_weightit <- lm(malaria_risk ~ grant,
                           data = grant_data_weightit, weights = ipw)
tidy(model_grant_weightit)

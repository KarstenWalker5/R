# Load libraries
library(tidyverse)  # ggplot, dplyr, and friends
library(ggridges)   # Ridge plots
library(ggstance)   
library(patchwork)  
library(scales)     
library(broom)  
library(broom.mixed)
library(infer)
library(rstan)
library(brms) 
library(ggplot2movies) 

# Reference post: https://www.andrewheiss.com/blog/2019/01/29/diff-means-half-dozen-ways/

###### Data ######
# Clean up data
set.seed(7)

movies_clean <- movies %>% 
  # Make a binary column for genre
  select(title, year, rating, Action, Comedy) %>% 
  filter(!(Action == 1 & Comedy == 1)) %>% 
  mutate(genre = case_when(Action == 1 ~ "Action",
                           Comedy == 1 ~ "Comedy",
                           TRUE ~ "Neither")) %>%
  filter(genre != "Neither") %>%
  # Make a numeric version of genre, where action = 1, comedy = 2
  mutate(genre_numeric = as.numeric(factor(genre))) %>% 
  # Make genre a factor
  mutate(genre = factor(genre)) %>% 
  select(-Action, -Comedy) %>% 
  # Randomly select 200 movies in each genre
  group_by(genre) %>% 
  sample_n(200) %>% 
  ungroup()

###### Visualize Data ######
# Theme
theme_fancy <- function() {
  theme_minimal(base_family = "Asap Condensed") +
    theme(panel.grid.minor = element_blank())
}

eda_boxplot <- ggplot(movies_clean, aes(x = genre, y = rating, fill = genre)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#0288b7", "#a90010"), guide = FALSE) + 
  scale_y_continuous(breaks = seq(1, 10, 1)) +
  labs(x = NULL, y = "Rating") +
  theme_fancy()

eda_histogram <- ggplot(movies_clean, aes(x = rating, fill = genre)) +
  geom_histogram(binwidth = 1, color = "white") +
  scale_fill_manual(values = c("#0288b7", "#a90010"), guide = FALSE) + 
  scale_x_continuous(breaks = seq(1, 10, 1)) +
  labs(y = "Count", x = "Rating") +
  facet_wrap(~ genre, nrow = 2) +
  theme_fancy() +
  theme(panel.grid.major.x = element_blank())

eda_ridges <- ggplot(movies_clean, aes(x = rating, y = fct_rev(genre), fill = genre)) +
  stat_density_ridges(quantile_lines = TRUE, quantiles = 2, scale = 3, color = "white") + 
  scale_fill_manual(values = c("#0288b7", "#a90010"), guide = FALSE) + 
  scale_x_continuous(breaks = seq(0, 10, 2)) +
  labs(x = "Rating", y = NULL,
       subtitle = "White line shows median rating") +
  theme_fancy()

(eda_boxplot | eda_histogram) / 
  eda_ridges + 
  plot_annotation(title = "Do comedies get higher ratings than action movies?",
                  subtitle = "Sample of 400 movies from IMDB",
                  theme = theme(text = element_text(family = "Asap Condensed"),
                                plot.title = element_text(face = "bold",
                                                          size = rel(1.5))))

group_diffs <- movies_clean %>% 
  group_by(genre) %>% 
  summarize(avg_rating = mean(rating, na.rm = TRUE)) %>% 
  mutate(diff_means = avg_rating - lead(avg_rating))

group_diffs

###### Frequentist T-test ######
t_test_eq <- t.test(rating ~ genre, data = movies_clean, var.equal = TRUE)

t_test_eq

t_test_eq_tidy <- tidy(t_test_eq) %>% 
  # Calculate difference in means, since t.test() doesn't actually do that
  mutate(estimate = estimate1 - estimate2) %>%
  select(starts_with("estimate"), everything())

t_test_eq_tidy

## Check if groups have equal variance

# Bartlett's test, homogeneity of variances based on the mean
bartlett.test(rating ~ genre, data = movies_clean)

# Levene's test, homogeneity of variances based on the median, more robust to outliers
car::leveneTest(rating ~ genre, data = movies_clean)

# Fligner-Killeen test: Check homogeneity of variances based on the median, so it’s more robust to outliers
fligner.test(rating ~ genre, data = movies_clean)

# Kruskal-Wallis test: Check homogeneity of distributions nonparametrically
kruskal.test(rating ~ genre, data = movies_clean)

# t-test, assuming unequal variance

t_test_uneq <- t.test(rating ~ genre, data = movies_clean)

t_test_uneq_tidy <- tidy(t_test_uneq) %>% 
  mutate(estimate = estimate1 - estimate2) %>% 
  select(starts_with("estimate"), everything())

t_test_uneq_tidy

###### Simulation-based tests ######
# Instead of dealing with all the assumptions of the data and finding the exact statistical test written by some dude 
# in the 1940s, we can use the power of bootstrapping, permutation, and simulation to construct a null distribution 
# and calculate confidence intervals. According to Allen Downey, there is actually only one statistical test and that 
# at their core, all statistical tests follow the same universal pattern:
# Step 1: Calculate a sample statistic, or . This is the main measure you care about: the difference in means, the 
#       average, the median, the proportion, the difference in proportions, the chi-squared value, etc.
# Step 2: Use simulation to invent a world where  is null. Simulate what the world would look like if there was 
#       no difference between two groups, or if there was no difference in proportions, or where the average value is a specific number.
# Step 3: Look at  in the null world. Put the sample statistic in the null world and see if it fits well.
# Step 4: Calculate the probability that  could exist in null world. This is the p-value, or the probability 
#       that you’d see a  at least that high in a world where there’s no difference.
# Step 5: Decide if  is statistically significant. Choose some evidentiary standard or threshold (like 0.05) 
#       for deciding if there’s sufficient proof for rejecting the null world.

# Calculate the difference in means
diff_means <- movies_clean %>% 
  specify(rating ~ genre) %>%
  # Order here means we subtract comedy from action (Action - Comedy)
  calculate("diff in means", order = c("Action", "Comedy"))

diff_means

# Generate a bootstrapped distribution of the difference in means based on our sample and calculate 
# the confidence interval:
  
boot_means <- movies_clean %>% 
  specify(rating ~ genre) %>% 
  generate(reps = 1000, type = "bootstrap") %>% 
  calculate("diff in means", order = c("Action", "Comedy"))

boostrapped_confint <- boot_means %>% get_confidence_interval()

boot_means %>% 
  visualize() + 
  shade_confidence_interval(boostrapped_confint,
                            color = "#8bc5ed", fill = "#85d9d2") +
  geom_vline(xintercept = diff_means$stat, size = 1, color = "#77002c") +
  labs(title = "Bootstrapped distribution of differences in means",
       x = "Action − Comedy", y = "Count",
       subtitle = "Red line shows observed difference; shaded area shows 95% confidence interval") +
  theme_fancy()


# Step 2: Invent a world where δ is null
genre_diffs_null <- movies_clean %>% 
  specify(rating ~ genre) %>% 
  hypothesize(null = "independence") %>% 
  generate(reps = 5000, type = "permute") %>% 
  calculate("diff in means", order = c("Action", "Comedy"))

# Step 3: Put actual observed δ in the null world and see if it fits
genre_diffs_null %>% 
  visualize() + 
  geom_vline(xintercept = diff_means$stat, size = 1, color = "#77002c") +
  scale_y_continuous(labels = comma) +
  labs(x = "Simulated difference in average ratings (Action − Comedy)", y = "Count",
       title = "Simulation-based null distribution of differences in means",
       subtitle = "Red line shows observed difference") +
  theme_fancy()

# Step 4: Calculate probability that observed δ could exist in null world
genre_diffs_null %>% 
  get_p_value(obs_stat = diff_means, direction = "both") %>% 
  mutate(p_value_clean = pvalue(p_value))

# KS test
ks_res <- sales_data %>%
  select(treat, sales) %>%
  mutate(treat = factor(treat)) %>%
  { ks.test(.$sales[.$treat == 1],
            .$sales[.$treat == 0]) }

ks_res

ecdf_df <- sales_data %>%
  mutate(treat = ifelse(treat == 1, "Treated", "Control")) %>%
  group_by(treat) %>%
  arrange(sales) %>%
  mutate(ecdf = ecdf(sales)(sales))

## Plotting KS stat
# Create a grid of unique sales values
grid <- sales_data %>% 
  distinct(sales) %>% 
  arrange(sales)

# Compute ECDF for both groups on the same grid
grid <- grid %>%
  mutate(
    ecdf_treated = ecdf(sales_data$sales[sales_data$treat == 1])(sales),
    ecdf_control = ecdf(sales_data$sales[sales_data$treat == 0])(sales),
    diff = abs(ecdf_treated - ecdf_control)
  )

# Find KS point
ks_point <- grid %>% 
  slice_max(diff, n = 1)

ggplot() +
  geom_step(data = ecdf_df,
            aes(x = sales, y = ecdf, color = treat),
            size = 1) +
  geom_segment(
    data = ks_point,
    aes(x = sales, xend = sales,
        y = ecdf_treated, yend = ecdf_control),
    color = "black", linetype = "dashed", linewidth = 1
  ) +
  annotate("text",
           x = ks_point$sales,
           y = max(ks_point$ecdf_treated, ks_point$ecdf_control),
           label = paste0("KS = ", round(ks_point$diff, 3)),
           vjust = -0.5, hjust = -0.2,
           size = 5) +
  labs(
    title = "Kolmogorov–Smirnov Test Visualization",
    subtitle = "ECDF Comparison: Treated vs Control",
    x = "Sales",
    y = "ECDF",
    color = "Group"
  ) +
  theme_fancy() +
  scale_color_manual(values = c("Control" = "#1f78b4",
                                "Treated" = "#e31a1c"))

###### Bayesian Regression ######

# Frequentist null hypothesis significance testing (NHST) determines the probability of the data given 
# a null hypothesis (i.e. , yielding results that are often unwieldy, phrased as the probability of rejecting 
# the null if it is true (hence all that talk of “null worlds” earlier). In contrast, Bayesian analysis 
# determines the probability of a hypothesis given the data (i.e.), resulting in probabilities that are directly 
# interpretable.

# Set Monte Carlo params:
CHAINS <- 4

ITER <- 2000

WARMUP <- 1000

BAYES_SEED <- 1234

options(mc.cores = parallel::detectCores())  

# $un the model. We set a couple priors for the intercept (normal, with a mean of 0 and standard deviation of 5) 
# and the beta/slope (normal, with a mean of 0 and a standard deviation of 1). 
# In this model, rating ~ genre, we assume equal variances between the groups.
# We use tidyMCMC from broom to calculate the medians of the posterior distributions and create confidence intervals.

brms_eq <- brm(
  # bf() is an alias for brmsformula() and lets you specify model formulas
  bf(rating ~ genre), 
  # Reverse the levels of genre so that comedy is the base case
  data = mutate(movies_clean, genre = fct_rev(genre)),
  prior = c(set_prior("normal(0, 5)", class = "Intercept"),
            set_prior("normal(0, 1)", class = "b")),
  chains = CHAINS, iter = ITER, warmup = WARMUP, seed = BAYES_SEED
)

brms_eq_tidy <- 
  tidyMCMC(brms_eq, conf.int = TRUE, conf.level = 0.95, 
           estimate.method = "median", conf.method = "HPDinterval")

brms_eq_tidy

# The intercept represents the mean comedy score, while the coefficient for action represents the difference
# from that mean, or the effect that we care about. These findings match what we got with the frequentist 
# work from earlier, only now we can change the interpretation: We’re 95% confident that the true population-level 
# difference in rating is between -0.968 and -0.374, with a median of -0.666.

# Potential Priors
get_prior(bf(rating ~ genre), data = movies_clean)

# To specify the prior distribution for the intercept, you can use set_prior("normal(0, 5)", class = "Intercept"), 
# which matches the table. More complicated formulas will have values in the dpar and nlpar columns, 
# and you can use those to drill down to specific priors for those terms (e.g. set_prior("cauchy(0, 1)", 
# class = "b", dpar = "sigma")).

# Also, you can plot the different distributions to get a sense of their shapes with either base R or with ggplot:

# Normal distribution: normal(0, 5)
curve(expr = dnorm(x, mean = 0, sd = 5), from = -20, to = 20)

# Cauchy distribution: cauchy(0, 1)
curve(expr = dcauchy(x, location = 0, scale = 1), from = -5, to = 5)

# ggplot
norm_ggplot <- ggplot(data = tibble(x = c(-20, 20)), aes(x = x)) +
  stat_function(fun = dnorm, n = 500, args = list(mean = 0, sd = 5)) +
  labs(title = "normal(0, 5)") +
  theme_fancy()

cauchy_ggplot <- ggplot(data = tibble(x = c(-5, 5)), aes(x = x)) +
  stat_function(fun = dcauchy, n = 500, args = list(location = 0, scale = 1)) +
  labs(title = "cauchy(0, 1)") +
  theme_fancy()

norm_ggplot / cauchy_ggplot

## Regression, assuming unequal variances
brms_uneq <- brm(
  bf(rating ~ genre, sigma ~ genre), 
  data = mutate(movies_clean, genre = fct_rev(genre)),
  prior = c(set_prior("normal(0, 5)", class = "Intercept"),
            set_prior("normal(0, 1)", class = "b"),
            set_prior("cauchy(0, 1)", class = "b", dpar = "sigma")),
  chains = CHAINS, iter = ITER, warmup = WARMUP, seed = BAYES_SEED
)

brms_uneq_tidy <- 
  tidyMCMC(brms_uneq, conf.int = TRUE, conf.level = 0.95, 
           estimate.method = "median", conf.method = "HPDinterval")

# sigma terms are on a log scale, so we need to exponentiate them back to the scale of the data.

brms_uneq_tidy %>% 
  mutate_at(vars(estimate, std.error, conf.low, conf.high),
            funs(ifelse(str_detect(term, "sigma"), exp(.), .)))

brms_uneq_tidy

# Regression, Bayesian Estimation Supersedes the t Test (BEST) 
# Post: https://psycnet.apa.org/doi/10.1037/a0029146
# Explanation: https://vuorre.netlify.com/post/2017/01/02/how-to-compare-two-groups-with-robust-bayesian-estimation-using-r-stan-and-brms/#robust-bayesian-estimation

# In the most simplest terms, the only difference between BEST and the unequal variance regression 
# above is that we model the data with a t distribution, which means we have a new parameter, (nu), 
# that changes the normality of the distribution (i.e. the degrees of freedom parameter in a t distribution). 
# Kruschke uses an exponential prior with a rate of 1/29 in his paper, so we do too. It looks like this:

ggplot(data = tibble(x = c(0, 100)), aes(x = x)) +
  stat_function(fun = dexp, n = 500, args = list(rate = 1/29)) +
  labs(title = "exponential(1/29)") +
  theme_fancy()

# Add family = student and set the prior for nu and we’re ready to go:
  
brms_uneq_robust <- brm(
    bf(rating ~ genre, sigma ~ genre), 
    family = student,
    data = mutate(movies_clean, genre = fct_rev(genre)),
    prior = c(set_prior("normal(0, 5)", class = "Intercept"),
              set_prior("normal(0, 1)", class = "b"),
              set_prior("cauchy(0, 1)", class = "b", dpar = "sigma"),
              set_prior("exponential(1.0/29)", class = "nu")),
    chains = CHAINS, iter = ITER, warmup = WARMUP, seed = BAYES_SEED
  )
  
brms_uneq_robust_tidy <- 
    tidyMCMC(brms_uneq_robust, conf.int = TRUE, conf.level = 0.95, 
             estimate.method = "median", conf.method = "HPDinterval") %>% 
    # Rescale sigmas
    mutate_at(vars(estimate, std.error, conf.low, conf.high),
              funs(ifelse(str_detect(term, "sigma"), exp(.), .)))

  brms_uneq_robust_tidy
  
# When specifying a model as y ~ a + bx as we did above with rating ~ genre, the coefficient for genre (b) 
# is actually the difference in means and doesn’t directly reflect the rating itself. In theory, 
# if we’re thinking about two groups with two different variances, we should model the distribution of each group, 
# not the distribution of the differences in groups. Analyzing the distributions of the two groups separately 
# and then calculating the difference should yield more transparent results. 

# Here’s an intercept-free version of the brms-based BEST regression from earlier. Note how we’re modeling both 
# the rating and the group sigmas without intercepts (rating ~ 0 + genre, sigma ~ 0 + genre), and that we no 
# longer specify a prior for the intercept (if we do, brm() yells at us). Also note that instead of modeling 
# the beta coefficient as a normal distribution centered around zero (since that represented the difference in means), 
# we specify the distribution of action and comedy ratings themselves. Because we’re dealing with actual ratings, 
# we can make them fairly well informed and constrained. For instance, no movie is rated below 1 or above 10, 
# and I’m guessing from past experience looking at ratings on Amazon and IMDB and elsewhere that people tend to 
# inflate their ratings. I’d guess that the distribution of ratings looks something like this: normally distributed 
# with a mean of 6, standard deviation of 2, and truncated at 1 and 10.
  
# the msm library has rtnorm, dtnorm, and family
# which let you plot and draw from a truncated normal distribution
# msm::rtnorm(1000, mean = 6, sd = 2, lower = 1, upper = 10)
ggplot(data = tibble(x = c(1, 10)), aes(x = x)) +
    stat_function(fun = dnorm, n = 500, args = list(mean = 6, sd = 2)) +
    labs(title = "normal(6, 2); truncated at 1 and 10") +
    theme_fancy()

# Run model as normal
brms_uneq_robust_groups <- brm(
  bf(rating ~ 0 + genre, sigma ~ 0 + genre), 
  family = student,
  data = mutate(movies_clean, genre = fct_rev(genre)),
  prior = c(
    # Set group mean prior
    set_prior("normal(6, 2)", class = "b", lb = 1, ub = 10),
    # Ser group variance priors. We keep the less informative cauchy(0, 1).
    set_prior("cauchy(0, 1)", class = "b", dpar = "sigma"),
    set_prior("exponential(1.0/29)", class = "nu")),
  chains = CHAINS, iter = ITER, warmup = WARMUP, seed = BAYES_SEED
)

# Exponentiate Sigma Coefficients
brms_uneq_robust_groups_tidy <- 
  tidyMCMC(brms_uneq_robust_groups, conf.int = TRUE, conf.level = 0.95, 
           estimate.method = "median", conf.method = "HPDinterval") %>% 
  # Rescale sigmas
  mutate_at(vars(estimate, std.error, conf.low, conf.high),
            funs(ifelse(str_detect(term, "sigma"), exp(.), .)))

# Because we calculated the group means themselves, we need to do an extra few steps to get the difference in means. 
# It’s fairly easy: we extract the posterior samples for each of the groups, subtract them from each other, 
# and then calculate the credible interval.

brms_uneq_robust_groups_post <- posterior_samples(brms_uneq_robust_groups) %>% 
  # We can exponentiate here!
  mutate_at(vars(contains("sigma")), funs(exp)) %>% 
  # For whatever reason, we need to log nu?
  mutate(nu = log10(nu)) %>% 
  mutate(diff_means = b_genreAction - b_genreComedy,
         diff_sigma = b_sigma_genreAction - b_sigma_genreComedy) %>% 
  # Calculate effect sizes, just for fun
  mutate(cohen_d = diff_means / sqrt((b_sigma_genreAction + b_sigma_genreComedy)/2),
         cles = dnorm(diff_means / sqrt((b_sigma_genreAction + b_sigma_genreComedy)), 0, 1))

brms_uneq_robust_groups_tidy_fixed <- 
  tidyMCMC(brms_uneq_robust_groups_post, conf.int = TRUE, conf.level = 0.95, 
           estimate.method = "median", conf.method = "HPDinterval")

brms_uneq_robust_groups_tidy_fixed

brms_uneq_robust_groups_tidy

# Compare all methods on one plot

# Make a bunch of data frames that have three columns: 
# estimate, conf.low, and conf.high

# Extract t-test results
t_test_eq_small <- t_test_eq_tidy %>% 
  select(estimate, conf.low, conf.high)

t_test_uneq_small <- t_test_uneq_tidy %>% 
  select(estimate, conf.low, conf.high)

# Extract simulation results
infer_simulation <- tibble(estimate = diff_means$stat,
                           conf.low = boostrapped_confint$`2.5%`,
                           conf.high = boostrapped_confint$`97.5%`)

# Extract brms regression results
brms_eq_small <- brms_eq_tidy %>% 
  filter(term == "b_genreAction") %>% 
  select(estimate, conf.low, conf.high)

brms_uneq_small <- brms_uneq_tidy %>% 
  filter(term == "b_genreAction") %>% 
  select(estimate, conf.low, conf.high)

brms_uneq_robust_small <- brms_uneq_robust_tidy %>% 
  filter(term == "b_genreAction") %>% 
  select(estimate, conf.low, conf.high)

brms_uneq_robust_groups_small <- brms_uneq_robust_groups_tidy_fixed %>% 
  filter(term == "diff_means") %>% 
  select(estimate, conf.low, conf.high)

# Put all these mini dataframes into a list column, then unnest
meta_diffs <- tribble(
  ~package, ~method, ~results,
  "t-test", "equal variances", t_test_eq_small,
  "t-test", "unequal variances", t_test_uneq_small,
  "infer", "simulation", infer_simulation,
  "brms", "equal variances", brms_eq_small,
  "brms", "unequal variances", brms_uneq_small,
  "brms", "BEST", brms_uneq_robust_small,
  "brms", "BEST; group means", brms_uneq_robust_groups_small
) %>% 
  unnest(results) %>% 
  mutate(method = paste0(package, ": ", method)) %>% 
  mutate(method = fct_inorder(method))

ggplot(meta_diffs, aes(x = estimate, y = fct_rev(method), color = package)) +
  geom_pointrangeh(aes(xmin = conf.low, xmax = conf.high), size = 1) +
  geom_vline(xintercept = 0, size = 1) +
  scale_color_viridis_d(option = "plasma", end = 0.9, guide = FALSE) +
  labs(x = "Mean rating for action movies − mean rating for comedies",
       y = NULL, caption = "Sample of 400 movies from IMDB",
       title = "Comedies get higher ratings than action movies",
       subtitle = "Effect is roughly the same regardless of method used") +
  expand_limits(x = 0) +
  theme_fancy() +
  theme(plot.title = element_text(face = "bold", size = rel(1.5)))

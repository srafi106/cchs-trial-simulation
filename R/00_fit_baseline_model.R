# Fits the Bayesian multilevel model for the probability of a cardiovascular
# condition (heart disease, heart attack or effects of stroke), ages 35+.
# The model is used in 01_build_virtual_population.R to give each covariate
# cell a baseline risk. Needs brms and rstan. The fit is cached in fit_dir.

library(brms)
library(dplyr)
if (!exists("fit_dir")) source("config.R")
source("R/functions_data.R")

d = recode_pumf(read_pumf(pumf_path))
dat = analytic_sample(d)
dat$w = dat$WTS_M / mean(dat$WTS_M)
cat("Analytic sample:", nrow(dat), "respondents\n")

# The weighted Bernoulli likelihood only depends on the summed weights of
# cases and non-cases in each cell, so the model is fitted on cells.
agg = dat %>%
  group_by(hr, prov, age_grp, sex, inc_q) %>%
  summarise(w1 = sum(w * cvd), w0 = sum(w * (1 - cvd)), .groups = "drop") %>%
  tidyr::pivot_longer(c(w1, w0), names_to = "y", values_to = "w") %>%
  mutate(cvd = ifelse(y == "w1", 1L, 0L)) %>%
  filter(w > 0) %>%
  as.data.frame()
stopifnot(isTRUE(all.equal(sum(agg$w), sum(dat$w))))
cat("Cell records:", nrow(agg), "\n")

priors = c(
  prior(normal(-2, 1.5), class = Intercept),
  prior(normal(0, 1.5), class = b),
  prior(exponential(1), class = sd)
)

fit = brm(
  cvd | weights(w) ~ age_grp * sex + inc_q + (1 | prov) + (1 | hr),
  data = agg, family = bernoulli(), prior = priors,
  chains = 4, iter = 2000, warmup = 1000, cores = 4, seed = 2026,
  control = list(adapt_delta = 0.95), backend = "rstan",
  file = file.path(fit_dir, fit_name), file_refit = "never"
)
print(summary(fit))

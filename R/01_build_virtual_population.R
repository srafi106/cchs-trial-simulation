# Builds the virtual patient population for the simulations.
# Target population: adults aged 50 or older. Each cell (health region x age
# group x sex x income quintile) gets a weighted population size and a baseline
# probability of a cardiovascular condition.
#
# Output (in output/simulation):
#   virtual_population.csv  one row per cell
#   external_info.csv       weighted prevalence and Kish effective sample size

library(dplyr)
if (!exists("sim_dir")) source("config.R")
source("R/functions_data.R")

d = recode_pumf(read_pumf(pumf_path))
dat = analytic_sample(d)
dat50 = dat[dat$DHHGAGE >= 4, ]
cat("Respondents aged 50+ in the analytic sample:", nrow(dat50), "\n")

frame = dat50 %>%
  group_by(hr, prov, age_grp, sex, inc_q) %>%
  summarise(N = sum(WTS_M), n_resp = n(), .groups = "drop") %>%
  as.data.frame()

fit_path = file.path(fit_dir, paste0(fit_name, ".rds"))
if (file.exists(fit_path) && requireNamespace("brms", quietly = TRUE)) {
  message("Cell probabilities from the fitted Bayesian model")
  fit = readRDS(fit_path)
  nd = transform(frame, w = 1, cvd = 0L)
  P = brms::posterior_epred(fit, newdata = nd, allow_new_levels = FALSE)
  frame$p_prev = colMeans(P)
} else {
  message("No fitted model found; using a survey-weighted logistic regression")
  dat50$w = dat50$WTS_M / mean(dat50$WTS_M)
  m = glm(cvd ~ age_grp * sex + inc_q + prov, data = dat50,
          family = quasibinomial(), weights = w)
  frame$p_prev = predict(m, newdata = frame, type = "response")
}

w = dat50$WTS_M
ext = data.frame(
  p_ext_weighted = weighted.mean(dat50$cvd, w),
  n_eff_kish = sum(w)^2 / sum(w^2),
  n_resp = nrow(dat50),
  pop_total = sum(w)
)
print(ext)
cat("Weighted mean of the cell probabilities:",
    round(weighted.mean(frame$p_prev, frame$N), 4), "\n")

write.csv(frame, file.path(sim_dir, "virtual_population.csv"), row.names = FALSE)
write.csv(ext, file.path(sim_dir, "external_info.csv"), row.names = FALSE)

# Part B: Bayesian borrowing of external control information.
#
# The external information is the weighted prevalence among respondents aged
# 50 or older. The control event probability gets a robust mixture prior:
# weight w0 on Beta(a1, b1) with mean p_ext and prior sample size ess_ext,
# and weight 1 - w0 on Beta(1, 1). A trial is a success when
# P(treatment event rate < control event rate | data) exceeds thr.
# Drift multiplies the control odds relative to the external rate.

if (!exists("sim_dir")) source("config.R")
if (!exists("called_from_run_all")) set.seed(seed)
source("R/functions_borrowing.R")

ext = read.csv(file.path(sim_dir, "external_info.csv"))
p_ext = ext$p_ext_weighted
a1 = p_ext * ess_ext
b1 = (1 - p_ext) * ess_ext
cat("\nPart B:", n_sim_B, "simulated trials per scenario\n")
cat("External control rate:", round(p_ext, 4), "| prior sample size:", ess_ext,
    "| Kish effective sample size:", round(ext$n_eff_kish), "\n")

designs = data.frame(
  design = c("Frequentist 1:1", "Bayesian 1:1, vague priors",
             "Bayesian 1:1, CCHS-informed control prior",
             "Bayesian 2:1, vague priors",
             "Bayesian 2:1, CCHS-informed control prior"),
  r_t = c(0.5, 0.5, 0.5, 2 / 3, 2 / 3),
  prior = c("freq", "vague", "ext", "vague", "ext"),
  stringsAsFactors = FALSE)

resB = list()
k = 0
for (drift in drifts) for (OR in c(1, OR_alt)) {
  pc = plogis(qlogis(p_ext) + log(drift))
  pt = plogis(qlogis(pc) + log(OR))
  for (alloc in unique(designs$r_t)) {
    nt = round(N_B * alloc)
    nc = N_B - nt
    yc = rbinom(n_sim_B, nc, pc)
    yt = rbinom(n_sim_B, nt, pt)
    for (j in which(designs$r_t == alloc)) {
      pr = designs$prior[j]
      if (pr == "freq") {
        ok = freq_success(yc, nc, yt, nt, alpha_one)
        ok[is.na(ok)] = FALSE
        out = c(success = mean(ok), bias_rd = NA, rmse_rd = NA, w_inf = NA)
      } else {
        po = post_success(yc, nc, yt, nt, if (pr == "ext") w0 else 0, a1, b1)
        err = po$rd - (pt - pc)
        out = c(success = mean(po$p_success > thr), bias_rd = mean(err),
                rmse_rd = sqrt(mean(err^2)),
                w_inf = if (pr == "ext") mean(po$wpost) else NA)
      }
      k = k + 1
      resB[[k]] = data.frame(design = designs$design[j], drift = drift,
                             true_OR = OR, control_rate = pc,
                             n_control = nc, n_treatment = nt, t(out))
    }
  }
  cat("  drift", drift, "OR", OR, "done\n")
}
sumB = do.call(rbind, resB)
names(sumB)[names(sumB) == "success"] = "p_success"
print(sumB, digits = 3)
write.csv(sumB, file.path(sim_dir, "partB_borrowing_summary.csv"),
          row.names = FALSE)
saveRDS(designs, file.path(sim_dir, "partB_designs.rds"))

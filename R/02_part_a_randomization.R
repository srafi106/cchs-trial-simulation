# Part A: comparison of four randomization schemes.
#
# Each simulated trial draws N_A participants from the virtual population,
# randomizes them 1:1 under each scheme, and analyzes the same simulated
# outcomes with an unadjusted and an adjusted logistic regression.
# Under the null the true odds ratio is 1; under the alternative it is OR_alt.
# One-sided Wald test at level alpha_one; benefit is a negative log odds ratio.

if (!exists("sim_dir")) source("config.R")
if (!exists("called_from_run_all")) set.seed(seed)
source("R/functions_randomization.R")

pop = read.csv(file.path(sim_dir, "virtual_population.csv"),
               stringsAsFactors = FALSE)
P = prepare_population(pop)
np = P$np
cat("Cells:", nrow(pop), "| provinces:", np,
    "| mean control event probability:",
    round(sum(P$cell_p * pop$p_prev), 4), "\n")

cat("\nPart A:", n_sim_A, "simulated trials of", N_A, "participants\n")
vars = c("max_smd", "abs_dp0", "abs_dn",
         "rej_un_null", "rej_un_alt", "rej_ad_null", "rej_ad_alt",
         "est_un_alt", "est_ad_alt")
A = lapply(schemes, function(s) {
  matrix(NA_real_, n_sim_A, length(vars), dimnames = list(NULL, vars))
})
names(A) = schemes
beta_alt = log(OR_alt)
zcrit = qnorm(alpha_one)

t0 = Sys.time()
for (r in seq_len(n_sim_A)) {
  p = draw_participants(P, N_A)
  B = make_B(p, np)
  Xad = make_Xad(p, np)
  Fm = make_F(p, np)
  u = runif(N_A)
  for (s in schemes) {
    arm = assign_arm(s, p, Fm, block_size, p_bc)
    y0 = as.integer(u < p$p0)
    y1 = as.integer(u < plogis(qlogis(p$p0) + beta_alt * arm))
    un0 = fit_arm(y0, arm)
    un1 = fit_arm(y1, arm)
    ad0 = fit_arm(y0, arm, Xad)
    ad1 = fit_arm(y1, arm, Xad)
    A[[s]][r, ] = c(
      max_smd(B, arm),
      100 * abs(mean(p$p0[arm == 1]) - mean(p$p0[arm == 0])),
      abs(sum(arm) - sum(1 - arm)),
      un0[1] / un0[2] < zcrit, un1[1] / un1[2] < zcrit,
      ad0[1] / ad0[2] < zcrit, ad1[1] / ad1[2] < zcrit,
      un1[1], ad1[1])
  }
  if (r %% 100 == 0) {
    cat("  trial", r, "of", n_sim_A, "-",
        round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n")
  }
}
saveRDS(A, file.path(sim_dir, "partA_replicates.rds"))

sumA = do.call(rbind, lapply(schemes, function(s) {
  m = A[[s]]
  data.frame(
    scheme = s,
    mean_max_SMD = mean(m[, "max_smd"]),
    p95_max_SMD = unname(quantile(m[, "max_smd"], 0.95)),
    mean_abs_prognostic_diff_pp = mean(m[, "abs_dp0"]),
    mean_abs_arm_size_diff = mean(m[, "abs_dn"]),
    type1_unadjusted = mean(m[, "rej_un_null"], na.rm = TRUE),
    type1_adjusted = mean(m[, "rej_ad_null"], na.rm = TRUE),
    power_unadjusted = mean(m[, "rej_un_alt"], na.rm = TRUE),
    power_adjusted = mean(m[, "rej_ad_alt"], na.rm = TRUE),
    bias_logOR_unadj = mean(m[, "est_un_alt"], na.rm = TRUE) - beta_alt,
    bias_logOR_adj = mean(m[, "est_ad_alt"], na.rm = TRUE) - beta_alt)
}))
sumA$mc_se_type1 = sqrt(0.025 * 0.975 / n_sim_A)
sumA$mc_se_power = sqrt(0.8 * 0.2 / n_sim_A)
print(round(sumA[, -1], 3))
write.csv(sumA, file.path(sim_dir, "partA_randomization_summary.csv"),
          row.names = FALSE)

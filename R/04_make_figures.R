# Figures from the saved simulation results.

if (!exists("sim_dir")) source("config.R")

# Figure 1: covariate balance by randomization scheme
A = readRDS(file.path(sim_dir, "partA_replicates.rds"))
png(file.path(fig_dir, "fig1_balance.png"), 1500, 1000, res = 200)
par(mar = c(10, 4.5, 2, 1))
boxplot(lapply(A, function(m) m[, "max_smd"]), las = 2,
        ylab = "Maximum standardized difference", col = "grey80",
        main = "Covariate balance by randomization scheme", cex.axis = 0.7)
dev.off()

# Figure 2: probability of success against control-rate drift
sumB = read.csv(file.path(sim_dir, "partB_borrowing_summary.csv"))
designs = readRDS(file.path(sim_dir, "partB_designs.rds"))
cols = c("black", "grey55", "firebrick", "steelblue", "darkorange")

png(file.path(fig_dir, "fig2_borrowing.png"), 2000, 1100, res = 200)
layout(matrix(c(1, 2, 3, 3), nrow = 2, byrow = TRUE), heights = c(4.2, 1))
for (OR in c(1, OR_alt)) {
  par(mar = c(4.5, 4.5, 3, 1))
  d = sumB[sumB$true_OR == OR, ]
  plot(NA, xlim = range(drifts), ylim = c(0, max(d$p_success, 0.05) * 1.05),
       xlab = "Odds multiplier: true control rate vs external",
       ylab = "Probability of success",
       main = if (OR == 1) "No treatment effect (type I error)" else
         paste0("True OR = ", OR, " (power)"))
  for (j in seq_len(nrow(designs))) {
    dd = d[d$design == designs$design[j], ]
    lines(dd$drift, dd$p_success, type = "b", col = cols[j], pch = 16, lwd = 2)
  }
  if (OR == 1) abline(h = alpha_one, lty = 2)
}
par(mar = c(0, 0, 0, 0))
plot.new()
legend("center", designs$design, col = cols, pch = 16, lwd = 2, ncol = 3,
       cex = 0.8, bty = "n")
dev.off()

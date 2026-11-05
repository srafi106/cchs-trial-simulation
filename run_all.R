# Runs the whole project from the project root:  Rscript run_all.R
# Set TRIAL_QUICK=1 for a short test run.
#
# The random seed is set once before Part A and not reset before Part B, so
# running this file reproduces the reported numbers. Running Part B on its own
# gives the same design comparisons with different Monte Carlo error.

source("config.R")

if (requireNamespace("brms", quietly = TRUE)) {
  source("R/00_fit_baseline_model.R")
} else {
  message("brms is not installed; the baseline risk will come from a ",
          "survey-weighted logistic regression")
}
source("R/01_build_virtual_population.R")

called_from_run_all = TRUE
set.seed(seed)
source("R/02_part_a_randomization.R")
source("R/03_part_b_borrowing.R")
source("R/04_make_figures.R")

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
cat("\nFinished. Results are in", normalizePath(out_dir), "\n")

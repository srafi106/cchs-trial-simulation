# Settings for the CCHS trial simulation project.
# Set the environment variable TRIAL_QUICK=1 for a short test run.

pumf_path = "data/pumf_cchs.csv"
out_dir = "output"
fit_dir = file.path(out_dir, "fits")
sim_dir = file.path(out_dir, "simulation")
fig_dir = file.path(out_dir, "figures")
fit_name = "m2_plus_income_agg"

quick = identical(Sys.getenv("TRIAL_QUICK"), "1")
seed = 20261101

# Part A: randomization schemes
N_A = 2400
n_sim_A = if (quick) 100 else 2000
block_size = 4
p_bc = 0.8

# Part B: borrowing of external control information
N_B = 2000
n_sim_B = if (quick) 1000 else 10000
ess_ext = 300
w0 = 0.5
thr = 0.975
drifts = c(0.75, 1, 1.25, 1.5)

# Shared by both parts
OR_alt = 0.65
alpha_one = 0.025

for (dd in c(out_dir, fit_dir, sim_dir, fig_dir)) {
  dir.create(dd, showWarnings = FALSE, recursive = TRUE)
}

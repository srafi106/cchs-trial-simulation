# cchs-trial-simulation

Simulation study of randomization schemes and Bayesian borrowing of external
control information in a hypothetical two-arm cardiovascular prevention trial.
The virtual participants are drawn from the 2022 Canadian Community Health
Survey (CCHS) public use microdata file (PUMF).

The CCHS is a cross-sectional survey. It supplies the covariate distribution
(province, age group, sex, household income quintile) and a baseline-risk
gradient for adults aged 50 or older. It has no treatment or incidence data, so
the trial outcome and the treatment effect are simulated, and the survey
prevalence of a cardiovascular condition is used as the control-arm event
probability.

## What is compared

**Part A, randomization schemes** (2,400 participants, 2,000 simulated trials)

- simple randomization
- permuted blocks of 4 within province
- permuted blocks of 4 within province x age group x sex
- Pocock-Simon minimization on province, age group, sex and income quintile

Each scheme is analyzed with an unadjusted and a covariate-adjusted logistic
regression. Outputs are covariate balance, type I error, power and bias.

**Part B, borrowing of external control information** (2,000 participants,
10,000 simulated trials per scenario)

Five designs: frequentist 1:1, Bayesian 1:1 and 2:1 with vague priors, and
Bayesian 1:1 and 2:1 with a robust mixture prior for the control event rate
built from the survey prevalence. The true control odds are multiplied by 0.75,
1, 1.25 or 1.5 relative to the external rate to create prior-data conflict.
Outputs are type I error, power, bias and RMSE of the risk difference, and the
posterior weight on the informative prior component.

## Data

The data are not included in this repository. Download the CCHS 2022 PUMF from
Statistics Canada (Open Licence applies) and save the main file as
`data/pumf_cchs.csv`. The code uses these columns:

| column | content |
|---|---|
| `DHHGAGE` | age group |
| `DHH_SEX` | sex |
| `INCDGRCA` | national household income quintile |
| `GEODGHR4` | health region |
| `GEOGPRV` | province |
| `CCCDGCAR` | heart disease, heart attack or effects of stroke (ages 35+) |
| `WTS_M` | sampling weight |

The bootstrap weights file is not needed for the simulations.

## Repository layout

```
config.R                      settings and paths
run_all.R                     runs everything in order
R/
  functions_data.R            reading and recoding the PUMF
  functions_randomization.R   participants, randomization schemes, analysis
  functions_borrowing.R       Bayesian analysis with a mixture prior
  00_fit_baseline_model.R     multilevel model for baseline risk (needs brms)
  01_build_virtual_population.R
  02_part_a_randomization.R
  03_part_b_borrowing.R
  04_make_figures.R
data/                         put pumf_cchs.csv here (not tracked)
output/                       created when the code runs (not tracked)
```

## Requirements

R 4.3 or later (developed under 4.5.1). Packages: `data.table`, `dplyr`.
Step 00 also needs `tidyr`, `brms` and `rstan`. The simulations themselves use
base R only.

If `brms` is not installed, step 01 gets the baseline risk from a
survey-weighted logistic regression instead of the Bayesian model. The results
will differ slightly.

## Running

From the project root:

```
Rscript run_all.R
```

A short test run (100 trials in Part A, 1,000 in Part B) takes about a minute:

```
TRIAL_QUICK=1 Rscript run_all.R
```

On a laptop, Part A takes about 16 minutes at 2,000 trials and Part B a few
minutes. Step 00 takes a few minutes with 4 cores. Its fit is cached in
`output/fits`, so later runs reuse it.

## Output

- `output/simulation/virtual_population.csv`, `external_info.csv`
- `output/simulation/partA_randomization_summary.csv`, `partA_replicates.rds`
- `output/simulation/partB_borrowing_summary.csv`
- `output/figures/fig1_balance.png`, `fig2_borrowing.png`
- `output/sessionInfo.txt`

## Results

Full-size run (2,000 trials in Part A, 10,000 per scenario in Part B).

Part A, true odds ratio 0.65, one-sided test at 0.025:

| scheme | mean max standardized difference | type I error, adjusted | power, adjusted |
|---|---|---|---|
| Simple | 0.086 | 3.1% | 87.7% |
| Province blocks | 0.069 | 2.7% | 89.0% |
| Province x age x sex blocks | 0.062 | 2.7% | 87.6% |
| Minimization | 0.030 | 2.6% | 87.7% |

Minimization gave the best balance. Type I error and power were similar across
schemes, and covariate adjustment added about 1.5 percentage points of power.

Part B, success probability:

| odds multiplier | true OR | Frequentist 1:1 | Bayesian 2:1, vague | Bayesian 2:1, informed prior |
|---|---|---|---|---|
| 0.75 | 1 (type I error) | 2.7% | 2.5% | 6.1% |
| 1.00 | 1 (type I error) | 2.8% | 2.5% | 1.4% |
| 1.00 | 0.65 (power) | 83.9% | 80.2% | 88.2% |
| 1.50 | 0.65 (power) | 92.8% | 90.3% | 78.2% |

With an accurate external rate, the 2:1 design with borrowing used 667 control
participants and had higher power than the standard 1:1 design with 1,000. When
the true control rate was below the external rate, type I error was inflated,
and when it was above, power fell.

## Notes

- The treatment effect, event rate and sample sizes are assumptions. The
  results apply to the settings in `config.R`.
- Participants are sampled from covariate cells, so correlation among
  covariates within a cell is not represented.
- The success threshold in Part B is fixed at 0.975 and is not calibrated
  across the drift scenarios.
- To reproduce the numbers above exactly, use `run_all.R`. It sets the seed
  once before Part A and does not reset it before Part B.

## License

MIT, see `LICENSE`.

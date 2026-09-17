# Identifying Item Bias Without Conditioning: A Difference-in-Differences Approach

R code accompanying the manuscript **“Identifying Item Bias Without Conditioning: A Difference-in-Differences Approach.”**

The proposed difference-in-differences (DID) approach identifies average item bias for the focal group without conditioning on a matching variable, under the structural assumptions, anchor validity, and item response function (IRF) equivalence. The sensitivity analysis bounds identification error when IRF equivalence is relaxed.

All analyses are organized in one script: [`reproduce_did_study.R`](reproduce_did_study.R).

## Contents

- **Study 1:** DID and conventional uniform logistic-regression DIF under IRF equivalence, varying matching-variable validity and reliability.
- **Study 2:** Unadjusted and sensitivity-adjusted confidence intervals under violations of IRF equivalence, using the distribution-free and normal-distribution bounds.
- **Empirical illustration:** The Korean sample from the 2023 Programme for the International Assessment of Adult Competencies (PIAAC).

The repository contains the R script, this README, and `.gitignore`. Output folders are created when the script runs. Raw PIAAC data and generated results are not included.

## Requirements

The simulations, simulation summaries, and base-R numerical checks require **R 3.6.0 or later**. The empirical illustration additionally requires **dplyr**, **sandwich**, and the Korean PIAAC public-use CSV.

Install the empirical-analysis packages once in R:

```r
install.packages(c("dplyr", "sandwich"))
```

The supplied empirical results were obtained with R 4.5.0, dplyr 1.1.4, and sandwich 3.1-1. The script records the R version, package versions, and run settings in its output folders. Package versions are recorded, not automatically installed or pinned.

## Run the code

Download this repository and set your working directory to the folder containing `reproduce_did_study.R`. Relative input and output paths are resolved from that directory.

### From a terminal

First run the numerical checks and a short simulation:

```sh
Rscript reproduce_did_study.R checks
Rscript reproduce_did_study.R smoke
```

Then select the required analysis:

| Command | Tasks |
| --- | --- |
| `Rscript reproduce_did_study.R all` | Checks, both simulation studies and their outputs, and the PIAAC illustration |
| `Rscript reproduce_did_study.R simulations` | Checks and both simulation studies, including tables and figures |
| `Rscript reproduce_did_study.R empirical` | PIAAC illustration only |
| `Rscript reproduce_did_study.R outputs` | Regenerate simulation tables and figures from saved `simulation.rds` files |
| `Rscript reproduce_did_study.R checks` | Numerical implementation checks only |
| `Rscript reproduce_did_study.R smoke` | Checks and both simulations with 20 replications per condition |

Running the script without an argument uses the `RUN_MODE` set at the top of the file, initially `"all"`. A command-line argument overrides that setting.

**Smoke mode is an execution check, not a reproduction of the manuscript's numerical results.** It writes to `results/smoke/`, separately from the full-run results. Full reproduction uses **5,000 replications per condition** and can take substantial time; runtime depends on the computer.

### From RStudio

1. Open `reproduce_did_study.R`.
2. Set the working directory to the script's folder, for example with **Session > Set Working Directory > To Source File Location**.
3. Set `RUN_MODE` near the top of the script to the desired mode, such as `"checks"`, `"smoke"`, `"simulations"`, or `"empirical"`.
4. Click **Source** to run the entire file.

For individual tasks, set `RUN_MODE <- "custom"` and edit the six `RUN_*` switches. For example, to run Study 2 alone:

```r
RUN_MODE <- "custom"
RUN_STUDY1_SIMULATION <- FALSE
RUN_STUDY1_OUTPUTS    <- FALSE
RUN_STUDY2_SIMULATION <- TRUE
RUN_STUDY2_OUTPUTS    <- TRUE
RUN_CHECKS           <- TRUE
RUN_EMPIRICAL        <- FALSE
```

Edit these settings **inside the script before sourcing it**. Changing variables only in the console will not override assignments in the file.

## PIAAC data

Download and extract the Korean public-use CSV, `prgkorp2.csv`, from the [OECD PIAAC 2nd Cycle Database](https://www.oecd.org/en/data/datasets/PIAAC-2nd-Cycle-Database.html). Place it beside the R script, or set `EMPIRICAL_DATA_FILE` to its local path.

```r
EMPIRICAL_DATA_FILE <- "prgkorp2.csv"
```

The reader expects the OECD semicolon-delimited CSV with `.` as the decimal separator. The configured missing-value codes are `.`, `.n`, and `.v`. CSV files do not retain the value labels available in labeled statistical-data formats.

The analysis uses:

| Purpose | Variable / selection |
| --- | --- |
| Test item | `E320004S` |
| Anchor item | `E320003S` |
| Age grouping | `AGEG10LFS`: code 2 for the reference group, code 4 for the focal group |
| Reference group | Adults aged 25–34 |
| Focal group | Adults aged 45–54 |
| Literacy plausible values | `PVLIT1` through `PVLIT10` |
| Analytic sample | Observed responses to both items; 285 reference and 251 focal respondents, N = 536 |

Observed item scores are checked to be 0 or 1. The script also checks the expected sample sizes, required columns, and empirical logistic-regression fits. If a check fails, verify the downloaded data and coding rather than removing the check to force execution.

The conventional DIF analysis fits `Y_T ~ M + G` separately for each plausible value and combines the group coefficients and model-based variances using multiple-imputation combining rules. DID is estimated by `lm((Y_T - Y_A) ~ G)` with an HC3 standard error. Both analyses use large-sample normal inference.

All empirical analyses are **unweighted**. The illustration does not provide a population-representative estimate for Korean adults. The sensitivity inputs are treated as fixed; their estimation uncertainty is not incorporated into the intervals.

Modes `all` and `empirical` require the local CSV and the two packages. Missing prerequisites stop execution; choose `simulations`, `checks`, or `smoke` to run without the PIAAC data.

## Simulation settings

| Setting | Study 1 | Study 2 |
| --- | --- | --- |
| Sample size, N | 500, 1,000, 2,000 | 500, 1,000, 2,000 |
| Item-bias estimand, tau | 0, −.05, −.10 | 0, −.05, −.10 |
| Matching-variable invalidity, delta | 0, .25, .50 | Not used |
| Matching-variable reliability, rho_M | 1.00, .80, .60 | Not used |
| Maximum absolute IRF difference, eta | 0 | 0, .05, .10, .15 |
| Conditions | 81 | 36 |
| Replications per condition | 5,000 | 5,000 |
| Base seed | 2026 | 2027 |

Both studies use reference and focal latent-trait distributions N(0, 1) and N(−.5, 1), respectively. Item responses are generated independently conditional on group membership and the latent trait. The focal-group logit coefficient is calibrated to the target value of tau.

Study 2 keeps both test-item IRFs fixed across eta and varies only the anchor IRF:

```r
P_T_reference <- plogis(1.5 * theta)
P_T_focal     <- plogis(1.5 * theta + gamma)
P_A           <- eta + (1 - 2 * eta) * plogis(1.5 * theta)
```

The bounds are:

```r
b_df     <- 2 * eta
b_normal <- eta * (4 * pnorm(abs(mu) / 2) - 2)
```

For an unadjusted confidence interval `[L, U]`, a bound `b` gives the sensitivity-adjusted interval `[L - b, U + b]`. The intervals are not clipped to [−1, 1]. Study 2 uses the data-generating values of eta and mu as fixed sensitivity inputs. Its main interval summaries report coverage of **tau**, not coverage of the population DID.

## Output files

Full-run results are written under `results/` by default.

| Manuscript output / purpose | File |
| --- | --- |
| Study 1 DID performance, Table S1 | `results/study1/table_study1_did.csv` |
| Study 1 logistic-regression DIF, Table S2 | `results/study1/table_study1_logistic.csv` |
| Study 1, Figure 4 | `results/study1/figure4.pdf` |
| Identification errors and bounds, Table S3 | `results/study2/tableS3_absolute_id_error_bounds.csv` |
| Study 2 intervals, tau = 0, Table S4 | `results/study2/tableS4_intervals.csv` |
| Study 2 intervals, tau = −.05, Table S5 | `results/study2/tableS5_intervals.csv` |
| Study 2 intervals, tau = −.10, Table S6 | `results/study2/tableS6_intervals.csv` |
| Study 2, Figure 5 | `results/study2/figure5.pdf` |
| Empirical DIF and DID, Table 1 | `results/empirical/table1_empirical_results.csv` |
| Empirical sensitivity intervals | `results/empirical/illustration_sensitivity.csv` |

Each simulation folder also contains `simulation.rds`, design and calibration CSV files, detailed summaries, figure data and captions, `metadata.txt`, and `sessionInfo.txt`. The empirical folder includes plausible-value results and run information. CSV files retain unrounded values; manuscript tables may round and rearrange columns.

Study 1 DID results are pooled from raw replications across matching-variable conditions. Logistic-regression summaries report invalid fits and warnings; rejection rates use valid fits as their denominators.

`study2_population_coverage_diagnostic.csv` is an additional check of coverage of the **population DID**. It is distinct from the coverage-of-tau results in Tables S4–S6.

### Expected empirical results

These are reference values from the supplied empirical analysis, rounded as in the manuscript:

| Method | Estimate | SE | 95% CI | p |
| --- | ---: | ---: | --- | ---: |
| Logistic-regression DIF | −.638 | .229 | [−1.086, −.190] | .005 |
| DID | −.119 | .042 | [−.201, −.036] | .005 |

With eta = .033 and an absolute standardized group difference of approximately .718421:

| Sensitivity bound | Bound, rounded | Adjusted 95% CI |
| --- | ---: | --- |
| Distribution-free | .066 | [−.267, .030] |
| Normal-distribution | .019 | [−.220, −.018] |

## Reproducibility notes

- The script sets the random-number generator and uses `base_seed + 1009 * condition_id` for each simulation condition. Keep the seeds, condition order, replication count, and random-draw order unchanged to reproduce the supplied design.
- Output-only mode uses the settings stored with the saved simulations. Changing `N_REP` or `ALPHA` does not change existing saved replications.
- Study 2 checks its saved model identifier before generating outputs. Results from an earlier difficulty/discrimination design cannot be substituted for the current fixed-test design.
- Rerunning a task replaces that task's files in the selected output folder. Use a different `OUTPUT_ROOT` for exploratory analyses.
- `.gitignore` excludes local raw data and generated results from ordinary Git staging. Keep the PIAAC CSV local when uploading repository files through a browser as well.

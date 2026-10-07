# Identifying Item Bias Without Conditioning: A Difference-in-Differences Approach

This repository contains one self-contained R script for the manuscript's two simulation studies and PIAAC illustration: **`reproduce_did_study.R`**. It produces the numerical results, Figures 4 and 5, and supplementary Tables S1–S4 using base R.

The proposed difference-in-differences (DID) approach identifies average item bias under a valid anchor and parallel response differences. Its sensitivity analysis bounds identification error under item response function (IRF) nonequivalence.

## Requirements

- R 3.6.0 or later; no contributed R packages are required.
- Arial and a PDF graphics device: native Quartz on macOS or Cairo on other systems. Set `MAKE_FIGURES <- FALSE` to export results and tables without figures.
- The Korean PIAAC public-use CSV, obtained separately, for the empirical illustration only.

Run commands from the repository folder. In RStudio, set that folder as the working directory, edit the settings at the top of the script, and source the file.

## Update figures and tables without rerunning simulations

```sh
Rscript reproduce_did_study.R presentation
```

**`presentation` is the default mode.** It reads four existing CSV summaries and updates Figures 4 and 5 and Tables S1–S4. It does not generate samples, fit models, or recompute the summary statistics. It leaves the input CSV files, raw simulation files, calibration file, and original run records unchanged.

The summaries must be in these locations, relative to `OUTPUT_DIR` (default: `results/`):

| Required file | Used for |
|---|---|
| `study1/did_summary.csv` | Figure 4 and Table S1 |
| `study1/logistic_summary.csv` | Figure 4 and Table S2 |
| `study2/main_summary.csv` | Figure 5 and Table S3 |
| `study2/eta_summary.csv` | Table S4 |

Copy existing summaries to these paths, or change `OUTPUT_DIR` to the folder containing the `study1/` and `study2/` folders. Use summaries from the same simulation run and retain their original precision. The Study 2 summaries must contain `coverage`, `mean_length`, and `rejection`; the script reports missing columns rather than estimating missing results.

If the summaries are missing, this mode stops with instructions. It never starts simulations automatically. The repository contains code and documentation; saved results must be generated or supplied locally.

Presentation runs write `presentation_settings.txt` and `presentation_sessionInfo.txt`. The settings record includes checksums of the four input summaries. Changes to the data-generating process or sensitivity assumptions require analysis from raw results or a new simulation run; changing settings in presentation mode does not alter existing estimates.

## Other run modes

| Mode | Action |
|---|---|
| `smoke` | Run both studies with 20 replications per condition in `results_smoke/` |
| `simulations` | Run both studies with 5,000 replications per condition |
| `study1` | Run Study 1 only, or reuse `STUDY1_SAVED` if specified |
| `study2` | Run Study 2 only, or reuse `STUDY2_SAVED` if specified |
| `outputs` | Recompute summaries, figures, and tables from saved raw RDS files |
| `empirical` | Run the PIAAC illustration only |
| `all` | Run both studies and the PIAAC illustration |

For a first run, check execution before starting the full simulations:

```sh
Rscript reproduce_did_study.R smoke
Rscript reproduce_did_study.R simulations
```

The full Study 1 run fits 405,000 logistic regressions and can take substantial time. Progress is printed by condition. Smoke results check execution and output generation; they are not manuscript results.

To recompute results from saved replication-level estimates:

```sh
Rscript reproduce_did_study.R outputs
```

This reads `results/study1/simulation.rds` and `results/study2/simulation.rds` by default. It generates no new samples and does not overwrite the input RDS files. Unlike `presentation`, it recalculates the summaries.

For compatible files from earlier versions, including `simulation_revised.rds`, edit:

```r
STUDY1_SAVED <- "path/to/study1/simulation_revised.rds"
STUDY2_SAVED <- "path/to/study2/simulation_revised.rds"
```

Then use `outputs`, or use `study1` or `study2` to reanalyze only that saved study. Saved files must contain compatible designs and replication-level estimates. The loader checks both before reuse.

Fresh simulation runs replace the corresponding generated files under `OUTPUT_DIR`. Change this setting to retain an earlier run. `smoke` always uses a separate directory and ignores saved-input paths. `presentation` also ignores saved-input paths because it reads CSV summaries only.

## Study 1

Study 1 compares DID with uniform logistic-regression DIF analysis as matching-variable validity and reliability vary. The identifying conditions for DID remain satisfied.

The fully crossed design contains 81 conditions, each with 5,000 replications in a full run:

| Factor | Values |
|---|---|
| Sample size, `N` | 500, 1,000, 2,000 |
| True average item-bias effect, `tau` | 0, -0.05, -0.10 |
| Direct group effect on the matching variable, `delta` | 0, 0.25, 0.50 |
| Matching-variable reliability, `rho` | 1, 0.80, 0.60 |

Group membership has probability 0.5. The latent distributions are normal with means 0 and -0.5 in the reference and focal groups, respectively, and variance 1 in both groups. A calibrated logit shift gives the specified probability-scale effect `tau`.

DID inference uses HC3 standard errors. Logistic-regression DIF analysis uses a two-sided Wald test for the group coefficient. Rejection rates at the 0.05 level represent Type I error when `tau = 0` and power otherwise. Failed logistic fits are recorded and excluded from rejection-rate denominators; their datasets are not replaced.

DID performance is pooled across the nine matching-variable conditions because these conditions affect neither the DID estimator nor the distribution of its inputs. A full run therefore pools 45,000 replications for each sample-size/effect-size combination. Logistic results retain all 81 conditions.

- **Figure 4:** rejection rates under matching-variable invalidity (A) and unreliability (B).
- **Table S1:** pooled DID bias, RMSE, SE/SD, coverage, and rejection rates.
- **Table S2:** logistic-regression DIF rejection rates for the full crossed design.

The Study 1 data-generating process, condition order, random-number draws, seed rule, and Figure 4 design are retained from the previous repository version.

## Study 2

Study 2 compares unadjusted and sensitivity-adjusted DID intervals using coverage, mean interval length, and rejection rates. Coverage is the proportion of intervals containing the true average item-bias effect `tau`. Rejection is the proportion excluding zero, representing Type I error when `tau = 0` and power otherwise.

The fully crossed design contains 36 conditions, each with 5,000 replications in a full run:

| Factor | Values |
|---|---|
| Sample size, `N` | 500, 1,000, 2,000 |
| True average item-bias effect, `tau` | 0, -0.05, -0.10 |
| Actual IRF discrepancy, `eta_true` (η₀) | 0, 0.05, 0.10, 0.15 |

The reference-condition test IRF is `p0(theta) = plogis(1.5 * theta)`. The anchor IRF is:

```r
eta_true + (1 - 2 * eta_true) * p0(theta)
```

Thus, `eta_true` is the supremum of the absolute difference between the reference-condition test IRF and the anchor IRF. The anchor remains valid. The test IRFs and calibrated logit shifts follow Study 1. In both studies, item responses are generated independently conditional on group and latent trait.

Given the unadjusted nominal 95% DID confidence interval `[L, U]`, the two sensitivity-adjusted intervals are:

- **Distribution-free:** `[L - 2 * eta, U + 2 * eta]`.
- **Normal-distribution:** `[L - b, U + b]`, where `b = eta * (4 * pnorm(abs(mu) / 2) - 2)` under equal-variance normal latent distributions.

The specified sensitivity parameter `eta` is distinct from the actual discrepancy `eta_true`. The main analysis sets `eta = eta_true` and `mu = -0.5`. Sensitivity inputs are fixed, and intervals are not truncated to the possible range of `tau`.

- **Figure 5:** coverage (A), mean interval length (B), and rejection rates (C) at `tau = -0.05`. The dotted line in Panel A marks 0.95 coverage. No Monte Carlo error bars are plotted.
- **Table S3:** all three measures at `tau = 0` and `tau = -0.10`, with `eta = eta_true`.
- **Table S4:** all three measures across `eta / eta_true = 0, 0.10, 0.25, 0.50, 0.75, 1, 1.25`, at `tau = -0.05` and `eta_true = 0.15`.

All methods and sensitivity specifications use the same simulated samples within each condition. Under IRF nonequivalence, the population DID can differ from `tau`; a high rejection rate alone therefore does not establish good performance. Coverage and interval length should be considered alongside rejection rates.

## Generated files

Paths below are relative to `OUTPUT_DIR` (default: `results/`).

| Folder | Files |
|---|---|
| `study1/` | `simulation.rds`, `design.csv`, `did_summary.csv`, `logistic_summary.csv`, `figure4.pdf`, `table_S1_did.tex`, `table_S2_logistic.tex` |
| `study2/` | `simulation.rds`, `design.csv`, `main_summary.csv`, `eta_summary.csv`, `did_diagnostics.csv`, `length_comparison.csv`, `eta_thresholds.csv`, `figure5.pdf`, `table_S3_other_effects.tex`, `table_S4_eta_choice.tex` |
| `empirical/` | `empirical_results.csv`, `sensitivity.csv`, `pv_results.csv`, `analysis_notes.txt` |
| Root, analysis runs | `gamma_calibration.csv`, `settings.txt`, `sessionInfo.txt` |
| Root, presentation runs | `presentation_settings.txt`, `presentation_sessionInfo.txt` |

Only files relevant to the selected mode are written. Presentation mode replaces the two figures (when enabled), four LaTeX tables, and its own run records. RDS files are written only when new samples are generated.

The raw RDS files contain the design and replication-level estimates. The CSV summaries include valid/failed counts and Monte Carlo standard errors. Study 2 diagnostics also report coverage of the population DID, relative interval lengths, and the design-specific sensitivity ratios at which each bound equals the absolute identification error.

Tables use the manuscript's CUP LaTeX macros (`\TBL`, `\TCH`, `fntable`, and `\botrule`). LaTeX assigns table numbers; the exported tables do not reset counters. Tables S3 and S4 each contain Coverage, Length, and Rejection columns for all three intervals.

Figures are vector PDFs with shared Arial typography, panel labels, sample-size headings, and line styling. Figure 5 has three panel rows and is taller than Figure 4. Font sizes are scaled for insertion at a width of 144 mm; set `FINAL_FIGURE_WIDTH_MM` to the actual insertion width and insert both figures at that width. The Figure 5 coverage axis includes space above 1 to keep triangle markers visible. Check the exported PDFs for local font substitution.

## PIAAC illustration

Obtain the Korean PIAAC public-use CSV separately and set `PIAAC_FILE` to its path. The default is `prgkorp2.csv`. The script expects a semicolon-delimited file with these columns:

| Variable | Use |
|---|---|
| `AGEG10LFS` | Age-group code: 2 for ages 25–34; 4 for ages 45–54 |
| `E320004S` | Binary test-item response |
| `E320003S` | Binary anchor-item response |
| `PVLIT1`–`PVLIT10` | Literacy plausible values |

```sh
Rscript reproduce_did_study.R empirical
```

Respondents must have observed binary responses on both items and belong to one of the two age groups. The manuscript sample contains 536 respondents: 285 reference-group and 251 focal-group respondents. The script warns if the supplied data yield different counts.

Logistic-regression DIF analysis is fitted separately for each plausible value; estimates and variances are combined using multiple-imputation rules. DID uses the item-response difference and an HC3 standard error. The illustration is unweighted and does not incorporate the complex survey design.

For the normal-distribution sensitivity scenario, `mu` is the mean signed standardized group difference across the ten plausible values. The default sensitivity parameter is `eta = 0.033`. Both inputs are treated as fixed; the adjusted intervals do not incorporate calibration uncertainty. `analysis_notes.txt` records the input checksum and analytic sample counts.

## Reproducibility

Simulation runs use `RNGkind("Mersenne-Twister", "Inversion", "Rejection")` and condition-specific seeds:

- Study 1: `2026 + 1009 * condition_id`.
- Study 2: `2027 + 1009 * condition_id`.

Random-number generation is serial, and each condition starts from its own seed. Small numerical differences across R versions and platforms may affect fitting or values near decision thresholds. Reusing saved results avoids generating new samples.

Analysis runs check effect calibration, the HC3 variance formula, the population DID identity, interval nesting, and the increase in mean interval length implied by each bound. Saved-data loading checks design compatibility and replication identifiers. Presentation mode checks the required summary columns and condition keys; it cannot verify the original data-generating process from summaries alone.

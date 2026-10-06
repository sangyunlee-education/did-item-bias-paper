# Identifying Item Bias Without Conditioning: A Difference-in-Differences Approach

This repository contains the R code for the manuscript's two simulation studies and PIAAC illustration. The proposed difference-in-differences (DID) approach identifies the average item-bias effect under the stated identifying conditions, including parallel response differences. Its sensitivity analysis bounds identification error when IRF equivalence does not hold.

All calculations, figures, and supplementary tables are produced by **`reproduce_did_study.R`**. No additional R scripts or contributed R packages are required.

## Requirements

- R 3.6.0 or later.
- Arial for the manuscript's figure typography.
- Native Quartz PDF support on macOS, or Cairo PDF support on other systems. Set `MAKE_FIGURES <- FALSE` to run the analyses and export tables without a PDF device.
- The Korean PIAAC public-use CSV for the empirical illustration only.

Run commands from the repository folder. Alternatively, open the script in RStudio, set its folder as the working directory, edit the settings at the top, and source the file.

## Running the script

Start with a small run that exercises both simulation studies and their outputs:

```sh
Rscript reproduce_did_study.R smoke
```

This uses 20 replications per condition and writes to `results_smoke/`. These results are for checking execution, not for the manuscript.

Run the full simulation studies with 5,000 replications per condition:

```sh
Rscript reproduce_did_study.R simulations
```

The script prints progress by condition. Study 1 includes 405,000 logistic regressions, so the full run can take substantial time. The studies can also be run separately:

```sh
Rscript reproduce_did_study.R study1
Rscript reproduce_did_study.R study2
```

To regenerate summaries, figures, and tables from saved raw simulations:

```sh
Rscript reproduce_did_study.R outputs
```

By default, this reads `results/study1/simulation.rds` and `results/study2/simulation.rds`. No new samples are generated, and the input RDS files are not overwritten. Saved simulation files are generated locally; the script does not download them.

For files from earlier versions, including `simulation_revised.rds`, set their paths near the top of the script:

```r
STUDY1_SAVED <- "path/to/study1/simulation_revised.rds"
STUDY2_SAVED <- "path/to/study2/simulation_revised.rds"
```

Then run `outputs`. To reanalyze only one saved study, set its path and run `study1` or `study2`. The loader checks the saved design and accepts the original repository schema and the revised raw-result schema. CSV summaries alone are insufficient for this reanalysis.

Fresh simulation runs replace the corresponding generated files under `OUTPUT_DIR`. To retain an earlier run, change `OUTPUT_DIR` before starting. The smoke run always uses a separate directory.

## Study 1

Study 1 compares DID with uniform logistic-regression DIF analysis as matching-variable validity and reliability vary. The identifying conditions for DID remain satisfied.

The fully crossed design contains 81 conditions:

| Factor | Values |
|---|---|
| Sample size, `N` | 500, 1,000, 2,000 |
| True average item-bias effect, `tau` | 0, −.05, −.10 |
| Direct group effect on the matching variable, `delta` | 0, .25, .50 |
| Matching-variable reliability, `rho` | 1, .80, .60 |

Group membership has probability .5. The latent distributions are standard normal in the reference group and normal with mean −.5 and variance 1 in the focal group. A calibrated logit shift gives the specified probability-scale effect `tau`.

DID inference uses HC3 standard errors. Logistic-regression DIF analysis uses a two-sided Wald test for the group coefficient. At the .05 level, rejection rates represent Type I error when `tau = 0` and power otherwise. Failed logistic fits are recorded and excluded from rejection-rate denominators; their datasets are not replaced.

DID performance is summarized over the nine matching-variable conditions because these conditions affect neither the DID estimator nor the distribution of its inputs. The full run therefore pools 45,000 replications per sample-size/effect-size combination. Logistic results retain all 81 conditions.

The Study 1 data-generating process, condition order, random-number draws, and seed rule are retained from the original script.

## Study 2

Study 2 compares unadjusted and sensitivity-adjusted DID intervals using coverage of `tau` and mean interval length. Coverage refers to the true average item-bias effect, not the population DID, which can differ from `tau` under IRF nonequivalence.

The fully crossed design contains 36 conditions:

| Factor | Values |
|---|---|
| Sample size, `N` | 500, 1,000, 2,000 |
| True average item-bias effect, `tau` | 0, −.05, −.10 |
| Actual extent of IRF nonequivalence, `eta_true` (η₀) | 0, .05, .10, .15 |

The reference-condition test IRF is `p0(theta) = plogis(1.5 * theta)`. The anchor IRF is

```r
eta_true + (1 - 2 * eta_true) * p0(theta)
```

Thus, `eta_true` is the supremum of the absolute difference between these IRFs. The focal-group test IRF retains the calibrated logit shift used in Study 1. Item responses are conditionally independent given group and latent ability.

Let `[L, U]` be the unadjusted nominal 95% DID confidence interval. The sensitivity-adjusted intervals are:

- **Distribution-free:** `[L - 2 * eta, U + 2 * eta]`.
- **Normal-distribution:** `[L - b, U + b]`, where `b = eta * (4 * pnorm(abs(mu) / 2) - 2)` under the specified equal-variance normal latent distributions.

Here, `eta` is the analyst's sensitivity parameter, distinct from the actual discrepancy `eta_true`. The main analysis sets `eta = eta_true` and `mu = -.5`. Both inputs are treated as fixed. Intervals are not truncated to the possible range of `tau`.

Figure 5 presents coverage (A) and mean interval length (B) for `tau = -.05`. Its layout, typography, panel labels, and dimensions follow Figure 4. The dotted line marks .95 coverage. Monte Carlo error bars are not plotted.

The supplementary tables report:

- **Table S3:** results for `tau = 0` and `tau = -.10`, with `eta = eta_true`.
- **Table S4:** results for `eta / eta_true = 0, .10, .25, .50, .75, 1, 1.25`, at `tau = -.05` and `eta_true = .15`.

All sensitivity specifications use the same simulated samples within a condition. Specifying `eta < eta_true` can still produce conservative coverage in this design because the bounds can exceed the actual absolute identification error. This does not establish validity for arbitrary misspecification of `eta`.

## Generated files

The default output directory is `results/`.

| Location | File | Contents |
|---|---|---|
| `study1/` | `simulation.rds` | Design and replication-level DID/logistic results |
| `study1/` | `design.csv` | 81 simulation conditions |
| `study1/` | `did_summary.csv` | Pooled DID bias, RMSE, SE/SD, coverage, and rejection rates |
| `study1/` | `logistic_summary.csv` | Condition-specific rejection rates and fit counts |
| `study1/` | `figure4.pdf` | Matching-variable invalidity and unreliability |
| `study1/` | `table_S1_did.tex` | Supplementary DID performance table |
| `study1/` | `table_S2_logistic.tex` | Supplementary logistic-regression DIF table |
| `study2/` | `simulation.rds` | Design and replication-level DID results |
| `study2/` | `design.csv` | 36 conditions, population DID, identification error, and bounds |
| `study2/` | `main_summary.csv` | Interval performance with `eta = eta_true` |
| `study2/` | `eta_summary.csv` | Interval performance across sensitivity specifications |
| `study2/` | `did_diagnostics.csv` | DID means and unadjusted coverage of population DID |
| `study2/` | `length_comparison.csv` | Relative mean lengths of the two adjusted intervals |
| `study2/` | `eta_thresholds.csv` | Design-specific ratios at which bounds equal absolute identification error |
| `study2/` | `figure5.pdf` | Coverage and mean interval length at `tau = -.05` |
| `study2/` | `table_S3_other_effects.tex` | Results at additional effect sizes |
| `study2/` | `table_S4_eta_choice.tex` | Results across specified values of `eta` |
| Root | `gamma_calibration.csv` | Calibrated logit shifts and achieved effects |
| Root | `sessionInfo.txt` | R version and platform information |
| Root | `settings.txt` | Settings and saved-input paths for the latest run |

Tables use the manuscript's CUP LaTeX macros (`\TBL`, `\TCH`, `fntable`, and `\botrule`) and are intended for inclusion in that template. Table numbers follow the manuscript's counters. Supplementary Study 2 results are exported as tables, not additional figures.

Figures are vector PDFs. Font sizes are scaled for an inclusion width of 144 mm; change `FINAL_FIGURE_WIDTH_MM` if the manuscript uses another width. Insert Figures 4 and 5 at the same width. Font substitution depends on the local graphics device, so check that Arial is installed and inspect the exported PDFs before submission.

## PIAAC illustration

Obtain the Korean PIAAC public-use file separately and set `PIAAC_FILE` to its local path. The default filename is `prgkorp2.csv`. The script expects a semicolon-delimited file with these columns:

| Variable | Use |
|---|---|
| `AGEG10LFS` | Age-group code: 2 for the reference group; 4 for the focal group |
| `E320004S` | Binary test-item response |
| `E320003S` | Binary anchor-item response |
| `PVLIT1`–`PVLIT10` | Literacy plausible values |

Run the illustration with:

```sh
Rscript reproduce_did_study.R empirical
```

Or run both studies and the illustration with:

```sh
Rscript reproduce_did_study.R all
```

Respondents must have observed binary responses on both items and belong to one of the two age groups. The manuscript sample contains 536 respondents: 285 reference-group and 251 focal-group respondents. The script warns if the supplied data yield different counts.

Logistic-regression DIF analysis is fitted separately for each plausible value; estimates and variances are combined using multiple-imputation rules. DID uses the item-response difference and an HC3 standard error. Analyses are unweighted and model-based. The illustration does not use survey weights or replicate weights.

For the normal-distribution sensitivity scenario, `mu` is the mean signed standardized group difference across the ten plausible values. The default sensitivity parameter is `eta = .033`. These inputs are treated as fixed; the adjusted intervals do not incorporate uncertainty in their calibration.

Outputs under `results/empirical/` are `empirical_results.csv`, `sensitivity.csv`, `pv_results.csv`, and `analysis_notes.txt`. The last file records the input checksum and analytic sample counts.

## Reproducibility and checks

The script uses `RNGkind("Mersenne-Twister", "Inversion", "Rejection")` and condition-specific seeds:

- Study 1: `2026 + 1009 * condition_id`.
- Study 2: `2027 + 1009 * condition_id`.

Random-number generation is serial, and each condition starts from its own seed. Reusing saved raw simulations avoids differences in newly generated samples. Small numerical differences across R versions and platforms may still affect fitting or values near decision thresholds.

The script checks effect calibration, the HC3 variance formula, the population DID identity, nesting of adjusted intervals, and the exact increase in mean interval length implied by each bound. Saved-data loading checks design compatibility and replication identifiers. CSV summaries retain valid/failed counts and Monte Carlo standard errors. Reported table counts are taken from the actual results, including when saved files are reused.

Run `smoke` locally before a full run. A successful smoke run checks execution and output generation; its small number of replications does not validate the manuscript's numerical results.

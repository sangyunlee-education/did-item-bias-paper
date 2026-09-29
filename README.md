# Identifying Item Bias Without Conditioning: A Difference-in-Differences Approach

R code for the manuscript **“Identifying Item Bias Without Conditioning: A Difference-in-Differences Approach.”**

DID identifies average item bias for the focal population without conditioning on the latent trait, under the structural setup, anchor validity, and expected item difference equivalence. IRF equivalence is a stronger sufficient condition. The sensitivity analysis bounds identification error under bounded IRF nonequivalence.

All analyses are in [`reproduce_did_study.R`](reproduce_did_study.R): settings, shared calculations, simulations, figures/tables, checks, PIAAC, and execution.

## Requirements

- R 3.6.0 or later. Simulations and tables use base R.
- PDF figures use Arial, with native Quartz on macOS or Cairo support in R on other systems. Install Arial if it is unavailable; font substitution can change the appearance. `FIGURE_FONT` controls the family.
- The PIAAC illustration additionally requires `dplyr`, `sandwich`, and the local Korean public-use CSV:

```r
install.packages(c("dplyr", "sandwich"))
```

Raw data and generated results are not included. Each run records the R/package versions; packages are not automatically installed or pinned.

## Run

Set the working directory to the folder containing the script. In RStudio, use **Session > Set Working Directory > To Source File Location**, edit `RUN_MODE` near the top, and click **Source**.

```r
RUN_MODE <- "simulations"  # First run: 5,000 replications per condition
```

The default runs both simulations and generates Figures 4–5 and Tables S1–S2. It does not require PIAAC data. Full simulation can take substantial time.

For existing saved simulations, change only the mode:

```r
RUN_MODE <- "outputs"
```

This regenerates figures and appendix tables from `results/study1/simulation.rds` and `results/study2/simulation.rds`, without rerunning simulations. The script checks that the saved Study 2 results use the specified IRFs.

The same modes are available from a terminal:

| Command | Purpose |
| --- | --- |
| `Rscript reproduce_did_study.R simulations` | Both simulations, figures, and appendix tables |
| `Rscript reproduce_did_study.R outputs` | Figures and appendix tables from saved simulations |
| `Rscript reproduce_did_study.R empirical` | PIAAC illustration only |
| `Rscript reproduce_did_study.R all` | Simulations and PIAAC illustration |
| `Rscript reproduce_did_study.R checks` | Numerical checks only |
| `Rscript reproduce_did_study.R smoke` | Both simulations with 20 replications per condition |

Smoke results are written separately under `results/smoke/` and are not manuscript results. A terminal argument overrides `RUN_MODE`. Edit settings inside the script before sourcing it; assignments made only in the console are overwritten.

## Manuscript outputs

| Output | Path |
| --- | --- |
| Figure 4 | `results/study1/figure4.pdf` |
| Figure 4 data | `results/study1/figure4_data.csv` |
| Table S1: DID performance | `results/study1/table_study1_did.csv` and `.tex` |
| Table S2: logistic-regression DIF rejection rates | `results/study1/table_study1_logistic.csv` and `.tex` |
| Figure 5 | `results/study2/figure5.pdf` |
| Figure 5 data | `results/study2/figure5_data.csv` |
| Empirical Table 1 | `results/empirical/table1_empirical_results.csv` |
| Empirical sensitivity intervals | `results/empirical/illustration_sensitivity.csv` |

The `.tex` files use the manuscript's CUP `\TBL`/`fntable` format and `booktabs`; paste their contents into the appendix in S1–S2 order. Labels remain `tab:sim1_did` and `tab:sim1_logistic`. CSV values retain full precision; LaTeX values are rounded for presentation.

Study 1 DID results pool raw replications across the nine matching-variable conditions within each sample-size and item-bias combination: nine rows, each based on 45,000 replications in a full run. Empirical SD, RMSE, and SE/SD are calculated from the pooled replications, not by averaging condition-level ratios. Logistic-regression results retain the full crossed design and use valid fits as denominators; the CSV includes failure and warning counts. The `rho_M` field denotes the manuscript's reliability parameter rho.

Figure 5 omits eta = 0, where the three intervals coincide; these results remain in its CSV.

Simulation folders retain `simulation.rds`, `metadata.txt`, and `sessionInfo.txt` to support regeneration and reproducibility.

## Simulation design

| Setting | Study 1 | Study 2 |
| --- | --- | --- |
| N | 500, 1,000, 2,000 | 500, 1,000, 2,000 |
| tau | 0, −.05, −.10 | 0, −.05, −.10 |
| Matching-variable invalidity, delta | 0, .25, .50 | Not used |
| Matching-variable reliability, rho | 1.00, .80, .60 | Not used |
| IRF difference bound, eta | 0 | 0, .05, .10, .15 |
| Conditions | 81 | 36 |
| Replications per condition | 5,000 | 5,000 |
| Base seed | 2026 | 2027 |

Both studies use latent-trait distributions N(0, 1) and N(−.5, 1), with conditionally independent responses. The focal logit coefficient is calibrated to tau. Study 2 varies only the anchor IRF:

```r
P_T_reference <- plogis(1.5 * theta)
P_T_focal <- plogis(1.5 * theta + gamma)
P_A <- eta + (1 - 2 * eta) * plogis(1.5 * theta)

b_df <- 2 * eta
b_normal <- eta * (4 * pnorm(abs(mu) / 2) - 2)
```

A sensitivity-adjusted interval is `[L - b, U + b]`, where `[L, U]` is the unadjusted DID interval. It is not clipped to [−1, 1]. Rejection means exclusion of zero: Type I error at tau = 0 and power at tau != 0. Study 2 uses the data-generating eta and mu as fixed sensitivity parameters.

Condition seeds are `base_seed + 1009 * condition_id`, with Mersenne-Twister/Inversion/Rejection RNG settings. Output-only mode uses saved settings; manuscript outputs require alpha = .05. Rerunning overwrites generated files; change `OUTPUT_ROOT` for separate runs.

## PIAAC illustration

Download `prgkorp2.csv` from the [OECD PIAAC 2nd Cycle Database](https://www.oecd.org/en/data/datasets/PIAAC-2nd-Cycle-Database.html) and place it beside the script, or change `EMPIRICAL_DATA_FILE` to its local path. The reader expects the semicolon-delimited CSV, with missing-value codes `.`, `.n`, and `.v`.

The analysis uses test item `E320004S`, anchor item `E320003S`, and literacy plausible values `PVLIT1`–`PVLIT10`. `AGEG10LFS` codes 2 and 4 define reference ages 25–34 and focal ages 45–54. Complete responses to both items yield N = 536 (285 reference, 251 focal); the script checks these counts.

Uniform logistic-regression DIF is fitted separately for each plausible value and combined using multiple-imputation rules. DID uses `lm((Y_T - Y_A) ~ G)` with an HC3 standard error. Analyses are unweighted and use large-sample normal inference. The fixed sensitivity value is eta = .033; the standardized group mean difference is calculated from the plausible values. Uncertainty in the sensitivity parameters is not incorporated.

The unweighted results are illustrative, not population-representative estimates for Korean adults.

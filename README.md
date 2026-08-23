# Identifying Item Bias Without Conditioning: A Difference-in-Differences Approach

## Reproducibility Materials

This repository contains the R code used to reproduce the Monte Carlo simulations and empirical illustration reported in the manuscript **“Identifying Item Bias Without Conditioning: A Difference-in-Differences Approach.”**

The repository is anonymized for double-blind peer review. Author-identifying information is intentionally omitted and will be added after the review process, if appropriate.

## Repository Contents

### `simulation.R`

Reproduces the two Monte Carlo simulation studies.

**Simulation Study 1** evaluates the proposed difference-in-differences (DID) estimator when the matching variable used by conventional DIF analysis is invalid or unreliable. The simulation varies:

- sample size: \(N = 500, 1000, 2000\);
- average item-bias effect: \(\tau = 0, -.05, -.10\);
- matching-variable invalidity: \(\delta = 0, .25, .50\); and
- matching-variable reliability: \(\rho_M = 1.00, .80, .60\).

The conditioning-based comparator is uniform logistic-regression DIF. DID inference uses HC3 heteroskedasticity-robust standard errors.

**Simulation Study 2** evaluates robustness of DID to violations of item response function (IRF) equivalence. Difficulty and discrimination departures are calibrated to the same maximum reference-group IRF discrepancy,

\[
\epsilon
=
\sup_{\theta}
\left|
P_T(0,\theta)-P_A(\theta)
\right|,
\]

with \(\epsilon = 0, .05, .10\). The simulation separately evaluates sampling performance for the population DID and coverage of the target item-bias estimand \(\tau\).

The script writes condition-specific summaries and source data for the simulation tables and figures to the `results/` directory.

### `empirical_illustration.R`

Reproduces the empirical illustration using the Korean sample from the **2023 Programme for the International Assessment of Adult Competencies (PIAAC)**.

The analysis:

- treats adults aged 25–34 as the reference group and adults aged 45–54 as the focal group;
- examines Literacy item `E320004S` as the test item;
- uses `E320003S` as the anchor item;
- uses the ten PIAAC Literacy plausible values (`PVLIT1`–`PVLIT10`) as the matching variable for conventional uniform logistic-regression DIF;
- estimates the proposed DID on the response-probability scale using HC3 heteroskedasticity-robust standard errors; and
- conducts the distribution-free and equal-variance normal sensitivity analyses reported in the manuscript.

The analytic sample contains 536 respondents with observed responses to both items (285 reference-group and 251 focal-group respondents).

## Data

The PIAAC public-use microdata are **not redistributed in this repository**.

The empirical illustration requires the Korean 2023 PIAAC public-use data file, available from the OECD PIAAC data repository:

https://www.oecd.org/en/about/programmes/piaac/piaac-data.html

After obtaining the Korean public-use CSV file, place it in the working directory and set the following line in `empirical_illustration.R` to the appropriate local file name or path:

```r
DATA_FILE <- "prgkorp2.csv"
```

The OECD CSV is semicolon-delimited. The script imports the file using the delimiter and missing-value codes required for the public-use file.

## Software Requirements

The analyses were implemented in **R**. The scripts require the following packages:

```r
install.packages(c(
  "dplyr",
  "sandwich"
))
```

The scripts print `sessionInfo()` at completion to facilitate reproducibility.

## Running the Analyses

Clone or download this repository and set the repository root as the R working directory.

To reproduce the Monte Carlo simulations:

```r
source("simulation.R")
```

To reproduce the empirical illustration after obtaining the PIAAC data:

```r
source("empirical_illustration.R")
```

The simulation script uses 5,000 Monte Carlo replications per condition, as reported in the manuscript. Consequently, the full simulation may require substantial computation time.

Both scripts create a `results/` directory automatically and save the principal reproducibility outputs as CSV files.

## Expected Empirical Results

With the Korean 2023 PIAAC public-use data used in the manuscript, `empirical_illustration.R` should reproduce the following rounded results:

| Method | Estimate | SE | 95% CI | p |
|---|---:|---:|---:|---:|
| Logistic DIF | −.638 | .229 | [−1.086, −.190] | .005 |
| DID | −.119 | .042 | [−.201, −.036] | .005 |

The logistic DIF estimate is on the log-odds scale, whereas the DID estimate is on the response-probability scale.

For the sensitivity analysis, the manuscript uses \(\epsilon=.033\). The distribution-free sensitivity-adjusted 95% interval is approximately \([-.267,.030]\). Under the equal-variance normal specification, using the observed standardized Literacy difference of approximately \(|\mu|=.72\), the adjusted interval is approximately \([-.220,-.018]\).

Small differences in unrounded output may arise from software versions or numerical routines.

## Reproducibility Notes

The code is organized to correspond as closely as possible to the analyses reported in the manuscript. The simulation and empirical illustration are kept in separate scripts so that the computationally intensive Monte Carlo analyses can be run independently of the PIAAC application.

Random-number generation for the simulation is initialized with a fixed seed in `simulation.R`. Because the manuscript reports Monte Carlo results based on finite replications, exact reproduction of the reported simulation values requires the same code, seed, software environment, and random-number generation behavior.

The empirical illustration does not require access to any nonpublic or proprietary data beyond the publicly available OECD PIAAC public-use file.

## Double-Blind Review

This repository is prepared for double-blind peer review. It intentionally excludes author names, affiliations, personal contact information, and other direct author identifiers. Please do not cite or distribute the repository in a manner that could compromise anonymous review.

## License and Citation

Formal citation information and licensing details will be added following peer review. During the anonymous review period, please refer to the accompanying manuscript by title:

> *Identifying Item Bias Without Conditioning: A Difference-in-Differences Approach*


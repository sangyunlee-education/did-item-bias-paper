# Identifying Item Bias Without Conditioning: A Difference-in-Differences Approach

## Reproducibility Materials

This repository contains the R code used for the simulation studies and empirical illustration reported in the manuscript *Identifying Item Bias Without Conditioning: A Difference-in-Differences Approach*.

The repository is prepared for double-blind peer review and contains no author-identifying information.

## Files

### `simulation.R`

Reproduces Simulation Studies 1 and 2.

**Simulation Study 1** evaluates the proposed difference-in-differences (DID) approach when the matching variable used in conventional DIF analysis is invalid or unreliable.

The simulation varies:

- Sample size: `N = 500, 1000, 2000`
- Item-bias effect: `tau = 0, -.05, -.10`
- Matching-variable invalidity: `delta = 0, .25, .50`
- Matching-variable reliability: `rho_M = 1.00, .80, .60`

The comparison method is uniform logistic-regression DIF. DID inference uses HC3 heteroskedasticity-robust standard errors. Each condition uses 5,000 Monte Carlo replications.

**Simulation Study 2** evaluates DID under violations of item response function (IRF) equivalence. Difficulty and discrimination departures are calibrated to the same maximum IRF discrepancy:

```text
epsilon = sup_theta |P_T(0, theta) - P_A(theta)|
```

with `epsilon = 0, .05, .10`.

The script saves simulation summaries and source data for the reported tables and figures in the `results/` directory.

### `empirical_illustration.R`

Reproduces the empirical illustration using the Korean sample from the 2023 Programme for the International Assessment of Adult Competencies (PIAAC).

The analysis uses:

- Reference group: adults aged 25–34
- Focal group: adults aged 45–54
- Test item: `E320004S`
- Anchor item: `E320003S`
- Matching variable for logistic DIF: `PVLIT1`–`PVLIT10`
- Analytic sample: `N = 536` (285 reference; 251 focal)

The script estimates conventional uniform logistic-regression DIF, the proposed DID with HC3 heteroskedasticity-robust standard errors, and the sensitivity analyses reported in the manuscript.

## Data

The PIAAC public-use data are not included in this repository.

The Korean 2023 PIAAC public-use data can be obtained from the OECD PIAAC data repository:

https://www.oecd.org/en/about/programmes/piaac/piaac-data.html

After downloading the Korean CSV file, place it in the working directory and specify its path in `empirical_illustration.R`:

```r
DATA_FILE <- "prgkorp2.csv"
```

## Requirements

The analyses require R and the following packages:

```r
install.packages(c("dplyr", "sandwich"))
```

## Running the Code

Run the simulation studies with:

```r
source("simulation.R")
```

Run the empirical illustration with:

```r
source("empirical_illustration.R")
```

Both scripts create a `results/` directory for the main output files.

## Expected Empirical Results

Running `empirical_illustration.R` with the Korean 2023 PIAAC public-use data should reproduce the following rounded results:

| Method | Estimate | SE | 95% CI | p |
| --- | ---: | ---: | ---: | ---: |
| Logistic DIF | -.638 | .229 | [-1.086, -.190] | .005 |
| DID | -.119 | .042 | [-.201, -.036] | .005 |

The logistic DIF estimate is on the log-odds scale. The DID estimate is on the response-probability scale.

The sensitivity analysis uses `epsilon = .033`. The sensitivity-adjusted 95% intervals are:

- Distribution-free: `[-.267, .030]`
- Equal-variance normal specification with `|mu| ≈ .72`: `[-.220, -.018]`

## Reproducibility

`simulation.R` uses a fixed random seed and 5,000 replications per condition. Exact Monte Carlo results may depend on the R version, package versions, and random-number generation environment.

Both scripts print `sessionInfo()` at completion.

## Anonymous Review

These materials are provided for double-blind peer review. Author names, affiliations, and contact information are intentionally omitted.

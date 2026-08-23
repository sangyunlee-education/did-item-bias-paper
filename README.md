# Identifying Item Bias Without Conditioning: A Difference-in-Differences Approach

## Reproducibility Materials

This repository contains the R code for the simulation studies and empirical illustration reported in the manuscript *Identifying Item Bias Without Conditioning: A Difference-in-Differences Approach*.

The repository is prepared for double-blind peer review and contains no author-identifying information.

## Files

### `simulation.R`

Contains the code for Simulation Studies 1 and 2.

**Simulation Study 1** evaluates the proposed difference-in-differences (DID) approach when the matching variable used in conventional DIF analysis is invalid or unreliable. The simulation varies:

- Sample size: `N = 500, 1000, 2000`
- Item-bias effect: `tau = 0, -.05, -.10`
- Matching-variable invalidity: `delta = 0, .25, .50`
- Matching-variable reliability: `rho_M = 1.00, .80, .60`

The comparison method is uniform logistic-regression DIF. DID inference uses HC3 heteroskedasticity-robust standard errors. Each condition uses 5,000 Monte Carlo replications.

**Simulation Study 2** evaluates the robustness of DID to violations of item response function (IRF) equivalence. Difficulty and discrimination departures are calibrated to the same maximum reference-group IRF discrepancy:

```text
epsilon = sup_theta |P_T(0, theta) - P_A(theta)|
```

The simulation considers `epsilon = 0, .05, .10` and evaluates the consequences of IRF nonequivalence for identification and inference.

The script saves the principal simulation outputs to the `results/` directory.

### `empirical_illustration.R`

Contains the code for the empirical illustration using the Korean sample from the 2023 Programme for the International Assessment of Adult Competencies (PIAAC).

The analysis uses:

- Reference group: adults aged 25–34
- Focal group: adults aged 45–54
- Test item: `E320004S`
- Anchor item: `E320003S`
- Matching variable for logistic DIF: `PVLIT1`–`PVLIT10`
- Analytic sample: `N = 536` (285 reference; 251 focal)

The script conducts conventional uniform logistic-regression DIF analysis, estimates the proposed DID using HC3 heteroskedasticity-robust standard errors, and performs the sensitivity analyses reported in the manuscript.

## Data

The PIAAC public-use microdata are not included in this repository.

The Korean 2023 PIAAC public-use data are available from the OECD PIAAC data repository:

https://www.oecd.org/en/about/programmes/piaac/piaac-data.html

After downloading the Korean CSV file, place it in the working directory and specify its local path in `empirical_illustration.R`:

```r
DATA_FILE <- "prgkorp2.csv"
```

## Requirements

The analyses require R and the following packages:

```r
install.packages(c("dplyr", "sandwich"))
```

## Running the Code

To run the simulation studies:

```r
source("simulation.R")
```

To run the empirical illustration after obtaining the PIAAC data:

```r
source("empirical_illustration.R")
```

Both scripts create a `results/` directory for the principal output files and print `sessionInfo()` at completion.

## Reproducibility Notes

`simulation.R` uses a fixed random seed and 5,000 Monte Carlo replications per condition. Because the reported simulation results are based on finite Monte Carlo samples, exact numerical reproduction may depend on the R version, package versions, and random-number generation environment.

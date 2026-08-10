# Simulation Code

This directory contains the R code for the two Monte Carlo simulation studies in:

**Identifying Item Bias Without Conditioning: A Difference-in-Differences Approach**

## Files

- `00_helpers.R`: latent-trait generation, numerical calibration, and HC3 inference functions.
- `01_study1_exact_equivalence.R`: Study 1 under exact item response function equivalence.
- `02_study2_falsification_diagnostic.R`: Study 2 for the falsification diagnostic.
- `03_summarize_simulation_results.R`: creates the compact summaries reported in the main text.

## Requirements

The simulation code uses **base R only**. No additional R packages are required.

## Running the simulations

Run all scripts from the repository root in this order:

```r
source("R/simulation/01_study1_exact_equivalence.R")
source("R/simulation/02_study2_falsification_diagnostic.R")
source("R/simulation/03_summarize_simulation_results.R")
```

Each simulation condition uses 5,000 Monte Carlo replications. For a quick check before the full run, temporarily change `N_REP <- 5000L` to a smaller value such as `50L` in the Study 1 and Study 2 scripts.

Outputs are written automatically to:

```text
results/simulation/
```

## Reproducibility

Each condition uses a fixed random-number seed. The code therefore reproduces the same Monte Carlo results whenever it is run with the same R implementation and settings.

The numerical calibration reproduces the population targets described in the manuscript. In Study 1, `delta_T` is calibrated separately for each focal-group latent-trait distribution and target value of `tau`. In Study 2, `c(alpha)` is calibrated so that `E(Y_T - Y_A | G = 0) = alpha`.

## Note on the manuscript tables

The present scripts use explicit fixed seeds for reproducibility. If the originally reported simulations were generated under a different or undocumented seed, Monte Carlo summaries may differ slightly from the values currently printed in the manuscript while converging to the same population targets and operating characteristics.

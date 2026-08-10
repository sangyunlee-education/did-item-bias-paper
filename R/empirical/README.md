# Empirical Illustration Code

This directory contains the R code for the empirical illustration in:

**Identifying Item Bias Without Conditioning: A Difference-in-Differences Approach**

## Files

- `00_prepare_empirical_data.R`: prepares the Korean PIAAC Cycle 2 locator-test analysis data.
- `01_conventional_dif_screening.R`: screens the eight items using Mantel--Haenszel, logistic-regression DIF, SIBTEST, and a 2PL IRT likelihood-ratio test.
- `02_pairwise_falsification_did.R`: runs the pairwise falsification diagnostic and DID analysis using the screened candidate anchors.
- `03_summarize_empirical_results.R`: creates formatted versions of the screening and pairwise result tables and prints the key counts reported in the manuscript.

## Requirements

The analysis uses the following R packages:

```r
install.packages(c("difR", "mirt", "sandwich"))
```

## Data

The empirical illustration uses the Korean Public Use File from PIAAC Cycle 2.

The raw PIAAC data are not included in this repository. Place the Korea CSV Public Use File in the `data/` directory and edit the short configuration block at the top of:

```text
R/empirical/00_prepare_empirical_data.R
```

The manuscript uses the following eight Numeracy locator-test items:

```text
C601C06S
C815P001S
C815P002S
C832P001S
C832P002S
C813P001S
C833P001S
C833P002S
```

Correct responses are coded `1`; incorrect responses and omissions are coded `0`. Males are the reference group (`G = 0`) and females are the focal group (`G = 1`). The complete-case analytic sample reported in the manuscript is `N = 5,780`, with 2,715 reference-group and 3,065 focal-group respondents.

## Running the empirical illustration

Run the scripts from the repository root in this order:

```r
source("R/empirical/00_prepare_empirical_data.R")
source("R/empirical/01_conventional_dif_screening.R")
source("R/empirical/02_pairwise_falsification_did.R")
source("R/empirical/03_summarize_empirical_results.R")
```

Outputs are written automatically to:

```text
results/empirical/
```

## Analysis workflow

### Stage 1: Conventional DIF screening

Each of the eight items is treated as the studied item in turn, with the remaining seven items providing the basis for matching or linking.

- Mantel--Haenszel uses the sum of the remaining seven items as the matching score.
- Logistic-regression DIF uses the standardized sum of the remaining seven items and a likelihood-ratio test of the group and group-by-score terms.
- SIBTEST uses the remaining seven items as the matching set.
- The 2PL IRT likelihood-ratio test constrains the remaining seven items as linking items while allowing the focal-group latent mean and variance to differ.

An item is retained as a candidate anchor only when none of the four procedures rejects at the `.05` level. The p-values are unadjusted.

### Stage 2: Falsification diagnostic and DID

Each remaining test item is paired with every candidate anchor. For each pair,

```text
Y_T - Y_A
```

is regressed on group membership. The intercept estimates

```text
alpha = E(Y_T - Y_A | G = 0),
```

and the group coefficient estimates the DID.

The falsification diagnostic tests `H0: alpha = 0`. HC3 heteroskedasticity-robust standard errors, two-sided tests, and 95% confidence intervals are used. No adjustment is made for the candidate-anchor screening stage or for multiple anchor--test comparisons.

## Expected screening outcome

The manuscript reports three candidate anchors:

```text
C601C06S
C815P001S
C833P001S
```

Pairing these three candidate anchors with the five remaining test items yields 15 anchor--test pairs.

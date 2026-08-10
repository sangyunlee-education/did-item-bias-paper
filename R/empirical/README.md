# Empirical Illustration Code

This directory contains the R code for the empirical illustration in:

**Identifying Item Bias Without Conditioning: A Difference-in-Differences Approach**

## Files

- `00_prepare_empirical_data.R`: reads and recodes the Korean PIAAC Cycle 2 data.
- `01_conventional_dif_screening.R`: performs the four leave-one-item-out conventional DIF screening procedures.
- `02_pairwise_falsification_did.R`: performs the pairwise falsification diagnostic and DID analyses using HC3 inference.
- `03_summarize_empirical_results.R`: formats the manuscript tables and checks the reproduced results against the values reported in the manuscript.

## R packages

```r
install.packages(c(
  "dplyr",
  "tidyr",
  "purrr",
  "tibble",
  "mirt",
  "difR",
  "sandwich"
))
```

## Data

The raw PIAAC data are not included in this repository.

The scripts use the Korea PIAAC Cycle 2 file:

```text
prgkorp2.csv
```

The preferred location is:

```text
data/piaac/prgkorp2.csv
```

For convenience, `00_prepare_empirical_data.R` also accepts `prgkorp2.csv` in the repository root.

The analysis uses:

- `GENDER_R = 1`: male/reference group (`G = 0`)
- `GENDER_R = 2`: female/focal group (`G = 1`)
- item score `1`: correct
- item scores `0` and `7`: incorrect/omitted
- all other values: missing

The eight Numeracy locator-test items are:

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

Complete-case analysis should yield `N = 5,780`, with 2,715 reference-group and 3,065 focal-group respondents.

## Running the analysis

Run the scripts from the repository root in this order:

```r
source("R/empirical/00_prepare_empirical_data.R")
source("R/empirical/01_conventional_dif_screening.R")
source("R/empirical/02_pairwise_falsification_did.R")
source("R/empirical/03_summarize_empirical_results.R")
```

The last script checks the reproduced results against the manuscript. A successful reproduction should end with:

```text
OVERALL: PASS
```

Outputs are written to:

```text
results/empirical/
```

## Analysis details

For each studied item, the remaining seven items provide the matching or linking set.

- Mantel--Haenszel uses the sum of the other seven items as the matching score and reproduces the manuscript analysis with `correct = FALSE`.
- Logistic-regression DIF compares a reduced model using the standardized matching score with a full model including group and the score-by-group interaction.
- SIBTEST uses the remaining items as the matching set.
- The 2PL IRT likelihood-ratio test frees the studied item across groups while the remaining seven items link the group scales.

An item is retained as a candidate anchor only when none of the four unadjusted tests rejects at `.05`.

Each remaining test item is then paired with every candidate anchor. For each pair, `Y_T - Y_A` is regressed on group membership. The intercept is the falsification-diagnostic target `alpha`, and the group coefficient is the DID. HC3 heteroskedasticity-robust standard errors, two-sided tests, and 95% confidence intervals are used.

The main empirical table includes all 15 test--anchor pairs, while DID estimates from pairs for which the falsification diagnostic rejects are not given an item-bias interpretation.

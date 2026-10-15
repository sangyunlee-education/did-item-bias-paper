# Identifying Item Bias Without Conditioning: A Difference-in-Differences Approach

`reproduce_did_study.R` reproduces the two simulation studies, PIAAC illustration, Figures 4–5, and Tables S1–S4. It requires R 3.6 or later and no additional R packages. PDF figures use Arial and a Quartz or Cairo device.

## Run

1. Obtain the Korean 2023 PIAAC public-use CSV from the [OECD PIAAC database](https://www.oecd.org/en/data/datasets/PIAAC-2nd-Cycle-Database.html).
2. Set `PIAAC_FILE` at the top of the script to its location (default: `prgkorp2.csv`). Use the original semicolon-delimited CSV; the script handles item-response recoding.
3. Set the repository folder as the working directory and run the entire script using **Source** in RStudio, or:

```sh
Rscript reproduce_did_study.R
```

The script runs Study 1, Study 2, the empirical illustration, and figure and table generation in sequence. Both simulation studies use 5,000 replications per condition. Each full run generates new samples using fixed seeds and replaces outputs with the same filenames.

## Outputs

- `results/study1/`: Study 1 results, Figure 4, and Tables S1–S2.
- `results/study2/`: Study 2 results, Figure 5, and Tables S3–S4.
- `results/empirical/`: DIF and DID estimates, sensitivity bounds, and confidence intervals.
- `results/run_records/`: analysis settings, input-data checksum, and R session information.

LaTeX tables use the manuscript’s CUP macros. Insert Figures 4 and 5 at the same width; the default final width is 144 mm.

PIAAC analyses are unweighted and do not account for the complex sampling design. Sensitivity inputs are treated as fixed. Small numerical differences across platforms or R versions are possible.

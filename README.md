# Identifying Item Bias Without Conditioning: A Difference-in-Differences Approach

`reproduce_did_study.R` reproduces the two simulation studies, PIAAC illustration, Figures 4–5, and Tables S1–S4. It requires R 3.6 or later and no additional R packages. PDF figures use Arial and a Quartz or Cairo device.

## Run

Run from the repository folder. In RStudio, set this folder as the working directory, choose `TASK` at the top of the script, and click **Source**. Alternatively:

```sh
Rscript reproduce_did_study.R study1
Rscript reproduce_did_study.R study2
Rscript reproduce_did_study.R empirical
Rscript reproduce_did_study.R figure4
Rscript reproduce_did_study.R figure5
Rscript reproduce_did_study.R tables
```

For `study1` and `study2`, set `DATA_SOURCE` in the script:

- `"saved"` (default): recompute summaries from saved `simulation.rds` files.
- `"new"`: generate samples using `N_REP` replications per condition (default: 5,000).

To reproduce from scratch, run `study1` and `study2` with `DATA_SOURCE <- "new"`, then run `figure4`, `figure5`, and `tables`.

Figure and table tasks read CSV summaries without running simulations:

- Figure 4: `study1/did_summary.csv` and `study1/logistic_summary.csv`
- Figure 5: `study2/main_summary.csv`
- Tables S1–S4: all three files above and `study2/eta_mu_summary.csv`

All paths are relative to `OUTPUT_DIR` (default: `results/`). Saved simulation files are `study1/simulation.rds` and `study2/simulation.rds`; use `STUDY1_FILE` and `STUDY2_FILE` to specify other locations. Generated outputs overwrite files with the same names.

## PIAAC illustration

Obtain the Korean PIAAC public-use CSV separately and set `PIAAC_FILE`. The semicolon-delimited file must contain `AGEG10LFS`, `E320004S`, `E320003S`, and `PVLIT1`–`PVLIT10`.

Analyses are unweighted and do not account for the complex sampling design. Sensitivity inputs are treated as fixed.

## Outputs

Results are saved under `results/study1/`, `results/study2/`, and `results/empirical/`. Run settings and R session information are saved under `results/run_records/`.

Figure 4 reports rejection rates. Figure 5 plots coverage against mean interval length, with IRF discrepancy levels in rows and sample sizes in columns. Tables S3–S4 report coverage and mean interval length; S4 varies the specified values of both sensitivity parameters, eta and mu.

LaTeX tables use the manuscript's CUP macros. Insert Figures 4 and 5 at the same width; `FINAL_FIGURE_WIDTH_MM` defaults to 144 mm.

Simulation seeds are fixed in the script. Small numerical differences across platforms or R versions are possible; reusing saved simulation files avoids generating new samples.

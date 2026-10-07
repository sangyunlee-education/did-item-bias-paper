# Identifying Item Bias Without Conditioning: A Difference-in-Differences Approach

`reproduce_did_study.R` reproduces the two simulation studies, PIAAC illustration, Figures 4–5, and Tables S1–S4. It requires R 3.6 or later and no additional R packages. Figures use Arial and a Quartz or Cairo PDF device.

## Run

Run from the repository folder. In RStudio, set that folder as the working directory and choose `RUN_MODE` at the top of the script.

```sh
# Update figures and tables from existing CSV summaries (default).
Rscript reproduce_did_study.R presentation

# Run both simulation studies with 5,000 replications per condition.
Rscript reproduce_did_study.R simulations

# Recompute summaries and outputs from saved raw RDS files.
Rscript reproduce_did_study.R outputs

# Run the PIAAC illustration.
Rscript reproduce_did_study.R empirical
```

`presentation` requires these four files under `OUTPUT_DIR` (default: `results/`):

- `study1/did_summary.csv`
- `study1/logistic_summary.csv`
- `study2/main_summary.csv`
- `study2/eta_summary.csv`

It generates no samples and leaves existing numerical results and analysis records unchanged. Missing summaries cause an error, not a new simulation run. Study 2 summaries must include coverage, mean interval length, and rejection rates.

`outputs` reads `study1/simulation.rds` and `study2/simulation.rds` under `OUTPUT_DIR`. Set `STUDY1_SAVED` and `STUDY2_SAVED` for other compatible RDS paths. Fresh simulation runs replace generated results; change `OUTPUT_DIR` to retain an earlier run.

Other modes are `study1`, `study2`, `all` (both studies and the illustration), and `smoke` (20 replications per condition in `results_smoke/`).

## Data and outputs

For the illustration, obtain the Korean PIAAC public-use CSV separately and set `PIAAC_FILE`. The semicolon-delimited file must contain `AGEG10LFS`, `E320004S`, `E320003S`, and `PVLIT1`–`PVLIT10`. Analyses are unweighted; sensitivity inputs are treated as fixed.

Results are saved under `results/study1/`, `results/study2/`, and `results/empirical/`. Figure 5 shows coverage (A), mean interval length (B), and rejection rates (C). Tables S3 and S4 report all three measures. LaTeX tables use the manuscript's CUP macros.

Insert Figures 4 and 5 at the same width. The default final width is 144 mm; adjust `FINAL_FIGURE_WIDTH_MM` if needed. Set `MAKE_FIGURES <- FALSE` to skip PDF output.

Simulation seeds are `2026 + 1009 * condition_id` for Study 1 and `2027 + 1009 * condition_id` for Study 2. Small numerical differences across platforms or R versions are possible. Reusing saved results avoids generating new samples. Run settings and R session information are saved separately for analysis and presentation runs.

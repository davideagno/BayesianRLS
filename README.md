# BayesianRLS

Code to reproduce the results contained in the paper:

Agnoletto, D, Rigon, T, and Dunson, D. B. (2026+) **Martingale posteriors for flexible regression via
recursive least squares**, arXiv.

The shared functions are in `RLS_functions.R` and `RLS_functions.cpp`. Each folder corresponds to a
section of the paper or of the supplementary material. Its scripts are run from that folder, and
write their results (`.RData` files, figures and summary tables) to the subfolder `output/`.

| Folder | Section | Output | Scripts |
|---|---|---|---|
| `sec5.1_uq_gaussian/` | 5.1 | Table 1 | `sim_normal.R`, then `results.R` |
| `sec5.2_additive_models/` | 5.2 | Table 2, Table S3 | `1_sim_data.R`, `2_rls.R`, `2_rls_norm.R`, `3_trad.R`, `4_bart.R`, `5_results.R`, in this order |
| `sec6_application/` | 6, S3, S4 | figures and results of the application | `data.R` (downloads the NHANES data), then `analysis.R` |
| `secS2.1_illustration/` | S2.1 | Figures S1 and S2 | `student_t.R`, `skew_t.R` |
| `secS2.2_uq_student_t/` | S2.2 | Table S1 | `sim_t.R`, then `../sec5.1_uq_gaussian/results.R` with argument `t` |
| `secS2.3_frequentist_am/` | S2.3 | Table S2 | `freq_AM.R`, then `freq_AM_results.R`; both use the output of `sec5.2_additive_models/` |

The sample size of the simulations is passed as the first argument, as a multiple of the number
of covariates: for instance, `Rscript 2_rls.R 100` runs the scenario with n = 100p. Without
arguments, the first scenario of the paper is used. The summaries of Table 1 and Table S1 take a
second argument, `norm` or `t`, for instance `Rscript results.R 10 t`.

`run_all.sh` runs everything in the right order, or only the given folders
(e.g. `./run_all.sh sec5.2_additive_models`). On macOS it keeps the machine awake, and its log
`run_all.log` reports for each step the elapsed time and the number of system sleep events: the
computing times are reliable only for steps with no sleep events.

R packages: `Rcpp`, `RcppArmadillo`, `mgcv`, `coda`, `brms`, `BART`, `robustBLME`, `extraDistr`,
`parallel`, `ggplot2`, `ggthemes`, `ggpubr`, `ggcorrplot`, `latex2exp`, `gridExtra`, `haven`, `dplyr`.

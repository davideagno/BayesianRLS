# BayesianRLS

Code to reproduce the results contained in the paper:

Agnoletto, D, Rigon, T, and Dunson, D. B. (2026+) **Martingale posteriors for flexible regression via
recursive least squares**, arXiv.

The shared functions are in `RLS_functions.R` and `RLS_functions.cpp`. Each folder corresponds to a
section of the paper or of the supplementary material. Its scripts are run from that folder, and
write their results (`.RData` files, figures and summary tables) to the subfolder `output/`.
The figures of the paper are written to `../../PAPER/figures` when that folder exists (as in the
authors' project), and to `output/` otherwise.
The `.RData` files contain, besides the metrics of the tables, the posterior summaries of each
dataset (means, standard deviations, HPD intervals and quantiles at each test point or for each
coefficient, with the true values), so that the summary scripts can be modified without
rerunning the simulations. Each simulation study uses 30 datasets per scenario.

| Folder | Section | Output | Scripts, in order |
|---|---|---|---|
| `sec5.1_uq_gaussian/` | 5.1 | Table 1 | `simulate.R`, `summarize_tables.R` |
| `sec5.2_illustration/` | 5.2, S2.1 | Figure 1, Figures S1 and S2 | `student_t_errors.R`, `skew_t_errors.R`, `posteriors_figure.R` |
| `sec5.3_additive_models/` | 5.3, S2.4 | Table 2, Table S3 | `1_simulate_data.R`, `2_fit_rls_student_t.R`, `3_fit_rls_gaussian.R`, `4_fit_blm_brms.R`, `5_fit_bart.R`, `6_summarize_table.R` |
| `sec6_application/` | 6, S3, S4 | Figure 2, Figures S3 to S8, results of Section 6 | `1_download_data.R` (downloads the NHANES data into `appl_data.RData`, already included), `2_analysis.R` |
| `secS2.2_uq_student_t/` | S2.2 | Table S1 | `simulate.R`, then `../sec5.1_uq_gaussian/summarize_tables.R` with argument `t` |
| `secS2.3_frequentist_am/` | S2.3 | Table S2 | `1_fit_frequentist_am.R`, `2_summarize_table.R`; both use the output of `sec5.3_additive_models/` |

The sample size of the simulations is passed as the first argument, as a multiple of the number
of covariates: for instance, `Rscript 2_fit_rls_student_t.R 100` runs the scenario with n = 100p. Without
arguments, the first scenario of the paper is used. The summaries of Table 1 and Table S1 take a
second argument, `norm` or `t`, for instance `Rscript summarize_tables.R 10 t`.

`run_all.sh` runs everything in the right order, or only the given folders
(e.g. `./run_all.sh sec5.3_additive_models`). On macOS it keeps the machine awake, and its log
`run_all.log` reports for each step the elapsed time and the number of system sleep events: the
computing times are reliable only for steps with no sleep events.

In `sec6_application/2_analysis.R`, the traditional Bayesian model is also fitted with a vague
Gamma(0.001, 0.001) prior on the smoothing parameters, for the sensitivity analysis of Section 6;
the trace plot of a selected coefficient is saved in `output/` as a diagnostic.

R packages: `Rcpp`, `RcppArmadillo`, `mgcv`, `coda`, `brms` (with `rstan`), `BART`, `robustBLME`,
`extraDistr`, `parallel`, `ggplot2`, `ggthemes`, `ggpubr`, `ggcorrplot`, `latex2exp`, `scales`,
`haven`, `dplyr`, `sn`.

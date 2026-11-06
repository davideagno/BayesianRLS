# Summaries for Table S2; also uses the RLS-t results of
# ../sec5.3_additive_models/2_fit_rls_student_t.R.

rm(list = ls())

library(mgcv)
library(coda)
library(brms)

p <- 15
n_mult <- if (length(commandArgs(TRUE)) > 0) as.numeric(commandArgs(TRUE)[1]) else 30  # n = n_mult * p
n_obs <- n_mult*p
df_true <- 5
sigma_true <- 1
sink(file.path("output", paste0("summary_n", n_mult, "p.txt")), split = TRUE)

load(paste0("../sec5.3_additive_models/output/sim_p", p, "_n", n_obs, "_df", df_true, "_sigma", sigma_true, "_rls_alt.RData"))
load(file.path("output", paste0("sim_p", p, "_n", n_obs, "_df", df_true, "_sigma", sigma_true, "_AM.RData")))


round(c(mean(mse_rls), mean(mse_apx), mean(mse_amt), mean(mse_amn)), 3)
round(c(mean(rowMeans(in_ci_rls)), mean(rowMeans(in_ci_apx)), mean(rowMeans(in_ci_amt)), mean(rowMeans(in_ci_amn))), 3)
round(c(mean(width_ci_rls), mean(width_apx), mean(width_amt), mean(width_amn)), 3)
round(c(mean(time_rls), mean(time_apx), mean(time_amt), mean(time_amn)), 2)


round(c(sd(mse_rls), sd(mse_apx), sd(mse_amt), sd(mse_amn)), 3)
# Standard deviation across datasets, as for the other rows
round(c(sd(colMeans(in_ci_rls)), sd(colMeans(in_ci_apx)), sd(colMeans(in_ci_amt)), sd(colMeans(in_ci_amn))), 3)
round(c(sd(colMeans(width_ci_rls)), sd(width_apx), sd(width_amt), sd(width_amn)), 3)
round(c(sd(time_rls), sd(time_apx), sd(time_amt), sd(time_amn)), 2)

sink()

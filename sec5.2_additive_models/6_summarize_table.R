# Section 5.2, step 5: summaries for Table 2.

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

load(file.path("output", paste0("sim_p", p, "_n", n_obs, "_df", df_true, "_sigma", sigma_true, "_rls_alt.RData")))
idx_rls <- which(is.na(mse_rls) == FALSE)
length(idx_rls)

load(file.path("output", paste0("sim_p", p, "_n", n_obs, "_df", df_true, "_sigma", sigma_true, "_rls_norm_alt.RData")))
idx_rls_norm <- which(is.na(mse_rls_norm) == FALSE)
length(idx_rls_norm)

load(file.path("output", paste0("sim_p", p, "_n", n_obs, "_df", df_true, "_sigma", sigma_true, "_trad.RData")))
idx_trad <- which(is.na(mse_trad) == FALSE)
length(idx_trad)

load(file.path("output", paste0("sim_p", p, "_n", n_obs, "_df", df_true, "_sigma", sigma_true, "_bart.RData")))
idx_bart <- which(is.na(mse_bart) == FALSE)
length(idx_bart)


round(c(mean(mse_rls), mean(mse_rls_norm), mean(mse_trad[idx_trad]), mean(mse_bart)), 3)
round(c(mean(rowMeans(in_ci_rls)), mean(rowMeans(in_ci_rls_norm)), mean(rowMeans(in_ci_trad[,idx_trad])), mean(rowMeans(in_ci_bart))), 3)
round(c(mean(rowMeans(width_ci_rls)), mean(rowMeans(width_ci_rls_norm)), mean(rowMeans(width_ci_trad[,idx_trad])), mean(rowMeans(width_ci_bart))), 3)
round(c(mean(mse_rls_q), mean(mse_rls_q_norm), mean(mse_trad_q[idx_trad]), mean(mse_bart_q)), 3)
round(c(mean(rowMeans(in_ci_rls_q)), mean(rowMeans(in_ci_rls_q_norm)), mean(rowMeans(in_ci_trad_q[,idx_trad])), mean(rowMeans(in_ci_bart_q))), 3)
round(c(mean(rowMeans(width_rls_q)), mean(rowMeans(width_rls_q_norm)), mean(rowMeans(width_trad_q[,idx_trad])), mean(rowMeans(width_bart_q))), 3)
round(c(mean(time_est+time_opt), mean(time_opt_norm+time_est_norm), mean(time_trad[idx_trad]), mean(time_bart)), 2)
round(c(mean(time_rls), mean(time_rls_norm), mean(time_trad[idx_trad]), mean(time_bart)), 2)


round(c(sd(mse_rls), sd(mse_rls_norm), sd(mse_trad[idx_trad]), sd(mse_bart)), 3)
round(c(sd(rowMeans(in_ci_rls)), sd(rowMeans(in_ci_rls_norm)), sd(rowMeans(in_ci_trad[,idx_trad])), sd(rowMeans(in_ci_bart))), 3)
round(c(sd(rowMeans(width_ci_rls)), sd(rowMeans(width_ci_rls_norm)), sd(rowMeans(width_ci_trad[,idx_trad])), sd(rowMeans(width_ci_bart))), 3)
round(c(sd(mse_rls_q), sd(mse_rls_q_norm), sd(mse_trad_q[idx_trad]), sd(mse_bart_q)), 3)
round(c(sd(rowMeans(in_ci_rls_q)), sd(rowMeans(in_ci_rls_q_norm)), sd(rowMeans(in_ci_trad_q[,idx_trad])), sd(rowMeans(in_ci_bart_q))), 3)
round(c(sd(rowMeans(width_rls_q)), sd(rowMeans(width_rls_q_norm)), sd(rowMeans(width_trad_q[,idx_trad])), sd(rowMeans(width_bart_q))), 3)
round(c(sd(time_est+time_opt), sd(time_opt_norm+time_est_norm), sd(time_trad[idx_trad]), sd(time_bart)), 2)
round(c(sd(time_rls), sd(time_rls_norm), sd(time_trad[idx_trad]), sd(time_bart)), 2)

# Table S3: average ESS of the draws of f(x) across the test set (mean and sd over datasets).
# The predictive resampling draws are independent, so their ESS is B by construction; the
# values estimated by coda are reported only as a check.
round(c(mean(ess_rls), mean(ess_rls_norm), mean(ess_trad[,idx_trad]), mean(ess_bart)), 1)
round(c(sd(colMeans(ess_rls)), sd(colMeans(ess_rls_norm)), sd(colMeans(ess_trad[,idx_trad])), sd(colMeans(ess_bart))), 1)

sink()

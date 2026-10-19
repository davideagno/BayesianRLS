
rm(list = ls())

library(mgcv)
library(coda)
library(brms)

p <- 15
n_obs <- 30*p
df_true <- 5
sigma_true <- 1

load(paste0("sim_p", p, "_n", n_obs, "_df", df_true, "_sigma", sigma_true, "_rls_alt.RData"))
idx_rls <- which(is.na(mse_rls) == FALSE)
length(idx_rls)

load(paste0("sim_p", p, "_n", n_obs, "_df", df_true, "_sigma", sigma_true, "_rls_norm_alt.RData"))
idx_rls_norm <- which(is.na(mse_rls_norm) == FALSE)
length(idx_rls_norm)

load(paste0("sim_p", p, "_n", n_obs, "_df", df_true, "_sigma", sigma_true, "_trad.RData"))
idx_trad <- which(is.na(mse_trad) == FALSE)
length(idx_trad)

load(paste0("sim_p", p, "_n", n_obs, "_df", df_true, "_sigma", sigma_true, "_bart.RData"))
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




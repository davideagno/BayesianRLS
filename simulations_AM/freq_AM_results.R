
rm(list = ls())

library(mgcv)
library(coda)
library(brms)

p <- 15
n_obs <- 100*p
df_true <- 5
sigma_true <- 1

load(paste0("sim_p", p, "_n", n_obs, "_df", df_true, "_sigma", sigma_true, "_rls_alt.RData"))
load(paste0("sim_p", p, "_n", n_obs, "_df", df_true, "_sigma", sigma_true, "_AM.RData"))


round(c(mean(mse_rls), mean(mse_apx), mean(mse_amt), mean(mse_amn)), 3)
round(c(mean(rowMeans(in_ci_rls)), mean(rowMeans(in_ci_apx)), mean(rowMeans(in_ci_amt)), mean(rowMeans(in_ci_amn))), 3)
round(c(mean(width_ci_rls), mean(width_apx), mean(width_amt), mean(width_amn)), 3)
round(c(mean(time_rls), mean(time_apx), mean(time_amt), mean(time_amn)), 2)


round(c(sd(mse_rls), sd(mse_apx), sd(mse_amt), sd(mse_amn)), 3)
round(c(sd(rowMeans(in_ci_rls)), sd(rowMeans(in_ci_apx)), sd(rowMeans(in_ci_amt)), sd(rowMeans(in_ci_amn))), 3)
round(c(sd(width_ci_rls), sd(width_apx), sd(width_amt), sd(width_amn)), 3)
round(c(sd(time_rls), sd(time_apx), sd(time_amt), sd(time_amn)), 2)



# Section S2.1 of the supplement: uncertainty quantification with Student-t
# responses (Table S1). Run from this folder, with n_obs = 3*p or 10*p; summaries
# with ../sec5.1_uq_gaussian/summarize_tables.R.

rm(list = ls())

library(mgcv)
library(coda)
library(brms)
library(parallel)

Rcpp::sourceCpp("../RLS_functions.cpp")
source("../RLS_functions.R")

# Preliminary settings
p <- 50
n_mult <- if (length(commandArgs(TRUE)) > 0) as.numeric(commandArgs(TRUE)[1]) else 3  # n = n_mult * p
n_obs <- n_mult*p
dir.create("output", showWarnings = FALSE)
n_all <- 1e5
sigma_true <- 2
nSim <- 30
df_true <- 5
var_true <- (sigma_true^2)*(df_true / (df_true-2))

B_ps <- 10000
N_ps <- 9 * n_obs   # N = 10 n, with the Gaussian correction of Algorithms 1 and 2

# beta_true, covariates and responses
set.seed(1997); beta_true <- rnorm(p, sd = 2)

fr_Sn_A <- fr_Sn_B <- fr_Sn_1 <- fr_Sn_2 <- fr_Sn_Vn <- rep(NA, nSim)
fr_Vn_A <- fr_Vn_B <- fr_Vn_1 <- fr_Vn_2 <- rep(NA, nSim)
time_A <- time_B <- time_1 <- time_2 <- rep(NA, nSim)
mse_A <- mse_B <- mse_1 <- mse_2 <- rep(NA, nSim)
mse_ex_Vn <- mse_ex_Sn <- fr_Vn_ex <- fr_Sn_ex <- frac_omit <- rep(NA, nSim)
# For any further analysis: posterior mean, covariance and summaries (mean, sd, HPD,
# 2.5% and 97.5% quantiles) of each method, and the exact quantities beta_n, V_n,
# Sigma_n / sigma^2 and V_N
post_mean_1 <- post_mean_2 <- post_mean_A <- post_mean_B <- beta_n_all <- matrix(NA, nrow = p, ncol = nSim)
post_cov_1 <- post_cov_2 <- post_cov_A <- post_cov_B <- array(NA, dim = c(p, p, nSim))
post_sum_1 <- post_sum_2 <- post_sum_A <- post_sum_B <- array(NA, dim = c(p, 6, nSim))
V_n_all <- Sigma_n_all <- V_N_all <- array(NA, dim = c(p, p, nSim))
in_ci_A <- in_ci_B <- in_ci_1 <- in_ci_2 <- in_ci_Vn <- in_ci_Sn <- matrix(NA, nrow = p, ncol = nSim)
width_A <- width_B <- width_1 <- width_2 <- width_Vn <- width_Sn <- rep(NA, nSim)

for (s in 1:nSim) {
  
  print(s)
  
  # Generate data
  set.seed(s)
  X_all <- matrix(rnorm(n_all*p), nrow = n_all, ncol = p)
  X <- X_all[1:n_obs,]
  X_new <- X_all[(n_obs+1):n_all,]
  mu_true <- c(X %*% beta_true)
  eps <- rt(n_obs, df = df_true)
  y <- mu_true + sigma_true*eps
  
  # RLS
  b_0 <- rep(0, p)
  V_0 <- diag(p)*1000
  bb <- rls(y, X, b_0, V_0)
  
  # Fixed covariates for predictive resamplings
  set.seed(s*3)
  X_ps <- matrix(nrow = NROW(X_new), ncol = p)
  for (i in 1:NROW(X_ps)) {
    idx <- sample(1:n_obs, size = 1)
    X_ps[i,] <- X[idx,]
  }
  
  # Sigma_n
  # Exact covariances under the covariate distribution used by predictive resampling
  cc <- Sigma_n_lm(bb$P, X_ps)
  # Remainder beyond these 1e5 covariates added through sigma^2 V (Section S5 of the supplement)
  cc$Sigma <- cc$Sigma + cc$V
  # Predictive resampling uses only the first N - n of them
  X_ps <- X_ps[1:N_ps, ]
  V_N_ps <- Sigma_n_lm(bb$P, X_ps)$V
  lowtri <- c(lower.tri(cc$Sigma, diag = TRUE))
  in_ci_beta_Sn <- rep(0, p)
  width_ci_beta_Sn <- rep(NA, p)
  for (j in 1:p) {
    ci_j <- qnorm(c(0.025, 0.975), mean = bb$beta[j], sd = sqrt(var_true)*sqrt(cc$Sigma[j,j]))
    width_ci_beta_Sn[j] <- abs(ci_j[2] - ci_j[1])
    if (ci_j[1] <= beta_true[j] & beta_true[j] <= ci_j[2]) in_ci_beta_Sn[j] <- 1
  }
  in_ci_Sn[,s] <- in_ci_beta_Sn
  width_Sn[s] <- mean(width_ci_beta_Sn)
  
  # V_n
  fr_Sn_Vn[s] <- sqrt(sum((c(cc$Sigma)[lowtri] - c(cc$V_n)[lowtri])^2))
  in_ci_beta_Vn <- rep(0, p)
  width_ci_beta_Vn <- rep(NA, p)
  for (j in 1:p) {
    ci_j <- qnorm(c(0.025, 0.975), mean = bb$beta[j], sd = sqrt(var_true)*sqrt(cc$V_n[j,j]))
    width_ci_beta_Vn[j] <- abs(ci_j[2] - ci_j[1])
    if (ci_j[1] <= beta_true[j] & beta_true[j] <= ci_j[2]) in_ci_beta_Vn[j] <- 1
  }
  in_ci_Vn[,s] <- in_ci_beta_Vn
  width_Vn[s] <- mean(width_ci_beta_Vn)
  
  
  # Monte Carlo reference: B_ps Gaussian draws with the exact posterior mean and
  # covariance, for sigma_n^2 = sigma^2 / a_n (V_n) and sigma_n = sigma (Sigma_n)
  set.seed(s*5)
  Z <- matrix(rnorm(p*B_ps), nrow = p, ncol = B_ps)
  ex_Vn <- t(c(bb$beta) + sqrt(var_true) * t(chol(cc$V_n)) %*% Z)
  ex_Sn <- t(c(bb$beta) + sqrt(var_true) * t(chol(cc$Sigma)) %*% Z)
  mse_ex_Vn[s] <- mean((bb$beta - colMeans(ex_Vn))^2)
  mse_ex_Sn[s] <- mean((bb$beta - colMeans(ex_Sn))^2)
  fr_Vn_ex[s] <- sqrt(sum((c(cc$V_n)[lowtri] - c(cov(ex_Vn) / (var_true))[lowtri])^2))
  fr_Sn_ex[s] <- sqrt(sum((c(cc$Sigma)[lowtri] - c(cov(ex_Sn) / (var_true))[lowtri])^2))
  # Fraction of posterior variance omitted by stopping at N (Section 3.2), against n / N
  frac_omit[s] <- mean(diag(V_N_ps) / diag(cc$V_n))
  beta_n_all[,s] <- bb$beta
  V_n_all[,,s] <- cc$V_n
  Sigma_n_all[,,s] <- cc$Sigma
  V_N_all[,,s] <- V_N_ps
  
  # 1 := PS standard with \sigma_n^2 = \sigma^2
  t_1 <- system.time(pp <- ps_lm_par_t(bb$beta, bb$P, X, N_ps, B_ps, sigma_true, df_true, nthreads = 7, correct = TRUE))
  time_1[s] <- t_1[3]
  mse_1[s] <- mean((bb$beta - colMeans(pp))^2)
  cov_1 <- cov(pp) / var_true
  fr_Sn_1[s] <- sqrt(sum((c(cc$Sigma)[lowtri] - c(cov_1)[lowtri])^2))
  fr_Vn_1[s] <- sqrt(sum((c(cc$V_n)[lowtri] - c(cov_1)[lowtri])^2))
  in_ci_beta_1 <- rep(0, p)
  width_ci_beta_1 <- rep(NA, p)
  for (j in 1:p) {
    ci_j <- robustBLME::hpd(pp[,j])
    width_ci_beta_1[j] <- abs(ci_j[2] - ci_j[1])
    if (ci_j[1] <= beta_true[j] & beta_true[j] <= ci_j[2]) in_ci_beta_1[j] <- 1
  }
  in_ci_1[,s] <- in_ci_beta_1
  width_1[s] <- mean(width_ci_beta_1)
  post_mean_1[,s] <- colMeans(pp)
  post_cov_1[,,s] <- cov(pp)
  post_sum_1[,,s] <- summarize_draws(t(pp))
  
  
  # 2 := PS fixed with \sigma_n^2 = \sigma^2
  t_2 <- system.time(aa <- ps_lm_fixed_par_t(bb$beta, bb$P, X_ps, B_ps, sigma_true, df_true, nthreads = 7, correct = TRUE))
  time_2[s] <- t_2[3]
  mse_2[s] <- mean((bb$beta - colMeans(aa))^2)
  cov_2 <- cov(aa) / var_true
  fr_Sn_2[s] <- sqrt(sum((c(cc$Sigma)[lowtri] - c(cov_2)[lowtri])^2))
  fr_Vn_2[s] <- sqrt(sum((c(cc$V_n)[lowtri] - c(cov_2)[lowtri])^2))
  in_ci_beta_2 <- rep(0, p)
  width_ci_beta_2 <- rep(NA, p)
  for (j in 1:p) {
    ci_j <- robustBLME::hpd(aa[,j])
    width_ci_beta_2[j] <- abs(ci_j[2] - ci_j[1])
    if (ci_j[1] <= beta_true[j] & beta_true[j] <= ci_j[2]) in_ci_beta_2[j] <- 1
  }
  in_ci_2[,s] <- in_ci_beta_2
  width_2[s] <- mean(width_ci_beta_2)
  post_mean_2[,s] <- colMeans(aa)
  post_cov_2[,,s] <- cov(aa)
  post_sum_2[,,s] <- summarize_draws(t(aa))
  
  
  # A := PS standard with \sigma_n^2 = \sigma^2 / a_n
  t_A <- system.time(vv <- ps_lm_par_t_alt(bb$beta, bb$P, X, N_ps, B_ps, sigma_true, df_true, nthreads = 7, keep_sigma = FALSE, correct = TRUE))
  time_A[s] <- t_A[3]
  mse_A[s] <- mean((bb$beta - colMeans(vv))^2)
  cov_A <- cov(vv) / var_true
  fr_Sn_A[s] <- sqrt(sum((c(cc$Sigma)[lowtri] - c(cov_A)[lowtri])^2))
  fr_Vn_A[s] <- sqrt(sum((c(cc$V_n)[lowtri] - c(cov_A)[lowtri])^2))
  in_ci_beta_A <- rep(0, p)
  width_ci_beta_A <- rep(NA, p)
  for (j in 1:p) {
    ci_j <- robustBLME::hpd(vv[,j])
    width_ci_beta_A[j] <- abs(ci_j[2] - ci_j[1])
    if (ci_j[1] <= beta_true[j] & beta_true[j] <= ci_j[2]) in_ci_beta_A[j] <- 1
  }
  in_ci_A[,s] <- in_ci_beta_A
  width_A[s] <- mean(width_ci_beta_A)
  post_mean_A[,s] <- colMeans(vv)
  post_cov_A[,,s] <- cov(vv)
  post_sum_A[,,s] <- summarize_draws(t(vv))
  
  
  # B := PS fixed with \sigma_n^2 = \sigma^2 / a_n
  t_B <- system.time(ll <- ps_lm_fixed_par_t_alt(bb$beta, bb$P, X_ps, B_ps, sigma_true, df_true, nthreads = 7, correct = TRUE))
  time_B[s] <- t_B[3]
  mse_B[s] <- mean((bb$beta - colMeans(ll))^2)
  cov_B <- cov(ll) / var_true
  fr_Sn_B[s] <- sqrt(sum((c(cc$Sigma)[lowtri] - c(cov_B)[lowtri])^2))
  fr_Vn_B[s] <- sqrt(sum((c(cc$V_n)[lowtri] - c(cov_B)[lowtri])^2))
  in_ci_beta_B <- rep(0, p)
  width_ci_beta_B <- rep(NA, p)
  for (j in 1:p) {
    ci_j <- robustBLME::hpd(ll[,j])
    width_ci_beta_B[j] <- abs(ci_j[2] - ci_j[1])
    if (ci_j[1] <= beta_true[j] & beta_true[j] <= ci_j[2]) in_ci_beta_B[j] <- 1
  }
  in_ci_B[,s] <- in_ci_beta_B
  width_B[s] <- mean(width_ci_beta_B)
  post_mean_B[,s] <- colMeans(ll)
  post_cov_B[,,s] <- cov(ll)
  post_sum_B[,,s] <- summarize_draws(t(ll))
  
  
  name_out <- file.path("output", paste0("sim_th_p", p, "_n", n_obs, "_t.RData"))
  save(fr_Sn_A, fr_Sn_B, fr_Sn_1, fr_Sn_2, fr_Sn_Vn, fr_Vn_A, fr_Vn_B, fr_Vn_1, fr_Vn_2,
       time_A, time_B, time_1, time_2, mse_A, mse_B, mse_1, mse_2,
       in_ci_A, in_ci_B, in_ci_1, in_ci_2, in_ci_Vn, in_ci_Sn,
       width_A, width_B, width_1, width_2, width_Vn, width_Sn,
       mse_ex_Vn, mse_ex_Sn, fr_Vn_ex, fr_Sn_ex, frac_omit,
       post_mean_1, post_mean_2, post_mean_A, post_mean_B, post_cov_1, post_cov_2, post_cov_A, post_cov_B,
       post_sum_1, post_sum_2, post_sum_A, post_sum_B, beta_n_all, V_n_all, Sigma_n_all, V_N_all, beta_true,
       file = name_out)
  
}


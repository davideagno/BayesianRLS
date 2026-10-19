
rm(list = ls())

library(mgcv)
library(coda)
library(brms)
library(parallel)

Rcpp::sourceCpp("RLS_functions.cpp")
source("RLS_functions.R")

# Preliminary settings
p <- 50
n_obs <- 10*p
n_all <- 1e5
sigma_true <- 2
nSim <- 30
df_true <- 5
var_true <- (sigma_true^2)*(df_true / (df_true-2))

B_ps <- 10000
N_ps <- n_all - n_obs

# beta_true, covariates and responses
set.seed(1997); beta_true <- rnorm(p, sd = 2)

fr_Sn_A <- fr_Sn_B <- fr_Sn_1 <- fr_Sn_2 <- fr_Sn_Vn <- rep(NA, nSim)
fr_Vn_A <- fr_Vn_B <- fr_Vn_1 <- fr_Vn_2 <- rep(NA, nSim)
time_A <- time_B <- time_1 <- time_2 <- rep(NA, nSim)
mse_A <- mse_B <- mse_1 <- mse_2 <- rep(NA, nSim)
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
  cc <- Sigma_n_lm(bb$P, X_new)
  lowtri <- c(lower.tri(cc$Sigma, diag = TRUE))
  in_ci_beta_Sn <- rep(0, p)
  width_ci_beta_Sn <- rep(NA, p)
  for (j in 1:p) {
    ci_j <- qnorm(c(0.025, 0.975), mean = bb$beta[j], sd = sigma_true*sqrt(cc$Sigma[j,j]))
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
    ci_j <- qnorm(c(0.025, 0.975), mean = bb$beta[j], sd = sigma_true*sqrt(cc$V_n[j,j]))
    width_ci_beta_Vn[j] <- abs(ci_j[2] - ci_j[1])
    if (ci_j[1] <= beta_true[j] & beta_true[j] <= ci_j[2]) in_ci_beta_Vn[j] <- 1
  }
  in_ci_Vn[,s] <- in_ci_beta_Vn
  width_Vn[s] <- mean(width_ci_beta_Vn)
  
  
  # 1 := PS standard with \sigma_n^2 = \sigma^2
  t_1 <- system.time(pp <- ps_lm_par_t(bb$beta, bb$P, X, N_ps, B_ps, sigma_true, df_true, nthreads = 7))
  time_1[s] <- t_1[3]
  mse_1[s] <- mean((beta_true - colMeans(pp))^2)
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
  
  
  # 2 := PS fixed with \sigma_n^2 = \sigma^2
  t_2 <- system.time(aa <- ps_lm_fixed_par_t(bb$beta, bb$P, X_ps, B_ps, sigma_true, df_true, nthreads = 7))
  time_2[s] <- t_2[3]
  mse_2[s] <- mean((beta_true - colMeans(aa))^2)
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
  
  
  # A := PS standard with \sigma_n^2 = \sigma^2 / a_n
  t_A <- system.time(vv <- ps_lm_par_t_alt(bb$beta, bb$P, X, N_ps, B_ps, sigma_true, df_true, nthreads = 7))
  time_A[s] <- t_A[3]
  mse_A[s] <- mean((beta_true - colMeans(vv))^2)
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
  
  
  # B := PS fixed with \sigma_n^2 = \sigma^2 / a_n
  t_B <- system.time(ll <- ps_lm_fixed_par_t_alt(bb$beta, bb$P, X_ps, B_ps, sigma_true, df_true, nthreads = 7))
  time_B[s] <- t_B[3]
  mse_B[s] <- mean((beta_true - colMeans(ll))^2)
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
  
  
  name_out <- paste0("sim_th_p", p, "_n", n_obs, "_t.RData")
  save(fr_Sn_A, fr_Sn_B, fr_Sn_1, fr_Sn_2, fr_Sn_Vn, fr_Vn_A, fr_Vn_B, fr_Vn_1, fr_Vn_2,
       time_A, time_B, time_1, time_2, mse_A, mse_B, mse_1, mse_2,
       in_ci_A, in_ci_B, in_ci_1, in_ci_2, in_ci_Vn, in_ci_Sn,
       width_A, width_B, width_1, width_2, width_Vn, width_Sn,
       file = name_out)
  
}


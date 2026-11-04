# Section 5.3: simulation with additive models (Table 2 and Table S3).
# Run scripts 1-6 of this folder in order, from this folder, with n_obs = 30*p and
# n_obs = 100*p. Step 1: simulate the training and test datasets.

rm(list = ls())

library(mgcv)
library(coda)
library(brms)

Rcpp::sourceCpp("../RLS_functions.cpp")
source("../RLS_functions.R")

# Preliminary settings
p <- 15
n_mult <- if (length(commandArgs(TRUE)) > 0) as.numeric(commandArgs(TRUE)[1]) else 30  # n = n_mult * p
n_obs <- n_mult*p
dir.create("output", showWarnings = FALSE)
nSim <- 30
n_test <- 1e4
df_val <- 5
sigma_val <- 1

# Generate data
dat_sim <- list()
for (s in 1:nSim) {
  dat_sim[[s]] <- generate_am_t(n = n_obs, p = p, df = df_val, sigma = sigma_val, seed = s)
}
dat_test_sim <- generate_am_t(n = n_test, p = p, df = df_val, seed = 2026*p)
name_file <- file.path("output", paste0("sim_p", p, "_n", n_obs, "_df", df_val, "_sigma", sigma_val, ".RData"))
save(dat_sim, dat_test_sim, file = name_file)


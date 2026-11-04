# Analysis of the NHANES serum triglycerides dataset (Section 6 of the paper and
# Sections S3-S4 of the supplementary material).
#
# Run from this folder (CODE/sec6_application/). Input: appl_data.RData (built by 1_download_data.R).
# Output: paper figures in ../../PAPER/figures (or output/), the trace plot and the
# posterior draws in output/.
#
# Steps:
#   1. data, train/test split and standardization
#   2. selection of the predictive distribution (Figures S3-S4)
#   3. spline design matrices
#   4. predictive Bayes: LOO-PLL, posterior mean, predictive resampling
#   5. traditional Bayesian additive model (brms) and frequentist additive models (mgcv)
#   6. out-of-sample MSE and computing times
#   7. figures

rm(list = ls())

Rcpp::sourceCpp("../RLS_functions.cpp")
source("../RLS_functions.R")

library(coda)
library(extraDistr)
library(brms)
library(mgcv)
library(ggplot2)
library(ggthemes)

out_dir <- "output"
dir.create(out_dir, showWarnings = FALSE)
# Paper figures are written directly to PAPER/figures when the script is run within the
# paper project, and to output/ otherwise
fig_dir <- if (dir.exists("../../PAPER/figures")) "../../PAPER/figures" else "output"

theme_appl <- theme_light(base_size = 9) +
  theme(
    strip.text = element_text(colour = "black"),
    strip.background = element_rect(fill = "gray82"),
    panel.grid.major = element_line(linewidth = 0.3, colour = "gray93"),
    panel.grid.minor = element_line(linewidth = 0.15, colour = "gray93")
  )


# 1. Data, train/test split and standardization --------------------------------

load("appl_data.RData")
contrasts(df$ethn) <- contr.treatment(levels(df$ethn), base = 3)

set.seed(2026)
size_train <- round(0.7 * NROW(df))
idx_train <- sort(sample(1:NROW(df), size = size_train, replace = FALSE))
idx_test <- sort(setdiff(1:NROW(df), idx_train))

# Continuous covariates: age, bmi, hdl, glucose, sbp, pir (waist is excluded)
cont_vars <- c(3:4, 6:9)

df_train <- df[idx_train, ]
y <- df_train$log_tg
covariates <- df_train[, cont_vars]
train_loc <- train_scale <- rep(NA, NCOL(covariates))
for (j in 1:NCOL(covariates)) {
  vec <- scale(covariates[, j])
  covariates[, j] <- vec
  train_loc[j] <- attr(vec, "scaled:center")
  train_scale[j] <- attr(vec, "scaled:scale")
}
mm <- model.matrix(~ethn, data = df_train)
covariates <- cbind(covariates, mm[, -1])
df_mgcv <- covariates
df_mgcv$y <- y

# Test set, standardized with the training location and scale
df_test_1 <- df[idx_test, ]
y_test <- df_test_1$log_tg
covariates_test <- df_test_1[, cont_vars]
for (j in 1:NCOL(covariates_test)) {
  covariates_test[, j] <- (covariates_test[, j] - train_loc[j]) / train_scale[j]
}
mm_test <- model.matrix(~ethn, data = df_test_1)
covariates_test <- cbind(covariates_test, mm_test[, -1])
df_test <- covariates_test
df_test$y <- y_test

am_formula <- y ~ s(age) + s(bmi) + s(hdl) + s(glucose) + s(sbp) + s(pir) +
  ethn1 + ethn2 + ethn4 + ethn6 + ethn7


# 2. Selection of the predictive distribution (supplement, Section S3) ---------

# Q-Q plot of the residuals of a Gaussian additive model (Figure S4)
fit_a <- mgcv::gam(am_formula, data = df_mgcv, method = "REML", family = gaussian())
res <- residuals(fit_a)
nus <- c(9, 10000)
n <- length(res)
df_qq <- do.call(rbind, lapply(nus, function(nu) {
  data.frame(nu = paste("df =", nu), teorici = qt(ppoints(n), df = nu), res = sort(res))
}))
# Reference line through the first and third quartiles, as in qqline()
df_line <- do.call(rbind, lapply(nus, function(nu) {
  y <- unname(quantile(res, c(0.25, 0.75)))
  x <- qt(c(0.25, 0.75), df = nu)
  slope <- diff(y) / diff(x)
  data.frame(nu = paste("df =", nu), slope = slope, intercept = y[1] - slope * x[1])
}))
livelli <- paste("df =", nus)
titoli <- c('"Student-t" ~ (nu == 9)', "Gaussian")
df_qq$nu <- factor(df_qq$nu, levels = livelli, labels = titoli)
df_line$nu <- factor(df_line$nu, levels = livelli, labels = titoli)

# Variance of the response within the quartiles of each covariate (Figure S3)
quartile_var <- function(x) {
  q <- quantile(x, c(0.25, 0.50, 0.75))
  c(
    var(y[x > 0 & x <= q[1]]), var(y[x > q[1] & x <= q[2]]),
    var(y[x > q[2] & x <= q[3]]), var(y[x > q[3] & x <= max(x)])
  )
}
v_age <- quartile_var(df_train$age)
v_bmi <- quartile_var(df_train$bmi)
v_hdl <- quartile_var(df_train$hdl)
v_glucose <- quartile_var(df_train$glucose)
v_sbp <- quartile_var(df_train$sbp)
v_pir <- quartile_var(df_train$pir)
v_pir[4] <- var(y[df_train$pir == 5.00]) # pir is top-coded at 5

variabili <- c("PIR", "BMI", "Age", "Glucose", "SBP", "HDL")
df_var <- data.frame(
  gruppo = rep(1:4, times = length(variabili)),
  varianza = c(v_pir, v_bmi, v_age, v_glucose, v_sbp, v_hdl),
  variabile = factor(rep(variabili, each = 4), levels = variabili)
)


# 3. Spline design matrices -----------------------------------------------------

# Columns of X: intercept, ethnicity dummies (unpenalized), then one thin plate
# regression spline basis per continuous covariate, with the centering constraint
# absorbed. Omega contains only the penalty matrices of the smooth terms.
n_disc <- NCOL(mm[, -1])
var_order <- c(6 + 1:n_disc, 1:6)
n_par <- 1 + n_disc

smooths <- lapply(var_order[-c(1:n_disc)], function(j) {
  spec <- eval(call("s", as.name(colnames(covariates)[j])))
  mgcv::smoothCon(spec, data = covariates, knots = NULL, absorb.cons = TRUE)[[1]]
})

Omega <- list()
X_splines <- rep(1, length(y))
for (j in 1:n_disc) {
  X_splines <- cbind(X_splines, covariates[, var_order[j]])
}
for (j in (n_disc + 1):length(var_order)) {
  X_splines <- cbind(X_splines, smooths[[j - n_disc]]$X)
  Omega[[j - n_disc]] <- smooths[[j - n_disc]]$S[[1]]
}
stopifnot(all(colSums(X_splines[, 2:n_par]) >= 2)) # each dummy needs at least two ones

X_splines_test <- rep(1, length(y_test))
for (j in 1:n_disc) {
  X_splines_test <- cbind(X_splines_test, covariates_test[, var_order[j]])
}
for (j in (n_disc + 1):length(var_order)) {
  X_splines_test <- cbind(X_splines_test, mgcv::PredictMat(smooths[[j - n_disc]], df_test))
}

b_0 <- rep(0, NCOL(X_splines))


# 4. Predictive Bayes -----------------------------------------------------------

# Skew Student-t P_n (supplement, Section S3): the skewness estimate is close to zero
opt_skt <- nlminb(
  start = c(log(5), log(1), 0, rep(0.5, length(Omega))),
  objective = npll_loo_alt_appl_skt,
  y = y, X = X_splines, beta_0 = b_0, Omega = Omega, n_par = n_par
)
c(df = exp(opt_skt$par[1]), sigma = exp(opt_skt$par[2]), alpha = opt_skt$par[3])

# Student-t P_n: hyperparameters selected by maximizing the LOO-PLL
t_opt <- system.time(opt <- nlminb(
  start = c(log(5), log(1), rep(0.5, length(Omega))),
  objective = npll_loo_alt_appl,
  y = y, X = X_splines, beta_0 = b_0, Omega = Omega, n_par = n_par
))
df_opt <- exp(opt$par[1])
sigma_opt <- exp(opt$par[2])
lambda_opt <- exp(opt$par[-c(1:2)])
c(df = df_opt, sigma = sigma_opt)

# Posterior mean
t_est <- system.time(bb <- rls_lambda_appl(y, X_splines, b_0, lambda_opt, Omega, n_par))

# Convergence of a single forward path (Figure S8)
set.seed(123)
c_ps <- ps_rs_conv_t_1_alt(X_splines, bb$beta, bb$P, 4e5, 1, df_opt, sigma_opt)

# Predictive resampling with fixed future covariates (Algorithm 2)
t_ps <- system.time(beta_ps <- ps_rs_parallel_t_fixed_alt(
  X_splines, bb$beta, bb$P,
  N = 3e5, B = 10000, df = df_opt, sigma = sigma_opt, seed = 1, nthreads = 7
))
beta_ps <- as.mcmc(beta_ps)


# 5. Competitors ----------------------------------------------------------------

# Traditional Bayesian additive model with Student-t likelihood (brms). The prior on
# the coefficients is Gaussian with precision blockdiag(lambda_j * Omega_j) / sigma^2,
# flat on the unpenalized columns; block j occupies b[strt[j]:(strt[j] + d[j] - 1)].
P <- length(Omega)
d <- vapply(Omega, NCOL, integer(1))
Dmax <- max(d)
strt <- as.integer(n_par + 1L + c(0L, head(cumsum(d), -1)))
stopifnot(strt[P] + d[P] - 1 == NCOL(X_splines))
ranks <- vapply(seq_len(P), function(j) {
  ev <- eigen(Omega[[j]], symmetric = TRUE, only.values = TRUE)$values
  sum(ev > max(ev) * 1e-10)
}, integer(1))
dat_splines <- data.frame(y = y, X_splines)
colnames(dat_splines) <- c("y", paste0("X", seq_len(ncol(X_splines))))
Omega_arr <- array(0, dim = c(P, Dmax, Dmax))
for (j in seq_len(P)) {
  Om <- as.matrix(Omega[[j]])
  Omega_arr[j, seq_len(d[j]), seq_len(d[j])] <- (Om + t(Om)) / 2
}
# Gamma(a, b) prior on the smoothing parameters
blm_stanvars <- function(a, b) {
  stanvar(as.integer(P), name = "P", scode = "  int<lower=1> P;") +
    stanvar(as.integer(Dmax), name = "Dmax", scode = "  int<lower=1> Dmax;") +
    stanvar(as.integer(d), name = "dblk", scode = "  array[P] int<lower=1> dblk;") +
    stanvar(strt, name = "strt", scode = "  array[P] int<lower=1> strt;") +
    stanvar(as.integer(ranks), name = "rk", scode = "  array[P] int<lower=1> rk;") +
    stanvar(Omega_arr, name = "Omega", scode = "  array[P] matrix[Dmax, Dmax] Omega;") +
    stanvar(
      scode = "  vector<lower=0>[P] lambda;  // smoothing parameters",
      block = "parameters"
    ) +
    stanvar(scode = paste0("
  target += gamma_lpdf(lambda | ", a, ", ", b, ");
  // beta ~ N(0, sigma^2 * S^{-1}),  S = blockdiag(lambda_j * Omega_j)
  // blocks may differ in size: block j occupies b[strt[j] : strt[j]+dblk[j]-1]
  for (j in 1:P) {
    vector[dblk[j]] beta_j = segment(b, strt[j], dblk[j]);
    matrix[dblk[j], dblk[j]] Om_j = block(Omega[j], 1, 1, dblk[j], dblk[j]);
    target += -0.5 * (lambda[j] / square(sigma)) * quad_form_sym(Om_j, beta_j)
              + 0.5 * rk[j] * log(lambda[j]) - rk[j] * log(sigma);
  }"), block = "model")
}
fit_blm <- function(stanvars) {
  brms::brm(
    formula = y ~ -1 + ., data = dat_splines, family = student(),
    stanvars = stanvars, prior = prior("", class = "b"),
    chains = 5, cores = 5, iter = 20000, warmup = 18000, seed = 1
  )
}

# Gamma(2, 0.5) prior, as in the simulation study of Section 5.3
t_blm <- system.time(blm <- fit_blm(blm_stanvars(2, 0.5)))
post_samp <- coda::as.mcmc(as.matrix(blm))
ess_blm <- mean(coda::effectiveSize(post_samp[, 1:NCOL(X_splines)]))

# Sensitivity to a vague Gamma(0.001, 0.001) prior on the smoothing parameters
t_blm_vague <- system.time(blm_vague <- fit_blm(blm_stanvars(0.001, 0.001)))
post_vague <- coda::as.mcmc(as.matrix(blm_vague))
ess_blm_vague <- mean(coda::effectiveSize(post_vague[, 1:NCOL(X_splines)]))

# Frequentist additive models (mgcv)
fit_t <- mgcv::gam(am_formula, data = df_mgcv, method = "REML", family = scat())
fit_n <- mgcv::gam(am_formula, data = df_mgcv, method = "REML", family = gaussian())


# 6. Out-of-sample MSE and computing times --------------------------------------

pm_blm <- colMeans(post_samp[, 1:NCOL(X_splines)])
mse_rls <- mean((df_test$y - c(X_splines_test %*% bb$beta))^2)
mse_blm <- mean((df_test$y - c(X_splines_test %*% pm_blm))^2)
mse_blm_vague <- mean((df_test$y - c(X_splines_test %*% colMeans(post_vague[, 1:NCOL(X_splines)])))^2)
mse_am_t <- mean((df_test$y - c(X_splines_test %*% fit_t$coefficients))^2)
mse_am_n <- mean((df_test$y - c(X_splines_test %*% fit_n$coefficients))^2)
c(rls = mse_rls, blm = mse_blm, blm_vague = mse_blm_vague, am_t = mse_am_t, am_n = mse_am_n)

c(
  rls_pm = t_opt[[3]] + t_est[[3]], rls_ps = t_ps[[3]],
  rls_tot = t_opt[[3]] + t_est[[3]] + t_ps[[3]], blm = t_blm[[3]], blm_vague = t_blm_vague[[3]]
)
c(ess_blm = ess_blm, ess_blm_vague = ess_blm_vague)
c(divergent = rstan::get_num_divergent(blm$fit), divergent_vague = rstan::get_num_divergent(blm_vague$fit))

save(beta_ps, post_samp, post_vague, c_ps, opt, opt_skt, bb, t_opt, t_est, t_ps, t_blm, t_blm_vague,
  file = file.path(out_dir, "appl_results.RData")
)


# 7. Figures --------------------------------------------------------------------

lev_v <- c("Age", "BMI", "HDL", "Glucose", "SBP", "PIR") # panel order
var_col <- c(Age = "age", BMI = "bmi", HDL = "hdl", Glucose = "glucose", SBP = "sbp", PIR = "pir")
sc_idx <- c(Age = 1, BMI = 2, HDL = 3, Glucose = 4, SBP = 5, PIR = 6)

# Partial residuals for the predictive Bayes panels: Gaussian additive model with the
# smoothing parameters selected by the LOO-PLL, whose fit coincides with bb$beta
fit_pl <- mgcv::gam(am_formula, data = df_mgcv, sp = pmin(lambda_opt, 1e6))
bb_chk <- rls_lambda_appl(y, X_splines, b_0, pmin(lambda_opt, 1e6), Omega, n_par = n_par)
max(abs(coef(fit_pl) - bb_chk$beta))
tmp <- tempfile()
png(tmp)
pd <- plot(fit_pl, residuals = TRUE)
dev.off()
unlink(tmp)
pres <- setNames(lapply(1:6, function(k) pd[[k]]$p.resid), lev_v)

# Evaluation grids (standardized scale) and corresponding spline bases
n_plot <- 500
grid <- lapply(var_col, function(v) seq(min(df_mgcv[[v]]), max(df_mgcv[[v]]), length = n_plot))
X_plot <- as.data.frame(grid)
names(X_plot) <- var_col
X_splines_plot <- cbind(1, do.call(cbind, lapply(smooths, function(sm) mgcv::PredictMat(sm, X_plot))))

# Columns of the coefficient vector for each smooth term (9 basis functions each)
idx <- setNames(lapply(0:5, function(k) n_par + 9 * k + 1:9), lev_v)

# Pointwise HPD bands from the posterior draws of each effect on the grid
hpd_rows <- function(draws) t(sapply(seq_len(NROW(draws)), function(i) robustBLME::hpd(draws[i, ])))
effects <- lapply(lev_v, function(v) {
  new_splines <- X_splines_plot[, idx[[v]] - n_disc]
  draws_ps <- tcrossprod(new_splines, beta_ps[, idx[[v]]])
  draws_blm <- tcrossprod(new_splines, post_samp[, idx[[v]]])
  list(
    new_splines = new_splines,
    ci_ps = hpd_rows(draws_ps),
    pm_blm = apply(draws_blm, 1, mean),
    ci_blm = hpd_rows(draws_blm)
  )
})
names(effects) <- lev_v

# Covariate values on the original scale
to_orig <- function(x, v) x * train_scale[sc_idx[v]] + train_loc[sc_idx[v]]

# Clip the bands to the range of the partial residuals of each covariate
clip_bands <- function(df_fit, df_pts) {
  rng <- do.call(rbind, lapply(lev_v, function(v) {
    r <- range(df_pts$pres[df_pts$var == v])
    data.frame(var = factor(v, levels = lev_v), ymin = r[1], ymax = r[2])
  }))
  i <- match(df_fit$var, rng$var)
  df_fit$low <- pmax(df_fit$low, rng$ymin[i])
  df_fit$upp <- pmin(df_fit$upp, rng$ymax[i])
  df_fit
}

# One colour per method, for both the line and its band: predictive Bayes in orange and
# traditional Bayes in blue (as in Figure S5), mgcv in grey
partial_effects_plot <- function(df_fit, df_pts, colour) {
  ggplot(df_fit, aes(x = x)) +
    geom_ribbon(aes(ymin = low, ymax = upp), fill = colour, alpha = 0.25) +
    geom_line(aes(y = fit), colour = colour, linewidth = 0.7) +
    geom_point(aes(y = pres), data = df_pts, size = 0.3, colour = "gray35", alpha = 0.6) +
    facet_wrap(~var, nrow = 2, ncol = 3, scales = "free") +
    labs(x = NULL, y = "Partial effect") +
    theme_appl
}

# Predictive Bayes (Figure 2)
df_fit_rls <- do.call(rbind, lapply(lev_v, function(v) {
  e <- effects[[v]]
  data.frame(
    x = to_orig(grid[[v]], v), fit = as.vector(e$new_splines %*% bb$beta[idx[[v]]]),
    low = e$ci_ps[, 1], upp = e$ci_ps[, 2], var = v
  )
}))
df_pts_rls <- do.call(rbind, lapply(lev_v, function(v) {
  data.frame(x = to_orig(df_mgcv[[var_col[v]]], v), pres = as.vector(pres[[v]]), var = v)
}))
df_fit_rls$var <- factor(df_fit_rls$var, levels = lev_v)
df_pts_rls$var <- factor(df_pts_rls$var, levels = lev_v)
df_fit_rls <- clip_bands(df_fit_rls, df_pts_rls)
pl_rls <- partial_effects_plot(df_fit_rls, df_pts_rls, "#fc7d0b")

# Traditional Bayes (Figure S6): posterior mean of each effect and partial residuals
# computed from the posterior mean of the coefficients
df_fit_blm <- do.call(rbind, lapply(lev_v, function(v) {
  e <- effects[[v]]
  data.frame(
    x = to_orig(grid[[v]], v), fit = e$pm_blm,
    low = e$ci_blm[, 1], upp = e$ci_blm[, 2], var = v
  )
}))
df_pts_blm <- do.call(rbind, lapply(lev_v, function(v) {
  pres_blm <- y - X_splines[, -idx[[v]]] %*% pm_blm[-idx[[v]]]
  data.frame(x = to_orig(df_mgcv[[var_col[v]]], v), pres = as.vector(pres_blm), var = v)
}))
df_fit_blm$var <- factor(df_fit_blm$var, levels = lev_v)
df_pts_blm$var <- factor(df_pts_blm$var, levels = lev_v)
df_fit_blm <- clip_bands(df_fit_blm, df_pts_blm)
pl_blm <- partial_effects_plot(df_fit_blm, df_pts_blm, "#1170aa")

# Frequentist additive model with Student-t likelihood (Figure S7): pointwise 95%
# confidence bands; partial residuals are response residuals plus the smooth term
z <- qnorm(0.975)
term_lab <- setNames(paste0("s(", var_col, ")"), lev_v)
df_fit_gam <- do.call(rbind, lapply(lev_v, function(v) {
  g <- grid[[v]]
  nd <- df_mgcv[rep(1, length(g)), , drop = FALSE]
  nd[[var_col[v]]] <- g
  pr <- predict(fit_t, newdata = nd, type = "terms", terms = term_lab[v], se.fit = TRUE)
  fit <- as.vector(pr$fit[, term_lab[v]])
  se <- as.vector(pr$se.fit[, term_lab[v]])
  data.frame(x = to_orig(g, v), fit = fit, low = fit - z * se, upp = fit + z * se, var = v)
}))
pt <- predict(fit_t, type = "terms")
r_resp <- residuals(fit_t, type = "response")
df_pts_gam <- do.call(rbind, lapply(lev_v, function(v) {
  data.frame(
    x = to_orig(df_mgcv[[var_col[v]]], v),
    pres = as.vector(r_resp + pt[, term_lab[v]]), var = v
  )
}))
df_fit_gam$var <- factor(df_fit_gam$var, levels = lev_v)
df_pts_gam$var <- factor(df_pts_gam$var, levels = lev_v)
df_fit_gam <- clip_bands(df_fit_gam, df_pts_gam)
pl_gam <- partial_effects_plot(df_fit_gam, df_pts_gam, "#57606c")

# Posterior densities of the intercept and of the ethnicity coefficients (Figure S5).
# Columns 2-6 are ethn1, ethn2, ethn4, ethn6 and ethn7 (codes of RIDRETH3); the baseline
# is non-Hispanic White (code 3).
titoli_d <- c(
  "Intercept", "Mexican American", "Other Hispanic", "Non-Hispanic Black",
  "Non-Hispanic Asian", "Other or multiracial"
)
df_d <- do.call(rbind, lapply(1:6, function(j) {
  data.frame(rls = as.numeric(beta_ps[, j]), blm = as.numeric(post_samp[, j]), title = titoli_d[j])
}))
df_d$title <- factor(df_d$title, levels = titoli_d)
pl_d <- ggplot(df_d) +
  facet_wrap(~title, nrow = 2, ncol = 3, scales = "free") +
  stat_density(aes(x = blm, col = "blm"), geom = "line", position = "identity", adjust = 1.5) +
  stat_density(aes(x = rls, col = "rls"), geom = "line", position = "identity", adjust = 1.5) +
  scale_color_tableau(
    name = "", palette = "Color Blind",
    labels = c("rls" = "Predictive Bayes", "blm" = "Traditional Bayes")
  ) +
  xlab("") +
  ylab("Density") +
  theme_appl +
  theme(
    legend.text.align = 0,
    axis.text.y = element_blank(),
    plot.title = element_text(hjust = 0.5),
    legend.position = "top"
  )

# Trace plots of a selected coefficient (diagnostic, saved in output/)
jdx <- 30
lab <- c("Traditional Bayes", "Predictive Bayes")
s1 <- as.numeric(post_samp[, jdx])
s2 <- as.numeric(beta_ps[, jdx])
df_tr <- rbind(
  data.frame(iter = seq_along(s1), value = s1, chain = lab[1]),
  data.frame(iter = seq_along(s2), value = s2, chain = lab[2])
)
df_tr$chain <- factor(df_tr$chain, levels = lab)
df_vl <- data.frame(chain = factor(lab[1], levels = lab), x = c(2, 4, 6, 8) * 1e3) # chain boundaries
pl_trace <- ggplot(df_tr, aes(x = iter, y = value)) +
  geom_line(linewidth = 0.3) +
  geom_vline(aes(xintercept = x), data = df_vl, linetype = "dashed", colour = "gray40", linewidth = 0.2) +
  facet_wrap(~chain, nrow = 1, scales = "free_x") +
  coord_cartesian(ylim = range(s1)) + # same y-limits in both panels
  labs(x = "Iteration", y = "Regression Coefficient") +
  theme_appl

# Selection of the predictive distribution (Figures S3-S4)
pl_qq <- ggplot(df_qq, aes(x = teorici, y = res)) +
  geom_point() +
  geom_abline(data = df_line, aes(slope = slope, intercept = intercept), colour = "red") +
  facet_wrap(~nu, scales = "free_x", labeller = label_parsed) +
  labs(x = "Theoretical quantiles", y = "Residuals") +
  theme_appl +
  theme(legend.position = "top")

pl_var <- ggplot(df_var, aes(x = factor(gruppo), y = varianza, colour = variabile, group = variabile)) +
  geom_line() +
  geom_point(size = 2) +
  scale_color_tableau(palette = "Tableau 10") +
  guides(colour = guide_legend(nrow = 1)) +
  scale_x_discrete(labels = paste0("Q", 1:4)) +
  coord_cartesian(ylim = c(0, 0.8)) +
  labs(x = "Covariates quartile", y = "Response variance", colour = "Covariate") +
  theme_appl +
  theme(legend.position = "top")

# Convergence of a single forward path (Figure S8). The path is plotted at every step for
# the first 1,000 forward samples and every 40 steps afterwards, to keep the pdf small.
keep <- unique(c(1:1000, seq(1000, NROW(c_ps), by = 40)))
df_conv <- data.frame(
  conv = c(c_ps[keep, ]), fs = rep(keep, NCOL(c_ps)),
  group = rep(1:NCOL(c_ps), each = length(keep)), title = "Convergence"
)
pl_conv <- ggplot(df_conv) +
  facet_wrap(~title) +
  geom_line(aes(x = fs, y = conv, group = group), alpha = 0.7) +
  labs(x = "Number of forward samples", y = expression(L[1] ~ "distance")) +
  ylim(c(0, max(c_ps))) +
  scale_x_continuous(labels = scales::label_comma()) +
  theme_appl +
  theme(
    legend.text.align = 0,
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5),
    legend.position = "top"
  )

ggsave(file.path(out_dir, "trace_plot.pdf"), pl_trace, width = 6.3, height = 2.4, device = cairo_pdf)
ggsave(file.path(fig_dir, "pl_rls.pdf"), pl_rls, width = 6.3, height = 3.8, device = cairo_pdf)
ggsave(file.path(fig_dir, "pl_d.pdf"), pl_d, width = 6.3, height = 3.4, device = cairo_pdf)
ggsave(file.path(fig_dir, "pl_blm.pdf"), pl_blm, width = 6.3, height = 3.8, device = cairo_pdf)
ggsave(file.path(fig_dir, "pl_mgcv.pdf"), pl_gam, width = 6.3, height = 3.8, device = cairo_pdf)
ggsave(file.path(fig_dir, "pl_appl_conv.pdf"), pl_conv, width = 6.3, height = 3, device = cairo_pdf)
ggsave(file.path(fig_dir, "pl_appl_var.pdf"), pl_var, width = 6.3, height = 3.2, device = cairo_pdf)
ggsave(file.path(fig_dir, "pl_appl_qq.pdf"), pl_qq, width = 6.3, height = 3, device = cairo_pdf)

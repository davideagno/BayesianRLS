
rm(list = ls())

Rcpp::sourceCpp("RLS_functions.cpp")
source("RLS_functions.R")

library(coda)
library(extraDistr)
library(brms)
library(mgcv)
library(ggplot2)

# Dataset
load("appl_data.RData")
round(cor(df[,c(3:9)]), 2)
set.seed(1)
size_train <- round(0.5*NROW(df))
idx_train <- sort(sample(1:NROW(df), size = size_train, replace = FALSE))
idx_test <- sort(setdiff(1:NROW(df), idx_train))

# Train set
df_train <- df[idx_train,]
y <- df_train$log_tg
covariates <- df_train[,c(3:4,6:9)]
train_loc <- train_scale <- rep(NA, NCOL(covariates))
for (j in 1:NCOL(covariates)) {
  vec <- scale(covariates[,j])
  covariates[,j] <- vec
  train_loc[j] <- attr(vec, "scaled:center")
  train_scale[j] <- attr(vec, "scaled:scale")
}
mm <- model.matrix(~ ethn, data = df_train)
covariates <- cbind(covariates, mm[,-1])
df_mgcv <- covariates
df_mgcv$y <- y

# Test set
df_test_1 <- df[idx_test,]
y_test <- df_test_1$log_tg
covariates_test <- df_test_1[,c(3:4,6:9)]
for (j in 1:NCOL(covariates_test)) {
  covariates_test[,j] <- (covariates_test[,j] - train_loc[j]) / train_scale[j]
}
mm_test <- model.matrix(~ ethn, data = df_test_1)
covariates_test <- cbind(covariates_test, mm_test[,-1])
df_test <- covariates_test
df_test$y <- y_test


# Building \Lambda and \beta_0
var_order <- c(7:11, 1:6)
n_disc <- NCOL(mm[,-1])
colnames(covariates)[var_order]
smooths <- lapply(var_order[-c(1:n_disc)], function(j) {
  spec <- eval(call("s", as.name(colnames(covariates)[j])))   # builds s(xj)
  mgcv::smoothCon(spec, data = covariates, knots = NULL, absorb.cons = TRUE)[[1]]
})
Omega <- Z <- list(); X_splines <- rep(1, length(y))
for (j in 1:n_disc) {
  Z[[j]] <- covariates[,var_order[j]]
  Omega[[j]] <- matrix(100, 1, 1)
  X_splines <- cbind(X_splines, Z[[j]])
}
for (j in (n_disc+1):length(var_order)) {
  Z[[j]] <- smooths[[j-n_disc]]$X
  Omega[[j]] <- smooths[[j-n_disc]]$S[[1]]
  X_splines <- cbind(X_splines, Z[[j]])
}
X_splines_test <- rep(1, length(y_test))
for (j in 1:n_disc) {
  X_splines_test <- cbind(X_splines_test, covariates_test[,var_order[j]])
}
for (j in (n_disc+1):length(var_order)) {
  X_splines_test <- cbind(X_splines_test, mgcv::PredictMat(smooths[[j-n_disc]], df_test))
}
b_0 <- rep(0, NCOL(X_splines))

# LOO-PLL optimization
t_opt <- system.time(opt <- nlminb(start = c(log(2), log(1), rep(0.5, length(var_order))), objective = npll_loo_alt,
                                   y = y, X = X_splines, beta_0 = b_0, Omega = Omega))
t_opt
opt
(df_opt <- exp(opt$par[1]))
(sigma_opt <- exp(opt$par[2]))
lambda_opt <- exp(opt$par[-c(1:2)])

# Posterior mean
bb <- rls_lambda(y, X_splines, b_0, lambda_opt, Omega)

# Convergence of predictive resampling
set.seed(1234)
c_ps <- ps_rs_conv_t_1_alt(X_splines, bb$beta, bb$P, 4e5, 1, df_opt, sigma_opt)
plot(c_ps, type = "l")

# Predictive resampling
beta_ps <- ps_rs_parallel_t_fixed_alt(X_splines, bb$beta, bb$P, N = 3e5, B = 10000,
                                      df = df_opt, sigma = sigma_opt, seed = 1, nthreads = 7)
beta_ps <- as.mcmc(beta_ps)
# save(beta_ps, file = "ps_appl_10k.RData")
# load("ps_appl_10k.RData")



# Traditional Bayesian additive model with Student-t likelihood
colnames(covariates)[var_order]
P <- length(var_order)
d <- vapply(Z, NCOL, integer(1))
stopifnot(length(d) == P, length(Omega) == P)
Dmax <- max(d)
# start index of each block within the coefficient vector b
strt <- as.integer(head(c(1L, cumsum(d) + 1L), P))
# rank of each penalty, using a scale-relative tolerance
ranks <- vapply(seq_len(P), function(j) {
  ev <- eigen(Omega[[j]], symmetric = TRUE, only.values = TRUE)$values
  sum(ev > max(ev) * 1e-10)
}, integer(1))
dat_splines <- data.frame(y = y, X_splines)
colnames(dat_splines) <- c("y", paste0("X", seq_len(ncol(X_splines))))
# Pack penalties into Omega_arr[j, , ], zero-padded out to Dmax
Omega_arr <- array(0, dim = c(P, Dmax, Dmax))
for (j in seq_len(P)) {
  Om <- as.matrix(Omega[[j]])
  Om <- (Om + t(Om)) / 2                  # enforce exact symmetry for quad_form_sym
  Omega_arr[j, seq_len(d[j]), seq_len(d[j])] <- Om
}
stanvars <-
  stanvar(as.integer(P),      name = "P",     scode = "  int<lower=1> P;") +
  stanvar(as.integer(Dmax),   name = "Dmax",  scode = "  int<lower=1> Dmax;") +
  stanvar(as.integer(d),      name = "dblk",  scode = "  array[P] int<lower=1> dblk;") +
  stanvar(strt,               name = "strt",  scode = "  array[P] int<lower=1> strt;") +
  stanvar(as.integer(ranks),  name = "rk",    scode = "  array[P] int<lower=1> rk;") +
  stanvar(Omega_arr,          name = "Omega", scode = "  array[P] matrix[Dmax, Dmax] Omega;") +
  stanvar(scode = "  vector<lower=0>[P] lambda;  // smoothing parameters",
          block = "parameters") +
  stanvar(scode = "
  // weakly-informative prior on the smoothing parameters
  target += gamma_lpdf(lambda | 2, 0.5);
  // beta ~ N(0, sigma^2 * S^{-1}),  S = blockdiag(lambda_j * Omega_j)
  // blocks may differ in size: block j occupies b[strt[j] : strt[j]+dblk[j]-1]
  for (j in 1:P) {
    vector[dblk[j]] beta_j = segment(b, strt[j], dblk[j]);
    matrix[dblk[j], dblk[j]] Om_j = block(Omega[j], 1, 1, dblk[j], dblk[j]);
    target += -0.5 * (lambda[j] / square(sigma)) * quad_form_sym(Om_j, beta_j)
              + 0.5 * rk[j] * log(lambda[j]) - rk[j] * log(sigma);
  }",
          block = "model")
blm <- brms::brm(formula = y ~ -1+., data = dat_splines, family = student(),
                 stanvars = stanvars, prior = prior("", class = "b"),
                 chains = 5, cores = 5, iter = 20000, warmup = 18000, seed = 1)
post_samp <- coda::as.mcmc(as.matrix(blm))
# save(blm, post_samp, file = "blm_appl_10k.RData")
# load("blm_appl_10k.RData")
coda::effectiveSize(post_samp)


# Frequentist Student-t additive model
fit_t <- mgcv::gam(y ~ s(age) + s(bmi) + s(hdl) + s(glucose) + s(sbp) + s(pir) +
                     ethn1 + ethn2 + ethn4 + ethn6 + ethn7,
                   data = df_mgcv, method = "REML", family = scat())


# Out of sample predictions for f(x)
pred_rls_test <- c(X_splines_test %*% bb$beta)
mse_rls <- mean((df_test$y - pred_rls_test)^2)
pm_blm <- colMeans(post_samp[,1:NCOL(X_splines)])
pred_blm_test <- c(X_splines_test %*% pm_blm)
mse_blm <- mean((df_test$y - pred_blm_test)^2)
pred_am_t_test <- c(X_splines_test %*% fit_t$coefficients)
mse_am_t <- mean((df_test$y - pred_am_t_test)^2)
mse_rls; mse_am_t; mse_blm










##########################################################################################
#                                        FIGURES                                         #
##########################################################################################

# Partial residuals for RLS-t
fit_pl <- mgcv::gam(y ~ s(age) + s(bmi) + s(hdl) + s(glucose) + s(sbp) + s(pir) +
                      ethn1 + ethn2 + ethn4 + ethn6 + ethn7,
                    data = df_mgcv, method = "REML", sp = lambda_opt,
                    paraPen = list(ethn1 = list(Omega[[1]]), ethn2 = list(Omega[[2]]),
                                   ethn4 = list(Omega[[3]]), ethn6 = list(Omega[[4]]),
                                   ethn7 = list(Omega[[5]])))
cbind(c(fit_pl$coefficients), c(bb$beta))
summary(fit_pl)
tmp <- tempfile(); png(tmp)
pd <- plot(fit_pl, residuals = TRUE)
dev.off(); unlink(tmp)

pres_age <- pd[[1]]$p.resid
pres_bmi <- pd[[2]]$p.resid
pres_hdl <- pd[[3]]$p.resid
pres_gluc <- pd[[4]]$p.resid
pres_sbp <- pd[[5]]$p.resid
pres_pir <- pd[[6]]$p.resid

n_plot <- 120
age_new <- seq(min(df_mgcv$age), max(df_mgcv$age), length = n_plot)
bmi_new <- seq(min(df_mgcv$bmi), max(df_mgcv$bmi), length = n_plot)
hdl_new <- seq(min(df_mgcv$hdl), max(df_mgcv$hdl), length = n_plot)
gluc_new <- seq(min(df_mgcv$gluc), max(df_mgcv$gluc), length = n_plot)
sbp_new <- seq(min(df_mgcv$sbp), max(df_mgcv$sbp), length = n_plot)
pir_new <- seq(min(df_mgcv$pir), max(df_mgcv$pir), length = n_plot)
X_plot <- as.data.frame(cbind(age = age_new, bmi = bmi_new, hdl = hdl_new, 
                              glucose = gluc_new, sbp = sbp_new, pir = pir_new))
X_splines_plot <- cbind(1, do.call(cbind, lapply(smooths, function(sm) mgcv::PredictMat(sm, X_plot))))


idx_age <- 7:15
age_new_splines <- X_splines_plot[,idx_age-n_disc]
age_splines_ps <- tcrossprod(age_new_splines, beta_ps[,idx_age])
age_splines_blm <- tcrossprod(age_new_splines, post_samp[,idx_age])
age_low <- age_upp <- rep(NA, length = NROW(age_new_splines))
age_low_blm <- age_upp_blm <- age_pm_blm <- rep(NA, length = NROW(age_new_splines))
for (i in 1:NROW(age_new_splines)) {
  ci <- robustBLME::hpd(age_splines_ps[i,])
  age_low[i] <- ci[1]
  age_upp[i] <- ci[2]
  age_pm_blm[i] <- mean(age_splines_blm[i,])
  ci_blm <- robustBLME::hpd(age_splines_blm[i,])
  age_low_blm[i] <- ci_blm[1]
  age_upp_blm[i] <- ci_blm[2]
}

idx_bmi <- 16:24
bmi_new_splines <- X_splines_plot[,idx_bmi-n_disc]
bmi_splines_ps <- tcrossprod(bmi_new_splines, beta_ps[,idx_bmi])
bmi_splines_blm <- tcrossprod(bmi_new_splines, post_samp[,idx_bmi])
bmi_low <- bmi_upp <- rep(NA, length = NROW(bmi_new_splines))
bmi_low_blm <- bmi_upp_blm <- bmi_pm_blm <- rep(NA, length = NROW(bmi_new_splines))
for (i in 1:NROW(bmi_new_splines)) {
  ci <- robustBLME::hpd(bmi_splines_ps[i,])
  bmi_low[i] <- ci[1]
  bmi_upp[i] <- ci[2]
  bmi_pm_blm[i] <- mean(bmi_splines_blm[i,])
  ci_blm <- robustBLME::hpd(bmi_splines_blm[i,])
  bmi_low_blm[i] <- ci_blm[1]
  bmi_upp_blm[i] <- ci_blm[2]
}

idx_hdl <- 25:33
hdl_new_splines <- X_splines_plot[,idx_hdl-n_disc]
hdl_splines_ps <- tcrossprod(hdl_new_splines, beta_ps[,idx_hdl])
hdl_splines_blm <- tcrossprod(hdl_new_splines, post_samp[,idx_hdl])
hdl_low <- hdl_upp <- rep(NA, length = NROW(hdl_new_splines))
hdl_low_blm <- hdl_upp_blm <- hdl_pm_blm <- rep(NA, length = NROW(hdl_new_splines))
for (i in 1:NROW(hdl_new_splines)) {
  ci <- robustBLME::hpd(hdl_splines_ps[i,])
  hdl_low[i] <- ci[1]
  hdl_upp[i] <- ci[2]
  hdl_pm_blm[i] <- mean(hdl_splines_blm[i,])
  ci_blm <- robustBLME::hpd(hdl_splines_blm[i,])
  hdl_low_blm[i] <- ci_blm[1]
  hdl_upp_blm[i] <- ci_blm[2]
}

idx_gluc <- 34:42
gluc_new_splines <- X_splines_plot[,idx_gluc-n_disc]
gluc_splines_ps <- tcrossprod(gluc_new_splines, beta_ps[,idx_gluc])
gluc_splines_blm <- tcrossprod(gluc_new_splines, post_samp[,idx_gluc])
gluc_low <- gluc_upp <- rep(NA, length = NROW(gluc_new_splines))
gluc_low_blm <- gluc_upp_blm <- gluc_pm_blm <- rep(NA, length = NROW(gluc_new_splines))
for (i in 1:NROW(gluc_new_splines)) {
  ci <- robustBLME::hpd(gluc_splines_ps[i,])
  gluc_low[i] <- ci[1]
  gluc_upp[i] <- ci[2]
  gluc_pm_blm[i] <- mean(gluc_splines_blm[i,])
  ci_blm <- robustBLME::hpd(gluc_splines_blm[i,])
  gluc_low_blm[i] <- ci_blm[1]
  gluc_upp_blm[i] <- ci_blm[2]
}

idx_sbp <- 43:51
sbp_new_splines <- X_splines_plot[,idx_sbp-n_disc]
sbp_splines_ps <- tcrossprod(sbp_new_splines, beta_ps[,idx_sbp])
sbp_splines_blm <- tcrossprod(sbp_new_splines, post_samp[,idx_sbp])
sbp_low <- sbp_upp <- rep(NA, length = NROW(sbp_new_splines))
sbp_low_blm <- sbp_upp_blm <- sbp_pm_blm <- rep(NA, length = NROW(sbp_new_splines))
for (i in 1:NROW(sbp_new_splines)) {
  ci <- robustBLME::hpd(sbp_splines_ps[i,])
  sbp_low[i] <- ci[1]
  sbp_upp[i] <- ci[2]
  sbp_pm_blm[i] <- mean(sbp_splines_blm[i,])
  ci_blm <- robustBLME::hpd(sbp_splines_blm[i,])
  sbp_low_blm[i] <- ci_blm[1]
  sbp_upp_blm[i] <- ci_blm[2]
}

idx_pir <- 52:60
pir_new_splines <- X_splines_plot[,idx_pir-n_disc]
pir_splines_ps <- tcrossprod(pir_new_splines, beta_ps[,idx_pir])
pir_splines_blm <- tcrossprod(pir_new_splines, post_samp[,idx_pir])
pir_low <- pir_upp <- rep(NA, length = NROW(pir_new_splines))
pir_low_blm <- pir_upp_blm <- pir_pm_blm <- rep(NA, length = NROW(pir_new_splines))
for (i in 1:NROW(pir_new_splines)) {
  ci <- robustBLME::hpd(pir_splines_ps[i,])
  pir_low[i] <- ci[1]
  pir_upp[i] <- ci[2]
  pir_pm_blm[i] <- mean(pir_splines_blm[i,])
  ci_blm <- robustBLME::hpd(pir_splines_blm[i,])
  pir_low_blm[i] <- ci_blm[1]
  pir_upp_blm[i] <- ci_blm[2]
}

# Settings
lev_v <- c("HDL", "Glucose", "SBP")   # ordine dei pannelli (colonne)
lev_m <- c("RLS", "BLM")              # ordine dei pannelli (righe)

# Original X (no standardization)
x_hdl  <- hdl_new  * train_scale[3] + train_loc[3]
x_gluc <- gluc_new * train_scale[4] + train_loc[4]
x_sbp  <- sbp_new  * train_scale[5] + train_loc[5]

# Partial residuals for BLM
pres_blm_hdl  <- y - X_splines[, -idx_hdl]  %*% pm_blm[-idx_hdl]
pres_blm_gluc <- y - X_splines[, -idx_gluc] %*% pm_blm[-idx_gluc]
pres_blm_sbp  <- y - X_splines[, -idx_sbp]  %*% pm_blm[-idx_sbp]

# Estimate curves + bands in long format
mk_fit <- function(x, fit, low, upp, var, method) {
  data.frame(x = as.vector(x), fit = as.vector(fit),
             low = as.vector(low), upp = as.vector(upp),
             var = var, method = method)
}

df_fit <- rbind(
  mk_fit(x_hdl,  hdl_new_splines  %*% c(bb$beta)[idx_hdl],
         smooth(hdl_low, twiceit = TRUE),  smooth(hdl_upp, twiceit = TRUE),  "HDL",     "RLS"),
  mk_fit(x_gluc, gluc_new_splines %*% c(bb$beta)[idx_gluc],
         gluc_low, gluc_upp, "Glucose", "RLS"),
  mk_fit(x_sbp,  sbp_new_splines  %*% c(bb$beta)[idx_sbp],
         sbp_low,  sbp_upp,  "SBP",     "RLS"),
  mk_fit(x_hdl,  hdl_pm_blm,  hdl_low_blm,  hdl_upp_blm,  "HDL",     "BLM"),
  mk_fit(x_gluc, gluc_pm_blm, gluc_low_blm, gluc_upp_blm, "Glucose", "BLM"),
  mk_fit(x_sbp,  sbp_pm_blm,  sbp_low_blm,  sbp_upp_blm,  "SBP",     "BLM")
)

# Partial residuals in long format
mk_pts <- function(x, pres, var, method) {
  data.frame(x = as.vector(x), pres = as.vector(pres),
             var = var, method = method)
}

px_hdl  <- df_mgcv$hdl     * train_scale[3] + train_loc[3]
px_gluc <- df_mgcv$glucose * train_scale[4] + train_loc[4]
px_sbp  <- df_mgcv$sbp     * train_scale[5] + train_loc[5]

df_pts <- rbind(
  mk_pts(px_hdl,  pres_hdl,      "HDL",     "RLS"),
  mk_pts(px_gluc, pres_gluc,     "Glucose", "RLS"),
  mk_pts(px_sbp,  pres_sbp,      "SBP",     "RLS"),
  mk_pts(px_hdl,  pres_blm_hdl,  "HDL",     "BLM"),
  mk_pts(px_gluc, pres_blm_gluc, "Glucose", "BLM"),
  mk_pts(px_sbp,  pres_blm_sbp,  "SBP",     "BLM")
)

# Factors with selected order
df_fit$var    <- factor(df_fit$var,    levels = lev_v)
df_fit$method <- factor(df_fit$method, levels = lev_m)
df_pts$var    <- factor(df_pts$var,    levels = lev_v)
df_pts$method <- factor(df_pts$method, levels = lev_m)

# 3. Limit on y common across models, covariate-specific
rng <- do.call(rbind, lapply(lev_v, function(v) {
  r <- range(df_pts$pres[df_pts$var == v])
  data.frame(var = factor(v, levels = lev_v), ymin = r[1], ymax = r[2])
}))

# match() invece di merge(): non riordina le righe e non perde i fattori
i <- match(df_fit$var, rng$var)
df_fit$low <- pmax(df_fit$low, rng$ymin[i]) 
df_fit$upp <- pmin(df_fit$upp, rng$ymax[i])

# un solo layer geom_blank per fissare i limiti in ogni pannello
df_blank <- rbind(
  data.frame(var = rng$var, y = rng$ymin),
  data.frame(var = rng$var, y = rng$ymax)
)
df_blank <- do.call(rbind, lapply(lev_m, function(m) {
  transform(df_blank, method = factor(m, levels = lev_m))
}))

# Plot
lab_m <- c(RLS = "Predictive Bayes", BLM = "Traditional Bayes")

pl_appl <- ggplot(df_fit, aes(x = x)) +
  geom_ribbon(aes(ymin = low, ymax = upp, fill = method, alpha = method)) +
  geom_line(aes(y = fit, colour = method), linewidth = 0.8) +
  geom_point(aes(y = pres), data = df_pts, size = 0.2, alpha = 0.3) +
  geom_blank(aes(y = y), data = df_blank, inherit.aes = FALSE) +
  facet_grid(method ~ var, scales = "free",
             labeller = labeller(method = lab_m)) +
  scale_fill_manual(values   = c(RLS = "blue",  BLM = "orange")) +
  scale_colour_manual(values = c(RLS = "red2",  BLM = "darkgreen")) +
  scale_alpha_manual(values  = c(RLS = 0.2,     BLM = 0.3)) +
  labs(x = NULL, y = "Partial effect") +
  theme_light() +
  theme(legend.position = "none",
        strip.text = element_text(size = 11, colour = "black"),
        strip.background = element_rect(fill = "gray82"),
        panel.grid.major = element_line(size = 0.3, colour = "gray93"),
        panel.grid.minor = element_line(size = 0.15, colour = "gray93"),)

pl_appl
# ggsave(filename = "pl_appl.pdf", plot = pl_appl, width = 10, height = 7)




# FIGURE 1 SUPPLEMENTARY                                          
lev_v <- c("Age", "BMI", "PIR") 
lev_m <- c("RLS", "BLM") 
x_age  <- age_new*train_scale[1] + train_loc[1]
x_bmi <- bmi_new*train_scale[2] + train_loc[2]
x_pir  <- pir_new*train_scale[6] + train_loc[6]
pres_blm_age  <- y - X_splines[,-idx_age] %*% pm_blm[-idx_age]
pres_blm_bmi <- y - X_splines[,-idx_bmi] %*% pm_blm[-idx_bmi]
pres_blm_pir  <- y - X_splines[,-idx_pir] %*% pm_blm[-idx_pir]
mk_fit <- function(x, fit, low, upp, var, method) {
  data.frame(x = as.vector(x), fit = as.vector(fit),
             low = as.vector(low), upp = as.vector(upp),
             var = var, method = method)
}
df_fit <- rbind(
  mk_fit(x_age,  age_new_splines  %*% c(bb$beta)[idx_age],
         age_low,  age_upp,  "Age",     "RLS"),
  mk_fit(x_bmi, bmi_new_splines %*% c(bb$beta)[idx_bmi],
         smooth(bmi_low, twiceit = TRUE), smooth(bmi_upp, twiceit = TRUE), "BMI", "RLS"),
  mk_fit(x_pir,  pir_new_splines  %*% c(bb$beta)[idx_pir],
         pir_low,  pir_upp,  "PIR",     "RLS"),
  mk_fit(x_age,  age_pm_blm,  age_low_blm,  age_upp_blm,  "Age",     "BLM"),
  mk_fit(x_bmi, bmi_pm_blm, bmi_low_blm, bmi_upp_blm, "BMI", "BLM"),
  mk_fit(x_pir,  pir_pm_blm,  pir_low_blm,  pir_upp_blm,  "PIR",     "BLM")
)
mk_pts <- function(x, pres, var, method) {
  data.frame(x = as.vector(x), pres = as.vector(pres),
             var = var, method = method)
}
px_age <- df_mgcv$age*train_scale[1] + train_loc[1]
px_bmi <- df_mgcv$bmi*train_scale[2] + train_loc[2]
px_pir <- df_mgcv$pir*train_scale[6] + train_loc[6]
df_pts <- rbind(
  mk_pts(px_age, pres_age, "Age", "RLS"),
  mk_pts(px_bmi, pres_bmi, "BMI", "RLS"),
  mk_pts(px_pir, pres_pir, "PIR", "RLS"),
  mk_pts(px_age, pres_blm_age, "Age", "BLM"),
  mk_pts(px_bmi, pres_blm_bmi, "BMI", "BLM"),
  mk_pts(px_pir, pres_blm_pir, "PIR", "BLM")
)
df_fit$var <- factor(df_fit$var, levels = lev_v)
df_fit$method <- factor(df_fit$method, levels = lev_m)
df_pts$var <- factor(df_pts$var, levels = lev_v)
df_pts$method <- factor(df_pts$method, levels = lev_m)
rng <- do.call(rbind, lapply(lev_v, function(v) {
  r <- range(df_pts$pres[df_pts$var == v])
  data.frame(var = factor(v, levels = lev_v), ymin = r[1], ymax = r[2])
}))
i <- match(df_fit$var, rng$var)
df_fit$low <- pmax(df_fit$low, rng$ymin[i])
df_fit$upp <- pmin(df_fit$upp, rng$ymax[i])
df_blank <- rbind(
  data.frame(var = rng$var, y = rng$ymin),
  data.frame(var = rng$var, y = rng$ymax)
)
df_blank <- do.call(rbind, lapply(lev_m, function(m) {
  transform(df_blank, method = factor(m, levels = lev_m))
}))
lab_m <- c(RLS = "Predictive Bayes", BLM = "Traditional Bayes")
pl_appl_2 <- ggplot(df_fit, aes(x = x)) +
  geom_ribbon(aes(ymin = low, ymax = upp, fill = method, alpha = method)) +
  geom_line(aes(y = fit, colour = method), linewidth = 0.8) +
  geom_point(aes(y = pres), data = df_pts, size = 0.2, alpha = 0.3) +
  geom_blank(aes(y = y), data = df_blank, inherit.aes = FALSE) +
  facet_grid(method ~ var, scales = "free",
             labeller = labeller(method = lab_m)) +
  scale_fill_manual(values = c(RLS = "blue",  BLM = "orange")) +
  scale_colour_manual(values = c(RLS = "red2",  BLM = "darkgreen")) +
  scale_alpha_manual(values = c(RLS = 0.2,     BLM = 0.3)) +
  labs(x = NULL, y = "Partial effect") +
  theme_light() +
  theme(legend.position = "none",
        strip.text = element_text(size = 11, colour = "black"),
        strip.background = element_rect(fill = "gray82"),
        panel.grid.major = element_line(size = 0.3, colour = "gray93"),
        panel.grid.minor = element_line(size = 0.15, colour = "gray93"),)

pl_appl_2
# ggsave(filename = "pl_appl_2.pdf", plot = pl_appl_2, width = 10, height = 7)







# FIGURE 2 SUPPLEMENTARY
grids <- list(age = age_new, bmi = bmi_new, hdl = hdl_new,
              glucose = gluc_new, sbp = sbp_new, pir = pir_new)

base <- df_mgcv[1, all.vars(formula(fit_t))[-1], drop = FALSE]

pres_t <- lapply(names(grids), function(v) {
  nd <- base[rep(1, length(grids[[v]])), , drop = FALSE]
  nd[[v]] <- grids[[v]]
  p    <- predict(fit_t, nd, type = "terms", se.fit = TRUE)
  term <- paste0("s(", v, ")")
  data.frame(var = v, x = grids[[v]],
             est = p$fit[, term],
             se  = p$se.fit[, term],
             lwr = p$fit[, term] - 1.96 * p$se.fit[, term],
             upr = p$fit[, term] + 1.96 * p$se.fit[, term])
})
age_low_t <- pres_t[[1]]$lwr
age_upp_t <- pres_t[[1]]$upr
age_est_t <- pres_t[[1]]$est
bmi_low_t <- pres_t[[2]]$lwr
bmi_upp_t <- pres_t[[2]]$upr
bmi_est_t <- pres_t[[2]]$est
hdl_low_t <- pres_t[[3]]$lwr
hdl_upp_t <- pres_t[[3]]$upr
hdl_est_t <- pres_t[[3]]$est
gluc_low_t <- pres_t[[4]]$lwr
gluc_upp_t <- pres_t[[4]]$upr
gluc_est_t <- pres_t[[4]]$est
sbp_low_t <- pres_t[[5]]$lwr
sbp_upp_t <- pres_t[[5]]$upr
sbp_est_t <- pres_t[[5]]$est
pir_low_t <- pres_t[[6]]$lwr
pir_upp_t <- pres_t[[6]]$upr
pir_est_t <- pres_t[[6]]$est
                                                 
lev_v <- c("HDL", "Glucose", "SBP") 
lev_m <- c("RLS", "AMt") 
x_hdl  <- hdl_new*train_scale[3] + train_loc[3]
x_gluc <- gluc_new*train_scale[4] + train_loc[4]
x_sbp  <- sbp_new*train_scale[5] + train_loc[5]
pres_t_hdl  <- y - X_splines[,-idx_hdl] %*% fit_t$coefficients[-idx_hdl]
pres_t_gluc <- y - X_splines[,-idx_gluc] %*% fit_t$coefficients[-idx_gluc]
pres_t_sbp  <- y - X_splines[,-idx_sbp] %*% fit_t$coefficients[-idx_sbp]
mk_fit <- function(x, fit, low, upp, var, method) {
  data.frame(x = as.vector(x), fit = as.vector(fit),
             low = as.vector(low), upp = as.vector(upp),
             var = var, method = method)
}
df_fit <- rbind(
  mk_fit(x_hdl,  hdl_new_splines  %*% c(bb$beta)[idx_hdl],
         smooth(hdl_low, twiceit = TRUE), smooth(hdl_upp, twiceit = TRUE),  "HDL",     "RLS"),
  mk_fit(x_gluc, gluc_new_splines %*% c(bb$beta)[idx_gluc],
         gluc_low, gluc_upp, "Glucose", "RLS"),
  mk_fit(x_sbp,  sbp_new_splines  %*% c(bb$beta)[idx_sbp],
         sbp_low,  sbp_upp,  "SBP",     "RLS"),
  mk_fit(x_hdl,  hdl_est_t,  hdl_low_t,  hdl_upp_t,  "HDL",     "AMt"),
  mk_fit(x_gluc, gluc_est_t, gluc_low_t, gluc_upp_t, "Glucose", "AMt"),
  mk_fit(x_sbp,  sbp_est_t,  sbp_low_t,  sbp_upp_t,  "SBP",     "AMt")
)
mk_pts <- function(x, pres, var, method) {
  data.frame(x = as.vector(x), pres = as.vector(pres),
             var = var, method = method)
}
px_hdl <- df_mgcv$hdl*train_scale[3] + train_loc[3]
px_gluc <- df_mgcv$gluc*train_scale[4] + train_loc[4]
px_sbp <- df_mgcv$sbp*train_scale[5] + train_loc[5]
df_pts <- rbind(
  mk_pts(px_hdl, pres_hdl, "HDL", "RLS"),
  mk_pts(px_gluc, pres_gluc, "Glucose", "RLS"),
  mk_pts(px_sbp, pres_sbp, "SBP", "RLS"),
  mk_pts(px_hdl, pres_t_hdl, "HDL", "AMt"),
  mk_pts(px_gluc, pres_t_gluc, "Glucose", "AMt"),
  mk_pts(px_sbp, pres_t_sbp, "SBP", "AMt")
)
df_fit$var <- factor(df_fit$var, levels = lev_v)
df_fit$method <- factor(df_fit$method, levels = lev_m)
df_pts$var <- factor(df_pts$var, levels = lev_v)
df_pts$method <- factor(df_pts$method, levels = lev_m)
rng <- do.call(rbind, lapply(lev_v, function(v) {
  r <- range(df_pts$pres[df_pts$var == v])
  data.frame(var = factor(v, levels = lev_v), ymin = r[1], ymax = r[2])
}))
i <- match(df_fit$var, rng$var)
df_fit$low <- pmax(df_fit$low, rng$ymin[i])
df_fit$upp <- pmin(df_fit$upp, rng$ymax[i])
df_blank <- rbind(
  data.frame(var = rng$var, y = rng$ymin),
  data.frame(var = rng$var, y = rng$ymax)
)
df_blank <- do.call(rbind, lapply(lev_m, function(m) {
  transform(df_blank, method = factor(m, levels = lev_m))
}))
lab_m <- c(RLS = "Predictive Bayes", AMt = "Frequentist Additive Model")
pl_appl_3 <- ggplot(df_fit, aes(x = x)) +
  geom_ribbon(aes(ymin = low, ymax = upp, fill = method, alpha = method)) +
  geom_line(aes(y = fit, colour = method), linewidth = 0.8) +
  geom_point(aes(y = pres), data = df_pts, size = 0.2, alpha = 0.3) +
  geom_blank(aes(y = y), data = df_blank, inherit.aes = FALSE) +
  facet_grid(method ~ var, scales = "free",
             labeller = labeller(method = lab_m)) +
  scale_fill_manual(values = c(RLS = "blue", AMt = "orange")) +
  scale_colour_manual(values = c(RLS = "red2", AMt = "darkgreen")) +
  scale_alpha_manual(values = c(RLS = 0.2, AMt = 0.3)) +
  labs(x = NULL, y = "Partial effect") +
  theme_light() +
  theme(legend.position = "none",
        strip.text = element_text(size = 11, colour = "black"),
        strip.background = element_rect(fill = "gray82"),
        panel.grid.major = element_line(size = 0.3, colour = "gray93"),
        panel.grid.minor = element_line(size = 0.15, colour = "gray93"),)
pl_appl_3
# ggsave(filename = "pl_appl_3.pdf", plot = pl_appl_4, width = 10, height = 7)


                                                 
# FIGURE 3 SUPPLEMENTARY
lev_v <- c("Age", "BMI", "PIR") 
lev_m <- c("RLS", "AMt") 
x_age  <- age_new*train_scale[1] + train_loc[1]
x_bmi <- bmi_new*train_scale[2] + train_loc[2]
x_pir  <- pir_new*train_scale[6] + train_loc[6]
pres_t_age  <- y - X_splines[,-idx_age] %*% fit_t$coefficients[-idx_age]
pres_t_bmi <- y - X_splines[,-idx_bmi] %*% fit_t$coefficients[-idx_bmi]
pres_t_pir  <- y - X_splines[,-idx_pir] %*% fit_t$coefficients[-idx_pir]
mk_fit <- function(x, fit, low, upp, var, method) {
  data.frame(x = as.vector(x), fit = as.vector(fit),
             low = as.vector(low), upp = as.vector(upp),
             var = var, method = method)
}
df_fit <- rbind(
  mk_fit(x_age,  age_new_splines  %*% c(bb$beta)[idx_age],
         age_low,  age_upp,  "Age",     "RLS"),
  mk_fit(x_bmi, bmi_new_splines %*% c(bb$beta)[idx_bmi],
         smooth(bmi_low, twiceit = TRUE), smooth(bmi_upp, twiceit = TRUE), "BMI", "RLS"),
  mk_fit(x_pir,  pir_new_splines  %*% c(bb$beta)[idx_pir],
         pir_low,  pir_upp,  "PIR",     "RLS"),
  mk_fit(x_age,  age_est_t,  age_low_t,  age_upp_t,  "Age",     "AMt"),
  mk_fit(x_bmi, bmi_est_t, bmi_low_t, bmi_upp_t, "BMI", "AMt"),
  mk_fit(x_pir,  pir_est_t,  pir_low_t,  pir_upp_t,  "PIR",     "AMt")
)
mk_pts <- function(x, pres, var, method) {
  data.frame(x = as.vector(x), pres = as.vector(pres),
             var = var, method = method)
}
px_age <- df_mgcv$age*train_scale[1] + train_loc[1]
px_bmi <- df_mgcv$bmi*train_scale[2] + train_loc[2]
px_pir <- df_mgcv$pir*train_scale[6] + train_loc[6]
df_pts <- rbind(
  mk_pts(px_age, pres_age, "Age", "RLS"),
  mk_pts(px_bmi, pres_bmi, "BMI", "RLS"),
  mk_pts(px_pir, pres_pir, "PIR", "RLS"),
  mk_pts(px_age, pres_t_age, "Age", "AMt"),
  mk_pts(px_bmi, pres_t_bmi, "BMI", "AMt"),
  mk_pts(px_pir, pres_t_pir, "PIR", "AMt")
)
df_fit$var <- factor(df_fit$var, levels = lev_v)
df_fit$method <- factor(df_fit$method, levels = lev_m)
df_pts$var <- factor(df_pts$var, levels = lev_v)
df_pts$method <- factor(df_pts$method, levels = lev_m)
rng <- do.call(rbind, lapply(lev_v, function(v) {
  r <- range(df_pts$pres[df_pts$var == v])
  data.frame(var = factor(v, levels = lev_v), ymin = r[1], ymax = r[2])
}))
i <- match(df_fit$var, rng$var)
df_fit$low <- pmax(df_fit$low, rng$ymin[i])
df_fit$upp <- pmin(df_fit$upp, rng$ymax[i])

df_blank <- rbind(
  data.frame(var = rng$var, y = rng$ymin),
  data.frame(var = rng$var, y = rng$ymax)
)
df_blank <- do.call(rbind, lapply(lev_m, function(m) {
  transform(df_blank, method = factor(m, levels = lev_m))
}))
lab_m <- c(RLS = "Predictive Bayes", AMt = "Frequentist Additive Model")
pl_appl_4 <- ggplot(df_fit, aes(x = x)) +
  geom_ribbon(aes(ymin = low, ymax = upp, fill = method, alpha = method)) +
  geom_line(aes(y = fit, colour = method), linewidth = 0.8) +
  geom_point(aes(y = pres), data = df_pts, size = 0.2, alpha = 0.3) +
  geom_blank(aes(y = y), data = df_blank, inherit.aes = FALSE) +
  facet_grid(method ~ var, scales = "free",
             labeller = labeller(method = lab_m)) +
  scale_fill_manual(values = c(RLS = "blue", AMt = "orange")) +
  scale_colour_manual(values = c(RLS = "red2", AMt = "darkgreen")) +
  scale_alpha_manual(values = c(RLS = 0.2, AMt = 0.3)) +
  labs(x = NULL, y = "Partial effect") +
  theme_light() +
  theme(legend.position = "none",
        strip.text = element_text(size = 11, colour = "black"),
        strip.background = element_rect(fill = "gray82"),
        panel.grid.major = element_line(size = 0.3, colour = "gray93"),
        panel.grid.minor = element_line(size = 0.15, colour = "gray93"),)
pl_appl_4
# ggsave(filename = "pl_appl_4.pdf", plot = pl_appl_5, width = 10, height = 7)

                                                 


# FIGURE CONVERGENCE
df_conv <- data.frame(conv = c(c_ps),
                      fs = rep(1:NROW(c_ps), NCOL(c_ps)),
                      group = rep(1:NCOL(c_ps), each = NROW(c_ps)),
                      title = "Convergence")
pl_conv <- ggplot(df_conv) +
  facet_wrap(~ title) +
  geom_line(aes(x = fs, y = conv, group = group), alpha = 0.7) +
  labs(x = "Number of forward samples", y = TeX("$L_1$ distance")) +
  ylim(c(0, max(c_ps))) +
  theme_light() +
  theme(legend.text.align = 0,
        strip.text = element_text(size = 12, colour = "black"),
        strip.background = element_rect(fill = "gray82"),
        panel.grid.major = element_line(size = 0.3, colour = "gray93"),
        panel.grid.minor = element_line(size = 0.15, colour = "gray93"),
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 8),
        axis.title = element_text(size = 12),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5),
        legend.position = "top")
pl_conv
# ggsave(filename = "pl_appl_conv.png", plot = pl_conv, width = 10, height = 7)





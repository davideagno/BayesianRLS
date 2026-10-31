
rm(list = ls())

Rcpp::sourceCpp("RLS_functions.cpp")
source("RLS_functions.R")

library(coda)
library(extraDistr)
library(brms)
library(mgcv)
library(ggplot2)


load("appl_data.RData")
contrasts(df$ethn) <- contr.treatment(levels(df$ethn), base = 3)

set.seed(2026)
size_train <- round(0.7*NROW(df))
idx_train <- sort(sample(1:NROW(df), size = size_train, replace = FALSE))
idx_test <- sort(setdiff(1:NROW(df), idx_train))

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


fit_a <- mgcv::gam(y ~ s(age) + s(bmi) + s(hdl) + s(glucose) + s(sbp) + s(pir) +
                     ethn1 + ethn2 + ethn4 + ethn6 + ethn7, data = df_mgcv, 
                   method = "REML", family = gaussian())
res <- residuals(fit_a)
nus <- c(9, 10000)
n   <- length(res)
# punti del Q-Q plot per ciascun nu
df_qq <- do.call(rbind, lapply(nus, function(nu) {
  data.frame(nu = paste("df =", nu),
             teorici = qt(ppoints(n), df = nu),
             res = sort(res))
}))
# retta come in qqline(): passa per il 1° e 3° quartile
df_line <- do.call(rbind, lapply(nus, function(nu) {
  y <- unname(quantile(res, c(0.25, 0.75)))
  x <- qt(c(0.25, 0.75), df = nu)
  slope <- diff(y) / diff(x)
  data.frame(nu = paste("df =", nu), slope = slope, intercept = y[1] - slope * x[1])
}))
livelli <- paste("df =", nus)
titoli <- c('"Student-t" ~ (nu == 9)', 'Gaussian')
# ...stesse due righe factor() di sopra...
df_qq$nu   <- factor(df_qq$nu,   levels = livelli, labels = titoli)
df_line$nu <- factor(df_line$nu, levels = livelli, labels = titoli)
pl_qq <- ggplot(df_qq, aes(x = teorici, y = res)) +
  geom_point() +
  geom_abline(data = df_line, aes(slope = slope, intercept = intercept),
              colour = "red") +
  facet_wrap(~ nu, scales = "free_x", labeller = label_parsed) +
  labs(x = "Theoretical quantiles", y = "Residuals") +
  theme_light() +
  theme(legend.position = "top",
        strip.text = element_text(size = 11, colour = "black"),
        strip.background = element_rect(fill = "gray82"),
        panel.grid.major = element_line(size = 0.3, colour = "gray93"),
        panel.grid.minor = element_line(size = 0.15, colour = "gray93"),)
pl_qq
# ggsave(filename = "pl_appl_qq.pdf", plot = pl_qq, width = 10, height = 7)

v_ethn <- rep(NA, 5)
v_ethn[1] <- var(y[which(df_train$ethn == 1)]); length(which(df_train$ethn == 1))
v_ethn[2] <- var(y[which(df_train$ethn == 2)]); length(which(df_train$ethn == 2))
v_ethn[3] <- var(y[which(df_train$ethn == 3)]); length(which(df_train$ethn == 3))
v_ethn[4] <- var(y[which(df_train$ethn == 4)]); length(which(df_train$ethn == 4))
v_ethn[5] <- var(y[which(df_train$ethn == 6 | df_train$ethn == 7)]); length(which(df_train$ethn == 6 | df_train$ethn == 7))

v_hdl <- rep(NA, 4)
v_hdl[1] <- var(y[which(df_train$hdl > 0 & df_train$hdl <= quantile(df_train$hdl, 0.25))])
v_hdl[2] <- var(y[which(df_train$hdl > quantile(df_train$hdl, 0.25) & df_train$hdl <= quantile(df_train$hdl, 0.50))])
v_hdl[3] <- var(y[which(df_train$hdl > quantile(df_train$hdl, 0.50) & df_train$hdl <= quantile(df_train$hdl, 0.75))])
v_hdl[4] <- var(y[which(df_train$hdl > quantile(df_train$hdl, 0.75) & df_train$hdl <= max(df_train$hdl))])

v_glucose <- rep(NA, 4)
v_glucose[1] <- var(y[which(df_train$glucose > 0 & df_train$glucose <= quantile(df_train$glucose, 0.25))])
v_glucose[2] <- var(y[which(df_train$glucose > quantile(df_train$glucose, 0.25) & df_train$glucose <= quantile(df_train$glucose, 0.50))])
v_glucose[3] <- var(y[which(df_train$glucose > quantile(df_train$glucose, 0.50) & df_train$glucose <= quantile(df_train$glucose, 0.75))])
v_glucose[4] <- var(y[which(df_train$glucose > quantile(df_train$glucose, 0.75) & df_train$glucose <= max(df_train$glucose))])

v_sbp <- rep(NA, 4)
v_sbp[1] <- var(y[which(df_train$sbp > 0 & df_train$sbp <= quantile(df_train$sbp, 0.25))])
v_sbp[2] <- var(y[which(df_train$sbp > quantile(df_train$sbp, 0.25) & df_train$sbp <= quantile(df_train$sbp, 0.50))])
v_sbp[3] <- var(y[which(df_train$sbp > quantile(df_train$sbp, 0.50) & df_train$sbp <= quantile(df_train$sbp, 0.75))])
v_sbp[4] <- var(y[which(df_train$sbp > quantile(df_train$sbp, 0.75) & df_train$sbp <= max(df_train$sbp))])

v_age <- rep(NA, 4)
v_age[1] <- var(y[which(df_train$age > 0 & df_train$age <= quantile(df_train$age, 0.25))])
v_age[2] <- var(y[which(df_train$age > quantile(df_train$age, 0.25) & df_train$age <= quantile(df_train$age, 0.50))])
v_age[3] <- var(y[which(df_train$age > quantile(df_train$age, 0.50) & df_train$age <= quantile(df_train$age, 0.75))])
v_age[4] <- var(y[which(df_train$age > quantile(df_train$age, 0.75) & df_train$age <= max(df_train$age))])

v_bmi <- rep(NA, 4)
v_bmi[1] <- var(y[which(df_train$bmi > 0 & df_train$bmi <= quantile(df_train$bmi, 0.25))])
v_bmi[2] <- var(y[which(df_train$bmi > quantile(df_train$bmi, 0.25) & df_train$bmi <= quantile(df_train$bmi, 0.50))])
v_bmi[3] <- var(y[which(df_train$bmi > quantile(df_train$bmi, 0.50) & df_train$bmi <= quantile(df_train$bmi, 0.75))])
v_bmi[4] <- var(y[which(df_train$bmi > quantile(df_train$bmi, 0.75) & df_train$bmi <= max(df_train$bmi))])

v_pir <- rep(NA, 4)
v_pir[1] <- var(y[which(df_train$pir > 0 & df_train$pir <= quantile(df_train$pir, 0.25))])
v_pir[2] <- var(y[which(df_train$pir > quantile(df_train$pir, 0.25) & df_train$pir <= quantile(df_train$pir, 0.50))])
v_pir[3] <- var(y[which(df_train$pir > quantile(df_train$pir, 0.50) & df_train$pir <= quantile(df_train$pir, 0.75))])
v_pir[4] <- var(y[which(df_train$pir == 5.00)])

variabili <- c("PIR", "BMI", "Age", "Glucose", "SBP", "Ethnicity", "HDL")
df_var <- data.frame(
  gruppo    = rep(1:4, times = length(variabili)),
  varianza  = c(v_pir, v_bmi, v_age, v_glucose, v_sbp, v_ethn[1:4], v_hdl),
  variabile = factor(rep(variabili, each = 4), levels = variabili)
)
colori <- setNames(palette()[c(1, 8, 4, 5, 6, 3, 2)], variabili)
pl_var <- ggplot(df_var, aes(x = factor(gruppo), y = varianza,
                             colour = variabile, group = variabile)) +
  geom_line() +
  geom_point(size = 2) +
  scale_colour_manual(values = colori) +
  scale_x_discrete(labels = paste0("Q", 1:4)) +
  coord_cartesian(ylim = c(0, 0.8)) +
  labs(x = "Covariates quartile", y = "Response variance", colour = "Covariate") +
  theme_light() +
  theme(legend.position = "top",
        strip.text = element_text(size = 11, colour = "black"),
        strip.background = element_rect(fill = "gray82"),
        panel.grid.major = element_line(size = 0.3, colour = "gray93"),
        panel.grid.minor = element_line(size = 0.15, colour = "gray93"),)
pl_var
# ggsave(filename = "pl_appl_var.pdf", plot = pl_var, width = 10, height = 7)



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


n_disc <- NCOL(mm[,-1])
var_order <- c(6 + 1:n_disc, 1:6)
n_par <- 1 + n_disc 
colnames(covariates)[var_order]

smooths <- lapply(var_order[-c(1:n_disc)], function(j) {
  spec <- eval(call("s", as.name(colnames(covariates)[j])))
  mgcv::smoothCon(spec, data = covariates, knots = NULL, absorb.cons = TRUE)[[1]]
})

# Omega contiene SOLO le penalita' delle smooth
Z <- Omega <- list(); X_splines <- rep(1, length(y))
for (j in 1:n_disc) {
  Z[[j]] <- covariates[,var_order[j]]
  X_splines <- cbind(X_splines, Z[[j]])
}
for (j in (n_disc+1):length(var_order)) {
  Z[[j]] <- smooths[[j-n_disc]]$X
  Omega[[j-n_disc]] <- smooths[[j-n_disc]]$S[[1]]
  X_splines <- cbind(X_splines, Z[[j]])
}
colSums(X_splines[, 2:n_par])           # ogni dummy deve avere almeno 2 uni
dim(X_splines)

X_splines_test <- rep(1, length(y_test))
for (j in 1:n_disc) {
  X_splines_test <- cbind(X_splines_test, covariates_test[,var_order[j]])
}
for (j in (n_disc+1):length(var_order)) {
  X_splines_test <- cbind(X_splines_test, mgcv::PredictMat(smooths[[j-n_disc]], df_test))
}

b_0 <- rep(0, NCOL(X_splines))

# Optimization
opt <- nlminb(start = c(log(5), log(1), 0, rep(0.5, length(Omega))), objective = npll_loo_alt_appl_skt,
              y = y, X = X_splines, beta_0 = b_0, Omega = Omega, n_par = n_par)
opt
(df_opt <- exp(opt$par[1]))
(sigma_opt <- exp(opt$par[2]))
(alpha_opt <- opt$par[3])
lambda_opt <- exp(opt$par[-c(1:3)])

t_opt <- system.time(opt <- nlminb(start = c(log(5), log(1), rep(0.5, length(Omega))), objective = npll_loo_alt_appl,
                                   y = y, X = X_splines, beta_0 = b_0, Omega = Omega, n_par = n_par))
t_opt
opt
(df_opt <- exp(opt$par[1]))
(sigma_opt <- exp(opt$par[2]))
lambda_opt <- exp(opt$par[-c(1:2)])

t_est <- system.time(bb <- rls_lambda_appl(y, X_splines, b_0, lambda_opt, Omega, n_par))

pred_rls_test <- c(X_splines_test %*% bb$beta)
mean((df_test$y - pred_rls_test)^2)

# Convergence
set.seed(123)
c_ps <- ps_rs_conv_t_1_alt(X_splines, bb$beta, bb$P, 4e5, 1, df_opt, sigma_opt)
plot(c_ps, type = "l")

# Predictive resampling
t_ps <- system.time(beta_ps <- ps_rs_parallel_t_fixed_alt(X_splines, bb$beta, bb$P, 
                                                          N = 3e5, B = 10000, 
                                                          df = df_opt, sigma = sigma_opt, 
                                                          seed = 1, nthreads = 7))
beta_ps <- as.mcmc(beta_ps)
# save(beta_ps, file = "ps_appl_10k_6.RData")
# load("ps_appl_10k_6.RData")



# Modello bayesiano tradizionale
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
  Om <- (Om + t(Om)) / 2
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
  target += gamma_lpdf(lambda | 0.001, 0.001);
  // beta ~ N(0, sigma^2 * S^{-1}),  S = blockdiag(lambda_j * Omega_j)
  // blocks may differ in size: block j occupies b[strt[j] : strt[j]+dblk[j]-1]
  for (j in 1:P) {
    vector[dblk[j]] beta_j = segment(b, strt[j], dblk[j]);
    matrix[dblk[j], dblk[j]] Om_j = block(Omega[j], 1, 1, dblk[j], dblk[j]);
    target += -0.5 * (lambda[j] / square(sigma)) * quad_form_sym(Om_j, beta_j)
              + 0.5 * rk[j] * log(lambda[j]) - rk[j] * log(sigma);
  }",
          block = "model")
t_blm <- system.time(blm <- brms::brm(formula = y ~ -1+., data = dat_splines, family = student(),
                                      stanvars = stanvars, prior = prior("", class = "b"),
                                      chains = 5, cores = 5, iter = 20000, warmup = 18000, seed = 1))
post_samp <- coda::as.mcmc(as.matrix(blm))
# save(blm, post_samp, file = "blm_appl_10k_6.RData")
# load("blm_appl_10k_6.RData")
coda::effectiveSize(post_samp)
mean(coda::effectiveSize(post_samp[,1:NCOL(X_splines)]))


jdx <- 30
coda::effectiveSize(post_samp[,jdx])
s1 <- as.numeric(post_samp[, jdx])
s2 <- as.numeric(beta_ps[, jdx])

# Effective sample size, riportata nelle etichette dei pannelli
ess <- c(coda::effectiveSize(s1), coda::effectiveSize(s2))
lab <- c("Traditional Bayes", "Predictive Bayes")

# Dati in formato lungo
df_tr <- rbind(
  data.frame(iter = seq_along(s1), value = s1, chain = lab[1]),
  data.frame(iter = seq_along(s2), value = s2, chain = lab[2])
)
df_tr$chain <- factor(df_tr$chain, levels = lab)

# Linee verticali tratteggiate solo nel primo pannello
df_vl <- data.frame(chain = factor(lab[1], levels = lab),
                    x     = c(2, 4, 6, 8) * 1e3)

pl_trace <- ggplot(df_tr, aes(x = iter, y = value)) +
  geom_line(linewidth = 0.3) +
  geom_vline(aes(xintercept = x), data = df_vl,
             linetype = "dashed", colour = "gray40", linewidth = 0.2) +
  facet_wrap(~ chain, nrow = 1, scales = "free_x") +
  # stessi limiti y del primo pannello anche nel secondo
  # (equivale a ylim = range(post_samp[, jdx]) del plot base)
  coord_cartesian(ylim = range(s1)) +
  labs(x = "Iteration", y = "Regression Coefficient") +
  theme_light() +
  theme(strip.text       = element_text(size = 11, colour = "black"),
        strip.background = element_rect(fill = "gray82"),
        panel.grid.major = element_line(linewidth = 0.3,  colour = "gray93"),
        panel.grid.minor = element_line(linewidth = 0.15, colour = "gray93"))

pl_trace
# ggsave(filename = "trace_plot.pdf", plot = pl_trace, width = 10, height = 3.5)


# save(t_opt, t_est, t_ps, t_blm, file = "time_appl_6.RData")
# load("time_appl_6.RData")
t_opt[3] + t_est[3] + t_ps[3]
t_blm[3]



fit_t <- mgcv::gam(y ~ s(age) + s(bmi) + s(hdl) + s(glucose) + s(sbp) + s(pir) +
                     ethn1 + ethn2 + ethn4 + ethn6 + ethn7, data = df_mgcv, 
                   method = "REML", family = scat())
fit_n <- mgcv::gam(y ~ s(age) + s(bmi) + s(hdl) + s(glucose) + s(sbp) + s(pir) +
                     ethn1 + ethn2 + ethn4 + ethn6 + ethn7, data = df_mgcv, 
                   method = "REML", family = gaussian())



pred_rls_test <- c(X_splines_test %*% bb$beta)
mse_rls <- mean((df_test$y - pred_rls_test)^2)
pm_blm <- colMeans(post_samp[,1:NCOL(X_splines)])
pred_blm_test <- c(X_splines_test %*% pm_blm)
mse_blm <- mean((df_test$y - pred_blm_test)^2)
pred_am_t_test <- c(X_splines_test %*% fit_t$coefficients)
mse_am_t <- mean((df_test$y - pred_am_t_test)^2)
pred_am_n_test <- c(X_splines_test %*% fit_n$coefficients)
mse_am_n <- mean((df_test$y - pred_am_n_test)^2)
mse_rls; mse_am_t; mse_am_n; mse_blm


fit_pl <- mgcv::gam(y ~ s(age) + s(bmi) + s(hdl) + s(glucose) + s(sbp) + s(pir) +
                      ethn1 + ethn2 + ethn4 + ethn6 + ethn7, data = df_mgcv, 
                    sp = pmin(lambda_opt, 1e6))

bb_chk <- rls_lambda_appl(y, X_splines, b_0, pmin(lambda_opt, 1e6), Omega, n_par = n_par)
max(abs(coef(fit_pl) - bb_chk$beta)) 

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

n_plot <- 500
age_new <- seq(min(df_mgcv$age), max(df_mgcv$age), length = n_plot)
bmi_new <- seq(min(df_mgcv$bmi), max(df_mgcv$bmi), length = n_plot)
hdl_new <- seq(min(df_mgcv$hdl), max(df_mgcv$hdl), length = n_plot)
gluc_new <- seq(min(df_mgcv$gluc), max(df_mgcv$gluc), length = n_plot)
sbp_new <- seq(min(df_mgcv$sbp), max(df_mgcv$sbp), length = n_plot)
pir_new <- seq(min(df_mgcv$pir), max(df_mgcv$pir), length = n_plot)
X_plot <- as.data.frame(cbind(age = age_new, bmi = bmi_new, hdl = hdl_new, 
                              glucose = gluc_new, sbp = sbp_new, pir = pir_new))
X_splines_plot <- cbind(1, do.call(cbind, lapply(smooths, function(sm) mgcv::PredictMat(sm, X_plot))))


idx_age <- n_par + 1:9
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


idx_bmi <- n_par + 10:18
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


idx_hdl <- n_par + 19:27
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


idx_gluc <- n_par + 28:36
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


idx_sbp <- n_par + 37:45
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


idx_pir <- n_par + 46:54
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


# Effetti parziali delle 6 covariate - solo il nostro modello (Predictive Bayes, RLS)
# Griglia 2 x 3

lev_v <- c("Age", "BMI", "HDL", "Glucose", "SBP", "PIR")   # ordine dei pannelli

# 1. Ascisse sulla scala originale (de-standardizzate)
x_age  <- age_new  * train_scale[1] + train_loc[1]
x_bmi  <- bmi_new  * train_scale[2] + train_loc[2]
x_hdl  <- hdl_new  * train_scale[3] + train_loc[3]
x_gluc <- gluc_new * train_scale[4] + train_loc[4]
x_sbp  <- sbp_new  * train_scale[5] + train_loc[5]
x_pir  <- pir_new  * train_scale[6] + train_loc[6]

px_age  <- df_mgcv$age     * train_scale[1] + train_loc[1]
px_bmi  <- df_mgcv$bmi     * train_scale[2] + train_loc[2]
px_hdl  <- df_mgcv$hdl     * train_scale[3] + train_loc[3]
px_gluc <- df_mgcv$glucose * train_scale[4] + train_loc[4]
px_sbp  <- df_mgcv$sbp     * train_scale[5] + train_loc[5]
px_pir  <- df_mgcv$pir     * train_scale[6] + train_loc[6]

beta_rls <- c(bb$beta)

# 2. Curve stimate + bande, in formato lungo
mk_fit <- function(x, fit, low, upp, var) {
  data.frame(x = as.vector(x), fit = as.vector(fit),
             low = as.vector(low), upp = as.vector(upp),
             var = var)
}

df_fit_rls <- rbind(
  mk_fit(x_age,  age_new_splines  %*% beta_rls[idx_age],
         age_low, age_upp, "Age"),
  mk_fit(x_bmi,  bmi_new_splines  %*% beta_rls[idx_bmi],
         smooth(bmi_low, twiceit = TRUE), smooth(bmi_upp, twiceit = TRUE), "BMI"),
  mk_fit(x_hdl,  hdl_new_splines  %*% beta_rls[idx_hdl],
         smooth(hdl_low, twiceit = TRUE), smooth(hdl_upp, twiceit = TRUE), "HDL"),
  mk_fit(x_gluc, gluc_new_splines %*% beta_rls[idx_gluc],
         gluc_low, gluc_upp, "Glucose"),
  mk_fit(x_sbp,  sbp_new_splines  %*% beta_rls[idx_sbp],
         sbp_low, sbp_upp, "SBP"),
  mk_fit(x_pir,  pir_new_splines  %*% beta_rls[idx_pir],
         pir_low, pir_upp, "PIR")
)

# 3. Residui parziali, in formato lungo
mk_pts <- function(x, pres, var) {
  data.frame(x = as.vector(x), pres = as.vector(pres), var = var)
}

df_pts_rls <- rbind(
  mk_pts(px_age,  pres_age,  "Age"),
  mk_pts(px_bmi,  pres_bmi,  "BMI"),
  mk_pts(px_hdl,  pres_hdl,  "HDL"),
  mk_pts(px_gluc, pres_gluc, "Glucose"),
  mk_pts(px_sbp,  pres_sbp,  "SBP"),
  mk_pts(px_pir,  pres_pir,  "PIR")
)

df_fit_rls$var <- factor(df_fit_rls$var, levels = lev_v)
df_pts_rls$var <- factor(df_pts_rls$var, levels = lev_v)

# 4. Ritaglio delle bande al range dei residui parziali di ogni covariata
#    (con un solo metodo per pannello non serve piu' geom_blank:
#     scales = "free" fissa da solo i limiti sui punti)
rng <- do.call(rbind, lapply(lev_v, function(v) {
  r <- range(df_pts_rls$pres[df_pts_rls$var == v])
  data.frame(var = factor(v, levels = lev_v), ymin = r[1], ymax = r[2])
}))

i <- match(df_fit_rls$var, rng$var)
df_fit_rls$low <- pmax(df_fit_rls$low, rng$ymin[i])
df_fit_rls$upp <- pmin(df_fit_rls$upp, rng$ymax[i])

# 5. Grafico
pl_rls <- ggplot(df_fit_rls, aes(x = x)) +
  geom_ribbon(aes(ymin = low, ymax = upp), fill = "blue", alpha = 0.2) +
  geom_line(aes(y = fit), colour = "red2", linewidth = 0.8) +
  geom_point(aes(y = pres), data = df_pts_rls, size = 0.2, alpha = 0.3) +
  facet_wrap(~ var, nrow = 2, ncol = 3, scales = "free") +
  labs(x = NULL, y = "Partial effect") +
  theme_light() +
  theme(strip.text       = element_text(size = 11, colour = "black"),
        strip.background = element_rect(fill = "gray82"),
        panel.grid.major = element_line(linewidth = 0.3,  colour = "gray93"),
        panel.grid.minor = element_line(linewidth = 0.15, colour = "gray93"))

pl_rls
# ggsave(filename = "pl_rls.pdf", plot = pl_rls, width = 10, height = 6.5)








# Effetti parziali delle 6 covariate - solo il modello BLM (Traditional Bayes)
# Griglia 2 x 3

# 2. Residui parziali BLM
pres_blm_age  <- y - X_splines[, -idx_age]  %*% pm_blm[-idx_age]
pres_blm_bmi  <- y - X_splines[, -idx_bmi]  %*% pm_blm[-idx_bmi]
pres_blm_hdl  <- y - X_splines[, -idx_hdl]  %*% pm_blm[-idx_hdl]
pres_blm_gluc <- y - X_splines[, -idx_gluc] %*% pm_blm[-idx_gluc]
pres_blm_sbp  <- y - X_splines[, -idx_sbp]  %*% pm_blm[-idx_sbp]
pres_blm_pir  <- y - X_splines[, -idx_pir]  %*% pm_blm[-idx_pir]

# 3. Curve stimate (media a posteriori) + bande, in formato lungo
mk_fit <- function(x, fit, low, upp, var) {
  data.frame(x = as.vector(x), fit = as.vector(fit),
             low = as.vector(low), upp = as.vector(upp),
             var = var)
}

df_fit_blm <- rbind(
  mk_fit(x_age,  age_pm_blm,  age_low_blm,  age_upp_blm,  "Age"),
  mk_fit(x_bmi,  bmi_pm_blm,  bmi_low_blm,  bmi_upp_blm,  "BMI"),
  mk_fit(x_hdl,  hdl_pm_blm,  hdl_low_blm,  hdl_upp_blm,  "HDL"),
  mk_fit(x_gluc, gluc_pm_blm, gluc_low_blm, gluc_upp_blm, "Glucose"),
  mk_fit(x_sbp,  sbp_pm_blm,  sbp_low_blm,  sbp_upp_blm,  "SBP"),
  mk_fit(x_pir,  pir_pm_blm,  pir_low_blm,  pir_upp_blm,  "PIR")
)

# 4. Residui parziali, in formato lungo
mk_pts <- function(x, pres, var) {
  data.frame(x = as.vector(x), pres = as.vector(pres), var = var)
}

df_pts_blm <- rbind(
  mk_pts(px_age,  pres_blm_age,  "Age"),
  mk_pts(px_bmi,  pres_blm_bmi,  "BMI"),
  mk_pts(px_hdl,  pres_blm_hdl,  "HDL"),
  mk_pts(px_gluc, pres_blm_gluc, "Glucose"),
  mk_pts(px_sbp,  pres_blm_sbp,  "SBP"),
  mk_pts(px_pir,  pres_blm_pir,  "PIR")
)

df_fit_blm$var <- factor(df_fit_blm$var, levels = lev_v)
df_pts_blm$var <- factor(df_pts_blm$var, levels = lev_v)

# 5. Ritaglio delle bande al range dei residui parziali di ogni covariata
rng <- do.call(rbind, lapply(lev_v, function(v) {
  r <- range(df_pts_blm$pres[df_pts_blm$var == v])
  data.frame(var = factor(v, levels = lev_v), ymin = r[1], ymax = r[2])
}))

i <- match(df_fit_blm$var, rng$var)
df_fit_blm$low <- pmax(df_fit_blm$low, rng$ymin[i])
df_fit_blm$upp <- pmin(df_fit_blm$upp, rng$ymax[i])

# 6. Grafico
pl_blm <- ggplot(df_fit_blm, aes(x = x)) +
  geom_ribbon(aes(ymin = low, ymax = upp), fill = "orange", alpha = 0.3) +
  geom_line(aes(y = fit), colour = "darkgreen", linewidth = 0.8) +
  geom_point(aes(y = pres), data = df_pts_blm, size = 0.2, alpha = 0.3) +
  facet_wrap(~ var, nrow = 2, ncol = 3, scales = "free") +
  labs(x = NULL, y = "Partial effect") +
  theme_light() +
  theme(strip.text       = element_text(size = 11, colour = "black"),
        strip.background = element_rect(fill = "gray82"),
        panel.grid.major = element_line(linewidth = 0.3,  colour = "gray93"),
        panel.grid.minor = element_line(linewidth = 0.15, colour = "gray93"))

pl_blm
# ggsave(filename = "pl_blm.pdf", plot = pl_blm, width = 10, height = 6.5)


library(ggthemes)

# sostituisci con i nomi dei tuoi coefficienti
titoli <- c("Intercept", "Mexican American", "Other Hispanic", "Non-Hispanic Black", "Other race")

df_d <- do.call(rbind, lapply(1:5, function(j) {
  data.frame(rls   = as.numeric(beta_ps[, j]),
             blm   = as.numeric(post_samp[, j]),
             title = titoli[j])
}))
df_d$title <- factor(df_d$title, levels = titoli)  # mantiene l'ordine dei pannelli

pl_d <- ggplot(df_d) + 
  facet_wrap(~ title, nrow = 2, ncol = 3, scales = "free") +
  stat_density(aes(x = blm, col = "blm"),
               geom = "line", position = "identity", adjust = 1.5) +
  stat_density(aes(x = rls, col = "rls"),
               geom = "line", position = "identity", adjust = 1.5) +
  scale_color_tableau(name = "", palette = "Color Blind",
                      labels = c("rls" = "Predictive Bayes", 
                                 "blm" = "Traditional Bayes")) +
  xlab("") + ylab("Density") +
  theme_light() +
  theme(legend.text.align = 0,
        strip.text = element_text(size = 11, colour = "black"),
        strip.background = element_rect(fill = "gray82"),
        panel.grid.major = element_line(size = 0.3, colour = "gray93"),
        panel.grid.minor = element_line(size = 0.15, colour = "gray93"),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 7),
        axis.title = element_text(size = 10),
        legend.title = element_text(size = 10),
        plot.title = element_text(hjust = 0.5),
        legend.position = "top")
pl_d
# ggsave(filename = "pl_d.pdf", plot = pl_d, width = 10, height = 6.5)





# corrispondenze: nome pannello -> colonna di df_mgcv, termine smooth, griglia, indice di scala
var_col  <- c(Age = "age", BMI = "bmi", HDL = "hdl",
              Glucose = "glucose", SBP = "sbp", PIR = "pir")
term_lab <- paste0("s(", var_col, ")"); names(term_lab) <- names(var_col)
grid     <- list(Age = age_new, BMI = bmi_new, HDL = hdl_new,
                 Glucose = gluc_new, SBP = sbp_new, PIR = pir_new)
sc_idx   <- c(Age = 1, BMI = 2, HDL = 3, Glucose = 4, SBP = 5, PIR = 6)

z <- qnorm(0.975)   # bande puntuali al 95%

# 1. Curve stimate + bande sulla griglia (scala standardizzata -> originale)
#    type = "terms" restituisce solo il contributo dello smooth richiesto,
#    quindi i valori delle altre covariate in newdata sono irrilevanti
df_fit_gam <- do.call(rbind, lapply(lev_v, function(v) {
  g  <- as.vector(grid[[v]])
  nd <- df_mgcv[rep(1, length(g)), , drop = FALSE]
  nd[[var_col[v]]] <- g
  pr  <- predict(fit_t, newdata = nd, type = "terms",
                 terms = term_lab[v], se.fit = TRUE)
  fit <- as.vector(pr$fit[, term_lab[v]])
  se  <- as.vector(pr$se.fit[, term_lab[v]])
  data.frame(x   = g * train_scale[sc_idx[v]] + train_loc[sc_idx[v]],
             fit = fit, low = fit - z * se, upp = fit + z * se,
             var = v)
}))

# 2. Residui parziali: residuo di risposta + contributo dello smooth
#    (link identita': equivale a y - tutti gli altri termini, come per BLM)
pt     <- predict(fit_t, type = "terms")
r_resp <- residuals(fit_t, type = "response")

df_pts_gam <- do.call(rbind, lapply(lev_v, function(v) {
  data.frame(x    = df_mgcv[[var_col[v]]] * train_scale[sc_idx[v]] + train_loc[sc_idx[v]],
             pres = as.vector(r_resp + pt[, term_lab[v]]),
             var  = v)
}))

df_fit_gam$var <- factor(df_fit_gam$var, levels = lev_v)
df_pts_gam$var <- factor(df_pts_gam$var, levels = lev_v)

# 3. Ritaglio delle bande al range dei residui parziali di ogni covariata
rng <- do.call(rbind, lapply(lev_v, function(v) {
  r <- range(df_pts_gam$pres[df_pts_gam$var == v])
  data.frame(var = factor(v, levels = lev_v), ymin = r[1], ymax = r[2])
}))

i <- match(df_fit_gam$var, rng$var)
df_fit_gam$low <- pmax(df_fit_gam$low, rng$ymin[i])
df_fit_gam$upp <- pmin(df_fit_gam$upp, rng$ymax[i])

# 4. Grafico
pl_gam <- ggplot(df_fit_gam, aes(x = x)) +
  geom_ribbon(aes(ymin = low, ymax = upp), fill = "pink3", alpha = 0.2) +
  geom_line(aes(y = fit), colour = "purple4", linewidth = 0.8) +
  geom_point(aes(y = pres), data = df_pts_gam, size = 0.2, alpha = 0.3) +
  facet_wrap(~ var, nrow = 2, ncol = 3, scales = "free") +
  labs(x = NULL, y = "Partial effect") +
  theme_light() +
  theme(strip.text       = element_text(size = 11, colour = "black"),
        strip.background = element_rect(fill = "gray82"),
        panel.grid.major = element_line(linewidth = 0.3,  colour = "gray93"),
        panel.grid.minor = element_line(linewidth = 0.15, colour = "gray93"))

pl_gam
# ggsave(filename = "pl_mgcv.pdf", plot = pl_gam, width = 10, height = 6.5)




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









# Section 5.2 (Figure 1) and Section S2.1 of the supplement: illustration with
# skew Student-t errors.
# Run from this folder. The errors are centered skew-t draws, and the skew-t P_n is
# centered as well, so that its mean is x' beta (Assumption 1 of the paper).

rm(list = ls())

library(sn)
library(ggplot2)
library(ggthemes)

Rcpp::sourceCpp("../RLS_functions.cpp")
source("../RLS_functions.R")

# Settings
p <- 5
n_obs <- 15
n_all <- 1e5
sigma_true <- 0.2   # scale of the errors
alpha_true <- 3     # skewness parameter
nu_true <- 3        # degrees of freedom
B_ps <- 10000
b_0 <- rep(0, p)
V_0 <- diag(p) * 1000

# Mean of the standard skew-t ST(0, 1, alpha, nu), nu > 1
st_mean <- function(alpha, nu) {
  d <- alpha / sqrt(1 + alpha^2)
  d * sqrt(nu / pi) * gamma((nu - 1) / 2) / gamma(nu / 2)
}

# Data: same covariates and coefficients as student_t_errors.R
set.seed(1997); beta_true <- rnorm(p, sd = 2)
set.seed(2)
X_all <- matrix(rnorm(n_all * p), nrow = n_all, ncol = p)
X <- X_all[1:n_obs, ]
X_new <- X_all[(n_obs + 1):n_all, ]
mu_true <- c(X %*% beta_true)

set.seed(133)
eps <- rst(n_obs, 0, 1, alpha_true, nu_true) - st_mean(alpha_true, nu_true)
y <- mu_true + sigma_true * eps

# Fixed covariates for predictive resampling
set.seed(5)
X_ps <- X[sample.int(n_obs, 9 * n_obs, replace = TRUE), ]   # N = 10 n, with the correction

# RLS
bb <- rls(y, X, b_0, V_0)

# Skew Student-t P_n: par = c(log(sigma), alpha, log(nu - 1))
cat("Selecting the hyperparameters of the skew-t P_n by LOO-PLL\n")
opt <- nlminb(start = c(log(0.3), 1, log(2)),
              objective = function(x) npll_loo_lm_st_alt(x, y, X, b_0, V_0))
sigma_opt <- exp(opt$par[1])
alpha_opt <- opt$par[2]
df_opt <- exp(opt$par[3]) + 1
cat("sigma =", sigma_opt, " alpha =", alpha_opt, " nu =", df_opt, "\n")
beta_ps <- ps_lm_fixed_par_st_alt(bb$beta, bb$P, X_ps, B_ps, sigma_opt, alpha_opt, df_opt,
                                  nthreads = 7, correct = TRUE)

# Gaussian P_n (Student-t with df = 1e5, as in student_t_errors.R)
cat("Selecting the scale of the Gaussian P_n by LOO-PLL\n")
opt_n <- nlminb(start = 1, objective = function(x) npll_loo_lm_t_alt(c(log(1e5), x), y, X, b_0, V_0))
sigma_opt_n <- exp(opt_n$par)
cat("sigma =", sigma_opt_n, "\n")
beta_ps_n <- ps_lm_fixed_par_t_alt(bb$beta, bb$P, X_ps, B_ps, sigma_opt_n, 1e5, nthreads = 7, correct = TRUE)

# Draws of the 4th coefficient, for the figure of Section 5.2 (posteriors_figure.R)
dir.create("output", showWarnings = FALSE)
draws <- data.frame(fit = beta_ps[, 4], gauss = beta_ps_n[, 4])
save(draws, file = file.path("output", "draws_skew_t.RData"))

# Figure
lab_distr <- c("fit" = "Skew Student-t", "gauss" = "Gaussian", "true" = "True")
lty <- c("fit" = "solid", "gauss" = "dashed", "true" = "dotdash")
cols <- c("fit" = "#fc7d0b", "gauss" = "#1170aa", "true" = "#a3acb9")
theme_ill <- theme_light(base_size = 9) +
  theme(legend.text.align = 0,
        strip.text = element_text(colour = "black"),
        strip.background = element_rect(fill = "gray82"),
        panel.grid.major = element_line(linewidth = 0.3, colour = "gray93"),
        panel.grid.minor = element_line(linewidth = 0.15, colour = "gray93"),
        axis.text.y = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "top")

# Data generating distribution and fitted P_n, both centered at zero
y_val <- seq(-1.2, 2.4, length = 2000)
df_plot <- data.frame(
  y     = rep(y_val, 3),
  dens  = c(dst(y_val, xi = -sigma_opt * st_mean(alpha_opt, df_opt), omega = sigma_opt,
                alpha = alpha_opt, nu = df_opt),
            dnorm(y_val, sd = sigma_opt_n),
            dst(y_val, xi = -sigma_true * st_mean(alpha_true, nu_true), omega = sigma_true,
                alpha = alpha_true, nu = nu_true)),
  distr = factor(rep(c("fit", "gauss", "true"), each = length(y_val)), levels = names(lab_distr)),
  title = "Data generating distribution")
pl_p <- ggplot(df_plot, aes(x = y, y = dens, colour = distr, linetype = distr)) +
  facet_wrap(~ title) + geom_line(linewidth = 0.5) +
  scale_colour_manual(name = "", labels = lab_distr, values = cols) +
  scale_linetype_manual(name = "", labels = lab_distr, values = lty) +
  labs(x = "Y", y = "Density") +
  theme_ill

# Posterior distribution of the 4th coefficient
idx <- 4
xlim_post <- range(sapply(list(beta_ps, beta_ps_n),
                          function(b) quantile(b[, idx], c(.003, .997))))
df_d1 <- data.frame(
  beta  = c(beta_ps[, idx], beta_ps_n[, idx]),
  distr = factor(rep(c("fit", "gauss"), each = B_ps), levels = names(lab_distr)),
  title = "Posterior distribution")
pl <- ggplot(df_d1, aes(x = beta, colour = distr, linetype = distr)) +
  facet_wrap(~ title) +
  stat_density(geom = "line", position = "identity", adjust = 1.8, linewidth = 0.5) +
  scale_colour_manual(name = "", labels = lab_distr, values = cols) +
  scale_linetype_manual(name = "", labels = lab_distr, values = lty) +
  coord_cartesian(xlim = xlim_post) +
  xlab("Beta") + ylab("Density") +
  theme_ill

pl_tot <- ggpubr::ggarrange(pl_p, pl)
# Fitted predictive distributions and posterior: diagnostic figure, not in the paper
ggsave(file.path("output", "pl_ill_skewt.pdf"), plot = pl_tot, width = 6.3, height = 3.5, device = cairo_pdf)

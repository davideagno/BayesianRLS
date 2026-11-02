# Section S2.1 of the supplement: illustration with skew Student-t errors
# (Figure S2). Preliminary: Figure S2 has not been regenerated with this script.

rm(list = ls())

library(ggplot2)
library(ggthemes)

Rcpp::sourceCpp("../RLS_functions.cpp")
source("../RLS_functions.R")

# Centered skew Student-t (Azzalini and Capitanio, 2003). If T ~ ST(0, 1, alpha, nu),
# P_n is the law of x' beta + phi_i * (T - mu), with mu = E(T), so that its mean is
# x' beta (Assumption 1). The data are generated with the same centered errors.

# Mean and variance of T ~ ST(0, 1, alpha, nu), nu > 2
skt_moments <- function(nu, alpha) {
  delta <- alpha / sqrt(1 + alpha^2)
  b_nu <- sqrt(nu / pi) * exp(lgamma((nu - 1) / 2) - lgamma(nu / 2))
  mu <- b_nu * delta
  list(mu = mu, var = nu / (nu - 2) - mu^2)
}

# Density of T - mu
dskt_c <- function(z, nu, alpha, log = FALSE) {
  t <- z + skt_moments(nu, alpha)$mu
  out <- log(2) + dt(t, df = nu, log = TRUE) +
    pt(alpha * t * sqrt((nu + 1) / (nu + t^2)), df = nu + 1, log.p = TRUE)
  if (log) out else exp(out)
}

# Draws of T - mu, via T = Z / sqrt(W / nu) with Z skew-normal and W ~ chi^2_nu
rskt_c <- function(n, nu, alpha) {
  delta <- alpha / sqrt(1 + alpha^2)
  z <- delta * abs(rnorm(n)) + sqrt(1 - delta^2) * rnorm(n)
  z / sqrt(rchisq(n, df = nu) / nu) - skt_moments(nu, alpha)$mu
}

# Negative LOO-PLL, as npll_loo_lm_t_alt (to which it reduces for alpha = 0)
# par = c(log(nu), log(phi), alpha)
npll_loo_lm_skt <- function(par, y, X, beta_0, V_0) {
  nu <- exp(par[1])
  phi <- exp(par[2])
  alpha <- par[3]
  if (!is.finite(nu) || nu <= 2 || !is.finite(phi)) return(1e100)
  V_n <- solve(solve(V_0) + crossprod(X))
  beta_n <- V_n %*% (solve(V_0, beta_0) + crossprod(X, y))
  h <- rowSums((X %*% V_n) * X)
  r <- c(y - X %*% beta_n) / (1 - h)
  phi_i <- phi / sqrt(1 - h)
  out <- sum(dskt_c(r / phi_i, nu, alpha, log = TRUE) - log(phi_i))
  if (!is.finite(out)) return(1e100)
  -out
}

# Predictive resampling with shared future covariates (Algorithm 2), as
# ps_lm_fixed_par_t_alt: draws beta_n + sum_i w_i * xi_i, w_i = phi_i a_i V_{i-1} x_i,
# with phi_i = phi * sqrt(1 + x_i' V_{i-1} x_i)
ps_lm_fixed_skt <- function(beta_n, V_n, X_future, B, phi, nu, alpha, chunk = 100) {
  N <- nrow(X_future)
  p <- ncol(X_future)
  W <- matrix(0, p, N)
  V <- V_n
  for (i in 1:N) {
    x <- X_future[i, ]
    Vx <- c(V %*% x)
    q <- 1 + sum(x * Vx)
    W[, i] <- phi * sqrt(q) / q * Vx
    V <- V - tcrossprod(Vx) / q
  }
  out <- matrix(NA_real_, B, p)
  for (start in seq(1, B, by = chunk)) {
    idx <- start:min(B, start + chunk - 1)
    xi <- matrix(rskt_c(N * length(idx), nu, alpha), N, length(idx))
    out[idx, ] <- t(c(beta_n) + W %*% xi)
  }
  out
}

p <- 5
n_obs <- 15
n_all <- 1e5
sigma_true <- 0.2
df_true <- 3
alpha_true <- 4

set.seed(1997); beta_true <- rnorm(p, sd = 2)
set.seed(2)
X_all <- matrix(rnorm(n_all*p), nrow = n_all, ncol = p)
X <- X_all[1:n_obs,]
X_new <- X_all[(n_obs+1):n_all,]
mu_true <- c(X %*% beta_true)
eps <- rskt_c(n_obs, nu = df_true, alpha = alpha_true)
y <- mu_true + sigma_true*eps

# RLS
b_0 <- rep(0, p)
V_0 <- diag(p)*1000
bb <- rls(y, X, b_0, V_0)

# Fixed covariates for predictive resamplings
set.seed(5)
X_ps <- matrix(nrow = NROW(X_new), ncol = p)
for (i in 1:NROW(X_ps)) {
  idx <- sample(1:n_obs, size = 1)
  X_ps[i,] <- X[idx,]
}

# Skew Student-t P_n
opt <- nlminb(start = c(2, 2, 0), objective = npll_loo_lm_skt,
              y = y, X = X, beta_0 = b_0, V_0 = V_0)
opt
(df_opt <- exp(opt$par[1]))
(sigma_opt <- exp(opt$par[2]))
(alpha_opt <- opt$par[3])
B_ps <- 10000
set.seed(10)
beta_ps <- ps_lm_fixed_skt(bb$beta, bb$P, X_ps, B_ps, sigma_opt, df_opt, alpha_opt)

# Gaussian P_n
opt_n <- nlminb(start = 1, objective = function(x) npll_loo_lm_t_alt(c(log(1e5), x), y, X, b_0, V_0))
(sigma_opt_n <- exp(opt_n$par))
beta_ps_n <-  ps_lm_fixed_par_t_alt(bb$beta, bb$P, X_ps, B_ps, sigma_opt_n, 1e5, nthreads = 7)

# Posterior
idx <- 4
df_d1 <- data.frame(rls = as.numeric(beta_ps[,idx]),
                    blm = as.numeric(beta_ps_n[,idx]),
                    title = "Posterior distribution")
pl <- ggplot(df_d1) +
  facet_wrap(~ title) +
  stat_density(aes(x = rls, col = "rls"),
               geom = "line", position = "identity", adjust = 1.8) +
  stat_density(aes(x = blm, col = "blm"),
               geom = "line", position = "identity", adjust = 1.8, linetype = "dashed") +
  scale_color_tableau( name = "", palette = "Color Blind",
                       labels = c("rls" = "Skew Student-t",
                                  "blm" = "Gaussian")) +
  xlab("Beta") + ylab("Density") +
  coord_cartesian(xlim = quantile(c(df_d1$rls, df_d1$blm), c(0.001, 0.999))) +
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
pl

# Plug-in predictive
y_val <- seq(-1.2, 2.4, length = 2000)
df_plot <- data.frame(
  y     = rep(y_val, 3),
  dens  = c(dnorm(y_val, sd = sigma_opt_n),
            dskt_c(y_val / sigma_opt, df_opt, alpha_opt) / sigma_opt,
            dskt_c(y_val / sigma_true, df_true, alpha_true) / sigma_true),
  distr = factor(rep(c("rls_n", "rls_t", "true"), each = length(y_val)),
                 levels = c("rls_n", "rls_t", "true")),
  title = "Data generating distribution"
)
lab_distr <- c("rls_t" = "Skew Student-t", "rls_n" = "Gaussian", "true" = "True")
pl_p <- ggplot(df_plot, aes(x = y, y = dens, colour = distr, linetype = distr)) +
  facet_wrap(~ title) +
  geom_line(linewidth = 0.5) +
  scale_color_tableau(name = "", palette = "Color Blind", labels = lab_distr) +
  scale_linetype_manual(name = "", labels = lab_distr,
                        values = c("rls_t" = "solid", "rls_n" = "dashed", "true" = "dotdash")) +
  labs(x = "Y", y = "Density") +
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
pl_p

pl_tot <- gridExtra::arrangeGrob(pl_p, pl, ncol = 2)
plot(pl_tot)
dir.create("output", showWarnings = FALSE)
ggsave(file.path("output", "pl_ill_skewt.pdf"), plot = pl_tot, width = 9, height = 5)

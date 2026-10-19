
generate_am_t <- function(n = 500, p = 4, df = 10, sigma = 1, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  # Generate covariates
  X <- matrix(runif(n * p, -2, 2), nrow = n, ncol = p)
  colnames(X) <- paste0("x", 1:p)
  
  # Define component functions
  f <- list(
    function(x) 0.5 * x^3 - x,              # cubic
    function(x) 3 * x * exp(-x^2),          # bump
    function(x) x^2 - 1,                    # quadratic
    function(x) 2 * sin(pi * x)             # smooth sinusoidal
  )
  
  # Evaluate each component
  F_mat <- matrix(0, nrow = n, ncol = p)
  for (j in 1:p) {
    fj <- f[[(j - 1) %% length(f) + 1]]  # recycle if p > length(f)
    F_mat[, j] <- fj(X[, j])
  }
  colnames(F_mat) <- paste0("f", 1:p)
  
  # Center each component (identifiability)
  F_mat <- scale(F_mat, center = TRUE, scale = FALSE)
  
  # Intercept + noise
  mu <- 1.0
  # eps <- rnorm(n, 0, sigma)
  eps <- rt(n, df)
  y <- mu + rowSums(F_mat) + sigma*eps
  
  list(
    y = y,
    X = X,
    F_true = F_mat,
    mu = mu,
    df = df,
    n = n,
    p = p,
    true_functions = f
  )
}

true_f_test <- function(obj, Xnew) {
  out <- matrix(nrow = NROW(Xnew), ncol = NCOL(Xnew))
  
  loc <- attr(obj$F_true, "scaled:center")
  f_list <- obj$true_functions
  
  for (j in 1:NCOL(Xnew)) {
    fj <- f_list[[(j - 1) %% length(f_list) + 1]]
    out[,j] <- fj(Xnew[,j])
    out[,j] <- out[,j] - loc[j]
  }
  colnames(out) <- paste0("f", 1:NCOL(Xnew))
  out
}

permute_data <- function(y, X, nPerm, seed = 1) {
  
  n_obs <- NROW(X)
  y_perm <- X_perm <- idx_perm <- list()
  set.seed(seed)
  for (k in 1:nPerm) {
    idx_perm[[k]] <- sample(1:n_obs, size = n_obs, replace = FALSE)
    y_perm[[k]] <- y[idx_perm[[k]]]
    X_perm[[k]] <- X[idx_perm[[k]],]
  }
  
  return(list(y = y_perm, X = X_perm, idx = idx_perm))
}
permute_idx <- function(n, nPerm, seed = 1) {
  
  n
  idx_perm <- list()
  set.seed(seed)
  for (k in 1:nPerm) {
    idx_perm[[k]] <- sample(1:n, size = n, replace = FALSE)
  }
  
  return(idx_perm)
}

est_rs <- function(X, y, beta_0, P_0, nPerm = 1, seed = 1) {
  require(mgcv)
  require(parallel)
  
  set.seed(seed)
  
  # Permutations of observations
  perm <- permute_idx(length(y), nPerm, seed)
  y_perm <- X_perm <- list()
  for (k in 1:nPerm) {
    idx_perm <- perm[[k]]
    y_perm[[k]] <- y[idx_perm]
    X_perm[[k]] <- list()
    for (j in 1:length(X)) {
      X_perm[[k]][[j]] <- X[[j]][idx_perm,]
    }
  }
  
  mclapply_function <- function(l) {
    set.seed(seed*l)
    out_par <- est_rs_cpp(X_perm[[l]], y_perm[[l]], beta_0, P_0)
    out_par
  }
  out <- mclapply(1:nPerm, mclapply_function)
  
  beta_out <- P_out <- list()
  intercept_out <- 0
  for (j in 1:length(X)) {
    dim_j <- length(out[[1]]$beta[[j]])
    beta_out[[j]] <- rep(0, dim_j)
    P_out[[j]] <- matrix(0, ncol = dim_j, nrow = dim_j)
    for (k in 1:nPerm) {
      beta_out[[j]] <- beta_out[[j]] + out[[k]]$beta[[j]]
      P_out[[j]] <- P_out[[j]] + out[[k]]$P[[j]]
      intercept_out <- intercept_out + out[[k]]$intercept
    }
    beta_out[[j]] <- beta_out[[j]] / nPerm
    P_out[[j]] <- P_out[[j]] / nPerm
    intercept_out <- intercept_out / nPerm
  }
  
  return(list(y = y, X = X, beta_0 = beta_0, P_0 = P_0, seed = seed,
              beta = beta_out, P = P_out, intercept = intercept_out))
}
est_rs_lambda <- function(y, X, beta_0, lambda, Omega, nPerm = 1, seed = 1) {
  require(mgcv)
  require(parallel)
  
  set.seed(seed)
  
  P_0 <- invert_lambda_Omega(lambda, Omega)
  
  # Permutations of observations
  perm <- permute_idx(length(y), nPerm, seed)
  y_perm <- X_perm <- list()
  for (k in 1:nPerm) {
    idx_perm <- perm[[k]]
    y_perm[[k]] <- y[idx_perm]
    X_perm[[k]] <- list()
    for (j in 1:length(X)) {
      X_perm[[k]][[j]] <- X[[j]][idx_perm,]
    }
  }
  
  mclapply_function <- function(l) {
    set.seed(seed*l)
    out_par <- est_rs_cpp(y_perm[[l]], X_perm[[l]], beta_0, P_0)
    out_par
  }
  out <- mclapply(1:nPerm, mclapply_function)
  
  beta_out <- P_out <- list()
  intercept_out <- 0
  for (j in 1:length(X)) {
    dim_j <- length(out[[1]]$beta[[j]])
    beta_out[[j]] <- rep(0, dim_j)
    P_out[[j]] <- matrix(0, ncol = dim_j, nrow = dim_j)
    for (k in 1:nPerm) {
      beta_out[[j]] <- beta_out[[j]] + out[[k]]$beta[[j]]
      P_out[[j]] <- P_out[[j]] + out[[k]]$P[[j]]
      intercept_out <- intercept_out + out[[k]]$intercept
    }
    beta_out[[j]] <- beta_out[[j]] / nPerm
    P_out[[j]] <- P_out[[j]] / nPerm
    intercept_out <- intercept_out / nPerm
  }
  
  return(list(y = y, X = X, beta_0 = beta_0, lambda = lambda, Omega = Omega, seed = seed,
              beta = beta_out, P = P_out, intercept = intercept_out))
}

pll_rs_t <- function(df, sigma, X, y, beta_0, P_0, nPerm = 1, seed = 1) {
  require(mgcv)
  require(parallel)
  
  set.seed(seed)
  
  # Permutations of observations
  perm <- permute_idx(length(y), nPerm, seed)
  y_perm <- X_perm <- list()
  for (k in 1:nPerm) {
    idx_perm <- perm[[k]]
    y_perm[[k]] <- y[idx_perm]
    X_perm[[k]] <- list()
    for (j in 1:length(X)) {
      X_perm[[k]][[j]] <- X[[j]][idx_perm,]
    }
  }
  
  mclapply_function <- function(l) {
    set.seed(seed*l)
    out_par <- pll_rs_t_cpp(df, sigma, X_perm[[l]], y_perm[[l]], beta_0, P_0)
    out_par
  }
  out <- mclapply(1:nPerm, mclapply_function)
  
  return(mean(unlist(out)))
}
pll_rs_t_lambda <- function(df, sigma, lambda, y, X, beta_0, Omega, nPerm = 1, seed = 1) {
  require(mgcv)
  require(parallel)
  
  set.seed(seed)
  
  # Permutations of observations
  perm <- permute_idx(length(y), nPerm, seed)
  y_perm <- X_perm <- list()
  for (k in 1:nPerm) {
    idx_perm <- perm[[k]]
    y_perm[[k]] <- y[idx_perm]
    X_perm[[k]] <- list()
    for (j in 1:length(X)) {
      X_perm[[k]][[j]] <- X[[j]][idx_perm,]
    }
  }
  P_0 <- invert_lambda_Omega(lambda, Omega)
  
  mclapply_function <- function(l) {
    set.seed(seed*l)
    out_par <- pll_rs_t_cpp(df, sigma, X_perm[[l]], y_perm[[l]], beta_0, P_0)
    out_par
  }
  out <- mclapply(1:nPerm, mclapply_function)
  
  return(mean(unlist(out)))
}

# funzioni con rls
ps_rs_conv_t_1 <- function(X, beta, P, N, B, df, sigma, seed = 1) {
  out <- matrix(nrow = N, ncol = B)
  n <- NROW(X)
  for (b in 1:B) {
    print(b)
    beta_b <- beta
    P_b <- P
    for (i in 1:N) {
      idx_new <- sample(1:n, size = 1)
      x_new <- X[idx_new,]
      y_new <- sum(x_new * beta_b) + rt(1, df)*sigma
      Px <- P_b %*% x_new
      a <- c(1 / (1 + crossprod(x_new, Px)))
      beta_b <- beta_b + Px * (a * (y_new - sum(x_new * beta_b)))
      P_b <- P_b - a * tcrossprod(Px)
      out[i,b] <- sum((unlist(beta) - unlist(beta_b))^2)
    }
  }
  out
}


# funzioni con rls
ps_rs_conv_t_1_alt <- function(X, beta, P, N, B, df, sigma, seed = 1) {
  out <- matrix(nrow = N, ncol = B)
  n <- NROW(X)
  for (b in 1:B) {
    print(b)
    beta_b <- beta
    P_b <- P
    for (i in 1:N) {
      idx_new <- sample(1:n, size = 1)
      x_new <- X[idx_new,]
      Px <- P_b %*% x_new
      a <- c(1 / (1 + crossprod(x_new, Px)))
      sigma_new <- sigma / sqrt(a)
      y_new <- sum(x_new * beta_b) + rt(1, df)*sigma_new
      beta_b <- beta_b + Px * (a * (y_new - sum(x_new * beta_b)))
      P_b <- P_b - a * tcrossprod(Px)
      out[i,b] <- sum((unlist(beta) - unlist(beta_b))^2)
    }
  }
  out
}


ps_rs_parallel_t_1 <- function(X, beta, P, N, B, df, sigma, seed = 1, 
                               mc.cores = parallel::detectCores() - 1) {
  
  # Check if parallel package is available
  if (!requireNamespace("parallel", quietly = TRUE)) {
    stop("Package 'parallel' is required but not installed.")
  }
  
  # On Windows, mclapply doesn't work, so fall back to lapply
  if (.Platform$OS.type == "windows") {
    warning("mclapply not available on Windows. Using sequential lapply instead.")
    mc.cores <- 1
  }
  
  p <- length(beta)
  
  set.seed(seed)
  # Function to run a single iteration with error handling
  run_iteration <- function(b) {
    tryCatch({
      # Each iteration gets a different seed
      iteration_seed <- seed + b - 1
      
      # Run single iteration
      beta_result <- ps_rs_single_iteration_t_1(
        X = X,
        beta = beta,
        P = P,
        N = N,
        sigma = sigma,
        df = df
      )
      
      return(list(success = TRUE, result = beta_result, iteration = b))
    }, error = function(e) {
      return(list(success = FALSE, error = as.character(e), iteration = b))
    })
  }
  
  # Test single iteration first
  message("Testing single iteration...")
  test_result <- run_iteration(1)
  if (!test_result$success) {
    stop("Single iteration test failed with error: ", test_result$error)
  }
  message("Single iteration test successful!")
  
  # Parallel execution
  message(sprintf("Running %d iterations on %d cores...", B, mc.cores))
  
  results <- parallel::mclapply(
    X = 1:B,
    FUN = run_iteration,
    mc.cores = mc.cores,
    mc.preschedule = TRUE
  )
  
  # Check for errors
  errors <- sapply(results, function(x) !x$success)
  if (any(errors)) {
    error_iterations <- which(errors)
    error_messages <- sapply(results[errors], function(x) x$error)
    stop(sprintf("Errors in iterations: %s\nFirst error: %s", 
                 paste(error_iterations, collapse = ", "),
                 error_messages[1]))
  }
  
  # Reorganize results into beta_ps format
  beta_ps <- vector("list", p)
  for (j in 1:p) {
    beta_ps[[j]] <- matrix(nrow = B, ncol = length(beta[[j]]))
    for (b in 1:B) {
      beta_ps[[j]][b, ] <- results[[b]]$result[[j]]
    }
  }
  
  # Return results in same format as original function
  return(list(
    beta_ps = beta_ps,
    X = X,
    beta = beta,
    P = P,
    B = B,
    N = N,
    sigma = sigma,
    df = df,
    seed = seed
  ))
}

ps_single_iteration_t_1 <- function(X, beta, P, N, sigma, df) {
  out <- beta
  P_ps <- P
  n <- nrow(X)
  idxs <- sample.int(n, size = N, replace = TRUE)
  t_vals <- rt(N, df)
  for (i in seq_len(N)) {
    x_new <- X[idxs[i], ]
    y_new <- sum(x_new * out) + t_vals[i] * sigma
    Px <- P_ps %*% x_new
    a <- c(1 / (1 + crossprod(x_new, Px)))
    out <- out + Px * (a * (y_new - sum(x_new * out)))
    P_ps <- P_ps - a * tcrossprod(Px)
  }
  out
}



estimate_uc_real_rate <- function(
    rate = c("1y", "10y"),
    N = 5,
    M = 20000,
    seed = 616,
    save_results = TRUE,
    results_dir = "results",
    verbose = TRUE
) {
  rate <- match.arg(rate)
  set.seed(seed)
  
  if (rate == "1y") {
    if (verbose) cat("Estimation for 1y real rate:\n")
    data <- loadShortRate()
    filename <- file.path(results_dir, "uc_estimates_1y.RData")
  } else {
    if (verbose) cat("Estimation for 10y real rate:\n")
    data <- loadLongRate()
    filename <- file.path(results_dir, "uc_estimates_10y.RData")
  }
  
  if (verbose) {
    cat("Annual sample starts in", min(data$year), "and ends in", max(data$year), "\n")
    cat("Filename for estimation results:", filename, "\n")
  }
  
  prior <- get_prior()
  
  T <- length(data$rr)
  phi <- matrix(NA_real_, N, M)
  sigv2 <- matrix(NA_real_, N, M)
  sigeta2 <- matrix(NA_real_, N, M)
  x <- array(NA_real_, c(N, M, 2, T))
  
  start <- Sys.time()
  
  for (i in seq_len(N)) {
    pars <- getRandomPars()
    
    if (verbose) {
      cat("Chain", i, "- initialized at:\n")
      cat("phi =", pars$phi, "sigv2 =", pars$sigv2, "sigeta2 =", pars$sigeta2, "\n")
    }
    
    for (j in seq_len(M)) {
      if (verbose && j %% 5000 == 0) {
        cat("iteration", j, "\n")
      }
      
      ## draw latent state variables
      mod <- ssm_uc(data$rr, pars, prior)
      sim <- KFAS::simulateSSM(mod, type = "states")
      pars$x <- t(sim[, , 1])
      
      ## AR(1) parameter
      pars$phi <- draw_regression(
        pars$x[2, -1],
        pars$x[2, -T],
        pars$sigv2,
        0,
        1 / prior$var_phi,
        restrict = TRUE
      )
      
      ## random walk innovation variance
      alpha1 <- prior$alpha_eta + (T - 2)
      delta1 <- prior$delta_eta + sum(diff(pars$x[1, ])^2)
      pars$sigeta2 <- 1 / rgamma(1, alpha1 / 2, delta1 / 2)
      
      ## AR(1) innovation variance
      G <- embed(pars$x[2, ], 2)
      ydat <- G[, 1]
      xdat <- G[, -1, drop = FALSE]
      resid <- ydat - xdat * pars$phi
      alpha1 <- prior$alpha_v + (T - 2)
      delta1 <- prior$delta_v + sum(resid^2)
      pars$sigv2 <- 1 / rgamma(1, alpha1 / 2, delta1 / 2)
      
      ## save draws
      phi[i, j] <- pars$phi
      sigv2[i, j] <- pars$sigv2
      sigeta2[i, j] <- pars$sigeta2
      x[i, j, , ] <- pars$x[1:2, ]
    }
    
    if (verbose) {
      cat("Duration so far:", format(Sys.time() - start), "\n")
    }
  }
  
  if (verbose) {
    cat(N * M, "iterations take\n")
    print(Sys.time() - start)
  }
  
  ## save full MCMC sample
  mcmc_save <- list(phi = phi, sigv2 = sigv2, sigeta2 = sigeta2)
  parnames <- c("phi", "sigv2", "sigeta2")
  n <- length(parnames)
  
  ## flatten multiple chains into one chain and drop burn-in sample (first half)
  ind <- rep((M / 2 + 1):M, N) + M * rep(0:(N - 1), each = M / 2)
  
  phi_vec <- as.numeric(mcmc_save$phi)[ind]
  sigv2_vec <- as.numeric(mcmc_save$sigv2)[ind]
  sigeta2_vec <- as.numeric(mcmc_save$sigeta2)[ind]
  X <- array(x, c(N * M, 2, T))[ind, , ]
  
  stopifnot(all(abs(phi_vec) < 1))
  
  lambda <- sqrt(sigeta2_vec / sigv2_vec) * abs(1 - phi_vec)
  
  if (verbose) {
    print(summary(lambda))
  }
  
  ## parameter estimates
  tbl <- matrix(NA_real_, n + 1, 4)
  rownames(tbl) <- c(parnames, "lambda (signal-to-noise)")
  colnames(tbl) <- c("MCMC mean", "median", "LB", "UB")
  
  theta <- cbind(phi_vec, sigv2_vec, sigeta2_vec, lambda)
  tbl[, 1] <- colMeans(theta, na.rm = TRUE)
  tbl[, 2:4] <- t(apply(theta, 2, quantile, c(.5, .025, .975), na.rm = TRUE))
  
  if (verbose) {
    print(round(tbl, 4))
    sigu <- sqrt(sigeta2_vec)
    print("Posterior distribution of sigma_u:")
    print(summary(sigu))
    print("Posterior distribution of phi:")
    print(summary(phi_vec))
    cat("# posterior correlations\n")
    print(cor(theta))
  }
  
  ## MCMC estimate of rstar
  data$rstar.mean <- apply(X[, 1, ], 2, mean)
  data$rstar.lb <- apply(X[, 1, ], 2, quantile, .025)
  data$rstar.ub <- apply(X[, 1, ], 2, quantile, .975)
  
  if (save_results) {
    dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
    cat("Final save path:", filename, "\n")
    save(data, X, sigeta2_vec, sigv2_vec, phi_vec, prior, file = filename)
  }
  
  invisible(list(data = data, filename = filename))
}
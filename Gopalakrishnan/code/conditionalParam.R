conditionalParamsForX_custom <- function(xVec,         # length q
                                         piVec,        # length K
                                         muMat,        # d x K
                                         SigmaList) {  # list of (d x d), [[k]]
  K <- length(piVec)
  q <- nrow(muMat)-1
  stopifnot(length(xVec) == q)

  make_posdef <- function(S, eps = 1e-6) {
    S_sym <- (S + t(S)) / 2
    ev <- eigen(S_sym, symmetric = TRUE)
    vals <- pmax(ev$values, eps)
    ev$vectors %*% diag(vals, length(vals)) %*% t(ev$vectors)
  }

  # Get Sigma_k
  getSigma <- if (is.array(SigmaList)) {
    function(k) {
      S <- SigmaList[, , k, drop = FALSE]
      matrix(S, nrow = q+1, ncol = q+1)
    }
  } else if (is.list(SigmaList)) {
    function(k) SigmaList[[k]]
  } else {
    stop("SigmaList must be a 3‑D array or a list of covariance matrices.")
  }

  # # Test if the sigma_xx is positive definite -> TRUE/FALSE
  # is_posdef <- function(S) {
  #   ev <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
  #   all(ev > 1e-8)
  # }
  #
  # Posterior weights p(k | x) ∝ pi_k * N_q(x | mu_xk, Sigma_xxk)
  logNumer <- numeric(K)
  for (k in seq_len(K)) {
    Sigmak      <- getSigma(k)
    Sigma_xx <- Sigmak[1:q, 1:q, drop = FALSE]
    Sigma_xx_safe <- make_posdef(Sigma_xx)
    mu_xk     <- muMat[1:q, k]

    logdens <- mvtnorm::dmvnorm(xVec, mean = mu_xk, sigma = Sigma_xx_safe, log = TRUE)
    if (!is.finite(logdens) || !is.finite(piVec[k]) || piVec[k] <= 0) {
      logNumer[k] <- -Inf
    } else {
      logNumer[k] <- log(piVec[k]) + logdens
    }
  }

  # log-sum-exp for stability
  finite_idx <- is.finite(logNumer)
  if (!any(finite_idx)) {
    post_w <- rep(1 / K, K)
  } else {
    m <- max(logNumer[finite_idx])
    numerator <- rep(0, K)
    numerator[finite_idx] <- exp(logNumer[finite_idx] - m)
    numerator_sum <- sum(numerator)
    if (!is.finite(numerator_sum) || numerator_sum <= 0) {
      post_w <- rep(1 / K, K)
    } else {
      post_w <- numerator / numerator_sum
    }
  }


  # Conditional z | x, k parameters
  cond_means <- numeric(K)
  cond_vars  <- numeric(K)

  for (k in seq_len(K)) {
    Sigmak       <- getSigma(k)
    Sigma_xx <- Sigmak[1:q, 1:q, drop = FALSE]     # Σ_xx
    Sigma_xx_safe <- make_posdef(Sigma_xx)
    Sigma_xz <- Sigmak[1:q, q+1, drop = FALSE]     # Σ_xz (q×1)
    Sigma_zx <- t(Sigma_xz)                        # Σ_zx (1×q)
    sigma_zz <- Sigmak[q+1, q+1, drop = FALSE]     # σ_zz

    mu_x <- muMat[1:q, k]
    mu_z <- muMat[q+1, k]

    # Solve y = Σ_xx^{-1}
    R <- chol(Sigma_xx_safe + diag(1e-10, q))        # ridge = 1e-10, for stability
    solve_Sxx <- function(v) backsolve(R, backsolve(R, v, transpose = TRUE))

    delta_x      <- xVec - mu_x                      # (x - μ_x) (q)
    Sxx_inv_dx   <- solve_Sxx(delta_x)               # Σ_xx^{-1}(x - μ_x) (q) - for mu
    Sxx_inv_xz   <- solve_Sxx(Sigma_xz)              # Σ_xx^{-1}Σ_xz (q×1) - for var

    cond_means[k] <- as.numeric(mu_z + Sigma_zx %*% Sxx_inv_dx)
    cond_vars[k] <- max(as.numeric(sigma_zz - Sigma_zx %*% Sxx_inv_xz), 1e-8)
  }

  list(
    post_weights = post_w,       # p(k | x)
    cond_means   = cond_means,   # E[z | x, k]
    cond_vars    = cond_vars     # Var[z | x, k]
  )
}

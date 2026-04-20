conditionalParamsForX_custom <- function(xVec,         # length q
                                         piVec,        # length K
                                         muMat,        # d x K
                                         SigmaList) {  # list of (d x d), [[k]]
  K <- length(piVec)
  q <- nrow(muMat)-1
  stopifnot(length(xVec) == q)
  if (anyNA(piVec) || any(!is.finite(piVec)) || any(piVec < 0) || sum(piVec) <= 0) {
    stop("piVec must be finite, non-negative, and sum to a positive value.", call. = FALSE)
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

  stabilize_cov <- function(S, eps = 1e-8, max_steps = 8L) {
    # Symmetrize then add adaptive ridge until Cholesky succeeds.
    S <- (S + t(S)) / 2
    q_loc <- nrow(S)
    d <- diag(S)
    scale <- mean(abs(d[is.finite(d)]))
    if (!is.finite(scale) || scale <= 0) scale <- 1

    for (step in 0:max_steps) {
      jitter <- if (step == 0) 0 else eps * scale * (10 ^ (step - 1))
      S_try <- if (jitter > 0) S + diag(jitter, q_loc) else S
      if (!inherits(try(chol(S_try), silent = TRUE), "try-error")) return(S_try)
    }

    ev <- eigen(S, symmetric = TRUE)
    floor_val <- eps * scale
    vals <- pmax(ev$values, floor_val)
    ev$vectors %*% diag(vals, q_loc) %*% t(ev$vectors)
  }

  # Posterior weights p(k | x) ∝ pi_k * N_q(x | mu_xk, Sigma_xxk)
  logNumer <- numeric(K)
  for (k in seq_len(K)) {
    Sigmak      <- getSigma(k)
    Sigma_xx <- Sigmak[1:q, 1:q, drop = FALSE]
    mu_xk     <- muMat[1:q, k]
    Sigma_xx_pd <- stabilize_cov(Sigma_xx)
    logdens <- mvtnorm::dmvnorm(xVec, mean = mu_xk, sigma = Sigma_xx_pd, log = TRUE)
    logNumer[k] <- log(piVec[k]) + logdens
  }

  # Handle +Inf robustly: if one or more components dominate with +Inf log-prob,
  # assign all mass to them uniformly.
  pos_inf <- is.infinite(logNumer) & (logNumer > 0)
  if (any(pos_inf)) {
    post_w <- numeric(K)
    post_w[pos_inf] <- 1 / sum(pos_inf)
  } else {
    finite_log <- is.finite(logNumer)
    if (!any(finite_log)) {
      stop(
        "conditionalParamsForX_custom: all logNumer are non-finite. ",
        "xVec = c(", paste(format(xVec, digits = 16, scientific = TRUE), collapse = ", "),
        "), piVec = c(", paste(format(piVec, digits = 16, scientific = TRUE), collapse = ", "),
        "), logNumer = c(", paste(format(logNumer, digits = 16, scientific = TRUE), collapse = ", "), ").",
        call. = FALSE
      )
    }

    m <- max(logNumer[finite_log])
    numerator <- exp(logNumer - m)
    numerator[!is.finite(numerator)] <- 0
    total <- sum(numerator)
    if (!is.finite(total) || total <= 0) {
      stop(
        "conditionalParamsForX_custom: invalid numerator sum = ",
        format(total, digits = 16, scientific = TRUE),
        ". logNumer = c(", paste(format(logNumer, digits = 16, scientific = TRUE), collapse = ", "), ").",
        call. = FALSE
      )
    }
    post_w <- numerator / total
  }


  # Conditional z | x, k parameters
  cond_means <- numeric(K)
  cond_vars  <- numeric(K)

  for (k in seq_len(K)) {
    Sigmak       <- getSigma(k)
    Sigma_xx <- stabilize_cov(Sigmak[1:q, 1:q, drop = FALSE])  # Σ_xx
    Sigma_xz <- Sigmak[1:q, q+1, drop = FALSE]     # Σ_xz (q×1)
    Sigma_zx <- t(Sigma_xz)                        # Σ_zx (1×q)
    sigma_zz <- Sigmak[q+1, q+1, drop = FALSE]     # σ_zz

    mu_x <- muMat[1:q, k]
    mu_z <- muMat[q+1, k]

    # Solve y = Σ_xx^{-1}
    R <- chol(Sigma_xx + diag(1e-10, q))             # ridge = 1e-10, for stability
    solve_Sxx <- function(v) backsolve(R, backsolve(R, v, transpose = TRUE))

    delta_x      <- xVec - mu_x                      # (x - μ_x) (q)
    Sxx_inv_dx   <- solve_Sxx(delta_x)               # Σ_xx^{-1}(x - μ_x) (q) - for mu
    Sxx_inv_xz   <- solve_Sxx(Sigma_xz)              # Σ_xx^{-1}Σ_xz (q×1) - for var

    cond_means[k] <- as.numeric(mu_z + Sigma_zx %*% Sxx_inv_dx)
    cond_vars[k] <- as.numeric(sigma_zz - Sigma_zx %*% Sxx_inv_xz)
  }

  list(
    post_weights = post_w,       # p(k | x)
    cond_means   = cond_means,   # E[z | x, k]
    cond_vars    = cond_vars     # Var[z | x, k]
  )
}

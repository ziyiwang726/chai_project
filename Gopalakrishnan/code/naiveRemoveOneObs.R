naiveRemoveOneObs <- function(i, df, fit, zMat, eps = 1e-8) {
  N <- nrow(df)  # total number of features
  K <- fit$G

  # identify X columns (except 'z')
  xcols <- setdiff(names(df), "z")

  # Drop the i-th row
  xAll <- as.matrix(df[-i, xcols, drop = FALSE])
  zAll <- df$z[-i]
  q <- ncol(xAll)

  # reweight weight matrix
  W <- zMat[-i, , drop=FALSE]             # Get the posterior weight without i-th
  nMinusOne <- colSums(W)
  piMinusOne <- nMinusOne / (N - 1)

  muMinusOne <- matrix(NA_real_, nrow = ncol(df), ncol = K)  # each column is the (x-mean, z-mean) for cluster k
  SigmaList  <- vector("list", length=K)   # cluster-specific cov for each k, before merging

  for (k in 1:K) {
    w  <- W[, k, drop = FALSE]                    # (N-1) x 1
    sw <- nMinusOne[k]

    if (sw < eps) {
      # In case there is no weight left for this component after removing i
      mu_xk <- colMeans(xAll)
      mu_zk <- mean(zAll)
      Sigma_k <- diag(q + 1)
    } else {
      # Weighted means
      mu_xk <- as.numeric(crossprod(w, xAll) / sw)   # length q
      mu_zk <- as.numeric(crossprod(w, zAll) / sw)   # scalar

      # Centered matrices
      xCenter <- sweep(xAll, 2, mu_xk, FUN = "-")          # (N-1) x q
      zCenter <- zAll - mu_zk                              # (N-1)
      D  <- cbind(xCenter, zCenter)                        # (N-1) x (q+1)

      # sqrt of the posterior for weighting
      w_sqrt <- as.numeric(sqrt(w))                    # (N-1) x 1
      Dw <- D * w_sqrt                                 # == (diag(w_sqrt) %*% D)   multiply each row by sqrt(w)
      Sigma_k <- crossprod(Dw) / sw                    # (q+1) x (q+1)
    }

    muMinusOne[, k] <- c(mu_xk, mu_zk)
    SigmaList[[k]]  <- Sigma_k
  }

  list(
    pi   = piMinusOne,      # length K
    mu   = muMinusOne,      # (q+1) x K; rows = x1..xq, z
    Sigma = SigmaList       # list of (q+1) x (q+1) covariance matrices
  )
}

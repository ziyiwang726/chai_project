#' Conditional Hypothesis testing using Auxiliary Information (chai) main function
#'
#' @param z A numeric vector that saving the z-statistics.
#' @param X A numeric matrix/vector (covariates/side-information) that corresponding to z.
#' @param K_vec An integer value or a range specifying the numbers of mixture components.
#' Suggest use 2 or above.
#' If there is a vector/range, the model will automatically select the "best" based on BIC.
#' The default is K_vec = 2:10.
#' @param B An integer value indicating the total number of samples the model generate to estimate the \eqn{\pi_0(x)}.
#' The default is B = 100. More details please see the original paper.
#'
#' @return A list with the following components:
#' \item{z}{The input z-statistics.}
#' \item{X}{The input X (side information).}
#' \item{K}{The optimal number of mixture components.}
#' \item{B}{The input B value.}
#' \item{clFDR}{The conditional local FDR  of every hypotheses.}
#' \item{pi0}{The estimated \eqn{\pi_0(x)} of every hypotheses.}
#' \item{post_w}{The posterior weight \eqn{w_{ik}} of every hypotheses (\eqn{i}) belong to which component (\eqn{k}).}
#' \item{ord}{The order of Hypothesis indices by increasing conditional local FDR (from smallest to largest).}
#' \item{clFDR_sorted}{The sorted conditional local FDR (from smallest to largest).}
#'
#' @seealso
#' \code{\link[clfdrselect]{clfdrselect}} for the rejected hypotheses identified,
#' \code{\link[performance]{performance}} for checking the model performance with known ground truth,
#' \code{\link[direction]{direction}} for the direction of the rejected hypotheses (binary outcomes only).
#'
#' @import mclust
#' @import admix
#'
#' @examples
#' \dontrun{
#' # Generate a simulation
#' set.seed(123)
#' n = 1000; n0 <- 950; n1 <- 50
#' z0 <- rnorm(n0, mean = 0, sd = 1)
#' x0 <- rnorm(n0, mean = 3, sd = 1)
#' z1 <- rnorm(n1, mean = 3, sd = 1)
#' x1 <- rnorm(n1, mean = 6, sd = 1)
#' z <- c(z0, z1)
#' X <- c(x0, x1)
#' gt <- seq((n0+1), n)
#'
#' # Fit the model
#' res <- chai(z, X, K_vec = 2:10, B = 100)
#'
#' # Check the rejection indices
#' clfdrselect(res$clFDR, q = 0.05)
#'
#' # Check performance with ground truth
#' performance(gt, clfdrselect(res$clFDR, q = 0.05))
#' }

#' @export

get_chai_seed <- function(default = 123L) {
  option_seed <- getOption("chai.seed", default)
  seed_value <- suppressWarnings(as.integer(option_seed))
  if (is.na(seed_value)) {
    default
  } else {
    seed_value
  }
}

fit_mclust_with_fallback <- function(df, K_vec = 2:10, timeout_sec = 30, jitter_sd = 1e-6) {
  fit <- NULL
  fit_error <- NULL

  if (.Platform$OS.type == "unix" && is.finite(timeout_sec) && timeout_sec > 0) {
    job <- parallel::mcparallel(mclust::Mclust(df, G = K_vec), silent = TRUE)
    collected <- parallel::mccollect(job, wait = FALSE, timeout = timeout_sec)

    if (length(collected)) {
      fit_candidate <- collected[[1]]
      if (inherits(fit_candidate, "try-error") || inherits(fit_candidate, "error")) {
        fit_error <- simpleError(as.character(fit_candidate))
      } else {
        fit <- fit_candidate
      }
    } else {
      tools::pskill(job$pid, tools::SIGKILL)
      parallel::mccollect(job, wait = FALSE)
      fit_error <- simpleError(sprintf("Mclust exceeded %s seconds.", timeout_sec))
    }
  } else {
    fit_candidate <- tryCatch(mclust::Mclust(df, G = K_vec), error = function(e) e)
    if (inherits(fit_candidate, "error")) {
      fit_error <- fit_candidate
    } else {
      fit <- fit_candidate
    }
  }

  if (!is.null(fit)) {
    return(fit)
  }

  warning(
    "Primary Mclust fit failed or timed out; retrying with scaled+jittered input. Error: ",
    conditionMessage(fit_error)
  )

  df_stable <- as.data.frame(lapply(df, function(col) {
    scaled <- as.numeric(scale(col))
    scaled[!is.finite(scaled)] <- 0
    scaled
  }))

  set.seed(get_chai_seed())
  for (col_name in names(df_stable)) {
    df_stable[[col_name]] <- df_stable[[col_name]] + stats::rnorm(nrow(df_stable), sd = jitter_sd)
  }

  mclust::Mclust(df_stable, G = K_vec)
}

chai <- function(z, X, K_vec = 2:10, B = 100) {
  # require(mclust); require(locfdr); require(admix); require(mvtnorm)

  df <- data.frame(as.data.frame(X))
  names(df) <- paste0("x", seq_len(ncol(df)))   # To avoid there is a 'z' column name inside the X
  df$z <- z
  xcols <- setdiff(names(df), "z")

  fit <- fit_mclust_with_fallback(df, K_vec = K_vec)   # modelNames = "VVV"
  zMat <- fit$z
  n <- nrow(df)

  lfdr_naive_all <- numeric(n)
  minRatioAll <- numeric(n)
  ratioZ_givenX_all <- numeric(n)
  post_w <- matrix(nrow = n, ncol = fit$G)

  for (i in 1:n) {
    np <- naiveRemoveOneObs(i, df, fit, zMat)
    xVec <- as.numeric(df[i, xcols])
    zVal <- df$z[i]
    cp <- conditionalParamsForX_custom(xVec, np$pi, np$mu, np$Sigma)
    post_w[i,] <- cp$post_weights

    set.seed(get_chai_seed())
    rMix1 <- rGaussianMix(n = B, cp$post_weights, cp$cond_means, sqrt(cp$cond_vars))

    admixMod <- admix::admix_model(knownComp_dist = "norm",
                                   knownComp_param = c("mean" = 0, "sd" = 1))

    result <- admix::admix_estim(samples = list(rMix1),
                                 admixMod = list(admixMod),
                                 est_method = 'PS')

    minRatioAll[i] <- 1 - get_mixing_weights(result)  # pi0

    ratioZ_givenX_all[i] <- ratio_z_given_x_custom(zVal, cp$post_weights, cp$cond_means, cp$cond_vars)

    lfdr_naive_all[i] <- minRatioAll[i] / ratioZ_givenX_all[i]
  }

  ord <- order(lfdr_naive_all)
  lFDR_sorted <- lfdr_naive_all[ord]
  avgFDR <- cumsum(lFDR_sorted) / seq_along(lFDR_sorted)

  return(list(
    z = z,
    X = X,
    K = fit$G,
    B = B,
    clFDR = lfdr_naive_all,
    pi0 = minRatioAll,
    post_w = post_w,
    ord = ord,
    clFDR_sorted = lFDR_sorted,
    avgFDR = avgFDR
  ))
}

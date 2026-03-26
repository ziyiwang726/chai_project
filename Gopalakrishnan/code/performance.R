#' Check model performance with known ground truth
#'
#' @param gt An index vector that saving ground truth indices.
#' @param sel A index vector that saving the indices of the selected hypotheses.
#' @param z A numeric vector that saving the z-statistics.
#' This only needs to be provided when Benjamini-Hochberg (BH) performance is also being evaluated.
#' The default is set to NULL.
#' @param bh_level The target FDR level for Benjamini-Hochberg (BH) performance. The default is bh_level = 0.05.
#' @param one.sided A Boolean operators indicating whether the p-value used to evaluate BH performance
#' are computed from the z-statistics using a one-sided or two-sided test.
#'
#' @return A list with the following components:
#' \item{FDP}{The false discovery proportion calculated based on the ground truth and selections.}
#' \item{Power}{The power calculated based on the ground truth and selections.}
#' \item{True_Signal}{The number or true signals from the ground truth.}
#' \item{Selected}{The number of selections.}
#' \item{True_Positive}{The number of selections that are true signals.}
#' \item{False_Positive}{The number of selections that are NOT true signals, only in selected.}
#' \item{False_Negative}{The number of true signals that were NOT selected, only in ground truth.}
#' \item{BH_FDP}{The false discovery proportion calculated based on the ground truth and BH selections.}
#' \item{BH_Power}{The power calculated based on the ground truth and BH selections.}
#' \item{BH_True_Positive}{The number of BH selections that are true signals.}
#' \item{BH_False_Positive}{The number of BH selections that are NOT true signals.}
#' \item{BH_False_Negative}{The number of true signals that were NOT selected by BH.}
#' \item{BH_Selected}{The number of BH selections.}
#'
#' @seealso
#' \code{\link[chai]{chai}} for the Conditional Hypothesis testing using Auxiliary Information (chai) function,
#' \code{\link[clfdrselect]{clfdrselect}} for the rejected hypotheses identified,
#' \code{\link[direction]{direction}} for the direction of the rejected hypotheses (binary outcomes only).
#'
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
#' res <- chai(z, X, K_vec = 2:6, R = 100)
#'
#' # Check the rejection indices
#' clfdrselect(res$clFDR, q = 0.05)
#'
#' # Check performance with ground truth
#' performance(gt, clfdrselect(res$clFDR, q = 0.05))
#' }

#' @export
performance <- function(gt, sel, z = NULL,
                        bh_level = 0.05, one.sided = TRUE) {
  gt <- as.vector(gt)
  sel <- as.vector(sel)

  # True positive
  TP <- length(intersect(gt, sel))
  # False positive
  FP <- length(setdiff(sel, gt))
  # False negative
  FN <- length(setdiff(gt, sel))

  FDP <- ifelse((TP + FP) > 0, FP / (TP + FP), 0)
  Power <- ifelse((TP + FN) > 0, TP / (TP + FN), 0)


  results <- list(
    FDP = FDP,
    Power = Power,
    True_Signal = length(gt),
    Selected = length(sel),
    True_Positive = TP,
    False_Positive = FP,
    False_Negative = FN
  )

  # BH performance, only if z is not NULL
  if (!is.null(z)) {
    if (one.sided) {
      p_values <- pnorm(z, lower.tail = FALSE)
    } else {
      p_values <- 2 * (1 - pnorm(abs(z)))
    }

    p_adjusted <- p.adjust(p_values, method = "BH")
    bh <- which(p_adjusted < bh_level)

    bh_TP <- length(intersect(gt, bh))
    bh_FP <- length(setdiff(bh, gt))
    bh_FN <- length(setdiff(gt, bh))
    bh_FDP <- ifelse((bh_TP + bh_FP) > 0, bh_FP / (bh_TP + bh_FP), 0)
    bh_Power <- ifelse((bh_TP + bh_FN) > 0, bh_TP / (bh_TP + bh_FN), 0)

    results$BH_FDP <- bh_FDP
    results$BH_Power <- bh_Power
    results$BH_True_Positive <- bh_TP
    results$BH_False_Positive <- bh_FP
    results$BH_False_Negative <- bh_FN
    results$BH_Selected <- length(bh)
  }

  return(results)
}

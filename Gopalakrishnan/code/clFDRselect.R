#' Conditional local FDR-based selection.
#' @description
#' This applies to chai output or estimated clFDR values.
#'
#' The procedure is as follows: we first sort the clFDR values in ascending order, and then define the number of rejection \eqn{l} as:
#'
#' \eqn{L = \max \left\{ l \in \{1, \dots, m\} \,:\, \frac{1}{l} \sum_{i=1}^{l} \text{clfdr}_{(i)} \leq q\right\}}
#'
#' where \eqn{q} is the desired FDR level.
#'
#' @param c A chai output or a numeric vector of conditional local FDR.
#' @param q The target FDR level. The default is q = 0.05.
#' @param max_clFDR The maximum allowable threshold when selecting rejections based on the target FDR level.
#' The default is max_clFDR = 1, meaning that no additional filtering is applied to the selected clFDR values.
#' In extreme cases, when there are many very small clFDR values, some relatively large clFDR values may also be rejected due to compensation in the average.
#' In such situations, it may be worth considering lowering max_clFDR.
#' @return An index vector indicating the indices of the rejected hypotheses, ordered by clFDR from smallest to largest.

#' @seealso
#' \code{\link[chai]{chai}} for the Conditional Hypothesis testing using Auxiliary Information (chai) function,
#' \code{\link[performance]{performance}} for checking the model performance with known ground truth.
#' \code{\link[direction]{direction}} for the direction of the rejected hypotheses (binary outcomes only).

#' @export

clfdrselect <- function(c, q = 0.05, max_clFDR = 1) {
  if (is.list(c) && !is.null(c$avgFDR)) {
    ord <- c$ord
    clFDR_sorted <- c$clFDR_sorted
    avgFDR <- c$avgFDR
  } else {
    clfdr <- as.numeric(c)
    ord <- order(clfdr)
    clFDR_sorted <- clfdr[ord]
    avgFDR <- cumsum(clFDR_sorted) / seq_along(clFDR_sorted)
  }
  valid <- which(avgFDR <= q & clFDR_sorted <= max_clFDR)
  if (length(valid) == 0) return(integer(0))
  ord[seq_len(max(valid))]
}

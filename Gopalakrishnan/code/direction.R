#' Inspect the direction of rejected hypotheses (binary outcomes only)
#' @description
#' This function is designed to quickly determine the direction of the hypotheses rejected by the model when the outcome is binary.
#' It is recommended for use together with \code{\link[chai]{chai}} and \code{\link[clfdrselect]{clfdrselect}}.
#' Based on the sign of the z-statistics, the function identifies each rejected feature as either higher or lower in level 1 of the outcome.
#'
#' @param outcome An binary factor outcome of length n.
#' Please ensure that its ordering matches the sample order used when the z-statistics were originally computed.
#' @param z_value A numberic vector of length m containing the z-statistics.
#' @param selected_idx An index vector containing the indices of the hypotheses rejected by the model.
#' Please ensure that the indices are aligned with the ordering of the z-statistics.
#'
#' @return A data frame with the following columns:
#' \item{feature}{The feature names of the selected rejected hypotheses, taken from the names of the `z_value` input.
#' If names are available, they are retained; otherwise, the indices are used.}
#' \item{z_value}{The corresponding z-statistic values.}
#' \item{direction}{The direction of each selected rejected hypothesis with respect to the binary outcome.}
#'
#' @seealso
#' \code{\link[chai]{chai}} for the Conditional Hypothesis testing using Auxiliary Information (chai) function,
#' \code{\link[clfdrselect]{clfdrselect}} for the rejected hypotheses identified,
#' \code{\link[performance]{performance}} for checking the model performance with known ground truth.

#' @export
direction <- function(outcome, z_value, selected_idx) {
  if (!is.factor(outcome) || length(levels(outcome)) != 2) {
    stop("'outcome' must be a binary factor (a factor with exactly 2 levels).")
  }

  first_level <- levels(outcome)[1]
  lab_hi <- paste0("higher_in_", first_level)
  lab_lo <- paste0("lower_in_", first_level)

  if (!is.null(names(z_value))) {
    feature_names <- names(z_value)[selected_idx]
  } else {
    # If no names exist, use the index numbers as strings
    feature_names <- as.character(selected_idx)
  }

  data.frame(
    feature   = feature_names,
    z_value   = unname(z_value[selected_idx]),
    direction = ifelse(z_value[selected_idx] > 0, lab_hi, lab_lo),
    stringsAsFactors = FALSE
  )
}

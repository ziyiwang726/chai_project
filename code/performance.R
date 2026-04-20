performance <- function(ground_truth, selected_features, z = NULL, bh_level = 0.05, one.sided = TRUE) {
  ground_truth <- as.vector(ground_truth)
  selected_features <- as.vector(selected_features)

  # True positive
  TP <- length(intersect(ground_truth, selected_features))
  # False positive
  FP <- length(setdiff(selected_features, ground_truth))
  # False negative
  FN <- length(setdiff(ground_truth, selected_features))

  FDP <- ifelse((TP + FP) > 0, FP / (TP + FP), 0)
  Power <- ifelse((TP + FN) > 0, TP / (TP + FN), 0)


  results <- list(
    FDP = FDP,
    Power = Power,
    'True Signal' = length(ground_truth),
    'Selected Signal' = length(selected_features),
    'True Positive' = TP,
    'False Positive (Only in selected)' = FP,
    'False Negative (Only in ground truth)' = FN
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

    bh_TP <- length(intersect(ground_truth, bh))
    bh_FP <- length(setdiff(bh, ground_truth))
    bh_FN <- length(setdiff(ground_truth, bh))
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

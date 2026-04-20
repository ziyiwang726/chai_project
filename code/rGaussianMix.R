# Function to simulate a mixture of univariate Gaussians
# Args:
#   n       : total number of samples
#   pi      : vector of mixing proportions (length K, summing to 1)
#   mu      : vector of means (length K)
#   sigma   : vector of standard deviations (length K)
# Returns:
#   A numeric vector of length n
rGaussianMix <- function(n, pi, mu, sigma, debug_label = NULL) {
  fmt <- function(x) paste(format(x, digits = 16, scientific = TRUE), collapse = ", ")
  label <- if (!is.null(debug_label) && nzchar(debug_label)) {
    paste0(" [", debug_label, "]")
  } else {
    ""
  }

  # Check inputs
  K <- length(pi)
  if (!is.numeric(pi) || anyNA(pi) || any(!is.finite(pi))) {
    stop("Invalid mixing proportions", label, ": pi has NA/NaN/Inf values. pi = c(",
         fmt(pi), ")", call. = FALSE)
  }

  sum_pi <- sum(pi)
  if (!is.finite(sum_pi) || abs(sum_pi - 1) > 1e-8) {
    stop("Mixing proportions must sum to 1", label, ". sum(pi) = ",
         format(sum_pi, digits = 16, scientific = TRUE),
         ", pi = c(", fmt(pi), ")", call. = FALSE)
  }

  if (isTRUE(getOption("chai.debug_pi", FALSE))) {
    message("rGaussianMix", label, " pi = c(", fmt(pi), "); sum(pi) = ",
            format(sum_pi, digits = 16, scientific = TRUE))
  }

  if (!(length(mu) == K && length(sigma) == K)) {
    stop("mu and sigma must have the same length as pi.")
  }

  # 1. Sample component labels for each observation
  comp <- sample(seq_len(K), size = n, replace = TRUE, prob = pi)

  # 2. Generate from the corresponding normal for each label
  x <- rnorm(n, mean = mu[comp], sd = sigma[comp])

  # 3. Return the simulated data
  return(x)
}

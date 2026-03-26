# Function to simulate a mixture of univariate Gaussians
# Args:
#   n       : total number of samples
#   pi      : vector of mixing proportions (length K, summing to 1)
#   mu      : vector of means (length K)
#   sigma   : vector of standard deviations (length K)
# Returns:
#   A numeric vector of length n
rGaussianMix <- function(n, pi, mu, sigma) {
  # Check inputs
  K <- length(pi)
  if (abs(sum(pi) - 1) > 1e-8) stop("Mixing proportions must sum to 1.")
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


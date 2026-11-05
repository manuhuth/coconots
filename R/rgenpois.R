#' @importFrom stats runif
#' @importFrom HMMpa dgenpois
rgenpois <- function(n, lambda, eta){
  u <- stats::runif(n)
  if (length(u) == 0) {
    return(numeric(0))
  }

  # Inversion of the Generalized Poisson CDF in C++ (rgenpoisInv), one
  # uniform per draw. lambda may be a scalar or a vector of per-draw means
  # (recycled to length(u)). The pmf reproduces HMMpa::dgenpois and the CDF is
  # accumulated from 0 upwards with the strict `last < u` stopping rule of the
  # former R implementation, so the drawn values are unchanged; for a constant
  # lambda the CDF table is built once and shared by all draws. If the
  # accumulated CDF saturates below u (rounding), the inversion returns the
  # count at saturation; a draw above 170 stops with the count-limit error
  # (see GPInvCdf in cfunctions.cpp).
  rgenpoisInv(u, rep_len(as.numeric(lambda), length(u)), as.numeric(eta))
}

#' @importFrom stats runif
#' @importFrom HMMpa dgenpois
rgenpois <- function(n, lambda, eta){
  u <- stats::runif(n)
  if (length(u) == 0) {
    return(numeric(0))
  }

  # lambda and eta are constant across draws, so the CDF table is built once
  # (up to the largest uniform) and all draws are inverted against it. The
  # accumulation order matches the former per-draw while-loop, so the
  # inversion indices are identical.
  maxu <- max(u)
  cdf <- HMMpa::dgenpois(0, lambda1 = lambda, lambda2 = eta)
  last <- cdf
  index <- 0
  while (last < maxu) {
    index <- index + 1
    last <- last + HMMpa::dgenpois(index, lambda1 = lambda, lambda2 = eta)
    cdf[index + 1] <- last
  }

  # smallest index j with cdf[j + 1] >= u; left.open reproduces the strict
  # `while (nain < unif)` stopping rule (ties stop the loop)
  as.numeric(findInterval(u, cdf, left.open = TRUE))
}

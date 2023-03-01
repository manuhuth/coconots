#' @title Compute the cumulative distribution function (CDF) of the GP1 convolution
#' @description Computes the CDF of the the GP2 convolution distribution with parameters lambda, alpha and eta, at a given quantile q, for a given value of y . 
#' @param q Value at which to evaluate the CDF
#' @param y non-negative integer 
#' @param par A vector of 3 parameters: lambda, alpha and eta
#' @return The CDF of the the GP1 convolution distribution at the given quantile
#' @export
pGP1 <- function(q, y, par) { # q is number we want to know the quantile of
  if ((round(y) != y) | (y < 0)) {
    stop("y must be a non-negative integer")
  }

  if (length(par) != 3) {
    stop("three parametes, lambda, alpha and eta, must be passed to the function")
  }

  q <- floor(q)

  lambda <- par[1]
  alpha <- par[2]
  eta <- par[3]

  if (lambda < 0) {
    stop("lambda must be greater than zero")
  }

  if ((alpha < 0) | (alpha > 1) | (eta < 0) | (eta > 1)) {
    stop("alpha and eta must be between or equal to zero and one")
  }

  CDF <- 0
  if (q >= 0) {
    for (j in 0:q) {
      CDF <- CDF + dGP1(j, y, par)
    }
  }


  return(CDF)
}

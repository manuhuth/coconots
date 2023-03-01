#' @title Compute the cumulative distribution function (CDF) of the random operator used in the GP2 model
#' @description Computes the CDF of the GP2 convolution distribution with parameters lambda,alpha1, alpha2, alpha3 and eta, at a given value q, for a given value of y and z . 
#' @param q The value at which to evaluate the CDF
#' @param y non-negative integer 
#' @param z non-negative integer 
#' @param par A vector of 5 parameters: lambda, alpha1, alpha2, alpha3 and eta
#' @return The CDF of the the GP2 convolution distribution at the given quantile
#' @export

pGP2 <- function(q, y, z, par) { # q is number we want to know the quantile of

  if ((round(y) != y) | (y < 0)) {
    stop("y must be a non-negative integer")
  }

  if ((round(z) != z) | (z < 0)) {
    stop("z must be a non-negative integer")
  }

  q <- floor(q)

  if (length(par) != 5) {
    stop("five parameters, lambda, alpha1, alpha2, alpha3 and eta, must be passed to the function")
  }

  lambda <- par[1]
  alpha1 <- par[2]
  alpha2 <- par[3]
  alpha3 <- par[4]
  eta <- par[5]

  if (lambda < 0) {
    stop("lambda must be greater than zero")
  }

  if ((alpha1 < 0) | (alpha1 > 1) | (eta < 0) | (eta > 1) | (alpha2 < 0) | (alpha2 > 1) | (alpha3 < 0) | (alpha3 > 1)) {
    stop("alpha and eta must be between or equal to zero and one")
  }



  CDF <- 0
  if (q >= 0) {
    for (j in 0:q) {
      CDF <- CDF + dGP2(j, y, z, par)
    }
  }
  return(CDF)
}

dGP2 <- function(x, y, z, par) {
  if ((round(y) != y) | (y < 0)) {
    stop("y must be a non-negative integer")
  }

  if ((round(z) != z) | (z < 0)) {
    stop("z must be a non-negative integer")
  }

  if (round(x) != x) {
    stop("x must be an integer")
  }

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



  if (x < 0) {
    result <- 0
  }

  if (x >= 0) {
    result <- dGP2h(x, y, z, lambda, alpha1, alpha2, alpha3, eta)
  }
  return(result)
}

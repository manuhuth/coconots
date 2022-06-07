rgenpois <- function(n, lambda, eta){
  vec <- numeric(n)

  for (j in 1:n) {
    index <- 0
    unif <- runif(1)
    nain <- dgenpois(0, lambda1 = lambda, lambda2 = eta)
    while(nain < unif) {
      index <- index + 1
      nain <- nain + dgenpois(index, lambda1 = lambda, lambda2 = eta)
    }
    vec[j] <- index
  }
  return(vec)
}

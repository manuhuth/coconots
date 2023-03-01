#' @title Generate Random Numbers from Generalized Poisson Distribution
#' @description Generates a vector of random numbers from a generalized Poisson distribution with parameters lambda and eta. The function uses a while loop to generate the numbers.
#' @param n Number of random numbers to generate
#' @param lambda Parameter for the generalized Poisson distribution
#' @param eta Parameter for the generalized Poisson distribution
#' @return A vector of random numbers from the generalized Poisson distribution
#' @export

rgenpois <- function(n, lambda, eta){
  vec <- numeric(n)

  for (j in 1:n) {
    index <- 0
    unif <- runif(1)
    nain <- dgenpois(0, lambda1 = lambda, lambda2 = eta)
    while(nain < unif) {
      index <- index + 1
      nain <- nain + HMMpa::dgenpois(index, lambda1 = lambda, lambda2 = eta)
    }
    vec[j] <- index
  }
  return(vec)
}

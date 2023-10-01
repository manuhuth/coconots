cocoSim_base <- function(type, order, par, size, seasonality = c(1, 2), init = NULL) {
  
  if (length(seasonality == 1)) {
    seasonality <- c(seasonality, seasonality + 1)
  }
  
  if ((type != "GP") & (type != "Poisson")) {
    stop("Option 'type' must be either Poisson or GP")
  }

  if ((order != 1) & order != (2)) {
    stop("Option 'order' must be 1 or 2")
  }

  if (!is.numeric(par)) {
    stop("Option 'par' must be a numeric vector")
  }

  if (par[1] <= 0) {
    stop("lambda must be a positive reel number")
  }

  if (!((all(par[-1] > 0)) & all(par[-1] < 1))) {
    stop(" The autoregressive parameters and eta must be between zero and one")
  }

  if (length(seasonality) == 2) {
    if (seasonality[1] >= seasonality[2]) {
      stop("The first parameter of 'seasonality' must not be smaller than the second parameter")
    }

    if ((seasonality[2] != round(seasonality[2])) | (seasonality[2] < 1)) {
      stop("The values of 'seasonality' must be positive integer values")
    }
  }

  if ((seasonality[1] != round(seasonality[1])) | (seasonality[1] < 1)) {
    stop("The values of 'seasonality' must be positive integer values")
  }

  if (!is.numeric(size)) {
    stop("The value of 'size' must be a positive integer value")
  }

  if ((size != round(size)) | (size < 1)) {
    stop("The value of 'size' must be a positive integer value")
  }

 



  start_time <- Sys.time()

  T <- size

  # 1 Poisson 1
  if ((order == 1) & (type == "Poisson")) {
    if (length(par) != 2) {
      stop(paste("Number of parameteres must equal 2 for the Poisson 1 model"))
    }

    lambda <- par[1]
    alpha <- par[2]
    eta <- 0
    psi <- eta * (1 - alpha) / lambda

    
    data <- stats::rpois(n = seasonality[1], lambda)
    if (!is.null(init) ) {
      data <- init[(length(init) - length(data)+1):(length(init))]
    } 
    
    N <- length(data)
    data <- c(data, rep(NaN, T - N))
    innovations <- stats::rpois(n = T, lambda)
    uniform <- stats::runif(n = T, 0, 1)

    data <- simGP1(20, lambda, alpha, eta, T, N, seasonality[1], data, uniform, innovations)


    end_time <- Sys.time()
    time <- end_time - start_time
    output <- list("time" = time, "data" = data)
  } # end Poisson 1



  # 2 GP1
  if ((order == 1) & (type == "GP")) {
    if (length(par) != 3) {
      stop(paste("Number of parameteres must equal 3 for the GP 1 model"))
    }
    lambda <- par[1]
    alpha <- par[2]
    eta <- par[3]

    data <- rgenpois(n = seasonality[1], lambda, eta)

    if (!is.null(init) ) {
      data <- init[(length(init) - length(data)+1):(length(init))]
    }

    N <- length(data)
    data <- c(data, rep(NaN, T - N))
    innovations <- rgenpois(n = T, lambda, eta)
    uniform <- stats::runif(n = T, 0, 1)

    data <- simGP1(20, lambda, alpha, eta, T, N, seasonality[1], data, uniform, innovations)

    end_time <- Sys.time()
    time <- end_time - start_time
    output <- list("time" = time, "data" = data)
  } # end GP1


  # 3 Poisson 2
  if ((order == 2) & (type == "Poisson")) {
    if (length(par) != 4) {
      stop(paste("Number of parameteres must equal 4 for the Poisson 2 model"))
    }
    lambda <- par[1]
    alpha1 <- par[2]
    alpha2 <- par[3]
    alpha3 <- par[4]
    eta <- 0

    data <- stats::rpois(n = seasonality[2], lambda)

    if (!is.null(init) ) {
      data <- init[(length(init) - length(data)+1):(length(init))]
    }

    N <- length(data)
    data <- c(data, rep(NaN, T - N))
    innovations <- stats::rpois(n = T, lambda)
    uniform <- stats::runif(n = T, 0, 1)
    data <- simGP2(
      20, lambda, alpha1, alpha2, alpha3, eta, T, N, seasonality[1],
      seasonality[2], data, uniform, innovations
    )

    end_time <- Sys.time()
    time <- end_time - start_time
    output <- list("time" = time, "data" = data)
  } # end P2



  # 4 GP2
  if ((order == 2) & (type == "GP")) {
    if (length(par) != 5) {
      stop(paste("Number of parameteres must equal 5 for the GP 2 model"))
    }
    lambda <- par[1]
    alpha1 <- par[2]
    alpha2 <- par[3]
    alpha3 <- par[4]
    eta <- par[5]

    data <- rgenpois(n = seasonality[2], lambda, eta)

    if (!is.null(init) ) {
      data <- init[(length(init) - length(data)+1):(length(init))]
    }

    N <- length(data)
    data <- c(data, rep(NaN, T - N))
    innovations <- rgenpois(n = T, lambda, eta)
    uniform <- stats::runif(n = T, 0, 1)
    data <- simGP2(
      20, lambda, alpha1, alpha2, alpha3, eta, T, N, seasonality[1],
      seasonality[2], data, uniform, innovations
    )

    end_time <- Sys.time()
    time <- end_time - start_time
    output <- list("time" = time, "data" = data)
  } # end GP2
  
  if (!is.null(init)){
    output$data <- output$data[(order+1):length(output$data)]
  }

  return(output)
}

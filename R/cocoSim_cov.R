cocoSim_cov <- function(type, order, par, size, xreg, seasonality = c(1, 2), init = NULL) {
  
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

  if (any(par[-c((length(par) - ncol(xreg)):length(par))] < 0) |
    any(par[-c((length(par) - ncol(xreg)):length(par))] > 1)) {
    stop("The autoregressive parameters and eta must be between zero and one")
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

  if (!(is.atomic(init)) | (!(is.numeric(init)) & !(is.null(init)))) {
    stop("The option 'init' must either be NULL or a numeric vector")
  }


  start_time <- Sys.time()

  T <- size


  # 1 Poisson 1
  if ((order == 1) & (type == "Poisson")) {
    if (length(par) != 1 + ncol(xreg)) {
      stop(paste("Number of parameters must equal", 1 + ncol(xreg), "for the Poisson 1 model with", ncol(xreg), "covariates"))
    }
    alpha <- par[1]
    eta <- 0


    vec_lambda <- c()
    for (j in 1:ncol(xreg)) {
      nam <- paste("lambda", j, sep = "")
      vec_lambda[j] <- assign(nam, par[j + 1])
    }

    data <- c()
    for (t in 1:seasonality[1]) {
      lambda_start1 <- exp(as.numeric(as.vector(xreg[t, ])) %*% vec_lambda)
      data[t] <- stats::rpois(n = 1, lambda_start1)
    }

    if (!is.null(init)) {
      data <- init
    }

    lambdas <- exp(xreg %*% vec_lambda)

    N <- length(data)
    data <- c(data, rep(NaN, T - N))
    innovations <- c()
    for (index in 1:T) {
      lambda <- lambdas[index]
      innovations[index] <- stats::rpois(n = 1, lambda)
    }
    uniform <- stats::runif(n = T, 0, 1)

    data <- simGP1cov(20, alpha, eta, vec_lambda, T, N, seasonality[1], data, xreg, uniform, innovations)


    end_time <- Sys.time()
    time <- end_time - start_time
    output <- list("time" = time, "data" = data)
  } # end Poisson 1



  # 2 GP1
  if ((order == 1) & (type == "GP")) {
    if (length(par) != 2 + ncol(xreg)) {
      stop(paste("Number of parameters must equal", 2 + ncol(xreg), "for the GP 1 model with", ncol(xreg), "covariates"))
    }
    alpha <- par[1]
    eta <- par[2]


    vec_lambda <- c()
    for (j in 1:ncol(xreg)) {
      nam <- paste("lambda", j, sep = "")
      vec_lambda[j] <- assign(nam, par[j + 2])
    }

    data <- c()
    for (t in 1:seasonality[1]) {
      lambda_start1 <- exp(as.numeric(as.vector(xreg[t, ])) %*% vec_lambda)
      data[t] <- rgenpois(n = 1, lambda_start1, eta)
    }

    if (!is.null(init)) {
      data <- init
    }

    lambdas <- exp(xreg %*% vec_lambda)

    N <- length(data)
    data <- c(data, rep(NaN, T - N))
    innovations <- c()
    for (index in 1:T) {
      lambda <- lambdas[index]
      innovations[index] <- rgenpois(n = 1, lambda, eta)
    }
    uniform <- stats::runif(n = T, 0, 1)

    data <- simGP1cov(20, alpha, eta, vec_lambda, T, N, seasonality[1], data, xreg, uniform, innovations)


    end_time <- Sys.time()
    time <- end_time - start_time
    output <- list("time" = time, "data" = data)
  } # end GP1


  # 3 Poisson 2
  if ((order == 2) & (type == "Poisson")) {
    if (length(par) != 3 + ncol(xreg)) {
      stop(paste("Number of parameters must equal", 3 + ncol(xreg), "for the Poisson 2 model with", ncol(xreg), "covariates"))
    }
    alpha1 <- par[1]
    alpha2 <- par[2]
    alpha3 <- par[3]
    eta <- 0
    U <- (1 - alpha1 - alpha2 - alpha3)^(-1)


    vec_lambda <- c()
    for (j in 1:ncol(xreg)) {
      nam <- paste("lambda", j, sep = "")
      vec_lambda[j] <- assign(nam, par[j + 3])
    }

    data <- c()
    for (t in 1:seasonality[2]) {
      lambda_start1 <- exp( as.numeric(as.vector(xreg[t, ])) %*% vec_lambda)
      data[t] <- stats::rpois(n = 1, lambda_start1)
    }

    if (!is.null(init)) {
      data <- init
    }

    lambdas <- exp(xreg %*% vec_lambda)

    N <- length(data)
    data <- c(data, rep(NaN, T - N))
    innovations <- c()
    for (index in 1:T) {
      lambda <- lambdas[index]
      innovations[index] <- stats::rpois(n = 1, lambda)
    }
    uniform <- stats::runif(n = T, 0, 1)

    data <- simGP2cov(20, alpha1, alpha2, alpha3, eta, vec_lambda, T, N, seasonality[1], seasonality[2], data, xreg, uniform, innovations)


    end_time <- Sys.time()
    time <- end_time - start_time
    output <- list("time" = time, "data" = data)
  } # end P2



  # 4 GP 2
  if ((order == 2) & (type == "GP")) {
    if (length(par) != 4 + ncol(xreg)) {
      stop(paste("Number of parameters must equal", 4 + ncol(xreg), "for the GP 2 model with", ncol(xreg), "covariates"))
    }

    alpha1 <- par[1]
    alpha2 <- par[2]
    alpha3 <- par[3]
    eta <- par[4]


    vec_lambda <- c()
    for (j in 1:ncol(xreg)) {
      nam <- paste("lambda", j, sep = "")
      vec_lambda[j] <- assign(nam, par[j + 4])
    }

    data <- c()
    for (t in 1:seasonality[2]) {
      lambda_start1 <- exp( as.numeric(as.vector(xreg[t, ])) %*% vec_lambda)
      data[t] <- rgenpois(n = 1, lambda_start1, eta)
    }

    if (!is.null(init)) {
      data <- init
    }

    lambdas <- exp( xreg %*% vec_lambda)

    N <- length(data)
    data <- c(data, rep(NaN, T - N))
    innovations <- c()
    for (index in 1:T) {
      lambda <- lambdas[index]
      innovations[index] <- rgenpois(n = 1, lambda, eta)
    }
    uniform <- stats::runif(n = T, 0, 1)

    data <- simGP2cov(20, alpha1, alpha2, alpha3, eta, vec_lambda, T, N, seasonality[1], seasonality[2], data, xreg, uniform, innovations)

    end_time <- Sys.time()
    time <- end_time - start_time
    output <- list("time" = time, "data" = data)
  } # end GP2

  end_time <- Sys.time()
  time <- end_time - start_time
  
  if (is.null(init)){
    warning("No burn-in period is specified using the init argument. Hence, the resulting simulated time series might not be stationary. You can add a custom burn-in period by passing it to the init argument. This could be, for example, done by simulating a burn-in period with appropriate covariates using cocoSim and passing the resulting time series to the init argument of a new cocoSim run.") 
  } else {
    output$data <- output$data[(length(init) + 1):length(output$data)]
  }
  
  return(output)
}

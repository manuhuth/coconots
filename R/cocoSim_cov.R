#' @importFrom stats rpois
#' @importFrom stats runif
cocoSim_cov <- function(type, order, par, size, xreg, seasonality = c(1, 2), init = NULL, link_function="log") {
  
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

  start_time <- Sys.time()

  T <- size


  # 1 Poisson 1
  if ((order == 1) & (type == "Poisson")) {
    if (length(par) != 1 + ncol(xreg)) {
      stop(paste("Number of parameters must equal", 1 + ncol(xreg), "for the Poisson 1 model with", ncol(xreg), "covariates"))
    }
    alpha <- par[1]
    eta <- 0


    vec_lambda <- par[-1]
    lambdas <- simCovLambdas(xreg, vec_lambda, link_function, T, init, order)

    data <- c()
    for (t in 1:seasonality[1]) {
      lambda_start_1 <- as.numeric(as.vector(xreg[t, ])) %*% vec_lambda
      data[t] <- stats::rpois(n = 1, applyLinkFunction(lambda_start_1, link_function))
    }
    
    if (!is.null(init) ) {
      data <- init[(length(init) - length(data)+1):(length(init))]
    } 

    N <- length(data)
    data <- c(data, rep(NaN, T - N))
    innovations <- stats::rpois(n = length(lambdas), lambda = lambdas)
    uniform <- stats::runif(n = T, 0, 1)

    rows <- simCovRows(T, N, init, xreg)
    data <- simGP1cov(alpha, eta, vec_lambda, T, N, seasonality[1], data,
                      xreg[rows, , drop = FALSE], uniform, innovations[rows],
                      link_function)

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


    vec_lambda <- par[-(1:2)]
    lambdas <- simCovLambdas(xreg, vec_lambda, link_function, T, init, order)

    data <- c()
    for (t in 1:seasonality[1]) {
      lambda_start1 <- as.numeric(as.vector(xreg[t, ])) %*% vec_lambda
      data[t] <- rgenpois(n = 1, applyLinkFunction(lambda_start1, link_function),
                          eta)
    }

    if (!is.null(init) ) {
      data <- init[(length(init) - length(data)+1):(length(init))]
    } 

    N <- length(data)
    data <- c(data, rep(NaN, T - N))
    # one uniform per time point, in time order (same random stream as the
    # former per-index rgenpois(n = 1, ...) loop)
    innovations <- rgenpois(n = length(lambdas), lambda = lambdas, eta = eta)
    uniform <- stats::runif(n = T, 0, 1)

    rows <- simCovRows(T, N, init, xreg)
    data <- simGP1cov(alpha, eta, vec_lambda, T, N, seasonality[1], data,
                      xreg[rows, , drop = FALSE], uniform, innovations[rows],
                      link_function)


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


    vec_lambda <- par[-(1:3)]
    lambdas <- simCovLambdas(xreg, vec_lambda, link_function, T, init, order)

    data <- c()
    for (t in 1:seasonality[2]) {
      if (nrow(xreg) > 1){
        lambda_start1 <- as.numeric(as.vector(xreg[t, ])) %*% vec_lambda
      } else{
        lambda_start1 <- c(0,0)
      }
      data[t] <- stats::rpois(n = 1, applyLinkFunction(lambda_start1, link_function))
    }

    if (!is.null(init) ) {
      data <- init[(length(init) - length(data)+1):(length(init))]
    } 

    N <- length(data)
    data <- c(data, rep(NaN, T - N))
    innovations <- stats::rpois(n = length(lambdas), lambda = lambdas)
    if (nrow(xreg) < 2){
      uniform <- stats::runif(n = T+2, 0, 1)
    } else{
      uniform <- stats::runif(n = T, 0, 1)
    }

    rows <- simCovRows(T, N, init, xreg)
    data <- simGP2cov(alpha1, alpha2, alpha3, eta, vec_lambda, T, N,
                      seasonality[1], seasonality[2], data,
                      xreg[rows, , drop = FALSE], uniform, innovations[rows],
                      link_function)


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


    vec_lambda <- par[-(1:4)]
    lambdas <- simCovLambdas(xreg, vec_lambda, link_function, T, init, order)

    data <- c()
    for (t in 1:seasonality[2]) {
      if (nrow(xreg) > 1){
        lambda_start1 <- as.numeric(as.vector(xreg[t, ])) %*% vec_lambda
      } else{
        lambda_start1 <- c(0,0)
      }
      # start values from the GP innovation distribution at lambda_t (as for
      # order 1 and on the Julia backend)
      data[t] <- rgenpois(n = 1, applyLinkFunction(lambda_start1, link_function),
                          eta)
    }

    if (!is.null(init) ) {
      data <- init[(length(init) - length(data)+1):(length(init))]
    } 

    N <- length(data)
    data <- c(data, rep(NaN, T - N))
    # one uniform per time point, in time order (same random stream as the
    # former per-index rgenpois(n = 1, ...) loop)
    innovations <- rgenpois(n = length(lambdas), lambda = lambdas, eta = eta)
    
    if (nrow(xreg) < 3){
      uniform <- stats::runif(n = T+2, 0, 1)
    } else{
      uniform <- stats::runif(n = T, 0, 1)
    }
    
    rows <- simCovRows(T, N, init, xreg)
    data <- simGP2cov(alpha1, alpha2, alpha3, eta, vec_lambda, T, N,
                      seasonality[1], seasonality[2], data,
                      xreg[rows, , drop = FALSE], uniform, innovations[rows],
                      link_function)

    end_time <- Sys.time()
    time <- end_time - start_time
    output <- list("time" = time, "data" = data)
  } # end GP2

  end_time <- Sys.time()
  time <- end_time - start_time
  
  if (is.null(init)){
    warning("No burn-in period is specified using the init argument. Hence, the resulting simulated time series might not be stationary. You can add a custom burn-in period by passing it to the init argument. This could be, for example, done by simulating a burn-in period with appropriate covariates using cocoSim and passing the resulting time series to the init argument of a new cocoSim run.") 
  } else {
    output$data <- output$data[(order + 1):length(output$data)]
  }
  
  output$link_function <- link_function
  
  return(output)
}


# lambda_t = link(x_t' beta) for every row of xreg (innovation means of the
# simulation, one per row). Simulation needs them finite and positive; the
# identity link gives lambda_t <= 0 for x_t' beta <= 0. Also checks that xreg
# has a row for every simulated observation (see simCovRows).
simCovLambdas <- function(xreg, vec_lambda, link_function, T, init, order) {
  n_needed <- if (is.null(init)) T else T - order
  if (nrow(xreg) < n_needed) {
    stop(paste("xreg must have one row per simulated observation, i.e. at least",
               n_needed, "rows"))
  }
  lambdas <- applyLinkFunctionVec(xreg %*% vec_lambda, link_function)
  bad <- which(!is.finite(lambdas) | lambdas <= 0)
  if (length(bad) > 0) {
    stop(sprintf(paste(
      "The innovation mean lambda_t = link(x_t' beta) must be finite and positive",
      "for every t, but lambda_t = %g at t = %d. With the identity link choose",
      "covariates and beta such that x_t' beta > 0."), lambdas[bad[1]], bad[1]))
  }
  lambdas
}

# Rows of xreg (and of the innovations, one per row) that drive the
# simulated positions N + 1, ..., T of `data`: simGP1cov / simGP2cov use row
# t - N of the xreg they are given (and innovation t - N) at position t.
# Observation t is driven by covariate row t, as in the likelihood:
#  - without init the positions are the observations 1..T (positions 1..N are
#    the start values, drawn at rows 1..N), so position t needs row t;
#  - with init, positions 1..N hold the last N values of init and position t
#    is observation t - N, which needs row t - N.
simCovRows <- function(T, N, init, xreg) {
  if (is.null(init)) {
    N + seq_len(max(T - N, 0))
  } else {
    seq_len(nrow(xreg))
  }
}

cocoForecastKStepsRCPP <- function(fit, k=1, number_simulations=500, covariates=NULL){

  init <- fit$ts[length(fit$ts)]
  if (fit$order == 2){
    init <- fit$ts[(length(fit$ts)-1):length(fit$ts)]
  }

  # All paths in one C++ call (forecastPathsBase / forecastPathsCov): each
  # path is the series the per-path cocoSim_base / cocoSim_cov(init = init)
  # call below returned for the same RNG state (same model, draw order and
  # covariate alignment), so the procedure is unchanged. The R loop remains
  # for inputs outside the regular case.
  out_matrix <- forecastPathsRCPP(fit, k, number_simulations, covariates, init)

  if (is.null(out_matrix)) {
    run_sim <- function(i){
      #cannot be cocosim as the init argument is not available in cocoSim
      if (!is.null(covariates)){
        output <- cocoSim_cov(fit$type, order = fit$order, par = fit$par,
                            size= k+length(init), xreg = covariates,
                            seasonality = c(1, 2), init = init, link_function=fit$link_function)$data
      } else {
        output <- cocoSim_base(fit$type, order = fit$order, par = fit$par,
                              size= k+length(init),
                              seasonality = c(1, 2), init = init)$data
      }
    }

    if (k > 1) {
      out_matrix <- t(sapply(1:number_simulations, run_sim))
    } else if (k==1){
      out_matrix <- matrix(sapply(1:number_simulations, run_sim), ncol=1)
    }
  }

  freq_table <- function(i){
    df <- as.data.frame(table(out_matrix[,i]) / number_simulations)
    colnames <- c("value", "frequency")
    colnames(df) <- colnames
    df["value"] <- as.numeric(levels(df[, "value"]))

    return(df)
  }

  return(lapply(1:ncol(out_matrix), freq_table))

}


# number_simulations x k matrix of simulated paths (column j = step j) from
# one C++ call, or NULL if the inputs are outside the regular case (then the
# caller uses the per-path cocoSim loop). Runs the parameter checks of
# cocoSim_base / cocoSim_cov once, with the same messages.
forecastPathsRCPP <- function(fit, k, number_simulations, covariates, init) {
  type <- fit$type
  order <- fit$order
  par <- fit$par

  regular <- is.numeric(k) && length(k) == 1 && is.finite(k) && k >= 1 && k == round(k) &&
    is.numeric(number_simulations) && length(number_simulations) == 1 &&
    is.finite(number_simulations) && number_simulations >= 1 &&
    number_simulations == round(number_simulations) &&
    is.character(type) && length(type) == 1 && type %in% c("GP", "Poisson") &&
    is.numeric(order) && length(order) == 1 && order %in% c(1, 2) &&
    is.numeric(par) && all(is.finite(par)) &&
    is.numeric(init) && length(init) == order && all(is.finite(init)) &&
    all(init >= 0) && all(init == round(init)) && all(init <= .Machine$integer.max)
  if (!regular) return(NULL)

  gp <- type == "GP"
  init <- as.integer(init)

  if (is.null(covariates)) {
    # checks of cocoSim_base
    if (par[1] <= 0) {
      stop("lambda must be a positive reel number")
    }
    if (!((all(par[-1] > 0)) & all(par[-1] < 1))) {
      stop(" The autoregressive parameters and eta must be between zero and one")
    }
    n_par <- (if (gp) 3 else 2) + 2 * (order - 1)
    if (length(par) != n_par) {
      stop(paste("Number of parameteres must equal", n_par, "for the", if (gp) "GP" else "Poisson",
                 order, "model"))
    }
    if (order == 1) {
      return(forecastPathsBase(1L, gp, par[1], par[2], 0, 0, if (gp) par[3] else 0,
                               init, as.integer(k), as.integer(number_simulations)))
    }
    return(forecastPathsBase(2L, gp, par[1], par[2], par[3], par[4], if (gp) par[5] else 0,
                             init, as.integer(k), as.integer(number_simulations)))
  }

  xreg <- covariates
  link_function <- fit$link_function
  if (!is.matrix(xreg) || !is.numeric(xreg) || nrow(xreg) < k || ncol(xreg) < 1 ||
      any(!is.finite(xreg)) || (order == 2 && nrow(xreg) <= 1) ||
      !is.character(link_function) || length(link_function) != 1) {
    return(NULL)
  }
  storage.mode(xreg) <- "double"

  # checks of cocoSim_cov
  if (any(par[-c((length(par) - ncol(xreg)):length(par))] < 0) |
    any(par[-c((length(par) - ncol(xreg)):length(par))] > 1)) {
    stop("The autoregressive parameters and eta must be between zero and one")
  }
  n_ar <- if (order == 1) 1 else 3
  n_par <- n_ar + gp + ncol(xreg)
  if (length(par) != n_par) {
    stop(paste("Number of parameters must equal", n_par, "for the", if (gp) "GP" else "Poisson",
               order, "model with", ncol(xreg), "covariates"))
  }

  alpha <- par[1:n_ar]
  eta <- if (gp) par[n_ar + 1] else 0
  vec_lambda <- par[-(1:(n_ar + gp))]
  init_lambdas <- numeric(order)
  for (t in 1:order) {
    lambda_start1 <- as.numeric(as.vector(xreg[t, ])) %*% vec_lambda
    init_lambdas[t] <- applyLinkFunction(lambda_start1, link_function)
  }
  lambdas <- applyLinkFunctionVec(xreg %*% vec_lambda, link_function)
  T <- k + order
  n_unif <- T
  if (order == 2 && nrow(xreg) < (if (gp) 3 else 2)) n_unif <- T + 2
  if (order == 1) alpha <- c(alpha, 0, 0)
  forecastPathsCov(as.integer(order), gp, alpha[1], alpha[2], alpha[3], eta,
                   vec_lambda, xreg, lambdas, init_lambdas, init, as.integer(k),
                   as.integer(n_unif), as.integer(number_simulations), link_function)
}

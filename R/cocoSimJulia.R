#' @importFrom JuliaConnectoR juliaCall
cocoSimJulia <- function(type, order, par, length, xreg, link_function){
  addJuliaFunctions()
  
  if (is.null(xreg)){
    par_use <- c(par[2:length(par)], par[1])
  } else {
    par_use <- par
  }
  
  # With covariates no burn-in, as on the RCPP backend (see ?cocoSim).
  n_burn_in <- if (is.null(xreg)) 200L else 0L
  return(JuliaConnectoR::juliaCall("cocoSim", type, order, par_use, length, xreg, link_function, n_burn_in))
}
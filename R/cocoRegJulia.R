#' @importFrom JuliaConnectoR juliaCall
cocoRegJulia <- function(type, order, data, xreg, start, link_function, lower_bound_covariates){
  addJuliaFunctions()
  # fit and drop the "optimization" entry (the solver object, not used in R) in one
  # round trip; returns a proxy of the result dict
  return(JuliaConnectoR::juliaLet(
    'fit = cocoReg(type, order, data, xreg, start, link, lb); delete!(fit, "optimization"); fit',
    type = type, order = order, data = data, xreg = xreg, start = start,
    link = link_function, lb = lower_bound_covariates))
}

#' Parameter, covariance matrix, SEs and log-likelihood of a Julia fit
#'
#' One round trip: the four entries are packed into one Float64 vector
#' (parameter, column-major covariance matrix, se, log_likelihood) instead of
#' materialising the whole result dict (data, covariates, optimizer object)
#' with juliaGet. Values are unchanged (Julia parameter ordering).
#' @noRd
juliaRegSummary <- function(fit_julia) {
  vals <- JuliaConnectoR::juliaLet(
    'convert(Vector{Float64}, vcat(fit["parameter"], vec(fit["covariance_matrix"]), fit["se"], fit["log_likelihood"]))',
    fit = fit_julia)
  vals <- as.numeric(vals)
  n_par <- as.integer(round(sqrt(length(vals)))) - 1L
  if ((n_par + 1L)^2 != length(vals)) {
    stop("Unexpected length of the Julia fit summary.")
  }
  list(parameter = vals[seq_len(n_par)],
       covariance_matrix = matrix(vals[n_par + seq_len(n_par * n_par)], nrow = n_par, ncol = n_par),
       se = vals[n_par + n_par * n_par + seq_len(n_par)],
       log_likelihood = vals[length(vals)])
}

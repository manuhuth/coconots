utils::globalVariables(c(
  "Time", "Observed", "Fitted", "Residuals", "Lag", "ACF", "x", "p"
))

#' Vectorised counterpart of the C++ applyLinkFunction()
#'
#' Applies the link to a whole vector at once instead of one scalar C++
#' round-trip per element. The branch semantics (including the treatment of
#' non-positive and NaN inputs) mirror applyLinkFunction() in cfunctions.cpp.
#'
#' @noRd
applyLinkFunctionVec <- function(x, link_function) {
  x <- as.numeric(x)
  if (link_function == "log") {
    exp(x)
  } else if (link_function == "identity") {
    x
  } else if (link_function == "relu") {
    out <- x
    out[!(!is.na(x) & x > 0)] <- 1e-10
    out
  } else if (link_function == "softplus") {
    out <- numeric(length(x))
    pos <- !is.na(x) & x > 0
    out[pos] <- x[pos] + log(1 + exp(-x[pos]))
    out[!pos] <- log(1 + exp(x[!pos]))
    out
  } else {
    stop("Unknown link function")
  }
}

#' Parameter checks matching the stops in dGP1()/pGP1() and dGP2()/pGP2(),
#' applied once instead of on every scalar density call.
#' @noRd
validate_gp1_par <- function(lambda, alpha, eta) {
  if (any(lambda < 0)) {
    stop("lambda must be greater than zero")
  }
  if ((alpha < 0) | (alpha > 1) | (eta < 0) | (eta > 1)) {
    stop("alpha and eta must be between or equal to zero and one")
  }
  invisible(NULL)
}

#' @noRd
validate_gp2_par <- function(lambda, alpha1, alpha2, alpha3, eta) {
  if (any(lambda < 0)) {
    stop("lambda must be greater than zero")
  }
  if ((alpha1 < 0) | (alpha1 > 1) | (eta < 0) | (eta > 1) | (alpha2 < 0) | (alpha2 > 1) | (alpha3 < 0) | (alpha3 > 1)) {
    stop("alpha and eta must be between or equal to zero and one")
  }
  invisible(NULL)
}

#' Extract a value from a materialised Julia Dict by key name
#'
#' `JuliaConnectoR::juliaGet()` materialises a Julia `Dict` into an R list with
#' parallel `keys` and `values` elements. The order of these elements follows
#' Julia's hash-determined `Dict` iteration order, which is stable for a given
#' Julia version but is NOT guaranteed across versions or if keys are
#' added/renamed on the Julia side. Indexing `values` by integer position is
#' therefore fragile. This helper looks the value up by key name instead.
#'
#' @param julia_dict A list returned by `JuliaConnectoR::juliaGet()` for a Julia
#'   `Dict`, i.e. containing parallel `keys` and `values` elements.
#' @param key The Julia dictionary key (character scalar) to retrieve.
#' @return The value stored under `key`.
#' @noRd
getJuliaValue <- function(julia_dict, key) {
  keys <- as.character(unlist(julia_dict[["keys"]]))
  idx <- match(key, keys)
  if (is.na(idx)) {
    stop(sprintf("Julia result is missing the expected key '%s'.", key))
  }
  julia_dict[["values"]][[idx]]
}
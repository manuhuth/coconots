utils::globalVariables(c(
  "Time", "Observed", "Fitted", "Residuals", "Lag", "ACF", "x", "p"
))

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
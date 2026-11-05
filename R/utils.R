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
#' Default JULIA_NUM_THREADS before the Julia session starts
#'
#' Coconots.jl parallelizes its likelihood evaluations, but a Julia process
#' starts single-threaded unless JULIA_NUM_THREADS is set when it is spawned.
#' JuliaConnectoR reads the environment variable at start-up, so this must run
#' before the first `juliaEval`/`juliaCall` of the session. Defaults to the
#' number of physical cores minus one (leaving one core free for R and the
#' system); a pre-existing JULIA_NUM_THREADS setting (e.g. from .Renviron or
#' the shell) is never overridden. Has no effect on an already-running Julia
#' session.
#'
#' @return Invisibly, the thread count that applies (or `NULL` if unknown).
#' @noRd
set_julia_threads <- function() {
  existing <- Sys.getenv("JULIA_NUM_THREADS")
  if (!identical(existing, "")) {
    return(invisible(existing))
  }
  cores <- tryCatch(parallel::detectCores(logical = FALSE),
                    error = function(e) NA_integer_)
  if (is.na(cores) || cores < 1L) {
    cores <- tryCatch(parallel::detectCores(logical = TRUE),
                      error = function(e) NA_integer_)
  }
  if (is.na(cores) || cores < 1L) {
    return(invisible(NULL))
  }
  threads <- max(1L, cores - 1L)
  Sys.setenv(JULIA_NUM_THREADS = as.character(threads))
  invisible(threads)
}

#' Default thread count for the RCPP likelihood kernels
#'
#' Number of physical cores minus one (leaving one core free for the system),
#' mirroring the set_julia_threads() convention. Capped at 2 when
#' _R_CHECK_LIMIT_CORES_ is set, so R CMD check on CRAN never sees more than
#' the two cores its policy allows. Falls back to 1 if the core count cannot
#' be determined.
#'
#' @noRd
coco_default_cores <- function() {
  cores <- tryCatch(parallel::detectCores(logical = FALSE),
                    error = function(e) NA_integer_)
  if (is.na(cores) || cores < 1L) {
    cores <- tryCatch(parallel::detectCores(logical = TRUE),
                      error = function(e) NA_integer_)
  }
  if (is.na(cores) || cores < 1L) {
    return(1L)
  }
  cores <- max(1L, cores - 1L)
  limit <- Sys.getenv("_R_CHECK_LIMIT_CORES_")
  if (nzchar(limit) && !identical(tolower(limit), "false")) {
    cores <- min(cores, 2L)
  }
  cores
}

#' Polish a constrOptim() fit with box-constrained L-BFGS-B
#'
#' Nelder-Mead inside constrOptim() stalls when the optimum lies on (or close
#' to) a bound such as eta = 0 or alpha_i = 0: the simplex collapses against
#' the infeasible region and the outer barrier loop stops early. This happens
#' e.g. for a GP model fitted to equidispersed counts, or an order-2 model
#' fitted to order-1 data. L-BFGS-B handles such bounds through its active
#' set, so a short L-BFGS-B run started at the constrOptim() solution reaches
#' the constrained optimum.
#'
#' Rows of `ui` with a single non-zero entry define the box bounds (moved
#' inward by `margin`); all rows of `ui`/`ci` are enforced by returning a
#' large finite value outside the feasible region. The polished point is
#' accepted only if it is feasible and lowers the objective `f` (the negative
#' log-likelihood); otherwise, and on any error, `fit` is returned unchanged.
#' `gr` (called like `f`) is the gradient of `f`; without it L-BFGS-B falls
#' back to finite differences.
#' @noRd
coco_polish_fit <- function(fit, f, ui, ci, gr = NULL, ..., margin = 1e-8) {
  dots <- list(...)
  f_eval <- function(p) {
    v <- tryCatch(do.call(f, c(list(p), dots)), error = function(e) NA_real_)
    if (length(v) != 1) NA_real_ else as.numeric(v)
  }
  par0 <- as.numeric(unlist(fit$par, use.names = FALSE))
  ui <- as.matrix(ui)
  n_par <- length(par0)
  if (ncol(ui) != n_par) {
    return(fit)
  }
  val0 <- f_eval(par0)
  if (!is.finite(val0)) {
    return(fit)
  }

  lower <- rep(-Inf, n_par)
  upper <- rep(Inf, n_par)
  for (j in seq_len(nrow(ui))) {
    nz <- which(ui[j, ] != 0)
    if (length(nz) == 1) {
      bound <- ci[j] / ui[j, nz]
      if (ui[j, nz] > 0) {
        lower[nz] <- max(lower[nz], bound + margin)
      } else {
        upper[nz] <- min(upper[nz], bound - margin)
      }
    }
  }
  if (any(lower > upper)) {
    return(fit)
  }

  feasible <- function(p) all(ui %*% p - ci >= 0)
  penalty <- abs(val0) * 1e3 + 1e6
  obj <- function(p) {
    if (!feasible(p)) {
      return(penalty)
    }
    v <- f_eval(p)
    if (is.finite(v)) v else penalty
  }
  # gradient of obj: zero where obj is the constant penalty (infeasible p;
  # a non-finite objective gives a non-finite gradient)
  obj_gr <- if (is.null(gr)) NULL else function(p) {
    zero <- rep(0, n_par)
    if (!feasible(p)) {
      return(zero)
    }
    g <- tryCatch(as.numeric(do.call(gr, c(list(p), dots))), error = function(e) NULL)
    if (length(g) != n_par || any(!is.finite(g))) zero else g
  }

  start <- pmin(pmax(par0, lower), upper)
  if (!feasible(start)) {
    return(fit)
  }

  res <- tryCatch(
    stats::optim(
      par = start, fn = obj, gr = obj_gr, method = "L-BFGS-B",
      lower = lower, upper = upper,
      control = list(maxit = 1000, factr = 1e5, ndeps = rep(1e-6, n_par))
    ),
    error = function(e) NULL
  )
  if (is.null(res)) {
    return(fit)
  }

  p <- as.numeric(res$par)
  if (any(!is.finite(p)) || !feasible(p)) {
    return(fit)
  }
  v <- f_eval(p)
  if (is.finite(v) && v < val0) {
    fit$par <- p
    fit$value <- v
  }
  fit
}

#' Negative log-likelihood with derivatives by automatic differentiation
#'
#' Evaluates the objective `mlf` of cocoReg (the negative log-likelihood of
#' the C++ likelihoods) together with its gradient (`deriv = 1`) or gradient
#' and Hessian (`deriv = 2`) by forward-mode automatic differentiation with
#' dual numbers (C++ `likelihoodAD()`). `par` is in the `coco$par` layout;
#' `xreg` is the covariate matrix or NULL. Returns `list(value, gradient)`
#' and, for `deriv = 2`, `hessian`.
#' @noRd
coco_ad_derivs <- function(par, data, type, order, seasonality, xreg = NULL,
                           link_function = "log", deriv = 2L) {
  likelihoodAD(as.numeric(par), as.integer(order), identical(type, "GP"), data,
               as.integer(seasonality[1]), as.integer(seasonality[2]), xreg,
               link_function, as.integer(deriv))
}

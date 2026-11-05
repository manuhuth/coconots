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

#' Largest supported count
#'
#' All transition densities are built from factorials of the counts, and
#' 171! overflows in double precision. Counts above this limit are outside
#' the supported range on both backends (the C++ factorial raises the same
#' error when a computation would need one of them).
#' @noRd
COCO_MAX_COUNT <- 170L

#' Count values of a data argument (vector, ts, matrix or data frame)
#' @noRd
coco_count_values <- function(data) {
  if (is.data.frame(data)) {
    data <- as.matrix(data)
  }
  as.numeric(data)
}

#' Error if any count exceeds COCO_MAX_COUNT; `what` names the count in the
#' message ("observed", "initial", "simulated").
#' @noRd
coco_check_count_limit <- function(data, what = "observed") {
  x <- coco_count_values(data)
  x <- x[!is.na(x)]
  if (length(x) > 0 && max(x) > COCO_MAX_COUNT) {
    stop(sprintf(paste0("coconots: counts above %d are not supported (factorial overflow ",
                        "in the C++ backend); %s count = %s"),
                 COCO_MAX_COUNT, what, format(max(x))), call. = FALSE)
  }
  invisible(NULL)
}

#' Error if the largest support point of a predictive pmf (argument `arg`,
#' e.g. max_x of cocoScore or max of predict) exceeds COCO_MAX_COUNT.
#' @noRd
coco_check_support_limit <- function(value, arg) {
  if (!is.null(value) && length(value) == 1 && !is.na(value) && value > COCO_MAX_COUNT) {
    stop(sprintf(paste0("coconots: '%s' = %s exceeds the largest supported count %d ",
                        "(counts above %d are not supported: factorial overflow in the ",
                        "C++ backend)"),
                 arg, format(value), COCO_MAX_COUNT, COCO_MAX_COUNT), call. = FALSE)
  }
  invisible(NULL)
}

#' Error if a series of length `n` is shorter than `min_n`; `fun` names the
#' function and `reason` states what the minimum is needed for.
#' @noRd
coco_check_length <- function(n, min_n, fun, reason) {
  if (n < min_n) {
    stop(sprintf("%s: the time series has %d observation(s), but at least %d are required (%s)",
                 fun, as.integer(n), as.integer(min_n), reason), call. = FALSE)
  }
  invisible(NULL)
}

#' Error if a vector of probabilities contains non-finite values (NaN or
#' Inf), which arise from floating-point overflow or underflow in the
#' density recursions (e.g. counts close to the supported maximum with a
#' large dispersion parameter eta).
#' @noRd
coco_check_finite_pmf <- function(p, what) {
  if (any(!is.finite(p))) {
    stop(sprintf(paste0("coconots: the %s is not finite (floating-point overflow in the ",
                        "density recursions for these counts and parameters; this can happen ",
                        "for counts close to the supported maximum of %d with a large eta)"),
                 what, COCO_MAX_COUNT), call. = FALSE)
  }
  invisible(NULL)
}

#' Error if the covariate matrix is rank deficient: the likelihood is then
#' not identified (the regression coefficients of aliased columns are not
#' determined) and the GLM starting values contain NA coefficients. The
#' aliased columns are those that the pivoted QR decomposition (tolerance of
#' stats::lm) moves behind the rank.
#' @noRd
coco_check_xreg_rank <- function(xreg) {
  X <- xreg
  if (is.data.frame(X)) {
    X <- as.matrix(X)
  }
  X <- as.matrix(X)
  if (!is.numeric(X) || ncol(X) == 0 || any(!is.finite(X))) {
    return(invisible(NULL))
  }
  qrx <- qr(X, tol = 1e-07)
  if (qrx$rank < ncol(X)) {
    names_x <- colnames(X)
    if (is.null(names_x)) {
      names_x <- paste0("X", seq_len(ncol(X)))
    }
    names_x[names_x == ""] <- paste0("X", which(names_x == ""))
    aliased <- qrx$pivot[(qrx$rank + 1):ncol(X)]
    stop(sprintf(paste0("xreg is rank deficient: column(s) %s are linear combinations of ",
                        "others (rank %d < %d columns); the model is not identified. ",
                        "Remove the aliased column(s)."),
                 paste(names_x[sort(aliased)], collapse = ", "), qrx$rank, ncol(X)),
         call. = FALSE)
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

#' Default (Nelder-Mead) constrained fit of cocoReg_cov
#'
#' Nelder-Mead inside constrOptim() with the strict outer loop of cocoReg_cov
#' (`outer.it` = 500, `outer.eps` = 1e-10) spends most of its likelihood
#' evaluations on slow final convergence. Here the Nelder-Mead run uses the
#' constrOptim() default outer loop (at most 100 iterations, outer.eps 1e-5),
#' which locates the optimum, and constrOptim() with BFGS and the exact
#' gradient `gr` (C++ AD) then converges from that point with the requested
#' `outer.it`/`outer.eps`; its log-barrier handles all linear constraints
#' `ui %*% par - ci >= 0`. The BFGS result is kept only if its objective is
#' finite and not larger. If BFGS fails, the Nelder-Mead run is continued
#' with the requested outer loop instead. `...` is passed to `f` and `gr`.
#' @noRd
coco_constr_nm_fit <- function(sta, f, gr, ui, ci, mu, outer.it, outer.eps,
                               control = list(), ...) {
  nm_fit <- function(theta, it, eps, ...) {
    stats::constrOptim(
      theta = theta, f = f, ui = ui, ci = ci, mu = mu, method = "Nelder-Mead",
      control = control, outer.iterations = it, outer.eps = eps,
      hessian = FALSE, ...
    )
  }
  fit <- nm_fit(sta, min(outer.it, 100), max(outer.eps, 1e-5), ...)
  fin <- tryCatch(
    stats::constrOptim(
      theta = fit$par, f = f, grad = gr, ui = ui, ci = ci, mu = mu,
      method = "BFGS", control = control, outer.iterations = outer.it,
      outer.eps = outer.eps, hessian = FALSE, ...
    ),
    error = function(e) NULL
  )
  if (is.null(fin) || !is.finite(fin$value) || any(!is.finite(fin$par))) {
    return(tryCatch(nm_fit(fit$par, outer.it, outer.eps, ...), error = function(e) fit))
  }
  if (isTRUE(fit$value < fin$value)) fit else fin
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
#' back to finite differences. `fg` (optional, called like `f`) returns
#' `list(value, gradient)` of `f` from one evaluation (the C++ AD value and
#' gradient); if given, L-BFGS-B uses it for both its objective and its
#' gradient: optim() always requests the gradient at the point where it has
#' just evaluated the objective, so the gradient is taken from that one AD
#' evaluation instead of a separate likelihood evaluation for the value.
#' The final acceptance test uses `f`.
#' @noRd
coco_polish_fit <- function(fit, f, ui, ci, gr = NULL, ..., fg = NULL,
                            margin = 1e-8) {
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
  if (!is.null(fg)) {
    # value and gradient from one AD evaluation; obj() stores the gradient
    # of the last point it evaluated, obj_gr() returns it for that point
    cache_p <- NULL
    cache_g <- NULL
    obj <- function(p) {
      cache_p <<- p
      cache_g <<- rep(0, n_par)
      if (!feasible(p)) {
        return(penalty)
      }
      r <- tryCatch(do.call(fg, c(list(p), dots)), error = function(e) NULL)
      v <- if (is.null(r) || length(r$value) != 1) NA_real_ else as.numeric(r$value)
      if (!is.finite(v)) {
        return(penalty)
      }
      g <- as.numeric(r$gradient)
      if (length(g) == n_par && all(is.finite(g))) cache_g <<- g
      v
    }
    obj_gr <- function(p) {
      if (!is.null(cache_p) && identical(p, cache_p)) {
        return(cache_g)
      }
      obj(p)
      cache_g
    }
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

#' Default (Nelder-Mead) constrained fit of cocoReg_base and cocoReg_cov
#'
#' In cocoReg_cov, `fallback` is coco_constr_nm_fit() followed by
#' coco_polish_fit(); the beta rows of `ui` (beta_j >= b.beta) are box
#' lower bounds without upper bound, which the box extraction below handles.
#'
#' constrOptim() with Nelder-Mead spends most of the time of a no-covariate
#' fit (and of cocoSoc(), which runs four such fits): every outer barrier
#' iteration restarts the simplex, and near a bound (eta = 0, alpha_i = 0 in
#' misspecified models) the simplex runs to its iteration limit before the
#' L-BFGS-B polish finds the bound. Here the L-BFGS-B stage of
#' coco_polish_fit() (exact AD gradient `gr`) is started at the starting
#' values `sta` directly, and its solution `p` is accepted if it passes a
#' local optimality check with the AD gradient g and Hessian H at `p`
#' (`fgh(p)`, a list with `value`, `gradient`, `hessian` of `f`, as
#' coco_ad_derivs(deriv = 2) returns it):
#'
#' * `p` is feasible and every general (non-box) row of `ui`/`ci` has slack
#'   above `tol_slack` (L-BFGS-B handles these rows only by a penalty);
#' * coordinates within `tol_bound` of a box bound with g pointing out of
#'   the box are fixed; on the remaining (free) coordinates H is positive
#'   definite and the Newton decrement 0.5 g' H^{-1} g (the decrease of `f`
#'   predicted by the local quadratic model) is at most `tol_dec`.
#'
#' Up to `newton_steps` Newton steps on the free coordinates (clipped to the
#' box, kept only if feasible and decreasing `f`) are taken to meet the
#' decrement bound. If the check fails, `fallback()` (the constrOptim() fit
#' followed by coco_polish_fit(), i.e. the former default) is run and the
#' better of the two fits is returned. `fg` (optional) is passed to
#' coco_polish_fit() (value and gradient of `f` from one AD evaluation for the
#' L-BFGS-B stage). Returns list(fit, derivs); `derivs` is
#' `fgh` at fit$par if the check passed (so the caller can reuse the Hessian
#' for the standard errors), else NULL. `...` is passed to `f` and `gr`.
#' @noRd
coco_fit_lbfgsb_checked <- function(sta, f, gr, fgh, ui, ci, fallback, ...,
                                    fg = NULL, margin = 1e-8, tol_dec = 1e-6,
                                    tol_bound = 1e-7, tol_slack = 1e-6,
                                    newton_steps = 2L) {
  dots <- list(...)
  f_eval <- function(p) {
    v <- tryCatch(do.call(f, c(list(p), dots)), error = function(e) NA_real_)
    if (length(v) != 1) NA_real_ else as.numeric(v)
  }
  # former default fit; the candidate `cand` (list(par, value)) is kept if it
  # is strictly better
  run_fallback <- function(cand = NULL) {
    fb <- fallback()
    fb_val <- f_eval(as.numeric(unlist(fb$par, use.names = FALSE)))
    if (!is.null(cand) && is.finite(cand$value) &&
        (!is.finite(fb_val) || cand$value < fb_val)) {
      return(list(fit = cand, derivs = NULL))
    }
    list(fit = fb, derivs = NULL)
  }

  sta <- as.numeric(sta)
  ui <- as.matrix(ui)
  n_par <- length(sta)
  if (ncol(ui) != n_par) {
    return(run_fallback())
  }
  feasible <- function(p) all(ui %*% p - ci >= 0)
  if (!feasible(sta)) {
    return(run_fallback())
  }
  val_sta <- f_eval(sta)
  if (!is.finite(val_sta)) {
    return(run_fallback())
  }

  # box bounds as in coco_polish_fit(); rows with several non-zero entries
  # are the general constraints
  lower <- rep(-Inf, n_par)
  upper <- rep(Inf, n_par)
  general <- logical(nrow(ui))
  for (j in seq_len(nrow(ui))) {
    nz <- which(ui[j, ] != 0)
    if (length(nz) == 1) {
      bound <- ci[j] / ui[j, nz]
      if (ui[j, nz] > 0) {
        lower[nz] <- max(lower[nz], bound + margin)
      } else {
        upper[nz] <- min(upper[nz], bound - margin)
      }
    } else {
      general[j] <- TRUE
    }
  }
  if (any(lower > upper)) {
    return(run_fallback())
  }

  fit <- coco_polish_fit(list(par = sta, value = val_sta), f, ui, ci, gr = gr,
                         ..., fg = fg, margin = margin)
  p <- as.numeric(unlist(fit$par, use.names = FALSE))
  fp <- f_eval(p)
  if (length(p) != n_par || any(!is.finite(p)) || !feasible(p) || !is.finite(fp)) {
    return(run_fallback())
  }

  for (step in 0:newton_steps) {
    d <- tryCatch(fgh(p), error = function(e) NULL)
    if (is.null(d)) break
    g <- as.numeric(d$gradient)
    H <- as.matrix(d$hessian)
    if (length(g) != n_par || any(!is.finite(g)) || any(dim(H) != n_par) ||
        any(!is.finite(H))) break
    if (any(general) &&
        any(ui[general, , drop = FALSE] %*% p - ci[general] <= tol_slack)) break

    fixed <- ((p - lower <= tol_bound) & (g > 0)) |
      ((upper - p <= tol_bound) & (g < 0))
    free <- which(!fixed)
    if (length(free) == 0) {
      return(list(fit = list(par = p, value = fp), derivs = d))
    }
    R <- tryCatch(chol(H[free, free, drop = FALSE]), error = function(e) NULL)
    if (is.null(R)) break
    dir <- backsolve(R, forwardsolve(t(R), g[free]))  # H_free^{-1} g_free
    dec <- 0.5 * sum(g[free] * dir)
    if (is.finite(dec) && dec <= tol_dec) {
      return(list(fit = list(par = p, value = fp), derivs = d))
    }
    if (step == newton_steps || !is.finite(dec)) break

    q <- p
    q[free] <- p[free] - dir
    q <- pmin(pmax(q, lower), upper)
    if (!feasible(q)) break
    fq <- f_eval(q)
    if (!is.finite(fq) || fq >= fp) break
    p <- q
    fp <- fq
  }
  run_fallback(list(par = p, value = fp))
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

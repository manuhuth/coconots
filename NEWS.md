# coconots 2.0.5

* `installJuliaPackages()` installs or updates Coconots.jl to version 0.2 or newer
  (an existing 0.1.x installation was previously kept).

## Bug fixes

* `cocoReg(julia = TRUE)` without covariates passed a user-supplied `start` to Julia
  in R order (`lambda` first) although the Julia backend expects `lambda` last. `start`
  is now reordered before the call. The `start` documentation states the expected order.
* `cocoSim()` with covariates takes `lambda_t` from row `t` of `xreg` (previously row
  `t - order`); GP2 start values account for `eta`.
* `cocoPit()` evaluates the PIT over `t = order + 1, ..., T` on both backends, uses the
  upper-tail chi-square p-value, and draws the confidence bands at the level `conf.alpha`
  (order 1 previously used the model `alpha`).
* Densities, the order-2 normalizer `P(y, z)` and the order-1 thinning term return the
  exact limits at zero rates and at `alpha = 0` and `alpha = 1` instead of `0 * Inf = NaN`.
* Inverse-transform samplers terminate at the count limit 170.
* Input checks: counts above 170 (factorial overflow), rank-deficient `xreg`, series too
  short for the model order, and `lambda_t <= 0` now raise informative errors.

## Performance

* The default RCPP fit uses L-BFGS-B with gradient and Hessian from C++ forward-mode
  automatic differentiation and an optimality check (KKT conditions, positive definite
  reduced Hessian, Newton decrement); `constrOptim` Nelder-Mead with a polishing step
  remains the fallback. Starting values for covariate coefficients come from a Poisson GLM.
* Faster likelihoods (sized pmf tables, recurrences, factored order-2 convolution),
  `cocoBoot()`, `cocoSim()` with covariates and `predict(k > 1)` (all replicates or paths
  in one C++ call), `cocoScore()` and `cocoResid()` (one pmf per distinct history), and
  fewer R-Julia round trips in `cocoReg(julia = TRUE)` and `cocoScore(julia = TRUE)`.

## Other changes

* `method.hessian` is deprecated and ignored with a warning; the Hessian is computed by
  automatic differentiation on both backends. The `numDeriv` dependency is removed.

# coconots 2.0.4

## Bug fixes

* Fixed a heap buffer over-read reported by CRAN's valgrind check. The post-hoc
  log-likelihood for Poisson models with covariates reused the Generalized-Poisson
  parameter slicing, which drops `eta`. Poisson models carry no `eta`, so the slicing
  also discarded the first covariate coefficient and passed a coefficient vector one
  element shorter than the number of covariate columns to `likelihoodGP1cov()` /
  `likelihoodGP2cov()`, which then read one `double` past the end of that vector on
  every observation. The stored `$likelihood` (and hence AIC/BIC) was wrong for these
  models; parameter estimates and standard errors were unaffected.
* `likelihoodGP2cov()` clamped an underflowing per-observation likelihood to
  `1/10^-12 = 1e+12`, rewarding the optimizer for entering the degenerate region.
  It now clamps to `1e-12`.
* The order-2 starting-value shrink loops assigned to a dead variable, so `alpha3`
  never shrank and the loop could fail to terminate.
* `cocoReg()` with a `data.frame` `xreg` overwrote `data` instead of coercing `xreg`.
* `cocoSoc()` now forwards its `julia` argument to the internal `cocoReg()` fits and
  to `cocoScore()`, instead of hardcoding `julia = TRUE` when scoring.
* Julia result dictionaries are now read by key name rather than by integer position,
  which depended on Julia's hash iteration order and could shift silently across
  Julia versions.

## Performance

* `cocoReg()` gained a live `cores` argument. The order-2 likelihood kernels of the
  RCPP backend are multithreaded via `std::thread`; accumulation stays serial in the
  original order, so results are bit-identical for any thread count. `cores` defaults
  to the number of physical cores minus one, capped at 2 while `_R_CHECK_LIMIT_CORES_`
  is set.
* Estimation kernels and the assessment tools (`cocoScore()`, `cocoPit()`,
  `cocoResid()`, `cocoBoot()`, `rgenpois()`) were rewritten around batched density
  evaluations, lookup tables and memoization. Numerical results are unchanged
  (order-1 paths exactly, order-2 paths to ~1 ulp).
* `JULIA_NUM_THREADS` is now set to the number of physical cores minus one before the
  first `JuliaConnectoR` call, unless the user already set it. A Julia process
  otherwise starts single-threaded, leaving `Coconots.jl`'s threaded likelihood idle.

# coconots 2.0.3

* Added a `softplus` link function for the conditional mean of the innovation rate
  in `cocoReg()` (covariate models). Uses a numerically stable formulation, is smooth
  everywhere, and guarantees a positive rate.
* Added S3 `summary` methods for forecast objects (`cocoForecast`,
  `cocoForecastCollection`) returning a data frame of point forecasts (mean, median,
  mode) and prediction intervals per forecast horizon.

# coconots 2.0.1

* Added a `NEWS.md` file to track changes to the package.
* Added new S3 methods for an object created with cocoReg.

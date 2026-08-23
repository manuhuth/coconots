# Changelog

## coconots 2.0.4

### Bug fixes

- Fixed a heap buffer over-read reported by CRAN’s valgrind check. The
  post-hoc log-likelihood for Poisson models with covariates reused the
  Generalized-Poisson parameter slicing, which drops `eta`. Poisson
  models carry no `eta`, so the slicing also discarded the first
  covariate coefficient and passed a coefficient vector one element
  shorter than the number of covariate columns to `likelihoodGP1cov()` /
  `likelihoodGP2cov()`, which then read one `double` past the end of
  that vector on every observation. The stored `$likelihood` (and hence
  AIC/BIC) was wrong for these models; parameter estimates and standard
  errors were unaffected.
- `likelihoodGP2cov()` clamped an underflowing per-observation
  likelihood to `1/10^-12 = 1e+12`, rewarding the optimizer for entering
  the degenerate region. It now clamps to `1e-12`.
- The order-2 starting-value shrink loops assigned to a dead variable,
  so `alpha3` never shrank and the loop could fail to terminate.
- [`cocoReg()`](https://manuhuth.github.io/coconots/reference/cocoReg.md)
  with a `data.frame` `xreg` overwrote `data` instead of coercing
  `xreg`.
- [`cocoSoc()`](https://manuhuth.github.io/coconots/reference/cocoSoc.md)
  now forwards its `julia` argument to the internal
  [`cocoReg()`](https://manuhuth.github.io/coconots/reference/cocoReg.md)
  fits and to
  [`cocoScore()`](https://manuhuth.github.io/coconots/reference/cocoScore.md),
  instead of hardcoding `julia = TRUE` when scoring.
- Julia result dictionaries are now read by key name rather than by
  integer position, which depended on Julia’s hash iteration order and
  could shift silently across Julia versions.

### Performance

- [`cocoReg()`](https://manuhuth.github.io/coconots/reference/cocoReg.md)
  gained a live `cores` argument. The order-2 likelihood kernels of the
  RCPP backend are multithreaded via `std::thread`; accumulation stays
  serial in the original order, so results are bit-identical for any
  thread count. `cores` defaults to the number of physical cores minus
  one, capped at 2 while `_R_CHECK_LIMIT_CORES_` is set.
- Estimation kernels and the assessment tools
  ([`cocoScore()`](https://manuhuth.github.io/coconots/reference/cocoScore.md),
  [`cocoPit()`](https://manuhuth.github.io/coconots/reference/cocoPit.md),
  [`cocoResid()`](https://manuhuth.github.io/coconots/reference/cocoResid.md),
  [`cocoBoot()`](https://manuhuth.github.io/coconots/reference/cocoBoot.md),
  `rgenpois()`) were rewritten around batched density evaluations,
  lookup tables and memoization. Numerical results are unchanged
  (order-1 paths exactly, order-2 paths to ~1 ulp).
- `JULIA_NUM_THREADS` is now set to the number of physical cores minus
  one before the first `JuliaConnectoR` call, unless the user already
  set it. A Julia process otherwise starts single-threaded, leaving
  `Coconots.jl`’s threaded likelihood idle.

## coconots 2.0.3

CRAN release: 2026-06-14

- Added a `softplus` link function for the conditional mean of the
  innovation rate in
  [`cocoReg()`](https://manuhuth.github.io/coconots/reference/cocoReg.md)
  (covariate models). Uses a numerically stable formulation, is smooth
  everywhere, and guarantees a positive rate.
- Added S3 `summary` methods for forecast objects (`cocoForecast`,
  `cocoForecastCollection`) returning a data frame of point forecasts
  (mean, median, mode) and prediction intervals per forecast horizon.

## coconots 2.0.1

CRAN release: 2025-07-24

- Added a `NEWS.md` file to track changes to the package.
- Added new S3 methods for an object created with cocoReg.

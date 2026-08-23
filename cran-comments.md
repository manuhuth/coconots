## Submission

This is a patch release (2.0.4). The previously published version is 2.0.3.

It addresses the valgrind issue reported for 2.0.3 under "Additional issues" on
<https://cran.r-project.org/web/checks/check_results_coconots.html>:

```
Invalid read of size 8
  at inner_product<double*, double*, double> (stl_numeric.h:198)
  by innerProduct(...)        src/cfunctions.cpp:36
  by likelihoodGP2cov(...)    src/cfunctions.cpp:528
Address ... is 0 bytes after a block of size 7,960 alloc'd
ERROR SUMMARY: 298 errors from 1 contexts
```

Cause: the post-hoc log-likelihood evaluation for Poisson models with covariates
(`R/cocoReg_cov.R`) reused the parameter slicing of the Generalized-Poisson models,
which drops the overdispersion parameter `eta`. Poisson models carry no `eta`, so the
slicing also discarded the first covariate coefficient and passed a coefficient vector
one element shorter than the number of covariate columns to the C++ kernel. The kernel
iterates over the covariate row and reads the corresponding coefficient, so it read one
`double` past the end of that vector once per observation (298 observations in the
affected test, matching the reported error count).

Fix: `pars[-c(1, 2)]` -> `pars[-1]` for the order-1 case and
`pars[-c(1:4)]` -> `pars[-c(1:3)]` for the order-2 case. Parameter estimates and
standard errors were unaffected; the stored `$likelihood` (and therefore AIC/BIC) for
Poisson models with covariates was wrong and is now correct. We verified with an
assertion on every covariate-kernel call that the coefficient vector and the number of
covariate columns now agree for all four model variants.

The release also contains further bug fixes and a substantial speed-up of the
estimation kernels and assessment tools, including an optional multithreaded order-2
likelihood; see NEWS.md. The multithreading uses `std::thread` (no new dependency) and
respects `_R_CHECK_LIMIT_CORES_` by capping the thread count at 2 while it is set.

## Test environments

* local macOS (darwin), R 4.5.2 (2025-10-31)
* GitHub Actions: ubuntu-latest, macOS-latest, windows-latest;
  R release, devel, oldrel-1; Julia 1.11

## R CMD check results

0 errors | 0 warnings | 2 notes

Both notes are specific to the local check environment and do not reflect package
issues:

* "Files 'README.md' or 'NEWS.md' cannot be checked without 'pandoc' being installed."
  (pandoc is not installed on the local machine; it is available on CRAN.)
* "Skipping checking HTML validation: 'tidy' doesn't look like recent enough HTML Tidy
  ... package 'V8' unavailable." (local tooling only.)

## Submission

This is a minor update (2.0.3). The previously published version is 2.0.2.

## Changes since 2.0.2

* New `softplus` link function for the conditional mean of the innovation rate in
  `cocoReg()`. It uses a numerically stable formulation, is smooth everywhere, and
  guarantees a positive rate.
* New S3 `summary` methods for forecast objects (`cocoForecast`,
  `cocoForecastCollection`) returning a data frame of point forecasts (mean, median,
  mode) and prediction intervals per forecast horizon.

These changes were made following reviewer feedback for the Journal of Statistical
Software submission, so that the CRAN version stays in sync with the accompanying
manuscript.

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

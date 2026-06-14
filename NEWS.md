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

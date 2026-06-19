# Simulation of Count Time Series

The function generates a time series of low counts from the (G)PAR model
class for a specified innovation distribution, sample size, lag order,
and parameter values.

## Usage

``` r
cocoSim(
  type,
  order,
  par,
  length,
  xreg = NULL,
  init = NULL,
  julia = FALSE,
  julia_seed = NULL,
  link_function = "log"
)
```

## Arguments

- type:

  character, either "Poisson" or "GP" indicating the type of the
  innovation distribution

- order:

  integer, either 1 or 2 indicating the order of the model

- par:

  numeric vector, the parameters of the model, the number of elements in
  the vector depends on the type and order specified.

- length:

  integer, the number of observations in the generated time series

- xreg:

  data frame of control variables (defaul: NULL)

- init:

  numeric vector, initial data to use (default: NULL). See details for
  more information on the usage.

- julia:

  If TRUE, the julia implementation is used. In this case, init is
  ignored but it might be faster (default: FALSE).

- julia_seed:

  Seed for the julia implementation. Only used if julia equals TRUE.

- link_function:

  Specifies the link function for the conditional mean of the innovation
  (\\\lambda\\). The default is `log`, but other available options
  include `identity` and `relu`. This parameter is applicable only when
  covariates are used. Note that using the `identity` link function may
  result in \\\lambda\\ becoming negative. To prevent this, ensure all
  covariates are positive and restrict the parameter \\\beta\\ to
  positive values.

## Value

a vector of the simulated time series

## Details

The function checks for valid input of the type, order, parameters, and
initial data before generating the time series.

The init parameter allows users to set a custom burn-in period for the
simulation. By default, when simulating with covariates, no burn-in
period is specified since there is no clear choice on the covariates.
However, the init argument gives users the flexibility to select an
appropriate burn-in period for the covariate case. One way to do this is
to simulate a time series using `cocoSim` with appropriate covariates
and pass the resulting time series to the init argument of a new
`cocoSim` run so that the first time series is used as the burn-in
period. If init is not specified for the covariate case, a warning will
be returned to prompt the user to specify a custom burn-in period. This
helps ensure that the simulation accurately captures the dynamics of the
system being modeled.

## Author

Manuel Huth

## Examples

``` r
#First Order Model
lambda <- 1
alpha <- 0.4
set.seed(12345)

# Simulate using the RCPP implementation
data_rcpp <- cocoSim(order = 1, type = "Poisson", par = c(lambda, alpha), length = 100)

# Second Order Model
lambda <- 1
alpha_1 <- 0.3
alpha_2 <- 0.1
alpha_3 <- 0.2
eta <- 0.2
data_j <- cocoSim(order = 2, type = "GP",
                  par = c(lambda, alpha_1, alpha_2,alpha_3,eta),
                  length = 100)
```

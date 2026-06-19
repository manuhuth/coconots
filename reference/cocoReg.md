# Fitting First- and Second Order (G)PAR Models

The function fits first- and second-order (Generalized) Poisson integer
autoregressive `(G)PAR` time series models for count data as discussed
in Jung and Tremayne (2011). Autoregressive dependence on past counts is
modeled using a special random operator that preserves integer values
and, through closure under convolution, ensures that the marginal
distribution remains within the same family as the innovations.

These models can be viewed as stationary finite-order Markov chains,
where the innovation distribution can be either Poisson or Generalized
Poisson, the latter accounting for overdispersion. Estimation is
performed via maximum likelihood, with an option to impose linear
constraints. Without constraints, parameters may fall outside the
theoretically feasible space, but optimization may be faster.

Method of moments estimators are used to initialize numerical
optimization, though custom starting values can be provided. If julia is
installed, users can opt to run the optimization in julia for
potentially faster computation and improved numerical stability via
automatic differentiation. See below for details on the julia
implementation.

## Usage

``` r
cocoReg(
  type,
  order,
  data,
  xreg = NULL,
  constrained.optim = TRUE,
  b.beta = -10,
  start = NULL,
  start.val.adjust = TRUE,
  method_optim = "Nelder-Mead",
  replace.start.val = 1e-05,
  iteration.start.val = 0.6,
  method.hessian = "Richardson",
  cores = 2,
  julia = FALSE,
  julia_installed = FALSE,
  link_function = "log"
)
```

## Arguments

- type:

  character, either "Poisson" or "GP" indicating the type of the
  innovation distribution

- order:

  integer, either 1 or 2 indicating the order of the model

- data:

  time series data to be used in the analysis

- xreg:

  optional matrix of explanatory variables (without constant term) for
  use in a regression model

- constrained.optim:

  logical indicating whether optimization should be constrained,
  currently only available in the R version

- b.beta:

  numeric value indicating the lower bound for the parameters of the
  explanatory variables for the optimization, currently only available
  in the R version

- start:

  optional numeric vector of starting values for the optimization

- start.val.adjust:

  logical indicating whether starting values should be adjusted,
  currently only available in the R version

- method_optim:

  character string indicating the optimization method to be used,
  currently only available in the R version. In the julia implementation
  this is by default the LBFGS algorithm

- replace.start.val:

  numeric value indicating the value to replace any invalid starting
  values, currently only available in the R version

- iteration.start.val:

  numeric value indicating the proportion of the interval to use as the
  new starting value, currently only available in the R version

- method.hessian:

  character string indicating the method to be used to approximate the
  Hessian matrix, currently only available in the R version

- cores:

  numeric indicating the number of cores to use, currently only
  available in the R version (default: 2)

- julia:

  if TRUE, the model is estimated with julia. This can improve
  computational speed significantly since julia makes use of derivatives
  using autodiff. In this case, only `type`, `order`, `data`, `xreg`,
  and `start` are used as other inputs (default: FALSE).

- julia_installed:

  if TRUE, the model R output will contain a julia compatible output
  element.

- link_function:

  Specifies the link function for the conditional mean of the innovation
  (\\\lambda\\). The default is `log`, but other available options
  include `identity`, `relu`, and `softplus`. This parameter is
  applicable only when covariates are used. Note that using the
  `identity` link function may result in \\\lambda\\ becoming negative.
  To prevent this, ensure all covariates are positive and restrict the
  parameter \\\beta\\ to positive values by setting `b.beta` to a small
  positive value.

## Value

an object of class coco. It contains the parameter estimates, standard
errors, the log-likelihood, and information on the model specifications.
If julia is used for parameter estimation or the julia installation
parameter is set to TRUE, the results contain an additional Julia
element that is called from the model julia assessment tools if they are
run with the julia implementation.

## Details

Let a time series of counts be \\\\X_t\\\\ and be \\R(\cdot)\\ a random
operator that differs between model specifications. For more details on
the random operator, see Jung and Tremayne (2011) and Joe (1996). The
general first-order model is of the form \$\$X_t = R(X\_{t-1}) +
W_t,\$\$ and the general second-order model of the form \$\$X_t =
R(X\_{t-1}, X\_{t-2}) + W_t,\$\$ where \\W_t\\ are i.i.d Poisson (\\W_t
\sim Po(\lambda_t)\\) or Generalized Poisson (\\W_t \sim GP(\lambda_t,
\eta)\\) innovations. Through closure under convolution the marginal
distributions of \\\\X_t\\\\ are therefore Poisson or Generalized
Poisson distributions, respectively.

If no covariates are used \\\lambda_t = \lambda\\ and if covariates are
used \$\$g(\lambda_t) = \left(\beta_0 + \sum\_{j = 1}^k \beta_j \cdot
z\_{t,j} \right),\$\$ whereby \\z\_{t,j}\\ is the \\j\\-th covariate at
time \\t\\ and \\g\\ is a link function. Currently supported link
functions are: the identity \\g(x) = x\\, the logarithmic link \\g(x) =
\ln x\\, relu \\g(x) = \max(x, \epsilon)\\ for small \\\epsilon \> 0\\,
and softplus \\g(x) = \log(\exp(x) - 1)\\. To ensure positivity of
\\\lambda\\ if the identity function is used, \\\beta_j, z\_{t,j} \> 0\\
must be enforced. Alternatively, computational values of \\\lambda \leq
0\\ can be set to a small positive value. This option is named 'relu',
due to its similarity to a ReLU function commonly used in machine
learning. The softplus link uses a numerically stable implementation of
its inverse \\\log(1 + \exp(x))\\, ensuring \\\lambda \> 0\\ and smooth
gradients everywhere.

Standard errors are computed by the square root of the diagonal elements
of the inverse Hessian.

This function is implemented in two versions. The default runs on RCPP.
An alternative version uses a julia implementation which can be chosen
by setting the argument julia to TRUE. In order to use this feature, a
running julia installation is required on the system. The RCPP
implementation uses the derivative-free Nelder-Mead optimizer to obtain
parameter estimates. The julia implementation makes use of julia's
automatic differentiation in order to obtain gradients such that it can
use the LBFGS algorithm for optimization. This enhances the numeric
stability of the optimization and yields an internal validation if both
methods yield qualitatively same parameter estimates. Furthermore, the
julia implementation can increase the computational speed significantly,
especially for large models.

The model assessment tools
[`cocoBoot`](https://manuhuth.github.io/coconots/reference/cocoBoot.md),
[`cocoPit`](https://manuhuth.github.io/coconots/reference/cocoPit.md),
and
[`cocoScore`](https://manuhuth.github.io/coconots/reference/cocoScore.md)
will use a julia implementation as well, if the `cocoReg` was run with
julia. Additionally, one can make the RCPP output of `cocoReg`
compatible with the julia model assessments by setting `julia_installed`
to true. In this case, the user can choose between the RCPP and the
julia implementation for model assessment.

## References

Jung, R. C. and Tremayne, A. R. (2011) Convolution-closed models for
count time series with applications. *Journal of Time Series Analysis*,
**32**, 268–280.

Joe, H. (1996) Time series models with univariate margins in the
convolution-closed infinitely divisible class. *Journal of Applied
Probability*, **33**, 664–677.

## Author

Manuel Huth

## Examples

``` r
## GP2 model without covariates
length <- 1000
par <- c(0.5,0.2,0.05,0.3,0.3)
data <- cocoSim(order = 2, type = "GP", par = par, length = length)
fit <- cocoReg(order = 2, type = "GP", data = data)
#> Warning: The initial value of alpha2 is not in the feasible region and was set to 1e-05 

##Poisson1 model with covariates
length <- 1000
period <- 12
sin <- sin(2*pi/period*(1:length))
cos <- cos(2*pi/period*(1:length))
cov <- cbind(sin, cos)
par <- c(0.2, 0.2, -0.2)
data <- cocoSim(order = 1, type = "Poisson", par = par, xreg = cov, length = length)
#> Warning: No burn-in period is specified using the init argument. Hence, the resulting simulated time series might not be stationary. You can add a custom burn-in period by passing it to the init argument. This could be, for example, done by simulating a burn-in period with appropriate covariates using cocoSim and passing the resulting time series to the init argument of a new cocoSim run.
fit <- cocoReg(order = 1, type = "Poisson", data = data, xreg = cov)
```

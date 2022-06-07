\name{cocoReg}
\alias{cocoReg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fitting (Generalized) Poisson Autoregressive models
}
\description{
The function \code{cocoReg} fits first and second order (Generalized) Poisson Autoregressive (G)PAR models presented in  (Jung and Tremayne, 2010). Lags enter the process via a random operator that accounts for the count structure of the data (Joe, 1996). The innovations can either follow a Poisson (PAR) or a Generalized Poisson (GPAR) distribution. The models can be thought of as stationary Markov chains of finite order. The distribution of the innovations can either be Poisson or Generalized Poisson.

Stochastic seasonality can be implemented by allowing arbitrary lags which are used for the random operator.

Estimation is done by maximum likelihood and can either include linear constraines or be unconstrained. If unconstrained estimation is chosen, it cannot be guaranteed that the parameters lie in the theoretically feasible parameter space. However the speed of the optimization can be increased.

Starting values for numerical optimization are obtained by method of moments estimators but can be set manually if desired.

}

possible extensions: identity as covariate link
\usage{
cocoReg(type = c("Poisson", "GP"), order = c(1,2), data, xreg = NULL, 
        seasonality = c(1, 2), mu = 1e-4,outer.it = 500, outer.eps = 1e-10,
        optim_control = FALSE, constrained.optim = TRUE, start = NULL,
        start.val.adjust = TRUE, replace.start.val = 1e-5, 
        iteration.start.val = 0.99)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{type}{distributon of the Innovation. Either "Poisson" or "GP" (Generalized Poisson).}
  \item{order}{either 1 or 2. Order of the autoregressive process.}
  \item{data}{vector or data frame. A univariate time series of counts.}
  \item{xreg}{Null, matrix or data frame. If NULL, no covariates are included. Otherwise, a the obejct must contain covariates as columns and the number of rows must equal the length of \code{data}.}
  \item{seasonality}{vector containing two non-negative integers for order two models and one non-negative integer for order one models. If two integers are used but a model of order one is fitted, only the first integer is used. }
  \item{mu}{small positive tuning paramater. See \code{constrOptim} for details.}
  \item{outer.it}{positive integer value. Controls the iterations of the barrier algorithm. See \code{constrOptim} for details.}
  \item{outer.eps}{small positive number. Relative convergence tolerance of the barrier algorithm. See \code{constrOptim} for details.}
  \item{constr}{bool. If TRUE, the optimization is performed unconstrained.}
  \item{start}{a vector with custom starting values. Number must equal the number of parameters. If NULL, moment conditions are used to compute starting values.}
  \item{start.adj}{bool. If TRUE, starting values are adjusted such that they fit the linear constraints.}
  \item{replst}{small positive value. If \code{start.adj} = TRUE, starting values smaller than zero are set to \code{replst} and starting values for the autoregressive parameters and the overdispersion parameter greater than 1 are set to 1-\code{replst}.}
  \item{it.sv}{float close to but smaller than one. Adjusts the linear constraints by multiplying any value iterative by \code{it.sv} until they are satisfied.}

}

\details{
\emph{General introduction to the methods with formulas (math-heaviest part of the vignette)}

%Let a time series of counts be \eqn{\{X_t\}}, a random operator be \eqn{R(\cdot)}, the relevant past history of \eqn{X_t} be \eqn{\mathcal{F}_{t-1}} and \eqn{W_t} be i.i.d discrete innovations. The general model is of the form
%\deqn{X_t = R(\mathcal{F}_{t-1}) + W_t.}
%For first-order and second-order models the relevant past history depends on the model's order and is defined by \eqn{\mathcal{F}_{t-1} = X_{t-s_1}} and \eqn{\mathcal{F}_{t-1} = (X_{t-s_1}, X_{t-s_2})}, respectively. \eqn{s_1} < \eqn{s_2} are the parameters indicating the degree of stochastic seasonality.

%The innovations follow either a Poisson or a generalized Poisson distribution. Such that either \eqn{W_t \sim Pois(\lambda_t)} or \eqn{W_t \sim GP(\lambda_t, \eta)}.
%If no covariates are used \eqn{\lambda_t = \lambda} and if covariates are used \deqn{\lambda_t = \exp{\left(\beta_0 + \sum_{j = 1}^p \beta_j \cdot y_{t,j} \right)},} whereby \eqn{y_{t,j}} is the \eqn{j}-th covariate at time \eqn{t}.

%Standard errors are computed by the square root of the diagonal elements of the inverse hessian.

}
\value{
Returns a named list of class "\code{coco.fit}"" or "\code{coco.fit.c}" containing the following elements:
\item{par}{fitted paramaters.}
\item{gradient}{numerical gradient evaluated at the fitted parameters.}
\item{hessian}{numerical hessian evaluated at the fitted parameters.}
\item{inv hessian}{inverse of the hessian.}
\item{se}{asymptotic standard erros. Obtained as the sqaure root of the diagonals of the inverse hessian.}
\item{ts}{time series used to fit the model.}
\item{type}{innovation distribution used to fit the model.}
\item{order}{order of autoregressive process that is used to fit the model.}
\item{seasonality}{stochastic seasonality.}
\item{likelihood}{value of the log-likelihood evaluated at the fitted parameters.}
\item{duration}{time it took to fit the model.}

}
\references{
Jung, R. C. and Tremayne, A. R. (2010) Convolution-closed models for count timeseries with applications. \emph{Journal of Time Series Analysis}, \bold{32}, 3, 268--280.

Joe, H. (1996) Time series models with univariate margins in the convolution-closed infinitely divisible class. \emph{Journal of Applied Probability}, 664--677.
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## GP2 model without covariates
size <- 1000
par <- c(0.5,0.2,0.05,0.3,0.3)
data.sim <- cocoSim(order = 2, type = "GP", par = par, size = size)
data <- data.sim$data
fit <- cocoReg(order = 2, type = "GP", data = data)

##Poisson1 model with covariates
size <- 1000
period <- 50
sin <- sin(2*pi/period*(1:size))
cos <- cos(2*pi/period*(1:size))
cov <- cbind(sin, cos)
par <- c(0.2, 0.1, 0.2, -0.2)
data.sim <- cocoSim(order = 1, type = "Poisson", par = par, xreg = cov, size = size)
data <- data.sim$data
fit <- cocoReg(order = 1, type = "Poisson", data = data, xreg = cov)
}


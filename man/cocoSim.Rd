\name{cocoSim}
\alias{cocoSim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulation PAR1, PAR2, GP1 and GP2 Models

}
\description{
Created with Rcpp (C++)

findings/to-do:

1) find theory for burn-in period. Currently: random (Generalized) Poisson variables;

2) if values get around 180, the values rise to around 170 and fall down in the next period approx. to lambda. (found out that only values up to (170!) can be computed and hence the function 'fact' produces NaN. (restrict the values with if function?)

possible extensions: ?
}
\usage{
cocoSim(type, order, par, length, xreg=NULL, seasonality = c(1,2), init = NULL)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{type}{distributon of the Innovation. Poisson or Generalized Poisson (GP).}
  \item{order}{integer that indicates the order of the autoregressive process.}
  \item{par}{vector containing the parameters used for the simulation. See details for the order of the parameters.}
  \item{xreg}{matrix containing covariates. Only in \code{coco.reg.c}. Number of rows must equal length of \code{data}.}
  \item{length}{positive integer value greater than the highest value of \code{seasonality}. Sample size of time series.}
  \item{seasonality}{vector containing two non-negative integers for order two models and one non-negative integer for order one models. If two integers are used but a model of order one is simulated, only the first integer is used.}
  \item{init}{vector of positive integer values. Can be used to pass inital values for the simulated time series. (Could also be another simulated time series with different model structrue)}
}
\details{
For the \code{coco.sim} function the order of the parameters is \eqn{\lambda, \alpha, \eta} and for the \code{coco.sim.c} function with \code{k} covariates the order of the parameters is \eqn{\alpha, \eta, \beta_0, \beta_1, \dots, \beta_k}. For models of order two \eqn{\alpha} is replaced by \eqn{\alpha_1, \alpha_2, \alpha_3} and for models with Poisson innovations \eqn{\eta} is left out.

}
\value{
a list containing the following elements:
\item{time}{time it needed to simulate the time series.}
\item{data}{simulated time series.}
}
\references{
%% ~put references to the literature/web site here ~
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
ssize <- 1000
par <- c(0.5,0.2,0.05,0.3,0.3)
data.sim <- cocoSim(order = 2, type = "GP", par = par, length = ssize)
data <- data.sim$data


##Poisson1 model with covariates
ssize <- 1000
period <- 50
sin <- sin(2*pi/period*(1:ssize))
cos <- cos(2*pi/period*(1:ssize))
cov <- cbind(sin, cos)
par <- c(0.2, 0.1, 0.2, -0.2)
data.sim <- cocoSim(order = 1, type = "Poisson", par = par, xreg = cov, length = ssize)
data <- data.sim$data


}

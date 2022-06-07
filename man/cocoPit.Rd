\name{cocoPit}
\alias{cocoPit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Probability Integral Tansformation (PIT)
}
\description{
The function provides the non-randomized PIT histogram of Czado et al. (2009) for model assessment.
}
\usage{
pit(coco, J = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{coco}{object of class \code{coco.fit} or \code{coco.fit.c}}
  \item{J}{positive integer. Number of bins used to compute the PIT.}
}
\details{
\emph{Usage description}

\emph{Brief math without formulas}
}
\value{
a list containing the following elements:
\item{type}{innovation distribution used to fit the model.}
\item{order}{order of autoregressive process that is used to fit the model.}
\item{ts}{time series used to fit the model.}
\item{cov}{covariates used to fit the model.}
\item{par}{fitted paramaters used to compute the PIT.}
\item{PIT}{plot containing the PIT.}
\item{PIT values}{values of the PIT.}
\item{duration}{time it took to compute the PIT.}
}
\references{
Czado, C., Gneiting, T. and Held, L. (2009) Predictive model assessment for count data. \emph{Biometrics} \bold{65}, 1254--61.

Jung, R. C. and Tremayne, A. R. (2010) Convolution-closed models for count time series with applications. \emph{Journal of Time Series Analysis}, \bold{32}, 3, 268--280.
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{cocoReg}} for fitting an object of class \code{coco.fit} or \code{coco.fit.c} that can be passed to the function.
}
\examples{
## GP2 model without covariates
size <- 1000
par <- c(0.5,0.2,0.05,0.3,0.3)
data.sim <- cocoSim(order = 2, type = "GP", par = par, size = size)
data <- data.sim$data
fit <- cocoReg(order = 2, type = "GP", data = data)
pit <- cocoPit(fit)

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
pit <- cocoPit(fit)
}

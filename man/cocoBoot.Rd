\name{cocoBoot}
\alias{cocoBoot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Parametric Bootstrap
}
\description{

The function provides the parametric bootstrap of Tsay (1992) for model assessment.



}
\usage{
cocoBoot(coco, numb.lags = 21, rep.Bootstrap = 400, confidence = 0.95)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{coco}{object of class \code{coco.fit} or \code{coco.fit.c}}
   \item{numb.lags}{positive integer. Number of lags that are examined within the bootstrap.}
   \item{rep.Bootstrap}{positive integer number. Number of parametric replications of the time series.}
   \item{confidence}{typically one of 0.90, 0.95, 0.99. Confidence level for the confidence intervals.}
}
\details{
\emph{Usage description}

\emph{Brief math without formulas?}
}
\value{
a list containing the following elements:
\item{type}{innovation distribution used to fit the model.}
\item{order}{order of autoregressive process that is used to fit the model.}
\item{ts}{time seris used to fit the model.}
\item{cov}{covariates used to fit the model.}
\item{PBT.plot}{plot containing the PBT results.}
\item{confidence}{Confidence level for the confidence intervals}
\item{duration}{time it took to compute the PBT.}

}
\references{
Tsay, R. S. (1992) Model checking via parametric bootstraps in time series analysis. \emph{Applied Statistics} \bold{41}, 1--15.
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
boot <- cocoBoot(fit)

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
boot <- coco.Boot(fit)
}


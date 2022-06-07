\name{cocoResid}
\alias{cocoResid}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
(Pearson) Residuals and Fitted Values
}
\description{

This function provides the (Pearson) residuals and fitted values for a model of class \code{coco.fit} or \code{coco.fit.c}.

}
\usage{
cocoResid(coco, val.num = 1e-5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{coco}{object of class \code{coco.fit} or \code{coco.fit.c}}
  \item{val.num}{small positive number. Stopping value for numerical approximation of expected values and variance.}
}
\details{
\emph{Usage description}

\emph{Brief math with formulas?}}
\value{
a list containing the following elements:
\item{cond.expec}{contains the conditional expectations.}
\item{cond.var}{contains the conditional variances.}
\item{residuals}{contains the residuals.}
\item{pe.resid}{contains the pearson residuals.}
\item{type}{innovation distribution used to fit the model.}
\item{order}{order of autoregressive process that is used to fit the model.}
\item{ts}{time seris used to fit the model.}
\item{cov}{covariates used to fit the model.}
\item{par}{fitted paramaters used to compute the residuals.}
\item{duration}{time it took to compute the residuals and fitted values.}


}
\references{
Jung, R. C. and Tremayne, A. R. (2010) Convolution-closed models for count timeseries with applications. \emph{Journal of Time Series Analysis}, \bold{32}, 3, 268--280.
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
data.sim <- cocosim(order = 2, type = "GP", par = par, size = size)
data <- data.sim$data
fit <- cocoReg(order = 2, type = "GP", data = data)
resid <- cocoResid(fit)

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
resid <- cocoResid(fit)
}

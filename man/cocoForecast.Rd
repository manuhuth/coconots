\name{cocoForecast}
\alias{cocoForecast}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
One step ahead forecast.
}
\description{
This function yields the probability mass of the one-step ahead forecast.
}
\usage{
cocoForecast(coco, max=10, xcast=NULL, plot = TRUE, title = "Probability mass",
             xlab="x", ylab="Probabilities", width_bars = 0.04,
             seasonality=c(1,2), decimals = 4)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{coco}{object of class \code{coco.fit} or \code{coco.fit.c}}
   \item{max}{positive integer. Maximum value for which the probability mass is computed.}
   \item{xcast}{vector. Vector of covariates used to compute \eqn{\lambda_{t+1}}.}
   \item{plot}{bool. If TRUE, a probability mass plot is plotted.}
   \item{title}{string. Title of the probability mass plot.}
   \item{xlab}{string. Label of the x-axis of the probability mass plot.}
   \item{ylab}{string. Label of the y-axis of the probability mass plot.}
   \item{width_bars}{positive float. Barwidth of the of the probability mass plot.}
   \item{decimals}{non-negative integer. Number of decimal floats within the axis of the probability mass plot.}
}
\details{
\emph{Usage description}

\emph{Brief math with formulas?}
}
\value{
A vector containing the probability mass. 
}
\references{
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
densities <- cocoForecast(fit)
densities$plot
}


\name{cocoSummary}
\alias{cocoSummary}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summary of Cocoreg Fit.
}
\description{
Summary function to display the regression results for a regression object created by \code{cocoreg}.
}
\usage{
cocoSummary(coco)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{coco}{object of class \code{coco.fit} or \code{coco.fit.c}. Usually created by \code{cocoreg}.}
}
\details{
The functions returns the parameter estimates and their standard errors, the value of the log-likelihod function evaluated at the parameters, and the logarithmic, quadratic as well as the ranked probability score.
}

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
size <- 1000
par <- c(0.5,0.2,0.05,0.3,0.3)
data.sim <- cocoSim(order = 2, type = "GP", par = par, size = size)
data <- data.sim$data
fit <- cocoReg(order = 2, type = "GP", data = data)
cocoSummary(fit)

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
cocoSummary(fit)
}

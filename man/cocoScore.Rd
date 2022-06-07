\name{cocoScore}
\alias{cocoScore}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Scoring Rules
}
\description{
This function provides the logarithmic, quadratic, and Ranked probability scores that can be used for model selection of fitted time series of counts models. 

}
\usage{
cocoScore(coco, val.num = 1e-05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{coco}{object of class "\code{coco.fit}"" or "\code{coco.fit.c}""}
  \item{val.num}{small positive number. Break value for numerical approximation of the scoring rules.}
}
\details{
Scoring rules are used as a model selection tool and are computed as averages over the set of in-sample predictions. They take the form \eqn{(T-p)^{-1} \sum_{t=p+1}^T s[F(x_t)]}, whereby \eqn{s[\cdot]} denotes a specific scoring rule and \eqn{x_t} the observed counts. The scoring rule should be minimized to improve the forecasting power of the model. Several scoring rules have been proposed by the literature and it is usually not clearly determined which of them should be applied and therefore it is recommended to use different scoring rules for model assessment.

Logarithmic score: \eqn{logs(F(x_t|\mathcal{F}_{t-1}))= -\log p(x_t|F_{t-1})}

Quadratic score: \eqn{qs(F(x_t|\mathcal{F}_{t-1}))= -2 p(x_t|F_{t-1}) + ||p||^2}

Ranked probability score: \eqn{rps(F(x_t|\mathcal{F}_{t-1}))= \sum_{j=0}^\infty \left[F(j) - 1(x_t \leq j)^2 \right] }
}

\value{
Returns a named list containing the following elements:
\item{log.score}{Logarithmic Score}
\item{quad.score}{Quadratic Score}
\item{rps.score}{Ranked probability Score}
}
\references{
Czado, C. and Gneitling, T. and Held, L. (2009) Predictive Model Assessment for Count Data. \emph{Biometrics}, \bold{65}, 4, 1254--1261.

Jung, R. C. and McCabe, B. P. and Tremayne, A. R. (2016) Model validation and diagnostics. \emph{Handbook of discrete-valued time series}, 189--218.

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
data.sim <- cocoSim(order = 2, type = "GP", par = par, size = size)
data <- data.sim$data
fit <- cocoReg(order = 2, type = "GP", data = data)
score <- cocoScore(fit)

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
score <- cocoScore(fit)
}


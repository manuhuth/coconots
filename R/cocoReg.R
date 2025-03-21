#' @title cocoReg
#' @description The function fits first- and second-order (Generalized) Poisson
#' integer autoregressive [(G)PAR] time series models for count data as discussed in Jung and Tremayne (2011). Autoregressive dependence on past counts is modeled using a special random operator that preserves integer values and, through closure under convolution, ensures that the marginal distribution remains within the same family as the innovations.
#'
#' These models can be viewed as stationary finite-order Markov chains, where the innovation distribution can be either Poisson or Generalized Poisson, the latter accounting for overdispersion. Estimation is performed via maximum likelihood, with an option to impose linear constraints. Without constraints, parameters may fall outside the theoretically feasible space, but optimization may be faster.
#'
#' Method of moments estimators are used to initialize numerical optimization, though custom starting values can be provided. If \proglang{julia} is installed, users can opt to run the optimization in \proglang{julia} for potentially faster computation and improved numerical stability via automatic differentiation. See below for details on the \proglang{julia} implementation.
#' @param type character, either "Poisson" or "GP" indicating the type of the innovation distribution
#' @param order integer, either 1 or 2 indicating the order of the model
#' @param data time series data to be used in the analysis
#' @param xreg optional matrix of explanatory variables (without constant term) for use in a regression model
#' @param constrained.optim logical indicating whether optimization should be constrained, currently only available in the R version
#' @param b.beta numeric value indicating the lower bound for the parameters of the explanatory variables for the optimization, currently only available in the \proglang{R} version
#' @param start optional numeric vector of starting values for the optimization
#' @param start.val.adjust logical indicating whether starting values should be adjusted, currently only available in the R version
#' @param method_optim character string indicating the optimization method to be used, currently only available in the R version. In the \proglang{julia} implementation this is by default the LBFGS algorithm
#' @param replace.start.val numeric value indicating the value to replace any invalid starting values, currently only available in the \proglang{R} version
#' @param iteration.start.val numeric value indicating the proportion of the interval to use as the new starting value, currently only available in the \proglang{R} version
#' @param method.hessian character string indicating the method to be used to approximate the Hessian matrix, currently only available in the \proglang{R} version
#' @param cores numeric indicating the number of cores to use, currently only available in the \proglang{R} version (default: 2)
#' @param julia if TRUE, the model is estimated with \proglang{julia}. This can improve computational speed significantly since \proglang{julia} makes use of derivatives using autodiff. In this case, only \code{type}, \code{order}, \code{data}, \code{xreg}, and \code{start} are used as other inputs (default: FALSE).
#' @param julia_installed if TRUE, the model \proglang{R} output will contain a \proglang{julia} compatible output element.
#' @param link_function Specifies the link function for the conditional mean of the innovation (\eqn{\lambda}). The default is `log`, but other available options include `identity` and `relu`. This parameter is applicable only when covariates are used. Note that using the `identity` link function may result in \eqn{\lambda} becoming negative. To prevent this, ensure all covariates are positive and restrict the parameter \eqn{\beta} to positive values by setting `b.beta` to a small positive value.
#' @importFrom JuliaConnectoR juliaGet
#' @author Manuel Huth
#' @return an object of class coco. It contains the parameter estimates, standard errors, the log-likelihood, 
#' and information on the model specifications. If \proglang{julia} is used for parameter estimation or the \proglang{julia} installation
#' parameter is set to TRUE, the results contain an additional Julia element that is called from the model \proglang{julia}
#' assessment tools if they are run with the \proglang{julia} implementation.
#' @details 
#' Let a time series of counts be \eqn{\{X_t\}} and be \eqn{R(\cdot)} a random operator that differs between model specifications.
#' For more details on the random operator, see Jung and Tremayne (2011) and Joe (1996).
#' The general first-order model is of the form
#' \deqn{X_t = R(X_{t-1}) + W_t,}
#' and the general second-order model of the form
#' \deqn{X_t = R(X_{t-1}, X_{t-2}) + W_t,}
#' where \eqn{W_t} are i.i.d Poisson (\eqn{W_t \sim Po(\lambda_t)}) or Generalized
#' Poisson (\eqn{W_t \sim GP(\lambda_t, \eta)}) innovations. Through closure under convolution
#' the marginal distributions of \eqn{\{X_t\}} are therefore Poisson or Generalized Poisson distributions, respectively.
#' 
#' If no covariates are used \eqn{\lambda_t = \lambda} and if covariates are used
#' \deqn{g(\lambda_t) = \left(\beta_0 + \sum_{j = 1}^k \beta_j \cdot z_{t,j} \right),}
#' whereby \eqn{z_{t,j}} is the \eqn{j}-th covariate at time \eqn{t} and \eqn{g} is a link function. 
#' Current supported link functions are the identity \eqn{g(x) = x} and a logarithmic link function 
#' \eqn{g(x) = \ln x}. To ensure positivity of \eqn{\lambda} if the identity function is used, \eqn{\beta_j, z_{t,j} > 0} must be enforced.
#' Alternatively, computational values of \eqn{\lambda \leq 0} can be set to a small positive value. 
#' This option is named 'relu', due to its similarity to a ReLu function commonly used in machine learning.
#' 
#' 
#' Standard errors are computed by the square root of the diagonal elements of the inverse Hessian.
#' 
#' This function is implemented in two versions. The default runs on RCPP. 
#' An alternative version uses a \proglang{julia} implementation which can be chosen by 
#' setting the argument \proglang{julia} to TRUE. In order to
#' use this feature, a running \proglang{julia} installation is required on the system.
#' The RCPP implementation uses the derivative-free Nelder-Mead optimizer
#' to obtain parameter estimates. The \proglang{julia} implementation makes use of \proglang{julia}'s
#' automatic differentiation in order to obtain gradients such that it can use the LBFGS algorithm for
#' optimization. This enhances the numeric stability of the optimization and yields 
#' an internal validation if both methods yield qualitatively same parameter estimates.
#' Furthermore, the \proglang{julia} implementation can increase the computational speed
#' significantly, especially for large models. 
#' 
#' The model assessment tools \code{\link{cocoBoot}}, \code{\link{cocoPit}}, and \code{\link{cocoScore}}
#' will use a \proglang{julia} implementation as well, if the \code{\link{cocoReg}} was run with \proglang{julia}.
#' Additionally, one can make the RCPP output of \code{\link{cocoReg}} compatible with the \proglang{julia}
#' model assessments by setting \code{julia_installed} to true. In this case, the user can choose
#' between the \pkg{RCPP} and the \proglang{julia} implementation for model assessment.
#' 
#' @references 
#' Jung, R. C. and Tremayne, A. R. (2011) Convolution-closed models for count time series with applications. \emph{Journal of Time Series Analysis}, \bold{32}, 268--280.
#' 
#' Joe, H. (1996) Time series models with univariate margins in the convolution-closed infinitely divisible class. \emph{Journal of Applied Probability}, \bold{33}, 664--677.
#' @examples
#' ## GP2 model without covariates
#' length <- 1000
#' par <- c(0.5,0.2,0.05,0.3,0.3)
#' data <- cocoSim(order = 2, type = "GP", par = par, length = length)
#' fit <- cocoReg(order = 2, type = "GP", data = data)
#' 
#' ##Poisson1 model with covariates
#' length <- 1000
#' period <- 12
#' sin <- sin(2*pi/period*(1:length))
#' cos <- cos(2*pi/period*(1:length))
#' cov <- cbind(sin, cos)
#' par <- c(0.2, 0.2, -0.2)
#' data <- cocoSim(order = 1, type = "Poisson", par = par, xreg = cov, length = length)
#' fit <- cocoReg(order = 1, type = "Poisson", data = data, xreg = cov)
#' @export

cocoReg <- function(type, order, data, xreg = NULL, 
                    constrained.optim = TRUE, b.beta = -10,
                    start = NULL, start.val.adjust = TRUE, method_optim = "Nelder-Mead",
                    replace.start.val = 1e-5, iteration.start.val = 0.6,
                    method.hessian = "Richardson", cores=2, julia=FALSE, 
                    julia_installed=FALSE, link_function="log") {
  seasonality <- c(1, 2) #will be used as argument in future versions
  
  if ((type != "GP") & (type != "Poisson")) {
    stop("Option 'type' must be either Poisson or GP")
  }
  
  if (julia){
    start_time <- Sys.time()
    fit_julia <- cocoRegJulia(type, order, data, xreg, start, link_function, b.beta)
    end_time <- Sys.time()
    fit_R <- JuliaConnectoR::juliaGet(fit_julia)
    julia_out <- transformJuliaRegOutputToR(xreg=xreg, pars=fit_R[["values"]][[8]],
                                            grad=NULL, hes=NULL,
                                            inv_hes=fit_R[["values"]][[7]],
                                            se=fit_R[["values"]][[11]],
                                            data=data,
                                            type=type,
                                            order=order,
                                            likelihood=fit_R[["values"]][[1]],
                                            end_time = end_time,
                                            start_time=start_time, 
                                            julia_reg=fit_julia)
    class(julia_out) <- "coco"
    julia_out$link_function <- link_function
    return(julia_out)
  }
  
  if (is.null(xreg)) {
    
    output <- cocoReg_base(
      type = type, order = order, data = data, seasonality = seasonality, 
      constrained.optim = constrained.optim, start = start,
      start.val.adjust = start.val.adjust, replace.start.val = replace.start.val, method_optim=method_optim,
      iteration.start.val = iteration.start.val, method.hessian = method.hessian, julia_installed=julia_installed, link_function=link_function
    )
  } else {
    output <- cocoReg_cov(
      type = type, order = order, data = data, xreg = xreg, seasonality = seasonality,
      constrained.optim = constrained.optim, b.beta = b.beta, start = start, method_optim=method_optim,
      start.val.adjust = start.val.adjust, replace.start.val = replace.start.val,
      iteration.start.val = iteration.start.val, method.hessian = method.hessian, julia_installed=julia_installed, link_function=link_function
    )
  }
  
  class(output) <- "coco"
  
  return(output)
}

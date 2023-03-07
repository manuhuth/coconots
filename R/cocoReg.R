#' @title cocoReg
#' @description The function fits first and second order (Generalized) Poisson Autoregressive (G)PAR
#' models presented in (Jung and Tremayne, 2010). Autoregressive dependence on past counts is modelled by a special random operator that not only preserve integer status, but also, via the property of closure under
#' convolution, ensure that the marginal distribution of the observed counts is from the same family as the innovations. 
#' The models can be thought of as stationary Markov chains of finite order, where the distribution of the innovations can either be Poisson or Generalized Poisson, where the latter can account for overdispersed data.
#' Maximum likelihood is used for estimation and the user can choose to include linear constraints or
#' not. If linear constraints are not included, it cannot be guaranteed that the parameters will lie in the
#' theoretically feasible parameter space, but the optimization process might be faster.
#' The function uses method of moments estimators to obtain starting values for the numerical optimization,
#' but the user can also specify their own starting values if desired. If Julia is installed, the user can choose whether the optimization is run in Julia
#' which might yield faster results, especially for large models with many observations, due to the use of automatic differentiation.
#' 
#' @param type character string indicating the type of model to be fitted
#' @param order integer vector indicating the order of the model
#' @param data time series data to be used in the analysis
#' @param xreg optional matrix of explanatory variables for use in a regression model
#' @param seasonality integer vector indicating the seasonal component of the model, currently only available in the R version
#' @param constrained.optim logical indicating whether optimization should be constrained, currently only available in the R version
#' @param b.beta numeric value indicating the lower bound for the parameters of the explanatory variables for the optimization, currently only available in the R version
#' @param start optional numeric vector of starting values for the optimization
#' @param start.val.adjust logical indicating whether starting values should be adjusted, currently only available in the R version
#' @param method_optim character string indicating the optimization method to be used, currently only available in the R version. In the julia implementation this is by default the LBFGS algorithm
#' @param replace.start.val numeric value indicating the value to replace any invalid starting values, currently only available in the R version
#' @param iteration.start.val numeric value indicating the proportion of the interval to use as the new starting value, currently only available in the R version
#' @param method.hessian character string indicating the method to be used to approximate the Hessian matrix, currently only available in the R version
#' @param cores numeric indicating the number of cores to use, currently only available in the R version
#' @param julia if TRUE, the model is estimated with Julia. This can improve the speed significantly since Julia makes use of derivatives using autodiff. In this case, only type, order, data, xreg, and start are used as other inputs.
#' @param julia_installed if TRUE, the model R output will contain a Julia compatible output element.
#' @author Manuel Huth
#' @return output of the regression analysis
#' 
#' @details 
#' Let a time series of counts be \eqn{\{X_t\}}, a random operator be \eqn{R(\cdot)}, the relevant past history of \eqn{X_t} be \eqn{\mathcal{F}_{t-1}} and \eqn{W_t} be i.i.d discrete innovations. The general model is of the form
#' \deqn{X_t = R(\mathcal{F}_{t-1}) + W_t.}
#' For first-order and second-order models the relevant past history depends on the model's order and is defined by \eqn{\mathcal{F}_{t-1} = X_{t-s_1}} and \eqn{\mathcal{F}_{t-1} = (X_{t-s_1}, X_{t-s_2})}, respectively. \eqn{s_1} < \eqn{s_2} are the parameters indicating the degree of stochastic seasonality.
#' The innovations follow either a Poisson or a generalized Poisson distribution. Such that either \eqn{W_t \sim Pois(\lambda_t)} or \eqn{W_t \sim GP(\lambda_t, \eta)}.
#' If no covariates are used \eqn{\lambda_t = \lambda} and if covariates are used \deqn{\lambda_t = \exp{\left(\beta_0 + \sum_{j = 1}^p \beta_j \cdot y_{t,j} \right)},} whereby \eqn{y_{t,j}} is the \eqn{j}-th covariate at time \eqn{t}.
#' Standard errors are computed by the square root of the diagonal elements of the inverse hessian.
#' Add more on julia usage
#' @references 
#' Jung, R. C. and Tremayne, A. R. (2010) Convolution-closed models for count timeseries with applications. \emph{Journal of Time Series Analysis}, \bold{32}, 3, 268--280.
#' 
#' Joe, H. (1996) Time series models with univariate margins in the convolution-closed infinitely divisible class. \emph{Journal of Applied Probability}, 664--677.
#' @examples
#' ## GP2 model without covariates
#' length <- 1000
#' par <- c(0.5,0.2,0.05,0.3,0.3)
#' data.sim <- cocoSim(order = 2, type = "GP", par = par, length = length)
#' data <- data.sim$data
#' fit <- cocoReg(order = 2, type = "GP", data = data)
#' 
#' ##Poisson1 model with covariates
#' length <- 1000
#' period <- 50
#' sin <- sin(2*pi/period*(1:length))
#' cos <- cos(2*pi/period*(1:length))
#' cov <- cbind(sin, cos)
#' par <- c(0.2, 0.2, -0.2)
#' data.sim <- cocoSim(order = 1, type = "Poisson", par = par, xreg = cov, length = length)
#' data <- data.sim$data
#' fit <- cocoReg(order = 1, type = "Poisson", data = data, xreg = cov)
#' @export

cocoReg <- function(type, order, data, xreg = NULL, seasonality = c(1, 2),
                    constrained.optim = TRUE, b.beta = -10,
                    start = NULL, start.val.adjust = TRUE, method_optim = "Nelder-Mead",
                    replace.start.val = 1e-5, iteration.start.val = 0.99,
                    method.hessian = "Richardson", cores=2, julia=FALSE, 
                    julia_installed=FALSE) {
  
  if (julia){
    start_time <- Sys.time()
    fit_julia <- cocoRegJulia(type, order, data, xreg, start)
    end_time <- Sys.time()
    fit_R <- JuliaConnectoR::juliaGet(fit_julia)
    return(transformJuliaRegOutputToR(xreg=xreg, pars=fit_R[["values"]][[8]],
                                      grad=NULL, hes=NULL,
                               inv_hes=fit_R[["values"]][[7]],
                               se=fit_R[["values"]][[11]],
                               data=data,
                               type=type,
                               order=order,
                               likelihood=-fit_R[["values"]][[1]],
                               end_time = end_time,
                               start_time=start_time, 
                               julia_reg=fit_julia))
  }
  
  if (is.null(xreg)) {
    
    output <- cocoReg_base(
      type = type, order = order, data = data, seasonality = seasonality, 
      constrained.optim = constrained.optim, start = start,
      start.val.adjust = start.val.adjust, replace.start.val = replace.start.val, method_optim=method_optim,
      iteration.start.val = iteration.start.val, method.hessian = method.hessian, julia_installed=julia_installed
    )
  } else {
    output <- cocoReg_cov(
      type = type, order = order, data = data, xreg = xreg, seasonality = seasonality,
      constrained.optim = constrained.optim, b.beta = -10, start = start, method_optim=method_optim,
      start.val.adjust = start.val.adjust, replace.start.val = replace.start.val,
      iteration.start.val = iteration.start.val, method.hessian = method.hessian, julia_installed=julia_installed
    )
  }
  
  class(output) <- "coco"

  return(output)
}

#' @title Fitting First- and Second Order (G)PAR Models 
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
    fit_julia <<- cocoRegJulia(type, order, data, xreg, start, link_function, b.beta)
    end_time <- Sys.time()
    
    juliaLet("global fit2 = x", x=fit_julia)
    fit_no_optimization_in_dict <- juliaEval('
               delete!(fit2, "optimization")
               ')
    
    fit_R <- JuliaConnectoR::juliaGet(fit_no_optimization_in_dict)

    julia_out <- transformJuliaRegOutputToR(xreg=xreg, pars=fit_R[["values"]][[8]],
                                            grad=NULL, hes=NULL,
                                            inv_hes=fit_R[["values"]][[7]],
                                            se=fit_R[["values"]][[10]],
                                            data=data,
                                            type=type,
                                            order=order,
                                            likelihood=fit_R[["values"]][[1]],
                                            end_time = end_time,
                                            start_time=start_time, 
                                            julia_reg=fit_no_optimization_in_dict)

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

#' @export
summary.coco <- function(object, ..., score = FALSE) {
  # Attach score flag to the object and set the class to "summary.coco"
  object$score <- score
  class(object) <- "summary.coco"
  return(object)
}

#' @export
print.summary.coco <- function(x, ...) {
  coco <- x
  
  # Create data vectors with rounded values
  estimates <- round(coco$par, 4)
  se_values <- round(coco$se, 4)
  t_values <- round(coco$par / coco$se, 4)
  
  # Combine into a data frame for display
  df <- data.frame("Estimate" = estimates,
                   "Std. Error" = se_values,
                   "t" = t_values)
  
  # Generate parameter names based on model order and type
  if (coco$order == 1) {
    if (coco$type == "Poisson") {
      param_names <- c("lambda", "alpha")
    } else if (coco$type == "GP") {
      param_names <- c("lambda", "alpha", "eta")
    } else {
      param_names <- character(0)
    }
  } else if (coco$order == 2) {
    if (coco$type == "Poisson") {
      param_names <- c("lambda", "alpha1", "alpha2", "alpha3")
    } else if (coco$type == "GP") {
      param_names <- c("lambda", "alpha1", "alpha2", "alpha3", "eta")
    } else {
      param_names <- character(0)
    }
  } else {
    param_names <- character(0)
  }
  
  # Append covariate names if they exist
  if (!is.null(coco$cov)) {
    if (is.null(colnames(coco$cov))) {
      colnames(coco$cov) <- paste0("X", seq_len(ncol(coco$cov)))
    }
    # Remove "lambda" from the initial parameters and add covariate names
    param_names <- c(param_names[param_names != "lambda"], colnames(coco$cov))
  }
  
  # Set the row names of the coefficients dataframe
  rownames(df) <- param_names
  
  cat("Coefficients:\n")
  print(df, print.gap = 3, quote = FALSE, na.print = "")
  
  # Check if Julia regression flag is set
  julia <- !is.null(coco$julia_reg)
  
  # Build the output string with basic model statistics
  output_lines <- paste0(
    "\nType: ", coco$type,
    "\nOrder: ", coco$order,
    "\n\nLog-likelihood: ", round(coco$likelihood, 4)
  )
  
  if (isTRUE(coco$score)) {
    # When score flag is true, compute the additional scores
    score_vals <- cocoScore(coco, julia = julia)
    output_lines <- paste0(
      output_lines,
      "\nLogarithmic score: ", round(score_vals$log.score, 4),
      "\nQuadratic score: ", round(score_vals$quad.score, 4),
      "\nRanked probability score: ", round(score_vals$rps.score, 4),
      "\nAIC: ", round(score_vals$aic, 4),
      "\nBIC: ", round(score_vals$bic, 4)
    )
  } else {
    # Calculate AIC and BIC using the default formulas
    aic_value <- round(2 * length(coco$par) - 2 * coco$likelihood, 4)
    bic_value <- round(length(coco$par) * log(length(coco$ts)) - 2 * coco$likelihood, 4)
    output_lines <- paste0(
      output_lines,
      "\nAIC: ", aic_value,
      "\nBIC: ", bic_value
    )
  }
  
  cat(output_lines, "\n")
}

#' @export
print.coco <- function(x, ...) {
  # Title and basic model information
  cat("Fitted coco Model\n")
  cat("====================================\n")
  cat("Type: ", x$type, "\n", sep = "")
  cat("Order: ", x$order, "\n", sep = "")
  cat("Link Function: ", x$link_function, "\n", sep = "")
  
  # Coefficient information (with names assigned via coef.coco)
  cat("\nCoefficients:\n")
  print(coef(x))
  
  # Log-likelihood
  cat("\nLog-likelihood: ", round(x$likelihood, 4), "\n", sep = "")
  
  # Duration of model estimation (if available)
  if(!is.null(x$duration)) {
    cat("Duration (seconds): ", as.numeric(x$duration), "\n", sep = "")
  }
  
  # Covariate information (if available)
  if(!is.null(x$cov)) {
    cov_names <- colnames(x$cov)
    if(is.null(cov_names)) {
      cov_names <- paste0("X", seq_len(ncol(x$cov)))
    }
    cat("Covariates: ", paste(cov_names, collapse = ", "), "\n")
  }
  
  invisible(x)
}

#' @export
fitted.coco <- function(object, ...) {
  # Compute the residuals (fitted values, residuals, and other diagnostics)
  res <- cocoResid(object, ...)
  # Return only the fitted values
  res$fitted
}

#' @export
residuals.coco <- function(object, ...) {
  res <- cocoResid(object, ...)
  res$residuals
}

#' @export
coef.coco <- function(object, ...) {
  params <- object$par
  names(params) <- get_coco_param_names(object)
  params
}


#' @export
vcov.coco <- function(object, ...) {
  if (!is.null(object[["inv hessian"]])) {
    V <- object[["inv hessian"]]
  } else {
    stop("Covariance matrix not available for this object")
  }
  
  param_names <- get_coco_param_names(object)
  rownames(V) <- param_names
  colnames(V) <- param_names
  V
}

#' @export
rstandard.coco <- function(object, ...) {
  res <- cocoResid(object, ...)
  res$pe.resid
}

#' @export
nobs.coco <- function(object, ...) {
  length(object$ts)
}

#' @export
extractAIC.coco <- function(fit, scale = 0, k = 2, ...) {
  df <- length(fit$par)
  aic_value <- -2 * fit$likelihood + k * df
  c(df, aic_value)
}

#' @export
logLik.coco <- function(object, ...) {
  object$likelihood
}

get_coco_param_names <- function(object) {
  # Define parameter names based on the model order and type.
  if (object$order == 1) {
    if (object$type == "Poisson") {
      param_names <- c("lambda", "alpha")
    } else if (object$type == "GP") {
      param_names <- c("lambda", "alpha", "eta")
    } else {
      param_names <- character(0)
    }
  } else if (object$order == 2) {
    if (object$type == "Poisson") {
      param_names <- c("lambda", "alpha1", "alpha2", "alpha3")
    } else if (object$type == "GP") {
      param_names <- c("lambda", "alpha1", "alpha2", "alpha3", "eta")
    } else {
      param_names <- character(0)
    }
  } else {
    param_names <- character(0)
  }
  
  # If covariates are included, remove "lambda" and append the covariate names.
  if (!is.null(object$cov)) {
    if (is.null(colnames(object$cov))) {
      colnames(object$cov) <- paste0("X", seq_len(ncol(object$cov)))
    }
    param_names <- c(param_names[param_names != "lambda"], colnames(object$cov))
  }
  
  param_names
}

#' @export
autoplot.coco <- function(object, which = "fitted", ...) {
  # Ensure ggplot2 namespace is loaded
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for this function to work")
  }
  
  # Calculate fitted values and residuals using your internal function.
  resid_out <- cocoResid(object, ...)
  ts_obs <- object$ts
  # Generate a time index; note that fitted values and residuals are computed for observations
  # starting at index (order + 1)
  time_index <- seq_along(ts_obs)
  valid_index <- time_index[(object$order + 1):length(ts_obs)]
  
  if (which == "fitted") {
    # Create a data frame for plotting the observed and fitted values.
    df <- data.frame(Time = valid_index,
                     Observed = ts_obs[(object$order + 1):length(ts_obs)],
                     Fitted = resid_out$fitted)
    
    p <- ggplot2::ggplot(df, ggplot2::aes(x = Time)) +
      ggplot2::geom_line(ggplot2::aes(y = Observed), color = "black") +
      ggplot2::geom_line(ggplot2::aes(y = Fitted), color = "steelblue", linetype = "dashed") +
      ggplot2::labs(title = "Observed vs. Fitted Values",
                    y = "Count") +
      ggplot2::theme_bw()
    
    return(p)
    
  } else if (which == "residuals") {
    # Plot the standardized (Pearson) residuals over time.
    df <- data.frame(Time = valid_index,
                     Residuals = resid_out$pe.resid)
    
    p <- ggplot2::ggplot(df, ggplot2::aes(x = Time, y = Residuals)) +
      ggplot2::geom_line(color = "black") +
      ggplot2::labs(title = "Pearson Residuals",
                    y = "Pearson Residual") +
      ggplot2::theme_bw()
    
    return(p)
    
  } else if (which == "acf") {
    # Create an ACF plot of the standardized residuals.
    acf_obj <- stats::acf(resid_out$pe.resid, plot = FALSE)
    df <- data.frame(Lag = as.numeric(acf_obj$lag),
                     ACF = acf_obj$acf)
    
    p <- ggplot2::ggplot(df, ggplot2::aes(x = Lag, y = ACF)) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
      ggplot2::geom_bar(stat = "identity", fill = "steelblue") +
      ggplot2::labs(title = "ACF of Pearson Residuals") +
      ggplot2::theme_bw()
    
    return(p)
    
  } else {
    stop("Invalid 'which' argument. Choose from 'fitted', 'residuals', or 'acf'.")
  }
}

#' @export
plot.coco <- function(x, interactive = TRUE, ...) {
  
  if (interactive) {
    # Display the standardized residuals plot first.
    p_resid <- autoplot.coco(x, which = "residuals", ...)
    print(p_resid)
    
    # Wait for user input before continuing.
    readline(prompt = "Press Enter to view the ACF plot...")
    
    # Now display the ACF plot of the standardized residuals.
    p_acf <- autoplot.coco(x, which = "acf", ...)
    print(p_acf)
    
  } else {
    # For non-interactive usage, default to a single plot.
    # For example, the observed vs. fitted plot.
    p <- autoplot.coco(x, which = "acf", ...)
    print(p)
  }
}
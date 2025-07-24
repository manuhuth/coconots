#' @title Residual Based Model Assessment Procedure
#' @description Calculates the (Pearson) residuals of a fitted model for model evaluation purposes.
#' @param coco An object of class "coco
#' @param val.num A non-negative real number that halts the calculation once the cumulative probability reaches 1-\code{val.num}
#' @author Manuel Huth
#' @importFrom stats var fitted residuals
#' @importFrom ggplot2 autoplot
#' @return a list that includes the (Pearson) residuals, conditional expectations, conditional variances,
#' and information on the model specifications.
#'@details The Pearson residuals are computed as the scaled
#'deviation of the observed count from its conditional expectation given the relevant 
#'past history, including covariates, if applicable. If a fitted model is correctly specified,
#'the Pearson residuals should exhibit mean zero, variance one, and no significant serial correlation. 
#'@export
cocoResid <- function(coco, val.num = 1e-10) {
  start.time <- Sys.time()

  if (val.num <= 0) {
    stop("Option val.num must be a non-negative real number")
  }

  seasonality <- c(1, 2) #will be used as argument in future versions coco$seasonality
  bool_covariates <- !is.null(coco$cov)
  par <- coco$par
  data <- coco$ts
  
  #handle covariates
  if (bool_covariates){
    # set up values for lambda
    xreg <- coco$cov
    vec_lambda <- par[(length(par)-ncol(xreg)+1):length(par)]
    lambdas <- c()
    for (j in 1:length(data)) {
      lambdas[j] <- applyLinkFunction(as.numeric(as.vector(xreg[j, ])) %*% vec_lambda, coco$link_function)
    }
  } else {
    lambda <- par[1]
  }
    

  ### GP/PAR1
  if (coco$order == 1) {

    alpha <- par[2-bool_covariates]
      
    if (coco$type == "Poisson"){
      eta <- 0
    } else {
      eta <- par[3-bool_covariates]
    }
      
    # Fitted Values
    fitted <- c(rep(mean(data), seasonality[1])) 
    varX <- c(rep(stats::var(data), seasonality[1])) 
      
    fitted_varX <- lapply((seasonality[1] + 1):(length(data)), function(j) {
      if (bool_covariates){
        lambda <- lambdas[j]
      }
      fitted_value <- alpha * data[j - seasonality[1]] + lambda / (1-eta)
      varX_value <- alpha * (1 - alpha) * data[j - seasonality[1]] + lambda / (1-eta)^3
      return(list(fitted = fitted_value, varX = varX_value))
    })
      
    fitted <- unlist(lapply(fitted_varX, function(x) x$fitted))
    varX <- unlist(lapply(fitted_varX, function(x) x$varX))
      
  } # end if GP/PAR1


  ### GP/PAR2
  if (coco$order == 2) {
      
    alpha1 <- par[2-bool_covariates]
    alpha2 <- par[3-bool_covariates]
    alpha3 <- par[4-bool_covariates]
      
    U <- 1 / (1 - alpha1 - alpha2 - alpha3)
      
    if (coco$type == "Poisson"){
      eta <- 0
    } else {
      eta <- par[5-bool_covariates]
    }
      
    # Fitted Values
    fitted <- c(rep(mean(data), seasonality[2])) 
    varX <- c(rep(stats::var(data), seasonality[2])) 
      
    fitted_varX <- lapply((seasonality[2] + 1):(length(data)), function(j) {
      if (bool_covariates){
        lambda <- lambdas[j]
      }
        
      y <- data[j - seasonality[1]]
      z <- data[j - seasonality[2]]
        
      # compute conditional mean of R
      meanR <- 0
      varR <- 0
      prob_accumulated <- 0
      probabilities <- c()
      r <- 0
      while (prob_accumulated < (1 - val.num)) {
        prob <- dR2(r, y, z, lambda, alpha1, alpha2, alpha3, eta)
        probabilities <- c(probabilities, prob)
        prob_accumulated <- prob_accumulated + prob
        r <- r + 1
      }
      
      #define objects for cond. mean and cond. variance computation
      n_relevant_probabilities <- r-1 #r-1 as +1 is added in the end within the while loop
      relevant_counts <- 0:(r-1)
      
      #compute conditional mean and conditional variance
      meanR <- sum(probabilities * relevant_counts)
      varR <-  sum(probabilities * ((relevant_counts - meanR)^2))
        
      fitted_value <- meanR + lambda / (1 - eta)
      varX_value <- varR + lambda / (1 - eta)^3

      return(list(fitted = fitted_value, varX = varX_value))
    })
      
    fitted <- unlist(lapply(fitted_varX, function(x) x$fitted))
    varX <- unlist(lapply(fitted_varX, function(x) x$varX))
      
  } # end if GP/PAR2


  # Residuals
  residuals <- data[(coco$order+1):length(data)] - fitted

  # Pearson Residuals
  peResid <- residuals / varX^0.5

  end.time <- Sys.time()

  time <- end.time - start.time

  list_out <- list(
    "fitted" = fitted, "residuals" = residuals, "pe.resid" = peResid,
    "cond.var" = varX, "type" = coco$type, "order" = coco$order, "ts" = coco$ts,
    "par" = par,"duration" = time, xreg = coco$cov
  )

  class(list_out) <- "cocoResid"

  return(list_out)
} # end function


#------------------S3 methods--------------------------------------------------
#' @importFrom ggplot2 autoplot theme_bw ggtitle
#' @export
ggplot2::autoplot
#' @exportS3Method
autoplot.cocoResid <- function(object, ...){
  forecast::ggAcf(object$pe.resid) + theme_bw() + ggtitle("Pearson Residuals") + ggplot2::theme(text = ggplot2::element_text(size = 20)) 
}

#' @exportS3Method
plot.cocoResid <- function(x, ...) {
  p <- autoplot(
    x,
    ...
  )
  print(p)
}

#'@exportS3Method
print.cocoResid <- function(x, ...) {
  print(autoplot(x, ...,))
  invisible(x)
}

#' @exportS3Method
print.cocoResid <- function(x, ...) {
  cat("cocoResid Object\n")
  cat("Model Type: ", x$type, "\n")
  cat("Model Order: ", x$order, "\n")
  cat("Number of Fitted Values: ", length(x$fitted), "\n")
  
  # Display covariate information
  if (is.null(x$xreg)) {
    cat("Covariates: None\n")
  } else {
    n_cov <- if (is.data.frame(x$xreg) || is.matrix(x$xreg)) ncol(as.matrix(x$xreg)) else 1
    cat("Covariates: Present (", n_cov, " variable(s))\n", sep = "")
  }
  
  cat("Calculation Duration: ", as.numeric(x$duration, units = "secs"), " seconds\n", sep = "")
  cat("\n")
  
  invisible(x)
}

#' @exportS3Method
summary.cocoResid <- function(object, ...) {
  # Build a summary list extracting key information from the cocoResid object
  res <- list(
    model_type = object$type,
    model_order = object$order,
    ts_length = length(object$ts),
    n_fitted = length(object$fitted),
    covariates = if (is.null(object$xreg)) {
      "None"
    } else if (is.data.frame(object$xreg) || is.matrix(object$xreg)) {
      paste("Present (", ncol(as.matrix(object$xreg)), " variable(s))", sep = "")
    } else {
      "Present (1 variable)"
    },
    residuals_summary = summary(object$residuals),
    pearson_residuals_summary = summary(object$pe.resid),
    cond_var_summary = summary(object$cond.var),
    duration = as.numeric(object$duration, units = "secs")
  )
  
  # Assign a class to the summary object so that print.summary.cocoResid is dispatched.
  class(res) <- "summary.cocoResid"
  return(res)
}

#' @exportS3Method
print.summary.cocoResid <- function(x, ...) {
  cat("---- Residual Analysis Summary ----\n")
  cat("Model Type: ", x$model_type, "\n")
  cat("Model Order: ", x$model_order, "\n")
  cat("Time Series Length: ", x$ts_length, "\n")
  cat("Number of Fitted Values: ", x$n_fitted, "\n")
  cat("Covariates: ", x$covariates, "\n\n")
  
  cat("--- Residuals ---\n")
  print(x$residuals_summary)
  
  cat("\n--- Pearson Residuals ---\n")
  print(x$pearson_residuals_summary)
  
  cat("\n--- Conditional Variance ---\n")
  print(x$cond_var_summary)
  
  cat("\nCalculation Duration: ", x$duration, " seconds\n", sep = "")
  cat("\n")
  invisible(x)
}


#' @exportS3Method
fitted.cocoResid <- function(object, ...) {
  object$fitted
}

#' @exportS3Method
residuals.cocoResid <- function(object, ...) {
  object$residuals
}


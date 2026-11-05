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

  coco_check_count_limit(coco$ts)
  coco_check_length(length(coco$ts), coco$order + 1, "cocoResid",
                    sprintf("an order-%d model has residuals for the observations after the first %d",
                            coco$order, coco$order))

  seasonality <- c(1, 2) #will be used as argument in future versions coco$seasonality
  bool_covariates <- !is.null(coco$cov)
  par <- coco$par
  data <- coco$ts
  
  #handle covariates
  if (bool_covariates){
    # set up values for lambda
    xreg <- coco$cov
    vec_lambda <- par[(length(par)-ncol(xreg)+1):length(par)]
    lambdas <- applyLinkFunctionVec(xreg %*% vec_lambda, coco$link_function)
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
      
    # Vectorised form of the former per-j loop
    #   fitted[j] <- alpha * data[j - 1] + lambda / (1-eta)
    #   varX[j]   <- alpha * (1 - alpha) * data[j - 1] + lambda / (1-eta)^3
    # with the same operations per element. The scalar operands are
    # replicated with rep() (which keeps their names), so every element also
    # gets the names the scalar expression gave it (names of the first
    # operand if it has any, else of the second).
    idx <- (seasonality[1] + 1):(length(data))
    n_idx <- length(idx)
    data_lag <- data[idx - seasonality[1]]
    lambda_t <- if (bool_covariates) lambdas[idx] else rep(lambda, n_idx)
    fitted <- rep(alpha, n_idx) * data_lag + lambda_t / rep(1 - eta, n_idx)
    varX <- rep(alpha * (1 - alpha), n_idx) * data_lag + lambda_t / rep((1 - eta)^3, n_idx)
      
  } # end if GP/PAR1


  ### GP/PAR2
  if (coco$order == 2) {
      
    alpha1 <- par[2-bool_covariates]
    alpha2 <- par[3-bool_covariates]
    alpha3 <- par[4-bool_covariates]
      
    if (coco$type == "Poisson"){
      eta <- 0
    } else {
      eta <- par[5-bool_covariates]
    }
      
    idx <- (seasonality[2] + 1):(length(data))
    n_idx <- length(idx)
    y_t <- data[idx - seasonality[1]]
    z_t <- data[idx - seasonality[2]]
    lambda_t <- if (bool_covariates) lambdas[idx] else rep(lambda, n_idx)

    # conditional mean and variance of the thinning part R given (y, z)
    momentsR <- function(y, z, lambda) {
      # P(R = r | y, z) for r = 0, 1, ... until 1 - val.num is covered,
      # computed in a single C++ call
      probabilities <- dR2_vec(y, z, lambda, alpha1, alpha2, alpha3, eta, val.num)
      coco_check_finite_pmf(probabilities, "conditional pmf of the thinning part R given (y, z)")
      relevant_counts <- 0:(length(probabilities) - 1)

      #compute conditional mean and conditional variance
      meanR <- sum(probabilities * relevant_counts)
      varR <-  sum(probabilities * ((relevant_counts - meanR)^2))
      c(meanR, varR)
    }

    if (bool_covariates) {
      mom <- vapply(seq_len(n_idx), function(i) momentsR(y_t[i], z_t[i], lambda_t[i]),
                    numeric(2))
    } else {
      # without covariates the moments depend on t only through (y, z), which
      # dR2_vec() receives as integers: one evaluation per distinct history
      key <- paste(as.integer(y_t), as.integer(z_t))
      first <- which(!duplicated(key))
      mom_u <- vapply(first, function(i) momentsR(y_t[i], z_t[i], lambda),
                      numeric(2))
      mom <- mom_u[, match(key, key[first]), drop = FALSE]
    }

    # former per-j values meanR + lambda / (1 - eta) and
    # varR + lambda / (1 - eta)^3 (with their names, as in the order-1 case)
    fitted <- mom[1, ] + lambda_t / rep(1 - eta, n_idx)
    varX <- mom[2, ] + lambda_t / rep((1 - eta)^3, n_idx)
      
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


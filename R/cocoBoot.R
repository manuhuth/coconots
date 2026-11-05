#' @title Bootstrap Based Model Assessment Procedure
#' @description Model checking procedure emphasising reproducibility in fitted models, as proposed by Tsay (1992).
#' @param coco An object of class coco
#' @param numb.lags Number of lags for which to compute sample autocorrelations (default: 21).
#' @param rep.Bootstrap Number of bootstrap replicates to use (default: 1000)
#' @param conf.alpha \eqn{100(1-\code{conf.alpha})\%} probability interval for the acceptance envelopes (default: 0.05)
#' @param julia  if TRUE, the bootstrap is run with \proglang{julia} (default: FALSE)
#' @param julia_seed Seed for the \proglang{julia} implementation. Only used if \proglang{julia} equals TRUE
#' @return an object of class cocoBoot. It contains the bootstrapped confidence intervals
#' of the autocorrelations and information on the model specifications.
#' @details Bootstrap-generated acceptance envelopes for the autocorrelation function provides an overall evaluation by comparing it with the sample autocorrelation function in a joint plot. 
#' @importFrom forecast Acf
#' @importFrom matrixStats rowQuantiles
#' @importFrom ggplot2 autoplot
#' @references 
#' Tsay, R. S. (1992) Model checking via parametric bootstraps in time series analysis. \emph{Applied Statistics} \bold{41}, 1--15.
#' @examples
#' lambda <- 1
#' alpha <- 0.4
#' set.seed(12345)
#' data <- cocoSim(order = 1, type = "Poisson", par = c(lambda, alpha), length = 100)
#' fit <- cocoReg(order = 1, type = "Poisson", data = data)
#'
#' # bootstrap model assessment - R implementation
#' boot_r <- cocoBoot(fit, rep.Bootstrap=400)
#' plot(boot_r)
#' @export

cocoBoot <- function(coco, numb.lags = 21, rep.Bootstrap = 1000,
                 conf.alpha = 0.05, julia = FALSE, julia_seed = NULL
                 ) {
  start.time <- Sys.time()
  confidence <- 1 - conf.alpha
  if ((confidence <= 0) | (confidence >= 1)) {
    stop("Option confidence must be a real number between 0 and 1")
  }

  if ((numb.lags != round(numb.lags)) | (numb.lags < 1)) {
    stop("The value of numb.lags must be a positive integer")
  }

  if ((rep.Bootstrap != round(rep.Bootstrap)) | (rep.Bootstrap < 1)) {
    stop("The value of rep.Bootstrap must be a positive integer")
  }
  
  data <- coco$ts
  
  if (!is.null(coco$julia_reg) & julia){
    if (!is.null(julia_seed)){
      setJuliaSeed(julia_seed)
    }
    addJuliaFunctions()
    # As on the RCPP path below: with covariates the replicates are simulated
    # without burn-in (cocoSim's documented default), without covariates with one.
    n_burn_in <- if (is.null(coco$cov)) 200L else 0L
    coco_boot <- JuliaConnectoR::juliaGet( JuliaConnectoR::juliaCall("cocoBoot", coco$julia_reg, 1:numb.lags, rep.Bootstrap, 1-confidence, n_burn_in))
    #acfdata <- coco_boot$values[[2]]
    ac <- t(getJuliaValue(coco_boot, "pacfs"))
    #confidence_bands <- data.frame(cbind(coco_boot$values[[3]], coco_boot$values[[4]]))
    #colnames(confidence_bands) <- c("upper", "lower")
  } else {
    
    seasonality <- c(1,2) #will be used as argument in future versions

    conf.alpha <- 1 - confidence

    par <- coco$par

    # Parametric Bootstrap
    T <- length(data)
    nlags <- numb.lags
    nB <- rep.Bootstrap
    ac <- matrix(NaN, nrow = nlags, ncol = nB)

    # All replicates are simulated and their sample ACFs computed in one C++
    # call (bootAcfBase / bootAcfCov). Each replicate is the series the former
    # per-replicate cocoSim() call returned (same model, burn-in, covariate
    # alignment and draw order of the random numbers) and its ACF is the
    # stats::acf estimate forecast::Acf returned, so the procedure is
    # unchanged; the thinning distributions are memoized across replicates.
    if ((is.null(coco$cov)) ) {
      # cocoSim(length = T + 10) with its burn-in of 200, first 10 values dropped
      bootCheckParBase(coco$type, coco$order, par)
      gp <- coco$type == "GP"
      if (coco$order == 1) {
        ac <- bootAcfBase(1L, gp, par[1], par[2], 0, 0, if (gp) par[3] else 0,
                          as.integer(T + 10 + 200), as.integer(T),
                          as.integer(nlags), as.integer(nB))
      } else {
        ac <- bootAcfBase(2L, gp, par[1], par[2], par[3], par[4], if (gp) par[5] else 0,
                          as.integer(T + 10 + 200), as.integer(T),
                          as.integer(nlags), as.integer(nB))
      }
    }


    if (!is.null(coco$cov)) {
      # cocoSim(length = T, xreg = xreg) without burn-in
      xreg <- as.matrix(coco$cov)
      link_function <- coco$link_function
      type <- coco$type
      order <- coco$order
      bootCheckParCov(type, order, par, xreg)
      gp <- type == "GP"
      n_ar <- if (order == 1) 1 else 3
      alpha <- par[1:n_ar]
      eta <- if (gp) par[n_ar + 1] else 0
      vec_lambda <- par[-(1:(n_ar + gp))]
      init_lambdas <- numeric(order)
      for (t in 1:order) {
        if (order == 2 && nrow(xreg) <= 1) {
          lambda_start1 <- c(0, 0)
        } else {
          lambda_start1 <- as.numeric(as.vector(xreg[t, ])) %*% vec_lambda
        }
        init_lambdas[t] <- applyLinkFunction(lambda_start1, link_function)
      }
      lambdas <- applyLinkFunctionVec(xreg %*% vec_lambda, link_function)
      n_unif <- T
      if (order == 2 && nrow(xreg) < (if (gp) 3 else 2)) n_unif <- T + 2
      if (order == 1) alpha <- c(alpha, 0, 0)
      ac <- bootAcfCov(as.integer(order), gp, alpha[1], alpha[2], alpha[3], eta,
                       vec_lambda, xreg, lambdas, init_lambdas, as.integer(T),
                       as.integer(n_unif), as.integer(nlags), as.integer(nB),
                       link_function)
    }


  } #end julia
  acfdata <- forecast::Acf(data, plot = FALSE, lag.max = numb.lags)$acf[2:(numb.lags + 1)]
  confidence_bands <- data.frame(matrixStats::rowQuantiles(ac, probs = c((1-confidence)/2, 1-(1-confidence)/2)))
  df_plot <- cbind(acfdata, 1:numb.lags, confidence_bands)
  colnames(df_plot) <- c("y", "x", "lower", "upper")
  

  end.time <- Sys.time()
  time <- end.time - start.time
  list_out <- list(
    "type" = coco$type, "order" = coco$order, "ts" = coco$ts, "cov" = coco$cov, means = acfdata,
    "confidence" = confidence_bands, "duration" = time, "df_plot" = df_plot, "rep.Bootstrap" = rep.Bootstrap
  )
  class(list_out) <- "cocoBoot"
  return(list_out)
} # end function

#' @exportS3Method
print.cocoBoot <- function(x, ...) {
  cat("cocoBoot Object")
  cat("\nModel Type: ", x$type, sep = "")
  cat("\nModel Order: ", x$order, sep = "")
  cat("\nBootstrap Draws: ", as.numeric(x$rep.Bootstrap), sep = "")
  cat("\nBootstrapping Duration: ", as.numeric(x$duration, units = "secs"), " seconds\n", sep = "")
  cat("\n")
  
  invisible(x)
}

#' @importFrom ggplot2 autoplot
#' @export
ggplot2::autoplot
#' @exportS3Method
autoplot.cocoBoot <- function(object, ...){
  
  pl <- ggplot2::ggplot(object$df_plot, ggplot2::aes_string("x", "y")) +
    ggplot2::geom_point()+
    ggplot2::geom_line()+
    ggplot2::geom_ribbon(data=object$df_plot,ggplot2::aes_string(ymin="lower",ymax="upper"),
                         fill="steelblue", alpha=0.3) +
    ggplot2::theme_bw() + ggplot2::xlab("Lags") + ggplot2::ylab("Autocorrelation") +
    ggplot2::ggtitle("Parametric Bootstrap") +
    ggplot2::theme(text = ggplot2::element_text(size = 20))
  pl
  
}

#' @exportS3Method
plot.cocoBoot <- function(x, ...) {
  p <- autoplot(
    x,
    ...
  )
  print(p)
}

#' @exportS3Method
summary.cocoBoot <- function(object, ...) {
  sum_obj <- list(
    model_type    = object$type,
    model_order   = object$order,
    ts_length     = length(object$ts),
    covariates    = if (is.null(object$cov)) {
      "None"
    } else {
      if (is.data.frame(object$cov) || is.matrix(object$cov)) {
        paste("Present (", ncol(as.matrix(object$cov)), " variable(s))", sep = "")
      } else {
        "Present (1 variable)"
      }
    },
    out_of_region = sum(object$df_plot$y < object$df_plot$lower | object$df_plot$y > object$df_plot$upper),
    total_rows    = nrow(object$df_plot),
    bootstrap_draws = as.numeric(object$rep.Bootstrap),
    boot_duration = as.numeric(object$duration, units = "secs")
  )
  
  class(sum_obj) <- "summary.cocoBoot"
  return(sum_obj)
}

#' @exportS3Method
print.summary.cocoBoot <- function(x, ...) {
  cat("---- Bootstrap Assessment Summary ----\n")
  cat("Model Type:         ", x$model_type, "\n")
  cat("Model Order:        ", x$model_order, "\n")
  cat("Time Series Length: ", x$ts_length, "\n")
  cat("Covariates:         ", x$covariates, "\n")
  cat("Bootstrap Draws:    ", x$bootstrap_draws, "\n")
  cat("Bootstrapping Duration: ", x$boot_duration, " seconds\n", sep = "")
  cat("ACF values outside bootstrapped confidence region: ", 
      x$out_of_region, " out of ", x$total_rows, "\n")
  cat("\n")
  invisible(x)
}


# Parameter checks of cocoSim / cocoSim_base (no covariates), with the same
# messages, run once for all bootstrap replicates.
bootCheckParBase <- function(type, order, par) {
  if (order == 2) {
    if (par[2] + par[4] > 1) {
      stop("The condition alpha_1 + alpha_3 < 1 ist not satisfied")
    }
    if ((par[2] < 0) | (par[3] < 0) | (par[4] < 0) | (par[2] >= 1) | (par[3] >= 1) | (par[4] >= 1)) {
      stop("The alpha parameters must be within the unit interval")
    }
    if ((type == "GP") & ((par[5] < 0) | (par[5] >= 1))) {
      stop("eta must be within the unit interval")
    }
  } else {
    if ((par[2] < 0) | (par[2] >= 1)) {
      stop("alpha must be within the unit interval")
    }
    if ((type == "GP") & ((par[3] < 0) | (par[3] >= 1))) {
      stop("eta must be within the unit interval")
    }
  }
  if ((type != "GP") & (type != "Poisson")) {
    stop("Option 'type' must be either Poisson or GP")
  }
  if ((order != 1) & order != (2)) {
    stop("Option 'order' must be 1 or 2")
  }
  if (!is.numeric(par)) {
    stop("Option 'par' must be a numeric vector")
  }
  if (par[1] <= 0) {
    stop("lambda must be a positive reel number")
  }
  if (!((all(par[-1] > 0)) & all(par[-1] < 1))) {
    stop(" The autoregressive parameters and eta must be between zero and one")
  }
  n_par <- (if (type == "GP") 3 else 2) + 2 * (order - 1)
  if (length(par) != n_par) {
    stop(paste("Number of parameteres must equal", n_par, "for the", type, order, "model"))
  }
  invisible(TRUE)
}

# Parameter checks of cocoSim / cocoSim_cov (covariates), with the same
# messages, run once for all bootstrap replicates.
bootCheckParCov <- function(type, order, par, xreg) {
  if (order == 2) {
    if (par[1] + par[3] > 1) {
      stop("The condition alpha_1 + alpha_3 < 1 ist not satisfied")
    }
    if ((par[1] < 0) | (par[2] < 0) | (par[3] < 0) | (par[1] >= 1) | (par[2] >= 1) | (par[3] >= 1)) {
      stop("The alpha parameters must be within the unit interval")
    }
    if ((type == "GP") & ((par[4] < 0) | (par[4] >= 1))) {
      stop("eta must be within the unit interval")
    }
  } else {
    if ((par[1] < 0) | (par[1] >= 1)) {
      stop("alpha must be within the unit interval")
    }
    if ((type == "GP") & ((par[2] < 0) | (par[2] >= 1))) {
      stop("eta must be within the unit interval")
    }
  }
  if ((type != "GP") & (type != "Poisson")) {
    stop("Option 'type' must be either Poisson or GP")
  }
  if ((order != 1) & order != (2)) {
    stop("Option 'order' must be 1 or 2")
  }
  if (!is.numeric(par)) {
    stop("Option 'par' must be a numeric vector")
  }
  if (any(par[-c((length(par) - ncol(xreg)):length(par))] < 0) |
    any(par[-c((length(par) - ncol(xreg)):length(par))] > 1)) {
    stop("The autoregressive parameters and eta must be between zero and one")
  }
  n_par <- (if (order == 1) 1 else 3) + (type == "GP") + ncol(xreg)
  if (length(par) != n_par) {
    stop(paste("Number of parameters must equal", n_par, "for the", type, order,
               "model with", ncol(xreg), "covariates"))
  }
  invisible(TRUE)
}

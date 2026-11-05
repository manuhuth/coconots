#' @title Probability Integral Transform Based Model Assessment Procedure
#' @description Computes the probability integral transform (PIT) and provides
#' the non-randomized PIT histogram for assessing absolute performance of a
#' fitted model as proposed by Czado et al. (2009).
#' @param coco An object of class coco
#' @param J Number of bins for the histogram (default: 10)
#' @param conf.alpha Significance level for the confidence intervals (default: 0.05)
#' @param julia  if TRUE, the PIT is computed with \proglang{julia} (default: FALSE)
#' @return an object of class cocoPit. It contains the probability integral
#' transform values,  p-value of the chi-square goodness of fit test and information on the model specifications.
#' @details The adequacy of a distributional assumption for a model is assessed by
#' checking the cumulative non-randomized PIT distribution for uniformity.
#' A useful graphical device is the PIT histogram, which displays this
#' distribution to J equally spaced bins. We supplement the graph by
#' incorporating approximately \eqn{100(1 - \alpha)\%} confidence intervals obtained
#' from a standard chi-square goodness-of-fit test of the null hypothesis that
#' the J bins of the histogram are drawn from a uniform distribution.
#' For details, see Jung, McCabe and Tremayne (2016).
#' 
#' @importFrom stats pchisq
#' @importFrom stats qchisq
#' @importFrom ggplot2 autoplot
#' 
#' @references 
#' Czado, C., Gneiting, T. and Held, L. (2009) Predictive model assessment for count data. \emph{Biometrics} \bold{65}, 1254--61.
#' 
#' Jung, R. C., McCabe, B.P.M. and Tremayne, A.R. (2016). Model validation and diagnostics. \emph{In Handbook of Discrete
#' Valued Time Series}. Edited by Davis, R.A., Holan, S.H., Lund, R. and Ravishanker, N.. Boca Raton: Chapman and
#' Hall, pp. 189--218.
#' 
#' Jung, R. C. and Tremayne, A. R. (2011) Convolution-closed models for count time series with applications. \emph{Journal of Time Series Analysis}, \bold{32}, 3, 268--280.
#' @author Manuel Huth
#' @examples
#' lambda <- 1
#' alpha <- 0.4
#' set.seed(12345)
#' data <- cocoSim(order = 1, type = "Poisson", par = c(lambda, alpha), length = 100)
#' fit <- cocoReg(order = 1, type = "Poisson", data = data)
#'
#' #PIT R implementation
#' pit_r <- cocoPit(fit)
#' plot(pit_r)
#' @export

cocoPit <- function(coco, J = 10, conf.alpha = 0.05, julia=FALSE) {
  start.time <- Sys.time()
  alpha <- conf.alpha

  if ((J != round(J)) | (J < 1)) {
    stop("The value of J must be a positive integer")
  }

  coco_check_count_limit(coco$ts)
  coco_check_length(length(coco$ts), coco$order + 1, "cocoPit",
                    sprintf("an order-%d model transforms the observations after the first %d",
                            coco$order, coco$order))
  
  if (!is.null(coco$julia_reg) & julia){
    data <- coco$ts
    addJuliaFunctions()
    coco_pit <- JuliaConnectoR::juliaGet( JuliaConnectoR::juliaCall("cocoPit", coco$julia_reg, J))
    J_test <- J
    u <- getJuliaValue(coco_pit, "bins")
    d <- getJuliaValue(coco_pit, "Pit_values")
  } else {

  data <- coco$ts
  seasonality <- c(1, 2) #will be used as argument in future versions coco$seasonality

  if ( is.null(coco$cov) ) {
    if ((coco$type == "GP") & (coco$order == 2)) {
      par <- coco$par
      lambda <- rep(par[1], length(data))
      alpha1 <- par[2]
      alpha2 <- par[3]
      alpha3 <- par[4]
      eta <- par[5]
      U <- 1 / (1 - alpha1 - alpha2 - alpha3)
      xreg <- NULL
    }

    if ((coco$type == "Poisson") & (coco$order == 2)) {
      par <- coco$par
      lambda <- rep(par[1], length(data))
      alpha1 <- par[2]
      alpha2 <- par[3]
      alpha3 <- par[4]
      U <- 1 / (1 - alpha1 - alpha2 - alpha3)
      xreg <- NULL
    }

    if ((coco$type == "GP") & (coco$order == 1)) {
      par <- coco$par
      lambda <- rep(par[1], length(data))
      alpha <- par[2]
      eta <- par[3]
      xreg <- NULL
    }

    if ((coco$type == "Poisson") & (coco$order == 1)) {
      par <- coco$par
      lambda <- rep(par[1], length(data))
      alpha <- par[2]
      xreg <- NULL
    }
  }



  if ( !is.null(coco$cov) ) {
    if ((coco$type == "Poisson") & (coco$order == 1)) {
      par <- coco$par
      alpha <- par[1]

      vec_lambda <- par[-(1)]
      xreg <- coco$cov


      lambda <- applyLinkFunctionVec(xreg %*% vec_lambda, coco$link_function)
    }

    if ((coco$type == "GP") & (coco$order == 1)) {
      par <- coco$par
      alpha <- par[1]
      eta <- par[2]

      vec_lambda <- par[-(1:2)]
      xreg <- coco$cov


      lambda <- applyLinkFunctionVec(xreg %*% vec_lambda, coco$link_function)
    }

    if ((coco$type == "Poisson") & (coco$order == 2)) {
      par <- coco$par
      alpha1 <- par[1]
      alpha2 <- par[2]
      alpha3 <- par[3]
      vec_lambda <- par[-(1:3)]
      xreg <- coco$cov

      U <- 1 / (1 - alpha1 - alpha2 - alpha3)

      lambda <- applyLinkFunctionVec(xreg %*% vec_lambda, coco$link_function)
    }

    if ((coco$type == "GP") & (coco$order == 2)) {
      par <- coco$par
      alpha1 <- par[1]
      alpha2 <- par[2]
      alpha3 <- par[3]
      eta <- par[4]

      vec_lambda <- par[-(1:4)]
      xreg <- coco$cov

      U <- 1 / (1 - alpha1 - alpha2 - alpha3)

      lambda <- applyLinkFunctionVec(xreg %*% vec_lambda, coco$link_function)
    }
  }





  ### PIT
  T <- length(data)
  J_test <- J
  J <- J + 1 # add 1 because the first bin is not used
  meanPIT <- matrix(0, nrow = T, ncol = J)
  u <- seq(0, 1, length = J)

  if (coco$type == "Poisson") {
    eta <- 0
  }
  if (coco$order == 1) {
    validate_gp1_par(lambda, alpha, eta)
  } else {
    validate_gp2_par(lambda, alpha1, alpha2, alpha3, eta)
  }
  if (any((round(data) != data) | (data < 0))) {
    stop("y must be a non-negative integer")
  }

  pit_row <- function(t) {
    # Px = P(X_t <= x_t), Px1 = P(X_t <= x_t - 1): both are prefixes of the
    # same cumulative sum of transition densities, computed in one C++ call
    if (coco$order == 1) {
      dvec <- dGP1h_vec(data[t], data[t - seasonality[1]], lambda[t], alpha, eta)
    } else {
      dvec <- dGP2h_vec(data[t], data[t - seasonality[1]], data[t - seasonality[2]],
                        lambda[t], alpha1, alpha2, alpha3, eta)
    }
    coco_check_finite_pmf(dvec, "predictive pmf")
    cdf <- cumsum(dvec)
    Px <- cdf[data[t] + 1]
    Px1 <- if (data[t] >= 1) cdf[data[t]] else 0

    temp <- (u - Px1) / (Px - Px1)
    temp[is.na(temp)] <- 1
    temp[temp < 0] <- 0
    temp[temp > 1] <- 1
    temp
  }

  # Every observation with a full history enters the PIT: t = order + 1, ..., T
  # (t = 2 for order 1, t = 3 for order 2), as in Coconots.jl's cocoPit
  first_t <- seasonality[coco$order] + 1
  if (is.null(coco$cov) && T >= first_t) {
    # Without covariates (lambda[t] constant) row t depends on t only through
    # the integer counts (x_t, x_{t-1}[, x_{t-2}]): one evaluation per
    # distinct triple, copied to all t where it occurs
    tt <- first_t:T
    key <- if (coco$order == 1) {
      paste(data[tt], data[tt - seasonality[1]])
    } else {
      paste(data[tt], data[tt - seasonality[1]], data[tt - seasonality[2]])
    }
    first <- which(!duplicated(key))
    rows_u <- vapply(tt[first], pit_row, numeric(J))
    meanPIT[tt, ] <- t(rows_u[, match(key, key[first]), drop = FALSE])
  } else {
    tt <- if (T >= first_t) first_t:T else integer(0)
    for (t in tt) {
      meanPIT[t, ] <- pit_row(t)
    } # end t
  }

  # mean over the time points used (rows t < first_t are not part of the PIT)
  PIT <- colMeans(meanPIT[tt, , drop = FALSE])
  d <- diff(PIT)
  u <- u[-c(1)]
  } #end julia
  
  # number of PIT time points (t = order + 1, ..., T); d sums to one, so
  # d * n_pit are the bin counts of the non-randomized PIT histogram
  n_pit <- length(data) - coco$order
  pval <- pitChisqPValue(d, n_pit, J_test)
  diff <- (sqrt(stats::qchisq(1-alpha, J_test-1)/J_test^2*n_pit) ) / n_pit
  confidence_band_values <- 1 / J_test + c(-1, 1) * diff 
  
  end.time <- Sys.time()

  time <- end.time - start.time
  list_out <- list(
    "type" = coco$type, "order" = coco$order, "ts" = coco$ts, "alpha" = alpha,
    "par" = par, "PIT values" = d, "duration" = time, "p_value" = pval,
    "confidence_bands" = confidence_band_values, "J" = J
  )
  
  class(list_out) <- "cocoPit"
  return(list_out)
}

# p-value of the chi-square goodness-of-fit test of uniformity of the J bins
# of a PIT histogram d (bin masses summing to one) over n_pit time points:
# upper tail P(chi^2_{J-1} >= stat) with stat = sum((n d - n/J)^2 / (n/J)).
# Used for the R and the Julia PIT alike.
pitChisqPValue <- function(d, n_pit, J) {
  stat <- sum((d * n_pit - n_pit / J)^2 / (n_pit / J))
  stats::pchisq(stat, J - 1, lower.tail = FALSE)
}

#' @importFrom ggplot2 autoplot
#' @export
ggplot2::autoplot
#' @exportS3Method
autoplot.cocoPit <- function(object, ...){
  
  df <- data.frame(object$`PIT values`, 1:rep(length(object$`PIT values`)),
                   rep(object$confidence_bands[1], length(object$`PIT values`)),
                   rep(object$confidence_bands[2], length(object$`PIT values`)))
  colnames(df) <- c("pit", "bins", "lower", "upper")
  df_bands <- df
  df_bands$bins[1] <- df$bins[1] - 0.25
  df_bands$bins[length(df$bins)] <- df$bins[length(df$bins)] + 0.25
  pl <- ggplot2::ggplot(data=df, mapping = ggplot2::aes_string(x = "bins", y = "pit")) +
    ggplot2::geom_bar(stat="identity", position="dodge", width=0.3) + 
    ggplot2::labs(title = "Pit histogram", x = "Bins", y = "") +
    ggplot2::scale_x_continuous(breaks=df$bins) +
    ggplot2::theme_bw() + ggplot2::theme(text = ggplot2::element_text(size = 20)) +
    ggplot2::geom_ribbon(data=df_bands, ggplot2::aes_string(ymin="lower",ymax="upper"),
                         fill="steelblue", alpha=0.3) +
    ggplot2::ylim(c(0, min(1.05, max(df$pit)*1.5)))
  pl
}

#' @exportS3Method
plot.cocoPit <- function(x, ...) {
  p <- autoplot(
    x,
    ...
  )
  print(p)
}

#'@exportS3Method
print.cocoPit <- function(x, ...) {
  print(autoplot(x, ...,))
  invisible(x)
}

#' @exportS3Method
print.cocoPit <- function(x, ...) {
  cat("cocoPit Object")
  cat("\nModel Type: ", x$type, sep = "")
  cat("\nModel Order: ", x$order, sep = "")
  cat("\nPIT p-value: ", x$p_value, sep = "")
  cat("\nNumber of bins: ", x$J-1, sep = "")
  cat("\nBootstrapping Duration: ", as.numeric(x$duration, units = "secs"), " seconds", sep = "")
  cat("\n")
  
  invisible(x)
}

#' @exportS3Method
summary.cocoPit <- function(object, ...) {
  sum_obj <- list(
    model_type          = object$type,
    model_order         = object$order,
    ts_length           = length(object$ts),
    significance_level  = object$alpha,
    pit_values_summary  = summary(object[["PIT values"]]),
    confidence_bands    = object$confidence_bands,
    n_bins              = object$J - 1,
    pit_pvalue          = object$p_value
  )
  
  class(sum_obj) <- "summary.cocoPit"
  return(sum_obj)
}

#' @exportS3Method
print.summary.cocoPit <- function(x, ...) {
  cat("---- PIT Assessment Summary ----\n")
  cat("Model Type:          ", x$model_type, "\n")
  cat("Model Order:         ", x$model_order, "\n")
  cat("Time Series Length:  ", x$ts_length, "\n")
  cat("Significance Level:  ", x$significance_level, "\n\n")
  
  cat("Confidence Bands:    [", x$confidence_bands[1], ", ", x$confidence_bands[2], "]\n", sep = "")
  cat("Number of bins:      ", x$n_bins, "\n\n")
  
  cat("PIT Values Summary:\n")
  print(x$pit_values_summary)
  cat("\nPIT p-value:         ", x$pit_pvalue, "\n")
  
  invisible(x)
}


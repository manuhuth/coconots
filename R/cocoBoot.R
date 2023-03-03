#' @title Bootstrap Confidence Intervals for Autocorrelations of a COCO Model
#' @description Computes bootstrap confidence intervals for the autocorrelations of a COCO model. The function can handle both Poisson and GP models, with order 1 or 2. The function also has options for the number of lags, the number of bootstrap replicates, the confidence level, and the labels for the plot.
#' @param coco An object of class coco.fit or coco.fit.c
#' @param numb.lags Number of lags for which to compute autocorrelations
#' @param rep.Bootstrap Number of bootstrap replicates to use
#' @param confidence Confidence level for the intervals
#' @param plot_main Plot title
#' @param xlab X-axis label for the plot
#' @param ylab Y-axis label for the plot
#' @return A plot of the autocorrelations with bootstrap confidence intervals
#' @references 
#' Tsay, R. S. (1992) Model checking via parametric bootstraps in time series analysis. \emph{Applied Statistics} \bold{41}, 1--15.
#' @export

cocoBoot <- function(coco, numb.lags = 21, rep.Bootstrap = 400,
                 confidence = 0.95, plot_main="Bootstrap", xlab = "Lag", 
                 ylab= "Autocorrelation", julia = FALSE
                 ) {
  start.time <- Sys.time()

  if ((confidence <= 0) | (confidence >= 1)) {
    stop("Option confidence must be a real number between 0 and 1")
  }

  if ((numb.lags != round(numb.lags)) | (numb.lags < 1)) {
    stop("The value of numb.lags must be a positive integer")
  }

  if ((rep.Bootstrap != round(rep.Bootstrap)) | (rep.Bootstrap < 1)) {
    stop("The value of rep.Bootstrap must be a positive integer")
  }
  
  if (!is.null(coco$julia_reg) & julia){
    coco_boot <- JuliaConnectoR::juliaGet( JuliaConnectoR::juliaCall("cocoBoot", coco$julia_reg, 1:numb.lags, rep.Bootstrap, 1-confidence))
    acfdata <- coco_boot$values[[2]]
    confidence <- data.frame(cbind(coco_boot$values[[3]], coco_boot$values[[4]]))
    colnames(confidence) <- c("lower", "upper")
  } else {
    data <- coco$ts
    seasonality <- coco$seasonality
  
    conf.alpha <- 1 - confidence
    
    if ( (methods::is(coco, "coco.fit")) ) {
      if ((coco$type == "GP") & (coco$order == 2)) {
        par <- coco$par
        lambda <- par[1]
        alpha1 <- par[2]
        alpha2 <- par[3]
        alpha3 <- par[4]
        eta <- par[5]
        U <- 1 / (1 - alpha1 - alpha2 - alpha3)
        xreg <- NULL
      }
  
      if ((coco$type == "Poisson") & (coco$order == 2)) {
        par <- coco$par
        lambda <- par[1]
        alpha1 <- par[2]
        alpha2 <- par[3]
        alpha3 <- par[4]
        U <- 1 / (1 - alpha1 - alpha2 - alpha3)
        xreg <- NULL
      }
  
      if ((coco$type == "GP") & (coco$order == 1)) {
        par <- coco$par
        lambda <- par[1]
        alpha <- par[2]
        eta <- par[3]
        xreg <- NULL
      }
  
      if ((coco$type == "Poisson") & (coco$order == 1)) {
        par <- coco$par
        lambda <- par[1]
        alpha <- par[2]
        xreg <- NULL
      }
    }
  
  
    if ( (methods::is(coco, "coco.fit.c")) )  {
      if ((coco$type == "Poisson") & (coco$order == 1)) {
        par <- coco$par
        alpha <- par[1]
  
        vec_lambda <- par[-(1)]
        xreg <- coco$cov
        data <- coco$ts
  
        # set up values for lambda
        lambda <- c()
        for (j in 1:length(data)) {
          lambda[j] <- exp(as.numeric(as.vector(xreg[j, ])) %*% vec_lambda)
        }
      }
  
      if ((coco$type == "GP") & (coco$order == 1)) {
        par <- coco$par
        alpha <- par[1]
        eta <- par[2]
        vec_lambda <- par[-(1:2)]
        xreg <- coco$cov
        data <- coco$ts
  
        lambda <- c()
        for (j in 1:length(data)) {
          lambda[j] <- exp(as.numeric(as.vector(xreg[j, ])) %*% vec_lambda)
        }
      }
  
      if ((coco$type == "Poisson") & (coco$order == 2)) {
        par <- coco$par
        alpha1 <- par[1]
        alpha2 <- par[2]
        alpha3 <- par[3]
        vec_lambda <- par[-(1:3)]
        xreg <- coco$cov
        data <- coco$ts
        U <- 1 / (1 - alpha1 - alpha2 - alpha3)
  
        lambda <- c()
        for (j in 1:length(data)) {
          lambda[j] <- exp(as.numeric(as.vector(xreg[j, ])) %*% vec_lambda)
        }
      }
  
      if ((coco$type == "GP") & (coco$order == 2)) {
        par <- coco$par
        alpha1 <- par[1]
        alpha2 <- par[2]
        alpha3 <- par[3]
        eta <- par[4]
  
        vec_lambda <- par[-(1:4)]
        xreg <- coco$cov
        data <- coco$ts
        U <- 1 / (1 - alpha1 - alpha2 - alpha3)
  
        lambda <- c()
        for (j in 1:length(data)) {
          lambda[j] <- exp(as.numeric(as.vector(xreg[j, ])) %*% vec_lambda)
        }
      }
    }
  
  
    # Parametric Bootstrap
    T <- length(data)
    nlags <- numb.lags
    nB <- rep.Bootstrap
    B <- matrix(NaN, nrow = T, ncol = nB)
    ac <- matrix(NaN, nrow = nlags, ncol = nB)
    conf <- conf.alpha
  
    if ((is.null(coco$cov)) ) {
      for (b in 1:nB) {
        help <- cocoSim(order = coco$order, type = coco$type, par = par, length = (T + 10), seasonality = seasonality)$data[11:(T + 10)]
        B[, b] <- help
        ac[, b] <- forecast::Acf(help, plot = FALSE, lag.max = nlags)$acf[2:(nlags + 1)]
      }
    }
  
  
    if (!is.null(coco$cov)) {
      xreg <- as.matrix(xreg)
      for (b in 1:nB) {
        help <- cocoSim(type = coco$type, order = coco$order, par = par, length = T, xreg = xreg, seasonality = seasonality)$data
        B[, b] <- help
        ac[, b] <- forecast::Acf(help, plot = FALSE, lag.max = nlags)$acf[2:(nlags + 1)]
      }
    }
  

    means <- rowMeans(ac)
    var <- apply(ac, 1, var)
  
    confidence <- matrix(NaN, nrow = nlags, ncol = 2)
    colnames(confidence) <- c("lower", "upper")
    for (j in 1:nlags) {
      upper <- stats::qnorm(1 - conf / 2, means[j], var[j]^0.5)
      lower <- stats::qnorm(conf / 2, means[j], var[j]^0.5)
      confidence[j, ] <- c(lower, upper)
    }
    acfdata <- forecast::Acf(data, plot = FALSE, lag.max = nlags)$acf[2:(nlags + 1)]
  } #end julia
  
  max <- max(c(acfdata, confidence[, "upper"])) + 0.5 * abs(max(c(acfdata, confidence[, "upper"])))
  
  if (max >= 1.1) {
    max <- 1.1
  }

  min <- min(c(acfdata, confidence[, "lower"])) - 0.5 * abs(min(c(acfdata, confidence[, "lower"])))
  if (max <= -1.1) {
    max <- -1.1
  }

  plot(acfdata, ylab = ylab, xlab = xlab, ylim = c(min, max), main = plot_main)
  graphics::points(confidence[, 1], pch = 3, col = c("red"))
  graphics::points(confidence[, 2], pch = 3, col = c("red"))
  q <- grDevices::recordPlot()

  end.time <- Sys.time()
  time <- end.time - start.time
  list <- list(
    "type" = coco$type, "order" = coco$order, "ts" = coco$ts, "cov" = coco$cov, means = acfdata,
    "PBT.plot" = q, "confidence" = confidence, "duration" = time
  )

  return(list)
} # end function

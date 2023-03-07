#' @title Probability integral transform plot for coco
#' @description Computes the probability integral transform (PIT) and provides the non-randomized PIT histogram for assessing absolute performance of a fitted model as proposed by Czado et al. (2009).
#' @param coco An object of class coco
#' @param J Number of bins for the histogram (default: 10)
#' @param ylab Label for the y-axis (default: "Relative frequency")
#' @param xlab Label for the x-axis (default: "Probability integral transform")
#' @param plot_main Title for the plot (default: "PIT")
#' @param julia  if TRUE, the PIT is computed with Julia.
#' @return A plot of the probability integral transform for the coco object.
#' @references 
#' Czado, C., Gneiting, T. and Held, L. (2009) Predictive model assessment for count data. \emph{Biometrics} \bold{65}, 1254--61.
#' 
#' Jung, R. C. and Tremayne, A. R. (2010) Convolution-closed models for count time series with applications. \emph{Journal of Time Series Analysis}, \bold{32}, 3, 268--280.
#' @author Manuel Huth
#' @examples
#' lambda <- 1
#' alpha <- 0.4
#' set.seed(12345)
#' data <- cocoSim(order = 1, type = "Poisson", par = c(lambda, alpha), length = 100)$data
#' #julia_installed = TRUE ensures that the fit object
#' #is compatible with the julia cocoPit implementation 
#' fit <- cocoReg(order = 1, type = "Poisson", data = data)
#'
#' #PIT R implementation
#' pit_r <- cocoPit(fit)
#' @export

cocoPit <- function(coco, J = 10, ylab = "Relative frequency",
                    xlab = "Probability integral transform", plot_main = "PIT", julia=FALSE) {
  start.time <- Sys.time()

  

  if ((J != round(J)) | (J < 1)) {
    stop("The value of J must be a positive integer")
  }
  
  if (!is.null(coco$julia_reg) & julia){
    addJuliaFunctions()
    coco_pit <- JuliaConnectoR::juliaGet( JuliaConnectoR::juliaCall("cocoPit", coco$julia_reg,J))
    coco_pit$keys
    u <- coco_pit$values[[2]]
    d <- coco_pit$values[[1]]
  } else {

  data <- coco$ts
  seasonality <- coco$seasonality

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

      U <- 1 / (1 - alpha1 - alpha2 - alpha3)

      lambda <- c()
      for (j in 1:length(data)) {
        lambda[j] <- exp(as.numeric(as.vector(xreg[j, ])) %*% vec_lambda)
      }
    }
  }





  ### PIT
  T <- length(data)
  J <- J + 1 # add 1 because the first bin is not used
  meanPIT <- matrix(0, nrow = T, ncol = J)

  for (t in (seasonality[2] + 1):T) {
    if ((coco$type == "Poisson") & (coco$order == 1)) {
      part <- c(lambda[t], alpha, 0)
      Px <- pGP1(data[t], data[t - seasonality[1]], part)
      if (data[t] >= 1) {
        Px1 <- pGP1(data[t] - 1, data[t - seasonality[1]], part)
      }
      if (data[t] < 1) {
        Px1 <- 0
      }
    }

    if ((coco$type == "GP") & (coco$order == 1)) {
      part <- c(lambda[t], alpha, eta)
      Px <- pGP1(data[t], data[t - seasonality[1]], part)
      if (data[t] >= 1) {
        Px1 <- pGP1(data[t] - 1, data[t - seasonality[1]], part)
      }
      if (data[t] < 1) {
        Px1 <- 0
      }
    }

    if ((coco$type == "Poisson") & (coco$order == 2)) {
      part <- c(lambda[t], alpha1, alpha2, alpha3, 0)
      Px <- pGP2(data[t], data[t - seasonality[1]], data[t - seasonality[2]], part)
      if (data[t] >= 1) {
        Px1 <- pGP2(data[t] - 1, data[t - seasonality[1]], data[t - seasonality[2]], part)
      }
      if (data[t] < 1) {
        Px1 <- 0
      }
    }

    if ((coco$type == "GP") & (coco$order == 2)) {
      part <- c(lambda[t], alpha1, alpha2, alpha3, eta)
      Px <- pGP2(data[t], data[t - seasonality[1]], data[t - seasonality[2]], part)
      if (data[t] >= 1) {
        Px1 <- pGP2(data[t] - 1, data[t - seasonality[1]], data[t - seasonality[2]], part)
      }
      if (data[t] < 1) {
        Px1 <- 0
      }
    }

    u <- seq(0, 1, length = J)
    N <- 1
    PIT <- matrix(0, nrow = N, ncol = J)
    for (i in 1:N) {
      temp <- numeric(J)
      for (s in 1:length(u)) {
        temp[s] <- (u[s] - Px1[i]) / (Px[i] - Px1[i])
        if (is.na(temp[s])) temp[s] <- 1
        if (temp[s] < 0) temp[s] <- 0
        if (temp[s] > 1) temp[s] <- 1
      }
      PIT[i, ] <- temp
    }
    meanPIT[t, ] <- colMeans(PIT)
  } # end t

  PIT <- colMeans(meanPIT)
  d <- diff(PIT)
  u <- u[-c(1)]
  } #end julia
  
  # pdf("PIT_mle_poisson.pdf", height = 5, width = 5)
  plot(u, d,
    type = "h", lend = 1, lwd = 18, ylab = ylab, xlab = xlab, main = plot_main,
    ylim = c(0, max(d) * 2), xaxt = "n", col = "navyblue", yaxs = "i"
  )
  graphics::axis(1, at = seq(0.1, 0.9, 0.1), labels = as.character(seq(0.1, 0.9, 0.1)))
  p <- grDevices::recordPlot()

  end.time <- Sys.time()

  time <- end.time - start.time
  list <- list(
    "type" = coco$type, "order" = coco$order, "ts" = coco$ts,
    "par" = par, "PIT" = p, "PIT values" = d, "duration" = time
  )

  return(list)
}

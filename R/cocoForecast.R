cocoForecast <- function(coco, max=10, xcast=NULL, plot = TRUE, title = "Probability mass",
             xlab="x", ylab="Probabilities", width_bars = 0.04, seasonality=c(1,2), decimals = 4) {

  if ((class(coco) != "coco.fit") & (class(coco) != "coco.fit.c")) {
    stop("The coco object must be from class coco.fit or coco.fit.c")
  }
  
  if (length(seasonality == 1)) {
    seasonality = c(seasonality, seasonality + 1)
  }
  
  data <- coco$ts
  y <- data[length(data) - seasonality[1] + 1]
  z <- data[length(data) - seasonality[2] + 1]
  parameter <-coco$par
  
  if (class(coco) == "coco.fit.c") { 
    number_covariates <- ncol(coco$cov) 
    betas <- parameter[(length(parameter)-number_covariates+1):length(parameter)]
    parameter <- head(parameter, -number_covariates)
    dot_product <- betas %*% c(xcast)
    lambda <- exp(dot_product)
    parameter <- c(lambda, parameter)
  }
  
  type <- coco$type
  order <- coco$order
  
  densities <- c()
  
  x <- 0:max
  
  if (order == 1){
    if (type == "Poisson"){parameter <- c(parameter,0)}
    for (number in x) { 
      densities[number+1] <- dGP1(x=number, y=y, par=parameter)
    }
  } 
  
  if (order == 2){
    if (type == "Poisson"){parameter <- c(parameter,0)}
    for (number in x) { 
      densities[number+1] <- dGP2(x=number, y=y, z=z, par=parameter)
    }
  }
  
  densities <- round(densities, decimals)
  
  if (isTRUE(plot)) {
  plot <- ggplot2::ggplot(mapping = ggplot2::aes(x = x, y = densities)) + ggplot2::geom_bar(stat="identity", position="dodge", width=width_bars) + 
    ggplot2::labs(title = title, x = xlab, y = ylab) + ggplot2::scale_x_continuous(breaks=x)
  }
  
  mode <- match(max(densities), densities) - 1
  
  distribution_function <- cumsum(densities)
  
  median <- min(which(distribution_function >= 0.5)) - 1
  
  out <- list("plot" = plot, "density" = densities, "mode" = mode, "median" = median)
  
  return(out) 
}

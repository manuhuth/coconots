cocoResid <- function(coco, val.num = 1e-5) {
  start.time <- Sys.time()

  if ((class(coco) != "coco.fit") & (class(coco) != "coco.fit.c")) {
    stop("The coco object must be from class coco.fit or coco.fit.c")
  }

  if (val.num <= 0) {
    stop("Option val.num must be a non-negative real number")
  }

  seasonality <- coco$seasonality

  if (class(coco) == "coco.fit") {
    xreg <- NULL
    ### PAR1
    if ((coco$type == "Poisson") & (coco$order == 1)) {
      par <- coco$par
      lambda <- par[1]
      alpha <- par[2]
      data <- coco$ts

      # Fitted Values
      fitted <- c(rep(mean(data), seasonality[1])) 
      varX <- c(rep(var(data), seasonality[1])) 
      for (j in (seasonality[1] + 1):(length(data) + 1)) {
        fitted[j] <- alpha * data[j - seasonality[1]] + lambda
        varX[j] <- alpha * (1 - alpha) * data[j - seasonality[1]] + lambda
      } # end for
    } # end if PAR1


    ### GP1
    if ((coco$type == "GP") & (coco$order == 1)) {
      par <- coco$par
      lambda <- par[1]
      alpha <- par[2]
      eta <- par[3]
      data <- coco$ts

      # Fitted Values
      fitted <- c(rep(mean(data), seasonality[1])) 
      varX <- c(rep(var(data), seasonality[1])) 
      for (j in (seasonality[1] + 1):(length(data) + 1)) {
        fitted[j] <- alpha * data[j - seasonality[1]] + lambda / (1 - eta)
        # compute conditional variance
        y <- data[j - seasonality[1]]
        meanR <- alpha * y # mean of random operator
        help <- 0
        for (r in 1:y) {
          help <- help + (r - meanR)^2 * (choose(y, r) * alpha * (1 - alpha) * (alpha + eta * (1 - alpha) / lambda * r)^(r - 1)
            * (1 - alpha + eta * (1 - alpha) / lambda * (y - r))^(y - r - 1)
            / (1 + eta * (1 - alpha) / lambda * y)^(y - 1))
        } # end for r

        varX[j] <- help + lambda / (1 - eta)^3
      } # end for j
    } # end if GP1



    ### PAR2
    if ((coco$type == "Poisson") & (coco$order == 2)) {
      par <- coco$par
      lambda <- par[1]
      alpha1 <- par[2]
      alpha2 <- par[3]
      alpha3 <- par[4]
      U <- 1 / (1 - alpha1 - alpha2 - alpha3)
      data <- coco$ts
      eta <- 0

      # Fitted Values
      fitted <- c(rep(mean(data), seasonality[2])) 
      varX <- c(rep(var(data), seasonality[2])) 
      
      for (j in (seasonality[2] + 1):(length(data) + 1)) {
        y <- data[j - seasonality[1]]
        z <- data[j - seasonality[2]]
        # compute conditional mean of R
        help <- 0
        prob <- val.num + 1
        r <- 0
        while (prob > val.num) {
          prob <- dR2(r, y, z, lambda, alpha1, alpha2, alpha3, eta)
          help <- help + r * prob
          r <- r + 1
        }
        
        meanR <- help
        fitted[j] <- meanR + lambda / (1 - eta)
        
        # compute conditional variance
        help <- 0
        prob <- val.num + 1
        r <- 0
        while (prob > val.num) {
          prob <- dR2(r, y, z, lambda, alpha1, alpha2, alpha3, eta)
          help <- help + (r - meanR)^2 * prob
          r <- r + 1
        }
        
        varX[j] <- help + lambda / (1 - eta)^3
      } # end for j
    } # end if PAR2

    ### GP2
    if ((coco$type == "GP") & (coco$order == 2)) {
      par <- coco$par
      lambda <- par[1]
      alpha1 <- par[2]
      alpha2 <- par[3]
      alpha3 <- par[4]
      eta <- par[5]
      U <- 1 / (1 - alpha1 - alpha2 - alpha3)
      data <- coco$ts

      # Fitted Values
      fitted <- c(rep(mean(data), seasonality[2])) 
      varX <- c(rep(var(data), seasonality[2])) 

      for (j in (seasonality[2] + 1):(length(data) + 1)) {
        y <- data[j - seasonality[1]]
        z <- data[j - seasonality[2]]
        # compute conditional mean of R
        help <- 0
        prob <- val.num + 1
        r <- 0
        while (prob > val.num) {
          prob <- dR2(r, y, z, lambda, alpha1, alpha2, alpha3, eta)
          help <- help + r * prob
          r <- r + 1
        }

        meanR <- help
        fitted[j] <- meanR + lambda / (1 - eta)

        # compute conditional variance
        help <- 0
        prob <- val.num + 1
        r <- 0
        while (prob > val.num) {
          prob <- dR2(r, y, z, lambda, alpha1, alpha2, alpha3, eta)
          help <- help + (r - meanR)^2 * prob
          r <- r + 1
        }

        varX[j] <- help + lambda / (1 - eta)^3
      } # end for j
    } # end if GP2

    predVal <- fitted[length(fitted)]
    predVar <- varX[length(fitted)]
    pred <- c(predVal, predVar)
    names(pred) <- c("one step ahead forecast", "forecast variance")
    fitted <- fitted[-length(fitted)]
    varX <- varX[-length(fitted)]

    # Residuals
    residuals <- data - fitted

    # Pearson Residuals
    peResid <- residuals / varX^0.5

    end.time <- Sys.time()

    time <- end.time - start.time

    list <- list(
      "fitted" = fitted, "resdiuals" = residuals, "pe.resid" = peResid,
      "cond.var" = varX, "type" = coco$type, "order" = coco$order, "ts" = coco$ts,
      "par" = par, "prediction" = pred, "duration" = time
    )
  } # end no covariates

  if (class(coco) == "coco.fit.c") {
    ### PAR1
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

      # Fitted Values
      fitted <- c(rep(mean(data), seasonality[1])) 
      varX <- c(rep(var(data), seasonality[1])) 
      for (j in (seasonality[1] + 1):(length(data))) {
        fitted[j] <- alpha * data[j - seasonality[1]] + lambda[j]
        varX[j] <- alpha * (1 - alpha) * data[j - seasonality[1]] + lambda[j]
      } # end for
    } # end if PAR1


    ### GP1
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

      # Fitted Values
      fitted <- c(mean(data)) 
      varX <- c(var(data)) 
      for (j in (seasonality[1] + 1):(length(data))) {
        fitted[j] <- alpha * data[j - seasonality[1]] + lambda[j] / (1 - eta)
        # compute conditional variance
        y <- data[j - seasonality[1]]
        meanR <- alpha * y # mean of random operator
        help <- 0
        for (r in 1:y) {
          help <- help + (r - meanR)^2 * (choose(y, r) * alpha * (1 - alpha) * (alpha + eta * (1 - alpha) / lambda[j] * r)^(r - 1)
            * (1 - alpha + eta * (1 - alpha) / lambda[j] * (y - r))^(y - r - 1)
            / (1 + eta * (1 - alpha) / lambda[j] * y)^(y - 1))
        } # end for r

        varX[j] <- help + lambda[j] / (1 - eta)^3
      } # end for j
    } # end if GP1



    ### PAR2
    if ((coco$type == "Poisson") & (coco$order == 2)) {
      par <- coco$par
      alpha1 <- par[1]
      alpha2 <- par[2]
      alpha3 <- par[3]

      vec_lambda <- par[-(1:3)]
      xreg <- coco$cov
      data <- coco$ts
      U <- 1 / (1 - alpha1 - alpha2 - alpha3)
      eta <- 0

      lambda <- c()
      for (j in 1:length(data)) {
        lambda[j] <- exp(as.numeric(as.vector(xreg[j, ])) %*% vec_lambda)
      }


      # Fitted Values
      fitted <- c(rep(mean(data), seasonality[2])) 
      varX <- c(rep(var(data), seasonality[2]))
      
      for (j in (seasonality[2] + 1):(length(data))) {
        y <- data[j - seasonality[1]]
        z <- data[j - seasonality[2]]
        # compute conditional mean of R
        help <- 0
        prob <- val.num + 1
        r <- 0
        while (prob > val.num) {
          prob <- dR2(r, y, z, lambda[j], alpha1, alpha2, alpha3, eta)
          help <- help + r * prob
          r <- r + 1
        }
        
        meanR <- help
        fitted[j] <- meanR + lambda[j] / (1 - eta)
        
        # compute conditional variance
        meanR <- fitted[j] - lambda[j] # mean of random operator
        help <- 0
        prob <- val.num + 1
        r <- 0
        while (prob > val.num) {
          prob <- dR2(r, y, z, lambda[j], alpha1, alpha2, alpha3, eta)
          help <- help + (r - meanR)^2 * prob
          r <- r + 1
        }
        
        varX[j] <- help + lambda[j] / (1 - eta)^3
      } # end for j
    } # end if PAR2

    ### GP2
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


      # Fitted Values
      fitted <- c(rep(mean(data), seasonality[2])) 
      varX <- c(rep(var(data), seasonality[2])) 

      for (j in (seasonality[2] + 1):(length(data))) {
        y <- data[j - seasonality[1]]
        z <- data[j - seasonality[2]]
        # compute conditional mean of R
        help <- 0
        prob <- val.num + 1
        r <- 0
        while (prob > val.num) {
          prob <- dR2(r, y, z, lambda[j], alpha1, alpha2, alpha3, eta)
          help <- help + r * prob
          r <- r + 1
        }

        meanR <- help
        fitted[j] <- meanR + lambda[j] / (1 - eta)

        # compute conditional variance
        meanR <- fitted[j] - lambda[j] # mean of random operator
        help <- 0
        prob <- val.num + 1
        r <- 0
        while (prob > val.num) {
          prob <- dR2(r, y, z, lambda[j], alpha1, alpha2, alpha3, eta)
          help <- help + (r - meanR)^2 * prob
          r <- r + 1
        }

        varX[j] <- help + lambda[j] / (1 - eta)^3
      } # end for j
    } # end if GP2


    # Residuals
    residuals <- data - fitted

    # Pearson Residuals
    peResid <- residuals / varX^0.5

    end.time <- Sys.time()
    time <- end.time - start.time

    list <- list(
      "cond.expec" = fitted, "resdiuals" = residuals, "pe.resid" = peResid,
      "cond.var" = varX, "type" = coco$type, "order" = coco$order, "ts" = coco$ts, "cov" = xreg,
      "par" = par, "duration" = time
    )
  } # end covariates


  return(list)
} # end function

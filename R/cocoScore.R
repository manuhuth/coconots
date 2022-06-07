cocoScore <- function(coco, val.num = 1e-10) {

  if ( (class(coco) != "coco.fit") & (class(coco) != "coco.fit.c")) {
    stop("The coco object must be from class coco.fit or coco.fit.c")
  }

  T <- length(coco$ts)
  par <- coco$par
  seas <- coco$seasonality
  data <- coco$ts

if (class(coco) == "coco.fit") {


  #Poisson 1
  if ( (coco$type == "Poisson") & (coco$order == 1) ){
    par <- c(par, 0)
    log.score <- 0
    quad.score <- 0
    rps.score <- 0

    for (t in (seas[1]+1):T) {
      #log score
      log.score <- log.score - log(dGP1(data[t], data[t-seas[1]], par))

      #quadratic score
      norm.p <- 0
      help.stop <- 0

      j <- 0
      while (help.stop < (1-val.num) )  {
        help <- dGP1(j,data[t-seas[1]], par)
        help.stop <- pGP1(j,data[t-seas[1]], par)
        norm.p <- norm.p + help^2
        j <- j+1
      }
      quad.score <- quad.score - 2 * dGP1(data[t], data[t-seas[1]], par) + norm.p

      #ranked probability score
      sum <- 0
      help <- 0
      j <- 0
      while (help < (1-val.num) ){
        ind <- 0
        if (j >= data[t]) {
          ind <- 1
        }
        help <- pGP1(j, data[t-seas[1]], par)
        sum <- sum + (help - ind)^2
        j <- j + 1
      }
      rps.score <- rps.score + sum
    }# end t

    log.score <- log.score / (T-seas[1])
    quad.score <- quad.score / (T-seas[1])
    rps.score <- rps.score / (T-seas[1])

  }# Poisson 1


  #GP1
  if ( (coco$type == "GP") & (coco$order == 1) ){

    log.score <- 0
    quad.score <- 0
    rps.score <- 0

    for (t in (seas[1]+1):T) {
      #log score
      log.score <- log.score -log(dGP1(data[t], data[t-seas[1]], par))

      #quadratic score
      norm.p <- 0
      help.stop <- 0

      j <- 0
      while (help.stop < (1-val.num) )  {
        help <- dGP1(j,data[t-seas[1]], par)
        help.stop <- pGP1(j,data[t-seas[1]], par)
        norm.p <- norm.p + help^2
        j <- j+1
      }
      quad.score <- quad.score - 2* dGP1(data[t], data[t-seas[1]], par) + norm.p

      #ranked probability score
      sum <- 0
      help <- 0
      j <- 0
      while (help < (1-val.num) ){
        ind <- 0
        if (j >= data[t]) {
          ind <- 1
        }
        help <- pGP1(j, data[t-seas[1]], par)
        sum <- sum + (help - ind)^2
        j <- j + 1
      }
      rps.score <- rps.score + sum
  }# end t
    log.score <- log.score / (T-seas[1])
    quad.score <- quad.score / (T-seas[1])
    rps.score <- rps.score / (T-seas[1])
  }# GP1

  #Poisson 2
  if ( (coco$type == "Poisson") & (coco$order == 2) ){
    par <- c(par, 0)
    log.score <- 0
    quad.score <- 0
    rps.score <- 0

    for (t in (seas[2]+1):T) {
      #log score
      log.score <- log.score -log(dGP2(data[t], data[t-seas[1]], data[t-seas[2]], par))

      #quadratic score
      norm.p <- 0
      help.stop <- 0

      j <- 0
      while (help.stop < (1-val.num) )  {
        help <- dGP2(j,data[t-seas[1]], data[t-seas[2]], par)
        help.stop <- pGP2(j,data[t-seas[1]], data[t-seas[2]], par)
        norm.p <- norm.p + help^2
        j <- j+1
      }
      quad.score <- quad.score - 2* dGP2(data[t], data[t-seas[1]], data[t-seas[2]], par) + norm.p

      #ranked probability score
      sum <- 0
      help <- 0
      j <- 0
      while (help < (1-val.num) ){
        ind <- 0
        if (j >= data[t]) {
          ind <- 1
        }
        help <- pGP2(j, data[t-seas[1]],data[t-seas[2]], par)
        sum <- sum + (help - ind)^2
        j <- j + 1
      }
      rps.score <- rps.score + sum
    }# end t
    log.score <- log.score / (T-seas[2])
    quad.score <- quad.score / (T-seas[2])
    rps.score <- rps.score / (T-seas[2])
  }# Poisson 2

  #GP2
  if ( (coco$type == "GP") & (coco$order == 2) ){

    log.score <- 0
    quad.score <- 0
    rps.score <- 0

    for (t in (seas[2]+1):T) {
      #log score
      log.score <- log.score -log(dGP2(data[t], data[t-seas[1]], data[t-seas[2]], par))

      #quadratic score
      norm.p <- 0
      help.stop <- 0

      j <- 0
      while (help.stop < (1-val.num) )  {
        help <- dGP2(j,data[t-seas[1]], data[t-seas[2]], par)
        help.stop <- pGP2(j,data[t-seas[1]], data[t-seas[2]], par)
        norm.p <- norm.p + help^2
        j <- j+1
      }
      quad.score <- quad.score - 2* dGP2(data[t], data[t-seas[1]], data[t-seas[2]], par) + norm.p

      #ranked probability score
      sum <- 0
      help <- 0
      j <- 0
      while (help < (1-val.num) ){
        ind <- 0
        if (j >= data[t]) {
          ind <- 1
        }
        help <- pGP2(j, data[t-seas[1]],data[t-seas[2]], par)
        sum <- sum + (help - ind)^2
        j <- j + 1
      }
      rps.score <- rps.score + sum
    }# end t
    log.score <- log.score / (T-seas[2])
    quad.score <- quad.score / (T-seas[2])
    rps.score <- rps.score / (T-seas[2])
  }# GP2

  list <- list("log.score" = log.score, "quad.score" = quad.score, "rps.score" = rps.score)
}


##covariates
  if (class(coco) == "coco.fit.c") {
    xreg <- cbind(rep(1,nrow(coco$cov)),coco$cov)


    #Poisson 1
    if ( (coco$type == "Poisson") & (coco$order == 1) ){
      betas <- par[-c(1)]
      lambdas <- exp(as.matrix(xreg) %*% betas)

      log.score <- 0
      quad.score <- 0
      rps.score <- 0

      for (t in (seas[1]+1):T) {
        #log score
        par <- c(lambdas[t], coco$par[1], 0)
        log.score <- log.score -log(dGP1(data[t], data[t-seas[1]], par))

        #quadratic score
        norm.p <- 0
        help.stop <- 0

        j <- 0
        while (help.stop < (1-val.num) )  {
          help <- dGP1(j,data[t-seas[1]], par)
          help.stop <- pGP1(j,data[t-seas[1]], par)
          norm.p <- norm.p + help^2
          j <- j+1
        }
        quad.score <- quad.score - 2* dGP1(data[t], data[t-seas[1]], par) + norm.p

        #ranked probability score
        sum <- 0
        help <- 0
        j <- 0
        while (help < (1-val.num) ){
          ind <- 0
          if (j >= data[t]) {
            ind <- 1
          }
          help <- pGP1(j, data[t-seas[1]], par)
          sum <- sum + (help - ind)^2
          j <- j + 1
        }
        rps.score <- rps.score + sum
      }# end t
      log.score <- log.score / (T-seas[1])
      quad.score <- quad.score / (T-seas[1])
      rps.score <- rps.score / (T-seas[1])
    }# Poisson 1


    #GP1

    if ( (coco$type == "GP") & (coco$order == 1) ){
      betas <- par[-c(1,2)]
      lambdas <- exp(as.matrix(xreg) %*% betas)
      log.score <- 0
      quad.score <- 0
      rps.score <- 0

      for (t in (seas[1]+1):T) {
        #log score
        par <- c(lambdas[t], coco$par[1], coco$par[2])
        log.score <- log.score -log(dGP1(data[t], data[t-seas[1]], par))

        #quadratic score
        norm.p <- 0
        help.stop <- 0

        j <- 0
        while (help.stop < (1-val.num) )  {

          help <- dGP1(j,data[t-seas[1]], par)
          help.stop <- pGP1(j,data[t-seas[1]], par)
          norm.p <- norm.p + help^2
          j <- j+1
        }
        quad.score <- quad.score - 2 * dGP1(data[t], data[t-seas[1]], par) + norm.p


        #ranked probability score
        sum <- 0
        help <- 0
        j <- 0
        while (help < (1-val.num) ){
          ind <- 0
          if (j >= data[t]) {
            ind <- 1
          }
          help <- pGP1(j, data[t-seas[1]], par)
          sum <- sum + (help - ind)^2
          j <- j + 1
        }
        rps.score <- rps.score + sum
      }# end t
      log.score <- log.score / (T-seas[1])
      quad.score <- quad.score / (T-seas[1])
      rps.score <- rps.score / (T-seas[1])
    }# GP1

    #Poisson 2
    if ( (coco$type == "Poisson") & (coco$order == 2) ){
      betas <- par[-c(1,2,3)]
      lambdas <- exp(as.matrix(xreg) %*% betas)
      log.score <- 0
      quad.score <- 0
      rps.score <- 0

      for (t in (seas[2]+1):T) {
        #log score
        par <- c(lambdas[t], coco$par[1], coco$par[2], coco$par[3], 0)
        log.score <- log.score -log(dGP2(data[t], data[t-seas[1]], data[t-seas[2]], par))

        #quadratic score
        norm.p <- 0
        help.stop <- 0

        j <- 0
        while (help.stop < (1-val.num) )  {
          help <- dGP2(j,data[t-seas[1]], data[t-seas[2]], par)
          help.stop <- pGP2(j,data[t-seas[1]], data[t-seas[2]], par)
          norm.p <- norm.p + help^2
          j <- j+1
        }
        quad.score <- quad.score - 2* dGP2(data[t], data[t-seas[1]], data[t-seas[2]], par) + norm.p

        #ranked probability score
        sum <- 0
        help <- 0
        j <- 0
        while (help < (1-val.num) ){
          ind <- 0
          if (j >= data[t]) {
            ind <- 1
          }
          help <- pGP2(j, data[t-seas[1]],data[t-seas[2]], par)
          sum <- sum + (help - ind)^2
          j <- j + 1
        }
        rps.score <- rps.score + sum
      }# end t
      log.score <- log.score / (T-seas[2])
      quad.score <- quad.score / (T-seas[2])
      rps.score <- rps.score / (T-seas[2])
    }# Poisson 2

    #GP2
    if ( (coco$type == "GP") & (coco$order == 2) ){
      betas <- par[-c(1,2,3,4)]
      lambdas <- exp(as.matrix(xreg) %*% betas)
      log.score <- 0
      quad.score <- 0
      rps.score <- 0

      for (t in (seas[2]+1):T) {
        #log score
        par <- c(lambdas[t], coco$par[1], coco$par[2], coco$par[3], coco$par[4])
        log.score <- log.score -log(dGP2(data[t], data[t-seas[1]], data[t-seas[2]], par))

        #quadratic score
        norm.p <- 0
        help.stop <- 0

        j <- 0
        while (help.stop < (1-val.num) )  {
          help <- dGP2(j,data[t-seas[1]], data[t-seas[2]], par)
          help.stop <- pGP2(j,data[t-seas[1]], data[t-seas[2]], par)
          norm.p <- norm.p + help^2
          j <- j+1
        }
        quad.score <- quad.score - 2* dGP2(data[t], data[t-seas[1]], data[t-seas[2]], par) + norm.p

        #ranked probability score
        sum <- 0
        help <- 0
        j <- 0
        while (help < (1-val.num) ){
          ind <- 0
          if (j >= data[t]) {
            ind <- 1
          }
          help <- pGP2(j, data[t-seas[1]],data[t-seas[2]], par)
          sum <- sum + (help - ind)^2
          j <- j + 1
        }
        rps.score <- rps.score + sum
      }# end t
      log.score <- log.score / (T-seas[2])
      quad.score <- quad.score / (T-seas[2])
      rps.score <- rps.score / (T-seas[2])
    }# GP2

    list <- list("log.score" = log.score, "quad.score" = quad.score, "rps.score" = rps.score)
  }


  return(list)

}# end function

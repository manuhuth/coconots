cocoReg_base <- function(type, order, data, seasonality = c(1, 2), #mu = 1e-4, outer.it = 500, outer.eps = 1e-10,
                         method_optim=method_optim,
                         optim_control = FALSE, constrained.optim = TRUE, start = NULL,
                         start.val.adjust = TRUE, replace.start.val = 1e-5,
                         iteration.start.val = 0.99, method.hessian = "Richardson", ...) {
  start_time <- Sys.time()

  if (replace.start.val <= 0) {
    stop("Option 'replace.start.val' must be a non-negative real number")
  }

  if ((type != "GP") & (type != "Poisson")) {
    stop("Option 'type' must be either Poisson or GP")
  }

  if ((constrained.optim != TRUE) & (constrained.optim != FALSE)) {
    stop("Option 'constrained.optim' must be either TRUE or FALSE")
  }

  if ((start.val.adjust != TRUE) & (start.val.adjust != FALSE)) {
    stop("Option 'start.val.adjust' must be either TRUE or False")
  }

  if ((type != "GP") & (type != "Poisson")) {
    stop("Option 'type' must be either Poisson or GP")
  }

  if (is.data.frame(data)) {
    data <- as.matrix(data)
  }

  if ((order != 1) & order != (2)) {
    stop("Option 'order' must be 1 or 2")
  }
  
  if (length(seasonality == 1)) {
    seasonality <- c(seasonality, seasonality + 1)
  }

  if (seasonality[1] >= seasonality[2]) {
    stop("The first parameter of 'seasonality' must not be bigger than the second parameter")
  }

  if ((seasonality[1] != round(seasonality[1])) | (seasonality[1] < 1) |
    (seasonality[2] != round(seasonality[2])) | (seasonality[2] < 1)) {
    stop("The values of 'seasonality' must be positive integer values")
  }

  if (isTRUE(optim_control)) {
    optim_control <- list(trace = 1)
  } else {
    optim_control <- list()
  }
  se_boundary <- 0.00001
  lag.max <- max(seasonality) + 1
  unconstrained.optim.lower <- -Inf
  unconstrained.optim.upper <- Inf 

  # PAR1
  if ((type == "Poisson") & (order == 1)) {
    PAR1 <- function(data) {
      mlf <- function(par, data) {
        lambda <- par[1]
        alpha <- par[2]
        eta <- 0

        T <- length(data)

        mlef <- likelihoodGP1(20, lambda, alpha, eta, T, seasonality[1], data)

        return(mlef)
      } # end function
      
      #unconstrained.optim.lower <- c(0,0)
      #unconstrained.optim.upper <- c(Inf, 1) 


      ## starting values
      if (length(start) == 0) {
        acf <- forecast::Acf(data, plot = FALSE, lag.max = lag.max)
        alpha_s <- acf$acf[seasonality[1] + 1]

        if (start.val.adjust == TRUE) {
          if (alpha_s < 0) {
            alpha_s <- replace.start.val
            warning(paste("The initial value of alpha is not in the feasible region and was set to", replace.start.val, ""))
          }
        }
        eta_s <- 0
        lambda_s <- mean(data) * (1 - eta_s) * (1 - alpha_s)

        if (start.val.adjust == TRUE) {
          if (lambda_s < 0) {
            lambda_s <- replace.start.val
            warning(paste("The initial value of lambda is not in the feasible region and was set to", replace.start.val, ""))
          }
        }

        sta <- c(lambda_s, alpha_s)
      } # end compute starting values

      if (length(start) != 0) {
        sta <- start
        if (length(start) != 2) {
          stop("Number of initial starting values must equal 2 for the Poisson 1 model")
        }
      } # end custom starting values


      if (constrained.optim == TRUE) {
        # constrained.optimaints
        ui <- rbind(c(1, 0), c(0, 1), c(0, -1))
        ci <- c(0, 0, -1)

        # constrained.optimained optimization process
        if (method_optim == "BFGS") {
          gradient <- function(par, data){return(numDeriv::grad(func=mlf, x=par, data=data))}
          fit <- stats::constrOptim(
            theta = sta, f = mlf, grad=gradient, ui = ui, ci = ci, data = data, #mu = mu,
            method = method_optim, control = optim_control,
            #outer.iterations = outer.it, outer.eps = outer.eps,
            hessian = FALSE, ...)          
        }else{
          fit <- stats::constrOptim(
            theta = sta, f = mlf, ui = ui, ci = ci, data = data,
            method = method_optim, control = optim_control,
            hessian = FALSE, ...
          )
        }
      } # end constrained.optim True

      if (constrained.optim == FALSE) {
        if ((length(unconstrained.optim.lower) != 2) & (unconstrained.optim.lower != -Inf)) {
          stop("Number of lower bounds must equal 2 for the Poisson 1 model")
        }

        if ((length(unconstrained.optim.upper) != 2) & (unconstrained.optim.upper != Inf)) {
          stop("Number of upper bounds must equal 2 for the Poisson 1 model")
        }
        
        fit <- stats::optim(
          par = sta, fn = mlf, gr = NULL, method = c(method_optim), data = data,
          lower = unconstrained.optim.lower, upper = unconstrained.optim.upper,
          hessian = FALSE, ...
        )
      } # end constrained.optim False

      pars <- as.vector(unlist(fit$par, use.names = FALSE))
      names(pars) <- c("lambda", "alpha")
      T <- length(data)
      likelihood <- -likelihoodGP1(20, pars[1], pars[2], 0, T, seasonality[1], data)

      gra <- numDeriv::grad(func=mlf, x=pars, method = method.hessian, data=data)
      hes <- numDeriv::hessian(func=mlf, x=pars, method = method.hessian, data=data)
      inv_hes <- solve(hes)
      se <- diag(inv_hes)^0.5
      names(se) <- c("lambda", "alpha")

      end_time <- Sys.time()
      if ((pars[2] < se_boundary) |  (pars[2] > 1 - se_boundary)) {
        warning("The estimate of alpha is close to the boundary. Standard errors might not be valid.")
      }

      list_func <- list(
        "par" = pars, "gradient" = gra, "hessian" = hes,
        "inv hessian" = inv_hes, "se" = se, "ts" = data, "type" = "Poisson",
        "seasonality" = seasonality, "order" = 1, "likelihood" = likelihood, "duration" = end_time - start_time,
        julia_reg = NULL
      )

      class(list_func) <- "coco.fit"
      return(list_func)
    } # end PAR1
    return(PAR1(data))
  } # end if PAR1



  # GP1
  if ((type == "GP") & (order == 1)) {
    GP1 <- function(data) {
      mlf <- function(par, data) {
        lambda <- par[1]
        alpha <- par[2]
        eta <- par[3]

        T <- length(data)

        mlef <- likelihoodGP1(20, lambda, alpha, eta, T, seasonality[1], data)

        return(mlef)
      } # end function

      #unconstrained.optim.lower <- c(0,0,0)
      #unconstrained.optim.upper <- c(Inf, 1, 1) 

      ## starting values
      if (length(start) == 0) {
        acf <- forecast::Acf(data, plot = FALSE, lag.max = lag.max)
        alpha_s <- acf$acf[seasonality[1] + 1]

        if (start.val.adjust == TRUE) {
          if (alpha_s < 0) {
            alpha_s <- replace.start.val
            warning(paste("The initial value of alpha is not in the feasible region and was set to", replace.start.val, ""))
          }
        }
        eta_s <- 1 - (mean(data) / stats::var(data))^0.5
        if (start.val.adjust == TRUE) {
          if (eta_s < 0) {
            eta_s <- replace.start.val
            warning(paste("The initial value of eta is not in the feasible region and was set to", replace.start.val, ""))
          }

          if (eta_s > 1) {
            eta_s <- 1 - replace.start.val
            warning(paste("The initial value of eta is not in the feasible region and was set to", 1 - replace.start.val, ""))
          }
        }


        lambda_s <- mean(data) * (1 - eta_s) * (1 - alpha_s)
        if (start.val.adjust == TRUE) {
          if (lambda_s < 0) {
            lambda_s <- replace.start.val
            warning(paste("The initial value of alpha2 is not in the feasible region and was set to", replace.start.val, ""))
          }
        }
        sta <- c(lambda_s, alpha_s, eta_s)
      }

      if (length(start) != 0) {
        sta <- start
        if (length(start) != 3) {
          stop("Number of initial starting values must equal 3 for the GP 1 model")
        }
      }

      if (constrained.optim == TRUE) {
        # constrained.optimaints
        ui <- rbind(c(1, 0, 0), c(0, 1, 0), c(0, -1, 0), c(0, 0, 1), c(0, 0, -1))
        ci <- c(0, 0, -1, 0, -1)

        # constrained.optimained optimization process
        if (method_optim == "BFGS") {
          gradient <- function(par, data){return(numDeriv::grad(func=mlf, x=par, data=data))}
          fit <- stats::constrOptim(
            theta = sta, f = mlf, grad=gradient, ui = ui, ci = ci, data = data,
            method = method_optim, control = optim_control,
            hessian = FALSE, ...)          
        }else{
          fit <- stats::constrOptim(
            theta = sta, f = mlf, ui = ui, ci = ci, data = data,
            method = method_optim, control = optim_control,

            hessian = FALSE, ...
          )
        }
      } # end constrained.optim TRUE

      if (constrained.optim == FALSE) {
        if ((length(unconstrained.optim.lower) != 3) & (unconstrained.optim.lower != -Inf)) {
          stop("Number of lower bounds must equal 3 for the GP 1 model")
        }

        if ((length(unconstrained.optim.upper) != 3) & (unconstrained.optim.upper != Inf)) {
          stop("Number of upper bounds must equal 3 for the GP 1 model")
        }

        fit <- stats::optim(
          par = sta, fn = mlf, gr = NULL, method = c(method_optim), data = data,
          lower = unconstrained.optim.lower, upper = unconstrained.optim.upper,
          hessian = FALSE, ...
        )
      } # end constrained.optim False



      pars <- as.vector(unlist(fit$par, use.names = FALSE))
      names(pars) <- c("lambda", "alpha", "eta")
      T <- length(data)
      likelihood <- -likelihoodGP1(20, pars[1], pars[2], pars[3], T, seasonality[1], data)


      gra <- numDeriv::grad(func=mlf, x=pars, method = method.hessian, data=data)
      hes <- numDeriv::hessian(func=mlf, x=pars, method = method.hessian, data=data)
      inv_hes <- solve(hes)
      se <- diag(inv_hes)^0.5
      names(se) <- c("lambda", "alpha", "eta")

      end_time <- Sys.time()
      
      if ((pars[2] < se_boundary) |  (pars[2] > 1 - se_boundary)) {
        warning("The estimate of alpha is close to the boundary. Standard errors might not be valid.")
      }
      if ((pars[3] < se_boundary) |  (pars[3] > 1 - se_boundary)) {
        warning("The estimate of eta is close to the boundary. Standard errors might not be valid.")
      }
      
      list_func <- list(
        "par" = pars, "gradient" = gra, "hessian" = hes,
        "inv hessian" = inv_hes, "se" = se, "ts" = data, "type" = "GP",
        "order" = 1, "seasonality" = seasonality, "likelihood" = likelihood, "duration" = end_time - start_time,
        julia_reg = NULL
      )
      class(list_func) <- "coco.fit"
      return(list_func)
    } # end GP1
    return(GP1(data))
  } # end if GP1


  ############ PAR2
  if ((type == "Poisson") & (order == 2)) {
    PAR2 <- function(data) {
      mlf <- function(par, data) {
        lambda <- par[1]
        alpha1 <- par[2]
        alpha2 <- par[3]
        alpha3 <- par[4]
        eta <- 0

        T <- length(data)

        mlef <- likelihoodGP2(20, lambda, alpha1, alpha2, alpha3, eta, T, seasonality[1], seasonality[2], data)

        return(mlef)
      } # end function
      
      #unconstrained.optim.lower <- c(0,0,0,0)
      #unconstrained.optim.upper <- c(Inf, 1, 1,1) 

      ## starting values

      if (length(start) == 0) {
        acf <- forecast::Acf(data, plot = FALSE, lag.max = lag.max)
        p1 <- acf$acf[seasonality[1] + 1]
        p2 <- acf$acf[seasonality[2] + 1]

        v <- c()
        for (t in (seasonality[2] + 1):length(data)) {
          v[t - seasonality[2]] <- (data[t] - mean(data)) * (data[t - seasonality[1]] - mean(data)) * (data[t - seasonality[2]] - mean(data))
        }
        mu111 <- sum(v) / length(data)


        eta_s <- 0

        alpha3_s <- mu111 * (1 - eta_s)^4 / (mean(data) * (1 + 2 * eta_s))

        if (start.val.adjust == TRUE) {
          if (alpha3_s < 0) {
            alpha3_s <- replace.start.val
            warning(paste("The initial value of alpha3 is not in the feasible region and was set to", replace.start.val, ""))
          }
          if (alpha3_s > 1) {
            alpha3_s <- 1 - replace.start.val
            warning(paste("The initial value of alpha3 is not in the feasible region and was set to", 1 - replace.start.val, ""))
          }
        }
        alpha2_s <- p2 - alpha3_s
        if (start.val.adjust == TRUE) {
          if (alpha2_s < 0) {
            alpha2_s <- replace.start.val
            warning(paste("The initial value of alpha2 is not in the feasible region and was set to", replace.start.val, ""))
          }
        }
        alpha1_s <- p1 - alpha3_s
        if (start.val.adjust == TRUE) {
          if (alpha1_s < 0) {
            alpha1_s <- replace.start.val
            warning(paste("The initial value of alpha1 is not in the feasible region and was set to", replace.start.val, ""))
          }

          if (alpha1_s + alpha2_s + alpha3_s > 1) {
            while (alpha1_s + alpha2_s + alpha3_s > 1) {
              alpha1_s <- alpha1_s * iteration.start.val
              alpha2_s <- alpha2_s * iteration.start.val
              alpha3_2 <- alpha3_s * iteration.start.val
            }
            warning(paste("The initial sum of 2*alpha1, alpha2 and alpha3 had not been in the feasible region and was adjusted", replace.start.val, ""))
          }


          if (2 * alpha1_s + alpha3_s > 1) {
            while (2 * alpha1_s + alpha3_s > 1) {
              alpha1_s <- alpha1_s * iteration.start.val
              alpha3_2 <- alpha3_s * iteration.start.val
            }
            warning(paste("The initial sum of alpha1 and alpha3 had not been in the feasible region and was adjusted", replace.start.val, ""))
          }
        }

        lambda_s <- mean(data) * (1 - eta_s) * (1 - alpha1_s - alpha2_s - alpha3_s)

        if (start.val.adjust == TRUE) {
          if (lambda_s < 0) {
            lambda_s <- replace.start.val
            warning(paste("The initial value of lambda is not in the feasible region and was set to", replace.start.val, ""))
          }
        }
        sta <- c(lambda_s, alpha1_s, alpha2_s, alpha3_s)
      } # end starting values


      if (length(start) != 0) {
        sta <- start
        if (length(start) != 4) {
          stop("Number of initial starting values must equal 4 for the Poisson 2 model")
        }
      }

      # constrained.optimaints
      if (constrained.optim == TRUE) {
        ui <- rbind(c(0, -1, -1, -1), c(0, -2, 0, -1), c(0, 1, 0, 0), c(0, 0, 1, 0), c(0, 0, 0, 1), c(1, 0, 0, 0))

        ci <- c(-1, -1, 0, 0, 0, 0)

        # constrained.optimained optimization process
        if (method_optim == "BFGS") {
          gradient <- function(par, data){return(numDeriv::grad(func=mlf, x=par, data=data))}
          fit <- stats::constrOptim(
            theta = sta, f = mlf, grad=gradient, ui = ui, ci = ci, data = data, 
            method = method_optim, control = optim_control,
            hessian = FALSE, ...)          
        }else{
          fit <- stats::constrOptim(
            theta = sta, f = mlf, ui = ui, ci = ci, data = data,
            method = method_optim, control = optim_control,
            hessian = FALSE, ...
          )
        }
      } # end constrained.optim True

      if (constrained.optim == FALSE) {
        if ((length(unconstrained.optim.lower) != 4) & (unconstrained.optim.lower != -Inf)) {
          stop("Number of lower bounds must equal 4 for the Poisson 4 model")
        }

        if ((length(unconstrained.optim.upper) != 4) & (unconstrained.optim.upper != Inf)) {
          stop("Number of upper bounds must equal 2 for the Poisson 4 model")
        }

        fit <- stats::optim(
          par = sta, fn = mlf, gr = NULL, method = c(method_optim), data = data,
          lower = unconstrained.optim.lower, upper = unconstrained.optim.upper,
          hessian = FALSE, ...
        )
      } # end constrained.optim False

      pars <- as.vector(unlist(fit$par, use.names = FALSE))
      names(pars) <- c("lambda", "alpha1", "alpha2", "alpha3")
      T <- length(data)
      likelihood <- -likelihoodGP2(20, pars[1], pars[2], pars[3], pars[4], 0, T, seasonality[1], seasonality[2], data)


      gra <- numDeriv::grad(func=mlf, x=pars, method = method.hessian, data=data)
      hes <- numDeriv::hessian(func=mlf, x=pars, method = method.hessian, data=data)
      inv_hes <- solve(hes)
      se <- diag(inv_hes)^0.5
      names(se) <- c("lambda", "alpha1", "alpha2", "alpha3")

      end_time <- Sys.time()
      
      if ((pars[2] < se_boundary) |  (pars[2] > 1 - se_boundary)) {
        warning("The estimate of alpha1 is close to the boundary. Standard errors might not be valid.")
      }
      if ((pars[3] < se_boundary) |  (pars[3] > 1 - se_boundary)) {
        warning("The estimate of alpha2 is close to the boundary. Standard errors might not be valid.")
      }
      if ((pars[4] < se_boundary) |  (pars[4] > 1 - se_boundary)) {
        warning("The estimate of alpha3 is close to the boundary. Standard errors might not be valid.")
      }
    
      
      list_func <- list(
        "par" = pars,
        "gradient" = gra, "hessian" = hes, "inv hessian" = inv_hes,
        "se" = se, "ts" = data, "type" = "Poisson", "order" = 2,
        "seasonality" = seasonality, "likelihood" = likelihood, "duration" = end_time - start_time,
        julia_reg = NULL
      )
      class(list_func) <- "coco.fit"
      return(list_func)
    } # end PAR2
    return(PAR2(data))
  } # end if PAR2


  ############ GP2
  if ((type == "GP") & (order == 2)) {
    GP2 <- function(data) {
      mlf <- function(par, data) {
        lambda <- par[1]
        alpha1 <- par[2]
        alpha2 <- par[3]
        alpha3 <- par[4]
        eta <- par[5]

        T <- length(data)

        mlef <- likelihoodGP2(20, lambda, alpha1, alpha2, alpha3, eta, T, seasonality[1], seasonality[2], data)

        return(mlef)
      } # end function
      
      #unconstrained.optim.lower <- c(0,0,0,0,0)
      #unconstrained.optim.upper <- c(Inf, 1, 1,1,1) 

      ## starting values
      if (length(start) == 0) {
        acf <- forecast::Acf(data, plot = FALSE, lag.max = lag.max)
        p1 <- acf$acf[[seasonality[1] + 1]]
        p2 <- acf$acf[[seasonality[2] + 1]]

        v <- c()
        for (t in (seasonality[2] + 1):length(data)) {
          v[t - seasonality[2]] <- (data[t] - mean(data)) * (data[t - seasonality[1]] - mean(data)) * (data[t - seasonality[2]] - mean(data))
        }
        mu111 <- sum(v) / length(data)

        eta_s <- 1 - (mean(data) / stats::var(data))^0.5

        if (start.val.adjust == TRUE) {
          if (eta_s < 0) {
            eta_s <- replace.start.val
            warning(paste("The initial value of eta is not in the feasible region and was set to", replace.start.val, ""))
          }
          if (eta_s > 1) {
            eta_s <- 1 - replace.start.val
            warning(paste("The initial value of eta is not in the feasible region and was set to", 1 - replace.start.val, ""))
          }
        }
        alpha3_s <- mu111 * (1 - eta_s)^4 / (mean(data) * (1 + 2 * eta_s))
        if (start.val.adjust == TRUE) {
          if (alpha3_s < 0) {
            alpha3_s <- replace.start.val
            warning(paste("The initial value of alpha3 is not in the feasible region and was set to", replace.start.val, ""))
          }
          if (alpha3_s > 1) {
            alpha3_s <- 1 - replace.start.val
            warning(paste("The initial value of alpha3 is not in the feasible region and was set to", 1 - replace.start.val, ""))
          }
        }
        alpha2_s <- p2 - alpha3_s
        if (start.val.adjust == TRUE) {
          if (alpha2_s < 0) {
            alpha2_s <- replace.start.val
            warning(paste("The initial value of alpha2 is not in the feasible region and was set to", replace.start.val, ""))
          }


          if (alpha2_s > 1) {
            alpha2_s <- 1 - replace.start.val
            warning(paste("The initial value of alpha2 is not in the feasible region and was set to", 1 - replace.start.val, ""))
          }
        }

        alpha1_s <- p1 - alpha3_s
        if (start.val.adjust == TRUE) {
          if (alpha1_s < 0) {
            alpha1_s <- replace.start.val
            warning(paste("The initial value of alpha1 is not in the feasible region and was set to", replace.start.val, ""))
          }


          if (alpha1_s + alpha2_s + alpha3_s > 1) {
            while (alpha1_s + alpha2_s + alpha3_s > 1) {
              alpha1_s <- alpha1_s * 0.99
              alpha2_s <- alpha2_s * 0.99
              alpha3_2 <- alpha3_s * 0.99
            }
            warning(paste("The initial sum of 2*alpha1, alpha2 and alpha3 had not been in the feasible region and was adjusted", ""))
          }


          if (2 * alpha1_s + alpha3_s > 1) {
            while (2 * alpha1_s + alpha3_s > 1) {
              alpha1_s <- alpha1_s * 0.99
              alpha3_2 <- alpha3_s * 0.99
            }
            warning(paste("The initial sum of alpha1 and alpha3 had not been in the feasible region and was adjusted", ""))
          }
        }

        lambda_s <- mean(data) * (1 - eta_s) * (1 - alpha1_s - alpha2_s - alpha3_s)
        if (start.val.adjust == TRUE) {
          if (lambda_s < 0) {
            lambda_s <- replace.start.val
            warning(paste("The initial value of lambda is not in the feasible region and was set to", replace.start.val, ""))
          }
        }
        sta <- c(lambda_s, alpha1_s, alpha2_s, alpha3_s, eta_s)
      }

      if (length(start) != 0) {
        sta <- start
        if (length(start) != 5) {
          stop("Number of initial starting values must equal 5 for the GP 2 model")
        }
      }

      if (constrained.optim == TRUE) {
        # constrained.optimaints
        ui <- rbind(
          c(0, -1, -1, -1, 0), c(0, -2, 0, -1, 0), c(0, 1, 0, 0, 0), c(0, 0, 1, 0, 0), c(0, 0, 0, 1, 0), c(1, 0, 0, 0, 0),
          c(0, 0, 0, 0, 1), c(0, 0, 0, 0, -1)
        )
        ci <- c(-1, -1, 0, 0, 0, 0, 0, -1)

        # constrained.optimained optimization process
        if (method_optim == "BFGS") {
          gradient <- function(par, data){return(numDeriv::grad(func=mlf, x=par, data=data))}
          fit <- stats::constrOptim(
            theta = sta, f = mlf, grad=gradient, ui = ui, ci = ci, data = data, 
            method = method_optim, control = optim_control,
            hessian = FALSE, ...)          
        }else{
          fit <- stats::constrOptim(
            theta = sta, f = mlf, ui = ui, ci = ci, data = data,
            method = method_optim, control = optim_control,
            hessian = FALSE, ...
          )
        }
      } # end constrained.optim True

      if (constrained.optim == FALSE) {
        if ((length(unconstrained.optim.lower) != 5) & (unconstrained.optim.lower != -Inf)) {
          stop("Number of lower bounds must equal 5 for the GP 2 model")
        }

        if ((length(unconstrained.optim.upper) != 5) & (unconstrained.optim.upper != Inf)) {
          stop("Number of upper bounds must equal 5 for the GP 2 model")
        }

        fit <- stats::optim(
          par = sta, fn = mlf, gr = NULL, method = c(method_optim), data = data,
          lower = unconstrained.optim.lower, upper = unconstrained.optim.upper,
          hessian = FALSE, ...
        )
      } # end constrained.optim False


      pars <- as.vector(unlist(fit$par, use.names = FALSE))
      names(pars) <- c("lambda", "alpha1", "alpha2", "alpha3", "eta")
      T <- length(data)
      likelihood <- -likelihoodGP2(20, pars[1], pars[2], pars[3], pars[4], pars[5], T, seasonality[1], seasonality[2], data)


      gra <- numDeriv::grad(func=mlf, x=pars, method = method.hessian, data=data)
      hes <- numDeriv::hessian(func=mlf, x=pars, method = method.hessian, data=data)
      inv_hes <- solve(hes)
      se <- diag(inv_hes)^0.5
      names(se) <- c("lambda", "alpha1", "alpha2", "alpha3", "eta")
      end_time <- Sys.time()
      se_boundary <- 0.0001

      if ((pars[2] < se_boundary) |  (pars[2] > 1 - se_boundary)) {
        warning("The estimate of alpha1 is close to the boundary. Standard errors might not be valid.")
      }
      if ((pars[3] < se_boundary) |  (pars[3] > 1 - se_boundary)) {
        warning("The estimate of alpha2 is close to the boundary. Standard errors might not be valid.")
      }
      if ((pars[4] < se_boundary) |  (pars[4] > 1 - se_boundary)) {
        warning("The estimate of alpha3 is close to the boundary. Standard errors might not be valid.")
      }
      
      if ((pars[5] < se_boundary) |  (pars[5] > 1 - se_boundary)) {
        warning("The estimate of eta is close to the boundary. Standard errors might not be valid.")
      }
    
      list_func <- list(
        "par" = pars, "gradient" = gra,
        "hessian" = hes, "inv hessian" = inv_hes, "se" = se, "ts" = data,
        "type" = "GP", "order" = 2, "seasonality" = seasonality, "likelihood" = likelihood,
        "duration" = end_time - start_time, julia_reg = NULL
      )
      class(list_func) <- "coco.fit"
      return(list_func)
    } # end GP2
    return(GP2(data))
  } # end if GP2
}

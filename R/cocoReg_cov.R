#' @importFrom forecast Acf
#' @importFrom stats glm lm constrOptim optim var na.omit
#' @importFrom JuliaConnectoR juliaCall
cocoReg_cov <- function(type, order, data, xreg, seasonality = c(1, 2), mu = 1e-4,
                        outer.it = 500, outer.eps = 1e-10, optim_control = FALSE, constrained.optim = TRUE, b.beta = -10,
                        start = NULL, start.val.adjust = TRUE, method_optim= "Nelder-Mead",
                        replace.start.val = 1e-5, iteration.start.val = 0.99,
                        method.hessian = NULL, julia_installed=FALSE, link_function="log", cores = 1, ...) {
  start_time <- Sys.time()

  if (is.data.frame(data)) {
    data <- as.matrix(data)
  }

  if (is.data.frame(xreg)) {
    xreg <- as.matrix(xreg)
  }

  if ((outer.it <= 0) | (outer.it != round(outer.it))) {
    stop("Option 'outer.it' must be a positive integer value")
  }

  if (outer.eps <= 0) {
    stop("Option 'outer.eps' must be a non-negative real number")
  }

  if (mu <= 0) {
    stop("Option 'mu' must be a non-negative real number")
  }

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
  ncov <- ncol(xreg)
  lag.max <- max(seasonality) + 1
  unconstrained.optim.lower <- -Inf
  unconstrained.optim.upper <- Inf
  
  df_covariates <- as.data.frame(cbind(data,xreg))
  
  if (link_function=="log"){
    covariates_glm_fit <- stats::glm(data ~ -1 + ., family="poisson", data=df_covariates)
    starting_values_covariates <- stats::na.omit(covariates_glm_fit$coefficients)
  } else if (link_function %in% c("identity", "relu", "softplus")) {
    covariates_glm_fit <- stats::lm(data ~ -1 + ., data=df_covariates)
    starting_values_covariates <- stats::na.omit(covariates_glm_fit$coefficients)
    starting_values_covariates[starting_values_covariates <= 0] <- 1e-5
  } else {
    stop("Choose either identity, log, relu, or softplus as link function.")
  }

  # PAR1
  if ((type == "Poisson") & (order == 1)) {
    PAR1 <- function(data) {
      xreg_matrix <- data.matrix(xreg)
      mlf <- function(par, data) {
        alpha <- par[1]
        eta <- 0
        vec_lambda <- par[-1]

        T <- length(data)
        mlef <- likelihoodGP1cov(alpha, eta, vec_lambda, T, seasonality[1], data, xreg_matrix, link_function)
        return(mlef)
      } # end function


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

        sta <- c(alpha_s, covariateStartingValues(data, xreg, link_function, 1 - alpha_s,
                                                  starting_values_covariates, b.beta))

      }

      if (length(start) != 0) {
        sta <- start
        if (length(start) != 1 + ncol(xreg)) {
          stop(paste("Number of initial starting values must equal", 1 + ncol(xreg), "for the Poisson 1 model with", ncol(xreg), "covariates"))
        }
      }

      # derivatives of mlf by automatic differentiation (C++ dual numbers)
      mlf_derivs <- function(par, data, deriv = 2L) {
        coco_ad_derivs(par, data, type, 1, seasonality, xreg_matrix, link_function, deriv)
      }
      gradient <- function(par, data) mlf_derivs(par, data, 1L)$gradient

      derivs_opt <- NULL
      if (constrained.optim == TRUE) {
        # constrained.optimaints
        ui <- rbind(c(1), c(-1))

        matrix1 <- matrix(0, 2, ncol(xreg) ) # for additional columns

        matrix2 <- matrix(0, ncol(xreg) , 1) # for additional rows
        einheitsmatrix <- diag(ncol(xreg) )
        lower_part <- cbind(matrix2, einheitsmatrix)

        ui <- cbind(ui, matrix1) # right part of ne matrix
        ui <- rbind(ui, lower_part)

        ci <- c(0, -1, rep(b.beta, ncol(xreg)))

        # constrained.optimained optimization process
        if (method_optim == "BFGS") {
          fit <- stats::constrOptim(
            theta = sta, f = mlf, grad=gradient, ui = ui, ci = ci, data = data, #mu = mu,
            method = method_optim, control = optim_control,
            #outer.iterations = outer.it, outer.eps = outer.eps,
            hessian = FALSE)          
          fit <- coco_polish_fit(fit, mlf, ui, ci, gr = gradient, data = data)
        }else{
          constr_fit <- function() {
            fit <- coco_constr_nm_fit(
              sta, mlf, gradient, ui, ci, mu = mu, outer.it = outer.it,
              outer.eps = outer.eps, control = optim_control, data = data
            )
            coco_polish_fit(fit, mlf, ui, ci, gr = gradient, data = data)
          }
          if (method_optim == "Nelder-Mead") {
            # L-BFGS-B with AD gradient from the starting values, accepted
            # after an AD optimality check; else the constrOptim fit above
            opt <- coco_fit_lbfgsb_checked(
              sta, mlf, gradient, function(p) mlf_derivs(p, data, 2L), ui, ci,
              fallback = constr_fit, data = data
            )
            fit <- opt$fit
            derivs_opt <- opt$derivs
          } else {
            fit <- constr_fit()
          }
        }
      } # end constrained.optim True

      if (constrained.optim == FALSE) {
        if (length(unconstrained.optim.lower) != 1 + ncol(xreg) & (unconstrained.optim.lower != -Inf)) {
          stop(paste("Number of lower bounds must equal", 1 + ncol(xreg), "for the Poisson 1 model with", ncol(xreg), "covariates"))
        }

        if (length(unconstrained.optim.upper) != 1 + ncol(xreg) & (unconstrained.optim.upper != Inf)) {
          stop(paste("Number of upper bounds must equal", 1 + ncol(xreg), "for the Poisson 1 model with", ncol(xreg), "covariates"))
        }

        fit <- stats::optim(
          par = sta, fn = mlf, gr = NULL, method = c("Nelder-Mead"), data = data,
          lower = unconstrained.optim.lower, upper = unconstrained.optim.upper,
          control = optim_control, hessian = FALSE
        )
      } # end constrained.optim False

      pars <- as.vector(unlist(fit$par, use.names = FALSE))
      lambs <- paste0("beta", 1:ncov)
      names(pars) <- c("alpha", lambs)
      T <- length(data)
      likelihood <- -likelihoodGP1cov(pars[1], 0, pars[-1], T, seasonality[1], data, xreg_matrix, link_function)

      derivs <- if (!is.null(derivs_opt)) derivs_opt else mlf_derivs(pars, data, 2L)
      gra <- derivs$gradient
      hes <- derivs$hessian
      inv_hes <- solve(hes)
      se <- diag(inv_hes)^0.5
      names(se) <- c("alpha", lambs)

      end_time <- Sys.time()
      if ((pars[1] < se_boundary) |  (pars[1] > 1 - se_boundary)) {
        warning("The estimate of alpha is close to the boundary. Standard errors might not be valid.")
      }
      
      if (julia_installed) {
        addJuliaFunctions()
        julia_reg <- JuliaConnectoR::juliaCall("Coconots.create_julia_dict", 
                                               list("parameter", "covariance_matrix", "log_likelihood",
                                                    "type", "order", "data", "covariates",
                                                    "link", "starting_values", "optimizer", 
                                                    "lower_bounds", "upper_bounds", "optimization", 
                                                    "max_loop"),
                                               list(pars, inv_hes, likelihood,
                                                    type, order, data, xreg,
                                                    link_function, NULL, NULL, NULL, NULL,
                                                    NULL, NULL))
      } else {julia_reg = NULL}
      
      
      list_func <- list(
        "par" = pars, "gradient" = gra,
        "hessian" = hes, "inv hessian" = inv_hes, "se" = se,
        "ts" = data, "cov" = xreg, "type" = "Poisson", "order" = 1,
        "seasonality" = seasonality, "likelihood" = likelihood,
        "duration" = end_time - start_time, julia_reg = julia_reg, "link_function"=link_function
      )
      
      return(list_func)
    } # end PAR1
    return(PAR1(data))
  } # end if PAR1



  #GP1
  if ((type == "GP") & (order == 1)) {
    GP1 <- function(data) {
      xreg_matrix <- data.matrix(xreg)
      mlf <- function(par, data) {
        alpha <- par[1]
        eta <- par[2]
        vec_lambda <- par[-(1:2)]

        T <- length(data)
        mlef <- likelihoodGP1cov(alpha, eta, vec_lambda, T, seasonality[1], data, xreg_matrix, link_function)
        return(mlef)
      } # end function


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

        sta <- c(alpha_s, eta_s, covariateStartingValues(data, xreg, link_function,
                                                         (1 - alpha_s) * (1 - eta_s),
                                                         starting_values_covariates, b.beta))

      }

      if (length(start) != 0) {
        sta <- start
        if (length(start) != 2 + ncol(xreg)) {
          stop(paste("Number of initial starting values must equal", 2 + ncol(xreg), "for the GP 1 model with", ncol(xreg), "covariates"))
        }
      }

      # derivatives of mlf by automatic differentiation (C++ dual numbers)
      mlf_derivs <- function(par, data, deriv = 2L) {
        coco_ad_derivs(par, data, type, 1, seasonality, xreg_matrix, link_function, deriv)
      }
      gradient <- function(par, data) mlf_derivs(par, data, 1L)$gradient

      derivs_opt <- NULL
      if (constrained.optim == TRUE) {
        # constrained.optimaints
        ui <- rbind(c(1, 0), c(-1, 0), c(0, 1), c(0, -1))

        matrix1 <- matrix(0, 4, ncol(xreg) ) # for additional columns

        matrix2 <- matrix(0, ncol(xreg) , 2) # for additional rows
        einheitsmatrix <- diag(ncol(xreg) )
        lower_part <- cbind(matrix2, einheitsmatrix)

        ui <- cbind(ui, matrix1) # right part of the matrix
        ui <- rbind(ui, lower_part)

        ci <- c(0, -1, 0, -1, rep(b.beta, ncol(xreg) ))

        # constrained.optimained optimization process
        if (method_optim == "BFGS") {
          fit <- stats::constrOptim(
            theta = sta, f = mlf, grad=gradient, ui = ui, ci = ci, data = data, #mu = mu,
            method = method_optim, control = optim_control,
            #outer.iterations = outer.it, outer.eps = outer.eps,
            hessian = FALSE)          
          fit <- coco_polish_fit(fit, mlf, ui, ci, gr = gradient, data = data)
        }else{
          constr_fit <- function() {
            fit <- coco_constr_nm_fit(
              sta, mlf, gradient, ui, ci, mu = mu, outer.it = outer.it,
              outer.eps = outer.eps, control = optim_control, data = data
            )
            coco_polish_fit(fit, mlf, ui, ci, gr = gradient, data = data)
          }
          if (method_optim == "Nelder-Mead") {
            # L-BFGS-B with AD gradient from the starting values, accepted
            # after an AD optimality check; else the constrOptim fit above
            opt <- coco_fit_lbfgsb_checked(
              sta, mlf, gradient, function(p) mlf_derivs(p, data, 2L), ui, ci,
              fallback = constr_fit, data = data
            )
            fit <- opt$fit
            derivs_opt <- opt$derivs
          } else {
            fit <- constr_fit()
          }
        }
      } # constrained.optim True

      if (constrained.optim == FALSE) {
        if (length(unconstrained.optim.lower) != 2 + ncol(xreg) & (unconstrained.optim.lower != -Inf)) {
          stop(paste("Number of lower bounds must equal", 2 + ncol(xreg), "for the GP 1 model with", ncol(xreg), "covariates"))
        }

        if (length(unconstrained.optim.upper) != 2 + ncol(xreg) & (unconstrained.optim.upper != Inf)) {
          stop(paste("Number of upper bounds must equal", 2 + ncol(xreg), "for the GP 1 model with", ncol(xreg), "covariates"))
        }

        fit <- stats::optim(
          par = sta, fn = mlf, gr = NULL, method = c("Nelder-Mead"), data = data,
          lower = unconstrained.optim.lower, upper = unconstrained.optim.upper,
          control = optim_control, hessian = FALSE
        )
      } # end constrained.optim False

      pars <- as.vector(unlist(fit$par, use.names = FALSE))
      lambs <- paste0("beta", 1:ncov)
      names(pars) <- c("alpha", "eta", lambs)
      T <- length(data)
      likelihood <- -likelihoodGP1cov(pars[1], pars[2], pars[-c(1, 2)], T, seasonality[1], data, xreg_matrix, link_function)

      derivs <- if (!is.null(derivs_opt)) derivs_opt else mlf_derivs(pars, data, 2L)
      gra <- derivs$gradient
      hes <- derivs$hessian
      inv_hes <- solve(hes)
      se <- diag(inv_hes)^0.5
      names(se) <- c("alpha", "eta", lambs)

      end_time <- Sys.time()
      
      if ((pars[1] < se_boundary) |  (pars[1] > 1 - se_boundary)) {
        warning("The estimate of alpha is close to the boundary. Standard errors might not be valid.")
      }
      if ((pars[2] < se_boundary) |  (pars[2] > 1 - se_boundary)) {
        warning("The estimate of eta is close to the boundary. Standard errors might not be valid.")
      }
      
      if (julia_installed) {
        addJuliaFunctions()
        julia_reg <- JuliaConnectoR::juliaCall("Coconots.create_julia_dict", 
                                               list("parameter", "covariance_matrix", "log_likelihood",
                                                    "type", "order", "data", "covariates",
                                                    "link", "starting_values", "optimizer", 
                                                    "lower_bounds", "upper_bounds", "optimization", 
                                                    "max_loop"),
                                               list(pars, inv_hes, likelihood,
                                                    type, order, data, xreg,
                                                    link_function, NULL, NULL, NULL, NULL,
                                                    NULL, NULL))
      } else {julia_reg = NULL}
      
      list_func <- list(
        "par" = pars, "gradient" = gra,
        "hessian" = hes, "inv hessian" = inv_hes, "se" = se,
        "ts" = data, "cov" = xreg, "type" = "GP", "order" = 1,
        "seasonality" = seasonality, "likelihood" = likelihood, "duration" = end_time - start_time,
        julia_reg = julia_reg, "link_function"=link_function
      )
 
      return(list_func)
    } # end PGP1
    return(GP1(data))
  } # end if GP1


  #PAR2
  if ((type == "Poisson") & (order == 2)) {
    PAR2 <- function(data) {
      xreg_matrix <- data.matrix(xreg)
      mlf <- function(par, data) {
        alpha1 <- par[1]
        alpha2 <- par[2]
        alpha3 <- par[3]
        eta <- 0
        vec_lambda <- par[-(1:3)]

        T <- length(data)
        mlef <- likelihoodGP2cov(alpha1, alpha2, alpha3, eta, vec_lambda, T, seasonality[1], seasonality[2], data, xreg_matrix, link_function, cores)

        return(mlef)
      } # end function

      ## starting values
      if (length(start) == 0) {
        acf <- forecast::Acf(data, plot = FALSE, lag.max = lag.max)
        p1 <- acf$acf[seasonality[1] + 1]
        p2 <- acf$acf[seasonality[2] + 1]

        mean_data <- mean(data)
        idx <- (seasonality[2] + 1):length(data)
        v <- (data[idx] - mean_data) * (data[idx - seasonality[1]] - mean_data) * (data[idx - seasonality[2]] - mean_data)
        mu111 <- sum(v) / length(data)


        eta_s <- 0
        alpha3_s <- mu111 * (1 - eta_s)^4 / (mean(data) * (1 + 2 * eta_s))

        if (start.val.adjust == TRUE) {
          if (alpha3_s < 0) {
            alpha3_s <- replace.start.val
            warning(paste("The initial value of alpha3 is not in the feasible region and was set to", replace.start.val, ""))
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
              alpha3_s <- alpha3_s * iteration.start.val
            }
            warning(paste("The initial sum of 2*alpha1, alpha2 and alpha3 had not been in the feasible region and was adjusted", replace.start.val, ""))
          }


          if (2 * alpha1_s + alpha3_s > 1) {
            while (2 * alpha1_s + alpha3_s > 1) {
              alpha1_s <- alpha1_s * iteration.start.val
              alpha3_s <- alpha3_s * iteration.start.val
            }
            warning(paste("The initial sum of alpha1 and alpha3 had not been in the feasible region and was adjusted", replace.start.val, ""))
          }
        }
        #starting values covariates

        sta <- c(alpha1_s, alpha2_s, alpha3_s,
                 covariateStartingValues(data, xreg, link_function,
                                         1 - alpha1_s - alpha2_s - alpha3_s,
                                         starting_values_covariates, b.beta))
      }


      if (length(start) != 0) {
        sta <- start
        if (length(start) != 3 + ncol(xreg)) {
          stop(paste("Number of initial starting values must equal", 3 + ncol(xreg), "for the Poisson 2 model with", ncol(xreg), "covariates"))
        }
      }

      # derivatives of mlf by automatic differentiation (C++ dual numbers)
      mlf_derivs <- function(par, data, deriv = 2L) {
        coco_ad_derivs(par, data, type, 2, seasonality, xreg_matrix, link_function, deriv)
      }
      gradient <- function(par, data) mlf_derivs(par, data, 1L)$gradient

      derivs_opt <- NULL
      if (constrained.optim == TRUE) {
        # constrained.optimaints
        ui <- rbind(c(-1, -1, -1), c(-2, 0, -1), c(1, 0, 0), c(0, 1, 0), c(0, 0, 1))

        matrix1 <- matrix(0, 5, ncol(xreg) ) # for additional columns

        matrix2 <- matrix(0, ncol(xreg) , 3) # for additional rows
        einheitsmatrix <- diag(ncol(xreg) )
        lower_part <- cbind(matrix2, einheitsmatrix)

        ui <- cbind(ui, matrix1) # right part of ne matrix
        ui <- rbind(ui, lower_part)

        ci <- c(-1, -1, 0, 0, 0, rep(b.beta, ncol(xreg) ))

        # constrained.optimained optimization process
        if (method_optim == "BFGS") {
          fit <- stats::constrOptim(
            theta = sta, f = mlf, grad=gradient, ui = ui, ci = ci, data = data, #mu = mu,
            method = method_optim, control = optim_control,
            #outer.iterations = outer.it, outer.eps = outer.eps,
            hessian = FALSE)          
          fit <- coco_polish_fit(fit, mlf, ui, ci, gr = gradient, data = data)
        }else{
          constr_fit <- function() {
            fit <- coco_constr_nm_fit(
              sta, mlf, gradient, ui, ci, mu = mu, outer.it = outer.it,
              outer.eps = outer.eps, control = optim_control, data = data
            )
            coco_polish_fit(fit, mlf, ui, ci, gr = gradient, data = data)
          }
          if (method_optim == "Nelder-Mead") {
            # L-BFGS-B with AD gradient from the starting values, accepted
            # after an AD optimality check; else the constrOptim fit above
            opt <- coco_fit_lbfgsb_checked(
              sta, mlf, gradient, function(p) mlf_derivs(p, data, 2L), ui, ci,
              fallback = constr_fit, data = data
            )
            fit <- opt$fit
            derivs_opt <- opt$derivs
          } else {
            fit <- constr_fit()
          }
        }
      } # end constrained.optim

      if (constrained.optim == FALSE) {
        if (length(unconstrained.optim.lower) != 3 + ncol(xreg) & (unconstrained.optim.lower != -Inf)) {
          stop(paste("Number of lower bounds must equal", 4 + ncol(xreg), "for the Poisson 2 model with", ncol(xreg), "covariates"))
        }

        if (length(unconstrained.optim.upper) != 3 + ncol(xreg) & (unconstrained.optim.upper != Inf)) {
          stop(paste("Number of upper bounds must equal", 3 + ncol(xreg), "for the Poisson 2 model with", ncol(xreg), "covariates"))
        }

        fit <- stats::optim(
          par = sta, fn = mlf, gr = NULL, method = c("Nelder-Mead"), data = data,
          lower = unconstrained.optim.lower, upper = unconstrained.optim.upper,
          control = optim_control, hessian = FALSE
        )
      } # end constrained.optim False

      pars <- as.vector(unlist(fit$par, use.names = FALSE))
      lambs <- paste0("beta", 1:ncov)
      names(pars) <- c("alpha1", "alpha2", "alpha3", lambs)
      T <- length(data)
      likelihood <- -likelihoodGP2cov(pars[1], pars[2], pars[3], 0, pars[-c(1:3)], T, seasonality[1], seasonality[2], data, xreg_matrix, link_function)


      derivs <- if (!is.null(derivs_opt)) derivs_opt else mlf_derivs(pars, data, 2L)
      gra <- derivs$gradient
      hes <- derivs$hessian
      inv_hes <- solve(hes)
      se <- diag(inv_hes)^0.5
      names(se) <- c("alpha1", "alpha2", "alpha3", lambs)

      end_time <- Sys.time()
      
      if ((pars[1] < se_boundary) |  (pars[1] > 1 - se_boundary)) {
        warning("The estimate of alpha1 is close to the boundary. Standard errors might not be valid.")
      }
      if ((pars[2] < se_boundary) |  (pars[2] > 1 - se_boundary)) {
        warning("The estimate of alpha2 is close to the boundary. Standard errors might not be valid.")
      }
      if ((pars[3] < se_boundary) |  (pars[3] > 1 - se_boundary)) {
        warning("The estimate of alpha3 is close to the boundary. Standard errors might not be valid.")
      }
      
      if (julia_installed) {
        addJuliaFunctions()
        julia_reg <- JuliaConnectoR::juliaCall("Coconots.create_julia_dict", 
                                               list("parameter", "covariance_matrix", "log_likelihood",
                                                    "type", "order", "data", "covariates",
                                                    "link", "starting_values", "optimizer", 
                                                    "lower_bounds", "upper_bounds", "optimization", 
                                                    "max_loop"),
                                               list(pars, inv_hes, likelihood,
                                                    type, order, data, xreg,
                                                    link_function, NULL, NULL, NULL, NULL,
                                                    NULL, NULL))
      } else {julia_reg = NULL}
      
      list_func <- list(
        "par" = pars,
        "gradient" = gra, "hessian" = hes, "inv hessian" = inv_hes,
        "se" = se, "ts" = data, "cov" = xreg, "type" = "Poisson",
        "order" = 2, "seasonality" = seasonality, "likelihood" = likelihood,
        "duration" = end_time - start_time, julia_reg = julia_reg, "link_function"=link_function
      )

      return(list_func)
    } # end PAR2
    return(PAR2(data))
  } # end if PAR2


  # GP2
  if ((type == "GP") & (order == 2)) {
    GP2 <- function(data) {
      xreg_matrix <- data.matrix(xreg)
      mlf <- function(par, data) {
        alpha1 <- par[1]
        alpha2 <- par[2]
        alpha3 <- par[3]
        eta <- par[4]
        vec_lambda <- par[-(1:4)]

        T <- length(data)
        mlef <- likelihoodGP2cov(alpha1, alpha2, alpha3, eta, vec_lambda, T, seasonality[1], seasonality[2], data, xreg_matrix, link_function, cores)

        return(mlef)
      } # end function

      ## starting values
      if (length(start) == 0) {
        acf <- forecast::Acf(data, plot = FALSE, lag.max = lag.max)
        p1 <- acf$acf[[seasonality[1] + 1]]
        p2 <- acf$acf[[seasonality[2] + 1]]

        mean_data <- mean(data)
        idx <- (seasonality[2] + 1):length(data)
        v <- (data[idx] - mean_data) * (data[idx - seasonality[1]] - mean_data) * (data[idx - seasonality[2]] - mean_data)
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
              alpha1_s <- alpha1_s * iteration.start.val
              alpha2_s <- alpha2_s * iteration.start.val
              alpha3_s <- alpha3_s * iteration.start.val
            }
            warning(paste("The initial sum of 2*alpha1, alpha2 and alpha3 had not been in the feasible region and was adjusted", ""))
          }


          if (2 * alpha1_s + alpha3_s > 1) {
            while (2 * alpha1_s + alpha3_s > 1) {
              alpha1_s <- alpha1_s * iteration.start.val
              alpha3_s <- alpha3_s * iteration.start.val
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

        sta <- c(alpha1_s, alpha2_s, alpha3_s, eta_s,
                 covariateStartingValues(data, xreg, link_function,
                                         (1 - alpha1_s - alpha2_s - alpha3_s) * (1 - eta_s),
                                         starting_values_covariates, b.beta))

      }

      if (length(start) != 0) {
        sta <- start
        if (length(start) != 4 + ncol(xreg)) {
          stop(paste("Number of initial starting values must equal", 4 + ncol(xreg), "for the GP 2 model with", ncol(xreg), "covariates"))
        }
      }


      # derivatives of mlf by automatic differentiation (C++ dual numbers)
      mlf_derivs <- function(par, data, deriv = 2L) {
        coco_ad_derivs(par, data, type, 2, seasonality, xreg_matrix, link_function, deriv)
      }
      gradient <- function(par, data) mlf_derivs(par, data, 1L)$gradient

      derivs_opt <- NULL
      if (constrained.optim == TRUE) {
        # constrained.optimaints
        ui <- rbind(
          c(-1, -1, -1, 0), c(-2, 0, -1, 0), c(1, 0, 0, 0), c(0, 1, 0, 0), c(0, 0, 1, 0),
          c(0, 0, 0, 1), c(0, 0, 0, -1)
        )
        matrix1 <- matrix(0, 7, ncol(xreg)) # for additional columns

        matrix2 <- matrix(0, ncol(xreg) , 4) # for additional rows
        einheitsmatrix <- diag(ncol(xreg) )
        lower_part <- cbind(matrix2, einheitsmatrix)

        ui <- cbind(ui, matrix1) # right part of the matrix
        ui <- rbind(ui, lower_part)

        ci <- c(-1, -1, 0, 0, 0, 0, -1, rep(b.beta, ncol(xreg) ))

        # constrained.optimained optimization process
        if (method_optim == "BFGS") {
          fit <- stats::constrOptim(
            theta = sta, f = mlf, grad=gradient, ui = ui, ci = ci, data = data, #mu = mu,
            method = method_optim, control = optim_control,
            #outer.iterations = outer.it, outer.eps = outer.eps,
            hessian = FALSE)          
          fit <- coco_polish_fit(fit, mlf, ui, ci, gr = gradient, data = data)
        }else{
          constr_fit <- function() {
            fit <- coco_constr_nm_fit(
              sta, mlf, gradient, ui, ci, mu = mu, outer.it = outer.it,
              outer.eps = outer.eps, control = optim_control, data = data
            )
            coco_polish_fit(fit, mlf, ui, ci, gr = gradient, data = data)
          }
          if (method_optim == "Nelder-Mead") {
            # L-BFGS-B with AD gradient from the starting values, accepted
            # after an AD optimality check; else the constrOptim fit above
            opt <- coco_fit_lbfgsb_checked(
              sta, mlf, gradient, function(p) mlf_derivs(p, data, 2L), ui, ci,
              fallback = constr_fit, data = data
            )
            fit <- opt$fit
            derivs_opt <- opt$derivs
          } else {
            fit <- constr_fit()
          }
        }
      } # end constrained.optim True

      if (constrained.optim == FALSE) {
        if (length(unconstrained.optim.lower) != 4 + ncol(xreg) & (unconstrained.optim.lower != -Inf)) {
          stop(paste("Number of lower bounds must equal", 4 + ncol(xreg), "for the GP 2 model with", ncol(xreg), "covariates"))
        }

        if (length(unconstrained.optim.upper) != 4 + ncol(xreg) & (unconstrained.optim.upper != Inf)) {
          stop(paste("Number of upper bounds must equal", 4 + ncol(xreg), "for the GP 2 model with", ncol(xreg), "covariates"))
        }

        fit <- stats::optim(
          par = sta, fn = mlf, gr = NULL, method = c("Nelder-Mead"), data = data,
          lower = unconstrained.optim.lower, upper = unconstrained.optim.upper,
          control = optim_control, hessian = FALSE
        )
      } # end constrained.optim False

      pars <- as.vector(unlist(fit$par, use.names = FALSE))
      lambs <- paste0("beta", 1:ncov)
      names(pars) <- c("alpha1", "alpha2", "alpha3", "eta", lambs)
      T <- length(data)
      likelihood <- -likelihoodGP2cov(pars[1], pars[2], pars[3], pars[4], pars[-c(1:4)], T, seasonality[1], seasonality[2], data, xreg_matrix, link_function)


      derivs <- if (!is.null(derivs_opt)) derivs_opt else mlf_derivs(pars, data, 2L)
      gra <- derivs$gradient
      hes <- derivs$hessian
      inv_hes <- solve(hes)
      se <- diag(inv_hes)^0.5
      names(se) <- c("alpha1", "alpha2", "alpha3", "eta", lambs)

      end_time <- Sys.time()
      
      if ((pars[1] < se_boundary) |  (pars[1] > 1 - se_boundary)) {
        warning("The estimate of alpha1 is close to the boundary. Standard errors might not be valid.")
      }
      if ((pars[2] < se_boundary) |  (pars[2] > 1 - se_boundary)) {
        warning("The estimate of alpha2 is close to the boundary. Standard errors might not be valid.")
      }
      if ((pars[3] < se_boundary) |  (pars[3] > 1 - se_boundary)) {
        warning("The estimate of alpha3 is close to the boundary. Standard errors might not be valid.")
      }
      
      if ((pars[4] < se_boundary) |  (pars[4] > 1 - se_boundary)) {
        warning("The estimate of eta is close to the boundary. Standard errors might not be valid.")
      }
      
      if (julia_installed) {
        addJuliaFunctions()
        julia_reg <- JuliaConnectoR::juliaCall("Coconots.create_julia_dict", 
                                               list("parameter", "covariance_matrix", "log_likelihood",
                                                    "type", "order", "data", "covariates",
                                                    "link", "starting_values", "optimizer", 
                                                    "lower_bounds", "upper_bounds", "optimization", 
                                                    "max_loop"),
                                               list(pars, inv_hes, likelihood,
                                                    type, order, data, xreg,
                                                    link_function, NULL, NULL, NULL, NULL,
                                                    NULL, NULL))
      } else {julia_reg = NULL}
      
      list_func <- list(
        "par" = pars,
        "gradient" = gra, "hessian" = hes, "inv hessian" = inv_hes,
        "se" = se, "ts" = data, "cov" = xreg, "type" = "GP", "order" = 2,
        "seasonality" = seasonality, "likelihood" = likelihood,
        "duration" = end_time - start_time, julia_reg = julia_reg, "link_function"=link_function
      )

      return(list_func)
    } # end GP2
    return(GP2(data))
  } # end if GP2
}

# Default starting values of the covariate coefficients beta.
# The (GP-)INAR models have E[X_t] = lambda / ((1 - sum(alpha)) (1 - eta)), so for
# slowly varying lambda_t the innovation rate is lambda_t = scale * mu_t with
# scale = (1 - sum(alpha)) (1 - eta) evaluated at the alpha/eta starts and mu_t the
# marginal mean of the counts. Per link:
#   log:            Poisson GLM of the counts on xreg with offset -log(scale)
#                   (with an intercept column: intercept shifted by log(scale));
#   identity, relu: scale times the least-squares fit of the counts on xreg;
#   softplus:       least-squares fit of softplus^{-1}(lambda_t) on xreg with
#                   lambda_t = exp(xreg beta_log) from the log-link GLM above.
# Returns `fallback` (the marginal fit) if scale is not in (0, 1] or the result is
# not admissible (non-finite, beta <= b.beta, or lambda_t <= 0 for some row).
# Mirrored in Coconots.jl `compute_covariate_starting_values`.
covariateStartingValues <- function(data, xreg, link_function, scale, fallback, b.beta) {
  if (!is.finite(scale) || scale <= 0 || scale > 1) {
    return(fallback)
  }
  y <- as.numeric(data)
  X <- data.matrix(xreg)
  beta <- tryCatch({
    if (link_function %in% c("identity", "relu")) {
      scale * stats::lm.fit(X, y)$coefficients
    } else {
      beta_log <- stats::glm.fit(X, y, family = stats::poisson(),
                                 offset = rep(-log(scale), length(y)))$coefficients
      if (link_function == "log") {
        beta_log
      } else {
        lambda <- exp(drop(X %*% beta_log))
        stats::lm.fit(X, lambda + log(-expm1(-lambda)))$coefficients
      }
    }
  }, error = function(e) NULL)
  if (is.null(beta) || length(beta) != length(fallback) ||
      any(!is.finite(beta)) || any(beta <= b.beta)) {
    return(fallback)
  }
  lin <- drop(X %*% beta)
  if (any(!is.finite(lin)) ||
      (link_function %in% c("identity", "relu") && any(lin <= 0)) ||
      (link_function == "log" && any(exp(lin) <= 0))) {
    return(fallback)
  }
  beta <- as.numeric(beta)
  names(beta) <- names(fallback)
  return(beta)
}

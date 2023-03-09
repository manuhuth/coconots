transformJuliaRegOutputToR <- function(xreg, pars, grad, hes, inv_hes, se,
                                       data, type, order, likelihood, end_time, start_time, julia_reg){
  if (is.null(xreg)){
    par_use <- c(pars[length(pars)], pars[1:(length(pars)-1)])
    se_use <- c(se[length(se)], se[1:(length(se)-1)])
  } else {
    par_use <- pars
    se_use <- se
  }
  
  return(list(
    "par" = par_use, "gradient" = grad, cov=xreg,
    "hessian" = hes, "inv hessian" = inv_hes, "se" = se_use, "ts" = data,
    "type" = type, "order" = order, "seasonality" = c(1,2), "likelihood" = likelihood,
    "duration" = end_time - start_time, julia_reg = julia_reg
  ))
}
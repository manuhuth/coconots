cocoSim <- function(type, order, par, length, xreg = NULL, seasonality = c(1, 2), init = NULL) {
  

  
  
  if (is.null(xreg)) {
    if (order == 2){
      
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
      if ((par[2] < 0) | (par[2] >= 1))  {
        stop("alpha must be within the unit interval")
      }
      if ((type == "GP") & ((par[3] < 0) | (par[3] >= 1))){
        stop("eta must be within the unit interval")
      }
    }
    
    if (is.null(init)) {
      length_burn_in <- 200
    } else {length_burn_in <- 0}
    
    size <- length + length_burn_in
    output <- cocoSim_base(type = type, order = order, par = par, size = size, seasonality = seasonality, init = init)
    output$data <- output$data[(length_burn_in+1):(length+length_burn_in)]
  } else {
    if (order == 2){
      
      if (par[1] + par[3] > 1) {
        stop("The condition alpha_1 + alpha_3 < 1 ist not satisfied")
      }
      
      if ((par[1] < 0) | (par[2] < 0) | (par[3] < 0) | (par[1] >= 1) | (par[2] >= 1) | (par[3] >= 1)) {
        stop("The alpha parameters must be within the unit interval")
      }

      if ((type == "GP") & ((par[4] < 0) | (par[4] >= 1))){
        stop("eta must be within the unit interval")
      }
      
    } else {
      if ((par[1] < 0) | (par[1] >= 1))  {
        stop("alpha must be within the unit interval")
      }
      if ((type == "GP") & ((par[2] < 0) | (par[2] >= 1))){
        stop("eta must be within the unit interval")
      }
    }
    output <- cocoSim_cov(
      type = type, order = order, par = par, size = length, xreg = xreg, seasonality = seasonality, init = init
    )
  }

  return(output)
}

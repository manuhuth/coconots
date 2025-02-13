#' @title Computes Scores for Various Models Maintaining a Common Sample 
#'
#' @description This function computes log, quadrtic and ranked probability scores for Poisson and Generalized Poisson models.
#'
#' @param data A numeric vector containing the data to be used for modeling
#' @param models A character string specifying which models to use. Default is `"all"`, which uses both Poisson and GP models.
#' @param print.progress A logical value indicating whether to print progress messages (Default: `TRUE`).
#' @param max_x_score An integer which is used as the maximum count for the computation
#'  of the score (defaul: `50`)
#' @param julia if TRUE, \code{cocoSoc} is run with \proglang{julia} (default: FALSE)
#' @param ... Additional arguments to be passed to the `cocoReg` function.
#' 
#' @author Manuel Huth
#' 
#' @details Supports model selection by computing score over a range of models while maintaining a common sample and 
#' a common specification.
#' @return A list of class `"cocoSoc"` containing:
#' \describe{
#'   \item{fits}{A list of fitted model objects.}
#'   \item{scores_list}{A list of score objects for each model.}
#'   \item{scores_df}{A data frame containing the logarithmic, quadratic, and ranked probability scores for each model.}
#' }
#'
#'@examples
#' pars <- c(1.3, 0.25, 0.03, 0.2, 0.3)
#' set.seed(12345)
#' data <- cocoSim(order = 2, type = "GP", par = pars, length = 500)
#' soc <- cocoSoc(data, julia=T)
#' summary(soc)
#'
#' @export
cocoSoc <- function(data, models = "all", print.progress=TRUE,
                       max_x_score=50, julia = FALSE, ...
          ) {
  fits <- list()
  scores <- list()
  
  n_models <- ifelse(models == "all", 2, 4)
  
  if (models == "all"){
    models_use <- c("Poisson", "GP")
  } else {
    models_use <- c(models)
  } 
  
  #exclude first observation for order 1 models
  index <- 1
  for (type in models_use){
    for (order in 1:2){
      data_use <- data[(order*-1 + 3):length(data)]
      if (type == "Poisson"){name <- "PAR"} else {name <- "GP"}
      
      if (print.progress){
        print(paste0("Computing ", name, order, " model."))
      }
      
      name_mod <- paste0(name, order)
      
      fits[[name_mod]] <- cocoReg(type, order, data_use, ...)

      scores[[name_mod]] <- cocoScore(fits[[name_mod]], max_x = max_x_score, julia=T)
      index <- index + 1
    }
  }
  
  log_score <- c()
  quad_score <- c()
  rps_score <- c()
  for (model in names(scores)){
    log_score <- c(log_score, scores[[model]]$log.score)
    quad_score <- c(quad_score, scores[[model]]$quad.score)
    rps_score <- c(rps_score, scores[[model]]$rps.score)
  }
  
  df <- data.frame(names(scores))
  colnames(df) <- c(" ")
  df[c("Logarithmic Score", "Quadratic Score", "Ranked Probability Score")] <- cbind(log_score, quad_score, rps_score)

  list_out <- list("fits" = fits, "scores_list" = scores, "scores_df" = df)
  class(list_out) <- "cocoSoc"
  
  return(list_out)

  } # end function

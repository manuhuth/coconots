#' @title Compute Scores for Various Models
#'
#' @description This function computes and returns scores for Poisson and Generalized Poisson models. It resembles Stata's `varsoc` function.
#'
#' @param data A numeric vector containing the data to be used for modeling.
#' @param models A character string specifying which models to use. Default is `"all"`, which uses both Poisson and GP models.
#' @param print.progress A logical value indicating whether to print progress messages. Default is `TRUE`.
#' @param val.num.score A numeric value specifying the validation number score. Default is `1e-10`.
#' @param ... Additional arguments to be passed to the `cocoReg` function.
#'
#' @return A list of class `"cocoVarsoc"` containing:
#' \describe{
#'   \item{fits}{A list of fitted model objects.}
#'   \item{scores_list}{A list of score objects for each model.}
#'   \item{scores_df}{A data frame containing the logarithmic, quadratic, and ranked probability scores for each model.}
#' }
#'
#'
#' @export
cocoVarsoc <- function(data, models = "all", print.progress=TRUE,
                       val.num.score=1e-10, ...
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

      scores[[name_mod]] <- cocoScore(fits[[name_mod]], val.num = val.num.score, julia=T)
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
  class(list_out) <- "cocoVarsoc"
  
  return(list_out)

  } # end function

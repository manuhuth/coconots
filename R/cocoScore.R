#' @title Scoring Rule Based Model Assessment Procedure
#' @description The function calculates the log, quadratic and ranked probability scores for assessing relative performance of a fitted model as proposed by Czado et al. (2009).
#' @param coco An object of class coco
#' @param max_x An integer which is used as the maximum count for the computation
#'  of the score. The default value is `50`
#' @param julia if TRUE, the scores are computed with Julia.
#' @return a list containing the log score, quadratic score and ranked probability score.
#' @details Scoring rules assign a numerical score based on the predictive
#'  distribution and the observed data  to measure the quality of probabilistic predictions.
#' They are provided here as a model selection tool and are computed as
#'  averages over the relevant set of (in-sample) predictions. Scoring rules are, generally, negatively oriented
#' penalties that one seeks to minimize. The literature has developed a large number of scoring
#' rules and, unless there is a unique and clearly defined underlying decision problem,
#' there is no automatic choice of a (proper) scoring rule to be used in any given situation.
#' Therefore, the use of a variety of scoring rules may be appropriate to take advantage of
#' specific emphases and strengths. Three proper scoring rules
#' (for a definition of the concept of propriety see Gneiting and Raftery, 2007)
#' which Jung, McCabe and Tremayne (2016) found to be particularly useful are implemented.
#' For more information see the references listed below.
#' @references 
#' Czado, C. and Gneitling, T. and Held, L. (2009) Predictive Model Assessment for Count Data. \emph{Biometrics}, \bold{65}, 4, 1254--1261.
#'
#' Gneiting, T. and Raftery, A. E. (2007) Strictly proper scoring rules, prediction, and estimation. \emph{Journal
#' of the American Statistical Association}, 102:359-378.
#'
#' Jung, Robert C., Brendan P. M. McCabe, and Andrew R. Tremayne. (2016). Model validation and diagnostics. \emph{In Handbook of Discrete
#' Valued Time Series}. Edited by Richard A. Davis, Scott H. Holan, Robert Lund and Nalini Ravishanker. Boca Raton: Chapman and
#' Hall, pp. 189--218.
#'
#' Jung, R. C. and Tremayne, A. R. (2011) Convolution-closed models for count timeseries with applications. \emph{Journal of Time Series Analysis}, \bold{32}, 3, 268--280.
#' @author Manuel Huth
#' @examples
#' lambda <- 1
#' alpha <- 0.4
#' set.seed(12345)
#' data <- cocoSim(order = 1, type = "Poisson", par = c(lambda, alpha), length = 100)
#' #julia_installed = TRUE ensures that the fit object
#' #is compatible with the julia cocoScore implementation 
#' fit <- cocoReg(order = 1, type = "Poisson", data = data)
#'
#' #assessment using scoring rules - R implementation
#' score_r <- cocoScore(fit)
#' @export

cocoScore <- function(coco, max_x = 50, julia=FALSE) {
  
  
  if (!is.null(coco$julia_reg) & julia){
    addJuliaFunctions()
    coco_score <- JuliaConnectoR::juliaGet( JuliaConnectoR::juliaCall("compute_scores", coco$julia_reg, max_x))
    return(list("log.score" = coco_score$values[[1]], "quad.score" = coco_score$values[[2]],
                 "rps.score" = coco_score$values[[3]]
                #"aic" = 2*length(coco$par) - 2 * coco$likelihood,
                #"bic" = length(coco$par) * log(length(coco$ts)) - 2 * coco$likelihood
                ))
  } else {
  T <- length(coco$ts)
  par <- coco$par
  seas <- c(1, 2) #will be used as argument in future versions
  data <- coco$ts

  if (coco$type == "Poisson"){
    par <- c(par, 0)
  }
  
  if (!is.null(coco$cov)){
    xreg <- coco$cov
    betas <- par[(length(par)-ncol(xreg)+1):length(par)]
    lambdas <- exp(as.matrix(xreg) %*% betas)
    par_original <- par
    
  }
  
  #Poisson/GP 1
  if  (coco$order == 1){
    calculate_scores <- function(t) {
      if (!is.null(coco$cov)){
        par <- c(lambdas[t], par_original[1:(length(par_original) - ncol(xreg))]) 
      }
      
      # log score
      log_score <- -log(dGP1(data[t], data[t - seas[1]], par))
      
      # quadratic score
      norm_p <- sum(sapply(0:max_x, function(j) dGP1(j, data[t - seas[1]], par)^2))
      quad_score <- -2 * dGP1(data[t], data[t - seas[1]], par) + norm_p
      
      # ranked probability score
      sum_sq_diff <- sum(sapply(0:max_x, function(j) (pGP1(j, data[t - seas[1]], par) - (j >= data[t]))^2))
      rps_score <- sum_sq_diff
      
      return(c(log_score, quad_score, rps_score))
    }
    
    # Use sapply for the outer loop
    scores <- sapply((seas[1] + 1):T, calculate_scores)
    
    # Summarize scores
    log.score <- mean(scores[1, ])
    quad.score <- mean(scores[2, ])
    rps.score <- mean(scores[3, ])

  }# Poisson/GP 1


  #Poisson/GP 2
  if (coco$order == 2){
    calculate_scores <- function(t) {
      if (!is.null(coco$cov)){
        par <- c(lambdas[t], par_original[1:(length(par_original) - ncol(xreg))]) 
      }
      
      # log score
      log_score <- -log(dGP2(data[t], data[t - seas[1]], data[t - seas[2]], par))
      
      # quadratic score
      norm_p <- sum(sapply(0:max_x, function(j) dGP2(j, data[t - seas[1]], data[t - seas[2]], par)^2))
      quad_score <- -2 * dGP2(data[t], data[t - seas[1]], data[t - seas[2]], par) + norm_p
      
      # ranked probability score
      sum_sq_diff <- sum(sapply(0:max_x, function(j) (pGP2(j, data[t - seas[1]], data[t - seas[2]], par) - (j >= data[t]))^2))
      rps_score <- sum_sq_diff
      
      return(c(log_score, quad_score, rps_score))
    }
    
    # Use sapply for the outer loop
    scores <- sapply((seas[2] + 1):T, calculate_scores)
    
    # Summarize scores
    log.score <- mean(scores[1, ])
    quad.score <- mean(scores[2, ])
    rps.score <- mean(scores[3, ])
  }# Poisson/GP 2

  list_out <- list("log.score" = log.score, "quad.score" = quad.score, "rps.score" = rps.score)

  }#end julia if
  
  #list_out["bic"] <- length(coco$par) * log(length(coco$ts)) - 2 * coco$likelihood
  #list_out["aic"] <- 2*length(coco$par) - 2 * coco$likelihood
  return(list_out)

}# end function

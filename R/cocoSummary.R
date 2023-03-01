#' @title cocoSummary
#' @description Summarizes the results of a cocoReg analysis
#' @param coco an object of class cocoReg
#' @param score logical indicating whether to include score statistics in the summary
#' @return summary of the cocoReg analysis
#' @author Manuel Huth
#' @export

cocoSummary <- function(coco, score=FALSE) {
  df <- data.frame(cbind(round(coco$par,4), round(coco$se,4)) )
  colnames(df) <- c("Estimate", "Std. Error")
  
  cat("Coefficients:\n")
  print(df, print.gap=3, quote=FALSE, na.print="")
  
  
  if (isTRUE(score)) { 
    sc <- cocoScore(coco)
    cat("\nType:", coco$type,    "\nOrder:", coco$order,
        "\n\nLog-likelihood:", round(coco$likelihood,4),
        "\nLogarithmic score:", round(sc$log.score,4),
        "\nQuadratic score:", round(sc$quad.score,4),
        "\nRanked probability score", round(sc$rps.score,4)
    )
  } else{
    cat("\nType:", coco$type,    "\nOrder:", coco$order,
        "\n\nLog-likelihood:", round(coco$likelihood,4)
    )
  }


}
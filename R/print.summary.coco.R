#' @export
print.summary.coco <- function(x, ...) {
  coco <- x
  df <- data.frame(cbind(round(coco$par,4), round(coco$se,4), round(coco$par/coco$se,4)) )
  colnames(df) <- c("Estimate", "Std. Error", "t")
  
  if ((coco$order == 1) & (coco$type == "Poisson")){
    rows <- c("lambda", "alpha")
  } else if ((coco$order == 1) & (coco$type == "GP")){
    rows <- c("lambda", "alpha", "eta")
  } else if ((coco$order == 2) & (coco$type == "Poisson")){
    rows <- c("lambda", "alpha1", "alpha2", "alpha3")
  } else if ((coco$order == 2) & (coco$type == "GP")){
    rows <- c("lambda", "alpha1", "alpha2", "alpha3", "eta")
  }
  
  if (!is.null(coco$cov)){
    if (is.null(colnames(coco$cov))){
      colnames(coco$cov) <- paste0("X", 1:ncol(coco$cov))
    }
    rows <- c(rows[rows != "lambda"], colnames(coco$cov) )
  }
  
  rownames(df) <- rows
  
  cat("Coefficients:\n")
  print(df, print.gap=3, quote=FALSE, na.print="")
  if (is.null(coco$julia_reg)){
    julia = FALSE
  } else{
    julia = TRUE
  }
  
  if (isTRUE(coco$score)) { 
    sc <- cocoScore(coco, julia=julia)
    cat("\nType:", coco$type,    "\nOrder:", coco$order,
        "\n\nLog-likelihood:", round(coco$likelihood,4),
        "\nLogarithmic score:", round(sc$log.score,4),
        "\nQuadratic score:", round(sc$quad.score,4),
        "\nRanked probability score", round(sc$rps.score,4)
        #"\nAIC:", round(sc$aic,4),
        #"\nBIC:", round(sc$bic,4)
    )
  } else{
    cat("\nType:", coco$type,    "\nOrder:", coco$order,
        "\n\nLog-likelihood:", round(coco$likelihood,4)
        #"\nAIC:", round(2*length(coco$par) - 2 * coco$likelihood,4),
        #"\nBIC:", round(length(coco$par) * log(length(coco$ts)) - 2 * coco$likelihood,4)
    )
  }
}
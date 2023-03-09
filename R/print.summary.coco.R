#' @export
print.summary.coco <- function(x, ...) {
  coco <- x
  df <- data.frame(cbind(round(coco$par,4), round(coco$se,4), round(coco$par/coco$se,4)) )
  colnames(df) <- c("Estimate", "Std. Error", "t")
  
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
    )
  } else{
    cat("\nType:", coco$type,    "\nOrder:", coco$order,
        "\n\nLog-likelihood:", round(coco$likelihood,4)
    )
  }
}
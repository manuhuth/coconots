#' @importFrom ggplot2 autoplot theme_bw ggtitle
#' @export
ggplot2::autoplot
#' @export
autoplot.cocoResid <- function(object, ...){
  forecast::ggAcf(object$pe.resid) + theme_bw() + ggtitle("Pearson Residuals") + ggplot2::theme(text = ggplot2::element_text(size = 20)) 
}

#' @export
plot.cocoResid <- function(x, ...) {
  p <- autoplot(
    x,
    ...
  )
  print(p)
}

#'@export
print.cocoResid <- function(x, ...) {
  print(autoplot(x, ...,))
  invisible(x)
}
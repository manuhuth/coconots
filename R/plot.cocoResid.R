#' @importFrom ggplot2 autoplot theme_bw ggtitle
#' @export
ggplot2::autoplot
#' @export
autoplot.cocoResid <- function(object, ...){
  forecast::ggAcf(object$pe.resid) + theme_bw() + ggtitle("Pearson Residuals")
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
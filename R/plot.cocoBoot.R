#' @importFrom ggplot2 autoplot
#' @export
ggplot2::autoplot
#' @export
autoplot.cocoBoot <- function(object, ...){
  
  pl <- ggplot2::ggplot(object$df_plot, ggplot2::aes_string("x", "y")) +
      ggplot2::geom_point()+
      ggplot2::geom_line()+
      ggplot2::geom_ribbon(data=object$df_plot,ggplot2::aes_string(ymin="lower",ymax="upper"),
                           fill="steelblue", alpha=0.3) +
    ggplot2::theme_bw() + ggplot2::xlab("Lags") + ggplot2::ylab("Autocorrelation") +
    ggplot2::ggtitle("Parametric Bootstrap") +
    ggplot2::theme(text = ggplot2::element_text(size = 20))
  pl
  
}

#' @export
plot.cocoBoot <- function(x, ...) {
  p <- autoplot(
    x,
    ...
  )
  print(p)
}

#'@export
print.cocoBoot <- function(x, ...) {
  print(autoplot(x, ...,))
  invisible(x)
}
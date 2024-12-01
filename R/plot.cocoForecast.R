#' @importFrom ggplot2 autoplot
#' @export
ggplot2::autoplot
#' @export
autoplot.cocoForecast <- function(object, breaks=NULL, width=0.1, ...){
  
  pl <- ggplot2::ggplot(mapping = ggplot2::aes(x = object$x, y = object$densities_plot)) +
      ggplot2::geom_bar(stat="identity", position="dodge", width=width) + 
      ggplot2::labs(title = "Probability mass forecast", x = "Support", y = "Probability mass") +
    ggplot2::xlim(c(0, max(object$x))) +
      ggplot2::theme_bw() + ggplot2::theme(text = ggplot2::element_text(size = 20)) 

  pl
}

#' @export
plot.cocoForecast <- function(x, ...) {
  p <- autoplot(
    x,
    ...
  )
  suppressWarnings({print(p)})
}

#' @export
plot.cocoForecast <- function(x, ...) {
  p <- autoplot(
    x,
    ...
  )
  suppressWarnings({print(p)})
}
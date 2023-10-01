#' @importFrom ggplot2 autoplot
#' @export
ggplot2::autoplot
#' @export
autoplot.cocoForecast <- function(object, ...){
  
  pl <- ggplot2::ggplot(mapping = ggplot2::aes(x = object$x, y = object$densities_plot)) +
      ggplot2::geom_bar(stat="identity", position="dodge", width=0.04) + 
      ggplot2::labs(title = "Probability mass forecast", x = "Support", y = "Probability mass") +
      ggplot2::scale_x_continuous(breaks=object$x) +
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
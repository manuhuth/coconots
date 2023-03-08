#' @importFrom ggplot2 autoplot
#' @export
ggplot2::autoplot
#' @export
autoplot.cocoPit <- function(object, ...){
  
  df <- data.frame(object$`PIT values`, 1:rep(length(object$`PIT values`)),
                   rep(object$confidence_bands[1], length(object$`PIT values`)),
                   rep(object$confidence_bands[2], length(object$`PIT values`)))
  colnames(df) <- c("pit", "bins", "lower", "upper")
  df_bands <- df
  df_bands$bins[1] <- df$bins[1] - 0.25
  df_bands$bins[length(df$bins)] <- df$bins[length(df$bins)] + 0.25
  pl <- ggplot2::ggplot(data=df, mapping = ggplot2::aes(x = bins, y = pit)) +
      ggplot2::geom_bar(stat="identity", position="dodge", width=0.3) + 
      ggplot2::labs(title = "Pit histogram", x = "Bins", y = "") +
      ggplot2::scale_x_continuous(breaks=df$bins) +
      ggplot2::theme_bw() + ggplot2::theme(text = ggplot2::element_text(size = 20)) +
      ggplot2::geom_ribbon(data=df_bands, ggplot2::aes(ymin=lower,ymax=upper),
                         fill="steelblue", alpha=0.3) +
      ggplot2::ylim(c(0, min(1.05, max(df$pit)*1.5)))
  pl
}

#' @export
plot.cocoPit <- function(x, ...) {
  p <- autoplot(
    x,
    ...
  )
  print(p)
}

#'@export
print.cocoPit <- function(x, ...) {
  print(autoplot(x, ...,))
  invisible(x)
}
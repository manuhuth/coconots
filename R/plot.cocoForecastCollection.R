#' @importFrom ggplot2 autoplot
#' @export
ggplot2::autoplot
#' @export
autoplot.cocoForecastCollection <- function(object, forecast_type="mode",
                                            n_lags_plot=30, ...){
  return_forecast <- function(i){
    L <- object[[i]]
    c(df_data[nrow(df_data), "t"]+i, NA, L[[forecast_type]], L[["lower"]], L[["upper"]])
  }
  
  relevant_data <- object[[1]]$data[max(1,(length(object[[1]]$data)-n_lags_plot+1)):length(object[[1]]$data)]
  df_data <- cbind(max(1,(length(object[[1]]$data)-n_lags_plot+1)):length(object[[1]]$data), relevant_data, NA, NA, NA)
  colnames(df_data) <- c("t", "x", "forecast", "lower", "upper")
  
  df <- as.data.frame(t(sapply(1:length(object), return_forecast)))
  colnames(df) <- c("t", "x", "forecast", "lower", "upper")
  
  df_plot <- rbind(df_data, df)
  
  df_transition <- df_plot[(nrow(df_plot)-length(object)):nrow(df_plot) ,]
  df_transition[1, "forecast"] <- df_transition[1, "x"] 
  
  pl <- ggplot2::ggplot(mapping = ggplot2::aes(x = df_plot[, "t"], y = df_plot[, "x"])) +
      ggplot2::geom_point(ggplot2::aes(color="Data")) + ggplot2::geom_line() + 
      ggplot2::geom_point(mapping = ggplot2::aes(x = df[, "t"], y = df[, "forecast"],
                                                 color = "Prediction")) +
      ggplot2::geom_line(mapping = ggplot2::aes(x = df_transition[, "t"],
                                                y = df_transition[, "forecast"]),
                         color = "steelblue", linetype = "dotted") +
      ggplot2::geom_errorbar(mapping = ggplot2::aes(
                                                    ymin = df_plot[, "lower"],
                                                    ymax = df_plot[, "upper"]), color = "steelblue") +
      ggplot2::labs(title = paste0(forecast_type, " forecast"), x = "Time", y = "Count (x)") +
      ggplot2::theme_bw() + ggplot2::theme(text = ggplot2::element_text(size = 20)) +
    ggplot2::scale_color_manual(name='',
                       breaks=c('Data', 'Prediction'),
                       values=c('Data'='black', 'Prediction'='steelblue'))

  pl
}

#' @export
plot.cocoForecastCollection <- function(x, ...) {
  p <- autoplot(
    x,
    ...
  )
  suppressWarnings({print(p)})
}

#'@export
print.cocoForecastCollection <- function(x, ...) {
  suppressWarnings({print(autoplot(x, ...,))})
  invisible(x)
}
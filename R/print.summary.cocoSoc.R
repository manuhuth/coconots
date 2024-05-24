#' @export
print.summary.cocoSoc <- function(x, ...) {

  print(x$scores_df, print.gap=3, quote=FALSE, na.print="")

}
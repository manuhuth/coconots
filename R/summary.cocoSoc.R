#' @export
summary.cocoSoc <- function(object, ...) {
  class(object) <- "summary.cocoSoc"
  return(object)
}
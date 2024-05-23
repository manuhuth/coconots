#' @export
summary.cocoVarsoc <- function(object, ...) {
  class(object) <- "summary.cocoVarsoc"
  return(object)
}
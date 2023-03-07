#' @export
summary.coco <- function(object, ...) {
  class(object) <- "summary.coco"
  return(object)
}
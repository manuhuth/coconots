#' @export
summary.coco <- function(object, ..., score=FALSE) {
  object$score <- score
  class(object) <- "summary.coco"
  return(object)
}
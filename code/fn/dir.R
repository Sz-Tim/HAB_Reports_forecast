#' Shortcut for dir(..., full.names=T)
#'
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
dirf <- function(...) {
  dir(..., full.names=T)
}



#' Shortcut for dir(..., full.names=T, recursive=T)
#'
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
dirrf <- function(...) {
  dir(..., full.names=T, recursive=T)
}



#' Shortcut for dir(..., recursive=T)
#'
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
dirr <- function(...) {
  dir(..., recursive=T)
}

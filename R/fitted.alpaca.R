#' @title
#' Extract \code{feglm} Fitted Values 
#' @description
#' \code{fitted.alpaca} is a generic function which extract fitted values from
#' an object returned by \code{feglm}.
#' @param 
#' object an object of class \code{alpaca}.
#' @param 
#' ... other arguments.
#' @return
#' The function \code{fitted.alpaca} returns a vector of fitted values.
#' @seealso
#' \code{\link{feglm}}
#' @export
fitted.alpaca <- function(object, ...) predict(object, type = "response")
#' @title
#' Extract estimated Structural Parameters
#' @description
#' \code{coef.alpaca} is a generic function which extracts
#' estimated structural parameters from objects returned by \code{feglm}.
#' @param 
#' object an object of class \code{alpaca}.
#' @param 
#' ... other arguments.
#' @return
#' The function \code{coef.alpaca} returns a named vector of estimated 
#' structural parameters.
#' @seealso
#' \code{\link{feglm}}
#' @export
coef.alpaca <- function(object, ...) object[["coefficients"]]

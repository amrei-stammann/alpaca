#' @title
#' Extract Matrix of estimated Structural Parameters
#' @description
#' \code{coef.summary.alpaca} is a generic function which extracts
#' the matrix of estimated structural parameters from objects returned by
#' \code{feglm}.
#' @param 
#' object an object of class \code{summary.alpaca}.
#' @param 
#' ... other arguments.
#' @return
#' The function \code{coef.summary.alpaca} returns a named matrix of estimates 
#' related to the structural parameters.
#' @seealso
#' \code{\link{feglm}}
#' @export
coef.summary.alpaca <- function(object, ...) object[["cm"]]
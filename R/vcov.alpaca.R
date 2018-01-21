#' @title
#' Extract estimated Covariance Matrix of the estimated Structural Parameters
#' @description
#' \code{vcov.alpaca} extracts the estimated covariance matrix of the estimated
#' structural parameters from objects returned by \code{feglm}. The estimation
#' is equal to inverse of Hessian after convergence.
#' @param 
#' object an object of class \code{alpaca}.
#' @param 
#' ... other arguments.
#' @return
#' The function \code{vcov.alpaca} returns a named matrix.
#' @seealso
#' \code{\link{feglm}}
#' @export
vcov.alpaca <- function(object, ...) {
  # Compute eigenvalues to check if the Hessian is invertible and compute
  # its inverse.
  H <- - object[["Hessian"]]
  ev <- abs(eigen(H, symmetric = TRUE, only.values = TRUE)[["values"]])
  if (min(ev) > .Machine[["double.eps"]] * max(ev) * 10.0) {
    V <- solve(H)
  } else {
    V <- matrix(Inf, nrow = nrow(H), ncol = ncol(H))
  }
  
  # Return covariance estimate.
  V
}
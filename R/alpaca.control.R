#' @title
#' Set \code{alpaca} Control Parameters
#' @description
#' Set and change \code{alpaca} control parameters.
#' @param
#' step.tol tolerance level for one of the stopping conditions in the IRLS
#' algorithm. This specific stopping condition is based on the euclidean norm of
#' the step size in iteration \eqn{r} and can be expressed as follows: 
#' \eqn{||\boldsymbol{\beta}_{r} - \boldsymbol{\beta}_{r - 1}||_{2} < 
#' \text{tol}}{||step|| < tol}. Default is \code{1.0e-06}.
#' @param
#' grad.tol tolerance level for one of the stopping conditions in the IRLS
#' algorithm. This specific stopping condition is based on the euclidean norm of
#' the gradient in iteration \eqn{r} and can be expressed as follows: 
#' \eqn{||\mathbf{g}_{r}||_{2} < \text{tol}}{||g|| < tol}.
#' Default is \code{1.0e-05}.
#' @param
#' dev.tol tolerance level for one of the stopping conditions in the IRLS
#' algorithm. This specific stopping condition is based on the relative change
#' of the deviance in iteration \eqn{r} and can be expressed as follows: 
#' \eqn{\Delta \text{dev}_{r} / \text{dev}_{r} < \text{tol}}{\Delta dev / dev <
#'  tol}. Default is \code{1.0e-08}.
#' @param
#' pseudo.tol tolerance level for the stopping condition of the pseudo demeaning
#' algorithm. The stopping condition is based on the euclidean norm of the step
#' size in iteration \eqn{i} and can be expressed as follows: 
#' \eqn{||\mathbf{v}_{i} - \mathbf{v}_{i - 1}||_{2} < \text{tol}}{||step|| <
#'  tol}. Default is \code{1.0e-08}.
#' @param
#' rho.tol tolerance level for the stephalving in iteration \eqn{r}. Stephalving
#' only takes place if the deviance in iteration \eqn{r} is worse than the one
#' in the previous iteration. The stopping condition can be expressed as
#' follows: \eqn{\rho < \text{tol}}{\rho < tol}. Default is \code{1.0e-03}.
#' @param
#' iter.max unsigned integer indicating the maximum number of iterations of the
#' IRLS algorithm.
#' @param
#' trace unsigned integer indicating if output should be produced for each
#' iteration. Default is \code{0}. See \code{Details}.
#' @param
#' drop.pc logical indicating to drop observations that are perfectly classified
#' and hence do not contribute to the log-likelihood. See \code{Details}.
#' Default is \code{TRUE}.
#' @param
#' ... other arguments.
#' @details
#' ...
#' @return
#' The function \code{alpaca.control} returns a named list of control 
#' parameters.
#' @seealso
#' \code{\link{feglm}}
#' @export
alpaca.control <- function(step.tol = 1.0e-06,
                           grad.tol = 1.0e-05,
                           dev.tol = 1.0e-08,
                           pseudo.tol = 1.0e-08,
                           rho.tol = 1.0e-03,
                           iter.max = 100L,
                           trace = 0L,
                           drop.pc = TRUE,
                           ...) {
  # Check validity of tolerance parameters.
  if (step.tol <= 0.0 || grad.tol <= 0.0 || dev.tol <= 0.0 || 
      pseudo.tol <= 0.0 || rho.tol <= 0.0) {
    stop("All tolerance paramerters should be greater than zero.")
  }
  
  # Check validity of 'iter.max'.
  if (iter.max <= 1L) {
    stop("Maximum number of iterations should be at least one.")
  }
  
  # Return list with control parameters.
  list(step.tol = step.tol,
       grad.tol = grad.tol,
       dev.tol = dev.tol,
       pseudo.tol = pseudo.tol,
       rho.tol = rho.tol,
       iter.max = as.integer(iter.max),
       trace = as.integer(trace),
       drop.pc = as.logical(drop.pc))
}
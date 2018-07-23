#' @title
#' Set \code{feglm} Control Parameters
#' @description
#' Set and change parameters used for fitting \code{feglm}.
#' @param
#' dev.tol tolerance level for the first stopping condition of the maximization routine. The 
#' stopping condition is based on the relative change of the deviance in iteration \eqn{r}
#' and can be expressed as follows: \eqn{(\text{dev}_{r - 1} \text{dev}_{r}) / \text{dev}_{r} < 
#' \text{tol}}{\Delta dev / dev < tol}. Default is \code{1.0e-08}.
#' @param
#' step.tol tolerance level for the second stopping condition of the maximization routine. The
#' stopping condition is based on the euclidean norm of the step size in iteration \eqn{r}
#' and can be expressed as follows: \eqn{\lvert\boldsymbol{\beta}_{r} - 
#' \boldsymbol{\beta}_{r - 1}\rvert_{2} < \text{tol}}{||\Delta \beta|| < tol}. Default is
#' \code{1.0e-08}.
#' @param
#' pseudo.tol tolerance level for the stopping condition of the \dQuote{pseudo demeaning} algorithm.
#' The stopping condition is based on the relative change of euclidean norm in iteration \eqn{i} and
#' can be expressed as follows: \eqn{\lvert\mathbf{v}_{i} - \mathbf{v}_{i - 1}\rvert_{2} < 
#' \text{tol} \lvert\mathbf{v}_{i - 1}\rvert}{||\Delta v|| / ||v_old|| < tol}. Default is
#' \code{1.0e-05}.
#' @param
#' rho.tol tolerance level for the stephalving in the maximization routine. Stephalving only takes
#' place if the deviance in iteration \eqn{r} is larger than the one of the previous iteration. If 
#' this is the case, 
#' \eqn{\lvert\boldsymbol{\beta}_{r} - \boldsymbol{\beta}_{r - 1}\rvert_{2}}{||\Delta \beta||} is 
#' halfed until the deviance is less or equal compared to the deviance of the previous iteration. 
#' Stephalving fails if the the following condition holds: \eqn{\rho < \text{tol}}{\rho < tol}, 
#' where \eqn{\rho}{\rho} is the stepcorrection factor. If stephalving fails the maximization
#' routine is canceled. Default is \code{1.0e-04}.
#' @param
#' iter.max unsigned integer indicating the maximum number of iterations in the maximization
#' routine.
#' @param
#' trace unsigned integer indicating if output should be produced in each iteration. Default is
#' \code{0}.
#' @param
#' drop.pc logical indicating to drop observations that are perfectly classified and hence do not
#' contribute to the log-likelihood. This option is useful to reduce the computational costs
#' of the maximization problem, since it reduces the number of observations and does not change the
#' estimates. Default is \code{TRUE}.
#' @return
#' The function \code{feglm.control} returns a named list of control 
#' parameters.
#' @seealso
#' \code{\link{feglm}}
#' @export
feglm.control <- function(dev.tol    = 1.0e-08,
                          step.tol   = 1.0e-08,
                          pseudo.tol = 1.0e-05,
                          rho.tol    = 1.0e-04,
                          iter.max   = 100L,
                          trace      = 0L,
                          drop.pc    = TRUE) {
  # Check validity of tolerance parameters
  if (step.tol <= 0.0 || dev.tol <= 0.0 || pseudo.tol <= 0.0 || rho.tol <= 0.0) {
    stop("All tolerance paramerters should be greater than zero.")
  }
  
  # Check validity of 'iter.max'
  if (iter.max < 1L) {
    stop("Maximum number of iterations should be at least one.")
  }
  
  # Return list with control parameters
  list(dev.tol    = dev.tol,
       step.tol   = step.tol,
       pseudo.tol = pseudo.tol,
       rho.tol    = rho.tol,
       iter.max   = as.integer(iter.max),
       trace      = as.integer(trace),
       drop.pc    = as.logical(drop.pc))
}
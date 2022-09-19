#' @title
#' Set \code{feglm} Control Parameters
#' @description
#' Set and change parameters used for fitting \code{\link{feglm}}.
#' 
#' \strong{Note:} \code{\link{feglm.control}} is deprecated and will be removed soon.
#' @param
#' dev.tol tolerance level for the first stopping condition of the maximization routine. The 
#' stopping condition is based on the relative change of the deviance in iteration \eqn{r}
#' and can be expressed as follows:
#' \eqn{|dev_{r} - dev_{r - 1}| / (0.1 + |dev_{r}|) < tol}{|dev - devold| / (0.1 + |dev|) < tol}.
#' Default is \code{1.0e-08}.
#' @param
#' center.tol tolerance level for the stopping condition of the centering algorithm.
#' The stopping condition is based on the relative change of the centered variable similar to
#' \link[lfe]{felm}. Default is \code{1.0e-08}.
#' @param
#' iter.max unsigned integer indicating the maximum number of iterations in the maximization
#' routine. Default is \code{25L}.
#' @param
#' limit unsigned integer indicating the maximum number of iterations of 
#' \code{\link[MASS]{theta.ml}}. Default is \code{10L}.
#' @param
#' trace logical indicating if output should be produced in each iteration. Default is \code{FALSE}.
#' @param
#' drop.pc logical indicating to drop observations that are perfectly classified/separated and hence 
#' do not contribute to the log-likelihood. This option is useful to reduce the computational costs of
#' the maximization problem and improves the numerical stability of the algorithm. Note that dropping
#' perfectly separated observations does not affect the estimates. Default is \code{TRUE}.
#' @param
#' keep.mx logical indicating if the centered regressor matrix should be stored. The centered regressor
#' matrix is required for some covariance estimators, bias corrections, and average partial effects. This
#' option saves some computation time at the cost of memory. Default is \code{TRUE}.
#' @param
#' conv.tol,rho.tol deprecated; step-halving is now similar to \code{glm.fit2}.
#' @param
#' pseudo.tol deprecated; use \code{center.tol} instead.
#' @param
#' step.tol deprecated; termination conditions is now similar to \code{\link[stats]{glm}}.
#' @param
#' ... arguments passed to the deprecated function \code{\link{feglm.control}}.
#' @return
#' The function \code{\link{feglmControl}} returns a named list of control 
#' parameters.
#' @seealso
#' \code{\link{feglm}}
#' @export
feglmControl <- function(
  dev.tol    = 1.0e-08,
  center.tol = 1.0e-08,
  iter.max   = 25L,
  limit      = 10L,
  trace      = FALSE,
  drop.pc    = TRUE,
  keep.mx    = TRUE,
  conv.tol   = NULL,
  rho.tol    = NULL,
  pseudo.tol = NULL,
  step.tol   = NULL
  ) {
  # 'conv.tol' is deprecated
  if (!is.null(conv.tol)) {
    warning("'conv.tol' is deprecated;", call. = FALSE)
  }
  
  # 'rho.tol' is deprecated
  if (!is.null(rho.tol)) {
    warning("'rho.tol' is deprecated;", call. = FALSE)
  }
  
  # 'pseudo.tol' is deprecated
  if (!is.null(pseudo.tol)) {
    warning("'pseudo.tol' is deprecated; please use 'center.tol' instead.", call. = FALSE)
    center.tol <- pseudo.tol
  }
  
  # 'step.tol' is deprecated
  if (!is.null(step.tol)) {
    warning("'step.tol' is deprecated;", call. = FALSE)
  }
  
  # Check validity of tolerance parameters
  if (dev.tol <= 0.0 || center.tol <= 0.0) {
    stop("All tolerance paramerters should be greater than zero.", call. = FALSE)
  }
  
  # Check validity of 'iter.max'
  iter.max <- as.integer(iter.max)
  if (iter.max < 1L) {
    stop("Maximum number of iterations should be at least one.", call. = FALSE)
  }
  
  # Check validity of 'limit'
  limit <- as.integer(limit)
  if (limit < 1L) {
    stop("Maximum number of iterations should be at least one.", call. = FALSE)
  }
  
  # Return list with control parameters
  list(
    dev.tol    = dev.tol,
    center.tol = center.tol,
    iter.max   = iter.max,
    limit      = limit,
    trace      = as.logical(trace),
    drop.pc    = as.logical(drop.pc),
    keep.mx    = as.logical(keep.mx)
    )
}


### Deprecated functions

#' @rdname feglmControl
#' @aliases feglmControl
#' @export
feglm.control <- function(...) {
  .Deprecated("feglmControl")
  do.call(feglmControl, list(...))
}

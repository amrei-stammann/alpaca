#' @title
#' Efficiently recover estimates of the fixed effects after fitting \code{feglm}
#' @description
#' Recover estimates of the fixed effects by alternating between the normal equations of the fixed 
#' effects as shown by Stammann (2018).
#' 
#' \strong{Remark}: The system might not have a unique solution since we do not take collinearity
#' into account. If the solution is not unique, an estimable function has to be applied to our
#' solution to get meaningful estimates of the fixed effects. See Gaure (n. d.) for an extensive
#' treatment of this issue.
#' @param
#' object an object of class \code{"feglm"}.
#' @param 
#' alpha.tol tolerance level for the stopping condition. The algorithm is stopped in iteration
#' \eqn{i} if \eqn{\lvert\boldsymbol{\alpha}_{i} - \boldsymbol{\alpha}_{i - 1}\rvert_{2} < 
#' \text{tol} \lvert\boldsymbol{\alpha}_{i - 1}\rvert_{2}}{||\Delta \alpha|| < tol ||\alpha_old||}.
#' Default is \code{1.0e-08}.
#' @return
#' The function \code{getFEs} returns a named vector of estimated fixed effects.
#' @references
#' Gaure, S. (n. d.). "Multicollinearity, identification, and estimable 
#' functions". Unpublished.
#' @references 
#' Stammann, A. (2018). "Fast and Feasible Estimation of Generalized Linear Models with 
#' High-Dimensional k-way Fixed Effects". Working Paper.
#' @seealso
#' \code{\link{feglm}}
#' @export
getFEs <- function(object, alpha.tol = 1.0e-08) {
  # Check validity of 'object'
  if(!inherits(object, "feglm")) {
    stop("'getFEs' called on a non-'feglm' object.")
  }
  
  # Construct auxiliary matrix to flatten the fixed effects
  lvls.k <- object[["lvls.k"]]
  k.vars <- names(lvls.k)
  fe <- object[["data"]][, k.vars, with = FALSE]
  fe[, (k.vars) := lapply(.SD, function(x) as.integer(factor(x)) - 1L)]
  A <- as.matrix(fe)
  dimnames(A) <- NULL
  B <- apply(A, 2L, order) - 1L
  rm(fe)
  
  # Recover fixed effects by alternating between normal equations
  alpha <- as.vector(get.alpha(object[["D.alpha"]], lvls.k, A, B, alpha.tol))
  names(alpha) <- object[["nms.fe"]]
  
  # Return estimates of the fixed effects
  alpha
}
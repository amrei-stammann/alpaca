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
#' \eqn{i} if \eqn{||\boldsymbol{\alpha}_{i} - \boldsymbol{\alpha}_{i - 1}||_{2} < 
#' tol ||\boldsymbol{\alpha}_{i - 1}||_{2}}{||\Delta \alpha|| < tol ||\alpha_old||}.
#' Default is \code{1.0e-08}.
#' @return
#' The function \code{\link{getFEs}} returns a named list containing named vectors of estimated 
#' fixed effects.
#' @references
#' Gaure, S. (n. d.). "Multicollinearity, identification, and estimable 
#' functions". Unpublished.
#' @references 
#' Stammann, A. (2018). "Fast and Feasible Estimation of Generalized Linear Models with 
#' High-Dimensional k-way Fixed Effects". ArXiv e-prints.
#' @seealso
#' \code{\link{feglm}}
#' @export
getFEs <- function(object = NULL, alpha.tol = 1.0e-08) {
  # Check validity of 'object'
  if (is.null(object)) {
    stop("'object' has to be specified.", call. = FALSE)
  } else if (!inherits(object, "feglm")) {
    stop("'getFEs' called on a non-'feglm' object.", call. = FALSE)
  }
  
  # Extract regressor matrix
  X <- model.matrix(object[["formula"]], object[["data"]], rhs = 1L)[, - 1L, drop = FALSE]
  nms.sp <- attr(X, "dimnames")[[2L]]
  attr(X, "dimnames") <- NULL
  
  # Construct auxilliary matrix to flatten the fixed effects
  lvls.k <- object[["lvls.k"]]
  k.vars <- names(lvls.k)
  fe <- object[["data"]][, k.vars, with = FALSE]
  fe[, (k.vars) := lapply(.SD, as.integer)]
  A <- as.matrix(fe) - 1L
  dimnames(A) <- NULL
  rm(fe)
  B <- apply(A, 2L, order) - 1L
  
  # Recover fixed effects by alternating between the solution of normal equations
  pi <- object[["eta"]] - as.vector(X %*% object[["coefficients"]])
  alpha <- as.vector(getAlpha(pi, lvls.k, A, B, alpha.tol))
  
  # Generate a list with one vector for each fixed effects category
  k <- length(lvls.k)
  tmp.var <- c(0L, cumsum(lvls.k))
  fe <- vector("list", k)
  for (i in seq.int(k)) {
    fe[[i]] <- alpha[seq(tmp.var[[i]] + 1L, tmp.var[[i + 1L]])]
    names(fe[[i]]) <- object[["nms.fe"]][[i]]
  }
  names(fe) <- k.vars
  
  # Return list of fixed effects estimates
  fe
}
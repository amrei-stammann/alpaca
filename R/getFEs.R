#' @title
#' getFEs: A function to efficiently recover the estimated fixed effects after
#' fitting \code{feglm}
#' @description
#' Recover the estimated fixed effects using the Kaczmarz method 
#' (Kaczmarz (1937)). This method is an iterative solver of a sparse system of
#' equations based on the method of alternating projections.
#' 
#' \strong{Remark}: The system might not have a unique solution since we do not
#' take collinearity into account. If the solution is not unique, an estimable
#' function has to be applied to our solution to get meaningful estimates of the
#' fixed effects . See Gaure (n. d.) for an extensive treatment of this issue.
#' 
#' \strong{Caution}: This may take a while depending on the number of
#' observations, the number of categories \eqn{k}, and the structure of the
#' underlying graph.
#' @param
#' obj an object of class \code{"alpaca"}.
#' @param 
#' alpha.tol tolerance level for the stopping condition. The algorithm is
#' stopped in iteration \eqn{i} if 
#' \eqn{||\boldsymbol{\alpha}_{i} - \boldsymbol{\alpha}_{i - 1}||_{2} < 
#' \text{tol}}{||step|| < tol}. Default is \code{1.0e-08}.
#' @param 
#' trace unsigned integer indicating if output should be produced for each
#' iteration. Default is \code{0}. See \code{details}.
#' @details 
#' ...
#' @return
#' The function \code{getFEs} returns a named vector of estimated fixed effects.
#' @references 
#' Kaczmarz, S. (1937). "Angenaeherte Aufloesung von Systemen linearer
#' Gleichungen". Bulletin International de l'Academie Polonaise des Sciences et
#' des Lettres. Classe des Sciences Mathematiques et Naturelles. Serie A,
#' Sciences Mathematiques. 35.
#' @references
#' Gaure, S. (n. d.). "Multicollinearity, identification, and estimable 
#' functions". Unpublished.
#' @seealso
#' \code{\link{feglm}}
#' @export
getFEs <- function(obj,
                   alpha.tol = 1.0e-08,
                   trace = 0L) {
  # Check validity of 'obj'.
  if(!inherits(obj, "alpaca")) {
    stop("'getFEs' called on a non-'alpaca' object.")
  }
  
  # Construct auxiliary matrix (needed for Kaczmarz).
  lvls.k <- c(0L, obj[["lvls.k"]])
  A <- obj[["data"]][["D"]] - 1L
  for (k in seq(ncol(A))) {
    A[, k] <- A[, k] + sum(lvls.k[seq(k)])
  }
  
  # Recover fixed effects using Kaczmarz and return.
  alpha <- as.vector(.kaczmarz(obj[["b"]], obj[["w.tilde"]], A, sum(lvls.k),
                               alpha.tol, trace))
  names(alpha) <- obj[["nms.fe"]]
  alpha
}
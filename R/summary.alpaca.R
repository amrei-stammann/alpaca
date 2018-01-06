#' @title
#' Summarizing Models of Class \code{alpaca}
#' @description
#' Summary statistics for objects of class \code{alpaca}.
#' @param 
#' object an object of class \code{alpaca}.
#' @param 
#' ... other arguments.
#' @return
#' Returns an object of class \code{summary.alpaca} which is a list of summary
#' statistics of \code{object}.
#' @seealso
#' \code{\link{feglm}}
#' @export
summary.alpaca <- function(object, ...) {
  # Compute coefficent matrix.
  cm.header <- c("Estimate", "Std. error", "z value", "Pr(> |z|)")
  beta.hat <- coef(object)
  se.beta <- sqrt(diag(vcov(object)))
  z.score <- beta.hat / se.beta
  p.value <- 2.0 * pnorm(- abs(z.score))
  cm <- cbind(beta.hat, se.beta, z.score, p.value)  
  rownames(cm) <- names(beta.hat)
  colnames(cm) <- cm.header
  
  # Return list.
  structure(list(cm = cm, 
                 maximum = object[["maximum"]],
                 deviance = object[["deviance"]],
                 nobs = object[["nobs"]],
                 lvls.k = object[["lvls.k"]],
                 nobs.na = object[["nobs.na"]],
                 nobs.pc = object[["nobs.pc"]],
                 formula = object[["formula"]],
                 family = object[["family"]]),
            class = "summary.alpaca")
}
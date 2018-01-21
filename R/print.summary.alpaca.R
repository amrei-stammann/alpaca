#' @title
#' Print \code{summary.alpaca}
#' @description
#' \code{print.summary.alpaca} is a generic function which displays summary
#' statistics from objects returned by \code{summary.alpaca}.
#' @param 
#' x an object of class \code{summary.alpaca}.
#' @param 
#' digits unsigned integer indicating the number of decimal places. Default is 
#' \code{max(3L, getOption("digits") - 3L)}.
#' @param 
#' ... other arguments.
#' @seealso
#' \code{\link{feglm}}
#' @export
print.summary.alpaca <- function(x,
                                 digits = max(3L, getOption("digits") - 3L),
                                 ...) {
  cat(x[["family"]], "\n\n")
  print(x[["formula"]])
  cat("\nn= ", x[["nobs"]],
      ", l= [", paste0(x[["lvls.k"]], collapse = ", "), "],\n",
      "log-likelihood= ", round(x[["maximum"]], digits = digits), ",\n",
      "deviance= ", round(x[["deviance"]], digits = digits), ",\n",
      sep = "")
  cat("\nStructural parameter(s):\n\n")
  printCoefmat(x[["cm"]], P.values = TRUE, has.Pvalue = TRUE,
               digits = digits)
  if(x[["nobs.na"]] > 0L) {
    cat("(", x[["nobs.na"]], "observation(s) deleted due to missingness )\n")
  }
  if(x[["nobs.pc"]] > 0L) {
    cat("(", x[["nobs.pc"]],
        "observation(s) deleted due to perfect classification )\n")
  }
}
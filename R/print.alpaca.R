#' @title
#' Print \code{alpaca}
#' @description
#' \code{print.alpaca} is a generic function which displays some minimal
#' information from objects returned by \code{feglm}.
#' @param 
#' x an object of class \code{alpaca}.
#' @param 
#' digits unsigned integer indicating the number of decimal places. Default is 
#' \code{max(3L, getOption("digits") - 3L)}.
#' @param 
#' ... other arguments.
#' @seealso
#' \code{\link{feglm}}
#' @export
print.alpaca <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat(x[["family"]], ", l= ", paste0(x[["lvls.k"]], collapse = ", ")
      , "\n\n", sep = "")
  print(coef(x), digits = digits)
}
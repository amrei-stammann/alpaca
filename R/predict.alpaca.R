#' @title
#' Predict Method for \code{feglm} Fits
#' @description
#' \code{predict.alpaca} is a generic function which obtains predictions from
#' an object returned by \code{feglm}.
#' @param 
#' object an object of class \code{alpaca}.
#' @param
#' type the type of prediction required. \code{"link"} is on the scale to the
#' linear predictor whereas \code{"response"} is on the scale of the response
#' variable. Default is \code{"link"}.
#' @param 
#' ... other arguments.
#' @return
#' The function \code{predict.alpaca} returns a vector of predictions
#' @seealso
#' \code{\link{feglm}}
#' @export
predict.alpaca <- function(object, type = c("link", "response"), ...) {
  # Check validity of 'type'.
  type <- match.arg(type)
  
  # Compute requested type of prediction.
  x <- as.vector(object[["data"]][["X"]] %*% coef(object) + object[["D.alpha"]])
  if (type == "response") {
    if (object[["family"]] == "logit") {
      x <- plogis(x)
    } else if (object[["family"]] == "poisson") {
      x <- exp(x)
    }
  }
  
  # Return prediction.
  x
}
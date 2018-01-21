#' @title
#' Generate a Simulated Data Set
#' @description
#' Constructs a simulated data set with \eqn{n} cross-sectional units observed
#' for \eqn{t} time periods according to the following data generating process:
#' 
#' \deqn{y_{it} = \mathbf{1}[\mathbf{x}_{it}^{\prime} \boldsymbol{\beta} +
#'  \alpha_{i} + \gamma_{t} + \epsilon_{it} > 0]\, ,}{y = 1[X \beta + \alpha + 
#'  \gamma + \epsilon],}
#' where \eqn{\mathbf{X}}{X} consists of 3 independent standard normally
#' distributed regressors and \eqn{\epsilon_{it}}{\epsilon} is an iid. logistic
#' error term with location zero and scale one. \eqn{\alpha_{i}}{\alpha} and 
#' \eqn{\gamma_{t}}{\gamma} are generated as iid. standard normal and 
#' \eqn{\boldsymbol{\beta} = [1, - 1, 1]^{\prime}}{\beta = [1, - 1, 1]'}.
#' @param
#' n number of cross-sectional units.
#' @param
#' t number of time periods.
#' @param
#' seed a seed to ensure reproducibility.
#' @return
#' The function \code{simData} returns a data.frame with 6 variables.
#' @export
simData <- function(n, t, seed) {
  set.seed(seed)
  beta <- c(1.0, - 1.0, 1.0)
  X.it <- matrix(rnorm(3 * n * t), nrow = n * t, ncol = 3)
  colnames(X.it) <- c("x1", "x2", "x3")
  alpha.i <- rnorm(n)
  gamma.t <- rnorm(t)
  epsilon.it <- rlogis(n * t)
  y.it <- X.it %*% beta + rep(alpha.i, each = t) + rep(gamma.t, times = n) + 
    epsilon.it > 0.0
  data.frame(y = as.integer(y.it), X.it, i = rep(seq(n), each = t),
             t = rep(seq(t), times = n))
}
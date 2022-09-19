#' @title
#' Generate an artificial data set for some GLM's with two-way fixed effects
#' @description
#' Constructs an artificial data set with \eqn{n} cross-sectional units observed for \eqn{t} time
#' periods for logit, poisson, or gamma models. The \dQuote{true} linear predictor
#' (\eqn{\boldsymbol{\eta}}{\eta}) is generated as follows:
#' \deqn{\eta_{it} = \mathbf{x}_{it}^{\prime} \boldsymbol{\beta} +
#'  \alpha_{i} + \gamma_{t} \, ,}{\eta = X \beta + \alpha + 
#'  \gamma,}
#' where \eqn{\mathbf{X}}{X} consists of three independent standard normally distributed regressors.
#' Both parameter referring to the unobserved heterogeneity (\eqn{\alpha_{i}}{\alpha} and 
#' \eqn{\gamma_{t}}{\gamma}) are generated as iid. standard normal and the structural parameters are
#' set to \eqn{\boldsymbol{\beta} = [1, - 1, 1]^{\prime}}{\beta = [1, - 1, 1]'}.
#' 
#' \strong{Note:} The poisson and gamma model are based on the logarithmic link function.
#' @param
#' n a strictly positive integer equal to the number of cross-sectional units.
#' @param
#' t a strictly positive integer equal to the number of time periods.
#' @param
#' seed a seed to ensure reproducibility.
#' @param
#' model a string equal to \code{"logit"}, \code{"poisson"}, or \code{"gamma"}. Default is 
#' \code{"logit"}.
#' @return
#' The function \code{\link{simGLM}} returns a data.frame with 6 variables.
#' @seealso
#' \code{\link{feglm}}
#' @export
simGLM <- function(
  n     = NULL,
  t     = NULL,
  seed  = NULL,
  model = c("logit", "poisson", "gamma")
  ) {
  # Validity check 'n'
  if (is.null(n)) {
    stop("'n' has to be specified.")
  } else {
    if (n < 1L) {
      stop("Number of cross-sectional units should be at least one.", call. = FALSE)
    }
    n <- as.integer(n)
  }
  
  # Validity check 't'
  if (is.null(t)) {
    stop("'t' has to be specified.")
  } else {
    if (t < 1L) {
      stop("Number of time periods should be at least one.", call. = FALSE)
    }
    t <- as.integer(t)
  }
  
  # Validity check 'seed'
  if (is.null(seed)) {
    stop("'seed' has to be specified.", call. = FALSE)
  } else {
    if (seed < 0L) {
      stop("'seed' has to be a positive integer.", call. = FALSE)
    }
    seed <- as.integer(seed)
  }
  
  # Match 'model'
  model <- match.arg(model)
  
  # Generate data
  set.seed(seed)
  beta <- c(1.0, - 1.0, 1.0)
  X <- matrix(rnorm(3L * n * t), n * t, 3L)
  colnames(X) <- c("x1", "x2", "x3")
  alpha <- rnorm(n)
  gamma <- rnorm(t)
  eta <- as.vector(X %*% beta + rep(alpha, each = t) + rep(gamma, n))
  if (model == "logit") {
    y <- as.integer(rlogis(n * t, eta) >= 0.5)
  } else if (model == "gamma") {
    y <- rgamma(n * t, 10.0, 10.0 / exp(eta))
  } else {
    y <- rpois(n * t, exp(eta))
  }
  
  # Return data.frame
  data.frame(
    i = rep(seq.int(n), each = t),
    t = rep.int(seq.int(t), n),
    y, X
    )
}
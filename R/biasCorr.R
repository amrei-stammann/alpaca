#' @title
#' Asymptotic bias-correction after fitting binary choice models with two-way error component
#' @description
#' \code{\link{biasCorr}} is a post-estimation routine that can be used to substantially reduce the 
#' incidental parameter bias problem (Neyman and Scott (1948)) present in non-linear fixed effects 
#' models (see Fernandez-Val and Weidner (2018) for an overview). The command applies the analytical 
#' bias-correction derived by Fernandez-Val and Weinder (2016) to obtain bias-corrected estimates of 
#' the structural parameters and is currently restricted to logit and probit models.
#' @param
#' object an object of class \code{"feglm"}; currently restricted to \code{\link[stats]{binomial}} 
#' with \code{"logit"} or \code{"probit"} link function.
#' @param
#' L unsigned integer indicating a bandwidth for the estimation of spectral densities proposed by 
#' Hahn and Kuersteiner (2011). Default is zero, which should be used if all regressors are 
#' assumed to be strictly exogenous. In the presence of weakly exogenous or predetermined 
#' regressors, Fernandez-Val and Weidner (2016, 2018) suggest to choose a bandwidth not higher than 
#' four.
#' @return
#' The function \code{\link{biasCorr}} returns a named list of classes \code{"biasCorr"} and 
#' \code{"feglm"}.
#' @references
#' Czarnowske, D. and Stammann, A. (2019). "Binary Choice Models with High-Dimensional Individual 
#' and Time Fixed Effects". ArXiv e-prints.
#' @references
#' Fernandez-Val, I. and Weidner, M. (2016). "Individual and time effects in nonlinear panel models 
#' with large N, T". Journal of Econometrics, 192(1), 291-312.
#' @references
#' Fernandez-Val, I. and Weidner M. (2018). "Fixed effects estimation of large-t panel data models". 
#' Annual Review of Economics, 10, 109-138.
#' @references
#' Hahn, J. and Kuersteiner, G. (2011). "Bias reduction for dynamic nonlinear panel models with 
#' fixed effects". Econometric Theory, 27(6), 1152-1191.
#' @references
#' Neyman, J. and Scott, E. L. (1948). "Consistent estimates based on partially consistent 
#' observations". Econometrica, 16(1), 1-32.
#' @seealso
#' \code{\link{feglm}}
#' @examples 
#' \donttest{
#' # Generate an artificial data set for logit models
#' library(alpaca)
#' data <- simGLM(1000L, 20L, 1805L, model = "logit")
#' 
#' # Fit 'feglm()'
#' mod <- feglm(y ~ x1 + x2 + x3 | i + t, data)
#' 
#' # Apply analytical bias-correction
#' mod.bc <- biasCorr(mod)
#' summary(mod.bc)
#' 
#' # Compute bias-corrected average partial effects
#' mod.ape <- getAPEs(mod.bc)
#' summary(mod.ape)
#' }
#' @export
biasCorr <- function(object = NULL, L = 0L) {
  # Check validity of 'object'
  if (is.null(object)) {
    stop("'object' has to be specified.", call. = FALSE)
  } else if (!inherits(object, "feglm")) {
    stop("'biasCorr' called on a non-'feglm' object.", call. = FALSE)
  }
  
  # Check provided object
  # TODO: Add further bias-corrections
  if (object[["family"]][["family"]] != "binomial" |
      !(object[["family"]][["link"]] %in% c("logit", "probit")) |
      length(object[["lvls.k"]]) != 2L) {
    stop(paste("'biasCorr' currently only supports logit and probit models with 2-way error",
               "component. Further corrections will be added in the future."), call. = FALSE)
  }
  
  # Extract model information
  formula <- object[["formula"]]
  family <- object[["family"]]
  lvls.k <- object[["lvls.k"]]
  control <- object[["control"]]
  
  # Extract model response and regressor matrix
  y <- object[["data"]][[1L]]
  if (!all(y %in% c(0L, 1L))) {
    stop("'biasCorr' currently only supports 'true' binary choice models.", call. = FALSE)
  }
  X <- model.matrix(formula, object[["data"]], rhs = 1L)[, - 1L, drop = FALSE]
  nms.sp <- attr(X, "dimnames")[[2L]] # Saves memory
  attr(X, "dimnames") <- NULL
  
  # Construct auxiliary matrix to flatten the fixed effects
  k.vars <- names(lvls.k)
  fe <- model.part(formula, object[["data"]], rhs = 2L)
  fe[, (k.vars) := lapply(.SD, as.integer)]
  A <- as.matrix(fe) - 1L
  dimnames(A) <- NULL
  rm(fe)
  B <- apply(A, 2L, order) - 1L
  
  # Compute derivatives and weights
  eta <- object[["eta"]]
  mu <- family[["linkinv"]](eta)
  mu.eta <- family[["mu.eta"]](eta)
  if (family[["link"]] == "logit") {
    v <- y - mu
    w <- mu.eta
    z <- w * (1.0 - 2.0 * mu)
  } else {
    w <- mu.eta / family[["variance"]](mu)
    v <- w * (y - mu)
    w <- w * mu.eta
    z <- - eta * w
  }
  
  # Centering variables
  MX <- object[["Score"]] / v
  
  # Compute \hat{B}
  b <- as.vector(groupSums(MX * z, w, A[, 1L], B[, 1L])) / 2.0
  if (L > 0L) {
    b <- b + as.vector(groupSumsSpectral(MX * w, v, w, L, A[, 1L], B[, 1L]))
  }
  
  # Compute \hat{D}
  d <- as.vector(groupSums(MX * z, w, A[, 2L], B[, 2L])) / 2.0
  
  # Compute bias-corrected structural parameters
  beta <- object[["coefficients"]] - solve(object[["Hessian"]], - (b + d))
  names(beta) <- nms.sp
 
  # Update result list
  object[["coefficients"]] <- beta
  object[["bandwidth"]] <- L
  
  # Add additional class to result list
  attr(object, "class") <- c("feglm", "biasCorr")
  
  # Return updated list
  object
}


### Internal function (not exported)
# Higher-order partial derivatives
partialMuEta <- function(eta, family, order) {
  f <- family[["mu.eta"]](eta)
  if (order == 2L) {
    if (family[["link"]] == "logit") {
      f * (1.0 - 2.0 * family[["linkinv"]](eta))
    } else {
      - eta * f
    }
  } else {
    if (family[["link"]] == "logit") {
      f * ((1.0 - 2.0 * family[["linkinv"]](eta))^2 - 2.0 * f)
    } else {
      (eta^2 - 1.0) * f
    }
  }
}
#' @title
#' Asymptotic bias correction after fitting binary choice models with a one-/two-/three-way error
#' component
#' @description
#' \code{\link{biasCorr}} is a post-estimation routine that can be used to substantially reduce the 
#' incidental parameter bias problem (Neyman and Scott (1948)) present in nonlinear fixed effects 
#' models (see Fernández-Val and Weidner (2018) for an overview). The command applies the analytical 
#' bias correction derived by Fernández-Val and Weidner (2016) and Hinz, Stammann, and Wanner (2020) 
#' to obtain bias-corrected estimates of the structural parameters and is currently restricted to 
#' \code{\link[stats]{binomial}} with one-, two-, and three-way fixed effects.
#' @param
#' object an object of class \code{"feglm"}; currently restricted to \code{\link[stats]{binomial}}.
#' @param
#' L unsigned integer indicating a bandwidth for the estimation of spectral densities proposed by 
#' Hahn and Kuersteiner (2011). Default is zero, which should be used if all regressors are 
#' assumed to be strictly exogenous with respect to the idiosyncratic error term. In the presence of 
#' weakly exogenous regressors, e.g. lagged outcome variables, Fernández-Val and Weidner (2016, 2018) 
#' suggest to choose a bandwidth between one and four. Note that the order of factors to be partialed 
#' out is important for bandwidths larger than zero (see vignette for details).
#' @param
#' panel.structure a string equal to \code{"classic"} or \code{"network"} which determines the 
#' structure of the panel used. \code{"classic"} denotes panel structures where for example the same
#' cross-sectional units are observed several times (this includes pseudo panels). 
#' \code{"network"} denotes panel structures where for example bilateral trade flows are observed 
#' for several time periods. Default is \code{"classic"}.
#' @return
#' The function \code{\link{biasCorr}} returns a named list of classes \code{"biasCorr"} and 
#' \code{"feglm"}.
#' @references
#' Czarnowske, D. and A. Stammann (2020). "Fixed Effects Binary Choice Models: Estimation and Inference
#' with Long Panels". ArXiv e-prints.
#' @references
#' Fernández-Val, I. and M. Weidner (2016). "Individual and time effects in nonlinear panel models 
#' with large N, T". Journal of Econometrics, 192(1), 291-312.
#' @references
#' Fernández-Val, I. and M. Weidner (2018). "Fixed effects estimation of large-t panel data 
#' models". Annual Review of Economics, 10, 109-138.
#' @references
#' Hahn, J. and G. Kuersteiner (2011). "Bias reduction for dynamic nonlinear panel models with 
#' fixed effects". Econometric Theory, 27(6), 1152-1191.
#' @references
#' Hinz, J., A. Stammann, and J. Wanner (2020). "State Dependence and Unobserved Heterogeneity
#' in the Extensive Margin of Trade". ArXiv e-prints.
#' @references
#' Neyman, J. and E. L. Scott (1948). "Consistent estimates based on partially consistent 
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
#' # Apply analytical bias correction
#' mod.bc <- biasCorr(mod)
#' summary(mod.bc)
#' }
#' @export
biasCorr <- function(
    object          = NULL,
    L               = 0L,
    panel.structure = c("classic", "network")
  ) {
  # Check validity of 'object'
  if (is.null(object)) {
    stop("'object' has to be specified.", call. = FALSE)
  } else if (!inherits(object, "feglm")) {
    stop("'biasCorr' called on a non-'feglm' object.", call. = FALSE)
  }
  
  # Check validity of 'panel.structure'
  panel.structure <- match.arg(panel.structure)
  
  # Extract model information
  beta.uncorr <- object[["coefficients"]]
  control <- object[["control"]]
  data <- object[["data"]]
  eps <- .Machine[["double.eps"]]
  family <- object[["family"]]
  formula <- object[["formula"]]
  lvls.k <- object[["lvls.k"]]
  nms.sp <- names(beta.uncorr)
  nt <- object[["nobs"]][["nobs"]]
  k.vars <- names(lvls.k)
  k <- length(lvls.k)
  
  # Check if binary choice model
  if (family[["family"]] != "binomial") {
    stop("'biasCorr' currently only supports binary choice models.", call. = FALSE)
  }
  
  # Check if provided object matches requested panel structure
  if (panel.structure == "classic") {
    if (!(k %in% c(1L, 2L))) {
      stop("panel.structure == 'classic' expects a one- or two-way fixed effects model.", call. = FALSE)
    }
  } else {
    if (!(k %in% c(2L, 3L))) {
      stop("panel.structure == 'network' expects a two- or three-way fixed effects model.", call. = FALSE)
    }
  }
  
  # Extract model response, regressor matrix, and weights
  y <- data[[1L]]
  X <- model.matrix(formula, data, rhs = 1L)[, - 1L, drop = FALSE]
  attr(X, "dimnames") <- NULL
  wt <- object[["weights"]]
  
  # Generate auxiliary list of indexes for different sub panels
  k.list <- getIndexList(k.vars, data)
  
  # Compute derivatives and weights
  eta <- object[["eta"]]
  mu <- family[["linkinv"]](eta)
  mu.eta <- family[["mu.eta"]](eta)
  v <- wt * (y - mu)
  w <- wt * mu.eta
  z <- wt * partialMuEta(eta, family, 2L)
  if (family[["link"]] != "logit") {
    h <- mu.eta / family[["variance"]](mu)
    v <- h * v
    w <- h * w
    z <- h * z
    rm(h)
  }
  
  # Center regressor matrix (if required)
  if (control[["keep.mx"]]) {
    MX <- object[["MX"]]
  } else {
    MX <- centerVariables(X, w, k.list, control[["center.tol"]])
  }
  
  # Compute bias terms for requested bias correction
  if (panel.structure == "classic") {
    # Compute \hat{B} and \hat{D}
    b <- as.vector(groupSums(MX * z, w, k.list[[1L]])) / 2.0 / nt
    if (k > 1L) {
      b <- b + as.vector(groupSums(MX * z, w, k.list[[2L]])) / 2.0 / nt
    }
    
    # Compute spectral density part of \hat{B}
    if (L > 0L) {
      b <- b + as.vector(groupSumsSpectral(MX * w, v, w, L, k.list[[1L]])) / nt
    }
  } else {
    # Compute \hat{D}_{1}, \hat{D}_{2}, and \hat{B}
    b <- as.vector(groupSums(MX * z, w, k.list[[1L]])) / 2.0 / nt
    b <- b + as.vector(groupSums(MX * z, w, k.list[[2L]])) / 2.0 / nt
    if (k > 2L) {
      b <- b + as.vector(groupSums(MX * z, w, k.list[[3L]])) / 2.0 / nt
    }
    
    # Compute spectral density part of \hat{B}
    if (k > 2L && L > 0L) {
      b <- b + as.vector(groupSumsSpectral(MX * w, v, w, L, k.list[[3L]])) / nt
    }
  }
  
  # Compute bias-corrected structural parameters
  b <- solve(object[["Hessian"]] / nt, - b)
  beta <- beta.uncorr - b
  names(beta) <- nms.sp
  
  # Update \eta and first- and second-order derivatives
  eta <- feglmOffset(object, as.vector(X %*% beta))
  mu <- family[["linkinv"]](eta)
  mu.eta <- family[["mu.eta"]](eta)
  v <- wt * (y - mu)
  w <- wt * mu.eta
  if (family[["link"]] != "logit") {
    h <- mu.eta / family[["variance"]](mu)
    v <- h * v
    w <- h * w
    rm(h)
  }
  
  # Update centered regressor matrix
  MX <- centerVariables(X, w, k.list, control[["center.tol"]])
  colnames(MX) <- nms.sp
  
  # Update Hessian
  H <- crossprod(MX * sqrt(w))
  dimnames(H) <- list(nms.sp, nms.sp)
  
  # Update result list
  object[["coefficients"]] <- beta
  object[["eta"]] <- eta
  if (control[["keep.mx"]]) object[["MX"]] <- MX
  object[["Hessian"]] <- H
  object[["coefficients.uncorr"]] <- beta.uncorr
  object[["bias.term"]] <- b
  object[["panel.structure"]] <- panel.structure
  object[["bandwidth"]] <- L
  
  # Add additional class to result list
  attr(object, "class") <- c("feglm", "biasCorr")
  
  # Return updated list
  object
}
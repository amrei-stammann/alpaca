#' @title
#' Asymptotic bias correction after fitting binary choice models with a two-/three-way error
#' component
#' @description
#' \code{\link{biasCorr}} is a post-estimation routine that can be used to substantially reduce the 
#' incidental parameter bias problem (Neyman and Scott (1948)) present in non-linear fixed effects 
#' models (see Fernández-Val and Weidner (2018) for an overview). The command applies the analytical 
#' bias correction derived by Fernández-Val and Weinder (2016) and Hinz, Stammann, and Wanner (2020) 
#' to obtain bias-corrected estimates of the structural parameters and is currently restricted to 
#' logit and probit models with two- and three-way fixed effects.
#' @param
#' object an object of class \code{"feglm"}; currently restricted to \code{\link[stats]{binomial}} 
#' with \code{"logit"} or \code{"probit"} link function.
#' @param
#' L unsigned integer indicating a bandwidth for the estimation of spectral densities proposed by 
#' Hahn and Kuersteiner (2011). Default is zero, which should be used if all regressors are 
#' assumed to be strictly exogenous with respect to the error term. In the presence of weakly exogenous 
#' or predetermined regressors, Fernández-Val and Weidner (2016, 2018) suggest to choose a bandwidth 
#' zero and four. Note that the order of factors to be partialed out is important for bandwidths larger 
#' than zero (see vignette for details).
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
biasCorr <- function(object          = NULL,
                     L               = 0L,
                     panel.structure = c("classic", "network")) {
  # TODO: Add further bias corrections
  # Check validity of 'object'
  if (is.null(object)) {
    stop("'object' has to be specified.", call. = FALSE)
  } else if (!inherits(object, "feglm")) {
    stop("'biasCorr' called on a non-'feglm' object.", call. = FALSE)
  }
  
  # Check validity of 'panel.structure'
  panel.structure <- match.arg(panel.structure)
  
  # Extract model information
  control <- object[["control"]]
  data <- object[["data"]]
  family <- object[["family"]]
  formula <- object[["formula"]]
  lvls.k <- object[["lvls.k"]]
  k.vars <- names(lvls.k)
  k <- length(lvls.k)
  
  # Check if binary choice model
  if (family[["family"]] != "binomial" || !(family[["link"]] %in% c("logit", "probit"))) {
    stop(paste("'biasCorr' currently only supports logit and probit models."), call. = FALSE)
  }
  
  # Check if provided object matches requested panel structure
  if (panel.structure == "classic") {
    if (length(lvls.k) != 2L) {
      stop(paste("panel.structure == 'classic' expects a two-way fixed effects model."), call. = FALSE)
    }
  } else {
    if (!(length(lvls.k) %in% c(2L, 3L))) {
      stop(paste("panel.structure == 'network' expects a two- or three-way fixed effects model."), call. = FALSE)
    }
  }
  
  # Extract model response and regressor matrix
  y <- data[[1L]]
  if (!all(y %in% c(0L, 1L))) {
    stop("'biasCorr' currently only supports models with binary outcome.", call. = FALSE)
  }
  X <- model.matrix(formula, data, rhs = 1L)[, - 1L, drop = FALSE]
  nms.sp <- attr(X, "dimnames")[[2L]]
  attr(X, "dimnames") <- NULL
  
  # Generate auxiliary list of indexes for different sub panels
  k.list <- getIndexList(k.vars, data)
  
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
  
  # Compute bias terms for requested bias correction
  if (panel.structure == "classic") {
    # Compute \hat{B}
    b <- as.vector(groupSums(MX * z, w, k.list[[1L]])) / 2.0
    if (L > 0L) {
      b <- b + as.vector(groupSumsSpectral(MX * w, v, w, L, k.list[[1L]]))
    }
    
    # Compute \hat{D}
    b <- b + as.vector(groupSums(MX * z, w, k.list[[2L]])) / 2.0
  } else {
    # Compute \hat{D}_{1}
    b <- as.vector(groupSums(MX * z, w, k.list[[1L]])) / 2.0
    
    # Compute \hat{D}_{2}
    b <- b + as.vector(groupSums(MX * z, w, k.list[[2L]])) / 2.0
    
    # Compute \hat{B}
    if (k == 3L) {
      b <- b + as.vector(groupSums(MX * z, w, k.list[[3L]])) / 2.0
      if (L > 0L) {
        b <- b + as.vector(groupSumsSpectral(MX * w, v, w, L, k.list[[3L]]))
      }
    }
  }
  
  # Compute bias-corrected structural parameters
  beta <- object[["coefficients"]] - solve(object[["Hessian"]], - b)
  names(beta) <- nms.sp
  
  # Update result list
  object[["coefficients"]] <- beta
  object[["panel.structure"]] <- panel.structure
  object[["bandwidth"]] <- L
  
  # Add additional class to result list
  attr(object, "class") <- c("feglm", "biasCorr")
  
  # Return updated list
  object
}
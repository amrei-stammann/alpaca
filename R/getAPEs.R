#' @title
#' Compute average partial effects after fitting binary choice models with two-way error component
#' @description
#' \code{\link{getAPEs}} is a post-estimation routine that can be used to estimate average partial 
#' effects with respect to all covariates in the model and the corresponding covariance matrix. The 
#' estimation of the covariance is based on a linear approximation (delta method). Note that 
#' the command automatically determines which of the regressors are continuous or binary.
#' 
#' \strong{Remark:} The routine currently does not allow to compute average partial effects based 
#' on functional forms like interactions and polynomials.
#' @param
#' object an object of class \code{"biasCorr"} or \code{"feglm"}; currently restricted to 
#' \code{\link[stats]{binomial}} with \code{"logit"} or \code{"probit"} link function.
#' @param
#' n.pop unsigned integer indicating a finite population correction for the estimation of the 
#' covariance matrix of the average partial effects proposed by 
#' Cruz-Gonzalez, Fernandez-Val, and Weidner (2017). The correction factor is computed as follows: 
#' \eqn{(n^{\ast} - n) / (n^{\ast} - 1)}{(n.pop - n) / (n.pop - 1)}, 
#' where \eqn{n^{\ast}}{n.pop} and \eqn{n}{n} are the size of the entire 
#' population and the full sample size. Default is \code{NULL}, which refers to a factor of one and 
#' is equal to an infinitely large population.
#' @param
#' weak.exo logical indicating if some of the regressors are assumed to be weakly exogenous (e.g. 
#' predetermined). If object is of class \code{"biasCorr"}, the option will be automatically set to 
#' \code{TRUE} if the choosen bandwidth parameter is larger than zero. Note that this option only 
#' affects the estimation of the covariance matrix. Default is \code{FALSE}, which assumes that all 
#' regressors are strictly exogenous.
#' @return
#' The function \code{\link{getAPEs}} returns a named list of class \code{"APEs"}.
#' @references
#' Cruz-Gonzalez, M., Fernandez-Val, I., and Weidner, M. (2017). "Bias corrections for probit and 
#' logit models with two-way fixed effects". The Stata Journal, 17(3), 517-545.
#' @references
#' Czarnowske, D. and Stammann, A. (2019). "Binary Choice Models with High-Dimensional Individual 
#' and Time Fixed Effects". ArXiv e-prints.
#' @references
#' Fernandez-Val, I. and Weidner, M. (2016). "Individual and time effects in nonlinear panel models 
#' with large N, T". Journal of Econometrics, 192(1), 291-312.
#' @references
#' Fernandez-Val, I. and Weidner, M. (2018). "Fixed effects estimation of large-t panel data 
#' models". Annual Review of Economics, 10, 109-138.
#' @references
#' Neyman, J. and Scott, E. L. (1948). "Consistent estimates based on partially consistent 
#' observations". Econometrica, 16(1), 1-32.
#' @seealso
#' \code{\link{biasCorr}}, \code{\link{feglm}}
#' @examples 
#' \donttest{
#' # Generate an artificial data set for logit models
#' library(alpaca)
#' data <- simGLM(1000L, 20L, 1805L, model = "logit")
#' 
#' # Fit 'feglm()'
#' mod <- feglm(y ~ x1 + x2 + x3 | i + t, data)
#' 
#' # Compute average partial effects
#' mod.ape <- getAPEs(mod)
#' summary(mod.ape)
#' 
#' # Apply analytical bias-correction
#' mod.bc <- biasCorr(mod)
#' summary(mod.bc)
#' 
#' # Compute bias-corrected average partial effects
#' mod.ape.bc <- getAPEs(mod.bc)
#' summary(mod.ape.bc)
#' }
#' @export
getAPEs <- function(object = NULL, n.pop = NULL, weak.exo = FALSE) {
  # Check validity of 'object'
  if (is.null(object)) {
    stop("'object' has to be specified.", call. = FALSE)
  } else if (!inherits(object, "feglm")) {
    stop("'getAPEs' called on a non-'feglm' object.", call. = FALSE)
  }
  
  # Check if provided object is a two-way logit or probit
  # TODO: Add further average partial effects
  if (object[["family"]][["family"]] != "binomial" |
      !(object[["family"]][["link"]] %in% c("logit", "probit")) |
      length(object[["lvls.k"]]) != 2L) {
    stop(paste0("'getAPEs' currently only supports logit and probit models with 2-way error ", 
                "component. Further models will be added in the future."), call. = FALSE)
  }
  
  # Extract model information
  beta <- object[["coefficients"]]
  formula <- object[["formula"]]
  family <- object[["family"]]
  lvls.k <- object[["lvls.k"]]
  control <- object[["control"]]
  nt.full <- object[["nobs"]][["nobs.full"]]
  nt <- object[["nobs"]][["nobs"]]
  p <- length(beta)
  
  # Check validity of 'n.pop'
  if (!is.null(n.pop)) {
    n.pop <- as.integer(n.pop)
    if (n.pop < nt.full) {
      warning(paste("Size of the entire population is lower than the full sample size.",
                    "Correction factor set to zero."), call. = FALSE)
      adj <- 0.0
    } else {
      adj <- (n.pop - nt.full) / (n.pop - 1L)
    }
  } else {
    adj <- 1.0
  }
  
  # Extract model response and regressor matrix
  y <- object[["data"]][[1L]]
  X <- model.matrix(formula, object[["data"]], rhs = 1L)[, - 1L, drop = FALSE]
  nms.sp <- attr(X, "dimnames")[[2L]] # Saves memory
  attr(X, "dimnames") <- NULL
  
  # Determine which of the regressors are binary
  binary <- apply(X, 2L, function(x) all(x %in% c(0L, 1L)))
  
  # Construct auxiliary matrix to flatten the fixed effects
  k.vars <- names(lvls.k)
  fe <- model.part(formula, object[["data"]], rhs = 2L)
  fe[, (k.vars) := lapply(.SD, as.integer)]
  A <- as.matrix(fe) - 1L
  dimnames(A) <- NULL
  rm(fe)
  B <- apply(A, 2L, order) - 1L
  
  # Compute projection of the regressor matrix
  eta <- object[["eta"]]
  mu <- family[["linkinv"]](eta)
  mu.eta <- family[["mu.eta"]](eta)
  if (family[["link"]] == "logit") {
    v <- y - mu
  } else {
    v <- mu.eta / family[["variance"]](mu) * (y - mu)
  }
  PX <- X - object[["Score"]] / v
  
  # If object is of class 'biasCorr' update linear predictor and the corresponding probabilities
  if (inherits(object, "biasCorr")) {
    eta <- as.vector(X %*% beta)
    eta <- feglmOffset(object, eta)
    mu <- family[["linkinv"]](eta)
    mu.eta <- family[["mu.eta"]](eta)
  }
  
  # Compute required derivatives and weights
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
  
  # Compute average partial effects and Jacobian
  Delta <- matrix(NA_real_, nt, p)
  Delta1 <- matrix(NA_real_, nt, p)
  J <- matrix(NA_real_, p, p)
  Delta[, !binary] <- mu.eta
  Delta1[, !binary] <- partialMuEta(eta, family, 2L)
  for (i in seq.int(p)) {
    if (binary[[i]]) {
      eta0 <- eta - X[, i] * beta[[i]]
      f1 <- family[["mu.eta"]](eta0 + beta[[i]])
      Delta[, i] <- family[["linkinv"]](eta0 + beta[[i]]) - family[["linkinv"]](eta0)
      Delta1[, i] <- f1 - family[["mu.eta"]](eta0)
      J[i, ] <- - sum(PX[, i] * Delta1[, i])
      J[i, i] <- sum(f1) + J[i, i]
      J[i, - i] <- colSums(Delta1[, i] * X[, - i, drop = FALSE]) + J[i, - i]
      rm(eta0, f1)
    } else {
      Delta[, i] <- beta[[i]] * Delta[, i]
      Delta1[, i] <- beta[[i]] * Delta1[, i]
      J[i, ] <- colSums((X - PX[, i]) * Delta1[, i])
      J[i, i] <- sum(mu.eta) + J[i, i]
    }
  }
  delta <- colSums(Delta) / nt.full
  Delta <- t(t(Delta) - delta)
  J <- J / nt.full
  rm(mu.eta, PX)
  
  # Compute projection and residual projection of \Psi
  MPsi <- centerVariables(Delta1 / sqrt(w), sqrt(w), A, B, lvls.k, control[["center.tol"]]) / 
    sqrt(w)
  PPsi <- Delta1 / w - MPsi
  rm(Delta1)
  
  # Compute analytical bias-correction of average partial effects
  if (inherits(object, "biasCorr")) {
    # Check validity of 'L' and 'weak.exo'
    L <- object[["bandwidth"]]
    if (L == 0L && weak.exo) {
      warning("Inconsistent choice of 'weak.exo'; argument set to FALSE.", call. = FALSE)
      weak.exo <- FALSE
    }
    
    # Compute second-order partial derivatives
    Delta2 <- matrix(NA_real_, nt, p)
    Delta2[, !binary] <- partialMuEta(eta, family, 3L)
    for (i in seq.int(p)) {
      if (binary[[i]]) {
        eta0 <- eta - X[, i] * beta[[i]]
        Delta2[, i] <- partialMuEta(eta0 + beta[[i]], family, 2L) - partialMuEta(eta0, family, 2L)
        rm(eta0)
      } else {
        Delta2[, i] <- beta[[i]] * Delta2[, i]
      }
    }
    
    # Compute \hat{B}
    b <- as.vector(groupSums(- PPsi * z + Delta2, w, A[, 1L], B[, 1L])) / 2.0
    if (L > 0L) {
      weak.exo <- TRUE
      b <- b + as.vector(groupSumsSpectral(MPsi * w, v, w, L, A[, 1L], B[, 1L]))
    }
    
    # Compute \hat{D}
    d <- as.vector(groupSums(- PPsi * z + Delta2, w, A[, 2L], B[, 2L])) / 2.0
    rm(Delta2)
    
    # Compute bias-corrected average partial effects
    delta <- delta - (b + d) / nt
  }
  rm(eta, mu, MPsi)
  
  # Compute standard errors
  J <- J %*% solve(object[["Hessian"]] / nt)
  Gamma <- t(tcrossprod(J, object[["Score"]])) - PPsi * v
  V <- crossprod(Gamma)
  if (adj > 0.0) {
    V <- V + adj * (groupSumsVar(Delta, A[, 1L], B[, 1L]) + groupSumsVar(Delta, A[, 2L], B[, 2L]) - 
                      crossprod(Delta))
    if (weak.exo) {
      V <- V + adj * 2.0 * groupSumsCov(Delta, Gamma, A[, 1L], B[, 1L])
    }
  }
  V <- V / nt^2
  
  # Add names
  names(delta) <- nms.sp
  dimnames(V) <- list(nms.sp, nms.sp)
  
  # Return list
  structure(list(delta = delta, vcov = V), class = "APEs")
}
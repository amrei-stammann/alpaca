#' @title
#' Compute average partial effects after fitting binary choice models with a two-/three-way error 
#' component
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
#' Cruz-Gonzalez, Fernández-Val, and Weidner (2017). The correction factor is computed as follows: 
#' \eqn{(n^{\ast} - n) / (n^{\ast} - 1)}{(n.pop - n) / (n.pop - 1)}, 
#' where \eqn{n^{\ast}}{n.pop} and \eqn{n}{n} are the size of the entire 
#' population and the full sample size. Default is \code{NULL}, which refers to a factor of one and 
#' is equal to an infinitely large population.
#' @param
#' panel.structure a string equal to \code{"classic"} or \code{"network"} which determines the 
#' structure of the panel used. \code{"classic"} denotes panel structures where for example the same
#' cross-sectional units are observed several times (this includes pseudo panels). 
#' \code{"network"} denotes panel structures where for example bilateral trade flows are observed 
#' for several time periods. Default is \code{"classic"}.
#' @param
#' sampling.fe a string equal to \code{"independence"} or \code{"unrestricted"} which imposes 
#' sampling assumptions about the unobserved effects. \code{"independence"} imposes that all 
#' unobserved effects are mutually independent sequences. \code{"unrestricted"} does not impose any
#' sampling assumptions. Note that this option only affects the estimation of the covariance. 
#' Default is \code{"independence"}.
#' @param
#' weak.exo logical indicating if some of the regressors are assumed to be weakly exogenous (e.g. 
#' predetermined). If object is of class \code{"biasCorr"}, the option will be automatically set to 
#' \code{TRUE} if the choosen bandwidth parameter is larger than zero. Note that this option only 
#' affects the estimation of the covariance matrix. Default is \code{FALSE}, which assumes that all 
#' regressors are strictly exogenous.
#' @return
#' The function \code{\link{getAPEs}} returns a named list of class \code{"APEs"}.
#' @references
#' Cruz-Gonzalez, M., I. Fernández-Val, and M. Weidner (2017). "Bias corrections for probit and 
#' logit models with two-way fixed effects". The Stata Journal, 17(3), 517-545.
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
#' Hinz, J., A. Stammann, and J. Wanner (2020). "State Dependence and Unobserved Heterogeneity
#' in the Extensive Margin of Trade". ArXiv e-prints.
#' @references
#' Neyman, J. and E. L. Scott (1948). "Consistent estimates based on partially consistent 
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
#' # Apply analytical bias correction
#' mod.bc <- biasCorr(mod)
#' summary(mod.bc)
#' 
#' # Compute bias-corrected average partial effects
#' mod.ape.bc <- getAPEs(mod.bc)
#' summary(mod.ape.bc)
#' }
#' @export
getAPEs <- function(object          = NULL,
                    n.pop           = NULL,
                    panel.structure = c("classic", "network"),
                    sampling.fe     = c("independence", "unrestricted"),
                    weak.exo        = FALSE) {
  # TODO: Add further average partial effects
  # Check validity of 'object'
  if (is.null(object)) {
    stop("'object' has to be specified.", call. = FALSE)
  } else if (!inherits(object, "feglm")) {
    stop("'getAPEs' called on a non-'feglm' object.", call. = FALSE)
  }
  
  # Extract prior information if available or check validity of 'panel.structure'
  biascorr <- inherits(object, "biasCorr")
  if (biascorr) {
    panel.structure <- object[["panel.structure"]]
    L <- object[["bandwidth"]]
    if (L > 0L) {
      weak.exo <- TRUE
    } else {
      weak.exo <- FALSE
    }
  } else {
    panel.structure <- match.arg(panel.structure)
  }
  
  # Check validity of 'sampling.fe'
  sampling.fe <- match.arg(sampling.fe)
  
  # Extract model information
  beta <- object[["coefficients"]]
  control <- object[["control"]]
  data <- object[["data"]]
  family <- object[["family"]]
  formula <- object[["formula"]]
  lvls.k <- object[["lvls.k"]]
  nt.full <- object[["nobs"]][["nobs.full"]]
  nt <- object[["nobs"]][["nobs"]]
  k <- length(lvls.k)
  k.vars <- names(lvls.k)
  p <- length(beta)
  
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
  y <- data[[1L]]
  X <- model.matrix(formula, data, rhs = 1L)[, - 1L, drop = FALSE]
  nms.sp <- attr(X, "dimnames")[[2L]]
  attr(X, "dimnames") <- NULL
  
  # Determine which of the regressors are binary
  binary <- apply(X, 2L, function(x) all(x %in% c(0.0, 1.0)))
  
  # Generate auxiliary list of indexes for different sub panels
  k.list <- getIndexList(k.vars, data)
  
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
  if (biascorr) {
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
      J[i, ] <- - colSums(PX * Delta1[, i])
      J[i, i] <- sum(f1) + J[i, i]
      J[i, - i] <- colSums(X[, - i, drop = FALSE] * Delta1[, i]) + J[i, - i]
      rm(eta0, f1)
    } else {
      Delta[, i] <- beta[[i]] * Delta[, i]
      Delta1[, i] <- beta[[i]] * Delta1[, i]
      J[i, ] <- colSums((X - PX) * Delta1[, i])
      J[i, i] <- sum(mu.eta) + J[i, i]
    }
  }
  delta <- colSums(Delta) / nt.full
  Delta <- t(t(Delta) - delta)
  rm(mu, mu.eta, PX)
  
  # Compute projection and residual projection of \Psi
  Psi <- - Delta1 / w
  rm(Delta1)
  MPsi <- centerVariables(Psi * sqrt(w), sqrt(w), k.list, control[["center.tol"]]) / sqrt(w)
  PPsi <- Psi - MPsi
  rm(Psi)
  
  # Compute analytical bias-correction of average partial effects
  if (biascorr) {
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
    
    # Compute bias terms for requested bias correction
    if (panel.structure == "classic") {
      # Compute \hat{B}
      b <- as.vector(groupSums(PPsi * z + Delta2, w, k.list[[1L]])) / 2.0
      if (weak.exo) {
        b <- b - as.vector(groupSumsSpectral(MPsi * w, v, w, L, k.list[[1L]]))
      }
      
      # Compute \hat{D}
      b <- b + as.vector(groupSums(PPsi * z + Delta2, w, k.list[[2L]])) / 2.0
    } else {
      # Compute \hat{D}_{1}
      b <- as.vector(groupSums(PPsi * z + Delta2, w, k.list[[1L]])) / 2.0
      
      # Compute \hat{D}_{2}
      b <- b + as.vector(groupSums(PPsi * z + Delta2, w, k.list[[2L]])) / 2.0
      
      # Compute \hat{B}
      if (k == 3L) {
        b <- b + as.vector(groupSums(PPsi * z + Delta2, w, k.list[[3L]])) / 2.0
        if (weak.exo) {
          b <- b - as.vector(groupSumsSpectral(MPsi * w, v, w, L, k.list[[3L]]))
        }
      }
    }
    rm(Delta2)
    
    # Compute bias-corrected average partial effects
    delta <- delta - b / nt.full
  }
  rm(eta, w, z, MPsi)
  
  # Compute covariance matrix
  J <- J %*% solve(object[["Hessian"]])
  Gamma <- tcrossprod(object[["Score"]], J) - PPsi * v
  rm(v, PPsi)
  V <- crossprod(Gamma)
  if (adj > 0.0) {
    # Simplify covariance if sampling assumptions are imposed
    if (sampling.fe == "independence") {
      V <- V + adj * (groupSumsVar(Delta, k.list[[1L]]) + 
                        groupSumsVar(Delta, k.list[[2L]]) - 
                        crossprod(Delta))
      if (panel.structure == "network") {
        if (k == 3L) {
          V <- V + adj * (groupSumsVar(Delta, k.list[[3L]]) - crossprod(Delta))
        }
      }
    } else {
      V <- V + adj * tcrossprod(colSums(Delta))
    }
    
    # Add covariance in case of weak exogeneity
    if (weak.exo) {
      if (panel.structure == "classic") {
        Cov <- groupSumsCov(Delta, Gamma, k.list[[1L]])
        V <- V + adj * (Cov + t(Cov))
        rm(Cov)
      } else {
        if (k == 3L) {
          Cov <- groupSumsCov(Delta, Gamma, k.list[[3L]])
          V <- V + adj * (Cov + t(Cov))
          rm(Cov)
        }
      }
    }
  }
  V <- V / nt.full^2
  
  # Add names
  names(delta) <- nms.sp
  dimnames(V) <- list(nms.sp, nms.sp)
  
  # Return list
  structure(list(delta           = delta,
                 vcov            = V,
                 panel.structure = panel.structure,
                 sampling.fe     = sampling.fe,
                 weak.exo        = weak.exo), class = "APEs")
}
#' @title
#' Efficiently fit glm's with high-dimensional \eqn{k}-way fixed effects
#' @description
#' \code{feglm} can be used to fit generalized linear models with many high-dimensional fixed
#' effects. The estimation procedure is based on unconditional maximum likelihood and can be
#' interpreted as a \dQuote{pseudo demeaning} approach that combines the work of Gaure (2013) and
#' Stammann et. al. (2016). For technical details see Stammann (2018). The routine is well suited
#' for large data sets that would be otherwise infeasible to use due to memory limitations.
#' 
#' \strong{Remark:} The term fixed effect is used in econometrician's sense of having intercepts for
#' each level in each category.
#' @param
#' formula an object of class \code{"formula"}: a symbolic description of the model to be fitted.
#' \code{formula} must be of type \code{y ~ x | k}, where the second part of the formula refers to
#' factors to be concentrated out.
#' @param
#' data an object of class \code{"data.frame"} containing the variables in the model.
#' @param
#' family a description of the error distribution and link function to be used in the model. 
#' Similiar to \code{glm.fit()} this has to be the result of a call to a family function. Default is
#' \code{binomial()}. See \code{\link{family}} for details of family functions.
#' @param
#' beta.start an optional vector of starting values for the structural parameters in the linear 
#' predictor. Default is \eqn{\boldsymbol{\beta} = \mathbf{0}}{\beta = 0}.
#' @param
#' D.alpha.start an optional vector of starting values for the fixed effects contribution in the
#' linear predictor.
#' @param
#' control a named list of parameters for controlling the fitting process. See 
#' \code{\link{feglm.control}} for details.
#' @details
#' If \code{feglm()} does not converge this is usually a sign of linear dependence between one or 
#' more regressors and a fixed effects category. In this case, you should carefully inspect your
#' model specification.
#' @return
#' The function \code{feglm} returns a named list of class \code{"feglm"}.
#' @references
#' Gaure, S. (2013). "OLS with Multiple High Dimensional Category Variables". Computational
#' Statistics and Data Analysis. 66.
#' @references 
#' Stammann, A., F. Heiss, and D. McFadden (2016). "Estimating Fixed Effects Logit Models with Large
#' Panel Data". Working paper.
#' @references 
#' Stammann, A. (2018). "Fast and Feasible Estimation of Generalized Linear Models with
#' High-Dimensional k-Way Fixed Effects". Working Paper.
#' @examples 
#' \dontrun{
#' # Generate artificial for logit models
#' library(alpaca)
#' data <- simGLM(1000L, 200L, 1805L, model = "logit")
#' 
#' # Fit 'feglm()'
#' mod <- feglm(y ~ x1 + x2 + x3 | i + j, data)
#' summary(mod)
#' 
#' # Recover estimates of fixed effects
#' alpha <- getFEs(mod)
#' }
#' @export
feglm <- function(formula       = NULL,
                  data          = NULL,
                  family        = binomial(),
                  beta.start    = NULL,
                  D.alpha.start = NULL,
                  control       = NULL) {
  # Validity of input argument (formula)
  if (is.null(formula)) {
    stop("'formula' has to be specified.")
  } else if (!inherits(formula, "formula")) {
    stop("'formula' has to be of class formula.")
  }
  
  # Validity of input argument (data)
  if (is.null(data)) {
    stop("'data' has to be specified.")
  } else if (!inherits(data, "data.frame")) {
    stop("'data' has to be of class data.frame.")
  }
  
  # Validity of input argument (family)
  # NOTE: Quasi families not supported since they are no maximum likelihood estimators
  if (!inherits(family, "family")) {
    stop("'family' has to be of class family")
  } else if (family[["family"]] %in% c("quasi", "quasipoisson", "quasibinomial")) {
    stop("Quasi-variants of 'family' are not supported.")
  } else if (family[["family"]] == "gaussian") {
    stop("Linear models are not supported. We recommend using 'lfe' for this purpose!")
  }
  
  # Validity of input argument (control)
  if (is.null(control)) {
    control <- list()
  } else if (!inherits(control, "list")) {
    stop("'control' has to be of class list.")
  }
  
  # Extract control list
  control <- do.call(feglm.control, control)
  
  # Update formula and do further validity check
  formula <- Formula(formula)
  lhs <- attr(formula, "lhs")[[1L]]
  if (length(formula)[[2L]] < 2L || length(formula)[[1L]] > 1L) {
    stop("'formula' uncorrectly specified.")
  }
  
  # Generate model.frame
  data <- as.data.table(data)
  mf <- as.data.table(suppressWarnings(model.frame(formula, data)))
  nobs.full <- nrow(mf)
  nobs.na <- length(attr(mf, "na.action"))
  if (nobs.na > 0L) {
    data <- data[- attr(mf, "na.action")]
  }
  
  # Ensure that model response is in line with the choosen model
  if (family[["family"]] == "binomial") {
    # Transform 'y' and check if the number of levels equals two
    mf[, as.character(lhs) := factor(eval(lhs))]
    if (mf[, length(levels(eval(lhs)))] != 2L) {
      stop("Model response has to be binary.")
    }
    
    # Ensure 'y' is 0-1 encoded
    mf[, as.character(lhs) := as.integer(eval(lhs)) - 1L]
  } else if (family[["family"]] %in% c("Gamma", "inverse.gaussian")) {
    # Check if 'y' is strictly positive
    if (mf[, any(eval(lhs) <= 0.0)]) {
      stop("Model response has to be strictly positive.")
    }
  } else {
    # Check if 'y' is positive
    if (mf[, any(eval(lhs) < 0.0)]) {
      stop("Model response has to be positive.")
    }
  }
  
  # Get name of the fixed effects variables
  k.vars <- attr(terms(formula, rhs = 2L), "term.labels")
  k <- length(k.vars)
  
  # Drop observations that do not contribute to the loglikelihood
  if (family[["family"]] %in% c("binomial", "poisson")) {
    if (control[["drop.pc"]] == TRUE) {
      zzz <- NULL # To suppress "no visible binding for global variable ‘zzz’"
      for (i in k.vars) {
        # Determine which observations too keep
        mf[, zzz := mean(eval(lhs)), by = eval(i)]
        if (family[["family"]] == "binomial") {
          keep <- mf[zzz > 0.0 & zzz < 1.0, which = TRUE]
        } else {
          keep <- mf[zzz > 0.0, which = TRUE]
        }
        
        # Subset model and data frame
        if (length(keep) > 0L) {
          mf <- mf[keep]
          data <- data[keep]
        }
        
        # Remove auxiliary variables
        mf[, zzz := NULL]
        rm(keep)
      }
    } 
  }
  
  # Determine number of dropped observations
  nobs <- nrow(mf)
  nobs.pc <- nobs.full - nobs
  
  # Validity of input argument (D.alpha.start)
  if (!is.null(D.alpha.start)) {
    if (length(D.alpha.start) != nobs) {
      stop("Length of 'D.alpha.start' has to be equal to the number of observations.")
    }
    D.alpha <- D.alpha.start
  } else {
    zzz <- NULL # To suppress "no visible binding for global variable ‘zzz’"
    D.alpha <- numeric(nobs)
    for (i in k.vars) {
      mf[, zzz := mean(eval(lhs)), by = eval(i)]
      if (family[["family"]] == "binomial") {
        D.alpha <- D.alpha + family[["linkfun"]]((mf[["zzz"]] + 0.5) / 2.0) / k
      } else if (family[["family"]] %in% c("Gamma", "inverse.gaussian")) {
        D.alpha <- D.alpha + family[["linkfun"]](mf[["zzz"]]) / k
      } else {
        D.alpha <- D.alpha + family[["linkfun"]](mf[["zzz"]] + 0.1) / k
      }
      mf[, zzz := NULL]
    }
  }
  rm(D.alpha.start)
  
  # Extract model response, regressor matrix, and fixed effects
  y <- mf[[lhs]]
  X <- model.matrix(formula, data, rhs = 1L)[, - 1L, drop = FALSE] # Does not work with 'mf'
  nms.sp <- attr(X, "dimnames")[[2L]] # Saves memory
  attr(X, "dimnames") <- NULL
  fe <- mf[, k.vars, with = FALSE]
  rm(mf)
  
  # Check for linear dependence in 'X'
  p <- ncol(X)
  if (qr(X)[["rank"]] < p) {
    stop("Linear dependent terms detected!")
  }
  
  # Validity of input argument (beta.start)
  if (!is.null(beta.start)) {
    if (length(beta.start) != p) {
      stop("Length of 'beta.start' has to be equal to the number of structural parameters.")
    }
    beta <- beta.start
  } else {
    beta <- numeric(p)
  }
  rm(beta.start)
  
  # Ensure factors are consecutive integers
  fe[, (k.vars) := lapply(.SD, factor)]
  nms.fe <- lapply(fe, levels)
  lvls.k <- sapply(nms.fe, length)
  nms.fe <- unlist(nms.fe)
  fe[, (k.vars) := lapply(.SD, as.integer)]
  
  # Generate auxiliary matrices for the pseudo demeaning
  A <- as.matrix(fe) - 1L
  B <- apply(A, 2L, order) - 1L
  rm(fe)
  
  # Compute inital quantities for the maximization routine
  eta <- as.vector(X %*% beta + D.alpha)
  rm(D.alpha)
  mu <- family[["linkinv"]](eta)
  wt <- rep(1.0, nobs)
  dev <- sum(family[["dev.resids"]](y, mu, wt))
  
  # Start maximization of the log-likelihood
  iter <- 1L
  repeat {
    # Progress information
    if (control[["trace"]] > 0L) {
      cat("--- Iteration=", iter, "---\n")
    }
    
    # Compute IWLS weights and dependent variable
    mu.eta <- family[["mu.eta"]](eta)
    w.tilde <- sqrt(mu.eta^2 / family[["variance"]](mu))
    nu <- as.vector((y - mu) / mu.eta)
    
    # Pseudo-demean variables
    M <- pseudo.demeaning(cbind(nu, X) * w.tilde, w.tilde, A, B, lvls.k, control[["pseudo.tol"]])
    
    # Compute update step using gradient and Hessian
    g <- t(M[, seq(2L, p + 1L)]) %*% M[, 1L]
    H <- crossprod(M[, seq(2L, p + 1L)])
    beta.upd <- as.vector(ginv(H) %*% g)
    # Or qr.solve ?
    
    # Update \eta
    eta.upd <- as.vector(nu - (M[, 1L] - M[, seq(2L, p + 1L), drop = FALSE] %*% beta.upd) / w.tilde)
    
    # Step correction based on residual deviance
    dev.old <- dev
    rho <- 1.0
    repeat {
      # Compute residual deviance
      mu <- family[["linkinv"]](eta + rho * eta.upd)
      dev <- sum(family[["dev.resids"]](y, mu, wt))
      
      # Check if deviance is not increasing
      if (is.finite(dev) && dev <= dev.old) {
        # Update \eta and \beta
        eta <- eta + rho * eta.upd
        beta <- beta + rho * beta.upd
        
        # Progress information
        if (control[["trace"]] > 1L) {
          cat("dev=", dev, "\n")
          cat("beta=", round(beta, 4L), "\n")
        }
        
        # Leave step-correction
        break
      }
      
      # Update \rho
      rho <- rho / 2.0
      
      # If \rho is to small throw error
      if (rho < control[["rho.tol"]]) {
        stop("Backtracking (step-halving) failed.")
      }
    }
    
    # Check termination condition
    tc1 <- (dev.old - dev) / dev < control[["dev.tol"]]
    tc2 <- sqrt(as.double(beta.upd %*% beta.upd)) < control[["step.tol"]]
    tc3 <- iter == control[["iter.max"]]
    if (tc1 || tc2 || tc3) {
      # Determine termination code
      if (tc1) {
        conv <- 0L
        if (control[["trace"]] > 0L) {
          cat("Convergence (small change in deviance)\n")
        }
      } else if (tc2) {
        conv <- 1L
        if (control[["trace"]] > 0L) {
          cat("Convergence (small step size)\n")
        }
      } else {
        if (control[["trace"]] > 0L) {
          cat("No convergence\n")
        }
        conv <- 9L
      }
      
      # Final computations and leave maximization routine
      D.alpha <- eta - X %*% beta
      G <- M[, seq(2L, p + 1L), drop = FALSE] * M[, 1L]
      break
    }
    
    # Increase counter
    iter <- iter + 1L
  }
  rm(eta, w.tilde, y, X, A, B)
  
  # Generate result list
  res <- list(coefficients  = beta,
              D.alpha       = D.alpha,
              Score         = G,
              Hessian       = H,
              deviance      = dev,
              conv          = conv,
              iter          = iter,
              nobs          = nobs,
              nobs.na       = nobs.na,
              nobs.pc       = nobs.pc,
              lvls.k        = lvls.k,
              nms.fe        = nms.fe,
              formula       = formula,
              data          = data,
              family        = family,
              control       = control)
  
  # Modify result list
  names(res[["coefficients"]]) <- nms.sp
  colnames(res[["Score"]]) <- nms.sp
  dimnames(res[["Hessian"]]) <- list(nms.sp, nms.sp)
  
  # Return result list
  structure(res, class = "feglm")
}
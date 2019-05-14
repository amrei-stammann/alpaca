#' @title
#' Efficiently fit glm's with high-dimensional \eqn{k}-way fixed effects
#' @description
#' \code{\link{feglm}} can be used to fit generalized linear models with many high-dimensional fixed
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
#' factors to be concentrated out. It is also possible to pass additional variables to 
#' \code{\link{feglm}} (e.g. to cluster standard errors). This can be done by specifying the third 
#' part of the formula: \code{y ~ x | k | add}.
#' @param
#' data an object of class \code{"data.frame"} containing the variables in the model.
#' @param
#' family a description of the error distribution and link function to be used in the model. 
#' Similiar to \code{\link[stats]{glm.fit}} this has to be the result of a call to a family 
#' function. Default is \code{binomial()}. See \code{\link[stats]{family}} for details of family 
#' functions.
#' @param
#' beta.start an optional vector of starting values for the structural parameters in the linear 
#' predictor. Default is \eqn{\boldsymbol{\beta} = \mathbf{0}}{\beta = 0}.
#' @param
#' eta.start an optional vector of starting values for the linear predictor.
#' @param
#' control a named list of parameters for controlling the fitting process. See 
#' \code{\link{feglmControl}} for details.
#' @details
#' If \code{\link{feglm}} does not converge this is usually a sign of linear dependence between 
#' one or more regressors and a fixed effects category. In this case, you should carefully inspect 
#' your model specification.
#' @return
#' The function \code{\link{feglm}} returns a named list of class \code{"feglm"}.
#' @references
#' Gaure, S. (2013). "OLS with Multiple High Dimensional Category Variables". Computational
#' Statistics and Data Analysis, 66.
#' @references 
#' Stammann, A., F. Heiss, and D. McFadden (2016). "Estimating Fixed Effects Logit Models with Large
#' Panel Data". Working paper.
#' @references 
#' Stammann, A. (2018). "Fast and Feasible Estimation of Generalized Linear Models with
#' High-Dimensional k-Way Fixed Effects". ArXiv e-prints.
#' @examples 
#' \donttest{
#' # Generate an artificial data set for logit models
#' library(alpaca)
#' data <- simGLM(1000L, 20L, 1805L, model = "logit")
#' 
#' # Fit 'feglm()'
#' mod <- feglm(y ~ x1 + x2 + x3 | i + t, data)
#' summary(mod)
#' }
#' @export
feglm <- function(formula       = NULL,
                  data          = NULL,
                  family        = binomial(),
                  beta.start    = NULL,
                  eta.start     = NULL,
                  control       = NULL) {
  # Validity of input argument (formula)
  if (is.null(formula)) {
    stop("'formula' has to be specified.", call. = FALSE)
  } else if (!inherits(formula, "formula")) {
    stop("'formula' has to be of class formula.", call. = FALSE)
  }
  
  # Validity of input argument (data)
  if (is.null(data)) {
    stop("'data' has to be specified.", call. = FALSE)
  } else if (!inherits(data, "data.frame")) {
    stop("'data' has to be of class data.frame.", call. = FALSE)
  }
  
  # Validity of input argument (family)
  # NOTE: Quasi families not supported since they are no maximum likelihood estimators
  if (!inherits(family, "family")) {
    stop("'family' has to be of class family", call. = FALSE)
  } else if (family[["family"]] %in% c("quasi", "quasipoisson", "quasibinomial")) {
    stop("Quasi-variants of 'family' are not supported.", call. = FALSE)
  } else if (family[["family"]] == "gaussian") {
    stop("Linear models are not supported. We recommend using 'lfe' for this purpose.",
         call. = FALSE)
  } else if (startsWith(family[["family"]], "Negative Binomial")) {
    stop("Please use 'feglm.nb' instead.", call. = FALSE)
  }
  
  # Validity of input argument (control)
  if (is.null(control)) {
    control <- list()
  } else if (!inherits(control, "list")) {
    stop("'control' has to be of class list.", call. = FALSE)
  }
  
  # Extract control list
  control <- do.call(feglmControl, control)
  
  # Update formula and do further validity check
  formula <- Formula(formula)
  lhs <- attr(formula, "lhs")[[1L]]
  if (length(formula)[[2L]] < 2L || length(formula)[[1L]] > 1L) {
    stop("'formula' uncorrectly specified.", call. = FALSE)
  }
  
  # Generate model.frame
  if (is.null(attr(data, "terms"))) {
    data <- suppressWarnings(model.frame(formula, data))
  }
  data <- as.data.table(data)
  nobs.full <- nrow(data)
  nobs.na <- length(attr(data, "na.action"))
  
  # Ensure that model response is in line with the choosen model
  if (family[["family"]] == "binomial") {
    # Check if 'y' is numeric
    if (data[, is.numeric(eval(lhs))]) {
      # Check if 'y' is in [0, 1]
      if (data[, any(eval(lhs) < 0.0 | eval(lhs) > 1.0)]) {
        stop("Model response has to be within the unit interval.", call. = FALSE)
      }
    } else {
      # Check if 'y' is factor and transform otherwise
      data[, (1L) := checkFactor(eval(lhs))]
      
      # Check if the number of levels equals two
      if (data[, length(levels(eval(lhs)))] != 2L) {
        stop("Model response has to be binary.", call. = FALSE)
      }
      
      # Ensure 'y' is 0-1 encoded
      data[, (1L) := as.integer(eval(lhs)) - 1L]
    }
  } else if (family[["family"]] %in% c("Gamma", "inverse.gaussian")) {
    # Check if 'y' is strictly positive
    if (data[, any(eval(lhs) <= 0.0)]) {
      stop("Model response has to be strictly positive.", call. = FALSE)
    }
  } else {
    # Check if 'y' is positive
    if (data[, any(eval(lhs) < 0.0)]) {
      stop("Model response has to be positive.", call. = FALSE)
    }
  }
  
  # Get names of the fixed effects variables
  k.vars <- attr(terms(formula, rhs = 2L), "term.labels")
  k <- length(k.vars)
  
  # Drop observations that do not contribute to the loglikelihood
  if (family[["family"]] %in% c("binomial", "poisson")) {
    if (control[["drop.pc"]]) {
      trms <- attr(data, "terms") # Store terms; required for model matrix
      tmp.var <- tempVar(data)
      for (i in k.vars) {
        data[, (tmp.var) := mean(eval(lhs)), by = eval(i)]
        if (family[["family"]] == "binomial") {
          data <- data[get(tmp.var) > 0.0 & get(tmp.var) < 1.0]
        } else {
          data <- data[get(tmp.var) > 0.0]
        }
        data[, (tmp.var) := NULL]
      }
      setattr(data, "terms", trms)
      rm(trms)
    }
  }
  
  # Transform fixed effects variables and potential cluster variables to factors
  data[, (k.vars) := lapply(.SD, checkFactor), .SDcols = k.vars]
  if (length(formula)[[2L]] > 2L) {
    add.vars <- attr(terms(formula, rhs = 3L), "term.labels")
    data[, (add.vars) := lapply(.SD, checkFactor), .SDcols = add.vars]
  }
  
  # Determine number of dropped observations
  nobs <- c(nobs.full = nobs.full,
            nobs.na   = nobs.na,
            nobs.pc   = nobs.full - nrow(data),
            nobs      = nrow(data))
  
  # Extract model response and regressor matrix
  y <- data[[1L]]
  X <- model.matrix(formula, data, rhs = 1L)[, - 1L, drop = FALSE]
  nms.sp <- attr(X, "dimnames")[[2L]] # Saves memory
  attr(X, "dimnames") <- NULL
  
  # Check for linear dependence in 'X'
  p <- ncol(X)
  if (qr(X)[["rank"]] < p) {
    stop("Linear dependent terms detected!", call. = FALSE)
  }
  
  # Compute and check starting guesses
  if (!is.null(beta.start) || !is.null(eta.start)) {
    # If both are specified, ignore eta.start
    if (!is.null(beta.start) && !is.null(eta.start)) {
      warning("'beta.start' and 'eta.start' are specified. Ignoring 'eta.start'.", call. = FALSE)
    }
    
    # Compute and check starting guesses
    if (!is.null(beta.start)) {
      # Validity of input argument (beta.start)
      if (length(beta.start) != p) {
        stop("Length of 'beta.start' has to be equal to the number of structural parameters.",
             call. = FALSE)
      }
      
      # Set starting guesses
      beta <- beta.start
      eta <- as.vector(X %*% beta)
    } else {
      # Validity of input argument (eta.start)
      if (length(eta.start) != nobs[[4L]]) {
        stop("Length of 'eta.start' has to be equal to the number of observations.", call. = FALSE)
      }
      
      # Set starting guesses
      beta <- numeric(p)
      eta <- eta.start
    }
  } else {
    # Compute starting guesses if not user specified
    beta <- numeric(p)
    tmp.var <- tempVar(data)
    eta <- numeric(nobs[[4L]])
    for (i in k.vars) {
      data[, (tmp.var) := mean(eval(lhs)), by = eval(i)]
      if (family[["family"]] == "binomial") {
        eta <- eta + family[["linkfun"]]((data[[tmp.var]] + 0.5) / 2.0) / k
      } else if (family[["family"]] %in% c("Gamma", "inverse.gaussian")) {
        eta <- eta + family[["linkfun"]](data[[tmp.var]]) / k
      } else {
        eta <- eta + family[["linkfun"]](data[[tmp.var]] + 0.1) / k
      }
      data[, (tmp.var) := NULL]
    }
  }
  rm(beta.start, eta.start)
  
  # Ensure factors are consecutive integers and generate auxiliary matrices to center variables
  fe <- model.part(formula, data, rhs = 2L)
  nms.fe <- lapply(fe, levels)
  fe[, (k.vars) := lapply(.SD, as.integer)]
  lvls.k <- sapply(fe, max)
  A <- as.matrix(fe) - 1L
  rm(fe)
  dimnames(A) <- NULL
  B <- apply(A, 2L, order) - 1L
  
  # Fit generalized linear model
  fit <- feglmFit(beta, eta, y, X, A, B, lvls.k, family, control)
  rm(y, X, eta, A, B)
  
  # Add names to \beta, Scores, and Hessian
  names(fit[["coefficients"]]) <- nms.sp
  colnames(fit[["Score"]]) <- nms.sp
  dimnames(fit[["Hessian"]]) <- list(nms.sp, nms.sp)
  
  # Return result list
  structure(c(fit, list(nobs    = nobs,
                        lvls.k  = lvls.k,
                        nms.fe  = nms.fe,
                        formula = formula,
                        data    = data,
                        family  = family,
                        control = control)), class = "feglm")
}


#' @title
#' Set \code{feglm} Control Parameters
#' @description
#' Set and change parameters used for fitting \code{\link{feglm}}.
#' 
#' \strong{Note:} \code{\link{feglm.control}} is deprecated and will be removed soon.
#' @param
#' dev.tol tolerance level for the first stopping condition of the maximization routine. The 
#' stopping condition is based on the relative change of the deviance in iteration \eqn{r}
#' and can be expressed as follows: \eqn{(dev_{r - 1} - dev_{r}) / (0.1 + dev_{r}) < 
#' tol}{\Delta dev / (0.1 + dev) < tol}. Default is \code{1.0e-08}.
#' @param
#' center.tol tolerance level for the stopping condition of the centering algorithm.
#' The stopping condition is based on the relative change of euclidean norm in iteration \eqn{i} and
#' can be expressed as follows: \eqn{||\mathbf{v}_{i} - \mathbf{v}_{i - 1}||_{2} < 
#' tol ||\mathbf{v}_{i - 1}||}{||\Delta v|| / ||v_old|| < tol}. Default is
#' \code{1.0e-05}.
#' @param
#' rho.tol tolerance level for the stephalving in the maximization routine. Stephalving only takes
#' place if the deviance in iteration \eqn{r} is larger than the one of the previous iteration. If 
#' this is the case, 
#' \eqn{||\boldsymbol{\beta}_{r} - \boldsymbol{\beta}_{r - 1}||_{2}}{||\Delta \beta||} is 
#' halfed until the deviance is less or numerically equal compared to the deviance of the previous
#' iteration. Stephalving fails if the the following condition holds: \eqn{\rho < tol}{\rho < tol}, 
#' where \eqn{\rho}{\rho} is the stepcorrection factor. If stephalving fails the maximization
#' routine is canceled. Default is \code{1.0e-04}.
#' @param
#' conv.tol tolerance level that accounts for rounding errors inside the stephalving routine when
#' comparing the deviance with the one of the previous iteration. Default is \code{1.0e-06}.
#' @param
#' iter.max unsigned integer indicating the maximum number of iterations in the maximization
#' routine. Default is \code{100L}.
#' @param
#' limit unsigned integer indicating the maximum number of iterations of 
#' \code{\link[MASS]{theta.ml}}. Default is \code{10L}.
#' @param
#' trace logical indicating if output should be produced in each iteration. Default is \code{FALSE}.
#' @param
#' drop.pc logical indicating to drop observations that are perfectly classified 
#' (perfectly seperated) and hence do not contribute to the log-likelihood. This option is useful to
#' reduce the computational costs of the maximization problem, since it reduces the number of
#' observations and does not affect the estimates. Default is \code{TRUE}.
#' @param
#' pseudo.tol deprecated; use \code{center.tol} instead.
#' @param
#' step.tol depreacted; termination conditions is now similar to \code{\link[stats]{glm}}.
#' @param
#' ... arguments passed to the deprecated function \code{\link{feglmControl}}.
#' @return
#' The function \code{\link{feglmControl}} returns a named list of control 
#' parameters.
#' @seealso
#' \code{\link{feglm}}
#' @export
feglmControl <- function(dev.tol    = 1.0e-08,
                         center.tol = 1.0e-05,
                         rho.tol    = 1.0e-04,
                         conv.tol   = 1.0e-06,
                         iter.max   = 100L,
                         limit      = 10L,
                         trace      = FALSE,
                         drop.pc    = TRUE,
                         pseudo.tol = NULL,
                         step.tol   = NULL) {
  # 'pseudo.tol' is deprecated
  if (!is.null(pseudo.tol)) {
    warning("'pseudo.tol' is deprecated; please use 'center.tol' instead.", call. = FALSE)
    center.tol <- pseudo.tol
  }
  
  # 'step.tol' is deprecated
  if (!is.null(step.tol)) {
    warning("'step.tol' is deprecated;", call. = FALSE)
  }
  
  # Check validity of tolerance parameters
  if (dev.tol <= 0.0 || center.tol <= 0.0 || rho.tol <= 0.0 || conv.tol <= 0.0) {
    stop("All tolerance paramerters should be greater than zero.", call. = FALSE)
  }
  
  # Check validity of 'iter.max'
  iter.max <- as.integer(iter.max)
  if (iter.max < 1L) {
    stop("Maximum number of iterations should be at least one.", call. = FALSE)
  }
  
  # Check validity of 'limit'
  limit <- as.integer(limit)
  if (limit < 1L) {
    stop("Maximum number of iterations should be at least one.", call. = FALSE)
  }
  
  # Return list with control parameters
  list(dev.tol    = dev.tol,
       center.tol = center.tol,
       rho.tol    = rho.tol,
       conv.tol   = conv.tol,
       iter.max   = iter.max,
       limit      = limit,
       trace      = as.logical(trace),
       drop.pc    = as.logical(drop.pc))
}


### Internal function (not exported)


# Fitting algorithm (similar to glm.fit)
feglmFit <- function(beta, eta, y, X, A, B, lvls.k, family, control) {
  # Compute initial quantities for the maximization routine
  mu <- family[["linkinv"]](eta)
  wt <- rep(1.0, length(y))
  dev <- sum(family[["dev.resids"]](y, mu, wt))
  null.dev <- sum(family[["dev.resids"]](y, mean(y), wt))
  
  # Start maximization of the log-likelihood
  for (iter in seq.int(control[["iter.max"]])) {
    # Compute weights and dependent variable
    mu.eta <- family[["mu.eta"]](eta)
    w.tilde <- sqrt(mu.eta^2 / family[["variance"]](mu))
    nu <- as.vector((y - mu) / mu.eta)
    
    # Centering variables
    M <- centerVariables(cbind(nu, X) * w.tilde, w.tilde, A, B, lvls.k, control[["center.tol"]])
    
    # Compute update step and update \eta
    beta.upd <- try(solveNR(M[, - 1L, drop = FALSE], M[, 1L]), silent = TRUE)
    if (inherits(beta.upd, "try-error")) {
      stop("Failure to solve Newton equations.", call. = FALSE)
    }
    eta.upd <- as.vector(nu - (M[, 1L] - M[, - 1L, drop = FALSE] %*% beta.upd) / w.tilde)
    
    # Stephalving based on residual deviance as common in glm's
    dev.old <- dev
    rho <- 1.0
    repeat {
      # Compute residual deviance
      mu <- family[["linkinv"]](eta + rho * eta.upd)
      dev <- sum(family[["dev.resids"]](y, mu, wt))
      
      # Check if deviance is not increasing
      if (is.finite(dev) && dev <= dev.old + control[["conv.tol"]] * dev) {
        # Update \eta, \beta, and leave stephalving
        eta <- eta + rho * eta.upd
        beta <- beta + rho * beta.upd
        break
      }
      
      # Update \rho
      rho <- rho / 2.0
      
      # If \rho is to small throw error
      if (rho < control[["rho.tol"]]) {
        # TODO: Add info that this might be due perfect colliniearity with fixed effects
        stop("Failure in stephalving.", call. = FALSE)
      }
    }
    
    # Progress information
    if (control[["trace"]]) {
      cat("Deviance=", round(dev, 2L), "- Iteration= ", iter, "\n")
    }
    
    # Check termination condition
    if ((dev.old - dev) / (0.1 + dev) < control[["dev.tol"]]) {
      if (control[["trace"]]) {
        cat("Convergence\n")
      }
      conv <- TRUE
      break
    }
  }
  
  # Compute Score and Hessian
  G <- M[, - 1L, drop = FALSE] * M[, 1L]
  H <- crossprod(M[, - 1L])
  
  # Return result list
  list(coefficients  = beta,
       eta           = eta,
       Score         = G,
       Hessian       = H,
       deviance      = dev,
       null.deviance = null.dev,
       conv          = conv,
       iter          = iter)
}


# Efficient offset algorithm to update the linear predictor
feglmOffset <- function(object, offset) {
  # Check validity of 'object'
  if (!inherits(object, "feglm")) {
    stop("'feglmOffset' called on a non-'feglm' object.")
  }
  
  # Extract required quantities from result list
  formula <- object[["formula"]]
  family <- object[["family"]]
  control <- object[["control"]]
  lvls.k <- object[["lvls.k"]]
  nobs <- object[["nobs"]][["nobs"]]
  lhs <- attr(formula, "lhs")[[1L]]
  k.vars <- names(lvls.k)
  k <- length(lvls.k)
  
  # Compute starting guess for \eta
  tmp.var <- tempVar(object[["data"]])
  eta <- numeric(nobs)
  for (i in k.vars) {
    object[["data"]][, (tmp.var) := mean(eval(lhs)), by = eval(i)]
    if (family[["family"]] == "binomial") {
      eta <- eta + family[["linkfun"]]((object[["data"]][[tmp.var]] + 0.5) / 2.0) / k
    } else if (family[["family"]] %in% c("Gamma", "inverse.gaussian")) {
      eta <- eta + family[["linkfun"]](object[["data"]][[tmp.var]]) / k
    } else {
      eta <- eta + family[["linkfun"]](object[["data"]][[tmp.var]] + 0.1) / k
    }
    object[["data"]][, (tmp.var) := NULL]
  }
  
  # Construct auxiliary matrix to flatten the fixed effects
  fe <- model.part(formula, object[["data"]], rhs = 2L)
  fe[, (k.vars) := lapply(.SD, as.integer)]
  A <- as.matrix(fe) - 1L
  dimnames(A) <- NULL
  rm(fe)
  B <- apply(A, 2L, order) - 1L

  # Extract dependent variable and compute initial quantities for the maximization routine
  y <- object[["data"]][[1L]]
  mu <- family[["linkinv"]](eta)
  wt <- rep(1.0, nobs)
  dev <- sum(family[["dev.resids"]](y, mu, wt))
  
  # Start maximization of the log-likelihood
  for (iter in seq.int(control[["iter.max"]])) {
    # Compute weights and dependent variable
    mu.eta <- family[["mu.eta"]](eta)
    w.tilde <- sqrt(mu.eta^2 / family[["variance"]](mu))
    q.tilde <- as.matrix(eta - offset + (y - mu) / mu.eta) * w.tilde
    
    # Centering dependent variable and compute \eta update
    Mq <- centerVariables(q.tilde, w.tilde, A, B, lvls.k, control[["center.tol"]])
    eta.upd <- as.vector((q.tilde - Mq)) / w.tilde + offset - eta
    
    # Stephalving based on residual deviance as common in glm's
    dev.old <- dev
    rho <- 1.0
    repeat {
      # Compute residual deviance
      mu <- family[["linkinv"]](eta + rho * eta.upd)
      dev <- sum(family[["dev.resids"]](y, mu, wt))
      
      # Check if deviance is not increasing
      if (is.finite(dev) && dev <= dev.old + control[["conv.tol"]] * dev) {
        # Update \eta and leave stephalving
        eta <- eta + rho * eta.upd
        break
      }
      
      # Update \rho
      rho <- rho / 2.0
      
      # If \rho is to small throw error
      if (rho < control[["rho.tol"]]) {
        # TODO: Add info that this might be due perfect colliniearity with fixed effects
        stop("Failure in stephalving.", call. = FALSE)
      }
    }
    
    # Check termination condition
    if ((dev.old - dev) / (0.1 + dev) < control[["dev.tol"]]) {
      break
    }
  }
  
  # Return \eta
  eta
}


### Deprecated functions

#' @rdname feglmControl
#' @aliases feglmControl
#' @export
feglm.control <- function(...) {
  .Deprecated("feglmControl")
  do.call(feglmControl, list(...))
}

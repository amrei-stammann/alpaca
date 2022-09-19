#' @title
#' Efficiently fit negative binomial glm's with high-dimensional \eqn{k}-way fixed effects
#' @description
#' \code{feglm.nb} can be used to fit negative binomial generalized linear models with many 
#' high-dimensional fixed effects (see \code{\link{feglm}}).
#' @param
#' formula,data,weights,beta.start,eta.start,control see \code{\link{feglm}}.
#' @param
#' init.theta an optional initial value for the theta parameter (see \code{\link[MASS]{glm.nb}}).
#' @param
#' link the link function. Must be one of \code{"log"}, \code{"sqrt"}, or \code{"identity"}.
#' @details
#' If \code{feglm.nb} does not converge this is usually a sign of linear dependence between one or 
#' more regressors and a fixed effects category. In this case, you should carefully inspect your
#' model specification.
#' @return
#' The function \code{feglm.nb} returns a named list of class \code{"feglm"}.
#' @references
#' Gaure, S. (2013). "OLS with Multiple High Dimensional Category Variables". Computational
#' Statistics and Data Analysis. 66.
#' @references 
#' Marschner, I. (2011). "glm2: Fitting generalized linear models with convergence problems". 
#' The R Journal, 3(2).
#' @references 
#' Stammann, A., F. Heiss, and D. McFadden (2016). "Estimating Fixed Effects Logit Models with Large
#' Panel Data". Working paper.
#' @references 
#' Stammann, A. (2018). "Fast and Feasible Estimation of Generalized Linear Models with
#' High-Dimensional k-Way Fixed Effects". ArXiv e-prints.
#' @seealso
#' \code{\link[MASS]{glm.nb}}, \code{\link{feglm}} 
#' @export
feglm.nb <- function(
  formula       = NULL,
  data          = NULL,
  weights       = NULL,
  beta.start    = NULL,
  eta.start     = NULL,
  init.theta    = NULL,
  link          = c("log", "identity", "sqrt"),
  control       = NULL
  ) {
  # Check validity of formula
  if (is.null(formula)) {
    stop("'formula' has to be specified.", call. = FALSE)
  } else if (!inherits(formula, "formula")) {
    stop("'formula' has to be of class formula.", call. = FALSE)
  }
  
  # Check validity of data
  if (is.null(data)) {
    stop("'data' has to be specified.", call. = FALSE)
  } else if (!inherits(data, "data.frame")) {
    stop("'data' has to be of class data.frame.", call. = FALSE)
  }
  
  # Check validity of link
  link <- match.arg(link)
  
  # Check validity of control
  if (is.null(control)) {
    control <- list()
  } else if (!inherits(control, "list")) {
    stop("'control' has to be a list.", call. = FALSE)
  }
  
  # Extract control list
  control <- do.call(feglmControl, control)
  
  # Update formula and do further validity check
  formula <- Formula(formula)
  if (length(formula)[[2L]] < 2L || length(formula)[[1L]] > 1L) {
    stop("'formula' uncorrectly specified.", call. = FALSE)
  }
  
  # Generate model.frame
  setDT(data)
  data <- data[, c(all.vars(formula), weights), with = FALSE]
  lhs <- names(data)[[1L]]
  nobs.full <- nrow(data)
  data <- na.omit(data)
  nobs.na <- nobs.full - nrow(data)
  nobs.full <- nrow(data)
  
  # Ensure that model response is in line with the chosen model
  if (data[, any(get(lhs) < 0.0)]) {
    stop("Model response has to be positive.", call. = FALSE)
  }
  
  # Get names of the fixed effects variables and sort
  k.vars <- attr(terms(formula, rhs = 2L), "term.labels")
  k <- length(k.vars)
  setkeyv(data, k.vars)
  
  # Generate temporary variable
  tmp.var <- tempVar(data)
  
  # Drop observations that do not contribute to the loglikelihood
  if (control[["drop.pc"]]) {
    for (j in k.vars) {
      data[, (tmp.var) := mean(get(lhs)), by = eval(j)]
      data <- data[get(tmp.var) > 0.0]
      data[, (tmp.var) := NULL]
    }
  }
  
  # Transform fixed effects variables and potential cluster variables to factors
  data[, (k.vars) := lapply(.SD, checkFactor), .SDcols = k.vars]
  if (length(formula)[[2L]] > 2L) {
    add.vars <- attr(terms(formula, rhs = 3L), "term.labels")
    data[, (add.vars) := lapply(.SD, checkFactor), .SDcols = add.vars]
  }
  
  # Determine number of dropped observations
  nt <- nrow(data)
  nobs <- c(
    nobs.full = nobs.full,
    nobs.na   = nobs.na,
    nobs.pc   = nobs.full - nt,
    nobs      = nt
  )
  
  # Extract model response and regressor matrix
  y <- data[[1L]]
  X <- model.matrix(formula, data, rhs = 1L)[, - 1L, drop = FALSE]
  nms.sp <- attr(X, "dimnames")[[2L]]
  attr(X, "dimnames") <- NULL
  p <- ncol(X)
  
  # Check for linear dependence in 'X'
  if (qr(X)[["rank"]] < p) {
    stop("Linear dependent terms detected.", call. = FALSE)
  }
  
  # Extract weights if required
  if (is.null(weights)) {
    wt <- rep(1.0, nt)
  } else {
    wt <- data[[weights]]
  }
  
  # Check validity of weights
  if (!is.numeric(wt)) {
    stop("weights must be numeric.", call. = FALSE)
  }
  if (any(wt < 0.0)) {
    stop("negative weights are not allowed.", call. = FALSE)
  }
  
  # Check starting guess of \theta
  if (is.null(init.theta)) {
    family <- poisson(link)
  } else {
    # Validity of input argument (beta.start)
    if (length(init.theta) != 1L) {
      stop("'init.theta' has to be a scalar.", call. = FALSE)
    } else if (init.theta <= 0.0) {
      stop("'init.theta' has to be strictly positive.", call. = FALSE)
    }
    family <- negative.binomial(init.theta, link)
  }
  rm(init.theta)
  
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
      if (length(eta.start) != nt) {
        stop("Length of 'eta.start' has to be equal to the number of observations.", call. = FALSE)
      }
      
      # Set starting guesses
      beta <- numeric(p)
      eta <- eta.start
    }
  } else {
    # Compute starting guesses if not user specified
    beta <- numeric(p)
    eta <- rep(family[["linkfun"]](sum(wt * (y + 0.1)) / sum(wt)), nt)
  }
  rm(beta.start, eta.start)
  
  # Get names and number of levels in each fixed effects category
  nms.fe <- lapply(data[, k.vars, with = FALSE], levels)
  lvls.k <- sapply(nms.fe, length)
  
  # Generate auxiliary list of indexes for different sub panels
  k.list <- getIndexList(k.vars, data)
  
  # Extract control arguments
  tol <- control[["dev.tol"]]
  limit <- control[["limit"]]
  iter.max <- control[["iter.max"]]
  trace <- control[["trace"]]
  
  # Initial negative binomial fit
  fit <- feglmFit(beta, eta, y, X, wt, k.list, family, control)
  beta <- fit[["coefficients"]]
  eta <- fit[["eta"]]
  dev <- fit[["deviance"]]
  theta <- suppressWarnings(
    theta.ml(
      y     = y,
      mu    = family[["linkinv"]](eta),
      n     = nt,
      limit = limit,
      trace = trace
    )
  )
  
  # Alternate between fitting glm and \theta
  conv <- FALSE
  for (iter in seq.int(iter.max)) {
    # Fit negative binomial model
    dev.old <- dev
    theta.old <- theta
    family <- negative.binomial(theta, link)
    fit <- feglmFit(beta, eta, y, X, wt, k.list, family, control)
    beta <- fit[["coefficients"]]
    eta <- fit[["eta"]]
    dev <- fit[["deviance"]]
    theta <- suppressWarnings(
      theta.ml(
        y     = y,
        mu    = family[["linkinv"]](eta),
        n     = nt,
        limit = limit,
        trace = trace
        )
      )
    
    # Progress information
    if (trace) {
      cat("Outer Iteration=", iter, "\n")
      cat("Deviance=", format(dev, digits = 5L, nsmall = 2L), "\n")
      cat("theta=", format(theta, digits = 5L, nsmall = 2L), "\n")
      cat("Estimates=", format(beta, digits = 3L, nsmall = 2L), "\n")
    }
    
    # Check termination condition
    dev.crit <- abs(dev - dev.old) / (0.1 + abs(dev))
    theta.crit <- abs(theta - theta.old) / (0.1 + abs(theta.old))
    if (dev.crit <= tol && theta.crit <= tol) {
      if (trace) {
        cat("Convergence\n")
      }
      conv <- TRUE
      break
    }
  }
  rm(y, X, eta)
  
  # Information if convergence failed
  if (!conv && trace) cat("Algorithm did not converge.\n") 
  
  # Add names to \beta, Hessian, and MX (if provided)
  names(fit[["coefficients"]]) <- nms.sp
  if (control[["keep.mx"]]) {
    colnames(fit[["MX"]]) <- nms.sp
  }
  dimnames(fit[["Hessian"]]) <- list(nms.sp, nms.sp)
  
  # Return result list
  structure(c(fit, list(
    theta      = theta,
    iter.outer = iter,
    conv.outer = conv,
    nobs       = nobs,
    lvls.k     = lvls.k,
    nms.fe     = nms.fe,
    formula    = formula,
    data       = data,
    family     = family,
    control    = control)
    ), class = c("feglm", "feglm.nb"))
}
#' @title
#' Efficiently fit negative binomial glm's with high-dimensional \eqn{k}-way fixed effects
#' @description
#' \code{feglm.nb} can be used to fit negative binomial generalized linear models with many 
#' high-dimensional fixed effects (see \code{\link{feglm}}).
#' @param
#' formula see \code{\link{feglm}}.
#' @param
#' data see \code{\link{feglm}}.
#' @param
#' beta.start see \code{\link{feglm}}.
#' @param
#' eta.start see \code{\link{feglm}}.
#' @param
#' init.theta an optional initial value for the theta parameter (see \code{\link[MASS]{glm.nb}}).
#' @param
#' link the link function. Must be one of \code{"log"}, \code{"sqrt"}, or \code{"identity"}.
#' @param
#' control see \code{\link{feglm}}.
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
#' Stammann, A., F. Heiss, and D. McFadden (2016). "Estimating Fixed Effects Logit Models with Large
#' Panel Data". Working paper.
#' @references 
#' Stammann, A. (2018). "Fast and Feasible Estimation of Generalized Linear Models with
#' High-Dimensional k-Way Fixed Effects". ArXiv e-prints.
#' @seealso
#' \code{\link[MASS]{glm.nb}}, \code{\link{feglm}} 
#' @export
feglm.nb <- function(formula       = NULL,
                     data          = NULL,
                     beta.start    = NULL,
                     eta.start     = NULL,
                     init.theta    = NULL,
                     link          = c("log", "identity", "sqrt"),
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
  
  # Validity of input argument (link)
  link <- match.arg(link)
  
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
  if (data[, any(eval(lhs) < 0.0)]) {
    stop("Model response has to be positive.", call. = FALSE)
  }
  
  # Get names of the fixed effects variables
  k.vars <- attr(terms(formula, rhs = 2L), "term.labels")
  k <- length(k.vars)
  
  # Drop observations that do not contribute to the loglikelihood
  if (control[["drop.pc"]]) {
    trms <- attr(data, "terms") # Store terms; required for model matrix
    tmp.var <- tempVar(data)
    for (i in k.vars) {
      data[, (tmp.var) := mean(eval(lhs)), by = eval(i)]
      data <- data[get(tmp.var) > 0.0]
      data[, (tmp.var) := NULL]
    }
    setattr(data, "terms", trms)
    rm(trms)
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
      eta <- eta + family[["linkfun"]](data[[tmp.var]] + 0.1) / k
      data[, (tmp.var) := NULL]
    }
  }
  rm(beta.start, eta.start, init.theta)
  
  # Ensure factors are consecutive integers and generate auxiliary matrices to center variables
  fe <- model.part(formula, data, rhs = 2L)
  nms.fe <- lapply(fe, levels)
  fe[, (k.vars) := lapply(.SD, as.integer)]
  lvls.k <- sapply(fe, max)
  A <- as.matrix(fe) - 1L
  rm(fe)
  dimnames(A) <- NULL
  B <- apply(A, 2L, order) - 1L
  
  # Initial negative binomial fit
  fit <- feglmFit(beta, eta, y, X, A, B, lvls.k, family, control)
  beta <- fit[["coefficients"]]
  eta <- fit[["eta"]]
  dev <- fit[["deviance"]]
  theta <- suppressWarnings(theta.ml(y, family[["linkinv"]](eta), n = nobs[[4L]],
                                     limit = control[["limit"]], trace = control[["trace"]]))
  
  
  # Alternate between fitting a glm and \theta
  for (iter in seq.int(control[["iter.max"]])) {
    # Fit negative binomial model
    dev.old <- dev
    family <- negative.binomial(theta, link)
    fit <- feglmFit(beta, eta, y, X, A, B, lvls.k, family, control)
    beta <- fit[["coefficients"]]
    eta <- fit[["eta"]]
    dev <- fit[["deviance"]]
    theta <- suppressWarnings(theta.ml(y, family[["linkinv"]](eta), n = nobs[[4L]],
                                       limit = control[["limit"]], trace = control[["trace"]]))
    
    # Progress information
    if (control[["trace"]]) {
      cat("Deviance =", round(dev, 2L), "Iterations -", iter, "\n")
    }
    
    # Check termination condition
    if ((dev.old - dev) / (0.1 + dev) < control[["dev.tol"]]) {
      if (control[["trace"]]) {
        cat("Convergence\n")
      }
      break
    }
  }
  rm(y, X, eta, A, B)
  
  # Add names to \beta, Scores, and Hessian
  names(fit[["coefficients"]]) <- nms.sp
  colnames(fit[["Score"]]) <- nms.sp
  dimnames(fit[["Hessian"]]) <- list(nms.sp, nms.sp)
  
  # Return result list
  structure(c(fit, list(theta   = theta,
                        nobs    = nobs,
                        lvls.k  = lvls.k,
                        nms.fe  = nms.fe,
                        formula = formula,
                        data    = data,
                        family  = family,
                        control = control)), class = c("feglm", "feglm.nb"))
}
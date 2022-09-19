#' @title
#' Efficiently fit glm's with high-dimensional \eqn{k}-way fixed effects
#' @description
#' \code{\link{feglm}} can be used to fit generalized linear models with many high-dimensional fixed
#' effects. The estimation procedure is based on unconditional maximum likelihood and can be
#' interpreted as a \dQuote{weighted demeaning} approach that combines the work of Gaure (2013) and
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
#' Similar to \code{\link[stats]{glm.fit}} this has to be the result of a call to a family 
#' function. Default is \code{binomial()}. See \code{\link[stats]{family}} for details of family 
#' functions.
#' @param
#' weights an optional string with the name of the 'prior weights' variable in \code{data}.
#' @param
#' beta.start an optional vector of starting values for the structural parameters in the linear 
#' predictor. Default is \eqn{\boldsymbol{\beta} = \mathbf{0}}{\beta = 0}.
#' @param
#' eta.start an optional vector of starting values for the linear predictor.
#' @param
#' control a named list of parameters for controlling the fitting process. See 
#' \code{\link{feglmControl}} for details.
#' @details
#' If \code{\link{feglm}} does not converge this is often a sign of linear dependence between 
#' one or more regressors and a fixed effects category. In this case, you should carefully inspect 
#' your model specification.
#' @return
#' The function \code{\link{feglm}} returns a named list of class \code{"feglm"}.
#' @references
#' Gaure, S. (2013). "OLS with Multiple High Dimensional Category Variables". Computational
#' Statistics and Data Analysis, 66.
#' @references 
#' Marschner, I. (2011). "glm2: Fitting generalized linear models with convergence problems". 
#' The R Journal, 3(2).
#' @references 
#' Stammann, A., F. Heiss, and D. McFadden (2016). "Estimating Fixed Effects Logit Models with 
#' Large Panel Data". Working paper.
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
feglm <- function(
    formula       = NULL,
    data          = NULL,
    family        = binomial(),
    weights       = NULL,
    beta.start    = NULL,
    eta.start     = NULL,
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
  
  # Check validity of family
  # NOTE: Quasi families not supported yet
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
  if (family[["family"]] == "binomial") {
    # Check if 'y' is numeric
    if (data[, is.numeric(get(lhs))]) {
      # Check if 'y' is in [0, 1]
      if (data[, any(get(lhs) < 0.0 | get(lhs) > 1.0)]) {
        stop("Model response has to be within the unit interval.", call. = FALSE)
      }
    } else {
      # Check if 'y' is factor and transform otherwise
      data[, (1L) := checkFactor(get(lhs))]
      
      # Check if the number of levels equals two
      if (data[, length(levels(get(lhs)))] != 2L) {
        stop("Model response has to be binary.", call. = FALSE)
      }
      
      # Ensure 'y' is 0-1 encoded
      data[, (1L) := as.numeric(get(lhs)) - 1.0]
    }
  } else if (family[["family"]] %in% c("Gamma", "inverse.gaussian")) {
    # Check if 'y' is strictly positive
    if (data[, any(get(lhs) <= 0.0)]) {
      stop("Model response has to be strictly positive.", call. = FALSE)
    }
  } else {
    # Check if 'y' is positive
    if (data[, any(get(lhs) < 0.0)]) {
      stop("Model response has to be positive.", call. = FALSE)
    }
  }
  
  # Get names of the fixed effects variables and sort
  k.vars <- attr(terms(formula, rhs = 2L), "term.labels")
  k <- length(k.vars)
  setkeyv(data, k.vars)
  
  # Generate temporary variable for operations on the data.table
  tmp.var <- tempVar(data)
  
  # Drop observations that do not contribute to the log likelihood if requested 
  if (family[["family"]] %in% c("binomial", "poisson")) {
    if (control[["drop.pc"]]) {
      repeat {
        # Drop observations that do not contribute to the log-likelihood
        ncheck <- nrow(data)
        for (j in k.vars) {
          data[, (tmp.var) := mean(get(lhs)), by = eval(j)]
          if (family[["family"]] == "binomial") {
            data <- data[get(tmp.var) > 0.0 & get(tmp.var) < 1.0]
          } else {
            data <- data[get(tmp.var) > 0.0]
          }
          data[, (tmp.var) := NULL]
        }
        
        # Check termination
        if (ncheck == nrow(data)) break
      }
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
        stop("Length of 'beta.start' has to be equal to the number of structural parameters.", call. = FALSE)
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
    if (family[["family"]] == "binomial") {
      eta <- rep(family[["linkfun"]](sum(wt * (y + 0.5) / 2.0) / sum(wt)), nt)
      # eta <- rep(mean(family[["linkfun"]]((y + 0.5) / 2.0)), nt)
      # eta <- family[["linkfun"]]((y + 0.5) / 2.0) # GLM
    } else if (family[["family"]] %in% c("Gamma", "inverse.gaussian")) {
      eta <- rep(family[["linkfun"]](sum(wt * y) / sum(wt)), nt)
      # eta <- rep(mean(family[["linkfun"]](y)), nt)
      # eta <- family[["linkfun"]](y) # GLM
    } else {
      eta <- rep(family[["linkfun"]](sum(wt * (y + 0.1)) / sum(wt)), nt)
      # eta <- rep(mean(family[["linkfun"]](y + 0.1)), nt)
      # eta <- family[["linkfun"]](y + 0.1) # GLM
    }
  }
  rm(beta.start, eta.start)
  
  # Get names and number of levels in each fixed effects category
  nms.fe <- lapply(data[, k.vars, with = FALSE], levels)
  lvls.k <- sapply(nms.fe, length)
  
  # Generate auxiliary list of indexes for different sub panels
  k.list <- getIndexList(k.vars, data)
  
  # Fit generalized linear model
  fit <- feglmFit(beta, eta, y, X, wt, k.list, family, control)
  rm(y, X, eta)
  
  # Add names to \beta, Hessian, and MX (if provided)
  names(fit[["coefficients"]]) <- nms.sp
  if (control[["keep.mx"]]) {
    colnames(fit[["MX"]]) <- nms.sp
  }
  dimnames(fit[["Hessian"]]) <- list(nms.sp, nms.sp)
  
  # Generate result list
  reslist <- c(
    fit, list(
      nobs    = nobs,
      lvls.k  = lvls.k,
      nms.fe  = nms.fe,
      formula = formula,
      data    = data,
      family  = family,
      control = control
      )
  )
  
  # Return result list
  structure(reslist, class = "feglm")
}
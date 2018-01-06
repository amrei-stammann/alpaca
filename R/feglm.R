#' @title
#' feglm: A function to efficiently estimate glm's with high-dimensional 
#' \eqn{k}-way fixed effects
#' @description
#' \strong{Caution}: Package version 0.1 preliminary!
#' 
#' \code{feglm} is used to fit fixed effects generalized linear models with many
#' high-dimensional fixed effects based on unconditional maximum likelihood.
#' The algorithm can be interpreted as a pseudo demeaning approach that combines
#' the work of Gaure (2013) and Stammann et. al. (2016). For technical details
#' see Stammann (2018).
#' 
#' \strong{Remark:} The term fixed effect is used in econometrician's sense of
#' having intercepts for each level in each category.
#' @param
#' formula an object of class \code{"formula"}: a symbolic description of the
#' model to be fitted. \code{formula} must be of type \code{y ~ x | k},
#' where the \eqn{k} refers to an identifier of different levels in category
#' \eqn{k}.
#' @param
#' data an object of class \code{"data.frame"} containing the variables in the
#' model.
#' @param
#' beta.start an optional vector of starting values used for the structural
#' parameters in the linear predictor. Default see \code{Details}.
#' @param
#' D.alpha.start an optional vector of starting values used for the fixed
#' effects contribution in the linear predictor. Default see \code{Details}.
#' @param
#' family the description of the error distribution and link function to be used
#' in the model. This has to be a character string naming the
#' model function. The value should be any of \code{"logit"} or
#' \code{"poisson"}. Default is \code{"logit"}.
#' @param
#' control a named list of parameters for controlling the fitting process.
#' Defaults see \code{\link{alpaca.control}}.
#' @details 
#' ...
#' @return
#' The function \code{feglm} returns a named list of class \code{"alpaca"}.
#' @references
#' Gaure, S. (2013). "OLS with Multiple High Dimensional Category Variables".
#' Computational Statistics and Data Analysis. 66.
#' @references 
#' Stammann, A., F. Heiss, and D. McFadden (2016). 
#' "Estimating Fixed Effects Logit Models with Large Panel Data". Working paper.
#' @references 
#' Stammann, A. (2018). "Fast and Feasible Estimation of Generalized Linear
#' Models with High-Dimensional K-Way Fixed Effects". Working Paper.
#' @export
feglm <- function(formula = NULL,
                  data = NULL,
                  beta.start = NULL,
                  D.alpha.start = NULL,
                  family = c("logit", "poisson"),
                  control = list()) {
  # TODO: Checks
  # Validity of input argument (formula).
  if (is.null(formula)) {
    stop("'formula' has to be specified.")
  } else if (!inherits(formula, "formula")) {
    stop("'formula' has to be of class formula.")
  }
  
  # Validity of input argument (data).
  if (is.null(data)) {
    stop("'data' has to be specified.")
  } else if (!inherits(data, "data.frame")) {
    stop("'data' has to be of class data.frame.")
  }
  
  # Validity of input argument (family).
  # TODO: ...
  family <- match.arg(family)
  switch (family, logit = fam <- 0L, poisson = fam <- 2L)
  
  # Validity of input argument (control).
  if (!inherits(control, "list")) {
    stop("'control' has to be of class list.")
  }
  
  # Extract control list.
  ctrl <- do.call(alpaca.control, control)
  
  # Update formula and drop missing values.
  formula <- Formula(formula)
  mf <- model.frame(formula = formula, data = data)
  nobs.na <- length(attr(mf, "na.action"))
  
  # Extract data.
  y <- model.response(mf)
  X <- as.matrix(model.matrix(update(formula, . ~ . - 1), data = mf, rhs = 1))
  D <- model.part(formula, data = mf, rhs = 2)
  
  # Extract names of structural parameters.
  nms.sp <- attr(X, "dimnames")[[2]]
  
  # Check validity of D.
  if(ncol(D) < 1L) {
    stop("Fixed effects uncorrectly specified.")
  }
  
  # Validity of input argument (beta.start).
  if (!is.null(beta.start)) {
    if (length(beta.start) != ncol(X)) {
      stop("Length of 'beta.start' has to be equal to the number of structural
           parameters.")
    }
  } else {
    if (family %in% c("logit", "probit")) {
      beta.start <- coef(lm(I(log((y + 0.5) / 2.0)) ~ X - 1))
    } else if (family == "poisson") {
      beta.start <- coef(lm(I(log(y + 0.1)) ~ X - 1))
    }
  }
  
  # Drop perfectly classified observations.
  if (ctrl[["drop.pc"]] == TRUE) {
    for (k in seq(ncol(D))) {
      # Compute group means.
      mean.tab <- aggregate(y ~ D[, k], FUN = mean)
      
      # Drop observations that do not contribute to the loglikelihood.
      if (family %in% c("logit", "probit")) {
        mean.tab <- mean.tab[(mean.tab[, 2] > 0 && mean.tab[, 2] < 1), ]
        idx <- mean.tab[, 1]
        idx <- D[, k] %in% idx
        y <- y[idx]
        X <- as.matrix(X[idx, ])
        D <- as.data.frame(D[idx, ])
      } else if (family == "poisson") {
        mean.tab <- mean.tab[mean.tab[, 2] > 0.0, ]
        idx <- mean.tab[, 1]
        idx <- D[, k] %in% idx
        y <- y[idx]
        X <- as.matrix(X[idx, ])
        D <- as.data.frame(D[idx, ])
      }
    }
  }
  nobs.pc <- nrow(mf) - length(y)
  
  # Ensure factors are consectuive integers.
  D <- sapply(D, function(x) as.integer(factor(x)))

  # Number of levels of k categories.
  lvls.k <- sapply(seq(ncol(D)), function(x) max(D[, x]))
  
  # Validity of input argument (D.alpha.start).
  if (!is.null(D.alpha.start)) {
    if (length(D.alpha.start) != length(y)) {
      stop("Length of 'D.alpha.start' has to be equal to the number of
           observations.")
    }
  } else {
    D.alpha.start <- numeric(length(y))
    if (family == "poisson") {
      D.alpha.start <- numeric(length(y))
      for (k in seq(ncol(D))) {
        D.alpha.start <- D.alpha.start + log(ave(y, D[, k]) + 0.1)
      }
    }
  }
  
  # Maximize maximum likelihood.
  mod <- .feglm(y, X, D - 1L, lvls.k, beta.start, D.alpha.start, fam, ctrl)
  
  # Modify 'mod'.
  mod[["coefficients"]] <- as.vector(mod[["coefficients"]])
  names(mod[["coefficients"]]) <- nms.sp
  mod[["D.alpha"]] <- as.vector(mod[["D.alpha"]])
  mod[["b"]] <- as.vector(mod[["b"]])
  mod[["w.tilde"]] <- as.vector(mod[["w.tilde"]])
  mod[["gradient"]] <- as.vector(mod[["gradient"]])
  names(mod[["gradient"]]) <- nms.sp
  colnames(mod[["gradient.cont"]]) <- nms.sp
  rownames(mod[["Hessian"]]) <- nms.sp
  colnames(mod[["Hessian"]]) <- nms.sp
  
  # Extend 'mod'.
  mod[["nobs"]] <- length(y)
  mod[["nobs.na"]] <- nobs.na
  mod[["nobs.pc"]] <- nobs.pc
  mod[["lvls.k"]] <- lvls.k
  mod[["formula"]] <- formula
  mod[["data"]] <- list(y = y, X = X, D = D)
  mod[["family"]] <- family
  mod[["control"]] <- ctrl
  
  # Return result list.
  structure(mod, class = "alpaca")
}
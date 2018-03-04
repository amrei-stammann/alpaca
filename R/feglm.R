#' @title
#' feglm: A function to efficiently estimate glm's with high-dimensional 
#' \eqn{k}-way fixed effects
#' @description
#' \strong{Caution}: Package version 0.1.2 preliminary!
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
#' model function. The value should be any of \code{"logit"}, \code{"probit"},
#' or \code{"poisson"}. Default is \code{"logit"}.
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
#' Models with High-Dimensional k-Way Fixed Effects". Working Paper.
#' @export
feglm <- function(formula = NULL,
                  data = NULL,
                  beta.start = NULL,
                  D.alpha.start = NULL,
                  family = c("logit", "probit", "poisson"),
                  control = list()) {
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
  family <- match.arg(family)
  switch (family,
          logit = family.uint <- 0L,
          probit = family.uint <- 1L,
          poisson = family.uint <- 2L)
  
  # Validity of input argument (control).
  if (!inherits(control, "list")) {
    stop("'control' has to be of class list.")
  }
  
  # Extract control list.
  control <- do.call(alpaca.control, control)
  
  # Update formula and drop missing values.
  formula <- Formula(formula)
  mf <- model.frame(formula = formula, data = data)
  nobs.na <- length(attr(mf, "na.action"))
  
  # Extract data.
  y <- model.response(mf)
  # Bugfix:
  # Needed to ensure that factor variables are correclty dummy encoded.
  # Old:
  # X <- as.matrix(model.matrix(update(formula, . ~ . - 1), data = mf, rhs = 1))
  # X <- as.matrix(model.part(formula, data = mf, rhs = 1L))
  X <- model.matrix(formula, data = mf, rhs = 1L)[, - 1L, drop = FALSE]
  D <- model.part(formula, data = mf, rhs = 2L)
  
  # Extract names of structural parameters.
  nms.sp <- attr(X, "dimnames")[[2]]
  
  # Check validity of D.
  K <- ncol(D)
  if(K < 1L) {
    stop("Fixed effects uncorrectly specified.")
  }
  
  # Validity of input argument (beta.start).
  if (!is.null(beta.start)) {
    if (length(beta.start) != ncol(X)) {
      stop("Length of 'beta.start' has to be equal to the number of structural parameters.")
    }
  } else {
    # Changed:
    # Zero seems to work better.
    # Old:
    # if (family %in% c("logit", "probit")) {
    #   beta.start <- coef(lm(I(log((y + 0.5) / 2.0)) ~ X - 1))
    # } else if (family == "poisson") {
    #   beta.start <- coef(lm(I(log(y + 0.1)) ~ X - 1))
    # }
    beta.start <- numeric(ncol(X))
  }
  
  # Drop perfectly classified observations.
  if (control[["drop.pc"]] == TRUE) {
    for (k in seq(K)) {
      # Compute group means.
      mean.tab <- aggregate(y ~ D[[k]], FUN = mean)
      
      # Drop observations that do not contribute to the loglikelihood.
      if (family %in% c("logit", "probit")) {
        # Bugfix:
        # Dropping perfectly classified observations for binomial models should
        # now work as intended. Thanks to jmboehm@github
        # Old:
        # mean.tab <- mean.tab[(mean.tab[, 2L] > 0.0 && mean.tab[, 2L] < 1.0), ]
        mean.tab <- mean.tab[(mean.tab[, 2L] > 0.0 & mean.tab[, 2L] < 1.0), ]
        idx <- mean.tab[, 1L]
        idx <- D[[k]] %in% idx
        y <- y[idx]
        X <- X[idx, , drop = FALSE]
        D <- D[idx, , drop = FALSE]
      } else if (family == "poisson") {
        mean.tab <- mean.tab[mean.tab[, 2L] > 0.0, ]
        idx <- mean.tab[, 1L]
        idx <- D[[k]] %in% idx
        y <- y[idx]
        X <- X[idx, , drop = FALSE]
        D <- D[idx, , drop = FALSE]
      }
    }
  }
  nobs.pc <- nrow(mf) - length(y)
  
  # Ensure factors are consectuive integers.
  # Changes:
  # Some performance optimization and extract the names of the fixed effects.
  # Old:
  # D <- sapply(D, function(x) as.integer(factor(x)))
  # lvls.k <- sapply(seq(ncol(D)), function(x) max(D[, x]))
  nms.k <- colnames(D)
  D[nms.k] <- lapply(D[nms.k], factor)
  nms.fe <- lapply(D[nms.k], levels)
  lvls.k <- sapply(nms.fe, length)
  nms.fe <- unlist(nms.fe)
  D <- sapply(D, as.integer)
  
  # Validity of input argument (D.alpha.start).
  if (!is.null(D.alpha.start)) {
    if (length(D.alpha.start) != length(y)) {
      stop("Length of 'D.alpha.start' has to be equal to the number of observations.")
    }
  } else {
    # Changed:
    # Better starting values improve convergence of pseudo-demeaning.
    D.alpha.start <- numeric(length(y))
    if (family %in% c("logit", "probit")) {
      for (k in seq(K)) {
        D.alpha.start <- D.alpha.start + log(ave(y, D[, k]) + 0.5 / 2.0) / K
      }
    } else if (family == "poisson") {
      for (k in seq(K)) {
        # D.alpha.start <- D.alpha.start + log(ave(y, D[, k]) + 0.1)
        D.alpha.start <- D.alpha.start + log(ave(y, D[, k]) + 0.1) / K
      }
    }
  }
  
  # Maximize maximum likelihood.
  mod <- .feglm(y, X, D - 1L, lvls.k, beta.start, D.alpha.start, family.uint, control)
  
  # Modify 'mod'.
  mod[["coefficients"]] <- as.vector(mod[["coefficients"]])
  names(mod[["coefficients"]]) <- nms.sp
  mod[["D.alpha"]] <- as.vector(mod[["D.alpha"]])
  mod[["b"]] <- as.vector(mod[["b"]])
  mod[["w.tilde"]] <- as.vector(mod[["w.tilde"]])
  mod[["gradient"]] <- as.vector(mod[["gradient"]])
  names(mod[["gradient"]]) <- nms.sp
  dimnames(mod[["gradient.cont"]]) <- list(seq(length(y)), nms.sp)
  dimnames(mod[["Hessian"]]) <- list(nms.sp, nms.sp)
  
  # Extend 'mod'.
  mod[["nobs"]] <- length(y)
  mod[["nobs.na"]] <- nobs.na
  mod[["nobs.pc"]] <- nobs.pc
  mod[["lvls.k"]] <- lvls.k
  mod[["nms.fe"]] <- nms.fe
  mod[["formula"]] <- formula
  mod[["data"]] <- list(y = y, X = X, D = D)
  mod[["family"]] <- family
  mod[["control"]] <- control
  
  # Return result list.
  structure(mod, class = "alpaca")
}
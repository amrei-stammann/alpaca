#' @title
#' Extract estimates of average partial effects
#' @description
#' \code{\link{coef.APEs}} is a generic function which extracts estimates of the average partial 
#' effects from objects returned by \code{\link{getAPEs}}.
#' @param 
#' object an object of class \code{"APEs"}.
#' @param 
#' ... other arguments.
#' @return
#' The function \code{\link{coef.APEs}} returns a named vector of estimates of the average partial 
#' effects.
#' @seealso
#' \code{\link{getAPEs}}
#' @export
coef.APEs <- function(object, ...) {
  object[["delta"]]
}
  

#' @title
#' Extract estimates of structural parameters
#' @description
#' \code{\link{coef.feglm}} is a generic function which extracts estimates of the structural parameters
#' from objects returned by \code{\link{feglm}}.
#' @param 
#' object an object of class \code{"feglm"}.
#' @param 
#' ... other arguments.
#' @return
#' The function \code{\link{coef.feglm}} returns a named vector of estimates of the structural 
#' parameters.
#' @seealso
#' \code{\link{feglm}}
#' @export
coef.feglm <- function(object, ...) {
  object[["coefficients"]]
}


#' @title
#' Extract coefficient matrix of average partial effects
#' @description
#' \code{\link{coef.summary.APEs}} is a generic function which extracts a coefficient matrix of 
#' average partial effects from objects returned by \code{\link{getAPEs}}.
#' @param 
#' object an object of class \code{"summary.APEs"}.
#' @param 
#' ... other arguments.
#' @return
#' The function \code{\link{coef.summary.APEs}} returns a named matrix of estimates related to the
#' average partial effects.
#' @seealso
#' \code{\link{getAPEs}}
#' @export
coef.summary.APEs <- function(object, ...) {
  object[["cm"]]
}


#' @title
#' Extract coefficient matrix of structural parameters
#' @description
#' \code{\link{coef.summary.feglm}} is a generic function which extracts a coefficient matrix of 
#' structural parameters from objects returned by \code{\link{feglm}}.
#' @param 
#' object an object of class \code{"summary.feglm"}.
#' @param 
#' ... other arguments.
#' @return
#' The function \code{\link{coef.summary.feglm}} returns a named matrix of estimates related to the
#' structural parameters.
#' @seealso
#' \code{\link{feglm}}
#' @export
coef.summary.feglm <- function(object, ...) {
  object[["cm"]]
}


#' @title
#' Extract \code{feglm} fitted values 
#' @description
#' \code{\link{fitted.feglm}} is a generic function which extracts fitted values from an object 
#' returned by \code{\link{feglm}}.
#' @param 
#' object an object of class \code{"feglm"}.
#' @param 
#' ... other arguments.
#' @return
#' The function \code{\link{fitted.feglm}} returns a vector of fitted values.
#' @seealso
#' \code{\link{feglm}}
#' @export
fitted.feglm <- function(object, ...) {
  object[["family"]][["linkinv"]](object[["eta"]])
}


#' @title
#' Predict method for \code{feglm} fits
#' @description
#' \code{\link{predict.feglm}} is a generic function which obtains predictions from an object 
#' returned by \code{\link{feglm}}.
#' @param 
#' object an object of class \code{"feglm"}.
#' @param
#' type the type of prediction required. \code{"link"} is on the scale of the linear predictor
#' whereas \code{"response"} is on the scale of the response variable. Default is \code{"link"}.
#' @param 
#' ... other arguments.
#' @return
#' The function \code{\link{predict.feglm}} returns a vector of predictions.
#' @seealso
#' \code{\link{feglm}}
#' @export
predict.feglm <- function(object, type = c("link", "response"), ...) {
  # Check validity of 'type'
  type <- match.arg(type)
  
  # Compute requested type of prediction
  x <- object[["eta"]]
  if (type == "response") {
    x <- object[["family"]][["linkinv"]](x)
  }
  
  # Return prediction
  x
}


#' @title
#' Print \code{APEs}
#' @description
#' \code{\link{print.APEs}} is a generic function which displays some minimal information from 
#' objects returned by \code{\link{getAPEs}}.
#' @param 
#' x an object of class \code{"APEs"}.
#' @param 
#' digits unsigned integer indicating the number of decimal places. Default is 
#' \code{max(3L, getOption("digits") - 3L)}.
#' @param 
#' ... other arguments.
#' @seealso
#' \code{\link{getAPEs}}
#' @export
print.APEs <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  print(x[["delta"]], digits = digits)
}


#' @title
#' Print \code{feglm}
#' @description
#' \code{\link{print.feglm}} is a generic function which displays some minimal information from 
#' objects returned by \code{\link{feglm}}.
#' @param 
#' x an object of class \code{"feglm"}.
#' @param 
#' digits unsigned integer indicating the number of decimal places. Default is 
#' \code{max(3L, getOption("digits") - 3L)}.
#' @param 
#' ... other arguments.
#' @seealso
#' \code{\link{feglm}}
#' @export
print.feglm <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat(sub("\\(.*\\)", "", x[["family"]][["family"]]), " - ",
      x[["family"]][["link"]], " link",
      ", l= [", paste0(x[["lvls.k"]], collapse = ", "), "]\n\n", sep = "")
  print(x[["coefficients"]], digits = digits)
}


#' @title
#' Print \code{summary.APEs}
#' @description
#' \code{\link{print.summary.APEs}} is a generic function which displays summary statistics from 
#' objects returned by \code{\link{summary.APEs}}.
#' @param 
#' x an object of class \code{"summary.APEs"}.
#' @param 
#' digits unsigned integer indicating the number of decimal places. Default is 
#' \code{max(3L, getOption("digits") - 3L)}.
#' @param 
#' ... other arguments.
#' @seealso
#' \code{\link{getAPEs}}
#' @export
print.summary.APEs <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("Estimates:\n")
  printCoefmat(x[["cm"]], P.values = TRUE, has.Pvalue = TRUE, digits = digits)
}


#' @title
#' Print \code{summary.feglm}
#' @description
#' \code{\link{print.summary.feglm}} is a generic function which displays summary statistics from 
#' objects returned by \code{\link{summary.feglm}}.
#' @param 
#' x an object of class \code{"summary.feglm"}.
#' @param 
#' digits unsigned integer indicating the number of decimal places. Default is 
#' \code{max(3L, getOption("digits") - 3L)}.
#' @param 
#' ... other arguments.
#' @seealso
#' \code{\link{feglm}}
#' @export
print.summary.feglm <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat(sub("\\(.*\\)", "", x[["family"]][["family"]]), " - ",
      x[["family"]][["link"]], " link\n\n", sep = "")
  print(x[["formula"]])
  cat("\nEstimates:\n")
  printCoefmat(x[["cm"]], P.values = TRUE, has.Pvalue = TRUE, digits = digits)
  cat("\nresidual deviance= ",
      format(x[["deviance"]], digits = max(5L, digits + 1L), nsmall = 2L),
      ",\n", sep = "")
  cat("null deviance= ",
      format(x[["null.deviance"]], digits = max(5L, digits + 1L), nsmall = 2L),
      ",\n", sep = "")
  cat("n= ", x[["nobs"]][["nobs"]],
      ", l= [", paste0(x[["lvls.k"]], collapse = ", "), "]\n", sep = "")
  if (x[["nobs"]][["nobs.na"]] > 0L | x[["nobs"]][["nobs.pc"]] > 0L) {
    cat("\n")
    if (x[["nobs"]][["nobs.na"]] > 0L) {
      cat("(", x[["nobs"]][["nobs.na"]], "observation(s) deleted due to missingness )\n")
    }
    if (x[["nobs"]][["nobs.pc"]] > 0L) {
      cat("(", x[["nobs"]][["nobs.pc"]], "observation(s) deleted due to perfect classification )\n")
    }
  }
  cat("\nNumber of Fisher Scoring Iterations:", x[["iter"]], "\n")
  if (!is.null(x[["theta"]])) {
    cat("\ntheta= ",
        format(x[["theta"]], digits = digits, nsmall = 2L),
        ", std. error= ",
        format(attr(x[["theta"]], "SE"), digits = digits, nsmall = 2L),
        "\n", sep = "")
  }
}


#' @title
#' Summarizing models of class \code{APEs}
#' @description
#' Summary statistics for objects of class \code{"APEs"}.
#' @param 
#' object an object of class \code{"APEs"}.
#' @param 
#' ... other arguments.
#' @return
#' Returns an object of class \code{"summary.APEs"} which is a list of summary statistics of 
#' \code{object}.
#' @seealso
#' \code{\link{getAPEs}}
#' @export
summary.APEs <- function(object, ...) {
  # Compute coefficent matrix
  est <- object[["delta"]]
  se <- sqrt(diag(object[["vcov"]]))
  z <- est / se
  p <- 2.0 * pnorm(- abs(z))
  cm <- cbind(est, se, z, p)  
  rownames(cm) <- names(est)
  colnames(cm) <- c("Estimate", "Std. error", "z value", "Pr(> |z|)")
  
  # Return coefficient matrix
  structure(list(cm = cm), class = "summary.APEs")
}


#' @title
#' Summarizing models of class \code{feglm}
#' @description
#' Summary statistics for objects of class \code{"feglm"}.
#' @param 
#' object an object of class \code{"feglm"}.
#' @param
#' type the type of covariance estimate required. \code{"hessian"} refers to the inverse
#' of the negative expected Hessian after convergence and is the default option. 
#' \code{"outer.product"} is the outer-product-of-the-gradient estimator, 
#' \code{"sandwich"} is the sandwich estimator (sometimes also refered as robust estimator), 
#' and \code{"clustered"} computes a clustered covariance matrix given some cluster variables.
#' @param
#' cluster a symbolic description indicating the clustering of observations.
#' @param
#' cluster.vars deprecated; use \code{cluster} instead.
#' @param 
#' ... other arguments.
#' @details
#' Multi-way clustering is done using the algorithm of Cameron, Gelbach, and Miller (2011). An 
#' example is provided in the vignette "Replicating an Empirical Example of International Trade".
#' @return
#' Returns an object of class \code{"summary.feglm"} which is a list of summary statistics of
#' \code{object}.
#' @references
#' Cameron, C., J. Gelbach, and D. Miller (2011). "Robust Inference With Multiway Clustering".
#' Journal of Business & Economic Statistics 29(2).
#' @seealso
#' \code{\link{feglm}}
#' @export
summary.feglm <- function(object,
                          type         = c("hessian", "outer.product", "sandwich", "clustered"),
                          cluster      = NULL,
                          cluster.vars = NULL,
                          ...) {
  # 'cluster.vars' is deprecated
  if (!is.null(cluster.vars)) {
    warning("'cluster.vars' is deprecated; please use 'cluster' instead.", call. = FALSE)
    if (!is.character(cluster.vars)) {
      stop("'cluster.vars' has to be of class character.", call. = FALSE)
    }
    cluster <- as.formula(paste0("~", paste0(cluster.vars, collapse = "+")))
  }
  
  # Compute coefficent matrix
  est <- object[["coefficients"]]
  se <- sqrt(diag(vcov(object, type, cluster)))
  z <- est / se
  p <- 2.0 * pnorm(- abs(z))
  cm <- cbind(est, se, z, p)  
  rownames(cm) <- names(est)
  colnames(cm) <- c("Estimate", "Std. error", "z value", "Pr(> |z|)")
  
  # Generate result list
  res <- list(cm            = cm, 
              deviance      = object[["deviance"]],
              null.deviance = object[["null.deviance"]],
              iter          = object[["iter"]],
              nobs          = object[["nobs"]],
              lvls.k        = object[["lvls.k"]],
              formula       = object[["formula"]],
              family        = object[["family"]])
  if (inherits(object, "feglm.nb")) {
    res[["theta"]] <- object[["theta"]]
  }
  
  # Return list
  structure(res, class = "summary.feglm")
}


#' @title
#' Extract estimates of the covariance matrix
#' @description
#' \code{\link{vcov.feglm}} computes an estimate of the covariance matrix of the estimator of the
#' structural parameters from objects returned by \code{\link{feglm}}. The estimate is obtained
#' using the Hessian, the scores, or a combination of boths after convergence.
#' @param 
#' object an object of class \code{"feglm"}.
#' @param
#' type the type of covariance estimate required. \code{"hessian"} refers to the inverse
#' of the negative expected Hessian after convergence and is the default option. 
#' \code{"outer.product"} is the outer-product-of-the-gradient estimator, 
#' \code{"sandwich"} is the sandwich estimator (sometimes also refered as robust estimator), 
#' and \code{"clustered"} computes a clustered covariance matrix given some cluster variables.
#' @param
#' cluster a symbolic description indicating the clustering of observations.
#' @param
#' cluster.vars deprecated; use \code{cluster} instead.
#' @param 
#' ... other arguments.
#' @details
#' Multi-way clustering is done using the algorithm of Cameron, Gelbach, and Miller (2011). An 
#' example is provided in the vignette "Replicating an Empirical Example of International Trade".
#' @return
#' The function \code{\link{vcov.feglm}} returns a named matrix of covariance estimates.
#' @references
#' Cameron, C., J. Gelbach, and D. Miller (2011). "Robust Inference With Multiway Clustering".
#' Journal of Business & Economic Statistics 29(2).
#' @seealso
#' \code{\link{feglm}}
#' @export
vcov.feglm <- function(object,
                       type         = c("hessian", "outer.product", "sandwich", "clustered"),
                       cluster      = NULL,
                       cluster.vars = NULL,
                       ...) {
  # Check validity of input argument 'type'
  type <- match.arg(type)
  
  # 'cluster.vars' is deprecated
  if (!is.null(cluster.vars)) {
    warning("'cluster.vars' is deprecated; please use 'cluster' instead.", call. = FALSE)
    if (!is.character(cluster.vars)) {
      stop("'cluster.vars' has to be a character.", call. = FALSE)
    }
    cluster <- as.formula(paste0("~", paste0(cluster.vars, collapse = "+")))
  }
  
  # Compute requested type of estimated covariance matrix
  p <- length(object[["coefficients"]])
  if (type == "hessian") {
    # Check if the Hessian is invertible and compute its inverse
    R <- try(chol(object[["Hessian"]]), silent = TRUE)
    if (inherits(R, "try-error")) {
      V <- matrix(Inf, p, p)
    } else {
      V <- chol2inv(R)
    }
  } else if (type == "outer.product") {
    # Check if the OPG is invertible and compute its inverse
    R <- try(chol(crossprod(object[["Score"]])), silent = TRUE)
    if (inherits(R, "try-error")) {
      V <- matrix(Inf, p, p)
    } else {
      V <- chol2inv(R)
    }
  } else {
    # Check if the Hessian is invertible and compute its inverse
    R <- try(chol(object[["Hessian"]]), silent = TRUE)
    if (inherits(R, "try-error")) {
      V <- matrix(Inf, p, p)
    } else {
      # Extract the score and compute the inverse of the empirical Hessian
      G <- object[["Score"]]
      A <- chol2inv(R)
      
      # Compute inner part of the sandwich formula
      if (type == "sandwich") {
        B <- crossprod(G)
      } else {
        # Check validity of input argument 'cluster'
        if (is.null(cluster)) {
          stop("'cluster' has to be specified.", call. = FALSE)
        } else if (!inherits(cluster, "formula")) {
          stop("'cluster' has to be of class formula.", call. = FALSE)
        }
        
        # Extract cluster variables
        cluster <- Formula(cluster)
        D <- try(object[["data"]][, all.vars(cluster), with = FALSE], silent = TRUE)
        if (inherits(D, "try-error")) {
          stop(paste("At least one cluster variable was not found.",
                     "Ensure to pass variables that are not part of the model itself, but are", 
                     "required to compute clustered standard errors, to 'feglm'.", 
                     "This can be done via 'formula'. See documentation for details."),
               call. = FALSE)
        }
        
        # Ensure cluster variables are factors
        cl.vars <- names(D)
        D[, (cl.vars) := lapply(.SD, checkFactor)]
        
        # Join cluster variables and scores
        sp.vars <- colnames(G)
        G <- cbind(D, G)
        rm(D)
        
        # Multiway clustering ala Cameron, Gelbach, and Miller (2011)
        setkeyv(G, cl.vars)
        B <- matrix(0.0, p, p)
        for (i in seq.int(length(cl.vars))) {
          # Compute outer product for all possible combinations
          cl.combn <- combn(cl.vars, i)
          B.r <- matrix(0.0, p, p)
          for (j in seq.int(ncol(cl.combn))) {
            cl <- cl.combn[, j]
            B.r <- B.r + 
              crossprod(as.matrix(
                G[, lapply(.SD, sum), by = cl, .SDcols = sp.vars][, sp.vars, with = FALSE]))
          }
          
          # Update outer product
          if (i %% 2L) {
            B <- B + B.r
          } else {
            B <- B - B.r
          }
        }
      }
      
      # Sandwich formula
      V <- A %*% B %*% A
    }
  }
  
  # Return covariance estimate
  V
}
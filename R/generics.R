#' @title
#' Extract estimates of structural parameters
#' @description
#' \code{coef.feglm} is a generic function which extracts estimates of the structural parameters
#' from objects returned by \code{feglm}.
#' @param 
#' object an object of class \code{feglm}.
#' @param 
#' ... other arguments.
#' @return
#' The function \code{coef.feglm} returns a named vector of estimates refering to the structural
#' parameters.
#' @seealso
#' \code{\link{feglm}}
#' @export
coef.feglm <- function(object, ...) object[["coefficients"]]


#' @title
#' Extract coefficient matrix of structural parameters
#' @description
#' \code{coef.summary.feglm} is a generic function which extracts a coefficient matrix of structural
#' parameters from objects returned by \code{feglm}.
#' @param 
#' object an object of class \code{summary.feglm}.
#' @param 
#' ... other arguments.
#' @return
#' The function \code{coef.summary.feglm} returns a named matrix of estimates related to the
#' structural parameters.
#' @seealso
#' \code{\link{feglm}}
#' @export
coef.summary.feglm <- function(object, ...) object[["cm"]]


#' @title
#' Extract \code{feglm} fitted values 
#' @description
#' \code{fitted.feglm} is a generic function which extracts fitted values from an object returned by
#' \code{feglm}.
#' @param 
#' object an object of class \code{feglm}.
#' @param 
#' ... other arguments.
#' @return
#' The function \code{fitted.feglm} returns a vector of fitted values.
#' @seealso
#' \code{\link{feglm}}
#' @export
fitted.feglm <- function(object, ...) predict(object, type = "response")


#' @title
#' Predict method for \code{feglm} fits
#' @description
#' \code{predict.feglm} is a generic function which obtains predictions from an object returned by
#' \code{feglm}.
#' @param 
#' object an object of class \code{feglm}.
#' @param
#' type the type of prediction required. \code{"link"} is on the scale to the linear predictor
#' whereas \code{"response"} is on the scale of the response variable. Default is \code{"link"}.
#' @param 
#' ... other arguments.
#' @return
#' The function \code{predict.feglm} returns a vector of predictions.
#' @seealso
#' \code{\link{feglm}}
#' @export
predict.feglm <- function(object, type = c("link", "response"), ...) {
  # Check validity of 'type'
  type <- match.arg(type)
  
  # Extract regressor matrix
  X <- model.matrix(object[["formula"]], object[["data"]], rhs = 1L)[, - 1L, drop = FALSE]
  dimnames(X) <- NULL
  
  # Compute requested type of prediction
  x <- as.vector(X %*% coef(object) + object[["D.alpha"]])
  if (type == "response") {
    x <- object[["family"]][["linkinv"]](x)
  }
  
  # Return prediction
  x
}


#' @title
#' Print \code{feglm}
#' @description
#' \code{print.feglm} is a generic function which displays some minimal information from objects
#' returned by \code{feglm}.
#' @param 
#' x an object of class \code{feglm}.
#' @param 
#' digits unsigned integer indicating the number of decimal places. Default is 
#' \code{max(3L, getOption("digits") - 3L)}.
#' @param 
#' ... other arguments.
#' @seealso
#' \code{\link{feglm}}
#' @export
print.feglm <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat(x[["family"]][["family"]],
      ", l= [", paste0(x[["lvls.k"]], collapse = ", "), "]\n\n", sep = "")
  print(coef(x), digits = digits)
}


#' @title
#' Print \code{summary.feglm}
#' @description
#' \code{print.summary.feglm} is a generic function which displays summary statistics from objects
#' returned by \code{summary.feglm}.
#' @param 
#' x an object of class \code{summary.feglm}.
#' @param 
#' digits unsigned integer indicating the number of decimal places. Default is 
#' \code{max(3L, getOption("digits") - 3L)}.
#' @param 
#' ... other arguments.
#' @seealso
#' \code{\link{feglm}}
#' @export
print.summary.feglm <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat(x[["family"]][["family"]], "\n\n")
  print(x[["formula"]])
  cat("\nl= [", paste0(x[["lvls.k"]], collapse = ", "), "], n= ", x[["nobs"]],
      ", deviance= ", round(x[["deviance"]], digits = digits), "\n", sep = "")
  cat("\nStructural parameter(s):\n\n")
  printCoefmat(x[["cm"]], P.values = TRUE, has.Pvalue = TRUE,
               digits = digits)
  if(x[["nobs.na"]] > 0L) {
    cat("(", x[["nobs.na"]], "observation(s) deleted due to missingness )\n")
  }
  if(x[["nobs.pc"]] > 0L) {
    cat("(", x[["nobs.pc"]],
        "observation(s) deleted due to perfect classification )\n")
  }
}


#' @title
#' Summarizing models of class \code{feglm}
#' @description
#' Summary statistics for objects of class \code{feglm}.
#' @param 
#' object an object of class \code{feglm}.
#' @param
#' type the type of covariance estimate required. \code{"empirical.hessian"} refers to the inverse
#' of the negative Hessian after convergence and is the default option. \code{"outer.product"} is
#' the outer-product-of-the-gradient estimator, \code{"sandwich"} is the sandwich estimator 
#' (sometimes also refered as robust estimator), and \code{"clustered"} computes a clustered
#' covariance matrix given some cluster variables.
#' @param
#' cluster.vars a character or vector of characters indicating the names of the cluster variables in
#' the data set.
#' @param 
#' ... other arguments.
#' @details
#' Multi-way clustering is done using the algorithm of Cameron, Gelbach, and Miller (2011).
#' @return
#' Returns an object of class \code{summary.feglm} which is a list of summary statistics of
#' \code{object}.
#' @references
#' Cameron, C., J. Gelbach, and D. Miller (2011). "Robust Inference With Multiway Clustering".
#' Journal of Business & Economic Statistics 29(2).
#' @seealso
#' \code{\link{feglm}}
#' @export
summary.feglm <- function(object,
                          type         = c("empirical.hessian", "outer.product",
                                           "sandwich", "clustered"),
                          cluster.vars = NULL,
                          ...) {
  # Compute coefficent matrix
  cm.header <- c("Estimate", "Std. error", "z value", "Pr(> |z|)")
  beta.hat <- coef(object)
  se.beta <- sqrt(diag(vcov(object, type, cluster.vars)))
  z.score <- beta.hat / se.beta
  p.value <- 2.0 * pnorm(- abs(z.score))
  cm <- cbind(beta.hat, se.beta, z.score, p.value)  
  rownames(cm) <- names(beta.hat)
  colnames(cm) <- cm.header
  
  # Return list
  structure(list(cm       = cm, 
                 maximum  = object[["maximum"]],
                 deviance = object[["deviance"]],
                 nobs     = object[["nobs"]],
                 lvls.k   = object[["lvls.k"]],
                 nobs.na  = object[["nobs.na"]],
                 nobs.pc  = object[["nobs.pc"]],
                 formula  = object[["formula"]],
                 family   = object[["family"]]),
            class = "summary.feglm")
}


#' @title
#' Extract estimates of the covariance matrix
#' @description
#' \code{vcov.feglm} computes an estimate of the covariance matrix of the estimator of the
#' structural parameters from objects returned by \code{\link{feglm}}. The estimate is obtained
#' using the Hessian, the scores, or a combination of boths after convergence.
#' @param 
#' object an object of class \code{feglm}.
#' @param
#' type the type of covariance estimate required. \code{"empirical.hessian"} refers to the inverse
#' of the negative Hessian after convergence and is the default option. \code{"outer.product"} is
#' the outer-product-of-the-gradient estimator, \code{"sandwich"} is the sandwich estimator 
#' (sometimes also refered as robust estimator), and \code{"clustered"} computes a clustered
#' covariance matrix given some cluster variables.
#' @param
#' cluster.vars a character or vector of characters indicating the names of the cluster variables in
#' the data set.
#' @param 
#' ... other arguments.
#' @details
#' Multi-way clustering is done using the algorithm of Cameron, Gelbach, and Miller (2011).
#' @return
#' The function \code{\link{vcov.feglm}} returns a named matrix of covariance estimates.
#' @references
#' Cameron, C., J. Gelbach, and D. Miller (2011). "Robust Inference With Multiway Clustering".
#' Journal of Business & Economic Statistics 29(2).
#' @seealso
#' \code{\link{feglm}}
#' @export
vcov.feglm <- function(object,
                       type         = c("empirical.hessian", "outer.product",
                                        "sandwich", "clustered"),
                       cluster.vars = NULL,
                       ...) {
  # Check validity of input argument 'type'
  type <- match.arg(type)
  
  # Compute requested type of estimated covariance matrix
  if (type == "empirical.hessian") {
    # Compute eigenvalues to check if the Hessian is invertible and compute its inverse
    H <- object[["Hessian"]]
    ev <- abs(eigen(H, symmetric = TRUE, only.values = TRUE)[["values"]])
    if (min(ev) > .Machine[["double.eps"]] * max(ev) * 10.0) {
      V <- solve(H)
    } else {
      p <- nrow(H)
      V <- matrix(Inf, p, p)
    }
  } else if (type == "outer.product") {
    # Compute eigenvalues to check if the OPG is invertible and compute its inverse
    G <- object[["Score"]]
    B <- crossprod(G)
    ev <- abs(eigen(B, symmetric = TRUE, only.values = TRUE)[["values"]])
    if (min(ev) > .Machine[["double.eps"]] * max(ev) * 10.0) {
      V <- solve(B)
    } else {
      p <- nrow(B)
      V <- matrix(Inf, p, p)
    }
  } else {
    # Extract the score and compute the inverse of the empirical Hessian
    G <- object[["Score"]]
    H <- object[["Hessian"]]
    ev <- abs(eigen(H, symmetric = TRUE, only.values = TRUE)[["values"]])
    if (min(ev) > .Machine[["double.eps"]] * max(ev) * 10.0) {
      A <- solve(H)
    } else {
      p <- nrow(H)
      A <- matrix(Inf, p, p)
    }
    
    # Compute inner part of the sandwich formula
    if (type == "sandwich") {
      B <- crossprod(G)
    } else {
      # Check validity of 'cluster.vars'
      if (is.null(cluster.vars)) {
        stop("'cluster.vars' has to be specified.")
      } else {
        if (any(!cluster.vars %in% names(object[["data"]]))){
          stop("At least one cluster variable was not found in the data.")
        }
      }
      
      # Extract cluster variables
      D <- object[["data"]][, cluster.vars, with = FALSE]
      
      # Join cluster variables and scores
      p <- ncol(G)
      sp.vars <- colnames(G)
      G <- cbind(D, G)
      
      # Start clustering
      cl <- length(cluster.vars)
      B <- matrix(0.0, p, p)
      for (i in seq(cl)) {
        # Compute inner part of the sandwich formula for all possible combinations
        combns <- combn(cluster.vars, i)
        B.r <- lapply(seq(ncol(combns)), function(x) {
          G.r <- G[, lapply(.SD, sum), by = eval(combns[, x]), .SDcols = sp.vars]
          crossprod(as.matrix(G.r[, sp.vars, with = FALSE]))
        })
        
        # Multiway clustering
        # See Cameron, Gelbach, and Miller (2011), "Robust Inference With Multiway Clustering"
        if (i %% 2L) {
          B <- B + Reduce("+", B.r)
        } else {
          B <- B - Reduce("+", B.r)
        }
      }
    }
    
    # Sandwich formula
    V <- A %*% B %*% A
  }
  
  # Return covariance estimate
  V
}


# Unload
.onUnload <- function(libpath) library.dynam.unload("alpaca", libpath)
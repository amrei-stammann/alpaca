### Internal functions (not exported)


# Checks if variable is a factor and transforms if necessary
checkFactor <- function(x) {
  if (is.factor(x)) {
    droplevels(x)
  } else {
    factor(x)
  }
}


# Fitting algorithm (similar to glm.fit)
feglmFit <- function(beta, eta, y, X, wt, k.list, family, control) {
  # Extract control arguments
  center.tol <- control[["center.tol"]]
  dev.tol <- control[["dev.tol"]]
  epsilon <- max(min(1.0e-07, dev.tol / 1000.0), .Machine[["double.eps"]])
  iter.max <- control[["iter.max"]]
  trace <- control[["trace"]]
  keep.mx <- control[["keep.mx"]]
  
  # Compute initial quantities for the maximization routine
  nt <- length(y)
  mu <- family[["linkinv"]](eta)
  dev <- sum(family[["dev.resids"]](y, mu, wt))
  null.dev <- sum(family[["dev.resids"]](y, mean(y), wt))
  
  # Generate temporary variables
  Mnu <- as.matrix(numeric(nt))
  MX <- X
  
  # Start maximization of the log-likelihood
  conv <- FALSE
  for (iter in seq.int(iter.max)) {
    # Store \eta, \beta, and deviance of the previous iteration
    eta.old <- eta
    beta.old <- beta
    dev.old <- dev
    
    # Compute weights and dependent variable
    mu.eta <- family[["mu.eta"]](eta)
    w <- (wt * mu.eta^2) / family[["variance"]](mu)
    w.tilde <- sqrt(w)
    nu <- (y - mu) / mu.eta
    
    # Centering variables
    Mnu <- centerVariables((Mnu + nu), w, k.list, center.tol)
    MX <- centerVariables(MX, w, k.list, center.tol)
    
    # Compute update step and update \eta
    beta.upd <- as.vector(qr.solve(MX * w.tilde, Mnu * w.tilde, epsilon))
    eta.upd <- nu - as.vector(Mnu - MX %*% beta.upd)
    
    # Step-halving with three checks
    # 1. finite deviance
    # 2. valid \eta and \mu
    # 3. improvement as in glm2
    rho <- 1.0
    for (inner.iter in seq.int(50L)) {
      eta <- eta.old + rho * eta.upd
      beta <- beta.old + rho * beta.upd
      mu <- family[["linkinv"]](eta)
      dev <- sum(family[["dev.resids"]](y, mu, wt))
      dev.crit <- is.finite(dev)
      val.crit <- family[["valideta"]](eta) && family[["validmu"]](mu)
      imp.crit <- (dev - dev.old) / (0.1 + abs(dev)) <= - dev.tol
      if (dev.crit && val.crit && imp.crit) break
      rho <- rho / 2.0
    }
    
    # Check if step-halving failed (deviance and invalid \eta or \mu)
    if (!dev.crit || !val.crit) {
      stop("Inner loop failed; cannot correct step size.", call. = FALSE)
    }
    
    # Stop if we do not improve
    if (!imp.crit) {
      eta <- eta.old
      beta <- beta.old
      dev <- dev.old
      mu <- family[["linkinv"]](eta)
    }
    
    # Progress information
    if (trace) {
      cat("Deviance=", format(dev, digits = 5L, nsmall = 2L), "Iterations -", iter, "\n")
      cat("Estimates=", format(beta, digits = 3L, nsmall = 2L), "\n")
    }
    
    # Check convergence
    dev.crit <- abs(dev - dev.old) / (0.1 + abs(dev))
    if (trace) cat("Stopping criterion=", dev.crit, "\n")
    if (dev.crit < dev.tol) {
      if (trace) cat("Convergence\n")
      conv <- TRUE
      break
    }
    
    # Update starting guesses for acceleration
    Mnu <- Mnu - nu
  }
  
  # Information if convergence failed
  if (!conv && trace) cat("Algorithm did not converge.\n") 
  
  # Update weights and dependent variable
  mu.eta <- family[["mu.eta"]](eta)
  w <- (wt * mu.eta^2) / family[["variance"]](mu)
  
  # Center variables
  MX <- centerVariables(X, w, k.list, center.tol)
  
  # Recompute Hessian
  H <- crossprod(MX * sqrt(w))
  
  # Generate result list
  reslist <- list(
    coefficients  = beta,
    eta           = eta,
    weights       = wt,
    Hessian       = H,
    deviance      = dev,
    null.deviance = null.dev,
    conv          = conv,
    iter          = iter
    )
  
  # Update result list
  if (keep.mx) reslist[["MX"]] <- MX
  
  # Return result list
  reslist
}


# Efficient offset algorithm to update the linear predictor
feglmOffset <- function(object, offset) {
  # Check validity of 'object'
  if (!inherits(object, "feglm")) {
    stop("'feglmOffset' called on a non-'feglm' object.")
  }

  # Extract required quantities from result list
  control <- object[["control"]]
  data <- object[["data"]]
  wt <- object[["weights"]]
  family <- object[["family"]]
  formula <- object[["formula"]]
  lvls.k <- object[["lvls.k"]]
  nt <- object[["nobs"]][["nobs"]]
  k.vars <- names(lvls.k)
  
  # Extract dependent variable
  y <- data[[1L]]
  
  # Extract control arguments
  center.tol <- control[["center.tol"]]
  dev.tol <- control[["dev.tol"]]
  iter.max <- control[["iter.max"]]
  
  # Generate auxiliary list of indexes to project out the fixed effects
  k.list <- getIndexList(k.vars, data)
  
  # Compute starting guess for \eta
  if (family[["family"]] == "binomial") {
    eta <- rep(family[["linkfun"]](sum(wt * (y + 0.5) / 2.0) / sum(wt)), nt)
    # eta <- rep(mean(family[["linkfun"]]((y + 0.5) / 2.0)), nt)
  } else if (family[["family"]] %in% c("Gamma", "inverse.gaussian")) {
    eta <- rep(family[["linkfun"]](sum(wt * y) / sum(wt)), nt)
    # eta <- rep(mean(family[["linkfun"]](y)), nt)
  } else {
    eta <- rep(family[["linkfun"]](sum(wt * (y + 0.1)) / sum(wt)), nt)
    # eta <- rep(mean(family[["linkfun"]](y + 0.1)), nt)
  }

  # Compute initial quantities for the maximization routine
  mu <- family[["linkinv"]](eta)
  dev <- sum(family[["dev.resids"]](y, mu, wt))
  Myadj <- as.matrix(numeric(nt))

  # Start maximization of the log-likelihood
  for (iter in seq.int(iter.max)) {
    # Store \eta, \beta, and deviance of the previous iteration
    eta.old <- eta
    dev.old <- dev

    # Compute weights and dependent variable
    mu.eta <- family[["mu.eta"]](eta)
    w <- (wt * mu.eta^2) / family[["variance"]](mu)
    yadj <- (y - mu) / mu.eta + eta - offset

    # Centering dependent variable and compute \eta update
    Myadj <- centerVariables((Myadj + yadj), w, k.list, center.tol)
    eta.upd <- yadj - as.vector(Myadj) + offset - eta

    # Step-halving with three checks
    # 1. finite deviance
    # 2. valid \eta and \mu
    # 3. improvement as in glm2
    rho <- 1.0
    for (inner.iter in seq.int(50L)) {
      eta <- eta.old + rho * eta.upd
      mu <- family[["linkinv"]](eta)
      dev <- sum(family[["dev.resids"]](y, mu, wt))
      dev.crit <- is.finite(dev)
      val.crit <- family[["valideta"]](eta) && family[["validmu"]](mu)
      imp.crit <- (dev - dev.old) / (0.1 + abs(dev)) <= - dev.tol
      if (dev.crit && val.crit && imp.crit) break
      rho <- rho / 2.0
    }

    # Check if step-halving failed
    if (!dev.crit || !val.crit) {
      stop("Inner loop failed; cannot correct step size.", call. = FALSE)
    }

    # Check termination condition
    if (abs(dev - dev.old) / (0.1 + abs(dev)) < dev.tol) break

    # Update starting guesses for acceleration
    Myadj <- Myadj - yadj
  }

  # Return \eta
  eta
}


# Generate auxiliary list of indexes for different sub panels
getIndexList <- function(k.vars, data) {
  indexes <- seq.int(0L, nrow(data) - 1L)
  lapply(k.vars, function(x, indexes, data) {
    split(indexes, data[[x]])
  }, indexes = indexes, data = data)
}


# Compute score matrix
getScoreMatrix <- function(object) {
  # Extract required quantities from result list
  control <- object[["control"]]
  data <- object[["data"]]
  eta <- object[["eta"]]
  wt <- object[["weights"]]
  family <- object[["family"]]
  
  # Update weights and dependent variable
  y <- data[[1L]]
  mu <- family[["linkinv"]](eta)
  mu.eta <- family[["mu.eta"]](eta)
  w <- (wt * mu.eta^2) / family[["variance"]](mu)
  nu <- (y - mu) / mu.eta
  
  # Center regressor matrix (if required)
  if (control[["keep.mx"]]) {
    MX <- object[["MX"]]
  } else {
    # Extract additional required quantities from result list
    formula <- object[["formula"]]
    k.vars <- names(object[["lvls.k"]])
    
    # Generate auxiliary list of indexes to project out the fixed effects
    k.list <- getIndexList(k.vars, data)
    
    # Extract regressor matrix
    X <- model.matrix(formula, data, rhs = 1L)[, - 1L, drop = FALSE]
    nms.sp <- attr(X, "dimnames")[[2L]]
    attr(X, "dimnames") <- NULL
    
    # Center variables
    MX <- centerVariables(X, w, k.list, control[["center.tol"]])
    colnames(MX) <- nms.sp
  }
  
  # Return score matrix
  MX * (nu * w)
}


# Higher-order partial derivatives for 'binomial()'
partialMuEta <- function(eta, family, order) {
  # Safeguard \eta if necessary
  if (family[["link"]] != "logit") {
    eta <- family[["linkfun"]](family[["linkinv"]](eta))
  }
  
  # Second- and third-order derivatives
  f <- family[["mu.eta"]](eta)
  if (order == 2L) {
    # Second-order derivative
    if (family[["link"]] == "logit") {
      f * (1.0 - 2.0 * family[["linkinv"]](eta))
    } else if (family[["link"]] == "probit") {
      - eta * f
    } else if (family[["link"]] == "cloglog") {
      f * (1.0 - exp(eta))
    } else {
      - 2.0 * eta / (1.0 + eta^2) * f
    }
  } else {
    # Third-order derivative
    if (family[["link"]] == "logit") {
      f * ((1.0 - 2.0 * family[["linkinv"]](eta))^2 - 2.0 * f)
    } else if (family[["link"]] == "probit") {
      (eta^2 - 1.0) * f
    } else if (family[["link"]] == "cloglog") {
      f * (1.0 - exp(eta)) * (2.0 - exp(eta)) - f
    } else {
      (6.0 * eta^2 - 2.0) / (1.0 + eta^2)^2 * f
    }
  }
}


# Returns suitable name for a temporary variable
tempVar <- function(data) {
  repeat {
    tmp.var <- paste0(sample(letters, 5L, replace = TRUE), collapse = "")
    if (!(tmp.var %in% colnames(data))) {
      break
    }
  }
  tmp.var
}


# Unload
.onUnload <- function(libpath) {
  library.dynam.unload("alpaca", libpath)
}
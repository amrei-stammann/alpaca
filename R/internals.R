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
feglmFit <- function(beta, eta, y, X, k.list, family, control) {
  # Extract control arguments
  center.tol <- control[["center.tol"]]
  dev.tol <- control[["dev.tol"]]
  epsilon <- min(1.0e-07, dev.tol / 1000.0)
  iter.max <- control[["iter.max"]]
  trace <- control[["trace"]]
  
  # Compute initial quantities for the maximization routine
  mu <- family[["linkinv"]](eta)
  wt <- rep(1.0, length(y))
  dev <- sum(family[["dev.resids"]](y, mu, wt))
  null.dev <- sum(family[["dev.resids"]](y, mean(y), wt))
  z <- as.matrix(numeric(length(y)))
  
  # Start maximization of the log-likelihood
  conv <- FALSE
  for (iter in seq.int(iter.max)) {
    # Compute weights and dependent variable
    mu.eta <- family[["mu.eta"]](eta)
    w.tilde <- sqrt(mu.eta^2 / family[["variance"]](mu))
    nu <- (y - mu) / mu.eta
    
    # Centering variables
    z <- centerVariables((z + nu) * w.tilde, w.tilde, k.list, center.tol)
    X <- centerVariables(X * w.tilde, w.tilde, k.list, center.tol)
    
    # Compute update step and update \eta
    beta.upd <- qr.solve(X, z, epsilon)
    eta.upd <- nu - as.vector(z - X %*% beta.upd) / w.tilde
    
    # Step-halving based on residual deviance
    dev.old <- dev
    rho <- 1.0
    for (inner.iter in seq.int(iter.max)) {
      # Compute residual deviance
      mu <- family[["linkinv"]](eta + rho * eta.upd)
      dev <- sum(family[["dev.resids"]](y, mu, wt))
      
      # Check if deviance is not increasing
      if (is.finite(dev) && (dev - dev.old) / (0.1 + abs(dev)) <= - dev.tol) {
        # Update \eta, \beta, and leave step-halving
        eta <- eta + rho * eta.upd
        beta <- beta + rho * beta.upd
        break
      }
      
      # Update \rho
      rho <- rho / 2.0
    }
    
    # Progress information
    if (trace) {
      cat("Iteration=", iter, "\n")
      cat("Deviance=", format(dev, digits = 5L, nsmall = 2L), "\n")
      cat("Estimates=", format(beta, digits = 3L, nsmall = 2L), "\n")
    }
    
    # Check termination condition
    if (abs(dev - dev.old) / (0.1 + abs(dev)) < dev.tol) {
      if (trace) {
        cat("Convergence\n")
      }
      conv <- TRUE
      break
    }
    
    # Update starting guesses for acceleration
    z <- z / w.tilde - nu
    X <- X / w.tilde
  }
  
  # Compute Score and Hessian
  G <- X * as.vector(z)
  H <- crossprod(X)
  
  # Return result list
  list(coefficients  = as.vector(beta),
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
  control <- object[["control"]]
  data <- object[["data"]]
  family <- object[["family"]]
  formula <- object[["formula"]]
  lvls.k <- object[["lvls.k"]]
  nt <- object[["nobs"]][["nobs"]]
  k.vars <- names(lvls.k)
  k <- length(lvls.k)
  lhs <- names(data)[[1L]]
  
  # Compute starting guess for \eta
  tmp.var <- tempVar(data)
  eta <- numeric(nt)
  for (i in k.vars) {
    data[, (tmp.var) := mean(get(lhs)), by = eval(i)]
    if (family[["family"]] == "binomial") {
      eta <- eta + family[["linkfun"]]((data[[tmp.var]] + 0.5) / 2.0) / k
    } else if (family[["family"]] %in% c("Gamma", "inverse.gaussian")) {
      eta <- eta + family[["linkfun"]](data[[tmp.var]]) / k
    } else {
      eta <- eta + family[["linkfun"]](data[[tmp.var]] + 0.1) / k
    }
    data[, (tmp.var) := NULL]
  }
  
  # Generate auxiliary list of indexes to project out the fixed effects
  indexes <- seq.int(0L, nrow(data) - 1L)
  k.list <- lapply(k.vars, function(x, indexes, data) {
    split(indexes, data[[x]])
  }, indexes = indexes, data = data)
  rm(indexes)
  
  # Extract control arguments
  center.tol <- control[["center.tol"]]
  dev.tol <- control[["dev.tol"]]
  iter.max <- control[["iter.max"]]
  
  # Extract dependent variable and compute initial quantities for the maximization routine
  y <- data[[1L]]
  mu <- family[["linkinv"]](eta)
  wt <- rep(1.0, nt)
  dev <- sum(family[["dev.resids"]](y, mu, wt))
  z <- as.matrix(numeric(length(y)))
  
  # Start maximization of the log-likelihood
  for (iter in seq.int(iter.max)) {
    # Compute weights and dependent variable
    mu.eta <- family[["mu.eta"]](eta)
    w.tilde <- sqrt(mu.eta^2 / family[["variance"]](mu))
    q <- eta - offset + (y - mu) / mu.eta
    
    # Centering dependent variable and compute \eta update
    z <- centerVariables((z + q) * w.tilde, w.tilde, k.list, center.tol)
    eta.upd <- q - as.vector(z / w.tilde) + offset - eta
    
    # Step-halving based on residual deviance as common in glm's
    dev.old <- dev
    rho <- 1.0
    for (inner.iter in seq.int(iter.max)) {
      # Compute residual deviance
      mu <- family[["linkinv"]](eta + rho * eta.upd)
      dev <- sum(family[["dev.resids"]](y, mu, wt))
      
      # Check if deviance is not increasing
      if (is.finite(dev) && (dev - dev.old) / (0.1 + abs(dev)) <= - dev.tol) {
        # Update \eta and leave step-halving
        eta <- eta + rho * eta.upd
        break
      }
      
      # Update \rho
      rho <- rho / 2.0
    }
    
    # Check termination condition
    if (abs(dev - dev.old) / (0.1 + abs(dev)) < dev.tol) {
      break
    }
    
    # Update starting guess for acceleration
    z <- z / w.tilde - q
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


# Higher-order partial derivatives
partialMuEta <- function(eta, family, order) {
  f <- family[["mu.eta"]](eta)
  if (order == 2L) {
    if (family[["link"]] == "logit") {
      f * (1.0 - 2.0 * family[["linkinv"]](eta))
    } else if (family[["link"]] == "probit") {
      - eta * f
    } else {
      f * (1.0 - exp(eta))
    }
  } else {
    if (family[["link"]] == "logit") {
      f * ((1.0 - 2.0 * family[["linkinv"]](eta))^2 - 2.0 * f)
    } else if (family[["link"]] == "probit") {
      (eta^2 - 1.0) * f
    } else {
      f * (1.0 - exp(eta)) * (2.0 - exp(eta)) - f
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
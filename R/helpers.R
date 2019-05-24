### Helper functions (not exported)


# Checks if variable is a factor and transforms if necessary
checkFactor <- function(x) {
  if (is.factor(x)) {
    droplevels(x)
  } else {
    factor(x)
  }
}


# Higher-order partial derivatives
partialMuEta <- function(eta, family, order) {
  f <- family[["mu.eta"]](eta)
  if (order == 2L) {
    if (family[["link"]] == "logit") {
      f * (1.0 - 2.0 * family[["linkinv"]](eta))
    } else {
      - eta * f
    }
  } else {
    if (family[["link"]] == "logit") {
      f * ((1.0 - 2.0 * family[["linkinv"]](eta))^2 - 2.0 * f)
    } else {
      (eta^2 - 1.0) * f
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
.onUnload <- function(libpath) library.dynam.unload("alpaca", libpath)
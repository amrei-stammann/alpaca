### Helper functions (not exported)


# Checks if variable is a factor and transforms if necessary
checkFactor <- function(x) {
  if (is.factor(x)) {
    droplevels(x)
  } else {
    factor(x)
  }
}


# Solve Newton-Rapshon equations using QR factorization
solveNR <- function(X, y) {
  QR <- qr(X)
  as.vector(backsolve(qr.R(QR), crossprod(qr.Q(QR), y)))
}


# Returns suitable name for a temporary variable
tempVar <- function(data) {
  repeat {
    tmpvar <- paste0(sample(letters, 5L, replace = TRUE), collapse = "")
    if (!(tmpvar %in% colnames(data))) break
  }
  tmpvar
}



# Unload
.onUnload <- function(libpath) library.dynam.unload("alpaca", libpath)
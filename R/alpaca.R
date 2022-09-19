#' @title
#' alpaca: A package for fitting glm's with high-dimensional \eqn{k}-way fixed effects
#' @description
#' Provides a routine to partial out factors with many levels during the optimization of the log-likelihood 
#' function of the corresponding generalized linear model (glm). The package is based on the algorithm 
#' described in Stammann (2018) and is restricted to glm's that are based on maximum likelihood estimation 
#' and nonlinear. It also offers an efficient algorithm to recover estimates of the fixed effects in a 
#' post-estimation routine and includes robust and multi-way clustered standard errors. Further the package 
#' provides analytical bias corrections for binary choice models derived by Fernández-Val and Weidner (2016) 
#' and Hinz, Stammann, and Wanner (2020).
#' 
#' \strong{Note:} Linear models are also beyond the scope of this package since there is already a
#' comprehensive procedure available \link[lfe]{felm}.
#' @name
#' alpaca-package
#' @docType
#' package
#' @importFrom
#' data.table setDT setkeyv := .SD
#' @importFrom
#' Formula Formula
#' @importFrom
#' MASS negative.binomial theta.ml
#' @importFrom
#' Rcpp evalCpp
#' @importFrom
#' stats as.formula binomial model.matrix na.omit poisson pnorm printCoefmat rgamma rlogis rnorm 
#' rpois terms vcov
#' @importFrom
#' utils combn
#' @useDynLib
#' alpaca, .registration = TRUE
NULL
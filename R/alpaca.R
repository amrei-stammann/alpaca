#' @title
#' alpaca: A package for fitting glm's with high-dimensional \eqn{k}-way fixed effects.
#' @description
#' Concentrates out factors with many levels during the optimization of the log-likelihood function
#' of the corresponding generalized linear model (glm). The package is restricted to glm's that are
#' based on maximum likelihood estimation. This excludes all quasi-variants of glm. The package also
#' offers an efficient algorithm to recover estimates of the fixed effects in a post-estimation
#' routine. It also includes robust and multi-way clustered standard errors.
#' 
#' \strong{Note:} Linear models are also beyond the scope of this package since there is already a
#' comprehensive procedure available \link[lfe]{felm}.
#' @name
#' alpaca
#' @docType
#' package
#' @importFrom
#' data.table as.data.table := .SD
#' @importFrom
#' Formula Formula
#' @importFrom
#' MASS ginv
#' @importFrom
#' Rcpp evalCpp
#' @importFrom
#' stats binomial coef model.frame model.matrix poisson pnorm predict printCoefmat rgamma rlogis
#' rnorm rpois terms vcov 
#' @importFrom
#' utils combn
#' @useDynLib
#' alpaca, .registration = TRUE
NULL
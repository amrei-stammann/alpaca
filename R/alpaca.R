#' @title
#' alpaca: A package for fitting glm's with high-dimensional \eqn{k}-way fixed
#' effects.
#' @description
#' "ALternating Projection Algorithm for glm's with
#' multiple high-dimensional CAtegorical variables" (alpaca).
#' 
#' \strong{Keywords}: Generalized linear models, alternating projections,  and 
#' high-dimensional \eqn{k}-way fixed effects.
#' @name
#' alpaca
#' @docType
#' package
#' @importFrom
#' Formula Formula model.part
#' @importFrom
#' Rcpp evalCpp
#' @importFrom
#' stats aggregate ave coef lm model.frame model.matrix model.response plogis pnorm predict
#'  printCoefmat rlogis rnorm vcov
#' @useDynLib
#' alpaca, .registration = TRUE
NULL
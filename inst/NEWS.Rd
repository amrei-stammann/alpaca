\name{NEWS}
\title{alpaca news}
\section{Changes in version 0.3.1}{
  \itemize{
    \item All routines now use \code{setDT()} instead of \code{as.data.table()} to avoid unnecessary copies (suggested in #6 @zauster).
    \item \code{feglm.nb()} now returns 'iter' and 'conv' based on iterations of the outer loop.
    \item Fixed a bug in \code{feglm()} that prevented to use \code{I()} for the dependent variable.
    \item Fixed an error in \code{getAPEs()} related to the covariance.
    \item The last line of \code{print.summary.feglm()} now ends with a line break (#6 @zauster).
    \item The internal function \code{feglmFit()} now correctly sets 'conv' if the algorithm does not converge (#5 @zauster).
    \item Fixed some typos in the vignette.
  }
}
\section{Changes in version 0.3}{
  \itemize{
    \item \code{feglm()} now allows to estimate binomial model with fractional response.
    \item Added \code{feglm.nb()} for negative binomial models.
    \item Added post-estimation routine \code{biasCorr()} for analytical bias-corrections (currently restricted to logit and probit models with two-way error component).
    \item Added post-estimation routine \code{getAPEs()} to estimate average partial effects and the corresponding standard errors (currently restricted to logit and probit models with two-way error component).
    \item \code{getFEs()} now returns a list of named vectors. Each vector refers to one fixed effects category (suggested in #4 @zauster).
    \item Changed stopping criteria to the one used by \code{glm()}.
    \item Changed least squares fit to QR (similar to \code{lsfit()} used by \code{glm()}).
    \item Source code and vignettes revised.
  }
}
\section{Changes in version 0.2}{
  \itemize{
    \item Initial release on CRAN.
  }
}
\name{NEWS}
\title{alpaca news}
\section{Changes in version 0.3.3}{
  \itemize{
    \item Stopping condition of \code{feglm.nb()} has been adjusted to better match that of \code{glm.nb()}.
    \item \code{feglm.nb()} now additionally returns 'iter.outer' and 'conv.iter' based on iterations of the outer loop. Previously 'iter' and 'conv' were overwritten.
    \item Step-halving in \code{feglmFit()} and \code{feglmOffset()} is now similar to \code{glm.fit2()}.
    \item Fixed an error in the covariance (influence function) of \code{getAPEs()}.
    \item Updated some references in the documentation and vignette.
    \item Fixed some typos in the documentation and vignette.
  }
}
\section{Changes in version 0.3.2}{
  \itemize{
    \item Added option 'panel.structure' to \code{biasCorr()} and \code{getAPEs()}. This option allows to choose between the two-way bias correction suggested by Fernández-Val and Weidner (2016) and the bias corrections for network data suggested by Hinz, Stammann, and Wanner (2020). Currently both corrections are restricted to probit and logit models.
    \item Added option 'sampling.fe' to \code{getAPEs()} to impose simplifying assumptions when estimating the covariance matrix.
    \item \code{feglm()} now permits to expand functions with \code{poly()} and \code{bs()} (#9 @tcovert).
    \item \code{feglm()} now uses an acceleration scheme suggested by Correia, Guimaraes, and Zylkin (2019) that uses smarter starting values for \code{centerVariables()}.
    \item Added an example of the three-way bias correction suggested by Hinz, Stammann, and Wanner (2020) to the vignette.
    \item The control parameter 'trace' now also returns the current parameter values as well as the residual deviance.
    \item Fixed an error in \code{getAPEs()} related to the estimation of the covariance.
    \item Fixed a bug in the internal function that is used to estimate spectral densities.
  }
}
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
    \item Changed stopping condition to the one used by \code{glm()}.
    \item Changed least squares fit to QR (similar to \code{lsfit()} used by \code{glm()}).
    \item Source code and vignettes revised.
  }
}
\section{Changes in version 0.2}{
  \itemize{
    \item Initial release on CRAN.
  }
}
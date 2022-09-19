# alpaca

## Info
An R-package for fitting glm's with high-dimensional k-way fixed effects.

Provides a routine to partial out factors with many levels during the optimization of the log-likelihood function of the corresponding generalized linear model (glm). The package is based on the algorithm described in [Stammann (2018)](https://arxiv.org/abs/1707.01815) and is restricted to glm's that are based on maximum likelihood estimation and non-linear. It also offers an efficient algorithm to recover estimates of the fixed effects in a post-estimation routine and includes robust and multi-way clustered standard errors. Further the package provides analytical bias corrections for binary choice models (logit and probit) derived by [Fernandez-Val and Weidner (2016)](https://www.sciencedirect.com/science/article/pii/S0304407615002997) and [Hinz, Stammann, and Wanner (2020)](https://arxiv.org/pdf/2004.12655.pdf).

If you have any suggestions for improvements or questions, feel free to [contact me](mailto:amrei.stammann@rub.de).

The package is also available on [CRAN](https://cran.r-project.org/package=alpaca).

## News

### alpaca v0.3.4 (Release Date: 2022-08-10)

Changes:

* Added `vcov.APEs()` generic to extract the covariance matrix after `getAPEs()`.
* Improved the finite sample performance of bias corrections for the average partial effects in case of perfectly classified observations.
* Bias corrections for the average partial effects, i.e. `getAPEs()` after `biasCorr()`, do not require an offset algorithm anymore.
* The default option 'n.pop' in `getAPEs()` has been changed. Now the estimated covariance consists of the delta method part only, i.e. correction factor = 0.
* Improved the numerical stability of the bias corrections.
* `biasCorr()` now also supports one-way fixed effects models.
* Added bias corrections for 'cloglog' and 'cauchit'.
* `feglm()` and `feglm.nb()` do not return a matrix of scores anymore. Instead they, optionally, return the centered regressor matrix. The corresponding option in `feglmControl()` is 'keep.mx'. Default is TRUE.
* Improved the numerical stability of the step-halving in `feglm()`.
* Changed the projection in the MAP algorithm.
* The default option 'center.tol' in `feglmControl()` has been lowered to better handle fitting problems that are not well-behaved.
* Added optional 'weights' argument to `feglm()` and `feglm.nb()`.
* Updated documentation.

### alpaca v0.3.3 (Release Date: 2020-10-30)

Changes:

* Stopping condition of `feglm.nb()` has been adjusted to better match that of `glm.nb()`.
* `feglm.nb()` now additionally returns 'iter.outer' and 'conv.iter' based on iterations of the outer loop. Previously 'iter' and 'conv' were overwritten.
* Step-halving in `feglmFit()` and `feglmOffset()` is now similar to `glm.fit2()`.
* Fixed an error in the covariance (influence function) of `getAPEs()`.
* Updated some references in the documentation and vignette.
* Fixed some typos in the documentation and vignette.

### alpaca v0.3.2 (Release Date: 2020-01-12)

Changes:

* Added option 'panel.structure' to `biasCorr()` and `getAPEs()`. This option allows to choose between the two-way bias correction suggested by Fernandez-Val and Weidner (2016) and the bias corrections for network data suggested by Hinz, Stammann, and Wanner (2020). Currently both corrections are restricted to probit and logit models.
* Added option 'sampling.fe' to `getAPEs()` to impose simplifying assumptions when estimating the covariance matrix.
* `feglm()` now permits to expand functions with `poly()` and `bs()` (#9 @tcovert).
* `feglm()` now uses an acceleration scheme suggested by Correia, Guimaraes, and Zylkin (2019) that uses smarter starting values for `centerVariables()`.
* Added an example of the three-way bias correction suggested by Hinz, Stammann, and Wanner (2019) to the vignette.
* The control parameter 'trace' now also returns the current parameter values as well as the residual deviance.
* Fixed an error in `getAPEs()` related to the estimation of the covariance.
* Fixed a bug in the internal function that is used to estimate spectral densities.

### alpaca v0.3.1 (Release Date: 2019-05-24)

Changes:

* All routines now use `setDT()` instead of `as.data.table()` to avoid unnecessary copies (suggested in #6 @zauster).
* `feglm.nb()` now returns 'iter' and 'conv' based on iterations of the outer loop.
* Fixed a bug in `feglm()` that prevented to use `I()` for the dependent variable.
* Fixed an error in `getAPEs()` related to the covariance.
* The last line of `print.summary.feglm()` now ends with a line break (#6 @zauster).
* The internal function `feglmFit()` now correctly sets 'conv' if the algorithm does not converge (#5 @zauster).

### alpaca v0.3 (Release Date: 2019-05-14)

Changes:

* Added `feglm.nb()` for negative binomial models.
* Added post-estimation routine `biasCorr()` for analytical bias-corrections (currently restricted to logit and probit models with two-way error component).
* Added post-estimation routine `getAPEs()` to estimate average partial effects and the corresponding standard errors (currently restricted to logit and probit models with two-way error component).
* `getFEs()` now returns a list of named vectors. Each vector refers to one fixed effects category.
* Changed stopping criteria to the one used by `glm()`.
* Vignettes revised.

### alpaca v0.2 (Release Date: 2018-07-23)

ATTENTION: Syntax changed slightly. Have a look at the vignettes or help files.

Changes:

* various improvements (glm architecture, clustered standard errors, speed improvements).
* Syntax now more similiar to `glm()`.

### alpaca v0.1.3 (Release Date: 2018-03-08)

Changes:

* added option "cloglog" to argument `family`.
* added checks and routines to ensure that the model response is correctly encoded.

### alpaca v0.1.2 (Release Date: 2018-03-04)

Changes:

* `factor()` should now work as intended.

### alpaca v0.1.1 (Release Date: 2018-01-21)

Changes:

* added option `"probit"` to argument `family`.
* some performance tweaks.
* extract names of the fixed effects and `getFEs()` returns a named vector.
* adjusted computation of starting values.
* computation of the update step (irls) made numerically more stable.


Bugfix:

* construction of the regressor matrix such that factor variables are correctly dummy encoded.
* dropping perfectly classified observations for binomial models should now work as intended (thanks to jmboehm@github).

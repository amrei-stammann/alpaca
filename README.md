# alpaca

## Info
An R-package for fitting glm's with high-dimensional k-way fixed effects.

Provides a routine to concentrate out factors with many levels during the optimization of the log-likelihood function of the corresponding generalized linear model (glm). The package is based on the algorithm proposed by Stammann (2018) [https://arxiv.org/abs/1707.01815] and is restricted to glm's that are based on maximum likelihood estimation and non-linear. It also offers an efficient algorithm to recover estimates of the fixed effects in a post-estimation routine and includes robust and multi-way clustered standard errors. Further the package provides an analytical bias-correction for binary choice models (logit and probit) derived by Fernandez-Val and Weidner (2016).

This package is well suited to estimate so called "pseudo poisson maximum likelihood" (PPML) models with high-dimensional fixed effects that are commonly used in the international trade literature (structural gravity models). See the empirical example in Stammann (2018) and the vignettes.

If you have any suggestions for improvements or questions, feel free to contact me (Amrei.Stammann@hhu.de).

The package is also available on CRAN https://cran.r-project.org/web/packages/alpaca/index.html.

## News

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

# alpaca

## Info
An R-package for fitting glm's with high-dimensional k-way fixed effects.

This is a preliminary R-package based on the working paper "Fast and Feasible Estimation of Generalized Linear Models with High-Dimensional k-way Fixed Effects" (Stammann, 2018), https://arxiv.org/abs/1707.01815. A special Newton-Raphson pseudo-demeaning algorithm is implemented, such that the estimation of glm's with high-dimensional fixed effects becomes feasible. 

This version of the package is able to estimate the entire glm family with high-dimensional k-way fixed effects (excluding quasi families and the linear model).

This package is well suited to estimate so called "pseudo poisson maximum likelihood" (PPML) models with high-dimensional fixed effects that are commonly used in the international trade literature (structural gravity models). See the empirical example in Stammann (2018) and the vignettes.

If you have any suggestions for improvements or questions, feel free to contact me (Amrei.Stammann@hhu.de).

The package can be installed by `devtools::install_github("amrei-stammann/alpaca", build_vignettes = TRUE)`.

## News

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

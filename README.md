# alpaca
An R-package for fitting glm's with high-dimensional k-way fixed effects.

This is a preliminary R-package based on the working paper "Fast and Feasible Estimation of Generalized Linear Models with High-Dimensional k-way Fixed Effects" (Stammann, 2018). A special iteratively re-weighted least squares pseudo-demeaning algorithm is implemented, such that the estimation of glm's with high-dimensional fixed effects becomes feasible.

This version of the package is able to estimate  poisson and logit models with high-dimensional k-way fixed effects. 
In near future the package will be extended to probit models and robust/clustered standard errors.

If you have any suggestions for improvements, feel free to contact me (Amrei.Stammann@hhu.de).

The package can be installed by `devtools::install_github("amrei-stammann/alpaca", build_vignettes = TRUE)`.

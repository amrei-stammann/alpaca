#' @importFrom generics augment
#' @export
#' @seealso [augment.lm()]
generics::augment

#' @importFrom generics tidy
#' @export
#' @seealso [tidy.lm()]
generics::tidy

#' @importFrom generics glance
#' @export
#' @seealso [glance.lm()]
generics::glance

#' @title Tidy a(n) feglm object
#' @description Tidy summarizes information about the components
#'   of a model. A model component might be a single term in a
#'   regression, a single hypothesis, a cluster, or a class.
#'   Exactly what tidy considers to be a model component varies
#'   across models but is usually self-evident. If a model has
#'   several distinct types of components, you will need to
#'   specify which components to return.
#'
#' @param x A `feglm` object returned from [alpaca::feglm()].
#' @param fe Logical indicating whether or not to include
#'   estimates of fixed effects. Defaults to `FALSE`.
#' @param type Character indicating the type of covariance
#'   estimate to be used for the standard errors. Possible values
#'   are `hessian`, `outer.product`, `sandwich`, `cluster`. See
#'   alpaca::summary.feglm for details. Defaults to `hessian`.
#' @param cluster a symbolic description indicating the
#'   clustering of observations.
#' @param conf.int Logical indicating whether or not to include a
#'   confidence interval in the tidied output. Defaults to
#'   `FALSE`. NOT IMPLEMENTED YET.
#' @param conf.level The confidence level to use for the
#'   confidence interval if `conf.int = TRUE`. Must be strictly
#'   greater than 0 and less than 1. Defaults to 0.95, which
#'   corresponds to a 95 percent confidence interval. NOT
#'   IMPLEMENTED YET.
#' @param ... Additional arguments. Not used. Needed to match generic
#'   signature only. **Cautionary note:** Misspelled arguments will be
#'   absorbed in `...`, where they will be ignored. If the misspelled
#'   argument has a default value, the default value will be used.
#'   For example, if you pass `conf.lvel = 0.9`, all computation will
#'   proceed using `conf.level = 0.95`. Additionally, if you pass
#'   `newdata = my_tibble` to an [augment()] method that does not
#'   accept a `newdata` argument, it will use the default value for
#'   the `data` argument.
#'
#' @return A [tibble::tibble()] with one row for each term in the
#'   regression. The tibble has columns:
#'
#'   \item{term}{The name of the regression term.}
#'   \item{estimate}{The estimated value of the regression term.}
#'   \item{std.error}{The standard error of the regression term.}
#'   \item{statistic}{The value of a statistic, almost always a
#'     T-statistic, to use in a hypothesis that the regression
#'     term is non-zero.}
#'   \item{p.value}{The two-sided p-value associated with the
#'     observed statistic.}
#'   \item{conf.low}{The low end of a confidence interval for the
#'     regression term. Included only if `conf.int = TRUE`.}
#'   \item{conf.high}{The high end of a confidence interval for
#'     the regression term. Included only if `conf.int = TRUE`.}
#'
#' @examples
#'
#' library(alpaca)
#'
#' dt <- simGLM(n = 40, t = 10, seed = 123, model = "poisson")
#'
#' res_feglm <- feglm(formula = y ~ x1 + x2 + x3 | i + t,
#'                    data = dt, family = poisson())
#' res_feglm
#' tidy(res_feglm)
#' tidy(res_feglm, type = "sandwich")
#' tidy(res_feglm, type = "clustered", cluster = ~ i + t)
#' augment(res_feglm)
#' glance(res_feglm)
#'
#' @aliases feglm_tidiers alpaca_tidiers
#' @family feglm tidiers
#' @seealso [alpaca::feglm()]
#' @export
tidy.feglm <- function(x,
                       fe = FALSE, type = "hessian",
                       cluster = NULL,
                       conf.int = NULL, conf.level = NULL,
                       ...) {
  if(!is.null(conf.int) & !is.null(conf.level))
    warning("Both conf.int and conf.level are not implemented yet. Ignoring now.")
  if(type != "clustered" & !is.null(cluster)) {
    warning("Specify clustered SE with type = 'clustered'. Ignoring the cluster expression now.")
    cluster <- NULL
  }

  xs <- summary(x, type = type, cluster = cluster)

  ## get the coefficient matrix and add the variable names as the
  ## first column
  ret <- data.frame(xs$cm)
  colnames(ret) <- c("estimate", "std.error", "statistic", "p.value")
  ## ret <- tibble::add_column(ret, term = rownames(xs$cm), .before = 1)
  ret$term <- rownames(xs$cm)
  ret <- ret[, c("term", "estimate", "std.error",
                 "statistic", "p.value")]

  ## tibble::as_tibble(ret)
  ret
}

#' @title Augment data with information from a(n) feglm
#'   object
#'
#' @description Augment accepts a model object and a dataset and
#'   adds information about each observation in the dataset. Most
#'   commonly, this includes predicted values in the `.fitted`
#'   column, residuals in the `.resid` column, and standard
#'   errors for the fitted values in a `.se.fit` column. New
#'   columns always begin with a `.` prefix to avoid overwriting
#'   columns in the original dataset.
#'
#'   Users may pass data to augment via either the `data`
#'   argument or the `newdata` argument. If the user passes data
#'   to the `data` argument, it **must** be exactly the data that
#'   was used to fit the model object. Pass datasets to `newdata`
#'   to augment data that was not used during model fitting. This
#'   still requires that all columns used to fit the model are
#'   present.
#'
#'   Augment will often behavior different depending on whether
#'   `data` or `newdata` is specified. This is because there is
#'   often information associated with training observations
#'   (such as influences or related) measures that is not
#'   meaningfully defined for new observations.
#'
#'   For convenience, many augment methods provide default `data`
#'   arguments, so that `augment(fit)` will return the augmented
#'   training data. In these cases augment tries to reconstruct
#'   the original data based on the model object, with some
#'   varying degrees of success.
#'
#'   The augmented dataset is always returned as a
#'   [tibble::tibble] with the **same number of rows** as the
#'   passed dataset. This means that the passed data must be
#'   coercible to a tibble. At this time, tibbles do not support
#'   matrix-columns. This means you should not specify a matrix
#'   of covariates in a model formula during the original model
#'   fitting process, and that [splines::ns()], [stats::poly()]
#'   and [survival::Surv()] objects are not supported in input
#'   data. If you encounter errors, try explicitly passing a
#'   tibble, or fitting the original model on data in a tibble.
#'
#'   We are in the process of defining behaviors for models fit
#'   with various [na.action] arguments, but make no guarantees
#'   about behavior when data is missing at this time.
#'
#' @inherit tidy.feglm params examples
#' @param data A [data.frame()] or [tibble::tibble()] containing
#'   the original data that was used to produce the object `x`.
#'   Defaults to `stats::model.frame(x)` so that
#'   `augment(my_fit)` returns the augmented original data. **Do
#'   not** pass new data to the `data` argument. Augment will
#'   report information such as influence and cooks distance for
#'   data passed to the `data` argument. These measures are only
#'   defined for the original training data.
#'
#' @return A [tibble::tibble()] containing the data passed to `augment`,
#'   and **additional** columns:
#'
#'   \item{.fitted}{The predicted response for that observation.}
#'   \item{.resid}{The residual for a particular point. Present
#'   only when data has been passed to `augment` via the `data`
#'   argument.}
#'
#' @family feglm tidiers
#' @seealso [alpaca::feglm()]
#' @importFrom stats fitted
#' @export
augment.feglm <- function(x, data = x$data, ...) {
  ## df <- tibble::as_tibble(x$data)
  df <- x$data
  fitted_vals <- as.vector(fitted(x))
  resid_vals <- as.vector(data[[1]]) - fitted_vals
  df$`.fitted` <- fitted_vals
  df$`.resid` <- resid_vals

  df
}

#' @title Glance at a(n) feglm object
#'
#' @description Glance accepts a model object and returns a
#'   [tibble::tibble()] with exactly one row of model summaries.
#'   The summaries are typically goodness of fit measures,
#'   p-values for hypothesis tests on residuals, or model
#'   convergence information.
#'
#'   Glance never returns information from the original call to
#'   the modelling function. This includes the name of the
#'   modelling function or any arguments passed to the modelling
#'   function.
#'
#'   Glance does not calculate summary measures. Rather, it farms
#'   out these computations to appropriate methods and gathers
#'   the results together. Sometimes a goodness of fit measure
#'   will be undefined. In these cases the measure will be
#'   reported as `NA`.
#'
#' @inherit tidy.feglm params examples
#'
#' @return A one-row [tibble::tibble] with columns:
#'
#'   \item{null.deviance}{Null deviance of the model.}
#'   \item{deviance}{Deviance of the model.}
#'   \item{nobs}{Number of observations of the model.}
#'
#' @export
glance.feglm <- function(x, ...) {
  ## TODO add AIC, BIC, logLik, ...
  ## ret <- tibble::tibble(
  ret <- data.frame(
                   null.deviance = x[["null.deviance"]],
                   deviance = x[["deviance"]],
                   nobs = x[["nobs"]]["nobs"]
                 )
  ret
}

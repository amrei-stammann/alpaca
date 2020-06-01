#' @templateVar class feglm
#' @template title_desc_tidy
#'
#' @param x A `feglm` object returned from [alpaca::feglm()].
#' @template param_confint
#' @param fe Logical indicating whether or not to include
#'   estimates of fixed effects. Defaults to `FALSE`.
#' @param type Character indicating the type of covariance
#'   estimate to be used for the standard errors. Possible values
#'   are `hessian`, `outer.product`, `sandwich`, `cluster`. See
#'   alpaca::summary.feglm for details. Defaults to `hessian`.
#' @param cluster a symbolic description indicating the
#'   clustering of observations.
#' @template param_unused_dots
#'
#' @evalRd return_tidy(regression = TRUE)
#'
#' @examples
#'
#' library(alpaca)
#'
#' dt <- simGLM(n = 40, t = 10, seed = 123, model = "poisson")
#' setDT(dt)
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
#' @export
#' @aliases feglm_tidiers alpaca_tidiers
#' @family feglm tidiers
#' @seealso [tidy()], [alpaca::feglm()]
tidy.feglm <- function(x, ## conf.int = FALSE, conf.level = .95,
                       fe = FALSE, type = "hessian",
                       cluster = NULL,
                       ...) {
  if(type != "clustered" & !is.null(cluster)) {
    warning("Specify clustered SE with type = 'clustered'. Ignoring the cluster expression now.")
    cluster <- NULL
  }

  xs <- summary(x, type = type, cluster = cluster)

  ## get the coefficient matrix and add the variable names as the
  ## first column
  ret <- data.frame(xs$cm)
  colnames(ret) <- c("estimate", "std.error", "statistic", "p.value")
  ret <- tibble::add_column(ret, term = rownames(xs$cm), .before = 1)

  tibble::as_tibble(ret)
}

#' @templateVar class feglm
#' @template title_desc_augment
#'
#' @inherit tidy.feglm params examples
#' @template param_data
#'
#' @evalRd return_augment()
#'
#' @export
#' @family feglm tidiers
#' @seealso [augment()], [alpaca::feglm()]
augment.feglm <- function(x, data = x$data, ...) {
  df <- tibble::as_tibble(x$data)
  fitted_vals <- as.vector(fitted(x))
  resid_vals <- as.vector(data[[1]]) - fitted_vals
  dplyr::mutate(df,
                .fitted = fitted_vals,
                .resid = resid_vals)
}

#' @templateVar class feglm
#' @template title_desc_glance
#'
#' @inherit tidy.feglm params examples
#'
#' @evalRd return_glance(
#'   "null.deviance",
#'   "deviance",
#'   "nobs"
#' )
#'
#' @export
glance.feglm <- function(x, ...) {
  ## TODO add AIC, BIC, logLik, ...
  ret <- tibble::tibble(
                   null.deviance = x[["null.deviance"]],
                   deviance = x[["deviance"]],
                   nobs = x[["nobs"]]["nobs"]
                 )
  ret
}

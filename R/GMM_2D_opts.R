#' Default configuration for GMM_2D
#'
#' A list with parameters customizing a GMM_2D. Each component of the
#' list is an effective argument for **nazwa funkcji głównej**.
#'
#' @param eps_change Criterion for stopping the EM algorithm. By default it is \code{1e-4}
#' @param max_iter Maximum number of iterations of EM algorithm. By default it is \code{max_iter = 5000}
#' @param SW Regularizing coefficient for covariance (nie ja)
#' @param max_var_ratio Maximum disimilarity between horizontal and vertical dispersion. By default it is \code{max_var_ratio = 5}
#' @param crit Information Criterion to select best number of components.
#' Possible "AIC","AICc", "BIC" (default) or "ICL-BIC".
#' @param res Parameter for modyifing resolution (>1 less spots). **1** is as a default value.
#' @param cov_type Type of covariance model. Possible "sphere","diag" or "full" (default).
#' @param show \code{True} if the figure should be displayed.
#' @param init_nb Number of random initial conditions. By default it is \code{init_nb = 10}.
#' @param KS Maximum number of GMM components. By default it is \code{KS = 5}.
#' @param D_thr Significance threshold for D statistic By default is 0.1.(Likelihood Ratio Test ?)
#' @param init_con Type of initial conditions. Could be "rand" (default) or "diag".
#'
#' @examples
#' # display all default settings
#' GMM_2D_opts
#'
#' # create a new settings object
#' custom.settings <- GMM_2D_opts
#' custom.settings$crit <- "BIC"
#' custom.settings
#'
#' @export
GMM_2D_opts <- list(
  eps_change = 1e-4,
  max_iter = 5000,
  SW = 10,            #regularizing coefficient for covariance
  max_var_ratio = 5,
  crit = 'BIC',
  res = 1,
  cov_type = 'full',
  show = TRUE,
  init_nb = 10,
  KS = 5,
  D_thr = 0.1,
  init_con = "rand"
)
class(GMM_2D_opts) <- "gmm2_opts"



GMM_1D_opts <- list(
  bla = 4,
  bla1 =34.5
)
class(GMM_1D_opts) <- "gmm2_opts"
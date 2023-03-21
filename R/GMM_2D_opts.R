#' Default configuration for 2D Gaussian Mixture decomposition
#'
#' A list with parameters customizing a GMM_2D. Each component of the
#' list is an effective argument for \code{\link{runGMM2D}}.
#'
#' @param eps_change Criterion for early stopping of EM (1e-7, by default).
#' @param max_iter Maximum number of iterations of EM algorithm. By default it is \code{max_iter = 50 000}.
#' @param SW Regularizing coefficient for covariance.
#' @param max_var_ratio Maximum dissimilarity between horizontal and vertical dispersion. By default it is \code{max_var_ratio = 5}
#' @param IC Information criterion used to select the number of model components.
#' Possible methods are "AIC","AICc", "BIC" (default), "ICL-BIC" or "LR".
#' @param cov_type Type of covariance defined for each model component. Possible "sphere","diag" or "full" (default).
#' @param init_nb Number of random initial conditions. By default it is \code{init_nb = 10}.
#' @param KS Maximum number of components of the model. By default it is \code{KS = 5}.
#' @param quick_stop Logical value. Determines if stop searching of the number of components earlier based on the Likelihood Ratio Test. Used to speed up the function (TRUE, by default).
#' @param signi Significance level set for Likelihood Ratio Test (0.05, by default).
#' @param init_con Type of initial conditions. Could be "rand" (default),"DP" or "diag".
#' @param fixed Logical value. Fit GMM for selected number of components given by KS (FALSE, by default).

#'
#' @examples
#' # display all default settings
#' GMM_2D_opts
#'
#' # create a new settings object
#' custom.settings <- GMM_2D_opts
#' custom.settings$IC <- "AIC"
#' custom.settings
#'
#' @export
GMM_2D_opts <- list(
  eps_change = 1e-7,
  max_iter = 50000,
  SW = 0.01,
  max_var_ratio = 5,
  IC = 'BIC',
  cov_type = 'full',
  init_nb = 10,
  KS = 5,
  quick_stop=FALSE,
  signi = 0.05,
  init_con = "rand",
  fixed=FALSE
)
class(GMM_2D_opts) <- "gmm2_opts"

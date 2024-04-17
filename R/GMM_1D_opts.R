#' Default configuration for 1D Gaussian Mixture decomposition
#'
#' A list with parameters customizing a GMM for 1D and binned data. Each component of the
#' list is an effective argument for \code{\link{runGMM}}.
#'
#' @param KS Maximum number of components of the model.
#' @param eps_change Criterion for early stopping of EM (1e-7, by default) given by the following formula:
#' \deqn{\sum{(|\alpha - \alpha_{old})|} + \frac{\sum{(\frac{|\sigma^2 - \sigma^2_{old}|}{\sigma^2})}}{length(\alpha)}}
#' @param max_iter Maximum number of iterations of EM algorithm. By default it is \code{max_iter = 10 000}.
#' @param SW Parameter for calculating minimum variance of each Gaussian component (0.25, by default) using the following formula:
#' \deqn{(\frac{SW*range(x)}{no.of.components)})^2}. Lower value means smaller component variance allowed.
#' @param IC Information criterion used to select the number of model components.
#' Possible methods are "AIC","AICc", "BIC" (default), "ICL-BIC" or "LR".
#' @param sigmas.dev Parameter used to define close GMM components that needs to be merged. For each component, standard deviation is multiplied by \code{sigmas.dev} to estimate the distance from component mean.
#' All other components within this distance are merged. By default it is \code{sigmas.dev = 1}. When \code{sigmas.dev = 0} no components are merged.
#' @param quick_stop Logical value. Determines if stop searching of the number of components earlier based on the Likelihood Ratio Test. Used to speed up the function (TRUE, by default).
#' @param signi Significance level set for Likelihood Ratio Test (0.05, by default).
#' @param fixed Logical value. Fit GMM for selected number of components given by KS (FALSE, by default).
#' @param plot Logical value. If TRUE (default), the figure visualizing GMM decomposition will be displayed.
#' @param col.pal Name of the RColorBrewer palette used in the figure. By default \code{"Blues"}.
#'
#' @examples
#' # display all default settings
#' GMM_1D_opts
#'
#' # create a new settings object
#' custom.settings <- GMM_1D_opts
#' custom.settings$IC <- "AIC"
#' custom.settings
#'
#' @export
GMM_1D_opts <- list(
  KS = 15,
  eps_change = 1e-7,
  max_iter = 10000,
  SW = 0.25,
  IC = 'BIC',
  sigmas.dev = 1,
  quick_stop = FALSE,
  signi = 0.05,
  fixed = FALSE,
  plot = TRUE,
  col.pal = "Blues"
)
class(GMM_1D_opts) <- "gmm2_opts"

#' Default configuration for GMM_2D
#'
#' A list with parameters customizing a GMM_2D embedding. Each component of the
#' list is an effective argument for **nazwa funkcji głównej**.
#'
#' argument1: bla bla bla
#'
#' argument2: tez bla bla bla
#'
#' @examples
#' # display all default settings
#' GMM_2D_opts
#'
#' # create a new settings object
#' custom.settings <- GMM_2D_opts
#' custom.settings$IC <- "BIC"
#' custom.settings
#'
#' @export
GMM_2D_opts <- list(

)
class(GMM2D_opts) <- "gmm2_opts"

#' Function to fit Gaussian Mixture Model (GMM) to 2D data
#'
#' Main function to perform GMM on 2D data. Function choose the optimal number of components of a 2D mixture normal distributions by minimizing the value of the information criterion.
#'
#' @param X Matrix of 2D data to decompose by GMM.
#' @param Y Vector of counts, with the same length as "X".
#' Applies only to binned data (Y = NULL, by default).
#' @param opts Parameters of run stored in \code{\link{GMM_2D_opts}} variable.
#'
#' @returns Function returns a \code{list} of GMM parameters for the estimated number of components: \describe{
#'  \item{model}{\describe{
#'  \item{alpha}{Weights (alpha) of each component.}
#'  \item{center}{Means of decomposition.}
#'  \item{covar}{Covariances of each component.}
#'  \item{KS}{Estimated number of components.}
#'  \item{logL}{Log-likelihood statistic for the estimated number of components.}
#'  \item{IC}{The value of the selected information criterion which was used to calculate the number of components.}
#'  \item{cls}{Assigment of point to the clusters.}}}
#'  \item{fig}{Plot of decomposition.}
#' }
#'
#' @examples
#' \dontrun{
#' data(example2D)
#' custom.settings <- GMM_2D_opts
#' custom.settings$fixed <- TRUE
#' custom.settings$KS <- 3
#' custom.settings$max_iter <- 5000
#' custom.settings$plot <- TRUE
#'
#' res <- runGMM2D(example2D[,1:2], example2D[,3], opts = custom.settings)
#' }
#'
#'
#' @export
runGMM2D <- function(X, Y = NULL, opts = NULL){

  if(is.null(opts)){opts = rGMMtest::GMM_2D_opts}

  # Check part
  if (!hasArg("X")){
    stop("No data.")}

  if (length(X) < 2){
    stop("Not enough data.")}


  # Main GMM run
  GModel <- gaussian_mixture_2D(X, Y, opts)

  # Plot part
  if (!is.null(Y)){
    pl <- plot_gmm_2D_binned(X, Y, GModel, opts)
  } else{
    pl <- plot_gmm_2D_orig(X, GModel, opts)
  }


  # Return
  mix_gmm <- list(model = GModel, fig = pl)


  # Print the plot
  if(opts$plot){
    print(pl)}

  return(mix_gmm)

}

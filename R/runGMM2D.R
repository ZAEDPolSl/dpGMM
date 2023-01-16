#' Function to perform Gaussian Mixture Model on 2D data
#'
#' Function to choose the optimal number of components of a 2D mixture normal distributions, minimising the value of the information criterion.
#'
#' @param X matrix of data to decompose by GMM.
#' @param Y Vector of counts, should be the same length as "data".
#' Applies only to binned data therefore the default is Y = NULL.
#' @param opts parameters of run saves in \code{\link{GMM_2D_opts}} variable
#' @param plot Logical value. If TRUE (default FALSE), the GMM figure will be displayed.
#'
#' @returns Function returns a \code{list} of GMM parameters for the optimal number of components: \describe{
#'  \item{alpha}{Weights (alpha) of each component}
#'  \item{center}{Means of decomposition}
#'  \item{covar}{Covariances of each component}
#'  \item{KS}{Optimal number of components}
#'  \item{logL}{Log-likelihood value for the optimal number of components}
#'  \item{IC}{Value of the selected information criterion which was used to calculate the optimal number of components}
#'  \item{fig}{Plot of decomposition}
#' }
#'
#' @examples
#' \dontrun{
#' data(example2D_1)
#' opts<-GMM_2D_opts
#' runGMM2D(example2D_1[,1:2], example2D_1[,3], opts, plot=T)
#' }
#'
#'
#' @export
runGMM2D <- function(X, Y=NULL, opts, plot=FALSE){
  # Check part
  if (!hasArg("X")){
    stop("No data.")}

  if (length(X) < 2){
    stop("Not enough data.")}


  # Main GMM run
  GModel<- gaussian_mixture_2D(X, Y, opts)

  # Plot part
  if (!is.null(Y)){
    pl<-plot_gmm_2D_binned(X, Y, GModel, opts)
  } else{
    pl<-plot_gmm_2D_orig(X, GModel, opts)
  }


  # return
  mix_gmm <- list(model = GModel, fig=pl)


  # Print the plot
  if(plot){
    print(pl)}

  return(mix_gmm)

}

#' Function to fit Gaussian Mixture Model (GMM) to 1D data
#'
#' Function fits GMM with initial conditions found using dynamic programming-based approach by using expectation-maximization (EM) algorithm.
#' The function works on original and binned (e.g. obtained by creating histogram on 1D data) data. Additionally, threshold values that allows to assign data to individual Gaussian components are provided.
#' Function allows to estimate the number of GMM components using five different information criteria and merging of similar components.
#'
#'
#' @param X Vector of 1D data for GMM decomposition.
#' @param KS Maximum number of components of the model.
#' @param opts Parameters of run saved in \code{\link{GMM_1D_opts}} variable.
#' @param Y Vector of counts, with the same length as "X".
#' Applies only to binned data (Y = NULL, by default).
#' @param plot Logical value. If TRUE (default), the figure visualizing GMM decomposition will be displayed.
#' @param col.pal Name of the RColorBrewer palette used in the figure. By default \code{"Blues"}.
#'
#' @returns Function returns a \code{list} which contains: \describe{
#'  \item{model}{A \code{list} of model component parameters - mean values (mu), standard deviations (sigma)
#'  and weights (alpha) for each component. Output of \code{\link{gaussian_mixture_vector}}.}
#'  \item{KS}{Estimaged number of model components.}
#'  \item{IC}{The value of the selected information criterion which was used to calculate the number of components.}
#'  \item{logLik}{Log-likelihood statistic for the estimated number of components.}
#'  \item{threshold}{Vector of thresholds between each component.}
#'  \item{cluster}{Assignment of original \code{X} values to individual components (clusters) by thresholds.}
#'  \item{fig}{ggplot object (output of the \code{\link{plot_gmm_1D}} function). It contains GMM decomposition together with a histogram of the data.}
#'  \item{QQplot}{ggplot object (output of the \code{\link{plot_QQplot}} function).
#'  It presents diagnostic Quantile-Quantile plot for a single normal distribution and fitted GMM.}
#' }
#'
#' @importFrom methods hasArg
#' @examples
#' \dontrun{
#' data(example)
#'
#' custom.settings <- GMM_1D_opts
#' custom.settings$sigmas.dev <- 1.5
#' custom.settings$max_iter <- 1000
#'
#' mix_test <- runGMM(example$Dist, KS = 10, opts = custom.settings)
#' mix_test$QQplot
#'
#' #example for binned data
#' data(binned)
#' custom.settings <- GMM_1D_opts
#' custom.settings$quick_stop <- TRUE
#'
#' binned_test <- runGMM(X = binned$V1, KS = 40, opts = custom.settings, Y = binned$V2, col.pal = "Dark2", plot = F)
#' binned_test$fig
#' }
#'
#' @seealso \code{\link{gaussian_mixture_vector}}, \code{\link{EM_iter}}
#'
#' @export
# runGMM <- function(X, KS, Y = NULL, fixed = FALSE , eps_change = 1e-7, max_iter = 50000, SW = 0.01, IC = "BIC", sigmas.dev = 2.5,
                   # plot = TRUE, col.pal = "Blues", quick_stop = TRUE, signi = 0.05) {
runGMM <- function(X, KS, opts = GMM_1D_opts, Y = NULL, plot = TRUE, col.pal = "Blues"){
  # Check part
  if (!hasArg("X")){
    stop("No data.")}

  if (length(X) < 2){
    stop("Not enough data.")}

  if (KS < 2){
    stop("KS (no of components) must be larger than 1.")}

  IC_list <- c("AIC","AICc", "BIC", "ICL-BIC", "LR")
  if (!opts$IC %in% IC_list) {
    stop("Criterion not implemented. Please use AIC, AICc, BIC, ICL-BIC or LR")
  }


  # Main GMM run
  # GModel <- gaussian_mixture_vector(X, KS, Y, fixed , eps_change, max_iter, SW, IC, quick_stop, signi)
  GModel <- gaussian_mixture_vector(X, KS, opts, Y)

  # GMM components merging
  if(opts$sigmas.dev > 0 & GModel$KS > 1){
    IC_tmp <- GModel$IC
    logL_tmp <- GModel$logL
    GModel <- rGMMtest:::gmm_merge(GModel$model, opts$sigmas.dev) ######################## corect to package name
    GModel$IC <- IC_tmp
    GModel$logL <- logL_tmp
  }

  # Generating distribution from model
  dist.plot <- generate_dist(X, GModel$model, 1e4)

  # Thresholds estimation
  if(GModel$KS > 1){
    thr <- find_thr_by_params(GModel$model, dist.plot, opts$sigmas.dev)
  } else {thr = NULL}

  # remove thresholds out of data range
  rem <- which(thr > max(X) | thr < min(X))
  if (length(rem) != 0){thr <- thr[-rem]}

  # Clusters assignment
  clust <- matrix(1, 1, length(X))
  for(i in 1:length(thr)){clust[X > thr[i]] <- i+1}

  # Plot generating
  pl <- plot_gmm_1D(X, dist.plot, Y, thr, pal = col.pal)

  # QQplot
  pl.qq <- plot_QQplot(X, GModel$model)

  # Output of function
  mix_gmm <- list(model = GModel$model, KS = nrow(GModel$model), IC = GModel$IC, logLik = GModel$logL,
                  threshold = thr, cluster = as.vector(clust), fig = pl, QQplot = pl.qq)
  names(mix_gmm)[3] <- opts$IC

  # Print the plot
  if(plot){
    p <- mix_gmm[["fig"]]
    print(p)
    }

  return(mix_gmm)
}

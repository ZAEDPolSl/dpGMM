#' Function to fit Gaussian Mixture Model (GMM) to 1D data
#'
#' Function fits GMM with initial conditions found using dynamic programming-based approach by using expectation-maximization (EM) algorithm.
#' The package works on original and binned (e.g. obtained by creating histogram on 1D data) data. Additionally, threshold values that allows to assign data to individual Gaussian components are provided.
#' Function allows to estimate the number of GMM components using five different information criteria and merging of similar components.
#'
#'
#' @param X Vector of 1D data for GMM decomposition.
#' @param KS Maximum number of components of the model.
#' @param Y Vector of counts, with the same length as "X".
#' Applies only to binned data (Y = NULL, by default).
#' @param fixed Fit GMM for selected number of components given by KS.
#' @param eps_change Criterion for early stopping of EM (1e-7, by default) given by the following formula:
#' \deqn{\sum{(|\alpha - \alpha_{old})|} + \frac{\sum{(\frac{|\sigma^2 - \sigma^2_{old}|}{\sigma^2})}}{length(\alpha)}}
#' @param max_iter Maximum number of iterations of EM algorithm. By default it is \code{max_iter = 50 000}
#' @param SW Parameter for calculating minimum variance of each Gaussian component (0.01, by default) using the following formula:
#' \deqn{\frac{SW*range(x)}{no.of.components)^2}}. Lower value means smaller component variance allowed.
#' @param IC Information criterion used to select the number of model components.
#' Possible methods are "AIC","AICc", "BIC" (default), "ICL-BIC" or "LR".
#' @param sigmas.dev Parameter used to define close GMM components that needs to be merged. For each component, standard deviation is multiplied by \code{sigmas.dev} to estimate the distance from component mean.
#' All other components within this distance are merged. By default it is \code{sigmas.dev = 2.5}. When \code{sigmas.dev = 0} no components are merged.
#' @param plot Logical value. If TRUE (default), the figure visualizing GMM decomposition will be displayed.
#' @param col.pal Name of the RColorBrewer palette used in the figure. By default \code{"Blues"}.
#' @param quick_stop Logical value. Determines if stop searching of the number of components earlier based on the Likelihood Ratio Test. Used to speed up the function (TRUE, by default).
#' @param signi Significance level set for Likelihood Ratio Test (0.05, by default).
#'
#' @returns Function returns a \code{list} which contains: \describe{
#'  \item{model}{A \code{list} of model component parameters - mean values (mu), standard deviations (sigma)
#'  and weights (alpha) for each component. Output of \code{gaussian_mixture_vector}.}
#'  \item{KS}{Estimaged number of model components.}
#'  \item{IC}{The value of the selected information criterion which was used to calculate the number of components.}
#'  \item{logLik}{Log-likelihood statistic for the estimated number of components.}
#'  \item{threshold}{Vector of thresholds between each component.}
#'  \item{cluster}{Assignment of original \code{X} values to individual components (clusters) by thresholds.}
#'  \item{fig}{ggplot object (output of the \code{plot_gmm_1D} function). It contains GMM decomposition together with a histogram of the data.}
#'  \item{QQplot}{ggplot object (output of the \code{plot_QQplot} function).
#'  It presents diagnostic Quantile-Quantile plot for a single normal distribution and fitted GMM.}
#' }
#'
#' @importFrom methods hasArg
#' @examples
#' \dontrun{
#' data(example)
#' mix_test1 <- runGMM(example$Dist, KS = 15, IC = "AICc", quick_stop = F, sigmas.dev = 0)
#' mix_test2 <- runGMM(example$Dist, KS = 10, IC = "BIC", sigmas.dev = 1.5)
#' mix_test2$QQplot
#'
#' data(binned)
#' binned_test <- runGMM(X = binned$V1, Y = binned$V2, KS = 40, col.pal = "Dark2", plot = F, quick_stop = T)
#' binned_test$fig
#' }
#'
#' @seealso \code{\link{gaussian_mixture_vector}}, \code{\link{EM_iter}}, \code{\link{generate_dist}}, \code{\link{find_thr_by_params}}
#'
#' @export
runGMM <- function(X, KS, Y = NULL, fixed=FALSE , eps_change = 1e-7, max_iter = 50000, SW=0.01, IC = "BIC", sigmas.dev = 2.5,
                   plot = TRUE, col.pal="Blues", quick_stop = TRUE, signi = 0.05) {
  # Check part
  if (!hasArg("X")){
    stop("No data.")}

  if (length(X) < 2){
    stop("Not enough data.")}

  if (KS < 2){
    stop("KS (no of components) must be larger than 1.")}

  IC_list <- c("AIC","AICc", "BIC", "ICL-BIC", "LR")
  if (!IC %in% IC_list) {
    stop("Criterion not implemented. Please use AIC, AICc, BIC, ICL-BIC or LR")
  }


  # Main GMM run
  GModel<- gaussian_mixture_vector(X, KS, Y, fixed , eps_change, max_iter, SW, IC, quick_stop, signi)

  # GMM components merging
  if(sigmas.dev > 0 & GModel$KS>1){
    IC_tmp <- GModel$IC
    logL_tmp <- GModel$logL
    GModel <- gmm_merge(GModel$model, sigmas.dev)
    GModel$IC <- IC_tmp
    GModel$logL <- logL_tmp
  }

  # Generating distribution from model
  dist.plot<-generate_dist(X, GModel$model, 1e4)

  # Thresholds estimation
  #thr <- find_thr_by_dist_VMM(dist.plot)
  if(GModel$KS>1){
    thr <- find_thr_by_params(GModel$model,dist.plot,sigmas.dev)
  } else {thr=NULL}
  # remove Thresholds out of data range
  thr<-thr[-which(thr>max(X) | thr<min(X))]

  # Clusters assignment
  clust <- matrix(1, 1, length(X))
  for(i in 1:length(thr)){clust[X>thr[i]] <- i+1}

  # Plot generating
  pl<-plot_gmm_1D(X, dist.plot, Y, thr, pal=col.pal)

  # QQplot
  pl.qq<-plot_QQplot(X,GModel$model)

  # Output of function
  mix_gmm <- list(model = GModel$model, KS = nrow(GModel$model), IC = GModel$IC, logLik = GModel$logL,
                  threshold = thr, cluster = as.vector(clust), fig=pl,QQplot=pl.qq)
  names(mix_gmm)[3] <- IC

  # Print the plot
  if(plot){
    p<-mix_gmm[["fig"]]
    print(p)}

  return(mix_gmm)
}

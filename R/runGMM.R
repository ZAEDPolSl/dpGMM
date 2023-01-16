#' Function to perform Gaussian Mixture Model
#'
#' Function runs GMM with stable initial points of EM algorithm and all possible control parameters within package.
#'
#' @param X Vector of data to decompose by GMM.
#' @param KS Maximum number of components to test.
#' @param Y Vector of counts, should be the same length as "X".
#' Applies only to binned data therefore the default is Y = NULL.
#' @param fixed perform GMM only for number of components given in KS.
#' @param change Stop of EM criterion (default Inf), value compared to following formula:
#' \deqn{\sum{(|\alpha - \alpha_{old})|} + \frac{\sum{(\frac{|\sigma^2 - \sigma^2_{old}|}{\sigma^2})}}{length(\alpha)}}
#' @param max_iter Maximum number of iterations of EM algorithm. By default it is \code{max_iter = 5000}
#' @param SW Minimum standard deviation of component (default 0.1). Values from 0 to 1 are treated according to their values, for SW >1 the following formula is applied:
#' \deqn{\frac{range(x)}{(SW*no.of.components))^2}}.
#' @param IC Information criterion to select best number of components.
#' Possible "AIC","AICc", "BIC" (default), "ICL-BIC" or "LR".
#' @param merge Logical value. If TRUE (default) overlapping components are merged at the distance defined in the \code{sigmas.dev} argument
#' @param sigmas.dev Defines the sigma distance to merge between the overlapping components of GMM. By default it is \code{sigma.dev = 2.5}
#' @param plot Logical value. If TRUE (default), the GMM figure will be displayed.
#' @param col.pal RColorBrewer palette name for coloring components on figure. By default \code{"Blues"}.
#' @param quick_stop Logical value. Determines to stop the EM algorithm when adding another component is no longer significant according to the Likelihood Ratio Test. Used to speed up the function (Default is TRUE).
#' @param signi Significance level for Likelihood Ratio Test. By default is 0.05.
#'
#' @returns Function returns a \code{list} which contains: \describe{
#'  \item{model}{A \code{list} of model component parameters - mean values (mu), standard deviations (sigma)
#'  and weights (alpha) for each component. Output of \code{gaussian_mixture_vector}.}
#'  \item{KS}{Optimal number of components.}
#'  \item{IC}{Value of the selected information criterion which was used to calculate the optimal number of components.}
#'  \item{logLik}{Log-likelihood value for the optimal number of components.}
#'  \item{threshold}{Vector of thresholds between each components.}
#'  \item{cluster}{Assignment of original \code{X} values to individual components (clusters).}
#'  \item{fig}{ggplot object (output of the \code{plot_gmm_1D} function). It contains decomposed distributions together with a histogram of the data.}
#'  \item{QQplot}{ggplot object (output of the \code{plot_QQplot} function).
#'  It contains fit diagnostic Quantile-Quantile plot for one normal distribution and fitted GMM.}
#' }
#'
#' @importFrom methods hasArg
#' @examples
#' \dontrun{
#' data(example)
#' mix_test1 <- runGMM(example$Dist, KS = 15, IC = "AICc", quick_stop = F, merge = F)
#' mix_test2 <- runGMM(example$Dist, KS = 10, IC = "BIC", merge = T, sigmas.dev = 1.5)
#' mix_test2$QQplot
#'
#' data(binned)
#' binned_test <- runGMM(X = binned$V1, Y = binned$V2, KS =40, col.pal ="Dark2", plot = F, quick_stop = T)
#' binned_test$fig
#' }
#'
#' @seealso \code{\link{gaussian_mixture_vector}}, \code{\link{EM_iter}}, \code{\link{generate_dist}}, \code{\link{find_thr_by_params}}
#'
#' @export
runGMM <- function(X, KS, Y = NULL, fixed=FALSE , change = Inf, max_iter = 5000, SW=0.1, IC = "BIC", merge = TRUE, sigmas.dev = 2.5,
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
    GModel<- gaussian_mixture_vector(X, KS, Y, fixed , change, max_iter, SW, IC, quick_stop, signi)

    if(merge & GModel$KS>1){
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

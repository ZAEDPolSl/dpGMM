#' Merging of overlapping components
#'
#' Function merges distributions of overlapping components at distances defined in the argument \code{sigmas.dev}.
#'
#' @param GModel A \code{data.frame} of model component parameters - the rows are, in turn:  mean values (mu), standard deviations (sigma)
#'  and weights (alpha) and each column corresponds to one component. Output of \code{EM_iter}
#' @param sigmas.dev Number of sigmas defining the distance to merge the overlapping components of GMM.
#' By default it is \code{sigma.dev = 2.5}
#'
#' @return Function returns a \code{list} which contains: \describe{
#'   \item{model}{\code{data.frame} of the GMM decomposition parameters after merging}
#'   \item{KS}{number of components after merge}
#' }
#'
#' @examples
#' GModel<-data.frame(alpha=c(0.45,0.5,0.05),
#'                   mu=c(-14,-2,5),
#'                   sigma=c(2,4,1.5))
#' gmm_merge(GModel, sigmas.dev=5)
#'
#'
#' @seealso \code{\link{runGMM}} and \code{\link{EM_iter}}
#'
#' @export
gmm_merge <- function(GModel, sigmas.dev=1.5){

  KS <- length(GModel$mu) #number of components

  merged <- list()
  merged$mu <- GModel$mu[1]
  merged$sigma <- GModel$sigma[1]
  merged$alpha <- GModel$alpha[1]

  mergedKS <- 1
  indx_merged <- 1

  for(jj in 2:KS){
    dd <- GModel$mu[jj] - merged$mu[mergedKS]
    delta <- min(GModel$sigma[jj], merged$sigma[mergedKS])*sigmas.dev

    if(dd < delta){
      pp_est <- c(GModel$alpha[jj], merged$alpha[mergedKS])
      mu_est <- c(GModel$mu[jj], merged$mu[mergedKS])
      sig_est <- c(GModel$sigma[jj], merged$sigma[mergedKS])
      ww_temp <- GModel$alpha[jj] + merged$alpha[mergedKS]
      merged$mu[mergedKS] <- (GModel$alpha[jj]*GModel$mu[jj] + merged$alpha[mergedKS]*merged$mu[mergedKS])/ww_temp
      merged$sigma[mergedKS] <- sqrt(sum(pp_est*(mu_est^2 + sig_est^2))/ww_temp - merged$mu[mergedKS]^2)
      merged$alpha[mergedKS] <- ww_temp
    }else{
      merged$mu <- c(merged$mu, GModel$mu[jj])
      merged$sigma <- c(merged$sigma, GModel$sigma[jj])
      merged$alpha = c(merged$alpha, GModel$alpha[jj])
      mergedKS <- mergedKS + 1
    }
    indx_merged <- c(indx_merged, mergedKS)
  }

  new_GModel <- list(model = as.data.frame(merged[1:3]), KS = length(merged$alpha))
  return(new_GModel)
}

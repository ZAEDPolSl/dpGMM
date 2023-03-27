#' Merging of overlapping components
#'
#' Function merges distributions of overlapping components at distances defined in the argument \code{sigmas.dev}.
#'
#' @param alpha Vector of alphas (weights) for each distribution.
#' @param mu Vector of means for each distribution.
#' @param sigma Vector of  sigmas for each distribution.
#' @param sigmas.dev Parameter used to define close GMM components that needs to be merged. For each component, standard deviation is multiplied by \code{sigmas.dev} to estimate the distance from component mean.
#' All other components within this distance are merged. By default it is \code{sigmas.dev = 2.5}. When \code{sigmas.dev = 0} no components are merged.
#'
#' @return Function returns a \code{list} which contains: \describe{
#'   \item{model}{\code{data.frame} of the GMM decomposition parameters after merging.}
#'   \item{KS}{number of components after merging.}
#' }
#'
#' @keywords internal
#'
#' @seealso \code{\link{runGMM}}
gmm_merge <- function(alpha, mu, sigma, sigmas.dev = 2.5){

  GModel <- data.frame(alpha = alpha,
                       mu = mu,
                       sigma = sigma)

  KS <- length(GModel$mu) #number of components

  merged <- list()
  merged$mu <- GModel$mu[1]
  merged$sigma <- GModel$sigma[1]
  merged$alpha <- GModel$alpha[1]

  mergedKS <- 1
  indx_merged <- 1

  for(jj in 2:KS){
    dd <- GModel$mu[jj] - merged$mu[mergedKS]
    delta <- min(GModel$sigma[jj], merged$sigma[mergedKS]) * sigmas.dev

    if(dd < delta){
      pp_est <- c(GModel$alpha[jj], merged$alpha[mergedKS])
      mu_est <- c(GModel$mu[jj], merged$mu[mergedKS])
      sig_est <- c(GModel$sigma[jj], merged$sigma[mergedKS])
      ww_temp <- GModel$alpha[jj] + merged$alpha[mergedKS]
      merged$mu[mergedKS] <- (GModel$alpha[jj] * GModel$mu[jj] + merged$alpha[mergedKS] * merged$mu[mergedKS])/ww_temp
      merged$sigma[mergedKS] <- sqrt(sum(pp_est * (mu_est^2 + sig_est^2))/ww_temp - merged$mu[mergedKS]^2)
      merged$alpha[mergedKS] <- ww_temp
    }else{
      merged$mu <- c(merged$mu, GModel$mu[jj])
      merged$sigma <- c(merged$sigma, GModel$sigma[jj])
      merged$alpha <- c(merged$alpha, GModel$alpha[jj])
      mergedKS <- mergedKS + 1
    }
    indx_merged <- c(indx_merged, mergedKS)
  }

  new_GModel <- list(model = as.data.frame(merged[1:3]), KS = length(merged$alpha))
  return(new_GModel)
}

#' QQplot for GMM decomposition
#'
#' Generator of mixed-normal distribution
#'
#' Function to generate mixed-normal distribution
#'
#' @param n Number of points to generate
#' @param alpha Vector of alphas (weights) of each GMM component
#' @param mu Vector of means of each GMM component
#' @param sigma Vector of  sigmas of each GMM component
#' @importFrom stats rnorm
#'
#' @returns List with following elements::\describe{
#'    \item{Dist}{Numeric vector with generated data}
#'    \item{Cls}{Numeric vector with labels of components to each point}
#' }
#'
#' @examples
#' data<-generate_norm1D(1000, alpha=c(0.2,0.4,0.4), mu=c(-15,0,15), sigma=c(1,2,3))
#'
#' @export
generate_norm1D <- function(n, alpha, mu, sigma){
  KS <- length(mu)
  dist <- numeric(n)
  pts.kl<-c()
  for(i in 1:n){
    pts.kl[i] <- sample.int(KS, 1L, prob=alpha)
    dist[i] <- rnorm(1, as.numeric(mu[pts.kl[i]]), as.numeric(sigma[pts.kl[i]]))
  }

  idx<-order(dist,decreasing = F)
  pts.kl<-pts.kl[idx]
  dist<-dist[idx]

  res<-list(Dist=dist,Cls=pts.kl)
  return(res)
}

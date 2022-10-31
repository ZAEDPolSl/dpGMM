#' Generating alpha for GMM distribution
#'
#' Function generating alpha values of mixture distribution. Sum of alpha equals 1.
#'
#' @param n Number of alpha sets
#' @param m Number of random alphas
#'
#' @importFrom stats runif
#' @importFrom Matrix rowSums
#'
#' @examples
#' alpha_rand(1,4)
#'
alpha_rand <- function(n,m) {
  ri <- matrix(stats::runif(m*n,0,1), ncol=m)
  ri<- sweep( ri, 1, Matrix::rowSums( ri), FUN="/")
  ri
}


#' Generation of GMM data with high precision
#'
#' Function to generate PDF of GMM distributions and its cumulative results with high lincespacing.
#'
#' @param data Vector of original data
#' @param GModel \code{data.frame} of GMM parameters i.e GModel$alpha, GModel$mu, GModel$sigma (correct \code{colnames} are obligatory)
#' @param precision Precision of point linespacing
#'
#' @importFrom stats dnorm
#' @importFrom pracma linspace
#' @importFrom Matrix rowSums
#'
#' @returns List with following elements::\describe{
#'    \item{x}{Numeric vector with equaliy spread data of given precison}
#'    \item{dist}{Matrix with PDF of each GMM component and cumulative distribution}
#' }
#'
#'
#' @seealso \code{\link{runGMM}} and \code{\link{generate_norm1D}}
#' @export
generate_dist<-function(data, GModel, precision){
  x_temp = pracma::linspace(min(data),max(data),precision)
  f_temp = matrix(0, precision, nrow(GModel))
  for(k in 1:nrow(GModel)){
    f_temp[,k] = GModel$alpha[k] * stats::dnorm(x_temp, mean = GModel$mu[k], sd =GModel$sigma[k])
  }

  f_temp <- as.data.frame(f_temp)
  f_temp$main <- Matrix::rowSums(f_temp)

  return(list(x=x_temp,dist=f_temp))
}

#' Generator of mixed-normal distributions
#'
#' Generator of mixed-normal distribution with given model parameters for certain points number.
#'
#' @param n Number of points to generate
#' @param alpha Vector of alphas (weights) for each distribution
#' @param mu Vector of means for each distribution
#' @param sigma Vector of  sigmas for each distribution
#'
#' @importFrom stats rnorm
#'
#' @returns List with following elements::\describe{
#'    \item{Dist}{Numeric vector with generated data}
#'    \item{Cls}{Numeric vector with classification of each point to particular mixed distribution}
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
    dist[i] <- stats::rnorm(1, as.numeric(mu[pts.kl[i]]), as.numeric(sigma[pts.kl[i]]))
  }

  idx<-order(dist,decreasing = F)
  pts.kl<-pts.kl[idx]
  dist<-dist[idx]

  res<-list(Dist=dist,Cls=pts.kl)
  return(res)
}

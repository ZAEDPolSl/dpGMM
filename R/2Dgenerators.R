#' Generator of 2D Gaussian mixture distribution
#'
#' Generator of 2D mixed normal distribution with given model parameters for certain points number.
#'
#' @param n Number of points to generate
#' @param alpha Vector of alphas (weights) for each distribution
#' @param mu Matrix of means for each distribution
#' @param cov Vector of covariances for each distribution.
#'
#' @importFrom mvtnorm rmvnorm
#'
#' @returns List with following elements::\describe{
#'    \item{Dist}{Numeric marix with generated data}
#'    \item{Cls}{Numeric vector with classification of each point to particular distribution}
#' }
#'
#' @examples
#' \dontrun{
#' data<-generate_norm2D(1500, alpha=c(0.2,0.4,0.4), mu=matrix(c(1,2,1,3,2,2),nrow=2),cov =c(0.01,0.02,0.03))
#' }
#' @seealso \code{\link{generate_dset2D}}
#' @export
generate_norm2D <- function(n, alpha, mu, cov){
  KS <- dim(mu)[2]
  dist <- matrix(NA,nrow=n,ncol=2)
  pts.kl <- numeric(n)
  for(i in 1:n){
    pts.kl[i] <- sample.int(KS, 1L, prob=alpha)
    dist[i,] <- mvtnorm::rmvnorm(1, mu[,pts.kl[i]], diag(2)*(cov[pts.kl[i]]), method="svd" )
  }
  res<-list(Dist=dist,Cls=pts.kl)
  return(res)
}

#' Create multiple random 2D Gaussian mixture datasets
#'
#' Generator of multiple 2D mixed normal distribution with given model parameters ranges.
#'
#' @param n Number of points to generate.
#' @param m Number of distribution to generate.
#' @param KS_range Range of possible number of components of generated distribution. Default \code{KS=2:8}.
#' @param mu_range Range of means of components of generated distribution. Default \code{-15:15}.
#' @param cov_range Range of means of components of generated distribution. Default \code{1:5}.
#'
#' @importFrom stats runif
#'
#' @returns List with 2D GMM distributions where each list contains elements of \code{\link{generate_norm2D}}

#'
#' @seealso \code{\link{generate_norm2D}}
#' @export
generate_dset2D <- function(n=1500,m=1500,KS_range=2:8,mu_range=c(-15,15),cov_range=c(1,5)){

  res<- list()
  for(i in 1:m){
    res_tmp <- list()

    # get no. of components
    res_tmp[["KS"]] <- sample(KS_range,1)

    # randomly generate components' parameters
    res_tmp[["mu"]] <- rbind(stats::runif(res_tmp[["KS"]],mu_range[1],mu_range[2]),
                             stats::runif(res_tmp[["KS"]],mu_range[1],mu_range[2]))
    res_tmp[["sigma"]] <- stats::runif(res_tmp[["KS"]],cov_range[1],cov_range[2])
    res_tmp[["alpha"]] <- stats::runif(res_tmp[["KS"]],0,1)
    res_tmp[["alpha"]] <- res_tmp[["alpha"]]/sum(res_tmp[["alpha"]])

    # generate data based on parameters
    tmp <- generate_norm2D(n, res_tmp[["alpha"]] , res_tmp[["mu"]] , res_tmp[["sigma"]] )
    res_tmp[["dist"]] <- tmp

    # store data
    res[[paste0("sim",i)]] <- res_tmp
  }
  return(res)
}

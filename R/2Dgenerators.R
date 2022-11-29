#' Generator of 2D Gaussian mixture distributions
#'
#' Generator of 2D mixed normal distribution with given model parameters for certain points number.
#'
#' @param n Number of points to generate
#' @param alpha Vector of alphas (weights) for each distribution
#' @param mu Vector of means for each distribution
#' @param sigma Vector of  sigmas for each distribution.
#'
#' @importFrom mvtnorm rmvnorm
#'
#' @returns List with following elements::\describe{
#'    \item{Dist}{Numeric vector with generated data}
#'    \item{Cls}{Numeric vector with classification of each point to particular distribution}
#' }
#'
#' @examples
#' \dontrun{
#' data<-generate_norm2D(1500, alpha=c(0.2,0.4,0.4), mu=matrix(c(1,2,1,3,2,2), nrow=2),
#' cov =c(0.01,0.02,0.03))
#' }
#' @seealso \code{\link{generate_dset2D}} and \code{\link{generate_norm1D}}
#' @export
generate_norm2D <- function(n, alpha, mu, cov){
  KS <- dim(mu)[2]
  dist <- matrix(NA,nrow=n,ncol=2)
  pts.kl <- numeric(n)
  for(i in 1:n){
    pts.kl[i] <- sample.int(KS, 1L, prob=alpha)
    dist[i,] <- mvtnorm::rmvnorm(1, mu[,pts.kl[i]], diag(2)*(cov[pts.kl[i]]), method="svd" )
  }

  # par(mfrow=c(1,3))
  # plot(dist[,1],dist[,2],main="Points")
  # plot(density(dist[,1]),main="Dens X")
  # plot(density(dist[,2]),main="Dens Y")

  res<-list(Dist=dist,Cls=pts.kl)
  return(res)
}

#' Create random 2D dataset
#' @export
generate_dset2D <- function(n=1500,m=1500,KS_range=2:8,mu_range=c(-15,15),sig_range=c(1,5)){

  res<- list()
  for(i in 1:m){
    res_tmp <- list()

    # get no. of components
    res_tmp[["KS"]] <- sample(KS_range,1)

    # randomly generate components' parameters
    res_tmp[["mu"]] <- rbind(runif(res_tmp[["KS"]],mu_range[1],mu_range[2]),
                             runif(res_tmp[["KS"]],mu_range[1],mu_range[2]))
    res_tmp[["sigma"]] <- runif(res_tmp[["KS"]],sig_range[1],sig_range[2])
    res_tmp[["alpha"]] <- runif(res_tmp[["KS"]],0,1)
    res_tmp[["alpha"]] <- res_tmp[["alpha"]]/sum(res_tmp[["alpha"]])

    # generate data based on parameters
    tmp <- generate_norm2D(n, res_tmp[["alpha"]] , res_tmp[["mu"]] , res_tmp[["sigma"]] )
    res_tmp[["dist"]] <- tmp

    # store data
    res[[paste0("sim",i)]] <- res_tmp
  }
  return(res)
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
# generate_dist2D<-function(data, GModel, precision){
#   x_temp = pracma::linspace(min(data),max(data),precision)
#   f_temp = matrix(0, precision, nrow(GModel))
#   for(k in 1:nrow(GModel)){
#     f_temp[,k] = GModel$alpha[k] * stats::dnorm(x_temp, mean = GModel$mu[k], sd =GModel$sigma[k])
#   }
#
#   f_temp <- as.data.frame(f_temp)
#   f_temp$main <- Matrix::rowSums(f_temp)
#
#   return(list(x=x_temp,dist=f_temp))
# }


#' Generation of GMM data with high precision
#'
#' Function to generate PDF of GMM distributions and its cumulative results with high lincespacing.
#'
#' @param X Vector of 1D data.
#' @param alpha Vector of alphas (weights) for each distribution.
#' @param mu Vector of means for each distribution.
#' @param sigma Vector of  sigmas for each distribution.
#' @param precision Precision of point linespacing.
#'
#' @importFrom stats dnorm
#' @importFrom pracma linspace
#' @importFrom Matrix rowSums
#'
#' @returns List with following elements:\describe{
#'    \item{x}{Numeric vector with equaliy spread data of given precison.}
#'    \item{dist}{Matrix with PDF of each GMM component and cumulative distribution.}
#' }
#'
#'
#' @seealso \code{\link{runGMM}} and \code{\link{generate_norm1D}}
#' @export
generate_dist <- function(X, alpha, mu, sigma, precision){

  GModel <- data.frame(alpha = alpha,
                       mu = mu,
                       sigma = sigma)

  x_temp <- pracma::linspace(min(X), max(X), precision)
  f_temp <- matrix(0, precision, nrow(GModel))

  for(k in 1:nrow(GModel)){
    f_temp[,k] <- GModel$alpha[k] * stats::dnorm(x_temp, mean = GModel$mu[k], sd =GModel$sigma[k])
  }

  f_temp <- as.data.frame(f_temp)
  f_temp$main <- Matrix::rowSums(f_temp)

  return(list(x = x_temp, dist = f_temp))
}

#' Generator of 1D mixed-normal distributions
#'
#' Generator of mixed-normal distribution with given model parameters for certain points number.
#'
#' @param n Number of points to generate.
#' @param alpha Vector of alphas (weights) for each distribution.
#' @param mu Vector of means for each distribution.
#' @param sigma Vector of  sigmas for each distribution.
#'
#' @importFrom stats rnorm
#'
#' @returns List with following elements:\describe{
#'    \item{Dist}{Numeric vector with generated data}
#'    \item{Cls}{Numeric vector with classification of each point to particular mixed distribution}
#' }
#'
#' @examples
#' \dontrun{
#' data <- generate_norm1D(1000, alpha = c(0.2, 0.4, 0.4), mu = c(-15, 0, 15), sigma = c(1, 2, 3))
#' }
#' @export
generate_norm1D <- function(n, alpha, mu, sigma){
  KS <- length(mu)
  dist <- numeric(n)
  pts.kl <- c()
  for(i in 1:n){
    pts.kl[i] <- sample.int(KS, 1L, prob = alpha)
    dist[i] <- stats::rnorm(1, as.numeric(mu[pts.kl[i]]), as.numeric(sigma[pts.kl[i]]))
  }

  idx <- order(dist, decreasing = F)
  pts.kl <- pts.kl[idx]
  dist <- dist[idx]

  res <- list(Dist = dist, Cls = pts.kl)
  return(res)
}

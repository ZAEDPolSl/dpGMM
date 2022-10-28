#' Generating alpha for GMM distribution
#'
#' Function generating alpha values of mixture distribution. Sum of alpha equals 1.
#'
#' @param n Number of alpha sets
#' @param m Number of random alphas
#'
#' @importFrom stats runif
#'
#' @examples
#' alpha_rand(1,4)
#'
#' @export
alpha_rand <- function(n,m) {
  ri <- matrix(runif(m*n,0,1), ncol=m)
  ri<- sweep( ri, 1, rowSums( ri), FUN="/")
  ri
}

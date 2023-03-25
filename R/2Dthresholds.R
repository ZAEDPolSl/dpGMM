#' Class assignment for 2D Gaussian Mixture Model data
#'
#' Function which assign each point of 2D matrix data to a cluster by maximum probability.
#'
#' @param X matrix of data to decompose by GMM.
#' @param gmm Results of \code{\link{gaussian_mixture_2D}} decomposition.
#'
#' @returns Return a vector of cluster assignment of each point of X matrix.
#'
#' @export
find_class_2D <- function(X, gmm){

  KS <- gmm$KS
  dist <- matrix(NA, nrow = dim(X)[1], ncol = KS+1)
  for (i in 1:KS){
      center <- gmm$center[i,]
      covariance <- gmm$covar[,,i]
      dist[,i] <- rGMMtest:::norm_pdf_2D(X, center, covariance)  ### PACKAGE NAME!!!!!!!!!!!!!!!!
  }
  dist[, KS+1] <- rowSums(dist[,1:KS])
  dist <- as.data.frame(dist)
  dist$cluster <- apply(dist[,1:KS], 1, function(x) which.max(x))

  return(dist$cluster)
}

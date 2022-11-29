#' Plot of GMM decomposition for 2D data
#'
#' Function plot the decomposed distribution together with histogram of data.
#' This plot is also return as regular output of \code{\link{runGMM}}.
#'
#' @param X Vector of data
#' @param dist Output of \code{generate_dist} function. Its list of following elements\describe{
#'    \item{x}{Vector of generated data}
#'    \item{dist}{Matrix of pdf of for each generated mixture model. The last column is sum of all previous ones}
#' }
#' @param Y Vector of X counts (dedicated to binned data). Default=NULL
#' @param threshold Vector with GMM cutoffs
#' @param pal RColorBrewer palette name
#'
#' @import ggplot2
#' @import RColorBrewer
#' @importFrom grDevices colorRampPalette
#' @importFrom reshape2 melt
#'
#' @seealso \code{\link{runGMM}}
#'
#' @examples
#' \dontrun{
#' data(example)
#' GModel<-data.frame(alpha=c(0.45,0.5,0.05),
#'                   mu=c(-14,-2,5),
#'                   sigma=c(2,4,1.5))
#' dist.plot<-generate_dist(example$Dist, GModel, 1e4)
#' thr <- find_thr_by_params(GModel,dist.plot)
#' plot_gmm_1D(example$Dist, dist.plot,Y=NULL,threshold=thr, pal="Dark2")
#' }
#'
#' @export
plot_gmm_2D_binned <- function(data, img, gmm, opts){
  #Plot 2D binned data versus GMM decomposition.
  m <- dim(img)[1];  n <- dim(img)[2]
  ploty <- matrix(0, m, n)

  for(c in 1:gmm$KS){
    tmp <- gmm$alpha[c]*(norm_pdf_2D(data, gmm$center[c,], gmm$covar[,,c]))
    ploty <- ploty + matrix(tmp, m, n)
  }
  scale <- sum(255-img)/sum(ploty)

  #plot z dodawaniem kszta??w w petli z plot_2DGauss na zdjecie
  cov_type <- opts$cov_type
  coors <- data.frame()
  for (a in 1:gmm$KS){
    center <- gmm$center[a,]- c(min(data[,1]) - 1, min(data[,2]) - 1)
    covariance <- pracma:::rot90(gmm$covar[,,a], 2)
    tmp <- plot_2DGauss(center, covariance, cov_type)
    tmp$KS <- rep(a, 100)
    coors <- rbind(coors, tmp)
  }

  p <- ggplot(coors, aes(x = ellipse_x_r, y = ellipse_y_r, group = KS)) + geom_path() +theme_bw()

  return(list(scale,p))
}

#' @export
plot_gmm_2D_orig <- function(X, gmm, opts){

  #Plot 2D data versus GMM decomposition.
  cov_type <- opts$cov_type
  coors <- data.frame()
  for (a in 1:gmm$KS){
    center <- gmm$center[a,]- c(min(X[,1]) - 1, min(X[,2]) - 1)
    covariance <- pracma::rot90(gmm$covar[,,a], 2)
    tmp <- plot_2DGauss(center, covariance, cov_type)
    tmp$KS <- rep(a, 100)
    coors <- rbind(coors, tmp)
  }

  p <- ggplot(coors, aes(x = ellipse_x_r, y = ellipse_y_r, group = KS)) + geom_path() +theme_bw()

  return(list(scale,p))
}

#' QQplot for GMM decomposition
#'
#' Function return ggplot object with fit diagnostic Quantile-Quantile plot for one normal distribution and fitted GMM.
#' This plot is also return as regular output of \code{\link{runGMM}}.
#'
#' @param data Vector of original data
#' @param GModel \code{data.frame} of GMM parameters i.e GModel$alpha, GModel$mu, GModel$sigma (correct \code{colnames} are obligatory)
#'
#' @import ggplot2
#' @importFrom  ggpubr ggarrange
#' @importFrom stats qqplot
#'
#' @examples
#' \dontrun{
#' data(example)
#' GModel<-data.frame(alpha=c(0.45,0.5,0.05),
#'                   mu=c(-14,-2,5),
#'                   sigma=c(2,4,1.5))
#' plot_QQplot(example$Dist,GModel)
#' }
#'
#' @seealso \code{\link{runGMM}} and \code{\link{gaussian_mixture_vector}}
#'
#' @export
plot_2DGauss <- function(center, covariance, cov_type){
  e <- eigen(covariance)
  eigenvec <- apply(e$vectors, 1, rev)
  eigenval <- diag(rev(e$values))

  #Get the index of the largest eigenvector
  largest_eigenvec_ind_c <- which(eigenval == max(eigenval), arr.ind = T)[1]
  largest_eigenvec <- eigenvec[,largest_eigenvec_ind_c]

  #Get the largest eigenvalue
  largest_eigenval <- max(eigenval)

  #Get the smallest eigenvector and eigenvalue
  if (largest_eigenvec_ind_c == 1){
    smallest_eigenval <- max(eigenval[,2])
    smallest_eigenvec <- eigenvec[,2]
  } else {
    smallest_eigenval <- max(eigenval[,1])
    smllest_eigenvec <- eigenvec[,1]
  }

  #Calculate the angle between the x-axis and the largest eigenvector
  angle <- atan2(largest_eigenvec[2], largest_eigenvec[1])

  # This angle is between -pi and pi.
  # Let's shift it such that the angle is between 0 and 2pi
  if (angle < 0){
    angle <- angle + 2*pi
  }

  #Get the confidence interval error ellipse
  #chisquare_val = 3.0349;   % 99%
  if (cov_type == "sphere"){
    chisquare_val <- sqrt(qchisq(.95, 1))
  } else{
    chisquare_val <- sqrt(qchisq(.95, 2))
  }

  theta_grid <- seq(from=0, to=2*pi, length.out=100)
  phi <- angle
  x0 <- center[1]
  y0 <- center[2]
  a <- chisquare_val*sqrt(largest_eigenval)
  b <- chisquare_val*sqrt(smallest_eigenval)

  #the ellipse in x and y coordinates
  ellipse_x_r <- a*cos(theta_grid)
  ellipse_y_r <- b*sin(theta_grid)

  #Define a rotation matrix
  R <- rbind(c(cos(phi), sin(phi)), c(-sin(phi), cos(phi)))

  #let's rotate the ellipse to some angle phi
  r_ellipse <- t(apply(cbind(ellipse_x_r, ellipse_y_r), 1, "*", rowSums(R)))
  r_ellipse[,1] <- r_ellipse[,1] + y0; r_ellipse[,2] <- r_ellipse[,2]+x0

  return(as.data.frame(r_ellipse))
}

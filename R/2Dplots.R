#' Plot of GMM decomposition for 2D binned data
#'
#' Function plot the decomposed distribution together with histogram of data.
#' This plot is also return as regular output of \code{\link{runGMM}}.
#'
#' @param X matrix of data to decompose by GMM.
#' @param Y vector of counts, should be the same length as "X".
#' @param gmm results of \code{\link{gaussian_mixture_2D}} decomposition
#' @param opts parameters of run saves in \code{\link{GMM_2D_opts}} variable
#'
#' @import ggplot2
#' @import RColorBrewer
#' @import patchwork
#' @importFrom grDevices colorRampPalette
#' @importFrom pracma rot90
#'
#'
#' @examples
#' \dontrun{
#' }
#'
#' @export
plot_gmm_2D_binned <- function(X, Y, gmm, opts){
   colnames(X)<-c("Coordinates_1","Coordinates_2")
  cov_type <- opts$cov_type
  crits<-c(0.15,0.35)
  elps<-list()
  for (j in 1:length(crits)){
    coors <- data.frame()
    for (a in 1:gmm$KS){
      center <- gmm$center[a,]#- c(min(X[,1]) - 1, min(X[,2]) - 1)
      covariance <- pracma::rot90(gmm$covar[,,a], 2)
      tmp <- ellips2D(center, covariance, cov_type, crits[j])
      tmp$KS <- rep(a, 100)
      coors <- rbind(coors, tmp)
    }
    elps[[j]]<-coors
  }


  col <- grDevices::colorRampPalette(brewer.pal(8,"Dark2"))(gmm$KS)


  p<-ggplot() +theme_bw()+
    geom_tile(aes(x = X$Coordinates_1, y = X$Coordinates_2, fill = Y),show.legend = F)+
    scale_fill_viridis_c()+
    geom_point(aes(x=gmm$center[,1],y=gmm$center[,2]),color="red",size=3)+
    xlab("X1")+ylab('X2')+
    geom_path(aes(x = elps[[1]]$V2, y = elps[[1]]$V1, group = coors$KS),color="black",show.legend = F,size=0.75,linetype="dashed")+
    geom_path(aes(x = elps[[2]]$V2, y = elps[[2]]$V1, group = coors$KS),color="black",show.legend = F,size=0.75,linetype="dashed")

    return(p)

}



#' Plot of GMM decomposition for 2D data
#'
#' Function plot the decomposed distribution together with histogram of data.
#' This plot is also return as regular output of \code{\link{runGMM}}.
#'
#' @param X Matrix of data
#' @param gmm results of \code{\link{gaussian_mixture_2D}} decomposition
#' @param opts parameters of run stored in \code{\link{GMM_2D_opts}} variable
#'
#' @import ggplot2
#' @import RColorBrewer
#' @importFrom grDevices colorRampPalette
#' @importFrom pracma rot90
#'
#'
#' @examples
#' \dontrun{
#' }
#'
#' @export
plot_gmm_2D_orig <- function(X, gmm, opts){
  X<-as.data.frame(X)
  colnames(X)<-c("X1","X2")
  X$cls<-gmm$cls
  #Plot 2D data versus GMM decomposition.
  cov_type <- opts$cov_type
  crits<-c(0.25,0.75,0.95)
  elps<-list()
    for (j in 1:length(crits)){
      coors <- data.frame()
      for (a in 1:gmm$KS){
        center <- gmm$center[a,]#- c(min(X[,1]) - 1, min(X[,2]) - 1)
        covariance <- pracma::rot90(gmm$covar[,,a], 2)
        tmp <- ellips2D(center, covariance, cov_type, crits[j])
        tmp$KS <- rep(a, 100)
        coors <- rbind(coors, tmp)
      }
    elps[[j]]<-coors
    }


  col <- grDevices::colorRampPalette(brewer.pal(8,"Dark2"))(gmm$KS)

  p<-ggplot()
  p<-p+geom_point(aes(x = X$X1, y = X$X2,alpha=0.75,color=factor(X$cls)),show.legend = F)+
    scale_color_manual(values=col)+
       geom_path(aes(x = elps[[1]]$V2, y = elps[[1]]$V1, group = coors$KS),color="grey35",show.legend = F,size=0.65,linetype="dashed")+ #color=factor(coors$KS)
       geom_path(aes(x = elps[[2]]$V2, y = elps[[2]]$V1, group = coors$KS),color="grey35",show.legend = F,size=0.65,linetype="dashed")+
       geom_path(aes(x = elps[[3]]$V2, y = elps[[3]]$V1, group = coors$KS),color="grey35",show.legend = F,size=0.65,linetype="dashed")
  p<-p+geom_point(aes(x=gmm$center[,1],y=gmm$center[,2]),color="black",size=3)+xlab("X1")+ylab('X2')+
    theme_bw()+theme(plot.title = element_text(hjust = 0.5))
  p
 # p<-ggMarginal(p, type = "density",color="#324376", xparams = list(size = 1),yparams = list(size = 1))
  return(p)
}

#' Ellipses for plot 2D GMM
#'
#' Function for defining ellipses in 2D GMM plot by the confidence interval.
#'
#' @param center Means of decomposition
#' @param covariance Covariances of each component
#' @param cov_type Type of covariance model same as in \code{\link{GMM_2D_opts}}. Possible "sphere","diag" or "full" (default).
#' @param crit Confidence interval level of ellipse. Default 0.95.
#'
#' @export
ellips2D <- function(center, covariance, cov_type, crit=0.95){
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
    chisquare_val <- sqrt(qchisq(crit, 1))
  } else{
    chisquare_val <- sqrt(qchisq(crit, 2))
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
  r_ellipse <- cbind(ellipse_x_r, ellipse_y_r)%*%R
  r_ellipse[,1] <- r_ellipse[,1]+ y0;
  r_ellipse[,2] <- r_ellipse[,2]+ x0

  return(as.data.frame(r_ellipse))
}

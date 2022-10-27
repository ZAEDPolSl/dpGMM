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
  
  theta_grid <- linspace(0, 2*pi)
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
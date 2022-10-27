plot_gmm_2D <- function(data, img, gmm, opts){
  #Plot 2D gel image versus GMM decomposition. 
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
    covariance <- rot90(gmm$covar[,,a], 2)
    tmp <- plot_2DGauss(center, covariance, cov_type)
    tmp$KS <- rep(a, 100)
    coors <- rbind(coors, tmp)
  }
  
  p <- ggplot(coors, aes(x = ellipse_x_r, y = ellipse_y_r, group = KS)) + geom_path() +theme_bw()
  
  return(list(scale,p))
}
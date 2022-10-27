diag_init_2D <- function(data, KS){
  init <- list()
  init$alpha <- matrix(1, 1, KS)/KS #equal mixing proportions
  
  n_row <- data[nrow(data),1]
  n_col <- data[nrow(data),2]
  a <- cbind(seq(1,n_row, len = (KS*2)+1), seq(1,n_col, len = (KS*2)+1))
  colnames(a) <- colnames(data)
  init$center <- matrix(round(a[seq(2, nrow(a), 2),]), ncol=2)
  init$center <- apply(init$center, c(1,2), as.integer)
  
  init$covar <- replicate(KS, diag(apply(as.matrix(data), 2, sd)/KS), simplify = "array")
  init$KS <- KS
  return(init)
}
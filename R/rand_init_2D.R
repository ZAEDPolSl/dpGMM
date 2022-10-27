rand_init_2D <- function(data, KS){
  # Random initialization of 2D Gaussian components.
  n <- dim(data)[1]
  init <- list()
  init$alpha <- matrix(1, 1, KS)/KS #equal mixing proportions
  init$center <- as.matrix(data[sample(n, KS),])
  rownames(init$center) <- NULL
  init$covar <- replicate(KS, diag(apply(as.matrix(data), 2, sd)/KS), simplify = "array")
  init$KS <- KS
  return(init)
}
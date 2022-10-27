generate_norm2D <- function(n, alpha, mu, cov){
  
  KS <- dim(mu)[2]
  dist <- matrix(NA,nrow=n,ncol=KS)
  for(i in 1:n){
    k <- sample.int(KS, 1L, prob=alpha)
    dist[i,] <- rmvnorm(1, as.numeric(mu[,k]), cov, method="svd")
  }
  return(dist)
}
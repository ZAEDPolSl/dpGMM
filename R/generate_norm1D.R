generate_norm1D <- function(n, alpha, mu, sigma){
  KS <- length(mu)
  dist <- numeric(n)
  pts.kl<-c()
  for(i in 1:n){
    pts.kl[i] <- sample.int(KS, 1L, prob=alpha)
    dist[i] <- rnorm(1, as.numeric(mu[pts.kl[i]]), as.numeric(sigma[pts.kl[i]]))
  }
  
  idx<-order(dist,decreasing = F)
  pts.kl<-pts.kl[idx]
  dist<-dist[idx]
  
  res<-list(Dist=dist,Cls=pts.kl)
  return(res)
}
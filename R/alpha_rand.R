alpha_rand <- function(n,m) {
  ri <- matrix(runif(m*n,0,1), ncol=m)
  ri<- sweep( ri, 1, rowSums( ri), FUN="/")
  ri
}
dyn_pr_split_w_aux <- function(data, ygreki){
  N <- length(data)
  #aux_mx
  aux_mx <- matrix(0, N, N)
  for (kk in 1:(N-1)){
    for (jj in (kk+1):N){
      aux_mx[kk,jj] <- my_qu_ix_w(data[kk:(jj-1)], ygreki[kk:(jj-1)])
    }
  }
  return(aux_mx)
}
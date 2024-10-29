#' Dynamic programming split.
#'
#' Supporting function for dynamic programming in 1D decomposition.
#'
#' @param xhist x coordinates of histogram.
#' @param yhist y coordinates of histogram.
#' @param K Number of components.
#' @param aux_mx Return of function \code{dyn_pr_split_w_aux}.
#'
#' @return List of parameters.
#'
#' @keywords internal
#'
dyn_pr_split_w <- function(xhist, yhist, K, aux_mx,s_corr){
  #initialize
  Q <- matrix(0, 1, K)
  N <- length(xhist)
  p_opt_idx <- matrix(0, 1, N)
  p_aux <- matrix(0, 1, N)
  opt_pals <- matrix(0, K, N)
  PAR_min_seg<- 5

  for (kk in 1:N){
    if (sum(yhist[kk:N])<PAR_min_seg){
      p_opt_idx[kk]<-Inf
    } else{
      p_opt_idx[kk] <- my_qu_ix_w(xhist[kk:N], yhist[kk:N],s_corr)
    }

  }

  #iterate
  for (kster in 1:K){
    for (kk in 1:(N-kster)){
      for (jj in (kk+1):(N-kster+1)){
        p_aux[jj] <- aux_mx[kk, jj] + p_opt_idx[jj]
      }
      mm <- min(p_aux[(kk+1):(N-kster+1)])
      ix <- which.min(p_aux[(kk+1):(N-kster+1)])
      p_opt_idx[kk] <- mm
      opt_pals[kster,kk] <- kk+ix[1]
    }
    Q[kster] <- p_opt_idx[1]
  }

  #restore optimal decisions
  opt_part <- matrix(0, 1, K)
  opt_part[1] <- opt_pals[K,1]

  if (K != 1){
    for(kster in rev(1:(K-1))){
      opt_part[K-kster+1] <- opt_pals[kster, opt_part[K-kster]]
    }
  }
  return(list(Q, opt_part))
}



#' Dynamic programming split of aux.
#'
#' Supporting function for dynamic programming in 1D decomposition.
#'
#' @param xhist x coordinates of histogram.
#' @param yhist y coordinates of histogram.
#'
#' @keywords internal
#'
dyn_pr_split_w_aux <- function(xhist, yhist,s_corr){
  N <- length(xhist)
  #aux_mx
  aux_mx <- matrix(0, N, N)
  for (kk in 1:(N-1)){
    for (jj in (kk+1):N){
      aux_mx[kk,jj] <- dpGMM:::my_qu_ix_w(xhist[kk:(jj-1)], yhist[kk:(jj-1)],s_corr)
    }
  }
  return(aux_mx)
}

#' Supporting function for spliters.
#'
#' Supporting function used in spliters \code{dyn_pr_split_w_aux} and \code{dyn_pr_split_w}.
#'
#' @param xinvec x coordinates of histogram.
#' @param yinvec y coordinates of histogram.
#'
#' @keywords internal
#'
my_qu_ix_w <- function(xinvec, yinvec, s_corr){
  PAR <- 1
  PAR_sig_min <- .1
  if ((xinvec[length(xinvec)] - xinvec[1]) <= PAR_sig_min || sum(yinvec) <= 1.0e-3){
    wyn <- Inf
  }else {
    wwec <- yinvec/sum(yinvec)
    var_bin = sum(((xinvec - sum(xinvec*wwec))^2) * wwec);
    if (var_bin>s_corr){
      wyn <- (PAR + sqrt(var_bin-s_corr))/(max(xinvec) - min(xinvec))
    } else{
      wyn<-Inf
    }
  }

  return(wyn)
}

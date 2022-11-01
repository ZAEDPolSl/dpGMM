#' NIE MAM ZIELONEGO POJĘCIA CO ROBIĄ TE FUNKCJE
#'

dyn_pr_split_w <- function(data, ygreki, K_gr, aux_mx){
  #initialize
  Q <- matrix(0, 1, K_gr)
  N <- length(data) #to jest length od x w hist
  p_opt_idx <- matrix(0, 1, N)
  p_aux <- matrix(0, 1, N)
  opt_pals <- matrix(0, K_gr, N)

  for (kk in 1:N){
    p_opt_idx[kk] <- my_qu_ix_w(data[kk:N], ygreki[kk:N])
  }

  #iterate
  for (kster in 1:K_gr){
    for (kk in 1:(N-kster)){
      for (jj in (kk+1):(N-kster+1)){
        p_aux[jj] <- aux_mx[kk, jj] + p_opt_idx[jj]
      }
      mm <- min(p_aux[(kk+1):(N-kster+1)])
      # ix <- order(p_aux[(kk+1):(N-kster+1)]==mm)
      ix <- which.min(p_aux[(kk+1):(N-kster+1)])
      # mm <- apply(p_aux[kk+1:N-kster+1], 2, min)
      # ix <- which(p_aux[kk+1:N-kster+1] == mm, arr.ind = T)[,"row"]
      p_opt_idx[kk] <- mm
      opt_pals[kster,kk] <- kk+ix[1]
    }
    Q[kster] <- p_opt_idx[1]
  }

  #restore optimal decisions
  opt_part <- matrix(0, 1, K_gr)
  opt_part[1] <- opt_pals[K_gr,1]

  if (K_gr != 1){
    for(kster in rev(1:(K_gr-1))){
      opt_part[K_gr-kster+1] <- opt_pals[kster, opt_part[K_gr-kster]]
    }
  }
  return(list(Q, opt_part))
}

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

my_qu_ix_w <- function(invec, yinwec){
  PAR = 1
  PAR_sig_min =.1
  if ((invec[length(invec)]-invec[1]) <= PAR_sig_min || sum(yinwec) <= 1.0e-3){
    wyn = Inf
  }else {
    wwec = yinwec/sum(yinwec)
    wyn1 = (PAR + sqrt(sum(((invec - sum(invec*wwec))^2)*wwec)))/(max(invec) - min(invec))
    wyn = wyn1
  }

  return(wyn)
}


#' Generating the probability density function.
#'
#' Function generating PDF values of a Gaussian distribution.
#'
#' @param x Numeric vecotr for which to calculate the PDF.
#' @param mu Mean value of the distribution.
#' @param sigma Standard deviation value of the distribution.
#'
#' @return Returns a vector of PDF value
#'
#' @example
#' data(example)
#' pdf_value <- norm_pdf(data, 2.5, 0.75)
#'
#' @export

norm_pdf <- function(x, mu, sigma){
  y<- exp(-0.5*(((x-mu)/sigma)^2))/(2.506628274631*sigma)
  return(y)
}

#' Random initialization of 2D GMM EM
#'
#' Function for generating random initial conditions of 2D GMM model.
#'
#' @param X matrix of data to decompose by GMM.
#' @param KS Number of components.
#'
#' @return Initial values for EM.
#'
#' @keywords internal
#'
rand_init_2D <- function(X, KS){
  # Random initialization of 2D Gaussian components.
  n <- dim(X)[1]
  init <- list()
  init$alpha <- matrix(1, 1, KS)/KS #equal mixing proportions
  init$center <- as.matrix(X[sample(n, KS),])
  rownames(init$center) <- NULL
  init$covar <- replicate(KS, diag(apply(as.matrix(X), 2, sd)/KS), simplify = "array")
  init$KS <- KS
  return(init)
}

#' Dynamic programming initialization of 2D GMM EM
#'
#' Function for generating initial conditions of 2D GMM model by dynamic programming.
#'
#' @param X Matrix of data to decompose by GMM.
#' @param Y Vector of counts, with the same length as "X".
#' @param K Number of components.
#'
#' @return Initial values for EM.
#'
#' @keywords internal
#'
DP_init_2D <- function(X, Y, KS){
  init <- list()

  # create boundary distributions
  if (sum(Y) > length(Y)){
    B_dist1 <- aggregate(x = list("y" = Y), by = list("x" = X[,1]), FUN = sum)
    B_dist2 <- aggregate(x = list("y" = Y), by = list("x" = X[,2]), FUN = sum)
  } else{
    h  <- hist(X[,1], breaks = seq(min(X[,1]), max(X[,1]), l=(min(max(20,round(sqrt(nrow(X)))), 100)+1)),plot = F)
    B_dist1 <- cbind(h$mids, h$counts)
    h  <- hist(X[,2], breaks = seq(min(X[,2]), max(X[,2]), l=(min(max(20,round(sqrt(nrow(X)))), 100)+1)),plot = F)
    B_dist2 <- cbind(h$mids, h$counts)
  }

  n1 <- nrow(B_dist1)
  n2 <- nrow(B_dist2)

  # use DP to find IC on boundary distributions
  aux_mx1 <- dpGMM:::dyn_pr_split_w_aux(B_dist1[,1],B_dist1[,2]) 
  aux_mx2 <- dpGMM:::dyn_pr_split_w_aux(B_dist2[,1],B_dist2[,2]) 

  tmp1 <- dpGMM:::dyn_pr_split_w(B_dist1[,1],B_dist1[,2], KS-1, aux_mx1)
  tmp2 <- dpGMM:::dyn_pr_split_w(B_dist2[,1],B_dist2[,2], KS-1, aux_mx2)

  opt_part1 <- tmp1[[2]]
  opt_part2 <- tmp2[[2]]

  part_cl1 <- c(1, opt_part1, n1+1)
  part_cl2 <- c(1, opt_part2, n2+1)
  pp_ini <- matrix(nrow=2, ncol=KS)
  mu_ini <- matrix(nrow=2, ncol=KS)
  sig_ini <- matrix(nrow=2, ncol=KS)

  for (k in 1:KS){
    invec <- B_dist1[(part_cl1[k]):(part_cl1[k+1]-1),1]
    yinwec <- B_dist1[(part_cl1[k]):(part_cl1[k+1]-1),2]
    wwec <- yinwec/sum(yinwec)
    pp_ini[1,k] <- sum(yinwec)/sum(B_dist1[,2])
    mu_ini[1,k] <- sum(invec*wwec)
    sig_ini[1,k] <- 0.5*(max(invec)-min(invec))

    invec <- B_dist2[(part_cl2[k]):(part_cl2[k+1]-1),1]
    yinwec <- B_dist2[(part_cl2[k]):(part_cl2[k+1]-1),2]
    wwec <- yinwec/sum(yinwec)
    pp_ini[2,k] <- sum(yinwec)/sum(B_dist2[,2])
    mu_ini[2,k] <- sum(invec*wwec)
    sig_ini[2,k] <- 0.5*(max(invec)-min(invec))
  }

  # aggregate 1D IC to 2D
  init$alpha <- matrix(nrow=1, ncol=KS*KS)
  init$center <- matrix(nrow=KS*KS, ncol=2)
  init$covar <- replicate(KS*KS, matrix(0,nrow=2,ncol=2), simplify = "array")
  k <- 1
  for (a in 1:KS){
    for (b in 1:KS){
      init$alpha[k] <- pp_ini[1,a] + pp_ini[2,b]
      init$center[k,] <- c(mu_ini[1,a], mu_ini[2,b])
      init$covar[,,k] <- rbind(c(sig_ini[1,a],0),c(0,sig_ini[2,b]))
      k <- k + 1
    }
  }
  init$alpha <- init$alpha/sum(init$alpha)
  init$KS <- length(init$alpha)

  # select KS initial conditions from KS*KS by maximizing log-likelihood
  fit_tmp <- calc_lLik2D(X,Y,init)$cum_pdf
  ind <- order(fit_tmp, decreasing=T)
  init$alpha <- init$alpha[ind[1:(2*KS)]]
  init$center <- init$center[ind[1:(2*KS)],]
  init$covar <- init$covar[,,ind[1:(2*KS)]]
  init$KS <- 2*KS

  while (init$KS > KS){
    fit_tmp <- numeric(init$KS)
    for (a in 1:init$KS){
      init_tmp <- init
      init_tmp$alpha <- init_tmp$alpha[-a]
      init_tmp$center <- init_tmp$center[-a,]
      init_tmp$covar <- init_tmp$covar[,,-a]
      init_tmp$KS <- init_tmp$KS-1
      fit_tmp[a] <- calc_lLik2D(X,Y,init_tmp)$lLik
    }
    ind <- which.min(abs(fit_tmp))
    init$alpha <- init$alpha[-ind]
    init$center <- init$center[-ind,]
    init$covar <- init$covar[,,-ind]
    init$KS <- init$KS - 1
  }
  init$alpha <- init$alpha/sum(init$alpha)

  return(init)
}


#' Diagonal initialization of 2D GMM
#'
#' Function for generating initial conditions of 2D GMM model by diagonal.
#'
#' @param X Matrix of 2D GMM data.
#' @param K Number of components.
#'
#' @return Initial values for EM.
#'
#' @keywords internal
#'
diag_init_2D <- function(X, KS){
  init <- list()
  init$alpha <- matrix(1, 1, KS)/KS #equal mixing proportions

  n_row <- X[nrow(X),1]
  n_col <- X[nrow(X),2]
  a <- cbind(seq(1,n_row, len = (KS*2)+1), seq(1,n_col, len = (KS*2)+1))
  colnames(a) <- colnames(X)
  init$center <- matrix(round(a[seq(2, nrow(a), 2),]), ncol=2)
  init$center <- apply(init$center, c(1,2), as.integer)

  init$covar <- replicate(KS, diag(apply(as.matrix(X), 2, sd)/KS), simplify = "array")
  init$KS <- KS
  return(init)
}

#' Probability distribution of 2D normal distribution
#'
#' Function for calculation PDF of 2D normal model.
#'
#' @param X matrix of 2D GMM data.
#' @param center centers of 2D distributions (means)
#' @param covar matrix of covariances
#'
#' @keywords internal
#'
norm_pdf_2D <- function(x, center, covar){
  den <- (6.283185307179585 * sqrt(det(covar)))
  x_centr <- sweep(as.matrix(x), 2, center, FUN = "-")
  tmp <- sweep(x_centr, 2, colSums(covar), FUN = "/") * x_centr
  num <- - 0.5 * rowSums(tmp)
  y <- exp(num)/den
  return(y)
}

#' 2D plot support
#'
#' Supporting function for 2D GMM ploting.
#'
#' @param X matrix of 2D GMM data.
#' @param Y y coordinates of histogram.
#'
#'
#' @keywords internal
#'
coords_to_img <- function(X,Y){

  x <- unique(round(X[,1]))
  y <- unique(round(X[,2]))
  img <- matrix(0,nrow=length(x),ncol=length(y))
  rownames(img) <- x
  colnames(img) <- y

  for (a in 1:nrow(X)){
    img[round(X[a,1]),round(X[a,2])] <- img[round(X[a,1]),round(X[a,2])] + Y[a]
  }
  img <- (max(img) - img)
  img <- 255 * img/max(img)

  return(img)
}

#' 2D plot support
#'
#' Transform image into coordinates data
#'
#' @param img image in 2D array.
#'
#'
#'
img_to_coords <- function(img){

  dim1 <- 1:nrow(img)
  dim2 <- 1:ncol(img)
  N = length(dim1)*length(dim2);

  data <- matrix(nrow=N, ncol=3);
  count <- 1
  for (a in 1:length(dim1)){
    for (b in 1:length(dim2)){
      data[count,1] <- dim1[a];
      data[count,2] <- dim2[b];
      data[count,3] <- img[a,b];
      count <- count + 1;
    }
  }

  data <- data.frame(data)
  colnames(data) <- c("Coordinates_1","Coordinates_2","Y")
  return(data)
}

#' Log-Likelihood for 2D Gaussian Mixture Model.
#'
#' Function calculate log-likelihood of 2D Gaussian distribution
#'
#' @param X matrix of 2D GMM data.
#' @param Y y coordinates of histogram.
#' @param gmm 2D GMM fit, output of \code{\link{gaussian_mixture_2D}} function
#'
#'@returns Function returns a \code{list} with: \describe{
#'  \item{lLik}{Log-Likelihood for 2D GMM}
#'  \item{cum_pdf}{Cumulative probability distibution}
#' }
#'
#'
#' @keywords internal
#'
calc_lLik2D <- function(X,Y,gmm){

  data <- rep(apply(X,1,as.list), times=Y)
  data <- lapply(data,function(x){unlist(x)})
  data <- do.call(rbind,data)

  #calculate density function
  f <- matrix(0, gmm$KS, nrow(data))
  for (a in 1:gmm$KS){
    f[a,] <- dpGMM:::norm_pdf_2D(data, gmm$center[a,], gmm$covar[,,a])
  }
  px <-  colSums(f * as.numeric(gmm$alpha))
  px[is.nan(px) | px==0] <- 1e-100

  res <- list()
  res[["lLik"]] <- sum(log(px)*Y)
  res[["cum_pdf"]] <- rowSums(f)
  return(res)
}

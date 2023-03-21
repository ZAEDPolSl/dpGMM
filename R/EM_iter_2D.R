#' Expectation–maximization algorithm for 2D data
#'
#' The function performs the EM algorithm to find the local maximum likelihood for the estimated Gaussian mixture parameters.
#'
#' @param X Matrix of 2D data to decompose by GMM.
#' @param Y Vector of counts, with the same length as "X".
#' Applies only to binned data (Y = NULL, by default).
#' @param init Vector of initial parameters for Gaussian components.
#' @param opts Parameters of run stored in \code{\link{GMM_2D_opts}} variable.
#'
#' @returns Function returns a \code{list} of GMM parameters for tested number of components: \describe{
#'  \item{alpha}{Weights (alpha) of each component.}
#'  \item{center}{Means of decomposition.}
#'  \item{covar}{Covariances of each component.}
#'  \item{KS}{Estimated number of components.}
#'  \item{logL}{Log-likelihood statistic for the estimated number of components.}
#'  \item{IC}{The value of the selected information criterion which was used to calculate the number of components.}
#' }
#'
#'
#' @examples
#' \dontrun{
#' }
#'
#' @seealso \code{\link{runGMM2D}}
#'
#' @export
EM_iter_2D <- function(X, Y, init, opts){

    Y <- as.vector(Y)
    N <- length(Y)
    change <- Inf
    L_new <- 0
    count <- 1
    KS <- init$KS
    bin_edge_sum <- sum(Y)
    alpha <- as.vector(init$alpha)
    center <- as.matrix(init$center)
    if (dim(center)[2]==KS) {center=t(center)}
    covar <- init$covar
    opts$SW <- diag(2)*opts$SW

    #MAIN LOOP
    while (change > opts$eps_change && count < opts$max_iter){
      L_old <- L_new

      #calculate density function
      f <- matrix(0, KS, N)
      for (a in 1:KS){
        f[a,] <- rGMMtest:::norm_pdf_2D(X, center[a,], covar[,,a]) ####### PACKAGE NAME!!!!!!!!!!!!!!!!
      }
      px <-  colSums(f * alpha)
      px[is.nan(px) | px==0] <- 1e-100

      for (a in 1:KS){
        pk <- (alpha[a]*f[a,]*Y)/px
        denom <- sum(pk)

        #maximization step
        alpha[a] <- denom/bin_edge_sum
        center[a,] <- colSums(apply(X, 2, "*", pk))/denom  #a[a,]

        if(opts$cov_type == "diag"){
          #diagonal covariance
          tmp <- sweep(as.matrix(X), 2, as.numeric(center[a,]), FUN = "-")^2
          covarnum <- colSums(apply(tmp, 2, "*", pk))
          sig2_min <- covarnum/denom
          covar[,,a] <- opts$SW + diag(sig2_min)
        }else if(opts$cov_type == "full"){
          #full covariance
          x_centr <- sweep(as.matrix(X), 2, as.numeric(center[a,]), FUN = "-")
          x_centr <- apply(x_centr, 2, "*", (sqrt(pk/denom)))
          covar[,,a] <- opts$SW + apply(x_centr, 2, "%*%", x_centr)
        }else if(opts$cov_type == "sphere"){
          #sphere covariance
          tmp <- sweep(as.matrix(X), 2, as.numeric(center[a,]), FUN = "-")^2
          covarnum <- colSums(apply(tmp, 2, "*", pk))
          sig2_tmp <- mean(covarnum/denom)
          covar[,,a] <- opts$SW + diag(2)*sig2_tmp
        }

        #check if singularity appeared
        if (det(covar[,,a]) <= 0.1){
          covar[,,a] <- covar[,,a]*diag(2)
          # alpha[a] <- 0
        }

        #check if shape is proper
        eig_tmp <- rev(eigen(covar[,,a])[[1]])
        if (max(eig_tmp/min(eig_tmp))>opts$max_var_ratio){
          if (eig_tmp[2]>eig_tmp[1]){
            if (opts$cov_type == "diag"){
              covar[1,1,a] <- opts$max_var_ratio * covar[2,2,a]
            }else{
              e <- eigen(covar[,,a])
              V <- apply(e$vectors, 1, rev); D <- diag(rev(e$values))
              D[1,1] <- opts$max_var_ratio * D[2,2]
              covar[,,a] <- V %*% D %*% t(V)
            }
          }
        }
      }
      L_new <- sum(log(px)*Y)
      change <- 100 * abs((L_new - L_old)/L_old)

      count <- count + 1
    }



    #RETURN RESULTS
    gmm <- list()
    gmm$alpha <- alpha
    gmm$center <- center
    gmm$covar <- covar
    gmm$KS <- KS
    gmm$logL <- L_new

    if(opts$IC == "ICL-BIC"){
      if(exists("pk")){
        pk[is.nan(pk) | pk==0] = 5e-324
        EN <- -sum(sum(pk*log(pk)))
      }else{
        EN <- Inf
      }
    }

    switch(opts$IC,
           "BIC" = crit_val <- abs(2*L_new) + (7*KS-1)*log(bin_edge_sum),
           "AIC" = crit_val <- abs(2*L_new) + 2*(7*KS-1),
           "AICc" = crit_val <- abs(2*L_new) + 2*(7*KS-1) * (bin_edge_sum/(bin_edge_sum - (7*KS-1) -1)),
           "ICL-BIC" = crit_val <- abs(2*L_new) + 2*EN + (7*KS-1)*log(bin_edge_sum)
    )

    gmm$IC <- crit_val
    return(gmm)
  }

gaussian_mixture_2D <- function(X, Y=NULL, opts){

  if (dim(X)[1] == 2){
    X <- t(X)
  }
  if (is.null(Y)){
    Y <- rep(1,nrow(X))
  }

  # initialize variables
  bin_edge_sum <- sum(Y)
  IC <- matrix(NaN, opts$init_nb, opts$KS)
  logL <- IC
  D <- matrix(NaN, 1, opts$KS)
  gmm <- list()

  if(opts$init_con != "rand"){
    opts$init_nb <- 1
  }

  # decomposition for 1 component
  gmm[[1]] <- EM_iter_2D(X, Y, rand_init_2D(X, 1), opts)
  logL[,1] <- matrix(gmm[[1]]$logL, opts$init_nb, 1)
  IC <- gmm[[1]]$IC
  IC[,1] <- matrix(IC, opts$init_nb, 1)
  # gmm[[1]]$IC <- IC

  # decomposition for >2 components
  stop <- 1
  k <- 2

  while(stop && k < opts$KS){
    # find initial conditions
    gmm_tmp <- list()
    IC_tmp <- matrix(NaN, opts$init_nb, 1)
    logL_tmp <- IC_tmp

    for(a in 1:opts$init_nb){
      # perform decomposition
      if(opts$init_con == "rand"){
        gmm_tmp[[a]] <- EM_iter_2D(X, Y, rand_init_2D(X, k), opts)
      } else if(opts$init_con == "diag"){
        gmm_tmp[[a]] <- EM_iter_2D(X, Y, diag_init_2D(X, k), opts)
      } else if(opts$init_con == "DP"){
        gmm_tmp[[a]] <- EM_iter_2D(X, Y, DP_init_2D(X, Y, k), opts)
      }
      logL_tmp[a,1] <- gmm_tmp[[a]]$logL
      IC_tmp[a,1] <- gmm_tmp[[a]]$IC
    }

    gmm[[k]] <- gmm_tmp
    logL[,k] <- logL_tmp
    IC[,k] <- IC_tmp

    # check convergence
    D[k] <- -2*median(logL[,k-1]) + 2*median(logL[,k])
    if(pchisq(D[k], 7, lower.tail = F) > opts$D_thr){stop <- 0}

    k <- k+1
  }

  cmp_nb <- which(apply(IC, 2, median) == min(apply(IC, 2, median), na.rm = T))

  if(cmp_nb >1){
    tmp <- min(abs(IC[,cmp_nb] - median(IC[,cmp_nb])))
    ind <- which(tmp == min(tmp))

    gmm_out <- gmm[[cmp_nb]][[ind]]
    gmm_out$IC <- IC[ind, cmp_nb]
  } else{
    ind <- 1

    gmm_out <- gmm[[cmp_nb]]
    gmm_out$IC <- IC[ind, cmp_nb]
  }

  return(gmm_out)
}

gaussian_mixture_2D <- function(data, counts, opts){
  
  # initialize variables
  bin_edge_sum <- sum(counts)
  BIC <- matrix(NaN, opts$init_nb, opts$KS)
  logL <- BIC
  D <- matrix(NaN, 1, opts$KS)
  gmm <- list()
  
  # decomposition for 1 component
  if(opts$init_con == "rand"){
    gmm[[1]] <- EM_iter_2D(data, counts, rand_init_2D(data, 1), opts)
  } else if(opts$init_con == "diag"){
    gmm[[1]] <- EM_iter_2D(data, counts, diag_init_2D(data, 1), opts)
  }
  logL[,1] <- matrix(gmm[[1]]$logL, opts$init_nb, 1)
  # IC <- abs(2*gmm[[1]]$logL) + (7*gmm[[1]]$KS-1)*log(bin_edge_sum) *opts$res
  IC <- gmm[[1]]$IC
  BIC[,1] <- matrix(IC, opts$init_nb, 1)
  # gmm[[1]]$BIC <- IC
  
  # decomposition for >2 components
  stop <- 1
  k <- 2
  
  while(stop && k < opts$KS){
    # find initial conditions randomly
    gmm_tmp <- list()
    BIC_tmp <- matrix(NaN, opts$init_nb, 1)
    logL_tmp <- BIC_tmp
    
    for(a in 1:opts$init_nb){
      # perform decomposition
      if(opts$init_con == "rand"){
        gmm_tmp[[a]] <- EM_iter_2D(data, counts, rand_init_2D(data, k), opts)
      } else if(opts$init_con == "diag"){
        gmm_tmp[[a]] <- EM_iter_2D(data, counts, diag_init_2D(data, k), opts)
      }
      logL_tmp[a,1] <- gmm_tmp[[a]]$logL
      # BIC_tmp[a,1] <- abs(2*gmm_tmp[[a]]$logL) + (7*gmm_tmp[[a]]$KS-1)*log(bin_edge_sum) *opts$res
      BIC_tmp[a,1] <- gmm_tmp[[a]]$IC
    }
    
    gmm[[k]] <- gmm_tmp
    logL[,k] <- logL_tmp
    BIC[,k] <- BIC_tmp
    
    # check convergence
    D[k] <- -2*median(logL[,k-1]) + 2*median(logL[,k])
    if(pchisq(D[k], 7, lower.tail = F) > opts$D_thr){stop <- 0}
    
    k <- k+1
  }
  
  cmp_nb <- which(apply(BIC, 2, median) == min(apply(BIC, 2, median), na.rm = T))
  
  if(cmp_nb >1){
    tmp <- min(abs(BIC[,cmp_nb] - median(BIC[,cmp_nb])))
    ind <- which(tmp == min(tmp))
    
    gmm_out <- gmm[[cmp_nb]][[ind]]
    gmm_out$IC <- BIC[ind, cmp_nb]
  } else{
    ind <- 1
    
    gmm_out <- gmm[[cmp_nb]]
    gmm_out$IC <- BIC[ind, cmp_nb]
  }
  
  return(gmm_out)
}
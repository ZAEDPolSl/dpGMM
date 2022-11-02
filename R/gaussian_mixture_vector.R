gaussian_mixture_vector <- function(data, KS, Y = NULL, change = Inf, max_iter = 5000, SW=NULL, IC = "BIC", quick_stop = TRUE, signi = 0.05){

  if (min(dim(as.matrix(data))) !=1){
    stop("data must be 1D signal.")
  }

  if(is.null(Y)){Y<-matrix(1, 1, length(data))}
  bin_edge_sum <- sum(Y)

  IC_list<-c("AIC","AICc","BIC", "ICL-BIC", "LR")

  N <- length(data)
  crit_vector <- matrix(NaN, KS, 1)
  logL <- crit_vector
  D <- matrix(NaN, 1, KS)
  alpha <- list()
  mu <- list()
  sigma <- list()

  #histogram of input data (for drawing and IC).. soon
  h  <- hist(data, breaks = seq(min(data), max(data), l=(max(min(20,round(sqrt(N))), 100)+1)),plot = F)
  y <- h$counts
  x <- h$mids
  #decomposition for 1 component
  rcpt <- EM_iter(data, 1, mean(data), sd(data), N, Y, change, max_iter, SW, IC)

  alpha[[1]] <- rcpt[[1]]
  mu[[1]] <- rcpt[[2]]
  sigma[[1]] <- rcpt[[3]]
  logL[1] <- rcpt[[4]]
  if(IC != "LR"){crit_vector[1] <- rcpt[[5]]}

  #decomposition for >2 components
  stop <- 1
  k <- 2
  Nb <- length(x)
  aux_mx <- rGMMtest:::dyn_pr_split_w_aux(x,y) #CORRECT TO FINAL NAME OF PACKAGE !!!!!!!!!!!!!!!!!

  while (stop && k < KS){
    tmp <- rGMMtest:::dyn_pr_split_w(x, y, k-1, aux_mx) #CORRECT TO FINAL NAME OF PACKAGE!!!!!!
    opt_part <- tmp[[2]]

    part_cl <- c(1, opt_part, Nb+1)
    pp_ini <- matrix(0, 1, k)
    mu_ini <- matrix(0, 1, k)
    sig_ini <- matrix(0, 1, k)

    for (kkps in 1:k){
      invec <- x[(part_cl[kkps]):(part_cl[kkps+1]-1)]
      yinwec <- y[(part_cl[kkps]):(part_cl[kkps+1]-1)]
      wwec <- yinwec/sum(yinwec)
      pp_ini[kkps] <- sum(yinwec)/sum(y)
      mu_ini[kkps] <- sum(invec*wwec)
      sig_ini[kkps] <- 0.5*(max(invec)-min(invec))
    }

    #perform decomposition
    rcpt1<- EM_iter(data, pp_ini, mu_ini, sig_ini, N, Y, change, max_iter, SW, IC)
    alpha[[k]] <- rcpt1[[1]]
    mu[[k]] <- rcpt1[[2]]
    sigma[[k]] <- rcpt1[[3]]
    logL[k] <- rcpt1[[4]]

    if(IC != "LR"){crit_vector[k] <- rcpt1[[5]]}

    if(quick_stop | IC == "LR"){D[k] <- -2*logL[k-1] + 2*logL[k]}

    if(quick_stop){
      if ((1 - pchisq(D[k],3)) > signi){stop <- 0}
      }

    if(IC == "LR"){crit_vector[k] <- D[k]}
    k <- k+1
  }

  if(IC == "LR"){
    LR_crit <- qchisq(1-signi, 3) # crit. val
    D <- D[c(1, which(D > LR_crit))]
    cmp_nb <- which.min(D)
    crit_est <- (1 - pchisq(D[cmp_nb], 3))
  }else{
    cmp_nb <- which(crit_vector == min(crit_vector, na.rm = T))
    crit_est <- crit_vector[cmp_nb]
    }

  pp_est <- alpha[[cmp_nb]]
  mu_est <- mu[[cmp_nb]]
  sig_est <- sigma[[cmp_nb]]
  logL_est <- logL[cmp_nb]

  GModel <- data.frame(mu =  mu_est, sigma = sig_est, alpha = pp_est)
  GModel <- GModel[order(GModel$mu),]
  res <- list(model=GModel, IC=crit_est, logL=logL_est, KS=length(pp_est))


  return(res)
}

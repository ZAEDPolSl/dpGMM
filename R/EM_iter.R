EM_iter <- function(x, alpha, mu, sig, N, Y = NULL, change = Inf, max_iter = 5000, SW = NULL, IC = "BIC"){
  
  if(is.null(Y)){Y<-matrix(1, 1, length(x))}
  bin_edge_sum <- sum(Y)
  
  alpha <- c(alpha)
  mu <- c(mu)
  sig <- c(sig)
  
  x <- sort(x)
  sig2 <- sig^2
  count <- 1
  eps_change <- 1e-7
  KS <- length(alpha)
  
  if (is.null(SW)){
    SW <- ((max(x)-min(x))/(5*KS))^2 #minimum variance
  }
  
  while (change > eps_change && count < max_iter){
    old_alpha <- alpha
    old_sig2 <- sig2
    #cat(alpha, "  ",mu, "   ",sig," \n")
    f <- matrix(0, KS, N)
    sig <- sqrt(sig2)
    
    tmp <- sig <= 0
    
    if (sum(tmp) > 1){
      change <- 0
    } else {
      if (sum(tmp)){
        sig[tmp] <- SW  
      }
      
      for(a in 1:KS){
        f[a,] <- norm_pdf(x, mu[a], sig[a])
      }
      
      px <- matrix(0, KS, N)
      for(i in 1:KS){
        px[i,] <- alpha[i] * f[i,]
      }
      px <- colSums(px)
      #px <- alpha*f
      px[is.nan(px) | px==0] = 5e-324
      
      for (a in 1:KS){
        pk <- ((alpha[a]*f[a,])*Y)/px
        denom <- sum(pk)
        mu[a] <- sum((pk*x)/denom)
        sig2num <- sum(pk*((x-mu[a])^2))
        sig2[a] <- max(SW, sig2num/denom)
        alpha[a] <- denom/bin_edge_sum
      }
      change <- sum(abs(alpha-old_alpha)) + sum(((abs(sig2 - old_sig2))/sig2))/(length(alpha))     
    }
    
    count <- count + 1
  }
  #RETURN RESULTS
  if(sum(tmp)>1){
    logL <- -Inf
  } else {
    logL <- sum(log(px)*Y)
  }
  mu_est <- sort(mu)
  ind <- order(mu)
  sig_est <- sqrt(sig2[ind])
  pp_est <- alpha[ind]
  
  #KS <- length(alpha) #będzie odpowiadać k
  #calculating the information criterion
  
  if(IC == "ICL-BIC"){
    if(change != 0){
      pk[is.nan(pk) | pk==0] = 5e-324
      EN <- -sum(sum(pk*log(pk)))
    } else {
      EN <- Inf
    }
  }
  
  switch(IC,
         "BIC" = crit_val <- -2*logL + (3*KS-1)*log(bin_edge_sum),
         "AIC" = crit_val <- -2*logL + 2*(3*KS-1),
         "AICc" = crit_val <- -2*logL + 2*(3*KS-1)*(bin_edge_sum/(bin_edge_sum-(3*KS-1)-1)),
         "ICL-BIC" = crit_val <- -2*logL + 2*EN + (3*KS-1)*log(bin_edge_sum),
         "LR" = crit_val <- NULL
  )
  # switch(IC,
  #        "BIC" = crit_val <- -2*logL + 2*log(bin_edge_sum),
  #        "AIC" = crit_val <- -2*logL + 2*(3*1-1),
  #        "AICc" = crit_val<- -2*logL + 2*(3*1-1)*(bin_edge_sum/(bin_edge_sum-(3*1-1)-1)), #dla k = 1
  # )
  
  return(list(pp_est, mu_est, sig_est, logL, crit_val))
}

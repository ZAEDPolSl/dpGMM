#' Expectation–maximization algorithm for 1D data
#'
#' The function performs the EM algorithm to find the local maximum likelihood for the estimated Gaussian mixture parameters.
#'
#' @param X Vector of 1D data for GMM decomposition.
#' @param alpha Vector containing the weights (alpha) for each component in the statistical model.
#' @param mu Vector containing the means (mu) for each component in the statistical model.
#' @param sig Vector containing the standard deviation (sigma) for each component in the statistical model.
#' @param Y Vector of counts, with the same length as "X".
#' Applies only to binned data (Y = NULL, by default).
#' @param opts Parameters of run saved in \code{\link{GMM_1D_opts}} variable.
#'
#' @return Returns a \code{list} of GMM parameter values that correspond to the local extremes for each component.\describe{
#'  \item{alpha}{Vector of optimal alpha (weights) values.}
#'  \item{mu}{Vector of optimal mu (means) values.}
#'  \item{sigma}{Vector of optimal sigma (standard devations) values.}
#'  \item{logLik}{Log-likelihood statistic for the estimated number of components.}
#'  \item{crit}{Value of the selected information criterion in local extreme of likelihood function.}
#' }
#'
#' @importFrom stats dnorm
#'
#' @seealso \code{\link{runGMM}} and \code{\link{gaussian_mixture_vector}}
#'
#' @export
# EM_iter <- function(X, alpha, mu, sig, Y = NULL, eps_change = 1e-7, max_iter = 50000, SW = 0.01, IC = "BIC"){
EM_iter <- function(X, alpha, mu, sig, Y = NULL, opts = NULL){

  if(is.null(opts)){opts <- dpGMM::GMM_1D_opts}

  if(is.null(Y)){Y <- matrix(1, 1, length(X))}
  bin_edge_sum <- sum(Y)

  alpha <- c(alpha)
  mu <- c(mu)
  sig <- c(sig)

  N <- length(X)
  X <- sort(X)
  sig2 <- sig^2
  count <- 1
  change <- Inf
  KS <- length(alpha)

  opts$SW <- (((max(X) - min(X)) * opts$SW)/KS)^2 #minimum variance


  while (change > opts$eps_change && count < opts$max_iter){
    old_alpha <- alpha
    old_sig2 <- sig2
    f <- matrix(0, KS, N)
    sig <- sqrt(sig2)

    tmp <- sig <= 0

    if (sum(tmp) > 1){
      change <- 0
    } else {
      if (sum(tmp)){
        sig[tmp] <- opts$SW
      }

      for(a in 1:KS){
        f[a,] <- dnorm(X, mu[a], sig[a])
      }

      px <- matrix(0, KS, N)
      for(i in 1:KS){
        px[i,] <- alpha[i] * f[i,]
      }
      px <- colSums(px)
      px[is.nan(px) | px==0] = 5e-324

      for (a in 1:KS){
        pk <- ((alpha[a]*f[a,])*Y)/px
        denom <- sum(pk)
        mu[a] <- sum((pk*X)/denom)
        sig2num <- sum(pk*((X-mu[a])^2))
        sig2[a] <- max(opts$SW, sig2num/denom)
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
    logL <- sum(log(px) * Y)
  }
  mu_est <- sort(mu)
  ind <- order(mu)
  sig_est <- sqrt(sig2[ind])
  pp_est <- alpha[ind]

  #calculating the information criterion
  if(opts$IC == "ICL-BIC"){
    if(change != 0){
      pk[is.nan(pk) | pk == 0] = 5e-324
      EN <- -sum(sum(pk * log(pk)))
    } else {
      EN <- Inf
    }
  }

  switch(opts$IC,
         "BIC" = crit_val <- -2 * logL + (3 * KS-1) * log(bin_edge_sum),
         "AIC" = crit_val <- -2 * logL + 2 * (3 * KS-1),
         "AICc" = crit_val <- -2 * logL + 2 * (3 * KS-1) * (bin_edge_sum/(bin_edge_sum-(3 * KS-1)-1)),
         "ICL-BIC" = crit_val <- -2 * logL + 2 * EN + (3 * KS-1) * log(bin_edge_sum),
         "LR" = crit_val <- NULL
  )

  return(list(alpha = pp_est, mu = mu_est, sigma = sig_est, logL = logL, crit = crit_val))
}

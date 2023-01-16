#' Gaussian mixture decomposition for a vector of data
#'
#'Function to choose the optimal number of components of a mixture normal distributions, minimising the value of the information criterion.
#'
#' @param X Vector of data to decompose by GMM.
#' @param KS Maximum number of components to test.
#' @param Y Vector of counts, should be the same length as "data".
#' Applies only to binned data therefore the default is Y = NULL.
#' @param change Stop of EM criterion (default Inf), value compared to following formula:
#' \deqn{\sum{(|\alpha - \alpha_{old})|} + \frac{\sum{(\frac{|\sigma^2 - \sigma^2_{old}|}{\sigma^2})}}{length(\alpha)}}
#' @param max_iter Maximum number of iterations of EM algorithm. By default it is \code{max_iter = 5000}
#' @param SW Minimum standard deviation of component (default 0.1). Values from 0 to 1 are treated according to their values, for SW >1 the following formula is applied:
#' \deqn{\frac{range(x)}{(SW*no.of.components))^2}}.
#' @param IC Information Criterion to select best number of components.
#' Possible "AIC","AICc", "BIC" (default), "ICL-BIC" or "LR".
#' @param quick_stop Logical value. Determines to stop the EM algorithm when adding another component is no longer significant according to the Likelihood Ratio Test. Used to speed up the function (Default is TRUE).
#' @param signi Significance level for Likelihood Ratio Test. By default is 0.05.
#'
#' @returns Function returns a \code{list} of GMM parameters for the optimal number of components: \describe{
#'  \item{model}{A \code{data.frame} of model component parameters - mean values (mu), standard deviations (sigma)
#'  and weights (alpha) for each component. Output of \code{EM_iter}}
#'  \item{IC}{Value of the selected information criterion which was used to calculate the optimal number of components}
#'  \item{logL}{Log-likelihood value for the optimal number of components}
#'  \item{KS}{Optimal number of components}
#' }
#'
#' @importFrom stats pchisq qchisq
#' @importFrom graphics hist
#'
#' @examples
#' \dontrun{
#' data <- generate_norm1D(1000, alpha=c(0.2,0.4,0.4), mu=c(-15,0,15), sigma=c(1,2,3))
#' exp <- gaussian_mixture_vector(data$Dist, KS = 10, IC = "AIC", quick_stop = FALSE)
#' }
#'
#' @seealso \code{\link{runGMM}} and \code{\link{EM_iter}}
#'
#' @export
gaussian_mixture_vector <- function(X, KS, Y = NULL, change = Inf, max_iter = 5000, SW=0.1, IC = "BIC", quick_stop = TRUE, signi = 0.05){

  if (min(dim(as.matrix(X))) !=1){
    stop("data must be 1D signal.")
  }

  if(is.null(Y)){Y<-matrix(1, 1, length(X))}
  bin_edge_sum <- sum(Y)

  IC_list<-c("AIC","AICc","BIC", "ICL-BIC", "LR")

  N <- length(X)
  crit_vector <- matrix(NaN, KS, 1)
  logL <- crit_vector
  D <- matrix(NaN, 1, KS)
  alpha <- list()
  mu <- list()
  sigma <- list()

  #histogram of input data (for drawing and IC).. soon
  h  <- hist(X, breaks = seq(min(X), max(X), l=(min(max(20,round(sqrt(N))), 100)+1)),plot = F)
  y <- h$counts
  x <- h$mids
  #decomposition for 1 component
  rcpt <- EM_iter(X, 1, mean(X), sd(X), N, Y, change, max_iter, SW, IC)

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
    rcpt1<- EM_iter(X, pp_ini, mu_ini, sig_ini, N, Y, change, max_iter, SW, IC)
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
    ind_tmp <- which(D > LR_crit)
    if (length(ind_tmp) > 0){
      cmp_nb <- which.min(D)
      crit_est <- (1 - pchisq(D[cmp_nb], 3))
    } else {
      cmp_nb <- 1
      crit_est <- NA
    }
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

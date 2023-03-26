#' Gaussian mixture decomposition for 1D data
#'
#' Function to estimate number of components of a mixture normal distributions, minimizing the value of the information criterion.
#'
#' @param X Vector of 1D data for GMM decomposition.
#' @param KS Maximum number of components of the model.
#' @param opts Parameters of run saved in \code{\link{GMM_1D_opts}} variable.
#' @param Y Vector of counts, with the same length as "X".
#' Applies only to binned data (Y = NULL, by default).
#'
#'
#' @returns Function returns a \code{list} of GMM parameters for the estimated number of components: \describe{
#'  \item{model}{A \code{list} of model component parameters - mean values (mu), standard deviations (sigma)
#'  and weights (alpha) for each component.}
#'  \item{IC}{The value of the selected information criterion which was used to calculate the number of components.}
#'  \item{logLik}{Log-likelihood statistic for the estimated number of components.}
#'  \item{KS}{Estimaged number of model components.}
#' }
#'
#' @importFrom stats pchisq qchisq
#' @importFrom graphics hist
#'
#' @examples
#' \dontrun{
#' data <- generate_norm1D(1000, alpha = c(0.2,0.4,0.4), mu = c(-15,0,15), sigma = c(1,2,3))
#'
#' custom.settings <- GMM_1D_opts
#' custom.settings$IC <- "AIC"
#'
#' exp <- gaussian_mixture_vector(data$Dist, KS = 10, opts = custom.settings)
#' }
#'
#' @seealso \code{\link{runGMM}} and \code{\link{generate_norm1D}}
#'
#' @export
# gaussian_mixture_vector <- function(X, KS, Y = NULL, fixed = FALSE , eps_change = 1e-7, max_iter = 50000, SW = 0.01, IC = "BIC", quick_stop = TRUE, signi = 0.05){
gaussian_mixture_vector <- function(X, KS, opts = GMM_1D_opts, Y = NULL){

  if (min(dim(as.matrix(X))) != 1){
    stop("data must be 1D signal.")
  }

  if(is.null(Y)){Y <- matrix(1, 1, length(X))}
  bin_edge_sum <- sum(Y)

  IC_list <- c("AIC", "AICc", "BIC", "ICL-BIC", "LR")

  N <- length(X)
  crit_vector <- matrix(NaN, KS, 1)
  logL <- crit_vector
  D <- matrix(NaN, 1, KS)
  alpha <- list()
  mu <- list()
  sigma <- list()

  #histogram of input data (for drawing and IC)
  h  <- hist(X, breaks = seq(min(X), max(X), l = (min(max(20, round(sqrt(N))), 100)  +1)), plot = F)
  y <- h$counts
  x <- h$mids
  #decomposition for 1 component
  # rcpt <- EM_iter(X, 1, mean(X), sd(X), Y, eps_change, max_iter, SW, IC)
  rcpt <- EM_iter(X, 1, mean(X), sd(X), Y, opts)

  alpha[[1]] <- rcpt[[1]]
  mu[[1]] <- rcpt[[2]]
  sigma[[1]] <- rcpt[[3]]
  logL[1] <- rcpt[[4]]
  if(opts$IC != "LR"){crit_vector[1] <- rcpt[[5]]}

  Nb <- length(x)
  aux_mx <- rGMMtest:::dyn_pr_split_w_aux(x, y) #CORRECT TO FINAL NAME OF PACKAGE !!!!!!!!!!!!!!!!!

  #decomposition for fixed KS number
  if (opts$fixed){
      k <- KS
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
        mu_ini[kkps] <- sum(invec * wwec)
        sig_ini[kkps] <- 0.5 * (max(invec)-min(invec))
      }

      # rcpt1<- EM_iter(X, pp_ini, mu_ini, sig_ini, Y, eps_change, max_iter, SW, IC)
      rcpt1<- EM_iter(X, pp_ini, mu_ini, sig_ini, Y, opts)

      pp_est <- rcpt1[[1]]
      mu_est <- rcpt1[[2]]
      sig_est <- rcpt1[[3]]
      logL_est <- rcpt1[[4]]
      crit_est <- rcpt1[[5]]


  } else{ #decomposition for range 1:KS

      #decomposition for >2 components
      stop <- 1
      k <- 2

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
        # rcpt1<- EM_iter(X, pp_ini, mu_ini, sig_ini, Y, eps_change, max_iter, SW, IC)
        rcpt1<- EM_iter(X, pp_ini, mu_ini, sig_ini, Y, opts)
        alpha[[k]] <- rcpt1[[1]]
        mu[[k]] <- rcpt1[[2]]
        sigma[[k]] <- rcpt1[[3]]
        logL[k] <- rcpt1[[4]]

        if(opts$IC != "LR"){crit_vector[k] <- rcpt1[[5]]}

        if(opts$quick_stop | opts$IC == "LR"){D[k] <- -2 * logL[k-1] + 2 * logL[k]}

        if(opts$quick_stop){
          if ((1 - pchisq(D[k], 3)) > opts$signi){stop <- 0}
          }

        if(opts$IC == "LR"){crit_vector[k] <- D[k]}
        k <- k+1
      }

      if(opts$IC == "LR"){
        LR_crit <- qchisq(1-opts$signi, 3) # crit. val
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

  }


  GModel <- data.frame(mu =  mu_est, sigma = sig_est, alpha = pp_est)
  GModel <- GModel[order(GModel$mu),]
  res <- list(model = GModel, IC = crit_est, logL = logL_est, KS = length(pp_est))


  return(res)
}

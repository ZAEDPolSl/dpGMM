#' Gaussian mixture decomposition for 2D data
#'
#' Function to choose the optimal number of components of a 2D mixture normal distributions, minimizing the value of the information criterion.
#'
#' @param X Matrix of 2D data to decompose by GMM.
#' @param Y Vector of counts, with the same length as "X".
#' Applies only to binned data (Y = NULL, by default).
#' @param opts Parameters of run saved in \code{\link{GMM_2D_opts}} variable.
#'
#' @returns Function returns a \code{list} of GMM parameters for the optimal number of components: \describe{
#'  \item{alpha}{Weights (alpha) of each component.}
#'  \item{center}{Means of decomposition.}
#'  \item{covar}{Covariances of each component.}
#'  \item{KS}{Estimated number of components.}
#'  \item{logL}{Log-likelihood statistic for the estimated number of components.}
#'  \item{IC}{The value of the selected information criterion which was used to calculate the number of components.}
#'  \item{cls}{Assigment of point to the clusters.}
#' }
#'
#' @importFrom stats pchisq
#'
#' @examples
#' \dontrun{
#' data(example2D)
#' custom.settings <- GMM_2D_opts
#' exp <- gaussian_mixture_2D(example2D[,1:2], example2D[,3], opts = custom.settings)
#' }
#'
#' @seealso \code{\link{runGMM2D}}, \code{\link{GMM_2D_opts}}
#'
#' @export
gaussian_mixture_2D <- function(X, Y = NULL, opts = NULL){

  if(is.null(opts)){opts <- rGMMtest::GMM_2D_opts}

  if (dim(X)[1] == 2){
    X <- t(X)
  }
  if (is.null(Y)){
    Y <- rep(1, nrow(X))
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
  gmm[[1]] <- EM_iter_2D(X, Y, rGMMtest:::rand_init_2D(X, 1), opts) ########## PACKAGE NAME
  logL[,1] <- matrix(gmm[[1]]$logL, opts$init_nb, 1)
  IC <- gmm[[1]]$IC
  IC <- matrix(NA, opts$init_nb, opts$KS+1)
  IC[,1] <- gmm[[1]]$IC



  if (opts$fixed){ # decomposition for fixed KS
    k <- opts$KS
    gmm_tmp <- list()
    IC_tmp <- matrix(NaN, opts$init_nb, 1)
    logL_tmp <- IC_tmp

    for(a in 1:opts$init_nb){
      # perform decomposition
      if(opts$init_con == "rand"){
        gmm_tmp[[a]] <- EM_iter_2D(X, Y, rGMMtest:::rand_init_2D(X, k), opts)  ########## PACKAGE NAME
      } else if(opts$init_con == "diag"){
        gmm_tmp[[a]] <- EM_iter_2D(X, Y, rGMMtest:::diag_init_2D(X, k), opts)  ########## PACKAGE NAME
      } else if(opts$init_con == "DP"){
        gmm_tmp[[a]] <- EM_iter_2D(X, Y, rGMMtest:::DP_init_2D(X, Y, k), opts)  ########## PACKAGE NAME
      }
      logL_tmp[a,1] <- gmm_tmp[[a]]$logL
      IC_tmp[a,1] <- gmm_tmp[[a]]$IC
    }

    gmm[[k]] <- gmm_tmp
    logL[,k] <- logL_tmp
    IC[,k] <- IC_tmp

  } else{ # decomposition for >2 components

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
        gmm_tmp[[a]] <- EM_iter_2D(X, Y, rGMMtest:::rand_init_2D(X, k), opts) ########## PACKAGE NAME
      } else if(opts$init_con == "diag"){
        gmm_tmp[[a]] <- EM_iter_2D(X, Y, rGMMtest:::diag_init_2D(X, k), opts) ########## PACKAGE NAME
      } else if(opts$init_con == "DP"){
        gmm_tmp[[a]] <- EM_iter_2D(X, Y, rGMMtest:::DP_init_2D(X, Y, k), opts) ########## PACKAGE NAME
      }
      logL_tmp[a,1] <- gmm_tmp[[a]]$logL
      IC_tmp[a,1] <- gmm_tmp[[a]]$IC
    }

    gmm[[k]] <- gmm_tmp
    logL[,k] <- logL_tmp
    IC[,k] <- IC_tmp


    # check convergence
    if(opts$quick_stop){
      D[k] <- -2 * median(logL[,k-1]) + 2 * median(logL[,k])
      if(stats::pchisq(D[k], 7, lower.tail = F) > opts$signi){stop <- 0}
    }

    k <- k+1
  }
  }


  cmp_nb <- which(apply(IC, 2, median) == min(apply(IC, 2, median), na.rm = T))

  if(cmp_nb > 1){
    tmp <- min(abs(IC[,cmp_nb] - median(IC[,cmp_nb])))
    ind <- which(tmp == min(tmp))

    gmm_out <- gmm[[cmp_nb]][[ind]]
    gmm_out$IC <- IC[ind, cmp_nb]
  } else{
    ind <- 1

    gmm_out <- gmm[[cmp_nb]]
    gmm_out$IC <- IC[ind, cmp_nb]
  }



  if (gmm_out$KS > 1){
  cls <- find_class_2D(X,gmm_out)
  } else {
  cls <- rep(1, dim(X)[1])
  }
  gmm_out$cls <- cls
  return(gmm_out)
}

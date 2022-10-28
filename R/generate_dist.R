#' Generation of GMM data with high precision
#'
#' Function to generate GMM distribution data with their PDF and cumulative results.
#'
#' @param data Vector of data
#' @param GModel List of GModel parameters i.e GModel$model$alpha, GModel$model$mu, GModel$model$sigma. Given as output of function gaussian_mixture_vector
#' @param precision Precision of point linespacing
#'
#' @importFrom stats dnorm
#' @importFrom pracma linspace
#'
#' @returns List with following elements::\describe{
#'    \item{x}{Numeric vector with equaliy spread data of given precison}
#'    \item{dist}{Matrix with PDF of each GMM component and cumulative distribution}
#' }
#'
#'
#' @examples
#' data<-generate_norm1D(1000, alpha=c(0.2,0.4,0.4), mu=c(-15,0,15), sigma=c(1,2,3))
#' GModel<-list(model=list(alpha=c(0.2,0.4,0.4), mu=c(-15,0,15), sigma=c(1,2,3)))
#' res<-generate_dist(data,GModel,1e4)
#'
#' @seealso \code{\link{runGMM}}
#'
#' @export
generate_dist<-function(data, GModel, precision){
x_temp = linspace(min(data),max(data),precision)
f_temp = matrix(0, precision, nrow(GModel))
for(k in 1:nrow(GModel)){
  f_temp[,k] = GModel$alpha[k] * dnorm(x_temp, mean = GModel$mu[k], sd =GModel$sigma[k])
}

f_temp <- as.data.frame(f_temp)
f_temp$main <- rowSums(f_temp)

return(list(x=x_temp,dist=f_temp))
}

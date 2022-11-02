#' Calculation of GMM thresholds
#'
#' Function to calculate cutoffs between each component of mixture normal distributions.
#'
#' @param GModel \code{data.frame} of GMM parameters i.e GModel$alpha, GModel$mu, GModel$sigma (correct \code{colnames} are obligatory)
#' @param input output of \code{generate_dist} function. Its necessary only if arithmetical approach fails in threshold estimation and \code{find_thr_by_dist} function is called.
#' It is a list with following elements:\describe{
#'    \item{x}{Numeric vector with equaliy spread data of given precison}
#'    \item{dist}{Matrix with PDF of each GMM component and cumulative distribution}
#' }
#'
#' @returns Return a vector of thresholds.
#'
#' @examples
#' data(example)
#' GModel<-data.frame(alpha=c(0.45,0.5,0.05),
#'                   mu=c(-14,-2,5),
#'                   sigma=c(2,4,1.5))
#' dist.plot<-generate_dist(example$Dist, GModel, 1e4)
#' thr <- find_thr_by_params(GModel,dist.plot)
#'
#' @seealso \code{\link{find_thr_by_dist}} and \code{\link{generate_dist}}
#' @export
find_thr_by_params <- function(GModel,input){
  tol = 1e-10
  thr2 <- c()
  for (i in 1 : (nrow(GModel)-1)){
    A = (1/(2*(GModel$sigma[i]^2)))-(1/(2*(GModel$sigma[i+1]^2)))
    B = GModel$mu[i+1]/(GModel$sigma[i+1]^2) - GModel$mu[i]/(GModel$sigma[i]^2)
    C = (GModel$mu[i]^2)/(2*GModel$sigma[i]^2) - (GModel$mu[i+1]^2)/(2*(GModel$sigma[i+1]^2)) - log((GModel$alpha[i]*GModel$sigma[i+1])/(GModel$alpha[i+1]*GModel$sigma[i]))

    if (abs(A) < tol){
      if(abs(B)<tol){
        print("Gaussians are the same!")
        x1 = NaN
        x2 = NaN
      }else{
        x1 = -C/B
        x2 = x1

        thr2 <- c(thr2, x1)
      }

    } else{
      delta = (B^2) - (4*A*C)

      if (delta<0){
        thr2<-c(thr2,find_thr_by_dist(input)[i])
      } else{

        x1 = (-B - sqrt(delta))/(2*A)
        x2 = (-B + sqrt(delta))/(2*A)

        if (x1 > GModel$mu[i] && x1 < GModel$mu[i+1]){
          thr2 <- c(thr2, x1)
        } else if (x2 > GModel$mu[i] && x2 < GModel$mu[i+1]){
          thr2 <- c(thr2, x2)
        } else {
          d1<-min(c(abs(x1-GModel$mu[i]),abs(x1-GModel$mu[i+1])))
          d2<-min(c(abs(x2-GModel$mu[i]),abs(x2-GModel$mu[i+1])))
          if (d1<d2){
            thr2<-c(thr2, x1)
          } else{
            thr2<-c(thr2, x2)
          }
        }
      }
    }
  }

  return(thr2)
}

#' Calculation of GMM thresholds from distribution
#'
#' Function to calculate cutoffs between each component of mixture normal distributions using probability distribution function.
#'
#' @param input output of \code{generate_dist} function. It is a list with following elements:\describe{
#'    \item{x}{Numeric vector with equaliy spread data of given precison}
#'    \item{dist}{Matrix with PDF of each GMM component and cumulative distribution}
#' }
#'
#' @returns Return a vector of thresholds.
#'
#' @examples
#' data(example)
#' GModel<-data.frame(alpha=c(0.45,0.5,0.05),
#'                   mu=c(-14,-2,5),
#'                   sigma=c(2,4,1.5))
#' dist.plot<-generate_dist(example$Dist, GModel, 1e4)
#' thr <- find_thr_by_dist(dist.plot)
#'
#' @seealso \code{\link{find_thr_by_params}} and \code{\link{generate_dist}}
#' @export
find_thr_by_dist<- function(input){
  x<-input$x
  KS<-(ncol(input$dist)-1)
  dist<-input$dist[,1:KS]

  if (KS==1){
    threshold<-NA
  } else{
    thr<-c()
    ind<-1
    for (i in 1:(KS-1)){
      if (i==1){ f1 = dist[,ind]} else { f1=rowSums(dist[,ind])}
      if (i==(KS-1)){ f2 = dist[,(-ind)]} else { f2=rowSums(dist[,-ind])}

      f_diff = abs(f1-f2)
      ix<-order(f_diff)
      a = 1
      thr_ind = ix[a]
      while (thr_ind == 1 || thr_ind == length(x)){
        a<-a+1
        thr_ind = ix[a]
      }

      thr[i] = x[thr_ind]
      ind<-c(ind,i+1)
    }
  }

  return(thr)
}
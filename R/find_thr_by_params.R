find_thr_by_params <- function(GModel,input){
  # function to find thr based on GMM params
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
          #print(paste0("For ",i, " im in dist select"))
      } else{
      
          x1 = (-B - sqrt(delta))/(2*A)
          x2 = (-B + sqrt(delta))/(2*A)
          
          if (x1 > GModel$mu[i] && x1 < GModel$mu[i+1]){
            thr2 <- c(thr2, x1)
          #  print(paste0("For ",i, " select x1"))
          } else if (x2 > GModel$mu[i] && x2 < GModel$mu[i+1]){
            thr2 <- c(thr2, x2)
          #  print(paste0("For ",i, " select x2"))
          } else {
          # print(paste0("For ",i, " in new else"))
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

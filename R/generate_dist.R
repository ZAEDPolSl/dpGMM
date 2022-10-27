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
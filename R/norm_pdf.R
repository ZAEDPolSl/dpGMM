norm_pdf <- function(x, mu, sigma){
  y = exp(-0.5*(((x-mu)/sigma)^2))/(2.506628274631*sigma)
  return(y)
}
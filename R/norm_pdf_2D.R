norm_pdf_2D <- function(x, center, covar){
  den <- (6.283185307179585 * sqrt(det(covar)))
  x_centr <- sweep(as.matrix(x), 2, center, FUN = "-") 
  tmp <- sweep(x_centr, 2, colSums(covar), FUN = "/") * x_centr
  num <- - 0.5 * rowSums(tmp)
  y <- exp(num)/den
  return(y)
}
my_qu_ix_w <- function(invec, yinwec){
  PAR = 1
  PAR_sig_min =.1
  if ((invec[length(invec)]-invec[1]) <= PAR_sig_min || sum(yinwec) <= 1.0e-3){
    wyn = Inf
  }else {
    wwec = yinwec/sum(yinwec)
    wyn1 = (PAR + sqrt(sum(((invec - sum(invec*wwec))^2)*wwec)))/(max(invec) - min(invec))
    wyn = wyn1
  }
  
  return(wyn)
}
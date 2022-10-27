find_thr_by_dist<- function(input){
  # function to find thr based on distribution
  # input <- output from generate_distribution 
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


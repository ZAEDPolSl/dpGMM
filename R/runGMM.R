runGMM <- function(X, KS, Y = NULL, change = Inf, max_iter = 5000, SW=NULL, IC = "BIC", merge = TRUE, sigmas.dev = 2.5,
                   precision=1e4, plot = TRUE, col.pal="Blues", quick_stop = TRUE, signi = 0.05) {
  # Check part
  if (!hasArg("X")){
    stop("No data.")}

  if (length(X) < 2){
    stop("Not enough data.")}

  IC_list <- c("AIC","AICc", "BIC", "ICL-BIC", "LR")
  if (!IC %in% IC_list) {
    stop("Criterion not implemented. Please use AIC, AICc, BIC, ICL-BIC, LR, gap or NbClust")
  }


  # Main GMM run
    GModel<- gaussian_mixture_vector(X, KS, Y, change, max_iter, SW, IC, quick_stop, signi)

    if(merge & GModel$KS>1){
      IC_tmp <- GModel$IC
      logL_tmp <- GModel$logL
      GModel <- gmm_merge(GModel$model, sigmas.dev)
      GModel$IC <- IC_tmp
      GModel$logL <- logL_tmp
    }

  # Generating distribution from model
    dist.plot<-generate_dist(X, GModel$model, precision)

  # Thresholds estimation
    #thr <- find_thr_by_dist_VMM(dist.plot)
    if(GModel$KS>1){
      thr <- find_thr_by_params(GModel$model,dist.plot)
    } else {thr=NULL}

  # Clusters assignment
    clust <- matrix(1, 1, length(X))
    for(i in 1:length(thr)){clust[X>thr[i]] <- i+1}

  # Plot generating
    pl<-plot_gmm_1D(X, dist.plot, Y, thr, pal=col.pal)

  # QQplot
    pl.qq<-plot_QQplot(X,GModel)

  # Output of function
    mix_gmm <- list(model = GModel$model, KS = nrow(GModel$model), IC = GModel$IC, logLik = GModel$logL,
                    threshold = thr, cluster = as.vector(clust), fig=pl,QQplot=pl.qq)
    names(mix_gmm)[3] <- IC

  # Print the plot
    if(plot){
      p<-mix_gmm[["fig"]]
      print(p)}

  return(mix_gmm)
}

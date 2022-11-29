test_2D_init <- function(){

  require(ggplot2)

  # generate 2D data
  m <- 1000
  data2D <- generate_dset2D(m=m)

  res_all <- list()
  for (a in 1:m){
    cat("Simulation no",a,"\n")

    X <- data2D[[a]]$dist$Dist
    Y <- rep(1,nrow(X))
    KS <- data2D[[a]]$KS

    opts <- list()
    opts$SW <- ((max(X)-min(X))/(100*KS))^2
    opts$eps_change <- 1e-5
    opts$max_iter <- 5000
    opts$cov_type <- "sphere"
    opts$max_var_ratio <- 10
    opts$crit <- "BIC"
    opts$res <- 1

    # run GMM with different IC
    gmm <- list()
    init <- rand_init_2D(X, KS)
    #gmm[["rand"]] <- init
    gmm[["rand"]] <- EM_iter_2D(X, Y, init , opts);

    init <- DP_init_2D(X, Y, KS)
    #gmm[["DP"]] <- init
    gmm[["DP"]] <- EM_iter_2D(X, Y, init , opts);

    init <- diag_init_2D(X, KS)
    #gmm[["diag"]] <- init
    gmm[["diag"]] <- EM_iter_2D(X, Y, init , opts);

    # calculate difference between truth and estimated values
    res_tmp <- data.frame(matrix(nrow=length(gmm),ncol=5))
    res_tmp[,1] <- paste0("sim",a)
    res_tmp[,2] <- names(gmm)
    for (b in 1:length(gmm)){
      gmm_tmp <- gmm[[b]]
      gt <- data2D[[a]]

      res_tmp2 <- matrix(nrow=gt$KS,ncol=3)
      for (cc in 1:gt$KS){
        tmp <- dist(rbind(t(gt$mu[,cc]),gmm_tmp$center), method = "euclidean")
        tmp_ind <- which.min(tmp[1:nrow(gmm_tmp$center)])

        # Defining euclidean distance between parameters
        res_tmp2[cc,1] <- tmp[tmp_ind]
        res_tmp2[cc,2] <- dist(c(gt$sigma[cc],(gmm_tmp$covar[1,1,tmp_ind])), method = "euclidean")[1]
        res_tmp2[cc,3] <- dist(c(gt$alpha[cc],gmm_tmp$alpha[tmp_ind]), method = "euclidean")[1]

        gmm_tmp$center[tmp_ind,] <- NA
      }
      res_tmp[b,3] <- mean(res_tmp2[,1])
      res_tmp[b,4] <- mean(res_tmp2[,2])
      res_tmp[b,5] <- mean(res_tmp2[,3])

    }
    colnames(res_tmp) <- c("Simulation","Method","Diff_mean","Diff_sig","Diff_alpha")
    res_all[[a]] <- res_tmp
  }
  res_all <- do.call(rbind, res_all)

  p1 <- ggplot(res_all, aes(x=Method, y=Diff_mean)) + geom_violin() +
    geom_boxplot(width=0.2) + theme_bw()
  p2 <- ggplot(res_all, aes(x=Method, y=Diff_sig)) + geom_violin() +
    geom_boxplot(width=0.2) + theme_bw()
  p3 <- ggplot(res_all, aes(x=Method, y=Diff_alpha)) + geom_violin() +
    geom_boxplot(width=0.2) + theme_bw()

  ml <- gridExtra::arrangeGrob(grobs=list(p1,p2,p3), nrow=1, ncol=3)
  ggsave("Test_2D.pdf",ml,height= 4, width = 6)
}

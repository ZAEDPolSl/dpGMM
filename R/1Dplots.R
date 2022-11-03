#' Plot of GMM decomposition
#'
#' Function plot the decomposed distribution together with histogram of data. Moreover the cut-off are marked.
#' This plot is also return as regular output of \code{\link{runGMM}}.
#'
#' @param X Vector of data
#' @param dist Output of \code{generate_dist} function. Its list of following elements\describe{
#'    \item{x}{Vector of generated data}
#'    \item{dist}{Matrix of pdf of for each generated mixture model. The last column is sum of all previous ones}
#' }
#' @param Y Vector of X counts (dedicated to binned data). Default=NULL
#' @param threshold Vector with GMM cutoffs
#' @param pal RColorBrewer palette name
#'
#' @import ggplot2
#' @import RColorBrewer
#' @importFrom grDevices colorRampPalette
#' @importFrom reshape2 melt
#'
#' @seealso \code{\link{runGMM}}
#'
#' @examples
#' \dontrun{
#' data(example)
#' GModel<-data.frame(alpha=c(0.45,0.5,0.05),
#'                   mu=c(-14,-2,5),
#'                   sigma=c(2,4,1.5))
#' dist.plot<-generate_dist(example$Dist, GModel, 1e4)
#' thr <- find_thr_by_params(GModel,dist.plot)
#' plot_gmm_1D(example$Dist, dist.plot,Y=NULL,threshold=thr, pal="Dark2")
#' }
#'
#' @export
plot_gmm_1D <- function(X, dist, Y=NULL, threshold=NA, pal=NULL){
  # estimating of bin width of histogram
  binwidth = (max(X)-min(X))/floor(sqrt(length(X)))
  #binwidth = 2.64*IQR(data$V1)*nrow(data)^(-1/3) # alternative

  # extract data from dist
  x_temp<-dist$x
  dist<-dist$dist

  # organize data for line plot
  colnames(dist)=paste('Comp:',1:ncol(dist),sep='')
  colnames(dist)[ncol(dist)]<-"Main"
  #dist<-dist#*binwidth*length(data) in new version not needed

  tmp <- reshape2::melt(dist,id.vars = NULL)
  tmp$xx <- rep(x_temp, ncol(dist))
  tmp$lin <- 1
  tmp$lin[which(tmp$variable == "Main")] <- 0
  if (ncol(dist)==2){
    col <- c("darkgreen", "grey25")
  } else{
    col <- grDevices::colorRampPalette(brewer.pal(8,pal))(ncol(dist))
    col <- c(col[2:ncol(dist)], "grey25")}


  p <- ggplot() + theme_bw()

  if (!is.null(Y)){
    binwidth<-(max(X)-min(X))/length(X)
    Y<-Y/(sum(Y)*binwidth)
    p<- p+geom_bar(aes(x=X,y=Y),stat="identity",color="grey65", fill="grey65",alpha=0.2)

  } else{
    p<- p+geom_histogram(aes(x=X,y=..density..),binwidth = binwidth,color="black", fill="grey",alpha=0.2)
  }

  p<-p+geom_line(aes(x=tmp$xx, y=tmp$value, group=tmp$variable, color=as.factor(tmp$variable), linetype=as.factor(tmp$lin)), size=1)+
    scale_color_manual(values=col,name="") +
    scale_linetype_manual(values=c("dashed","solid")) + xlab("x")+ylab("Density")+guides(linetype="none")

  if (sum(!is.na(threshold)))
    p<-p+geom_vline(xintercept = threshold, lty = "dashed", col = "red")

  return(p)
}

#' QQplot for GMM decomposition
#'
#' Function return ggplot object with fit diagnostic Quantile-Quantile plot for one normal distribution and fitted GMM.
#' This plot is also return as regular output of \code{\link{runGMM}}.
#'
#' @param data Vector of original data
#' @param GModel \code{data.frame} of GMM parameters i.e GModel$alpha, GModel$mu, GModel$sigma (correct \code{colnames} are obligatory)
#'
#' @import ggplot2
#' @importFrom  ggpubr ggarrange
#' @importFrom stats qqplot
#'
#' @examples
#' \dontrun{
#' data(example)
#' GModel<-data.frame(alpha=c(0.45,0.5,0.05),
#'                   mu=c(-14,-2,5),
#'                   sigma=c(2,4,1.5))
#' plot_QQplot(example$Dist,GModel)
#' }
#'
#' @seealso \code{\link{runGMM}} and \code{\link{gaussian_mixture_vector}}
#'
#' @export
plot_QQplot<-function(data,GModel){
  tor<-generate_norm1D(length(data), GModel$alpha, GModel$mu, GModel$sigma)
  quants<-stats::qqplot(data,tor$Dist,plot.it = F)
  tmp<-data.frame(data=quants$x,theor=quants$y)

  p2<-ggplot(tmp,aes(theor,data))+theme_bw()+
    geom_smooth(method='lm',formula=y~x,color="red",se=F)+
    geom_point(size=2.5,shape=18,color="#324376")+
    ylab("Data")+xlab("GMM fit")+
    ggtitle("QQ plot: GMM")+theme(plot.title = element_text(hjust = 0.5))

  p1 <- ggplot(tmp, aes(sample = data))+ stat_qq(size=2.5,shape=18,color="#324376") + stat_qq_line(color="red")+
    ylab("Data")+xlab("Normal distribution")+theme_bw()+
    ggtitle("QQ plot: one dist.")+theme(plot.title = element_text(hjust = 0.5))

  return(ggpubr::ggarrange(p1,p2,align="hv"))
}

#' Plot of GMM decomposition for data vector
#'
#' Function plot the decomposed distribution together with histogram of data. Moreover the cut-off are marked.
#'
#' @param X Vector of data
#' @param dist Output of generate_dist function. Its list of following elements\describe{
#'    \item{x}{Vector of generated data}
#'    \item{dist}{Matrix of pdf of x for generated mixture model. The last column is sum of all previous ones.}
#' }
#' @param Y Vector of X counts (binned data). Default=NULL
#' @param threshold Vector with GMM cutoffs
#' @param pal RColorBrewer palette name
#'
#' @importFrom ggplot2 ggplot aes aes_string geom_smooth geom_point ylab xlab ggtitle theme scale_color_manual geom_histogram geom_bar
#' @importFrom ggplot2 guides scale_linetype_manual geom_vline geom_line
#' @importFrom stats qqplot
#' @importFrom reshape2 melt
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette
#'
#' @seealso \code{\link{runGMM}} and \code{\link{generate_dist}}
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

  tmp <- melt(dist,id.vars = NULL)
  tmp$xx <- rep(x_temp, ncol(dist))
  tmp$lin <- 1
  tmp$lin[which(tmp$variable == "Main")] <- 0
  if (ncol(dist)==2){
    col <- c("darkgreen", "grey25")
  } else{
  col <- colorRampPalette(brewer.pal(8,pal))(ncol(dist))
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

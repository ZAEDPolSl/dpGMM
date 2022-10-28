#' QQplot for GMM decomposition
#'
#' Function return ggplot object with fit diagnostic Quantile-Quantile plot for one normal distribution and fitted GMM.
#' The QQplot is also return as regular output of runGMM.
#'
#' @param data Vector of data
#' @param GModel List of GModel parameters i.e GModel$model$alpha, GModel$model$mu, GModel$model$sigma. Given as output of function gaussian_mixture_vector.
#'
#' @importFrom ggpubr ggarrange
#' @importFrom ggplot2 ggplot aes aes_string geom_smooth geom_point ylab xlab ggtitle theme
#' @importFrom ggplot2 stat_qq stat_qq_lin
#' @importFrom stats qqplot
#'
#' @examples
#' data<-generate_norm1D(1000, alpha=c(0.2,0.4,0.4), mu=c(-15,0,15), sigma=c(1,2,3))
#' GModel<-list(model=list(alpha=c(0.2,0.4,0.4), mu=c(-15,0,15), sigma=c(1,2,3)))
#' plot_QQplot(data,GModel)
#' @seealso \code{\link{runGMM}} and \code{\link{gaussian_mixture_vector}}
#'
#' @export
plot_QQplot<-function(data,GModel){
tor<-generate_norm1D(length(data), GModel$model$alpha, GModel$model$mu, GModel$model$sigma)
quants<-qqplot(data,tor$Dist,plot.it = F)
tmp<-data.frame(data=quants$x,theor=quants$y)

p2<-ggplot(tmp,aes(theor,data))+theme_bw()+
  geom_smooth(method='lm',formula=y~x,color="red",se=F)+
  geom_point(size=2.5,shape=18,color="#324376")+
  ylab("Data")+xlab("GMM fit")+
  ggtitle("QQ plot: GMM")+theme(plot.title = element_text(hjust = 0.5))

p1 <- ggplot(tmp, aes(sample = data))+ stat_qq(size=2.5,shape=18,color="#324376") + stat_qq_line(color="red")+
  ylab("Data")+xlab("Normal distribution")+theme_bw()+
  ggtitle("QQ plot: one dist.")+theme(plot.title = element_text(hjust = 0.5))

return(ggarrange(p1,p2,align="hv"))
}

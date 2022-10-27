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
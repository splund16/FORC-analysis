xx<-c(-200,-170,-125,-110,-100,-(5*18:14),-60,-30,20,100,150,200)

p<-ggplot(data=fit1d[fit1d$rnd.H_r%in%c(xx),],
          aes(x=H,y=pred,group=rnd.H_r,colour=rnd.H_r))+
  geom_line(size=2)+
  scale_colour_gradientn(colours=rainbow(6),
                         guide="colourbar",
                         limits=c(-200,200),
                         name=expression(mu[0]*H["R"]*" (mT)"))+
  #  geom_point()+
  xlab(expression(mu[0]*"H (mT)"))+
  ylab(expression("Moment  (x10"^"-3"*"A-m"^2*"/kg)"))+
  #guides(color=guide_legend(title=expression("H"["R"]*"  (Oe)")))+
  #geom_smooth()+
  theme_bw()+
  theme(text=element_text(size=22))
print(p)


xx<-c(-200,-75,-60,-30,20)
t.d<-fit1d[fit1d$rnd.H_r%in%c(xx),]
t.d$H_r_lab<-factor(paste0("Hr",t.d$rnd.H_r),levels=paste0("Hr",xx))
levels(t.d$H_r_lab)<-paste("mu[0]*'H'[R]*' = '*",xx,"*' mT'")
p<-ggplot(data=t.d,
          aes(x=H,y=dmom.dH,group=H_r,
              colour=rnd.H_r))+
  facet_wrap(vars(H_r_lab),labeller = label_parsed,ncol=5)+
  scale_colour_gradientn(colours=rainbow(6),
                         guide="colourbar",
                         limits=c(-200,200),
                         name=expression(mu[0]*"H"["R"]*"  (mT)"))+
  #  geom_line(data=fit1d$surface[fit1d$surface$rnd.H_r%in%c(xx),],aes(x=H,y=pred,group=H_r),colour="red",size=2,alpha=.5)+
  geom_point()+
  ylab(expression(partialdiff*" Moment / "*partialdiff*" "*mu[0]*"H  (A-m"^2*"/kgT)"))+
  xlab(expression(mu[0]*"H (mT)"))+
  # geom_smooth(color="red")+
  theme_bw()+
  theme(text=element_text(size=20),
        axis.text.x=element_text(angle=90,hjust=1,vjust=.4),
        axis.title=element_text(face="bold"))
print(p)

xx<-c(20,60,80,100,120)
t.d<-fit1d[fit1d$rnd.H%in%c(xx),]
t.d$h_lab<-factor(paste0("H",t.d$rnd.H),levels=paste0("H",xx))
levels(t.d$h_lab)<-paste("mu[0]*'H = '*",xx,"*' mT'")

p<-ggplot(data=t.d,
          aes(x=H_r,y=dmom.dH,colour=H_r))+
  scale_colour_gradientn(colours=rainbow(6),
                         guide="colourbar",
                         limits=c(-200,200),
                         name=expression(mu[0]*"H"["R"]*" (mT)"))+
  facet_wrap(vars(h_lab),labeller = label_parsed,ncol=5)+
  geom_point()+
  ylab(expression(partialdiff*" Moment / "*partialdiff*" "*mu[0]*"H  (A-m"^2*"/kgT)"))+
  xlab(expression(mu[0]*"H"["R"]*"  (mT)"))+
  theme_bw()+
  theme(text=element_text(size=22),
        axis.text.x=element_text(angle=90,hjust=1,vjust=.4),
          axis.title=element_text(face="bold"))
  # theme(axis.text=element_text(size=12),
  #       strip.text.x = element_text(size = 12),
  #       axis.title=element_text(size=14,face="bold"))
print(p)


p<-ggplot(data=t.d,
          aes(x=H_r,y=dmom.dH.pred))+
  facet_wrap(vars(h_lab),labeller = label_parsed,ncol=5)+
  geom_line()+
  geom_line(aes(y=weighted.dmom.dH.pred),color="red")+
  ylab(expression(partialdiff*" Moment / "*partialdiff*" "*mu[0]*"H  (A-m"^2*"/kgT)"))+
  xlab(expression(mu[0]*"H"["R"]*"  (mT)"))+
  theme_bw()+
  theme(text=element_text(size=22),
        axis.text.x=element_text(angle=90,hjust=1,vjust=.4),
        axis.title=element_text(face="bold"))
# theme(axis.text=element_text(size=12),
#       strip.text.x = element_text(size = 12),
#       axis.title=element_text(size=14,face="bold"))
print(p)

t.dat<-dat[dat$rnd.H==80,]

plot(t.dat$H_r,t.dat$dmom.dH.pred,type="l")
lines(t.dat$H_r,t.dat$weighted.dmom.dH.pred,col=2)


big_res<-fit1d%>%
  group_by(H_r)%>%
  summarise(sd=sd((Moment-pred)/Moment_SE))%>%
  arrange(-sd)

xx<-sort(2.5*round(big_res$H_r[1:12]/2.5))
t.d<-fit1d[fit1d$rnd.H_r%in%c(xx),]
t.d$H_r_lab<-factor(paste0("Hr",t.d$rnd.H_r),levels=paste0("Hr",xx))
levels(t.d$H_r_lab)<-paste("mu[0]*'H'[R]*' = '*",xx,"*' mT'")

p<-ggplot(data=t.d,
          aes(x=H,y=(Moment-pred)/Moment_SE,group=H_r,colour=H_r))+
  scale_colour_gradientn(colours=rainbow(6),
                         guide="colourbar",
                         limits=c(-200,200),
                         name=expression(mu[0]*"H"["R"]*" (mT)"))+
  facet_wrap(vars(H_r_lab),labeller = label_parsed)+
  geom_point()+
  xlab(expression(mu[0]*"H (mT)"))+
  ylab("Standardized Residual")+
  theme_bw()+
  theme(text=element_text(size=22),
        axis.text.x=element_text(angle=90,hjust=1,vjust=.4),
          axis.title=element_text(face="bold"))
  # theme(axis.text=element_text(size=12),
  #       strip.text.x = element_text(size = 12),
  #       axis.title=element_text(size=14,face="bold"))
print(p)

p<-ggplot(data=t.d,
          aes(x=H,y=Moment,group=H_r,colour=H_r))+
  scale_colour_gradientn(colours=rainbow(6),
                         guide="colourbar",
                         limits=c(-200,200),
                         name=expression(mu[0]*"H"["R"]*" (mT)"))+
  facet_wrap(vars(H_r_lab),labeller = label_parsed)+
  geom_point()+
  geom_line(data=t.d,aes(x=H,y=pred),colour="black")+
  xlab(expression(mu[0]*"H (mT)"))+
  ylab(expression("Moment  (x10"^"-3"*"A-m"^2*"/kg)"))+
  theme_bw()+
  theme(text=element_text(size=22),
        axis.text.x=element_text(angle=90,hjust=1,vjust=.4),
          axis.title=element_text(face="bold"))
  # theme(axis.text=element_text(size=12),
  #       strip.text.x = element_text(size = 12),
  #       axis.title=element_text(size=14,face="bold"))
print(p)



big_res<-fit1d%>%
  group_by(rnd.H)%>%
  summarise(sd=sd((Moment-pred)/Moment_SE))%>%
  arrange(-sd)

p<-ggplot(data=big_res,aes(x=rnd.H,y=sd))+
  geom_point()+theme_bw()+ xlab(expression(mu[0]*"H (mT)"))+
  ylab("Standard Deviation of 
Standardized Residual")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))
print(p)

xx<-sort(2.5*round(big_res$rnd.H[1:12]/2.5))
t.d<-fit1d[fit1d$rnd.H%in%c(xx),]
t.d$h_lab<-factor(paste0("H",t.d$rnd.H),levels=paste0("H",xx))
levels(t.d$h_lab)<-paste("mu[0]*'H = '*",xx,"*' mT'")

p<-ggplot(data=t.d,
          aes(x=H_r,y=(Moment-pred)/Moment_SE,group=h_lab,colour=H_r))+
  scale_colour_gradientn(colours=rainbow(6),
                         guide="colourbar",
                         limits=c(-200,200),
                         name=expression(mu[0]*"H (mT)"))+
  facet_wrap(vars(h_lab),labeller = label_parsed)+
  geom_point()+
  xlab(expression(mu[0]*"H"["R"]*" (mT)"))+
  ylab("Standardized Residual")+
  theme_bw()+
  theme(text=element_text(size=22),
        axis.text.x=element_text(angle=90,hjust=1,vjust=.4),
          axis.title=element_text(face="bold"))
  # theme(axis.text=element_text(size=12),
  #       strip.text.x = element_text(size = 12),
  #       axis.title=element_text(size=14,face="bold"))
print(p)




library(tseries)
par(mfrow=c(1,1))

AC<-data.frame(hr=NULL,lag=NULL,cor=NULL,cor_sign=NULL)
for(hr in unique(fit1d$rnd.H_r)){
  AC<-rbind(AC,data.frame(hr=hr,lag=0:20,
                          cor=acf((fit1d$Moment-fit1d$pred)[fit1d$rnd.H_r==hr ]
                                  ,lag=20,pl=FALSE)$acf[,,1],
                          cor_sign=acf(  sign(fit1d$Moment-fit1d$pred)[fit1d$rnd.H_r==hr ]
                                         ,lag=20,pl=FALSE)$acf[,,1]))
}

lag.labs <- paste("Lag =",(0:20)*2.5," mT")
names(lag.labs) <- as.character((0:20))

p<-ggplot(data=AC[AC$lag>0&AC$lag<10,],
          aes(x=as.numeric(as.character(hr)),y=cor))+
  facet_wrap(vars(lag),labeller=labeller(lag=lag.labs))+
  xlab(expression(mu[0]*"H"["R"]*" (mT)"))+
  ylab("Correlation")+
  geom_point()+
  theme_bw()+
  theme(text=element_text(size=22),
        axis.text.x=element_text(angle=90,hjust=1,vjust=.4),
          axis.title=element_text(face="bold"))
  # theme(axis.text=element_text(size=12),
  #       strip.text.x = element_text(size = 12),
  #       axis.title=element_text(size=14,face="bold"))
print(p)

p<-ggplot(data=AC[AC$lag>0&AC$lag<10,],
          aes(x=as.numeric(as.character(hr)),y=cor_sign))+
  facet_wrap(vars(lag),labeller=labeller(lag=lag.labs))+
  xlab(expression(mu[0]*"H"["R"]*" (mT)"))+
  ylab("Correlation of Residual Sign")+
  geom_point()+
  theme_bw()+
  theme(text=element_text(size=22),
        axis.text.x=element_text(angle=90,hjust=1,vjust=.4),
        axis.title=element_text(face="bold"))
  # theme(axis.text=element_text(size=12),
  #       strip.text.x = element_text(size = 12),
  #       axis.title=element_text(size=14,face="bold"))
print(p)


xx<-10*(-1:10)
t.d<-fit1d[fit1d$rnd.H%in%c(xx),]
t.d$h_lab<-factor(paste0("H",t.d$rnd.H),levels=paste0("H",xx))
levels(t.d$h_lab)<-paste("mu[0]*'H = '*",xx,"*' mT'")

p<-ggplot(data=t.d,
          aes(x=H_r,y=d2mom.dH.dHr))+
  facet_wrap(vars(h_lab),labeller = label_parsed)+
  xlab(expression(mu[0]*"H"["R"]*" (mT)"))+
  ylab(expression(partialdiff^2*" Moment / "*partialdiff*" "*mu[0]*"H"["R"]*" "*partialdiff*" "*mu[0]*"H  (x10"^3*"A-m"^2*"/kgT"^2*")"))+
  guides(color=guide_legend(title="H  (Oe)"))+
  #  geom_line(data=fit1d$surface[fit1d$surface$rnd.H_r%in%c(xx),],aes(x=H,y=pred,group=H_r),colour="red",size=2,alpha=.5)+
  geom_point()+
  # geom_smooth(color="red")+
  theme_bw()+
  theme(text=element_text(size=22),
        axis.text.x=element_text(angle=90,hjust=1,vjust=.4),
        axis.title=element_text(face="bold"))
  # theme(axis.text=element_text(size=12),
  #       strip.text.x = element_text(size = 12),
  #       axis.title=element_text(size=14,face="bold"))
print(p)

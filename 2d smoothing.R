source("Plotting Functions.R")

df <- data.frame(H = ave.dat[,"H"],
                 H_r = ave.dat[,"H_r"],
                 rnd.H = ave.dat[,"rnd.H"],
                 z = ave.dat[,"Moment"])
df$H_b<-(df$H+df$H_r)/2
df$H_c<-(df$H-df$H_r)/2

require("locfit")
mod<-locfit(z~lp(H,H_r,deg=3,nn=30/nrow(df)),data=df)
pred<-fitted(mod)
deriv.mod<-locfit(z~lp(H,H_r,deg=3,nn=30/nrow(df)),data=df,deriv=c(1,2))
deriv<-fitted(deriv.mod)

H_deriv.mod<-locfit(z~lp(H,H_r,deg=3,nn=30/nrow(df)),data=df,deriv=c(1))

# plots(ave.dat)
# plots(ave.dat,pred=pred,title="Locfit")
# plots(ave.dat,deriv=deriv,title="Locfit",solid=FALSE)

if(0){
  require("mgcv")
  # b2 <- gam(z ~  ti(H)+ ti(H_r)+ti(H,H_r, bs=c("tp", "tp")),data=df) ## BAD FIT!!
  # plots(dat,pred=b2$fitted.values,title="gam")
  # b3 <- gam(z~s(H_r,bs="tp",k=20)+s(H,bs="tp",k=20)+s(H,H_r,bs="tp",k=20),data=df) ## BAD FIT!!
  # plots(dat,pred=b3$fitted.values,title="gam")
  # b4 <- gam(z~s(H,H_r,bs="tp",fx=FALSE),data=df) ## BAD FIT!!
  # plots(dat,pred=b4$fitted.values,title="gam")
  # b5 <- gam(z~s(H_r,bs="cr",k=50)+s(H,bs="cr",k=50),data=df) ## BAD FIT!!
  # plots(dat,pred=b5$fitted.values,title="gam")
  b6 <- gam(z ~ te(H_r, H,k=c(40,40)), data = df)
  # b7 <- gam(z~s(H,H_r,bs="ad",k=15,m=3),data=df,gamma=1.4) ## BAD FIT!!
  # plots(dat,pred=b7$fitted.values,title="gam")
  # pred2<-b2$fitted.values
  # deriv2<-numeric.2nd.deriv(df=df,fit=b2)
  # plots(dat,pred=pred2,title="gam",all.plots = TRUE)
  # plots(dat,deriv=deriv2,title="gam")
  
  ## Thin plate splines in fields package
  library(fields)
  samp<-sample(1:nrow(dat),2500)
  dat2<-dat[samp,]
  df2<-df[samp,]
  Tps.fit<-Tps(x=dat2[,c("H","H_r")],Y=dat2[,"Moment"],scale.type="unscaled")
  pred2<-Tps.fit$fitted.values
  plots(dat2,pred=pred2,title="TPS")
  deriv2<-numeric.2nd.deriv(df=df2,fit=Tps.fit)
  plots(dat2,deriv=deriv2,title="TPS")
  
  
  library(MBA)
  spline <- mba.surf(data.frame(x,y,z),100,100)
  ?predict.mba.surf
  ?mba.surf
}

### Analyze residuals for temporal trends
full.pred<-predict(mod,newdata=data.frame(H=dat[,"H"],H_r=dat[,"H_r"]))
dat<-cbind(dat,resid=dat[,"Moment"]-full.pred)


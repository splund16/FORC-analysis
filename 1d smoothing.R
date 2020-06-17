require(fields)
require(dplyr)
fit.1d<-function(dat,  CV=FALSE){
  
  H<-"H"
  Moment<-"Moment"
  path<-dat$path
  ### Try fitting a single path using cubic splines
  dat$weighted.d2mom.dH.dHr<-dat$weighted.dmom.dH.pred<-dat$weighted.dmom.dH.cv.pred<-
    dat$d2mom.dH.dHr<-dat$dmom.dH.pred<-dat$dmom.dH.cv.pred<-dat$dmom.dH.resid<-
    dat$weighted.dmom.dH<-dat$weighted.cv.pred<-dat$weighted.pred<-
    dat$dmom.dH<-dat$cv.pred<-dat$pred<-rep(NA,nrow(dat))
  
  library(fields); 
  
  ### Estimate dmom.dH
  uni_hr<-unique(dat$H_r)
  weighted.spline.fit<-spline.fit<-list()
  spline.der<-function(spline.fit,x,delta.x){
    hi<-predict(spline.fit,x=x+delta.x)
    low<-predict(spline.fit,x=x-delta.x)
    .5*(hi-low)/delta.x
  }
  
  
  for(hr_ind in 1:length(uni_hr)){
    ind<-dat$H_r==uni_hr[hr_ind]
    invisible(spline.fit[[hr_ind]]<-sreg(dat$H[ind],dat$Moment[ind]))
    dat$dmom.dH[ind]<-spline.der(spline.fit[[hr_ind]],dat$rnd.H[ind],.0001)
    dat$pred[ind]<-spline.fit[[hr_ind]]$fitted.values
    k<-sample(rep(1:10,1000)[1:sum(ind)])
    if(CV){
      for(jj in 1:10){
        invisible(t.fit<-sreg(dat$H[ind][k!=jj],dat$Moment[ind][k!=jj]))
        dat$cv.pred[ind][k==jj]<-predict(t.fit,x=dat$H[ind][k==jj])
      }
    }
  }
  
  ### Evaluate d2mom.dHdHr
  ####Fit spline dmom.dH across H_r initial time
  spline.fit2<-spline.fit3<-list()
  h_table<-table(dat$rnd.H)
  for(ii in 1:sum(h_table>30)){
    H<-as.numeric(names(h_table)[ h_table>30])[ii]
    ind<-dat$rnd.H==H
    invisible(spline.fit2[[ii]]<-sreg(dat$H_r[ind],dat$dmom.dH[ind]))
    dat$dmom.dH.resid[ind]<-spline.fit2[[ii]]$residuals
    dat$dmom.dH.pred[ind]<-spline.fit2[[ii]]$fitted.values
    dat$d2mom.dH.dHr[ind]<-spline.der(spline.fit2[[ii]],dat$rnd.H_r[ind],.001)
    k<-sample(rep(1:10,1000)[1:sum(ind)])
    for(jj in 1:10){
      invisible(t.fit<-sreg(dat$H_r[ind][k!=jj],dat$dmom.dH[ind][k!=jj]))
      dat$dmom.dH.cv.pred[ind][k==jj]<-predict(t.fit,x=dat$H_r[ind][k==jj])
    }
  }
  
  if(nrow(dat)>5000){
    require(locfit)
    
    ## Construct weights to use when fitting Moment vs H
    dat$res<-dat$Moment-dat$pred
    dat$abs.res<-abs(dat$res)
    smooth.abs.res<-locfit(abs.res~lp(H,H_r,deg=3,nn=.001),data=dat)
    dat$weights<-1/fitted(smooth.abs.res)^2
    
    for(hr_ind in 1:length(uni_hr)){
      ind<-dat$H_r==uni_hr[hr_ind]
      invisible(weighted.spline.fit[[hr_ind]]<-sreg(dat$H[ind],dat$Moment[ind],weights = dat$weights[ind]))
      dat$weighted.dmom.dH[ind]<-spline.der(weighted.spline.fit[[hr_ind]],dat$rnd.H[ind],.0001)
      dat$weighted.pred[ind]<-weighted.spline.fit[[hr_ind]]$fitted.values
      k<-sample(rep(1:10,1000)[1:sum(ind)])
      if(CV){
        for(jj in 1:10){
          invisible(t.fit<-sreg(dat$H[ind][k!=jj],dat$Moment[ind][k!=jj],weights = dat$weights[ind][k!=jj]))
          dat$weighted.cv.pred[ind][k==jj]<-predict(t.fit,x=dat$H[ind][k==jj])
        }
      }
    }
    
    ### Refit splines to dmom.dH using weighted estimates
    for(ii in 1:sum(h_table>30)){
      H<-as.numeric(names(h_table)[ h_table>30])[ii]
      ind<-dat$rnd.H==H
      invisible(spline.fit3[[ii]]<-sreg(dat$H_r[ind],dat$weighted.dmom.dH[ind]))
      dat$weighted.dmom.dH.pred[ind]<-spline.fit3[[ii]]$fitted.values
    }
    
    ## Construct weights to use when fitting dMoment/dH vs H_r
    dat$res2<-dat$weighted.dmom.dH-dat$weighted.dmom.dH.pred
    dat$abs.res2<-abs(dat$res2)
    smooth.abs.res2<-locfit(abs.res2~lp(H,H_r,deg=3,nn=.001),data=dat)
    dat$weights2[!is.na(dat$weighted.dmom.dH.pred)]<-1/fitted(smooth.abs.res2)^2
    dat$weights2<- dat$weights2/mean( dat$weights2,na.rm=TRUE)
    for(ii in 1:sum(h_table>30)){
      H<-as.numeric(names(h_table)[ h_table>30])[ii]
      ind<-dat$rnd.H==H
      invisible(spline.fit3[[ii]]<-sreg(dat$H_r[ind],dat$weighted.dmom.dH[ind],weights=dat$weights2[ind]))
      dat$weighted.dmom.dH.pred[ind]<-spline.fit3[[ii]]$fitted.values
      dat$weighted.d2mom.dH.dHr[ind]<-spline.der(spline.fit3[[ii]],dat$rnd.H_r[ind],.001)
      k<-sample(rep(1:10,1000)[1:sum(ind)])
      if(CV){
        for(jj in 1:10){
          invisible(t.fit<-sreg(dat$H_r[ind][k!=jj],dat$weighted.dmom.dH[ind][k!=jj],weights=dat$weights2[ind][k!=jj]))
          dat$weighted.dmom.dH.cv.pred[ind][k==jj]<-predict(t.fit,x=dat$H_r[ind][k==jj])
        }
      }
    }
  }
  dat
}




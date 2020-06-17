setwd("C:/Users/spl/Documents/Magnetic Nanoparticle with Cindi Dennis/")
library(akima)
files<- c("CD032620.dat", "CD040220.dat", "CD040320.dat", "CD040920.dat", "CD041020.dat",
          "CD041620.dat", "CD041720.dat", "CD042320.dat", "CD043020.dat", "CD050120.dat",
          "CD050720.dat", "CD050820.dat", "CD051420.dat", "CD051520.dat", "CD052120.dat",
          "CD052220.dat", "CD101117.dat", "CD121017.dat")
for(f in 2:length(files)){
  try({
    File.name <-gsub(".dat","",files[f])
    File.name.mod<-""
    ## Does magfield decrease within minor loops?
    BACKWARD<-0
    MagField<-"Magnetic Field (Oe)"
    Moment<-"Moment (emu)"
    MomSE<-"M. Std. Err. (emu)"
    
    #File.name="FORC022814"; BACKWARD<-1
    
    ### Load data
    source("LoadFORCdata2.R")
    
    ### Use average across repeated loops and repeated measurements to conduct 2d smoothing
    source("2d smoothing.R")
    
    ### Compare repeated measurements for given h_r, h within a given minor loop path
    if(max(n)>1) source("RepeatedMeasurementsWithinMinorLoop.R")
    
    ### Compare Results Across Repeated Minor Loops
    if(any(table(rnd.H_r)>1)) source("RepeatedMinorLoops_new.R")
    
    source("1d smoothing.R")
    fit1d<-fit.1d(dat=ave.dat)
    
    fit1d$compare2d_pred<-predict(mod,newdata=fit1d)
    fit1d$compare2d_dmom.dH<-predict(H_deriv.mod,newdata=fit1d)
    fit1d$compare2d_d2mom.dH.dHr<-predict(deriv.mod,newdata=fit1d)
    
    
    #image.plot(x=fit1d$rnd.H,y=fit1d$rnd.H_r,z=fit1d$d2mom.dH.dHr)
    
    ### Create a matrix with derivative values
    uni_x<-sort(unique(fit1d$rnd.H))
    uni_y<-sort(unique(fit1d$rnd.H_r))
    pdf(paste0("Results/",File.name,"/",File.name,
               "_Selected_Fitting_results.pdf"),height=7,width=12)   
    x<-1:length(uni_x);names(x)<-uni_x
    y<-1:length(uni_y);names(y)<-uni_y
    z<-matrix(NA,length(y),length(x))
    z[cbind(y[as.character(fit1d$rnd.H_r)],x[as.character(fit1d$rnd.H)])]<-fit1d$d2mom.dH.dHr
    image.plot(z,y=uni_x,x=uni_y)
    
    uni_x<-sort(unique(fit1d$H_b))
    uni_y<-sort(unique(fit1d$H_c))
    
    x<-1:length(uni_x);names(x)<-uni_x
    y<-1:length(uni_y);names(y)<-uni_y
    z<-matrix(NA,length(y),length(x))
    z[cbind(y[as.character(fit1d$H_c)],x[as.character(fit1d$H_b)])]<-fit1d$d2mom.dH.dHr
    image.plot(z,y=uni_x,x=uni_y)
    
    library(akima)
    surf<--fit1d$compare2d_d2mom.dHdHr
    zz<-akima::interp(x=fit1d$H_r,y=fit1d$H,z=surf,nx=200,ny=200)
    image.plot(zz$z,x=zz$x,y=zz$y,xlab="H_r",ylab="H",main="Locfit 2nd derivative")
    
    surf<--fit1d$compare2d_d2mom.dHdHr
    zz<-akima::interp(x=fit1d$H_b,y=fit1d$H_c,z=surf,nx=200,ny=200)
    image.plot(zz$z,x=zz$x,y=zz$y,xlab="H_b",ylab="H_c",main="Locfit 2nd derivative")
    
    surf<-(fit1d$Moment-fit1d$compare2d_pred)
    zz<-akima::interp(x=fit1d$H_r,y=fit1d$H,z=surf,nx=200,ny=200)
    image.plot(zz$z,x=zz$x,y=zz$y,xlab="H_r",ylab="H",main="Locfit residual")
    
    surf<-(fit1d$Moment-fit1d$compare2d_pred)/fit1d$Moment_SE
    zz<-akima::interp(x=fit1d$H_r,y=fit1d$H,z=surf,nx=200,ny=200)
    image.plot(zz$z,x=zz$x,y=zz$y,xlab="H_r",ylab="H",main="Locfit residual/SE(Moment)",zlab="Residual/SE(Moment)")
    
    z[cbind(y[as.character(fit1d$H_r)],x[as.character(fit1d$H)])]<- -fit1d$compare2d_smooth
    image.plot(z,y=uni_x,x=uni_y)
    
    image.plot(matrix(c(rep(0,6),1,1,1),3,3))
    
    ind<-abs(fit1d$H_c+20)<50&abs(fit1d$H_b-80)<40
    
    # ### First derivative wrt H according to 1-d cubic splines
    # plots(fit1d,pred=fit1d$dH,which.plots=1:3)
    # 
    # ### First derivative wrt H according to 2-d locfit
    # plots(fit1d,pred=fit1d$compare2d_dmom_dH,which.plots=2)
    # 
    # ### Difference in first derivative
    # plots(fit1d,pred=fit1d$compare2d_dmom_dH-fit1d$dH,which.plots=2)
    
    ###  See relationship in estimated first derivatives
    plot(fit1d$dH,fit1d$compare2d_dmom_dH)
    
    # ### Difference in smoothed surface
    # plots(fit1d,pred=fit1d$compare2d_pred-fit1d$pred,which.plots=2)
    # 
    # ###  See relationship in residuals for fitted "surface"
    plot(fit1d$Moment-fit1d$pred,fit1d$Moment-fit1d$compare2d_pred)
    
    # ### Residuals from 1-d smoothing
    # plots(fit1d,pred=fit1d$Moment-fit1d$pred,which.plots=2)
    # 
    # ### Residuals from 2-d smoothing
    # plots(fit1d,pred=fit1d$Moment-fit1d$compare2d_pred,which.plots=2)
    # 
    # ### Residuals from 1-d smoothing in terms of number of standard errors
    # plots(fit1d,pred=(fit1d$Moment-fit1d$pred)/fit1d$Moment_SE,which.plots=2)
    # 
    # ### Residuals from 2-d smoothing in terms of number of standard errors
    # plots(fit1d,pred=(fit1d$Moment-fit1d$compare2d_pred)/fit1d$Moment_SE,which.plots=2) 
    
    d1<-density((fit1d$Moment-fit1d$compare2d_pred)/fit1d$Moment_SE,from=-10,to=10)
    d2<-density((fit1d$Moment-fit1d$pred)/fit1d$Moment_SE,from=-10,to=10)
    plot(d1,ylim=c(0,max(c(d1$y,d2$y))),xlab="Moment Residual / Moment SE")
    lines(d2,col=2)
    legend("topright",legend=c("2d locfit","1d cubic spline"),col=1:2,lty=1)
    
    # plots(fit1d,deriv=fit1d$d2mom.dH.dHr,solid=FALSE,which.plots=5)
    # plots(fit1d,deriv=-fit1d$compare2d_d2mom.dH.dHr,solid=FALSE,which.plots=5)
    # plots(fit1d,deriv=fit1d$d2mom.dH.dHr+fit1d$compare2d_d2mom.dH.dHr,solid=FALSE,which.plots=5)
    
    
    ind<-abs(fit1d$H_c+20)<50&abs(fit1d$H_b-80)<40
    
    
    plot(fit1d$d2mom.dH.dHr[ind],fit1d$compare2d_d2mom.dH.dHr[ind])
    plot(fit1d$weighted.d2mom.dH.dHr,fit1d$compare2d_d2mom.dH.dHr,)
    
    #plots(fit1d[ind,],deriv=fit1d$d2mom.dH.dHr[ind],solid=FALSE,which.plots=4:5)
    #plots(fit1d[ind,],deriv=-fit1d$compare2d_d2mom.dH.dHr[ind],solid=FALSE,which.plots=4:5)
    
    
    ### Some characteristics of the overall experimental design and system behavior
    source("Experiment Overview.R")
    
    layout(matrix(1:6,2,3))

    
    pres_plots<-function(dat1,dat2=NULL,method="",exclude_rng=NULL){
      library(akima)
      if(is.null(dat2))dat2<-dat1
      CA<-CL<-1.7
      par(mfcol=c(2,3))
      par(mai=c(.7,.7,1,1.2))
      nx<-ny<-300
      
      zz<-akima::interp(x=dat1$H_r,y=dat1$H,z=dat1$Moment,nx=nx,ny=ny)
      fields::image.plot(zz$z,x=zz$x,y=zz$y,xlab=expression(mu[0]*H[R]),
                         ylab=expression(mu[0]*H),main="Moment",
                         zlim=range(ave.dat$Moment),cex.lab=CL,cex.axis=CA)
      
      zz<-akima::interp(x=dat1$H_r,y=dat1$H,z=dat1$pred,nx=nx,ny=ny)
      fields::image.plot(zz$z,x=zz$x,y=zz$y,
                         xlab=expression(mu[0]*H[R]*" (mT)"),ylab=expression(mu[0]*H*" (mT)"),
                         main=paste(method,"
Smoothed Moment"),
                         xlim=range(ave.dat$H_r),ylim=range(ave.dat$H),zlim=range(ave.dat$Moment),
                         cex.lab=CL,cex.axis=CA)
      dat1$rnd.H_r<-2.5*round(dat1$H_r/2.5)
      
      H<-sort(unique(dat1$rnd.H))
      H_r<-sort(unique(dat1$rnd.H_r))
      zz<-rep(NA,length(H)*length(H_r))
      names(zz)<-paste(rep(H_r,each=length(H)),rep(H,length(H_r)))
      zz[paste(dat1$rnd.H_r,dat1$rnd.H)]<-(dat1$Moment-dat1$pred)/dat1$Moment_SE
      zmat<-matrix(zz,length(H),length(H_r))
      if(!is.null(exclude_rng))zmat[H>exclude_rng[1]&H<exclude_rng[2],]<-NA
      fields::image.plot(t(zmat),x=H_r,y=H,
                         xlab=expression(mu[0]*H[R]*" (mT)"),
                         ylab=expression(mu[0]*H*" (mT)"),
                         main=paste(method,"Standardized Residual"),
                         xlim=range(H_r),
                         ylim=range(H),
                         zlim=quantile(zmat,c(.003,.997),na.rm=TRUE),
                         cex.lab=CL,cex.axis=CA)
      
      par(mai=c(.7,.7,1,.2))
      norm.res<-(dat1$Moment-dat1$pred)/dat1$Moment_SE
      d<-density(norm.res)
      Ylim=range(d$y)
      Xlim=range(c(-3,3,d$x[d$y>(Ylim[2]/100)]))
      plot(density(norm.res,from=Xlim[1]-.05*diff(Xlim),to=Xlim[2]+.05*diff(Xlim)),
           xlim=Xlim,
           ylim=Ylim,
           lwd=2,
           lty=2,
           xlab="Standardized Residual",
           main="",
           cex.lab=CL,cex.axis=CA)
      grid()
      xx<-seq(-10,10,.01)
      #lines(xx,dt(xx,24),col=4,lwd=2)
      #    text(sum(c(.15,.85)*Xlim),(sum(c(.05,.95)*(Ylim))),expression(sigma),cex=CA)
      #text(sum(c(.15,.85)*Xlim),(sum(c(.15,.85)*(Ylim))),1.04,col=4,cex=CA)
      #    text(sum(c(.15,.85)*Xlim),(sum(c(.15,.85)*(Ylim))),signif(sd(norm.res),3),cex=CA)
      if("cv.pred"%in%names(dat1)){
        if(sum(!is.na(dat1$cv.pred))>100){
          ind<-which(!is.na(dat1$cv.pred))
          norm.cv.res<-((dat1$Moment-dat1$cv.pred)/dat1$Moment_SE)[ind]
          lines(density(norm.cv.res,from=Xlim[1]-.05*diff(Xlim),to=Xlim[2]+.05*diff(Xlim)),lwd=2)
          legend("topleft",title= "CV (SD)",
                 legend=c(paste0("Yes (",signif(sd(norm.cv.res),3),")"),
                          paste0("No (",signif(sd(norm.res),3),")")),
                 lty=c(1,2),lwd=2,cex=1.3,bg="white")
        }
      } 
      par(mai=c(.7,.7,1,1.2))
      
      library(stats)
      library(lmenssp)
      dat3<-dat2[!is.na(dat2$dmomdhdhr),]
      zz<-akima::interp(x=dat3$H_r,y=dat3$H,z=-dat3$dmomdhdhr,nx=nx,ny=ny)
      fields::image.plot(zz$z,x=zz$x,y=zz$y,xlab=expression(mu[0]*H[R]*" (mT)"),
                         ylab=expression(mu[0]*H*" (mT)"),
                         main=bquote(.(method)~":" ~ rho),
                         xlim=range(dat3$H_r),ylim=range(dat3$H),
                         cex.lab=CL,cex.axis=CA)
      
      zz<-akima::interp(x=dat3$H_b,y=dat3$H_c,z=-dat3$dmomdhdhr,nx=nx,ny=ny)
      fields::image.plot(zz$z,x=zz$x,y=zz$y,xlab=expression(mu[0]*H[b]*" (mT)"),
                         ylab=expression(mu[0]*H[c]*" (mT)"),
                         main=bquote(.(method)~":" ~ rho),
                         xlim=range(dat3$H_b),ylim=range(dat3$H_c),
                         cex.lab=CL,cex.axis=CA)
      
      par(mfrow=c(1,3))
      par(mai=c(.7,.7,.6,.2))
      fields::image.plot(zz$z,x=zz$x,y=zz$y,
                         xlab=expression(mu[0]*H[b]*" (mT)"),ylab=expression(mu[0]*H[c]*" (mT)"),main=bquote(.(method)~":" ~ rho),
                         xlim=range(dat3$H_b),ylim=range(dat3$H_c),
                         cex.lab=CL,cex.axis=CA)
      plot(zz$x,rowSums(zz$z,na.rm=TRUE)*mean(diff(zz$x)),
           xlab=expression(mu[0]*H[b]*" (mT)"),ylab="Density",main=paste(method),type="l",
           xlim=range(dat3$H_b),cex.lab=CL,cex.axis=CA)
      plot(zz$y,colSums(zz$z,na.rm=TRUE)*mean(diff(zz$y)),xlab=expression(mu[0]*H[c]*" (mT)"),
           ylab="Density",main=paste(method),type="l",
           xlim=range(dat3$H_c),cex.lab=CL,cex.axis=CA)
      
      mat1<-cbind(zz$x,rowSums(zz$z,na.rm=TRUE)*mean(diff(zz$y)))
      colnames(mat1)<-c("H_b","Density")
      
      mat2<-cbind(zz$y,colSums(zz$z,na.rm=TRUE)*mean(diff(zz$x)))
      colnames(mat2)<-c("H_c","Density")
      
      write.csv(mat1,
                file=paste0("Results/",File.name,gsub(":","through",method),"_H_b.csv"),
                row.names=FALSE)
      write.csv(mat2,
                file=paste0("Results/",File.name,gsub(":","through",method),"_H_c.csv"),
                row.names=FALSE)
    }
    
    fit1d$dmomdhdhr<-fit1d$d2mom.dH.dHr
    pres_plots(dat1=fit1d,method="Unweighted Sequential 1d Splines")
    
    fit1d$dmomdhdhr<-fit1d$compare2d_d2mom.dH.dHr
    fit1d$pred<-fit1d$compare2d_pred
    pres_plots(dat1=fit1d,method="LocFit Quadratic")
    
    dev.off()
  })
}

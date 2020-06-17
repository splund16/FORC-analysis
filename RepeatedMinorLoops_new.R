if(!"min.loop.rep"%in%ls())min.loop.rep<-2
repeated<-as.numeric(names(table(rnd.H_r)[table(rnd.H_r)>=min.loop.rep]))

matching.paths<-matrix(NA,length(repeated),max(table(rnd.H_r)))
ob.ind1<-ob.ind2<-NULL
for(i in 1:length(repeated)){
  matching.paths[i,1:sum(rnd.H_r==repeated[i])]<-which(rnd.H_r==repeated[i])
  ob.ind1<-c(ob.ind1,which(path==matching.paths[i,1]))
  ob.ind2<-c(ob.ind2,which(path==matching.paths[i,2]))
}
rownames(matching.paths)<- repeated

pdf(file=paste0(getwd(),"/Results/",File.name,"/",File.name," Repeated Minor Loop Plots.pdf"),height=8,width=10)

H_r.list<-sort(unique(rnd.H_r))
for(r in 1:length(repeated)){
  t.dat3<-dat[dat[,"path"]%in%matching.paths[r,],]
  H_use<-unique(t.dat3[,"rnd.H"])
  if(length(H_use)>10) H_use<-H_use[c(1,round((1:9)*length(H_use)/9))]
  path.col<-colorRampPalette(c("purple","blue","light green", "orange", "red"))(sum(!is.na(matching.paths[r,])))
  field.col<-colorRampPalette(c("purple","blue","light green", "orange", "red"))(length(H_use))
  #field.ind<-round((H_use-min(H_use))/diff(range(H_use))*50)+1
  
  leg<-matching.paths[r,!is.na(matching.paths[r,])]
  prefit<-smooth.spline(y=t.dat3[,Moment],x=t.dat3[,"rnd.H"])
  fit<-smooth.spline(y=t.dat3[,Moment],x=t.dat3[,"H"],df=prefit$df)
  names(fit$y)<-fit$x
  resid<-t.dat3[,Moment]-fit$y[as.character(t.dat3[,"H"])]
  global.resid<-t.dat3[,"resid"]
  
  
  plot(.5*(resid+global.resid),resid-global.resid,
       main=bquote("Repeated Minor Loop: H"["R"] ~ "=" ~ .(repeated[r]) ~ "(Oe)"),
       xlab="Average of Global and Local Residuals",
       ylab="Local Residual - Global Residual")
  
  plot(1e6,1e6,ylim=range(t.dat3[,Moment])+c(0,diff(range(t.dat3[,Moment]))/3.5),xlim=range(t.dat3[,"H"]),
       ylab="Moment (emu)",
       xlab="Applied Field (Oe)",
       main=bquote("Repeated Minor Loop: H"["R"] ~ "=" ~ .(repeated[r]) ~ "(Oe)"))
  grid()
  for(P in 1:length(path.col)){
    p<-matching.paths[r,P]
    t.dat4<-t.dat3[t.dat3[,"path"]==matching.paths[r,P],]
    for(f in unique(t.dat4[,"H"])){
      yy<-t.dat4[t.dat4[,"H"]==f,Moment]
      lines(f+1:length(yy),yy,col=path.col[P])   }
  }
  lines(fit$x,fit$y,lwd=1,col=makeTransparent(1))
  legend("top",title="Path",legend=leg,ncol=8,fill=path.col)
  YLIM<-range(c(resid,global.resid))
  YLIM<-YLIM+c(0,diff(YLIM)/3.5)

  for(type in c("Local","Global")){
    if(type=="Local"){
      t.resid<-resid
      mod.descrip<-bquote("Moment Residual from Smoothing Spline for H "["R"] ~ "=" ~ .(repeated[r]) ~ "(Oe), (emu)")
    }
    if(type=="Global"){
      t.resid<-global.resid
      mod.descrip<-"Moment Residual from Global Fit, (emu)"
    }
  
    XLIM<-range(t.dat3[,"H"])
    plot(1e6,1e6,ylim=YLIM,xlim=XLIM,
         ylab=mod.descrip,
         xlab="Applied Field (Oe)",cex.lab=.9,
         main=bquote("Repeated Minor Loop: H"["R"] ~ "=" ~ .(repeated[r]) ~ "(Oe)"),
         sub="Line segments display repeated measurements at fixed field")
    grid()
    for(P in 1:length(path.col)){
      p<-matching.paths[r,P]
      t.dat4<-t.dat3[t.dat3[,"path"]==p,]
      for(f in H_use){
        yy<-t.resid[t.dat3[,"path"]==p][t.dat4[,"rnd.H"]==f]
        if(length(yy>0))
          lines(f+(1:length(yy))*diff(XLIM)/(70*length(yy)),yy,
                col=path.col[P])
      }
    }
    
    leg2<-leg
    t.col2<-c(path.col,"gray","black")
    for(side in c(-1,1)){
      t.H_r<-repeated[r]
      H_r.side<-H_r.list[which(H_r.list==t.H_r)+side]
      leg2<-c(leg2,paste0("H_r=",H_r.side))
      t.dat.side<-as.data.frame(dat[dat[,"rnd.H_r"]==H_r.side,c(Moment,"H")])
      t.dat.side$resid<-t.dat.side[,Moment]-predict(fit,x=t.dat.side[,"H"])$y
      t.dat.side.magfield<-unique(t.dat.side$H)
      for(f in H_use){
        m<-unique(t.dat.side$H)[which.min((unique(t.dat.side$H)-f)^2)]
        yy<-t.dat.side$resid[t.dat.side$H==m]
        lines(m+(1:length(yy))*diff(XLIM)/(70*length(yy)),yy,
              col=c("gray","black")[(side+3)/2])
      }
    }
    legend("top",title="Path",legend=leg2,ncol=5,fill=t.col2,bty="n")
  }
   
  for(type in c("Local","Global")){
    if(type=="Local"){
      t.resid<-resid
      mod.descrip<-bquote("Moment Residual from Smoothing Spline for H "["R"] ~ "=" ~ .(repeated[r]) ~ "(Oe), (emu)")
    }
    if(type=="Global"){
      t.resid<-global.resid
      mod.descrip<-"Moment Residual from Global Fit, (emu)"
    } 
    XLIM=range(t.dat3[,"H_r"]-repeated[r])
    XLIM=XLIM+c(-1,1)*diff(XLIM)/25
    plot(1e6,1e6,ylim=YLIM,xlim=XLIM,
         ylab=mod.descrip,
         xlab="Measured Reversal Field - Nominal Reversal Field (Oe)",
         main=bquote("Repeated Minor Loop: H"["R"] ~ "=" ~ .(repeated[r]) ~ "(Oe)"),
         sub="Line segments display repeated measurements at fixed field")
    grid()
    for(P in 1:length(path.col)){
      p<-matching.paths[r,P]
      t.dat4<-t.dat3[t.dat3[,"path"]==p,]
      for(ff in 1:length(H_use)){
        f<-H_use[ff]
        t.ind<-t.dat3[,"path"]==p&(t.dat3[,"rnd.H"]==f)
        yy<-t.resid[t.ind]
        xx<-t.dat3[t.ind,"H_r"]-repeated[r]
        lines(xx+(1:length(yy))*diff(XLIM)/(50*length(yy)),yy,col=field.col[ff])
      }
      text(mean(t.dat3[t.dat3[,"path"]==p,"H_r"]-repeated[r]),mean(resid[t.dat3[,"path"]==p]),p,pos=4)
    }
    
    legend("top",title="Field",legend=H_use,ncol=5,fill=field.col)
  }
  for(type in c("Local","Global")){
    if(type=="Local"){
      t.resid<-resid
      mod.descrip<-bquote("Moment Residual from Smoothing Spline for H "["R"] ~ "=" ~ .(repeated[r]) ~ "(Oe), (emu)")
    }
    if(type=="Global"){
      t.resid<-global.resid
      mod.descrip<-"Moment Residual from Global Fit, (emu)"
    }
    XLIM=range(t.dat3[,"path"])
    plot(1e6,1e6,ylim=YLIM,xlim=XLIM,
         ylab=mod.descrip,
         xlab="Path Index",
         main=bquote("Repeated Minor Loop: H"["R"] ~ "=" ~ .(repeated[r]) ~ "(Oe)"),
         sub="Line segments display repeated measurements at fixed field")
    grid()
    for(P in 1:length(path.col)){
      p<-matching.paths[r,P]
      t.dat4<-t.dat3[t.dat3[,"path"]==p,]
      for(ff in 1:length(H_use)){
        f<-H_use[ff]
        t.ind<-t.dat3[,"path"]==p&(t.dat3[,"rnd.H"]==f)
        yy<-t.resid[t.ind]
        xx<-t.dat3[t.ind,"path"]
        lines(xx+(1:length(yy))*diff(XLIM)/(75*length(yy)),yy,col=field.col[ff])
      }
    }
    legend("top",title="Field",legend=H_use,ncol=5,fill=field.col)
  }
}

dev.off()





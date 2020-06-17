pdf(file=paste0(DIR,"Results/",File.name,"/",File.name," Experiment Overview Plots.pdf"),height=12,width=12)
plot(H_r,as.numeric(table(path)),ylab="Number of Observations in Minor Loop",
     main=paste("Basic Applied Field Stats for",File.name), lwd=3,
     col=colorRampPalette(c("purple","blue","light green", "orange", "red"))(51)[round((unique(path)-1)/(max(path)-1)*50)+1],)
grid()

color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
  scale = (length(lut)-1)/(max-min)
  plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main='')
  mtext(title,side=3,line=0)
  axis(2, pretty(ticks), las=1)
  for (i in 1:(length(lut)-1)) {
    y = (i-1)/scale + min
    rect(0,y,10,y+1/scale, col=lut[i], border=NA)
  }
}
subplot(color.bar(colorRampPalette(c("purple","blue","light green",  "orange", "red"))(51),1,max(path),nticks=5,title="Path Index"),x=sum(range(H_r)*c(.8,.2)),y=sum(range(table(path))*c(.8,.2)),size=c(1,2))

legend("topright",legend=c(paste("Median Time to Last Measurement:",round(med.max.time,2),"(s)"),
                           paste("# of Loops:",max(path)),
                           paste("Minor Loop direction:",if(BACKWARD==1) "Decreasing" else "Increasing", "Field"),
                           paste("Range of App. Field: (", round(min(dat[,"H"])),", ",round(max(dat[,"H"])),")")))

tmp2<-subplot(hist(n,main="",
                   breaks=0:max(n)+.5,col=2,
                   xlab="",
                   ylab="",
                   yaxt="n",
                   cex.main=.7,cex.lab=.8),
              x=sum(range(H_r)*c(.1,.9)),
              y=sum(range(table(path))*c(.16,.84)),size=c(2,1))

op <- par(no.readonly=TRUE)
par(tmp2)
mtext("# of Observations at Field",1,line=2,cex=.7)
par(op)

tmp2<-subplot(hist(table(rnd.H_r),main="", col=2, xlab="", ylab="",
                   yaxt="n", cex.main=.7,cex.lab=.8),
              x=sum(range(H_r)*c(.1,.9)), y=sum(range(table(path))*c(.5,.5)),size=c(2,1))

op <- par(no.readonly=TRUE)
par(tmp2)
mtext("# of Times Loop Repeats",1,line=2,cex=.7)
par(op)


plot(dat[rev.pt==1,1],H_r,ylab="H_r",xlab="Minor Loop Starting Time")
grid()


## To control file size, down sample if needed    
if(nrow(dat)>5000) ind<-round(1:5000/5000*nrow(dat)) else ind<-1:nrow(dat)
plot(dat[ind,1],dat[ind,2],xlab="Time Stamp (sec)",ylab="Temperature (K)")
grid()

plot(dat[ind,"rep.ind"],dat[ind,"time_diff"],
     xlab="Measurement Repeat Index",
     ylab="Time Since First Measurement (s)")
grid()

sub<-sample(nrow(dat),1000)
plot(dat[sub,"Time"], dat[sub,"resid"],
     pch=16,col=makeTransparent(1,50),
     main="Large Time Scale Temporal Drift",
     ylab="Residual (emu)",
     xlab="Time Stamp (s)")
trend<-smooth.spline(x=dat[,"Time"],y=dat[,"resid"],df=5)
lines(trend$x,trend$y,col=2,lwd=2)
grid()

gap<-max(spot.mn.dat[spot.mn.dat[,"path_step"]<11,"del_mom"])/2
n.trace<-floor(max(spot.mn.dat[,"path"])/5)
YLIM<-range(spot.mn.dat[spot.mn.dat[,"path_step"]<11,"del_mom"])+c(0,n.trace*gap)
plot(NULL,NULL,xlim=c(0,10),#range(dat[,"H"]),
     ylim=YLIM,
     yaxt="n",
     main="Minor Loop Beginnings",
     ylab="Moment Trajectory at Initial Steps of Minor Loop",
     xlab="Index within minor loop")
for(i in 1:n.trace)
  abline(h=(i-1)*gap,col="lightgray",lty=2)
COL<-colorRampPalette(c("purple","blue","light green", "orange", "red"))(floor(max(path)/5+1))

for(i in 1:n.trace){
  t.H_r<-sort(unique(as.numeric(spot.mn.dat[,"H_r"])))[5*i]
  t.path<-unique(spot.mn.dat[spot.mn.dat[,"H_r"]==t.H_r,"path"])[1]
  t.dat<-spot.mn.dat[spot.mn.dat[,"path"]==t.path,]
  lines(t.dat[,"path_step"],t.dat[,"del_mom"]+(i-1)*gap,lwd=2,
        col=COL[round((t.dat[1,"H_r"]-min(dat[,"H_r"]))/(diff(range(dat[,"H_r"])))*floor(max(path)/5+1))]  )
}

subplot(color.bar(
  colorRampPalette(c("purple","blue","light green",  "orange", "red"))(floor(max(path)/5+1)),
  min(dat[,"H_r"]),max(dat[,"H_r"]),nticks=5,title="H_R"),
  x=8,y=sum(YLIM*c(.04,.96)),size=c(1,1))

plot(dat[sub,"path_step"], dat[sub,"resid"],
     pch=16,col=makeTransparent(1,50),
     main="Minor Loop Beginnings",
     ylab="Residual (emu)",
     xlab="Field Index within Minor Loop")
trend<-smooth.spline(x=dat[,"path_step"],y=dat[,"resid"],df=5)
lines(trend$x,trend$y,col=2,lwd=2)
grid()

COL<-colorRampPalette(c("purple","blue","light green", "orange", "red"))(21)

YLIM=range(dat[sub,"resid"])
plot(dat[sub,"H_r"], dat[sub,"resid"],
     pch=16,col=COL[round(21*(dat[sub,"H_r_lag1"]-min(dat[,"H_r"]))/diff(range(dat[,"H_r"])))],
     main="Minor Loop Beginnings",
     ylab="Residual (emu)",
     xlab="Reversal Field")
subplot(color.bar(COL,
  min(dat[,"H_r"]),max(dat[,"H_r"]),nticks=5,title="H_R Lag1"),
  x=8,y=sum(YLIM*c(.04,.96)),size=c(1,1))
grid()
dev.off()


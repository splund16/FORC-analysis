## How many repeats are present in at least 20% of fields
n.x <- table(dat[,"rep.ind"])
n.x<-max(which(n.x>=(.2*n.x[1])))

pdf(file=paste(paste(getwd(),"Results",File.name,sep="/"),"/",File.name," Repeated Measurement Plots.pdf",sep=""),height=12,width=12)
## To control file size, down sample if needed    
if(nrow(dat)>5000) ind<-round(1:5000/5000*nrow(dat)) else ind<-1:nrow(dat)

rep.diff.list<-list()
rep.ratio.list<-list()
for(i in 1:(n.x-1)){
  ind<-which(dat[,"rep.ind"]==(i+1))
  rep.diff.list[[i]]<-dat[ind,Moment]-dat[ind-1,Moment]
  rep.ratio.list[[i]]<-dat[ind,Moment]/dat[ind-1,Moment]
  names(rep.diff.list[[i]])<-names(rep.diff.list[[i]])<-ind
}

par(mai=c(1,1,1,.02))

plot(-1,-1,xlim=range(time.steps),ylim=range(c(sapply(rep.diff.list,function(x)mean(x>0)),.5)),ylab=expression("Proportion Mom"["i+1"]*">Mom"["i"]),xlab="~Observation Time i (s)",
     main="Sign test for temporal trends",cex.lab=2,cex.axis=2)
grid()
for(i in 1:(n.x-1)){
  n<-length(rep.diff.list[[i]])
  p<-mean(rep.diff.list[[i]]>0)
  points(time.steps[i],p)
  lines(rep(time.steps[i],2),p+c(-2,2)*sqrt(p*(1-p)/n))
}
abline(h=.5,col="orange")

plot(-1,-1,xlim=range(time.steps),
     ylim=quantile(log(unlist(rep.ratio.list)),na.rm=TRUE,c(.1,.9)),
     ylab=expression("log2(Mom"["i+1"]/"Mom"["i"]*")"),xlab="~Observation Time i (s)",
     main="Size of temporal trends",cex.lab=2,cex.axis=2)
grid()
for(i in 1:(n.x-1))
  boxplot(log2(rep.ratio.list[[i]]),add=TRUE,at=time.steps[i],yaxt="n")


plot(-1,-1,xlim=range(time.steps),
     ylim=quantile(unlist(rep.diff.list),na.rm=TRUE,c(.1,.9)),
     ylab=expression("Mom"["i+1"]-"Mom"["i"]),xlab="~Observation Time i (s)",
     main="Size of temporal trends",cex.lab=2,cex.axis=2)
grid()
for(i in 1:(n.x-1)){
  boxplot(rep.diff.list[[i]],add=TRUE,at=time.steps[i],yaxt="n")
}

## Which h-reverse have increased variance?
rowvar<-unlist(as.data.frame(dat)%>%group_by(spot)%>%summarise(var.mom=var(Moment))%>%dplyr::select(var.mom))
rowvar<-rowvar[!is.na(rowvar)]
## Which applied fields have increased variance?

## Examine mean-variance relationship
use<-which(spot.change==1)[!is.na(rowvar)]
use2<-1:length(use)

if(length(use2)>5000)  use2<-use2[round((1:5000)/5000*length(use2))]

use<-use[use2]
plot(dat[use,"H"],log10(rowvar[use2]),xlab="Applied field",ylab="log10(Moment Variance)",cex.lab=2,cex.axis=2,
     col=colorRampPalette(c("purple","blue","light green", "orange", "red"))(51)[round((dat[use,"H_r"]-min(H_r))/(max(H_r)-min(H_r))*50)+1],
     lwd=1.5)
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

subplot(color.bar(colorRampPalette(c("purple","blue","light green",  "orange", "red"))(51),min(H_r),max(H_r),
                  nticks=5,title=expression("H"["R"])),
        x=sum(range(dat[use,"H"])*c(.85,.15)),
        y=sum(range(log10(rowvar[use2]),na.rm=TRUE)*c(.15,.85)),size=c(.7,2))


## Fit spline to identify region of increased variability in loops with low reversal fields
var.fit<-smooth.spline(x=dat[use,"H"][dat[use,"H_r"]< -40],y=log(rowvar[use2][dat[use,"H_r"]< -40])) 
if (n.x > 49) {
  nc <- 7
  nr <- 7
} else {
  nc <- ceiling(sqrt(n.x+1))
  nr <- ceiling((n.x+1)/nc)
}
H_r.cut<- -50
### Examine which field points are slow to approach equilibrium    
for(Hr.filter in 1:3){
  for (i in 1:length(rep.diff.list)) {
    if (i %in% (1 + nr * nc * (0:10))) {
      nf <- layout(matrix(1:(nc * nr), nr, nc, byrow = TRUE), 
                   widths = c(1.1 + 0.05 * nc, rep(1, nc - 1)), 
                   heights = c(rep(1, nr - 1), 1.1 + 0.05 * nr))
    }
    b <- l <- t <- r <- 0
    ylab <- xlab <- ""
    yaxt <- xaxt <- "n"
    if (i %in% (1 + nc * (0:nr))) {
      l <- 1
      ylab = expression(Delta*"Mom. ("*10^-5*"emu)")
      yaxt = "s"
    }
    if (i > n.x-nc) {
      b <- 1
      xlab = "Applied Field"
      xaxt = "s"
    }
    
    if(i<= nc) t<-.2
    par(mai = c(b, l, t, r),bty="o")
    Col <- makeTransparent(2)
    use<-1:length(rep.diff.list[[i]])
    t.H_r<-dat[as.numeric(names(rep.diff.list[[i]])),"H_r"]
    t.H<-dat[as.numeric(names(rep.diff.list[[i]])),"H"]
    t.path<-dat[as.numeric(names(rep.diff.list[[i]])),"path"]
    if(Hr.filter%%3==2) use<-use[t.H_r[use]<H_r.cut]
    if(Hr.filter%%3==0) use<-use[t.H_r[use]>H_r.cut]
    if(length(use>2000))  use<-use[round(1:2000)/2000*length(use)]
    Xlim<-range(var.fit$x[var.fit$y>(.7*median(var.fit$y)+.3*max(var.fit$y))])
    Xlim<-Xlim+diff(Xlim)*c(-.15,.15)
    
    
    Ylim=quantile(unlist(rep.diff.list)*10^5,c(.001,.975),na.rm=TRUE)
    suppressWarnings(plot(t.H[use], (rep.diff.list[[i]][use])*10^5, xlab = xlab, xaxt = xaxt, ylab = ylab, 
                          yaxt = yaxt, ylim=Ylim,xlim=Xlim, main = "", cex.lab = 1.7, 
                          cex.axis = 1.7, col = makeTransparent(colorRampPalette(c("purple","blue","light green", "orange", "red"))(101)[round((t.H_r[use]-min(t.H_r[use]))/(max(t.H_r[use])-min(t.H_r[use]))*100)+1],200),
                          pch=(ceiling(t.path[use]/max(t.path)*5))
    ))
    grid()
    text(Xlim[1],Ylim[1],paste(time.steps[i],"s vs. ",time.steps[i+1],"s",sep=""),adj=c(0,0),cex=1.5)
  }
  
  if(i==length(rep.diff.list)){
    par(mai=c(.5*b,0,0,r),bty="n")
    plot(0,0,col="white",xaxt="n",yaxt="n",xlab="",ylab="",ylim=Ylim,xlim=Xlim)
    path.bin<-round((0:5)/5*max(t.path))
    legend("topright",title="Path",pch=1:5,legend=paste(path.bin[-length(path.bin)],path.bin[-1],sep="-"),cex=1.4,bty="n")
    subplot(color.bar(colorRampPalette(c("purple","blue","light green",  "orange", "red"))(101),
                      min(t.H_r[use]),max(t.H_r[use]),nticks=5,title=expression("H"["R"])),
            x=sum(Xlim*c(.75,.25)),y=sum(Ylim*c(.75,.25)),size=c(.6,1.5))
  }
  par(bty="o")
}

dev.off()
if(0){
  
  ### This plot is only intended for the file(s) in which two separate loops appeared
  par(mfrow=c(1,1),mar=c(5, 4, 4, 2) + 0.1)
  ind2<-(t.H_r<(-99))&(abs(t.applied-85)<25.5)
  t.mat2<-t.mat[ind2,]
  t.H_r2<-t.H_r[ind2]
  t.path2<-t.path[ind2]
  t.field2<-t.field[ind2]
  plot(t.path2,t.mat2[,2]-t.mat2[,1],
       main="Path Index Relationship to Inner vs. Outer Loop", xlab="Path Index",
       ylab=expression("Mom."[t[2]]-"Mom."[t[1]]*" ("*10^-5*"emu)"),
       pch=ceiling((t.field[ind2]-min(t.field[ind2]))/(max(t.field[ind2])-min(t.field[ind2]))*5),
       col = makeTransparent(colorRampPalette(c("purple","blue","light green", "orange", "red"))(101)[round((t.H_r2-min(t.H_r2))/(max(t.H_r2)-min(t.H_r2))*100)+1],200))
  subplot(color.bar(colorRampPalette(c("purple","blue","light green",  "orange", "red"))(101),min(t.H_r2),max(t.H_r2),nticks=5,title=expression("H"["R"]*"(Oe)")),
          x=.1*max(t.path2),
          y=sum(range(t.mat2[,2]-t.mat2[,1])*c(.9,.1)),size=c(.6,1.5))
  breaks<-round((min(t.field2)*(5:0)+max(t.field2)*(0:5))/5)
  legend("topleft",legend=paste(breaks[-6],breaks[-1],sep="-"),pch=1:5,title=expression(H[applied]*"(Oe)"),bg="white")
  text(mean(range(t.path2)),sum(range(t.mat2[,2]-t.mat2[,1])*c(.05,.95)),"This plot is primarily useful when two separate loops appear in preceding moment drift plots")
  
  grid()    
  spline.fit<-smooth.spline(x=t.path2,y=t.mat2[,2]-t.mat2[,1])
  lines(spline.fit,lwd=4)
  dev.off()
}
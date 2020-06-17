require(Hmisc)
require(dplyr)
if(!file.exists(paste("Results",File.name,sep='/'))) dir.create(paste("Results",File.name,sep='/'))
DIR<-getwd()
if(substr(DIR,nchar(DIR),nchar(DIR))!="/") DIR<-paste0(DIR,"/")
##Load data
### Find start of data
where<-readLines(paste(DIR,File.name,".dat",sep=""),n=80)
SKIP<-which(where=="[Data]")+1
names<-scan(paste(DIR,File.name,".dat",sep=""),sep=",",skip=SKIP-1,what="character",nlines=1,quote="")
dat<-scan(paste(DIR,File.name,".dat",sep=""),sep=",",skip=SKIP)
dat<-matrix(dat,ncol=length(names),byrow=TRUE)
colnames(dat)<-names
n.unique<-apply(dat,2,function(v) length(unique(v)))
dat<-dat[,n.unique>1]
colnames(dat)[which(colnames(dat)==Moment)]<-"Moment"
Moment<-"Moment"
colnames(dat)[which(colnames(dat)==MagField)]<-"H"
H<-"H"

colnames(dat)[which(colnames(dat)==MomSE)]<-"Moment_SE"
dat[,"H"]<-dat[,"H"]/10

### Identify distinct paths by reversal points
rev.pt<-c(1,dat[-1,H]<dat[-nrow(dat),H])
if(BACKWARD==1)rev.pt<-c(1,dat[-1,H]>dat[-nrow(dat),H])
H_r<-dat[rev.pt==1,H]
path<-cumsum(rev.pt)
dat<-cbind(dat,path=path,H_r=H_r[path],
           H_r_lag1=c(NA,H_r[-length(H_r)])[path],
           H_r_lag2=c(NA,NA,H_r[-(length(H_r)+1:0)])[path])


dat<-cbind(dat,rnd.H=round(dat[,H]/2.5)*2.5)
rnd.H<-"rnd.H"
rnd.H_r<-dat[rev.pt==1,rnd.H]
dat<-cbind(dat,rnd.H_r=rnd.H_r[path])

spot.lab<-paste(dat[,"path"],dat[,"rnd.H_r"],dat[,"rnd.H"])
n.rep<-table(spot.lab)
rep.ind<-NULL
spot<-0
for(i in 1:length(n.rep)){
  t.lab<-spot.lab[spot+1]
  rep.ind<-c(rep.ind,1:n.rep[t.lab])
  spot<-length(rep.ind)
}
#dat[1:30,c("path",H,"H_r","rep.ind")]
spot.change<-c(1,dat[-1,H]!=dat[-nrow(dat),H])
if(any((rep.ind==1)!=(spot.change==1))){
  rep.ind<-spot.change<-rep(NA,nrow(dat))
  warning("Replicate measurement index 'rep.ind' has failed.")
} 
dat<-cbind(dat,rep.ind=rep.ind,spot=cumsum(spot.change),spot.change=spot.change)

### Include time elapsed since first measurement replicate
colnames(dat)[1]<-"Time"
dat<-as.data.frame(dat)%>%
  group_by(path,rnd.H_r,rnd.H)%>%
  mutate(time_diff=Time-min(Time))%>%
  group_by(rep.ind)%>%
  mutate(med.time_diff=median(time_diff))
dat<-as.matrix(dat)

time.steps<-signif(sort(unique(dat[,'med.time_diff'])),3)

### median of time difference between first and last replicate measurements
med.max.time<-as.data.frame(dat)%>%
  group_by(path,rnd.H_r,rnd.H)%>%
  summarise(max_time=max(time_diff))
med.max.time<-median(med.max.time$max_time)

### Create index used repeated observations at given field within minor loops
#t.id<-paste(path,dat[,H],rep(H_r,table(path)))

### Estimate average time duration of measurements
#time.den<-density(diff(dat[,1]))
#hist(diff(dat[,1])[diff(dat[,1])<40])
#plot(sort(diff(dat[,1])),log="y")
#mode.time.step<-round(time.den$x[which.max(time.den$y)],1)

### How many repeated measurements are taken at a step
n<-dat[c(which(dat[,"spot.change"]==1)-1,nrow(dat)),"rep.ind"]
#un.id<-unique(t.id)

#par(mar=c(5, 4, 4, 2) + 0.1)
### Display basic experiment description


if(0){#max(n)>1){
  print("Reformatting data matrix so repeated observations from same field and path are contained in single row")
  repeated.measurements<-TRUE
  
  if(paste0(gsub(".dat","",File.name),"Reform.Rdata")%in%list.files(paste0(DIR,"Results/",File.name))){
    load(paste0(DIR,"Results/",File.name,"/",gsub(".dat","",File.name),"Reform.Rdata"))
  }else{
    ## Create matrix where repeated observations from same field and path are contained in single row
    t.mat<-matrix(NA,length(un.id),max(n))
    rownames(t.mat)<-un.id
    j<-0
    print("Progress out of:")
    for(i in 1:19) cat("X")
    cat(" Done!\n")
    for(i in 1:length(un.id)){
      if(i%in%round((1:19)*length(un.id)/20))cat("X")
      tt.id<-t.id[j+1]
      t.mat[i,1:n[tt.id]]<-dat[j+1:n[tt.id],Moment]
      j<-j+n[tt.id]
    }
    colnames(t.mat)<-paste(Moment,signif(time.steps,2))
    
    ## Track row indices in original data matrix (Useful to get values of other variables)
    step.ind<-cumsum(c(0,n[un.id]))+1
    step.ind<-step.ind[-length(step.ind)]
    #Make sure it's the correct mapping by comparing values of first observations
    if(any(dat[step.ind,Moment]!=t.mat[,1])) stop("Error in formatting repeated observations")
    ## Add values from all other variables (at first observation)
    reform.dat<-cbind(dat[step.ind,],t.mat)
    dim(reform.dat)
    head(reform.dat)
    ## Save the restructed dataset
    save(reform.dat,file=paste0(DIR,"Results/",File.name,"/",gsub(".dat","",File.name),"Reform.Rdata"))
  }
}
require(dplyr)

step<-as.data.frame(dat)%>%group_by(path)%>%mutate(path_step=spot-min(spot))
step<-step$path_step
dat<-cbind(dat,path_step=step)

###Average across repeated minor loops and replicate measurements
ave.dat<-as.data.frame(dat[,c(H,Moment,"rnd.H","rnd.H_r","H_r","Moment_SE")])%>%
  group_by(rnd.H_r)%>%
  mutate(H_r=mean(H_r))%>%
  group_by(H_r,rnd.H)%>%
  summarise(rnd.H_r=2.5*round(mean(H_r)/2.5),H=mean(H),Moment=mean(Moment),n=length(Moment),Moment_SE=mean(Moment_SE))%>%
  arrange(H_r,H)
mn.rev.pt<-c(1,ave.dat[-1,H]<ave.dat[-nrow(ave.dat),H])
if(BACKWARD==1) mn.rev.pt<-c(1,ave.dat[-1,H]>ave.dat[-nrow(ave.dat),H])
mn.path<-cumsum(mn.rev.pt)
ave.dat<-data.frame(ave.dat,rev.pt=mn.rev.pt,path=mn.path)
ave.dat$H_b=(ave.dat$H_r+ave.dat$H)/2
ave.dat$H_c=(ave.dat$H-ave.dat$H_r)/2

spot.mn.dat<-as.data.frame(dat[,c("Time",H,Moment,"rnd.H","rnd.H_r","path","H_r","path_step","H_r_lag1","H_r_lag2")])
spot.mn.dat<-spot.mn.dat%>%
  group_by(path,rnd.H_r,rnd.H)%>%
  summarise(Time=mean(Time),H=mean(H),Moment=mean(Moment),
            n=length(Moment),H_r=mean(H_r),
            H_r_lag1=mean(H_r_lag1),H_r_lag2=mean(H_r_lag2),
            path_step=mean(path_step))
spot.mn.dat<-spot.mn.dat%>%group_by(path)%>%mutate(del_mom=Moment-Moment[1])

mn.rev.pt<-c(1,spot.mn.dat[-1,H]<spot.mn.dat[-nrow(spot.mn.dat),H])
if(BACKWARD==1) mn.rev.pt<-c(1,spot.mn.dat[-1,H]>spot.mn.dat[-nrow(spot.mn.dat),H])
spot.mn.dat<-cbind(as.matrix(spot.mn.dat),rev.pt=mn.rev.pt)


makeTransparent<-function(someColor, alpha = 100) {
  newColor <- col2rgb(someColor)
  apply(newColor, 2, function(curcoldata) {
    rgb(red = curcoldata[1], green = curcoldata[2], blue = curcoldata[3], alpha = alpha, 
        maxColorValue = 255)
  })
}

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


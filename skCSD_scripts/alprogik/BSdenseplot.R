source("alprogik/Colors_BlueRed.R")
library(data.table)
library(RColorBrewer)
library(fields)

M<-512
lambda.all<-c("1e-05","1e-04","0.001","0.01","0.1" ) #c(0.1,0.01,0.001)#
R.all<-2^(3:7) #c(8,16,128)#
ScriptLocation<-"/home/csdori/ksCSD_2014/trunk/"
SimPath<-paste0(ScriptLocation,"simulation/")
ElNumb<-c(8,16, 128)
#ElNumb<-c(10,20,50,100)
locationsKernel<-c("BS_d50_el8_CosChanging",  "BS_d50_el16_CosChanging","BS_d50_el128_CosChanging")#paste0("BS_d50_el",ElNumb) 
#locationsKernel<-c(locationsKernel,paste0("Y_el",ElNumb,"_m200_600"))
locationsData<-locationsKernel #c("BS_d50_el8_GaussChanging190", "BS_d50_el128_GaussChanging190")#paste0("BS_d50_el",ElNumb,"_Gauss10")  #locationsKernel
outname1<-"/kernelOut_4"
#outname1<-"/kernelOut_poli_CP3"
#lapply(paste0(SimPath,locationsData, outname1), dir.create)


MinError <- numeric()
MinErrorR<-numeric()
MinErrorL<-numeric()



#ImageAllParCurrReg(SimPath,locationsKernel,locationsData, 0,c(8,32,128),c(0.1,1e-5))


#ImageAllParCurrReg<-function(SimPath,locationsKernel,locationsData,SNR,  R2Plot,lambda2Plot){



#png(paste0(SumPlotPlace,"BSnarrow.png"),width=1000, height=900, pointsize=22)
tiff(paste0(SumPlotPlace,"/BS_Cosine.tiff"),width=10, height=9, res=400, unit="in" )
plot.new()
par(mfrow=c(3,2),oma=c(4,5,4,4))#,mar=c(1,1,1,1))
#par(mfrow=c(length(lambda.all)+2,length(R.all)),oma=c(2,3,4,2))#,mar=c(1,1,1,1))

#layout(rbind(c(1,2,3,4,5)),widths=c(0.5,0.2,1,3,0.2),respect=FALSE) 
#mtext(c("oma[SOUTH<-1]", "oma[WEST<-2]", "oma[NORTH<-3]", "oma[EAST<-4]"),     c(SOUTH<-1, WEST<-2, NORTH<-3, EAST<-4),      line=0.4, cex=1.2, col="green", outer=TRUE)
#title(main="skCSD Reconstruction", outer=TRUE, cex=3)




for(Dlength in 1:length(locationsKernel)){
inname<-paste0(SimPath, locationsKernel[Dlength],"/")
outname<-paste0(SimPath, locationsData[Dlength])
OrigCurr<-as.matrix(read.table(paste(SimPath,"/",locationsData[Dlength],"/membcurr",sep="")))
 SameplaceAll<-readLines(paste0(inname,outname1,"/sameplace.txt"))
 coord<-as.matrix(read.table((paste0(inname,"/coordsmid_x_y_z"))))[105:156]
    seg.nb<-max(as.matrix(read.table(paste(inname,outname1,'/connection_touched_once',sep=''))))

seg.length<-as.matrix(read.table(paste(SimPath,"/",locationsData[Dlength],"/","seglength",sep="")))
funaramvonal<-function(x) x/seg.length
memb.currents.line1<-array(0, c(dim(OrigCurr)))
memb.currents.line1<-apply(OrigCurr,2,funaramvonal) 
#memb.currents.line1[which(memb.currents.line1<(0.5*min(memb.currents.line1)))]<-0.5*min(memb.currents.line1)
times2Plot<-1:10 #as.matrix(read.table(paste0(inname,"time")))







BesteEmatr<-as.matrix(read.table(paste(SimPath,"/",locationsData[Dlength],"/kernelOut_4/EL1ErrorVarR_SNR0",sep="")))
CurrMin<- min(BesteEmatr[26:30,])
IndecesMin<-which(BesteEmatr==CurrMin,arr.ind=TRUE)
MinErrorR<-c(MinErrorR,R.all[IndecesMin[1]-25])
MinErrorL<-c(MinErrorL,lambda.all[IndecesMin[2]])


R2Plot<- MinErrorR[Dlength]
lambda2Plot<- MinErrorL[Dlength]

#if(file.exists(paste(inname,outname1,"/params",sep=""))==FALSE) next
#params<-as.matrix(read.table(paste(inname,outname1,"/params",sep="")))

#memb.currents.line1[which(memb.currents.line1<(0.7*max(memb.currents.line1)))]<-0.5*max(memb.currents.line1)
#OrigCurrS<-as.matrix(read.table(paste(outname,outname1,"/membcurr_smoothed",sep="")))
sCSD<-as.matrix(read.table(paste0(outname,"/sCSDginv")))
sCSD<-sCSD/(700/ElNumb[Dlength])

#Select just few Rs and Ms for plotting
R.allplot<-R2Plot#seq(5,130,25)
M<-512
lambda.allplot<-lambda2Plot#c("1e-05","1e-04","0.001","0.01","0.1" )



sCSDCoord<-as.matrix(read.table(paste0(outname,"/sCSDCoord")))


if(Dlength==1){
 par(mfg=c(1,1))
col2<-ColoursDori(memb.currents.line1)[[1]]
ExtVal<-0.11 #ColoursDori(memb.currents.line1)[[2]]
par(mar=c(1,0.1,0.1,0.1))
SegNb<-dim(memb.currents.line1)[1]
image(times2Plot, sort(coord),t(memb.currents.line1[,times2Plot]),col=col2,zlim=c(-ExtVal,ExtVal),xaxt="n", ylab="z (um)")
mtext("Ground Truth",3)
mtext("z (um)",side=2,line=2.5)
axis(1,round(seq(0,max(times2Plot)-10,,5)),las=2)
#axis(4,round(seq(1,SegNb-10,,4)),las=2)
mtext("Segment ID",3,line=5)
col2<-ColoursDori(memb.currents.line1)[[1]]
ExtVal<-ColoursDori(memb.currents.line1)[[2]]

par(mfg=c(1,2),mar=c(0.3,0.1,0.1,19)) 
image.plot(times2Plot,1:SegNb,t(memb.currents.line1[,times2Plot]),col=col2,zlim=c(-ExtVal,ExtVal),xaxt="n", yaxt="n",axis.args=list( at=c(-ExtVal, ExtVal), labels=c("-0.11 nA/um", "0.11 nA/um")),
 legend.only=TRUE)

}
par(mar=c(0,8,0,8)) 
#if(R2Plot<4) par(mar=c(0,8,0,8)) 
#if(R2Plot>3) par(mar=c(0,4,0,4)) 


#if(CellType=="BS"){
par(mfg=c(Dlength+1,1))
par(mar=c(1,0.1,1,0.1))
col2<-ColoursDori(sCSD)[[1]]
ExtVal<-0.11#ColoursDori(sCSD)[[2]]
image(times2Plot,sCSDCoord[,3],t(sCSD[,times2Plot]),col=col2,zlim=c(-ExtVal,ExtVal),xaxt="n",ylab="z (um)")
mtext("sCSD",3)
mtext(paste0(ElNumb[Dlength]," Electrodes"),side=2,line=4)
mtext("z (um)",side=2,line=2.5)
axis(1,round(seq(0,max(times2Plot)-10,,5)),las=2)
if(Dlength==2){axis(1,at=c(1:10),paste(seq(10,190,20))) 
mtext("R",1,line=3)
mtext(c( "Segment ID"), 4, line=2.5)
}

#par(mar=c(0,3,0,5)) 
#}


#image(c(times2Plot),1:dim(memb.currents.line1)[1],array(0,c(length(times2Plot),dim(memb.currents.line1)[1])),xlab="Time",ylab="Segment #",col="WHITE")
# image(t(OrigCurrS),col=rainbow(500),main="Smoothed Ground Truth")
  par(mfg=c(Dlength+1,2))

    skCSD.all<-array(0,c(seg.nb, length(times2Plot)))
#for(m in 1:length(lambda.all)){
 # for(r in 1:length(R.all)){
    m<-1
    r<-1
    Lambda<-lambda.allplot[m]
    R<-R.allplot[r]
    #Reading in the LFP and adding noise
    #Ktilda_M512_R55lambda1e-04
currName<-paste0(outname,outname1,"/skCSDall_M",M,"_R",R,"lambda",Lambda)
    if(file.exists(currName)) {
      skCSD.all.part<-as.matrix(read.table(currName)) 
      

 for(i in 1: seg.nb){
      sameplace<-numeric() 
      sameplace<-as.integer(unlist(strsplit(SameplaceAll," ")[i]))
      #Plot the current related ti the different segments
    #  matplot(cbind(t(C.calc[c(sameplace),])),t='l', xlab="Time steps", ylab="ksCSD", main="Values of current at the same location")
      if(length(sameplace)==1) skCSD.all[i,]<- skCSD.all.part[sameplace,times2Plot]
      if(length(sameplace)>1) skCSD.all[i,]<-as.matrix(colSums( skCSD.all.part[c(sameplace),times2Plot]),nrow=1)
    }      
      
      
     
                               }else next
col2<-ColoursDori(skCSD.all)[[1]]
ExtVal<-0.11#ColoursDori(skCSD.all)[[2]]
par(mar=c(1,0.1,1,0.1))
image(times2Plot,1:SegNb,t(skCSD.all[,times2Plot]),col=col2,zlim=c(-ExtVal,ExtVal),xaxt="n",yaxt="n")
if(m==1) mtext(paste0("skCSD"),side=3)#mtext(paste0(ElNumb[Dlength]," Electrodes"),side=3)
#if(r==1) mtext(paste0("Lambda", Lambda),side=2)
if(r==length(R.allplot)) {axis(4,at=round(seq(1,SegNb-10,,4)),las=2)}
if(m==length(lambda.allplot)) axis(1,round(seq(0,max(times2Plot)-10,,5)),las=2)
mtext(c( "Segment ID"), 4, line=2.5)
if(Dlength==2){axis(1,at=c(1:10),paste(seq(10,190,20))) 

mtext("R",1,line=3)
}


#dev.off()
}


dev.off()



####################################
############################### Type 2
###########################################x




MinError <- numeric()
MinErrorR<-numeric()
MinErrorL<-numeric()



#ImageAllParCurrReg(SimPath,locationsKernel,locationsData, 0,c(8,32,128),c(0.1,1e-5))


#ImageAllParCurrReg<-function(SimPath,locationsKernel,locationsData,SNR,  R2Plot,lambda2Plot){



#png(paste0(SumPlotPlace,"BSnarrow.png"),width=1000, height=900, pointsize=22)
tiff(paste0(SumPlotPlace,"/BS_CosineMatplot.tiff"),width=10, height=9, res=400, unit="in" )
plot.new()
par(mfrow=c(5,1),oma=c(4,5,4,4))#,mar=c(1,1,1,1))
#par(mfrow=c(length(lambda.all)+2,length(R.all)),oma=c(2,3,4,2))#,mar=c(1,1,1,1))

#layout(rbind(c(1,2,3,4,5)),widths=c(0.5,0.2,1,3,0.2),respect=FALSE) 
#mtext(c("oma[SOUTH<-1]", "oma[WEST<-2]", "oma[NORTH<-3]", "oma[EAST<-4]"),     c(SOUTH<-1, WEST<-2, NORTH<-3, EAST<-4),      line=0.4, cex=1.2, col="green", outer=TRUE)
#title(main="skCSD Reconstruction", outer=TRUE, cex=3)




source(paste0(ScriptLocation,"alprogik/errorCalcFun.R"))

L1errors<-numeric()






for(Dlength in 1:length(locationsKernel)){
inname<-paste0(SimPath, locationsKernel[Dlength],"/")
outname<-paste0(SimPath, locationsData[Dlength])
OrigCurr<-as.matrix(read.table(paste(SimPath,"/",locationsData[Dlength],"/membcurr",sep="")))
 SameplaceAll<-readLines(paste0(inname,outname1,"/sameplace.txt"))
 coord<-as.matrix(read.table((paste0(inname,"/coordsmid_x_y_z"))))[105:156]
    seg.nb<-max(as.matrix(read.table(paste(inname,outname1,'/connection_touched_once',sep=''))))

seg.length<-as.matrix(read.table(paste(SimPath,"/",locationsData[Dlength],"/","seglength",sep="")))
funaramvonal<-function(x) x/seg.length
memb.currents.line1<-array(0, c(dim(OrigCurr)))
memb.currents.line1<-apply(OrigCurr,2,funaramvonal) 
#memb.currents.line1[which(memb.currents.line1<(0.5*min(memb.currents.line1)))]<-0.5*min(memb.currents.line1)
times2Plot<-1:10 #as.matrix(read.table(paste0(inname,"time")))

coordsMid<-matrix(as.matrix(read.table(paste0(outname,"/coordsmid_x_y_z"))),ncol=3)

BesteEmatr<-as.matrix(read.table(paste(SimPath,"/",locationsData[Dlength],"/kernelOut_4/EL1ErrorVarR_SNR0",sep="")))
CurrMin<- min(BesteEmatr[26:30,])
IndecesMin<-which(BesteEmatr==CurrMin,arr.ind=TRUE)
MinErrorR<-c(MinErrorR,R.all[IndecesMin[1]-25])
MinErrorL<-c(MinErrorL,lambda.all[IndecesMin[2]])


R2Plot<- MinErrorR[Dlength]
lambda2Plot<- MinErrorL[Dlength]



#Select just few Rs and Ms for plotting
R.allplot<-R2Plot#seq(5,130,25)
M<-512
lambda.allplot<-lambda2Plot#c("1e-05","1e-04","0.001","0.01","0.1" )



if(Dlength==1){
 par(mfg=c(1,1))
par(mar=c(1,0.1,1,0.1))

SegNb<-dim(memb.currents.line1)[1]

matplot(sort(coord),memb.currents.line1[c(2,1,3:52),times2Plot],t="l",ylim=c(-0.13,0.13),xaxt="n",main="Ground Truth",ylab="CSD (nA/um)")

}
#par(mar=c(0,8,0,8)) 
#if(R2Plot<4) par(mar=c(0,8,0,8)) 
#if(R2Plot>3) par(mar=c(0,4,0,4)) 


  par(mfg=c(Dlength+1,1))

    skCSD.all<-array(0,c(seg.nb, length(times2Plot)))
#for(m in 1:length(lambda.all)){
 # for(r in 1:length(R.all)){
    m<-1
    r<-1
    Lambda<-lambda.allplot[m]
    R<-R.allplot[r]
    #Reading in the LFP and adding noise
    #Ktilda_M512_R55lambda1e-04
currName<-paste0(outname,outname1,"/skCSDall_M",M,"_R",R,"lambda",Lambda)
    if(file.exists(currName)) {
      skCSD.all.part<-as.matrix(read.table(currName)) 
      

 for(i in 1: seg.nb){
      sameplace<-numeric() 
      sameplace<-as.integer(unlist(strsplit(SameplaceAll," ")[i]))
      #Plot the current related ti the different segments
    #  matplot(cbind(t(C.calc[c(sameplace),])),t='l', xlab="Time steps", ylab="ksCSD", main="Values of current at the same location")
      if(length(sameplace)==1) skCSD.all[i,]<- skCSD.all.part[sameplace,times2Plot]
      if(length(sameplace)>1) skCSD.all[i,]<-as.matrix(colSums( skCSD.all.part[c(sameplace),times2Plot]),nrow=1)
    }      
      
      
     
                               }else next


L1errors<-c(L1errors,L1Error(skCSD.all[,1:10], memb.currents.line1[,1:10])$lerrorInTime)

par(mar=c(1,0.1,1,0.1))
if(Dlength!=3) matplot(coordsMid[c(2,1,3:52),3],skCSD.all[c(2,1,3:52),times2Plot],t="l",ylim=c(-0.13,0.13),xaxt="n",main=paste(ElNumb[Dlength], "Electrodes"),ylab="CSD (nA/um)")

if(Dlength==3) matplot(coordsMid[c(2,1,3:52),3],skCSD.all[c(2,1,3:52),times2Plot],t="l",ylim=c(-0.13,0.13),xlab="x (um)",main=paste(ElNumb[Dlength], "Electrodes"),ylab="CSD (nA/um)")
}

L1errors<-matrix(L1errors,ncol=Dlength)
matplot(L1errors,t="l")
#dev.off()
}


dev.off()




source(paste0(ScriptLocation,"alprogik/errorCalcFun.R"))

t1<-L1Error(skCSD.all[,1:10], memb.currents.line1[,1:10])








tiff(paste0(SumPlotPlace,"/BS_CosineImagePl.tiff"),width=10, height=9, res=400, unit="in" )
plot.new()
par(mfrow=c(5,1),oma=c(4,5,4,4))#,mar=c(1,1,1,1))
#par(mfrow=c(length(lambda.all)+2,length(R.all)),oma=c(2,3,4,2))#,mar=c(1,1,1,1))

#layout(rbind(c(1,2,3,4,5)),widths=c(0.5,0.2,1,3,0.2),respect=FALSE) 
#mtext(c("oma[SOUTH<-1]", "oma[WEST<-2]", "oma[NORTH<-3]", "oma[EAST<-4]"),     c(SOUTH<-1, WEST<-2, NORTH<-3, EAST<-4),      line=0.4, cex=1.2, col="green", outer=TRUE)
#title(main="skCSD Reconstruction", outer=TRUE, cex=3)




source(paste0(ScriptLocation,"alprogik/errorCalcFun.R"))

L1errors<-numeric()






for(Dlength in 1:length(locationsKernel)){
  inname<-paste0(SimPath, locationsKernel[Dlength],"/")
  outname<-paste0(SimPath, locationsData[Dlength])
  OrigCurr<-as.matrix(read.table(paste(SimPath,"/",locationsData[Dlength],"/membcurr",sep="")))
  SameplaceAll<-readLines(paste0(inname,outname1,"/sameplace.txt"))
  coord<-as.matrix(read.table((paste0(inname,"/coordsmid_x_y_z"))))[105:156]
  seg.nb<-max(as.matrix(read.table(paste(inname,outname1,'/connection_touched_once',sep=''))))
  
  seg.length<-as.matrix(read.table(paste(SimPath,"/",locationsData[Dlength],"/","seglength",sep="")))
  funaramvonal<-function(x) x/seg.length
  memb.currents.line1<-array(0, c(dim(OrigCurr)))
  memb.currents.line1<-apply(OrigCurr,2,funaramvonal) 
  #memb.currents.line1[which(memb.currents.line1<(0.5*min(memb.currents.line1)))]<-0.5*min(memb.currents.line1)
  times2Plot<-seq(1,25,by=3)#1:10 #as.matrix(read.table(paste0(inname,"time")))
  
  coordsMid<-matrix(as.matrix(read.table(paste0(outname,"/coordsmid_x_y_z"))),ncol=3)
  
  BesteEmatr<-as.matrix(read.table(paste(SimPath,"/",locationsData[Dlength],"/kernelOut_4/EL1ErrorVarR_SNR0",sep="")))
  CurrMin<- min(BesteEmatr[26:30,])
  IndecesMin<-which(BesteEmatr==CurrMin,arr.ind=TRUE)
  MinErrorR<-c(MinErrorR,R.all[IndecesMin[1]-25])
  MinErrorL<-c(MinErrorL,lambda.all[IndecesMin[2]])
  
  
  R2Plot<- MinErrorR[Dlength]
  lambda2Plot<- MinErrorL[Dlength]
  
  
  
  #Select just few Rs and Ms for plotting
  R.allplot<-R2Plot#seq(5,130,25)
  M<-512
  lambda.allplot<-lambda2Plot#c("1e-05","1e-04","0.001","0.01","0.1" )
  
  
  
  if(Dlength==1){
    par(mfg=c(1,1))
    par(mar=c(1,0.1,1,0.1))
    
    SegNb<-dim(memb.currents.line1)[1]
    
    #matplot(sort(coord),memb.currents.line1[c(2,1,3:52),times2Plot],t="l",ylim=c(-0.13,0.13),xaxt="n",main="Ground Truth",ylab="CSD (nA/um)")
    
    col2<-ColoursDori(memb.currents.line1)[[1]]
    ExtVal<-0.11 #ColoursDori(memb.currents.line1)[[2]]
    par(mar=c(1,0.1,0.1,0.1))
    SegNb<-dim(memb.currents.line1)[1]
    image(times2Plot/2, sort(coord),t(memb.currents.line1[c(2,1,3:52),times2Plot]),col=col2,zlim=c(-ExtVal,ExtVal),xaxt="n", ylab="z (um)")
    mtext("Ground Truth",3)
    mtext("z (um)",side=2,line=2.5)
  #  axis(1,round(seq(0,max(times2Plot)-10,,5)),las=2)
    
    
    
  }
  #par(mar=c(0,8,0,8)) 
  #if(R2Plot<4) par(mar=c(0,8,0,8)) 
  #if(R2Plot>3) par(mar=c(0,4,0,4)) 
  
  
  par(mfg=c(Dlength+1,1))
  
  skCSD.all<-array(0,c(seg.nb, length(times2Plot)))
  #for(m in 1:length(lambda.all)){
  # for(r in 1:length(R.all)){
  m<-1
  r<-1
  Lambda<-lambda.allplot[m]
  R<-R.allplot[r]
  #Reading in the LFP and adding noise
  #Ktilda_M512_R55lambda1e-04
  currName<-paste0(outname,outname1,"/skCSDall_M",M,"_R",R,"lambda",Lambda)
  if(file.exists(currName)) {
    skCSD.all.part<-as.matrix(read.table(currName)) 
    
    
    for(i in 1: seg.nb){
      sameplace<-numeric() 
      sameplace<-as.integer(unlist(strsplit(SameplaceAll," ")[i]))
      #Plot the current related ti the different segments
      #  matplot(cbind(t(C.calc[c(sameplace),])),t='l', xlab="Time steps", ylab="ksCSD", main="Values of current at the same location")
      if(length(sameplace)==1) skCSD.all[i,]<- skCSD.all.part[sameplace,times2Plot]
      if(length(sameplace)>1) skCSD.all[i,]<-as.matrix(colSums( skCSD.all.part[c(sameplace),times2Plot]),nrow=1)
    }      
    
    
    
  }else next
  
  
  L1errors<-c(L1errors,L1Error(skCSD.all[,1:10], memb.currents.line1[,1:10])$lerrorInTime)
  
  par(mar=c(1,0.1,1,0.1))
  if(Dlength!=3) #matplot(coordsMid[c(2,1,3:52),3],skCSD.all[c(2,1,3:52),times2Plot],t="l",ylim=c(-0.13,0.13),xaxt="n",main=paste(ElNumb[Dlength], "Electrodes"),ylab="CSD (nA/um)")
  {
    col2<-ColoursDori(memb.currents.line1)[[1]]
    ExtVal<-0.11 #ColoursDori(memb.currents.line1)[[2]]
    par(mar=c(1,0.1,0.1,0.1))
    SegNb<-dim(memb.currents.line1)[1]
    image(times2Plot, sort(coord),t(skCSD.all[c(2,1,3:52),times2Plot]),col=col2,zlim=c(-ExtVal,ExtVal),xaxt="n", ylab="z (um)")
    mtext("Ground Truth",3)
    mtext("z (um)",side=2,line=2.5)
    axis(1,round(seq(0,max(times2Plot)-10,,5)),las=2)
    
    
  }
  if(Dlength==3) #matplot(coordsMid[c(2,1,3:52),3],skCSD.all[c(2,1,3:52),times2Plot],t="l",ylim=c(-0.13,0.13),xlab="x (um)",main=paste(ElNumb[Dlength], "Electrodes"),ylab="CSD (nA/um)")
  {
    col2<-ColoursDori(memb.currents.line1)[[1]]
    ExtVal<-0.11 #ColoursDori(memb.currents.line1)[[2]]
    par(mar=c(1,0.1,0.1,0.1))
    SegNb<-dim(memb.currents.line1)[1]
    image(times2Plot, sort(coord),t(skCSD.all[c(2,1,3:52),times2Plot]),col=col2,zlim=c(-ExtVal,ExtVal),xaxt="n", ylab="z (um)")
    image.plot(times2Plot, sort(coord),t(skCSD.all[c(2,1,3:52),times2Plot]),col=col2,zlim=c(-ExtVal,ExtVal),xaxt="n", ylab="z (um)",
      legen.only=TRUE,add = TRUE     )
    
    mtext("Ground Truth",3)
    mtext("z (um)",side=2,line=2.5)
    axis(1,round(seq(0,max(times2Plot)-10,,5)),las=2)
    
  }
    
    }

L1errors<-matrix(L1errors,ncol=Dlength)
matplot(L1errors,t="l")
#dev.off()
}


dev.off()
















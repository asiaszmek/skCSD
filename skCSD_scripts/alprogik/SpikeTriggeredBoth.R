library(data.table)
library(akima)
source("alprogik/Colors_BlueRed.R")
source("/home/csdori/ksCSD_2014/trunk/alprogik/Colors_BlueRed.R")
setwd("/home/csdori/ksCSD_2014/trunk/simulation/GanglionOsc_128Regular")
Vmem<-as.matrix(read.table("somav.txt"))
Time<-as.matrix(read.table("time"))
elcoord<-matrix(as.matrix(read.table(paste0("/home/csdori/ksCSD_2014/trunk/simulation/Ganglion_d17.5_128Regular/elcoord_x_y_z"))),ncol=3)
#Detecting local maximums

#which maximums are bigger than a treshold
LocMax<-diff(sign(diff(Vmem)))==-2
Sp1<-which(Vmem[ LocMax==TRUE]>-40,arr.ind=TRUE)
Sp2<-which(LocMax==TRUE,arr.ind=TRUE)
Peaks<-Sp2[Sp1,1]

plot(Time, Vmem, t="l")
abline(v=Time[Peaks],col="PURPLE")

PeaksBef<-Peaks[which(Time[Peaks]<390)]
PeaksAft<-Peaks[which(Time[Peaks]>390)]

BaseFile<-"/home/csdori/ksCSD_2014/trunk/simulation/GanglionOsc_256Regular/"
TargetFile<-"/home/csdori/ksCSD_2014/trunk/simulation/Gang128Osc_SpikeTRriggered/"
TargetFileBef<-"/home/csdori/ksCSD_2014/trunk/simulation/Gang128Osc_SpikeTRriggeredBefore/"
TargetFileAft<-"/home/csdori/ksCSD_2014/trunk/simulation/Gang128Osc_SpikeTRriggeredAfter/"
membcurr<-as.matrix(fread(paste0(BaseFile,"membcurr")))
examplelfp<-as.matrix(fread(paste0(BaseFile,"myLFP")))
segcoord<-as.matrix(read.table(paste0(BaseFile,"segcoordinates.txt")))





sCSDTCalc<-function(Elcord,CellCord,const){
  
  d1<-dim(Elcord)[1]
d2<-dim(CellCord)[1]
  TransMatr<-array(0, c(d1,d2))
  for(i in 1:d2){
    for(j in 1:d1){
      TransMatr[j,i]<-1/const*1/sqrt(sum((Elcord[j,]-CellCord[i,])^2))
    }
  }
  return(TransMatr)
}



TransMatrix<-sCSDTCalc(elcoord,segcoord,4*pi*0.5)
#sCSD<-solve(TransMatrix)%*%LFP
NewLFP<-TransMatrix%*%membcurr


#########################################################x


WindowHalf<-60
LFPSpiketr<-array(0,c(dim(NewLFP)[1],WindowHalf*2))
membBef<-array(0,c(dim(membcurr)[1],WindowHalf*2))
for (i in 1: length(PeaksBef)){
LFPSpiketr<-LFPSpiketr+NewLFP[,(PeaksBef[i]-WindowHalf):(PeaksBef[i]+WindowHalf-1)]
membBef<-membBef+membcurr[,(PeaksBef[i]-WindowHalf):(PeaksBef[i]+WindowHalf-1)]

}
LFPSpiketr<-LFPSpiketr/length(PeaksBef)
membBef<-membBef/length(PeaksBef)

LFPSpAft<-array(0,c(dim(NewLFP)[1],WindowHalf*2))
membAft<-array(0,c(dim(membcurr)[1],WindowHalf*2))

for (i in 1: length(PeaksAft)){
LFPSpAft<-LFPSpAft+NewLFP[,(PeaksAft[i]-WindowHalf):(PeaksAft[i]+WindowHalf-1)]
membAft<-membAft+membcurr[,(PeaksBef[i]-WindowHalf):(PeaksBef[i]+WindowHalf-1)]

}
LFPSpAft<-LFPSpAft/length(PeaksAft)
membAft<-membAft/length(PeaksAft)



#Save this LFP and run skCSD on it...

write.table(NewLFP,paste0(TargetFile,"myLFPwhole"),row.names = FALSE, col.names = FALSE)

write.table(LFPSpiketr,paste0(TargetFileBef,"myLFP"),row.names = FALSE, col.names = FALSE)
write.table(membBef,paste0(TargetFileBef,"membSpTrig"),row.names = FALSE, col.names = FALSE)

write.table(LFPSpAft,paste0(TargetFileAft,"myLFP"),row.names = FALSE, col.names = FALSE)
write.table(membAft,paste0(TargetFileAft,"membSpTrig"),row.names = FALSE, col.names = FALSE)



file.copy(paste0(BaseFile,"membcurr"), paste0(TargetFile,"membcurr"), overwrite = TRUE)
file.copy(paste0(BaseFile,"somav.txt"), paste0(TargetFile,"somav.txt"), overwrite = TRUE)
file.copy(paste0(BaseFile,"time"), paste0(TargetFile,"time"), overwrite = TRUE)






















#plot LFP
Plotwidth<-10
CellType<-"gangNew"
xA<-1
yA<-2
xAEl<-1 #1
yAEl<-2
PlotTitle<-"Ganglion Cell"
  outname1<-"/kernelOut_4/"
  ToPlot<--c(70)#-c(1:10,70)
DataFolder<-TargetFileBef #Bef #BaseFile#"/home/csdori/ksCSD_2014/trunk/simulation/gang_5x5_100/"
DataFolderAft<-TargetFileAft
segstart<-matrix(as.matrix(read.table(paste0(DataFolder,'coordsstart_x_y_z'), colClasses='numeric')),ncol=3)
segend<-matrix(as.matrix(read.table(paste0(DataFolder,'coordsend_x_y_z'), colClasses='numeric')),ncol=3)
segmid<-matrix(as.matrix(read.table(paste0(DataFolder,'coordsmid_x_y_z'), colClasses='numeric')),ncol=3)
segwidth1<-as.matrix(read.table(paste0(DataFolder,'segdiam_x_y_z'), colClasses='numeric'))
segwidth<-segwidth1/max(segwidth1)*30 # not realistic?

#membcurr<-as.matrix(read.table(paste0(DataFolder,'membSpTrig')))

 source(paste0(ScriptLocation,"alprogik/KernSmoothDistance.R"))

 membcurr<-KernSmoothDistanceSpTrig(30,DataFolder, DataFolder, outname1,"30")[ToPlot,]
membcurrAft<-KernSmoothDistanceSpTrig(30,DataFolderAft, DataFolderAft, outname1,"30")[ToPlot,]
LFP<-as.matrix(read.table(paste0(DataFolder,'myLFP')))


elec<-elcoord
Xrange<-range(segstart[,xA],elec[,xA])*1.05
Yrange<-range(segstart[,yA],elec[,yA])*1.05
XrangeEl<-range(segstart[,xAEl],elec[,xAEl])*1.05
YrangeEl<-range(segstart[,yAEl],elec[,yAEl])*1.05



 
  Error2Plot<-as.matrix(read.table(paste0(DataFolder,"ErrorL1Smoothed_Noise_SNR0")))#/kernelOut_4/ErrorCV_NoiseElec_SNR0"))) #ErrorL1Smoothed_Noise_SNR0")))
#32-es símításos
CurrMin<- min(Error2Plot)
MinIndex<-which(Error2Plot==CurrMin,arr.ind=TRUE)
  Rread<-R.all[MinIndex[1]]
Lread<-lambda.all[MinIndex[2]]

  Error2PlotAft<-as.matrix(read.table(paste0(DataFolderAft,"ErrorL1Smoothed_Noise_SNR0")))#/kernelOut_4/ErrorCV_NoiseElec_SNR0"))) #ErrorL1Smoothed_Noise_SNR0")))
#32-es símításos
CurrMinAft<- min(Error2PlotAft)
MinIndexAft<-which(Error2PlotAft==CurrMinAft,arr.ind=TRUE)
  RreadAft<-R.all[MinIndexAft[1]]
LreadAft<-lambda.all[MinIndexAft[2]]


#Rread<-32
C.calc<-as.matrix(read.table(paste0(DataFolder,outname1,"/ksCSD_Matr512_R",Rread,"lambda",Lread,"SNR0")))
C.calcAft<-as.matrix(read.table(paste0(DataFolderAft,outname1,"/ksCSD_Matr512_R",RreadAft,"lambda",LreadAft,"SNR0")))

seg.nb<-608

SameplaceAll<-readLines(paste0(BaseFile,outname1,"/sameplace.txt"))
skCSD.all<-array(0,c(seg.nb, dim(LFP)[2]))
skCSD.allAft<-array(0,c(seg.nb, dim(LFP)[2]))
    for(i in 1: seg.nb){
      sameplace<-numeric() 
      sameplace<-as.integer(unlist(strsplit(SameplaceAll," ")[i]))
      #Plot the current related ti the different segments
    #  matplot(cbind(t(C.calc[c(sameplace),])),t='l', xlab="Time steps", ylab="ksCSD", main="Values of current at the same location")
      if(length(sameplace)==1) {skCSD.all[i,]<-C.calc[sameplace,]
      skCSD.allAft[i,]<-C.calcAft[sameplace,]}
      if(length(sameplace)>1) {skCSD.all[i,]<-as.matrix(colSums(C.calc[c(sameplace),]),nrow=1)
      skCSD.allAft[i,]<-as.matrix(colSums(C.calcAft[c(sameplace),]),nrow=1)
      }
    }
    
   
skCSD<-skCSD.all[ToPlot,]
skCSDAft<-skCSD.allAft[ToPlot,]
 dimseg<-dim(skCSD)[1]
dimseg1<-608






dir.create(paste0(DataFolder, "/PlotsBoth"))
for(t2 in 1:dim(LFP)[2]){
TimeInstant<-t2
png(paste0(DataFolder,"/PlotsBoth/SpTrig_",t2,".png"))
par(mfrow=c(2,2))



plot(1,type='n',xlim=Xrange,ylim=Yrange,xlab='x (um)',ylab='y (um)',main="Ground Truth",sub='with synaptic excitations',asp=1)


rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray")
# col1<-ColoursDori(membcurr[,TimeInstant])[[1]]
# ExtVal<-ColoursDori(membcurr[,TimeInstant])[[2]]
col1<-ColoursDori(c(membcurr,membcurrAft))[[1]]
ExtVal<-ColoursDori(c(membcurr,membcurrAft))[[2]]

#highlight the interesting part of the plot
symbols(segmid[ToPlot,xA],segmid[ToPlot,yA], circles=c(diff(range(segstart[,xA]))/Plotwidth* abs(membcurr[,TimeInstant])/ExtVal), add=TRUE, bg=colors()[182], fg=colors()[182], inches=FALSE)
points(elec[,xA],elec[,yA], pch=8) #diff(Xrange)

#col.limits<-max(abs(range(skCSD, membrane.curr.line)))#,memb.currents.vonal[,min.indexes])))

szinskala2<-color.scale(c(membcurr[c(dimseg:1),TimeInstant]),col=col1,zlim=c(-ExtVal*1.01,ExtVal*1.01))
segments(segstart[c(dimseg:1)[ToPlot],xA], segstart[c(dimseg:1)[ToPlot],yA], segend[c(dimseg:1)[ToPlot],xA], segend[c(dimseg:1)[ToPlot],yA], col= szinskala2,lwd=segwidth[c(dimseg:1)[ToPlot]],lend=1)
#segments(segstart[1:dimseg,xA], segstart[1:dimseg,yA], segend[1:dimseg,xA], segend[1:dimseg,yA], col= szinskala2,lwd=segwidth[1:dimseg])
#par(mar=c(0.5,0.1,1,0.5)) 
#mtext("mA/um ",side=4,line=1,cex=0.7)

ColMatr<-as.matrix(c(1:256))
col2<-ColoursDori(ColMatr)[[1]]
# add.image(max(Xrange),min(Yrange)-25, t(ColMatr), adj.x = 0.0, adj.y = 0.0,  image.width = 1, image.height = 0.05, col = col2)
# #axis(4,at=c(min(Yrange),max(Yrange)),labels=c("sink","source"))
# axis(4,at=c(min(Yrange),max(Yrange)),labels=c(round(-ExtVal*1.01,4),round(ExtVal*1.01,4)))
add.image(max(XrangeEl-10),mean(YrangeEl), t(ColMatr), adj.x = 0.0, adj.y = 0.5,  image.width = 1, image.height = 0.05, col = col2)
#axis(4,at=c(min(Yrange),max(Yrange)),labels=c("sink","source"))
axis(4,at=c(min(YrangeEl),max(YrangeEl)),labels=c("sink","source"))




##############################################################################
# After

plot(1,type='n',xlim=Xrange,ylim=Yrange,xlab='x (um)',ylab='y (um)',main="Ground Truth",sub='without synaptic excitations',asp=1)


rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray")
# col1<-ColoursDori(membcurr[,TimeInstant])[[1]]
# ExtVal<-ColoursDori(membcurr[,TimeInstant])[[2]]


#highlight the interesting part of the plot
symbols(segmid[ToPlot,xA],segmid[ToPlot,yA], circles=c(diff(range(segstart[,xA]))/Plotwidth* abs(membcurrAft[,TimeInstant])/ExtVal), add=TRUE, bg=colors()[182], fg=colors()[182], inches=FALSE)
points(elec[,xA],elec[,yA], pch=8) #diff(Xrange)

#col.limits<-max(abs(range(skCSD, membrane.curr.line)))#,memb.currents.vonal[,min.indexes])))

szinskala2<-color.scale(c(membcurrAft[c(dimseg:1),TimeInstant]),col=col1,zlim=c(-ExtVal*1.01,ExtVal*1.01))
segments(segstart[c(dimseg:1)[ToPlot],xA], segstart[c(dimseg:1)[ToPlot],yA], segend[c(dimseg:1)[ToPlot],xA], segend[c(dimseg:1)[ToPlot],yA], col= szinskala2,lwd=segwidth[c(dimseg:1)[ToPlot]],lend=1)
#segments(segstart[1:dimseg,xA], segstart[1:dimseg,yA], segend[1:dimseg,xA], segend[1:dimseg,yA], col= szinskala2,lwd=segwidth[1:dimseg])
#par(mar=c(0.5,0.1,1,0.5)) 
#mtext("mA/um ",side=4,line=1,cex=0.7)

ColMatr<-as.matrix(c(1:256))
col2<-ColoursDori(ColMatr)[[1]]
# add.image(max(Xrange),min(Yrange)-25, t(ColMatr), adj.x = 0.0, adj.y = 0.0,  image.width = 1, image.height = 0.05, col = col2)
# #axis(4,at=c(min(Yrange),max(Yrange)),labels=c("sink","source"))
# axis(4,at=c(min(Yrange),max(Yrange)),labels=c(round(-ExtVal*1.01,4),round(ExtVal*1.01,4)))
add.image(max(XrangeEl-10),mean(YrangeEl), t(ColMatr), adj.x = 0.0, adj.y = 0.5,  image.width = 1, image.height = 0.05, col = col2)
#axis(4,at=c(min(Yrange),max(Yrange)),labels=c("sink","source"))
axis(4,at=c(min(YrangeEl),max(YrangeEl)),labels=c("sink","source"))



#Reconstruction
col1<-ColoursDori(c(skCSD))[[1]]
ExtVal<-ColoursDori(c(skCSD))[[2]]


#TimeInstant<-TimeInstant+1
plot(1,type='n',xlim=Xrange,ylim=Yrange,xlab='x (um)',ylab='y (um)',main="skCSD",sub='with synaptic excitations',asp=1)
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray")
# col1<-ColoursDori(skCSD[,TimeInstant])[[1]]
# ExtVal<-ColoursDori(skCSD[,TimeInstant])[[2]]
#col1<-ColoursDori(skCSD)[[1]]
#ExtVal<-ColoursDori(membcurr)[[2]] #ColoursDori(skCSD)[[2]]
symbols(segmid[ToPlot,xA],segmid[ToPlot,yA], circles=c(diff(range(segstart[,xA]))/Plotwidth* abs(skCSD[,TimeInstant])/ExtVal), add=TRUE, bg=colors()[182], fg=colors()[182], inches=FALSE)


#col.limits<-max(abs(range(skCSD, membrane.curr.line)))#,memb.currents.vonal[,min.indexes])))

szinskala2<-color.scale(c(skCSD[c(dimseg:1),TimeInstant]),col=col1,zlim=c(-ExtVal*1.01,ExtVal*1.01))
points(elec[,xA],elec[,yA], pch=8)
segments(segstart[c(dimseg:1)[ToPlot],xA], segstart[c(dimseg:1)[ToPlot],yA], segend[c(dimseg:1)[ToPlot],xA], segend[c(dimseg:1)[ToPlot],yA], col= szinskala2,lwd=segwidth[c(dimseg:1)[ToPlot]],lend=1)
#par(mar=c(0.5,0.1,1,0.5)) 
#mtext("mA/um ",side=4,line=1,cex=0.7)

ColMatr<-as.matrix(c(1:256))
col2<-ColoursDori(ColMatr)[[1]]
#add.image(max(Xrange),min(Yrange)-25, t(ColMatr), adj.x = 0.0, adj.y = 0.0,  image.width = 1, image.height = 0.05, col = col2)
# add.image(t(ColMatr),adj.x = 0.0, adj.y = 0.0,  image.width = 1, image.height = 0.05, col = col2)
# #axis(4,at=c(min(Yrange),max(Yrange)),labels=c("sink","source"))
# # axis(4,at=c(min(Yrange),max(Yrange)),labels=c(round(-ExtVal*1.01,4),round(ExtVal*1.01,4)))
# axis(4,labels=c(round(-ExtVal*1.01,4),round(ExtVal*1.01,4)))
add.image(max(XrangeEl-10),mean(YrangeEl), t(ColMatr), adj.x = 0.0, adj.y = 0.5,  image.width = 1, image.height = 0.05, col = col2)
#axis(4,at=c(min(Yrange),max(Yrange)),labels=c("sink","source"))
axis(4,at=c(min(YrangeEl),max(YrangeEl)),labels=c("sink","source"))

##############################################################################
col1<-ColoursDori(c(skCSDAft))[[1]]
ExtVal<-ColoursDori(c(skCSDAft))[[2]]

#TimeInstant<-TimeInstant+1
plot(1,type='n',xlim=Xrange,ylim=Yrange,xlab='x (um)',ylab='y (um)',main='skCSD',sub='with synaptic excitations',asp=1)
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray")
# col1<-ColoursDori(skCSD[,TimeInstant])[[1]]
# ExtVal<-ColoursDori(skCSD[,TimeInstant])[[2]]
#col1<-ColoursDori(skCSD)[[1]]
#ExtVal<-ColoursDori(membcurr)[[2]] #ColoursDori(skCSD)[[2]]
symbols(segmid[ToPlot,xA],segmid[ToPlot,yA], circles=c(diff(range(segstart[,xA]))/Plotwidth* abs(skCSDAft[,TimeInstant])/ExtVal), add=TRUE, bg=colors()[182], fg=colors()[182], inches=FALSE)


#col.limits<-max(abs(range(skCSD, membrane.curr.line)))#,memb.currents.vonal[,min.indexes])))

szinskala2<-color.scale(c(skCSDAft[c(dimseg:1),TimeInstant]),col=col1,zlim=c(-ExtVal*1.01,ExtVal*1.01))
points(elec[,xA],elec[,yA], pch=8)
segments(segstart[c(dimseg:1)[ToPlot],xA], segstart[c(dimseg:1)[ToPlot],yA], segend[c(dimseg:1)[ToPlot],xA], segend[c(dimseg:1)[ToPlot],yA], col= szinskala2,lwd=segwidth[c(dimseg:1)[ToPlot]],lend=1)
#par(mar=c(0.5,0.1,1,0.5)) 
#mtext("mA/um ",side=4,line=1,cex=0.7)

ColMatr<-as.matrix(c(1:256))
col2<-ColoursDori(ColMatr)[[1]]
#add.image(max(Xrange),min(Yrange)-25, t(ColMatr), adj.x = 0.0, adj.y = 0.0,  image.width = 1, image.height = 0.05, col = col2)
# add.image(t(ColMatr),adj.x = 0.0, adj.y = 0.0,  image.width = 1, image.height = 0.05, col = col2)
# #axis(4,at=c(min(Yrange),max(Yrange)),labels=c("sink","source"))
# # axis(4,at=c(min(Yrange),max(Yrange)),labels=c(round(-ExtVal*1.01,4),round(ExtVal*1.01,4)))
# axis(4,labels=c(round(-ExtVal*1.01,4),round(ExtVal*1.01,4)))
add.image(max(XrangeEl-10),mean(YrangeEl), t(ColMatr), adj.x = 0.0, adj.y = 0.5,  image.width = 1, image.height = 0.05, col = col2)
#axis(4,at=c(min(Yrange),max(Yrange)),labels=c("sink","source"))
axis(4,at=c(min(YrangeEl),max(YrangeEl)),labels=c("sink","source"))






#points(elec[,xA],elec[,yA])
#if(file.exists(elignorename))  points(elec[ El2Ignore$x,xA],elec[ El2Ignore$x,yA],pch=4,col="RED")


#plot(LFP[91,],t="l",ylab="LFP close to soma")
      #abline(v=TimeInstant,col="RED", lwd=2)

par(oma=c(0,0,2,0))
#title(paste("Time:", round(SimTime[TimeInstant],3), "ms") , outer=TRUE)

#title(paste(t2) , outer=TRUE)
dev.off()
}

#Plot to compare spike triggered averages
#with synaptic input

# only oscillatory

# constant current input




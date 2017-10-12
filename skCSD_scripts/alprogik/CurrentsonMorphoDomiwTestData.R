#source("/media/zoe/BA0ED4600ED416EB/agy/ksCSD_2014/trunk/alprogik/CurrentsonMorpho2.R")
#source("/media/BA0ED4600ED416EB/agy/ksCSD_2014/trunk/alprogik/Colors_BlueRed.R")
source("alprogik/Colors_BlueRed.R")
library(data.table)
library(MASS)
#image scale function for plotting color bars for image plots

DataFolder<-"/home/csdori/ksCSD_2014/trunk/simulation/cell_Cos3/"#cell_Cos2/"#Domi_20ms/"#cell_D/" #gang_5x5_100/"
KernelLocs<-"/home/csdori/ksCSD_2014/trunk/simulation/Domi14/kernelOut_4/"
sig<-0.3

#DataFolder<-"/home/csdori/ksCSD_2014/trunk/simulation/Domi100ms/"#cell_D/" #gang_5x5_100/"


#DataFolder<-#"/home/csdori/ksCSD_2014/trunk/simulation/Ganglion_d17.5_128Regular/"#BS_d50_el128/"
#DataFolder<-"/home/csdori/ksCSD_2014/trunk/simulation/Y_el4x16_RotS_symm_d50_ver0/"
#DataFolder<-"/media/zoe/BA0ED4600ED416EB/agy/ksCSD_2014/trunk/simulation/Y_el4x16_Rot_symm_d50_ver0/"
# DataFolder<-"/home/csdori/ksCSD_2014/trunk/simulation/GanglionOsc_128Regular/" #gang_5x5_100/"
 #ffmpeg  -framerate 10 -pattern_type glob -i 'CSD_Morpho*.png' CSD.mp4
#DataFolder<-"/media/zoe/BA0ED4600ED416EB/agy/ksCSD_2014/trunk/simulation/Hex_8x8_inter70_d30/"

CellType<-"Domi" #"gangNew"#"Y" # "gangNew" #"Y-rotated" #


if(CellType=="Domi"){
  xA<-2
  yA<-1
  xAEl<-2 #1
  yAEl<-1
  PlotTitle<-"Domi"
  outname1<-"/kernelOut_4/"
  ToPlot<-c(1:895)
}

inname<-DataFolder



PlotMorpho<-function(DataFolder,Plotwidth){
El2Ignore<-numeric()
segstart<--matrix(as.matrix(read.table(paste0(DataFolder,'coordsstart_x_y_z'), colClasses='numeric')),ncol=3)
segend<--matrix(as.matrix(read.table(paste0(DataFolder,'coordsend_x_y_z'), colClasses='numeric')),ncol=3)
segmid<--matrix(as.matrix(read.table(paste0(DataFolder,'coordsmid_x_y_z'), colClasses='numeric')),ncol=3)
coordsEnd<--as.matrix(read.table(paste0(DataFolder,"/coordsend_x_y_z")))
elec<--matrix(as.matrix(read.table(paste0(DataFolder,'elcoord_x_y_z'), colClasses='numeric')),ncol=3)
seg.nb<-dim(segmid)[1]

WhichTimePlot<-1:5#100:300
SimTime<-as.matrix(read.table(paste0(DataFolder,'/time')))
LFP<-(as.matrix(fread(paste0(DataFolder,'myLFP'))))
membcurr<-(as.matrix(fread(paste0(DataFolder,'membcurr'))))
LFP<-LFP[,WhichTimePlot]
SimTime<-SimTime[WhichTimePlot]
#Smooth LFP on 2D
library(akima)
#LFPinterpol<-interp(elec[,xAEl],elec[,yAEl],LFP[,2])
#image(interp(elec[,xAEl],elec[,yAEl],LFP[,99],seq(-200,200,,40),seq(-200,600,,40)))
#image(interp(elec[,xAEl],elec[,yAEl],LFP[,1],seq(-500,500,,40),seq(-400,600,,40)))




TimeMax<-length(SimTime)
#TimeInstant<-which(SimTime==45) #45
#TimeInstant<-TimeInstant+2
########!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#skCSD<-as.matrix(read.table(paste0(DataFolder,outname1,"/skCSD_plain")))#paste0(outname,outname1,"/skCSDall_M",256,"_R",55,"_SNR0_lambda0")))

if(CellType=="Y-rotated" | CellType=="Y" | CellType=="gangNew" | CellType=="Domi"){
  SNR<-0
  Timeframe=TRUE
  if (Timeframe==TRUE){
 Error2Plot<- as.matrix(read.table(paste0(DataFolder,outname1,"/ErrorL1Smoothed_Noise_SNR0")))#_1_1000")))#as.matrix(read.table(paste0(DataFolder,"/ErrorTable")))
#32-es símításos
CurrMin<- min(Error2Plot)
MinIndex<-which(Error2Plot==CurrMin,arr.ind=TRUE)
  
  } else {
  Error2Plot<- as.matrix(read.table(paste0(DataFolder,"/ErrorTable")))
#32-es símításos
CurrMin<- min(Error2Plot[(2*length(R.all)+1):(3*length(R.all)),])
MinIndex<-which(Error2Plot[(2*length(R.all)+1):(3*length(R.all)),]==CurrMin,arr.ind=TRUE)
  
  }
}
#MinIndex<-data.matrix(read.table(paste(DataFolder,outname1,"/MinIndex_Noise","_SNR",SNR,sep="")))
M<-512
lambda.all<-c("1e-05","1e-04","0.001","0.01","0.1" )
R.all<-2^(3:7)
Rread<-R.all[MinIndex[1]] #32
#MinIndex[2]<-2
Lread<-lambda.all[MinIndex[2]]
R<-Rread
lamb<-MinIndex[2]
plotFileName<-paste0("ComparingL1TIme",R,"_",Lread)
dir.create(paste0(DataFolder,plotFileName))



Kmatr<-as.matrix(read.table(paste(KernelLocs,"/K_M",M,"_R",R,"lambda",lambda.all[lamb],sep="")))
Ktildematr<-as.matrix(read.table(paste(KernelLocs,"/Ktilda_M",M,"_R",R,"lambda",lambda.all[lamb],sep="")))

Tmatr<-(4*pi*sig)*t(Ktildematr)%*%ginv(Kmatr)



skCSD.all.part<-Tmatr%*%LFP
write.table(skCSD.all.part, paste(outname,"/ksCSD_Matr",M,"_R",R,"lambda",
                                  lambda.all[lamb], "SNR",SNR,sep=""),row.names = FALSE,col.names = FALSE)



skCSD.all<-array(0,c(seg.nb, dim(LFP)[2]))

SameplaceAll<-readLines(paste0(KernelLocs,"/sameplace.txt"))
#seg.nb<-max(as.matrix(read.table(paste(inname,outname1,'/connection_touched_once',sep=''))))

for(i in 1: seg.nb){
  sameplace<-numeric() 
  sameplace<-as.integer(unlist(strsplit(SameplaceAll," ")[i]))
  #Plot the current related ti the different segments
  #  matplot(cbind(t(C.calc[c(sameplace),])),t='l', xlab="Time steps", ylab="ksCSD", main="Values of current at the same location")
  if(length(sameplace)==1) skCSD.all[i,]<- skCSD.all.part[sameplace,]
  if(length(sameplace)>1) skCSD.all[i,]<-as.matrix(colSums( skCSD.all.part[c(sameplace),]),nrow=1)
}
skCSD<-skCSD.all
##############!!!!!!!!!!!!!!!!!!!!!




seglength<-as.matrix(read.table(paste0(DataFolder,'/seglength')))
funaramvonal<-function(x) x/seglength

membcurrLine<-apply(membcurr,2,funaramvonal)


#Plotratio<-Xrange/max(segwidth)
segwidth1<-as.matrix(read.table(paste0(DataFolder,'segdiam_x_y_z'), colClasses='numeric'))
segwidth<-segwidth1/max(segwidth1)*30 # not realistic?
#membcurr<-membcurrO

for(Tcounter in seq(1,length(SimTime),1)){
  TimeInstant<-Tcounter
  
plotname<-paste(paste(DataFolder,plotFileName,"/CSDTest_MorphoCSD", "-", 
            format(formatC(Tcounter, digits = 4, flag = "0"),scientific=FALSE), sep = ""), "png", sep = ".")


if(CellType=="Domi"){
  png(plotname, height=1000, width=900, pointsize=20)
  #TimeInstant<-500
 # par(mfrow=c(1,2))
  
  
}


dimseg<-dim(segstart)[1]
Xrange<-range(segstart[,xA],elec[,xA])*1.05
Yrange<-range(segstart[,yA],elec[,yA])*1.05
XrangeEl<-range(segstart[,xAEl],elec[,xAEl])*1.05
YrangeEl<-range(segstart[,yAEl],elec[,yAEl])*1.05
#LFP 
#TimeInstant<-200
#col1<-ColoursDori(LFP[,TimeInstant])[[1]]
#ExtVal<-ColoursDori(LFP[,TimeInstant])[[2]]
#plot.new()
split.screen(rbind(c(0,1,0.71, 1), c(0, 0.49, 0, 0.7),c( 0.51,1, 0, 0.7) ))
screen(1)
par(mar =c(4, 3, 2, 2))#par(oma=c(0,0,0,0))

matplot(SimTime/20-10,t(LFP[14:10,]),t="l", xlab="Time (ms)",lwd=2, ylab="Potential (uV)",lty=1,col=2:6)
abline(v=SimTime[TimeInstant]/20-10, lwd=2)
legend("bottomleft",paste(c(1:5)),col=2:6,ncol=5, pch=20,title="Electrodes")




NonLinPar<-1
what2Plot<-atan(LFP*NonLinPar/max(abs(LFP)*(pi/2-0.1)))

col1<-ColoursDoriLFP(what2Plot)[[1]]
ExtVal<-ColoursDoriLFP(what2Plot)[[2]]
screen(2)

#BreaksIm<-tan(seq(-ExtVal, ExtVal,,201)/ExtVal*(pi/2-0.1)*NonLinPar)/tan(pi/2-0.1)

plot(1,type='n',xlim=XrangeEl,ylim=YrangeEl,xlab='x (um)',ylab='y (um)',main="Test Ground Truth",asp=1)
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray")
# col1<-ColoursDori(skCSD[,TimeInstant])[[1]]
# ExtVal<-ColoursDori(skCSD[,TimeInstant])[[2]]
what2Plot<-atan(membcurrLine*NonLinPar/max(abs(membcurrLine)*(pi/2-0.1)))

col1<-ColoursDori(what2Plot)[[1]]
ExtVal<-ColoursDori(what2Plot)[[2]]




#symbols(segmid[c(1:dimseg)[ToPlot],xA],segmid[c(1:dimseg)[ToPlot],yA], circles=c(diff(range(segstart[,xA]))/Plotwidth* abs(skCSD[,TimeInstant])/ExtVal), add=TRUE, bg=colors()[182], fg=colors()[182], inches=FALSE)


#col.limits<-max(abs(range(skCSD, membrane.curr.line)))#,memb.currents.vonal[,min.indexes])))

szinskala2<-color.scale(c(what2Plot[,TimeInstant]),col=col1,zlim=c(-ExtVal*1.01,ExtVal*1.01))
points(elec[,xA],elec[,yA], pch=8)
segments(segstart[c(1:dimseg)[ToPlot],xA], segstart[c(1:dimseg)[ToPlot],yA], segend[c(1:dimseg)[ToPlot],xA], segend[c(1:dimseg)[ToPlot],yA], col= szinskala2,lwd=segwidth[c(1:dimseg)[ToPlot]],lend=1)
#par(mar=c(0.5,0.1,1,0.5)) 
#mtext("mA/um ",side=4,line=1,cex=0.7)
ColMatr<-as.matrix(c(-100:100))

what2Plot<-atan(ColMatr*NonLinPar/max(abs(ColMatr)*(pi/2-0.1)))

col2<-ColoursDori(what2Plot)[[1]]

#add.image(max(Xrange),min(Yrange)-25, t(ColMatr), adj.x = 0.0, adj.y = 0.0,  image.width = 1, image.height = 0.05, col = col2)
# add.image(t(ColMatr),adj.x = 0.0, adj.y = 0.0,  image.width = 1, image.height = 0.05, col = col2)
# #axis(4,at=c(min(Yrange),max(Yrange)),labels=c("sink","source"))
# # axis(4,at=c(min(Yrange),max(Yrange)),labels=c(round(-ExtVal*1.01,4),round(ExtVal*1.01,4)))
# axis(4,labels=c(round(-ExtVal*1.01,4),round(ExtVal*1.01,4)))
add.image(max(XrangeEl-10),mean(YrangeEl), t(what2Plot), adj.x = 0.0, adj.y = 0.5,  image.width = 1, image.height = 0.05, col = col2)
#axis(4,at=c(min(Yrange),max(Yrange)),labels=c("sink","source"))
axis(4,at=c(min(YrangeEl),max(YrangeEl)),labels=c("sink","source"))

#Reconstruction
#################ksCSD

#Ground Truth
#TimeInstant<-TimeInstant+1
#######################x

#Reconstruction

screen(3)
#TimeInstant<-TimeInstant+1
plot(1,type='n',xlim=XrangeEl,ylim=YrangeEl,xlab='x (um)',ylab='y (um)',main="skCSD",asp=1)
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray")
# col1<-ColoursDori(skCSD[,TimeInstant])[[1]]
# ExtVal<-ColoursDori(skCSD[,TimeInstant])[[2]]
what2Plot<-atan(skCSD*NonLinPar/max(abs(skCSD)*(pi/2-0.1)))

col1<-ColoursDori(what2Plot)[[1]]
ExtVal<-ColoursDori(what2Plot)[[2]]




#symbols(segmid[c(1:dimseg)[ToPlot],xA],segmid[c(1:dimseg)[ToPlot],yA], circles=c(diff(range(segstart[,xA]))/Plotwidth* abs(skCSD[,TimeInstant])/ExtVal), add=TRUE, bg=colors()[182], fg=colors()[182], inches=FALSE)


#col.limits<-max(abs(range(skCSD, membrane.curr.line)))#,memb.currents.vonal[,min.indexes])))

szinskala2<-color.scale(c(what2Plot[,TimeInstant]),col=col1,zlim=c(-ExtVal*1.01,ExtVal*1.01))
points(elec[,xA],elec[,yA], pch=8)
segments(segstart[c(1:dimseg)[ToPlot],xA], segstart[c(1:dimseg)[ToPlot],yA], segend[c(1:dimseg)[ToPlot],xA], segend[c(1:dimseg)[ToPlot],yA], col= szinskala2,lwd=segwidth[c(1:dimseg)[ToPlot]],lend=1)
#par(mar=c(0.5,0.1,1,0.5)) 
#mtext("mA/um ",side=4,line=1,cex=0.7)
ColMatr<-as.matrix(c(-100:100))

what2Plot<-atan(ColMatr*NonLinPar/max(abs(ColMatr)*(pi/2-0.1)))

col2<-ColoursDori(what2Plot)[[1]]

#add.image(max(Xrange),min(Yrange)-25, t(ColMatr), adj.x = 0.0, adj.y = 0.0,  image.width = 1, image.height = 0.05, col = col2)
# add.image(t(ColMatr),adj.x = 0.0, adj.y = 0.0,  image.width = 1, image.height = 0.05, col = col2)
# #axis(4,at=c(min(Yrange),max(Yrange)),labels=c("sink","source"))
# # axis(4,at=c(min(Yrange),max(Yrange)),labels=c(round(-ExtVal*1.01,4),round(ExtVal*1.01,4)))
# axis(4,labels=c(round(-ExtVal*1.01,4),round(ExtVal*1.01,4)))
add.image(max(XrangeEl-10),mean(YrangeEl), t(what2Plot), adj.x = 0.0, adj.y = 0.5,  image.width = 1, image.height = 0.05, col = col2)
#axis(4,at=c(min(Yrange),max(Yrange)),labels=c("sink","source"))
axis(4,at=c(min(YrangeEl),max(YrangeEl)),labels=c("sink","source"))

#points(elec[,xA],elec[,yA])
#if(file.exists(elignorename))  points(elec[ El2Ignore$x,xA],elec[ El2Ignore$x,yA],pch=4,col="RED")
par(oma=c(0,0,2,0))
title(paste("Time:", format(round(SimTime[TimeInstant]/20-10,3),nsmall=2), "ms") , outer=TRUE)
dev.off()
#close.screen(all = TRUE)  
}
}
PlotMorpho(DataFolder,50)

#system("convert -delay 40 *.png Compare.gif") #taurin muxik

# cleaning up
#file.remove(list.files(pattern=".png"))

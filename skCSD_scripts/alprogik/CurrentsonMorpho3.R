#source("/media/zoe/BA0ED4600ED416EB/agy/ksCSD_2014/trunk/alprogik/CurrentsonMorpho2.R")
#source("/media/BA0ED4600ED416EB/agy/ksCSD_2014/trunk/alprogik/Colors_BlueRed.R")
source("alprogik/Colors_BlueRed.R")
library(data.table)
#image scale function for plotting color bars for image plots


#DataFolder<-#"/home/csdori/ksCSD_2014/trunk/simulation/Ganglion_d17.5_128Regular/"#BS_d50_el128/"
#DataFolder<-"/home/csdori/ksCSD_2014/trunk/simulation/Y_el4x16_RotS_symm_d50_ver0/"
#DataFolder<-"/media/zoe/BA0ED4600ED416EB/agy/ksCSD_2014/trunk/simulation/Y_el4x16_Rot_symm_d50_ver0/"
 DataFolder<-"/home/csdori/ksCSD_2014/trunk/simulation/GanglionOsc_128Regular/"
 #ffmpeg  -framerate 10 -pattern_type glob -i 'CSD_Morpho*.png' CSD.mp4
#DataFolder<-"/media/zoe/BA0ED4600ED416EB/agy/ksCSD_2014/trunk/simulation/Hex_8x8_inter70_d30/"
plotFileName<-"ComparingSmoothedL1"
dir.create(paste0(DataFolder,plotFileName))
CellType<-"gangNew"#"Y" # "gangNew" #"Y-rotated" #

06202831913
if(CellType=="gangNew"){
xA<-1
yA<-2
xAEl<-1 #1
yAEl<-2
PlotTitle<-"Ganglion Cell"
  outname1<-"/kernelOut_4/"
  ToPlot<--c(1:10,70)
TanhPar<-0.9
}

if(CellType=="BS"){
  xA<-1
yA<-3
xAEl<-2 #1
yAEl<-3
ToPlot<-1:52
PlotTitle<-"BS"
  outname1<-"/kernelOut_poli_CP3/"
}

if(CellType=="Y"){
xA<-1
yA<-3
xAEl<-1
yAEl<-3
ToPlot<-1:86
PlotTitle<-"Y-shaped"
  outname1<-"/kernelOut_4/"
}


if(CellType=="Y-rotated"){
xA<-1
yA<-3
xAEl<-2 #1
yAEl<-3
ToPlot<-1:86
PlotTitle<-"Y-shaped"
  outname1<-"/kernelOut_4/"
}
PlotMorpho<-function(DataFolder,Plotwidth){
El2Ignore<-numeric()
segstart<-matrix(as.matrix(read.table(paste0(DataFolder,'coordsstart_x_y_z'), colClasses='numeric')),ncol=3)
segend<-matrix(as.matrix(read.table(paste0(DataFolder,'coordsend_x_y_z'), colClasses='numeric')),ncol=3)
segmid<-matrix(as.matrix(read.table(paste0(DataFolder,'coordsmid_x_y_z'), colClasses='numeric')),ncol=3)
elec<-matrix(as.matrix(read.table(paste0(DataFolder,'elcoord_x_y_z'), colClasses='numeric')),ncol=3)
seg.nb<-dim(segmid)[1]
LFP<-as.matrix(fread(paste0(DataFolder,'myLFP')))
#Smooth LFP on 2D
library(akima)
LFPinterpol<-interp(elec[,xAEl],elec[,yAEl],LFP[,2])
#image(interp(elec[,xAEl],elec[,yAEl],LFP[,99],seq(-200,200,,40),seq(-200,600,,40)))
image(interp(elec[,xAEl],elec[,yAEl],LFP[,99],seq(-500,500,,40),seq(-400,600,,40)))
if(file.exists(paste0(DataFolder,"kCSDCurr2"))==TRUE){
KCSDcurr<- as.matrix(fread(paste0(DataFolder,"kCSDCurr2")))
xloc<- as.matrix(read.table(paste0(DataFolder,"csdxLoc")))
yloc<- as.matrix(read.table(paste0(DataFolder,"csdyLoc")))
xloc1<-unique(c(xloc))
yloc1<-unique(c(yloc))
}
KcsdFileName<-'Dori_kCSDBest.h5'
if(file.exists(paste0(DataFolder,KcsdFileName))==TRUE){
  
  library("rhdf5")
  EstCSD<-h5read(paste0(DataFolder,KcsdFileName), '/est_csd')
  EstPot<-h5read(paste0(DataFolder,KcsdFileName), '/est_pot')
  SpaceX<-c(h5read(paste0(DataFolder,KcsdFileName), '/space_X'))
  SpaceY<-c(h5read(paste0(DataFolder,KcsdFileName), '/space_Y'))
  KCSDPar<-h5read(paste0(DataFolder,KcsdFileName), '/kcsd_result')
  H5close()
  #Zlim<-range(EstCSD)
}
SimTime<-as.matrix(read.table(paste0(DataFolder,'/time')))
TimeMax<-length(SimTime)
TimeInstant<-which(SimTime==45) #45
TimeInstant<-TimeInstant+2
########!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#skCSD<-as.matrix(read.table(paste0(DataFolder,outname1,"/skCSD_plain")))#paste0(outname,outname1,"/skCSDall_M",256,"_R",55,"_SNR0_lambda0")))
if(CellType=="Y-rotated" | CellType=="Y" | CellType=="gangNew"){
  SNR<-0
  
  Error2Plot<-as.matrix(read.table(paste0(DataFolder,"/ErrorTable")))
#32-es símításos
CurrMin<- min(Error2Plot)
MinIndex<-which(Error2Plot==CurrMin,arr.ind=TRUE)
  
  
  
#MinIndex<-data.matrix(read.table(paste(DataFolder,outname1,"/MinIndex_Noise","_SNR",SNR,sep="")))
M<-512
lambda.all<-c("1e-05","1e-04","0.001","0.01","0.1" )
R.all<-2^(3:7)
Rread<-R.all[MinIndex[1]]
Lread<-lambda.all[MinIndex[2]]
C.calc<-as.matrix(read.table(paste0(DataFolder,outname1,"/skCSDall_M512_R",Rread,"lambda",Lread)))
} else {
  C.calc<-as.matrix(read.table(paste0(DataFolder,outname1,"/skCSDall_M512_R30lambda1e-04")))#30lambda1e-04")))#80lambda0.001")))
}
SameplaceAll<-readLines(paste0(DataFolder,outname1,"/sameplace.txt"))
skCSD.all<-array(0,c(seg.nb, dim(LFP)[2]))
    for(i in 1: seg.nb){
      sameplace<-numeric() 
      sameplace<-as.integer(unlist(strsplit(SameplaceAll," ")[i]))
      #Plot the current related ti the different segments
    #  matplot(cbind(t(C.calc[c(sameplace),])),t='l', xlab="Time steps", ylab="ksCSD", main="Values of current at the same location")
      if(length(sameplace)==1) skCSD.all[i,]<-C.calc[sameplace,]
      if(length(sameplace)>1) skCSD.all[i,]<-as.matrix(colSums(C.calc[c(sameplace),]),nrow=1)
    }
skCSD<-skCSD.all[ToPlot,]
##############!!!!!!!!!!!!!!!!!!!!!
coordsEnd<-as.matrix(read.table(paste0(DataFolder,"/coordsend_x_y_z")))
seglength<-as.matrix(read.table(paste0(DataFolder,'/seglength')))
funaramvonal<-function(x) x/seglength
membCurr<-as.matrix(read.table(paste0(DataFolder,'/membcurr')))
membcurr0<-apply(membCurr,2,funaramvonal)
somaPot<-as.matrix(read.table(paste0(DataFolder,'/somav.txt')))
#Plotratio<-Xrange/max(segwidth)
segwidth1<-as.matrix(read.table(paste0(DataFolder,'segdiam_x_y_z'), colClasses='numeric'))
segwidth<-segwidth1/max(segwidth1)*30 # not realistic?
#membcurr<-membcurrO
#membcurrO<-membcurr
membcurr<-membcurr0[ToPlot,] 
for(Tcounter in 1:(TimeMax/2)){
  TimeInstant<-1+Tcounter*2
  
plotname<-paste(paste(DataFolder,plotFileName,"/CSDL1_MorphoCSD", "-", 
            formatC(Tcounter, digits = 3, flag = "0"), sep = ""), "png", sep = ".")
            
            
            
            
            
            
if(CellType=="Y"  | CellType=="Y-rotated" | CellType=="BS" ){            
#png(plotname, height=500, width=1200, pointsize=20)
png(plotname,width=5, height=3,units='in',res=300, pointsize=)
#TimeInstant<-500
par(mfrow=c(1,4))
}
if(CellType=="gangNew"){
png(plotname, height=1300, width=900, pointsize=20)
#TimeInstant<-500
#par(mfrow=c(2,2))
split.screen(rbind(c(0,0.95,0.8, 0.95), c(0, 0.45, 0.4, 0.8),c( 0.55,1, 0.4, 0.8), c(0, 0.45, 0, 0.4),c( 0.55,1, 0, 0.4)))
screen(1)
par(mar = rep(1, 4))#par(oma=c(0,0,0,0))
 plot(SimTime,somaPot,t="l", xlab="Time (ms)",ylab=expression("V"["soma"]))
      abline(v=SimTime[TimeInstant])
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
screen(2)
LFP<-tanh(TanhPar*LFP/max(abs(LFP)))
col1<-ColoursDoriLFP(LFP)[[1]]
ExtVal<-ColoursDoriLFP(LFP)[[2]]
if(CellType!="BS") image(interp(elec[,xAEl],elec[,yAEl],LFP[,TimeInstant],seq(XrangeEl[1],XrangeEl[2],,40),seq(YrangeEl[1],YrangeEl[2],,40)),col=col1, 
      zlim=c(-ExtVal, ExtVal),xlim=XrangeEl,ylim=YrangeEl,main="LFP",xlab='x (um)',ylab='y (um)', asp=1)
if(CellType=="BS") image(as.matrix(LFP[,TimeInstant]),col=col1, 
      zlim=c(-ExtVal, ExtVal),xlim=XrangeEl,ylim=YrangeEl,main="LFP",xlab='x (um)',ylab='y (um)', asp=1)
segments(segstart[c(1:dimseg),xAEl], segstart[c(1:dimseg),yAEl], segend[c(1:dimseg),xAEl], segend[c(1:dimseg),yAEl],lwd=1)
#mtext("mV",side=4,line=1,cex=0.7)
points(elec[,xAEl],elec[,yAEl], pch=8)
ColMatr<-as.matrix(c(1:256))
col2<-ColoursDoriLFP(ColMatr)[[1]]
# add.image(max(XrangeEl-10),min(YrangeEl), t(ColMatr), adj.x = 0.0, adj.y = 0.0,  image.width = 1, image.height = 0.05, col = col2)
# #axis(4,at=c(min(Yrange),max(Yrange)),labels=c("sink","source"))
# axis(4,at=c(min(YrangeEl),max(YrangeEl)),labels=c(round(-ExtVal,4),round(ExtVal,4)))
add.image(max(XrangeEl-10),mean(YrangeEl), t(ColMatr), adj.x = 0.0, adj.y = 0.5,  image.width = 1, image.height = 0.05, col = col2)
#axis(4,at=c(min(Yrange),max(Yrange)),labels=c("sink","source"))
axis(4,at=c(min(YrangeEl),max(YrangeEl)),labels=c("sink","source"))#labels=c(round(-ExtVal,4),round(ExtVal,4)))
if(file.exists(paste0(DataFolder,KcsdFileName))==TRUE){
#Zlim<-range(EstCSD)
screen(3)
EstCSD<-tanh(TanhPar*EstCSD/max(abs(EstCSD)))
col1<-ColoursDori(data.matrix(EstCSD))[[1]]
ExtVal<-ColoursDori(data.matrix(EstCSD))[[2]]
image(unique(SpaceX)*1000, unique(SpaceY)*1000, t(data.matrix(EstCSD[TimeInstant,,])),xlim=XrangeEl,ylim=YrangeEl, col=col1,
,main="kCSD",xlab='x (um)',ylab='y (um)',asp=1,zlim=c(-ExtVal, ExtVal))
#image(unique(c(xloc)),unique(c(yloc)),matrix(KCSDcurr[,TimeInstant],nrow=length(yloc1)),col=col1,    xlim=XrangeEl,ylim=YrangeEl
segments(segstart[c(1:dimseg),xAEl], segstart[c(1:dimseg),yAEl], segend[c(1:dimseg),xAEl], segend[c(1:dimseg),yAEl],lwd=1)
points(elec[,xAEl],elec[,yAEl], pch=8)
#mtext("??",side=4,line=1,cex=0.7)
ColMatr<-as.matrix(c(1:256))
col2<-ColoursDori(ColMatr)[[1]]
# 
# add.image(max(XrangeEl)-10,min(YrangeEl), t(ColMatr), adj.x = 0.0, adj.y = 0.0,  image.width = 1, image.height = 0.05, col = col2)
# #axis(4,at=c(min(Yrange),max(Yrange)),labels=c("sink","source"))
# axis(4,at=c(min(YrangeEl),max(YrangeEl)),labels=c(format(-ExtVal, digits=3),format(ExtVal,digits=3)))
add.image(max(XrangeEl-10),mean(YrangeEl), t(ColMatr), adj.x = 0.0, adj.y = 0.5,  image.width = 1, image.height = 0.05, col = col2)
#axis(4,at=c(min(Yrange),max(Yrange)),labels=c("sink","source"))
axis(4,at=c(min(YrangeEl),max(YrangeEl)),labels=c("sink","source"))#labels=c(round(-ExtVal,4),round(ExtVal,4))
}
#Ground Truth
#TimeInstant<-TimeInstant+1
#######################x
screen(4)
plot(1,type='n',xlim=Xrange,ylim=Yrange,xlab='x (um)',ylab='y (um)',main="Ground Truth",asp=1)
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray")
# col1<-ColoursDori(membcurr[,TimeInstant])[[1]]
# ExtVal<-ColoursDori(membcurr[,TimeInstant])[[2]]
membcurr<-tanh(TanhPar*membcurr/max(abs(membcurr)))
col1<-ColoursDori(membcurr)[[1]]
ExtVal<-ColoursDori(membcurr)[[2]]
#highlight the interesting part of the plot
symbols(segmid[c(1:dimseg)[ToPlot],xA],segmid[c(1:dimseg)[ToPlot],yA], circles=c(diff(range(segstart[,xA]))/Plotwidth* abs(membcurr[,TimeInstant])/ExtVal), add=TRUE, bg=colors()[182], fg=colors()[182], inches=FALSE)
points(elec[,xA],elec[,yA], pch=8) #diff(Xrange)
#col.limits<-max(abs(range(skCSD, membrane.curr.line)))#,memb.currents.vonal[,min.indexes])))
szinskala2<-color.scale(c(membcurr[,TimeInstant]),col=col1,zlim=c(-ExtVal*1.01,ExtVal*1.01))
segments(segstart[c(1:dimseg)[ToPlot],xA], segstart[c(1:dimseg)[ToPlot],yA], segend[c(1:dimseg)[ToPlot],xA], segend[c(1:dimseg)[ToPlot],yA], col= szinskala2,lwd=segwidth[c(1:dimseg)[ToPlot]],lend=1)
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
screen(5)
#TimeInstant<-TimeInstant+1
plot(1,type='n',xlim=Xrange,ylim=Yrange,xlab='x (um)',ylab='y (um)',main="skCSD",asp=1)
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray")
# col1<-ColoursDori(skCSD[,TimeInstant])[[1]]
# ExtVal<-ColoursDori(skCSD[,TimeInstant])[[2]]
skCSD<-tanh(TanhPar*skCSD/max(abs(skCSD)))
col1<-ColoursDori(skCSD)[[1]]
ExtVal<-ColoursDori(membcurr)[[2]] #ColoursDori(skCSD)[[2]]
symbols(segmid[c(1:dimseg)[ToPlot],xA],segmid[c(1:dimseg)[ToPlot],yA], circles=c(diff(range(segstart[,xA]))/Plotwidth* abs(skCSD[,TimeInstant])/ExtVal), add=TRUE, bg=colors()[182], fg=colors()[182], inches=FALSE)
#col.limits<-max(abs(range(skCSD, membrane.curr.line)))#,memb.currents.vonal[,min.indexes])))
szinskala2<-color.scale(c(skCSD[,TimeInstant]),col=col1,zlim=c(-ExtVal*1.01,ExtVal*1.01))
points(elec[,xA],elec[,yA], pch=8)
segments(segstart[c(1:dimseg)[ToPlot],xA], segstart[c(1:dimseg)[ToPlot],yA], segend[c(1:dimseg)[ToPlot],xA], segend[c(1:dimseg)[ToPlot],yA], col= szinskala2,lwd=segwidth[c(1:dimseg)[ToPlot]],lend=1)
#par(mar=c(0.5,0.1,1,0.5)) 
#mtext("mA/um ",side=4,line=1,cex=0.7)
ColMatr<-as.matrix(c(1:256))
col2<-ColoursDori(ColMatr)[[1]]
add.image(max(XrangeEl-10),mean(YrangeEl), t(ColMatr), adj.x = 0.0, adj.y = 0.5,  image.width = 1, image.height = 0.05, col = col2)
#axis(4,at=c(min(Yrange),max(Yrange)),labels=c("sink","source"))
axis(4,at=c(min(YrangeEl),max(YrangeEl)),labels=c("sink","source"))
#points(elec[,xA],elec[,yA])
#if(file.exists(elignorename))  points(elec[ El2Ignore$x,xA],elec[ El2Ignore$x,yA],pch=4,col="RED")
par(oma=c(0,0,2,0))
title(paste("Time:", round(SimTime[TimeInstant],3), "ms") , outer=TRUE)
dev.off()
}
}
PlotMorpho(DataFolder,3)

#system("convert -delay 40 *.png Compare.gif") #taurin muxik

# cleaning up
#file.remove(list.files(pattern=".png"))

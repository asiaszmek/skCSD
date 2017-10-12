source("alprogik/Colors_BlueRed.R")
library(data.table)
#image scale function for plotting color bars for image plots


#DataFolder<-#"/home/csdori/ksCSD_2014/trunk/simulation/Ganglion_d17.5_128Regular/"#BS_d50_el128/"
#DataFolder<-"/home/csdori/ksCSD_2014/trunk/simulation/Y_el4x16_RotS_symm_d50_ver0/"
#DataFolder<-"/media/zoe/BA0ED4600ED416EB/agy/ksCSD_2014/trunk/simulation/Y_el4x16_Rot_symm_d50_ver0/"
DataFolder<-"/home/csdori/ksCSD_2014/trunk/simulation/gang_5x5_400/"
#ffmpeg  -framerate 10 -pattern_type glob -i 'CSD_Morpho*.png' CSD.mp4
#DataFolder<-"/media/zoe/BA0ED4600ED416EB/agy/ksCSD_2014/trunk/simulation/Hex_8x8_inter70_d30/"
plotFileName<-"ComparingSmoothed"
dir.create(paste0(DataFolder,plotFileName))
CellType<-"gangNew"#"Y" # "gangNew" #"Y-rotated" #


if(CellType=="gangNew"){
  xA<-1
  yA<-2
  xAEl<-1 #1
  yAEl<-2
  PlotTitle<-"Ganglion Cell"
  outname1<-"/kernelOut_4/"
  ToPlot<--c(1:10,70)
}

El2Ignore<-numeric()
segstart<-matrix(as.matrix(read.table(paste0(DataFolder,'coordsstart_x_y_z'), colClasses='numeric')),ncol=3)
segend<-matrix(as.matrix(read.table(paste0(DataFolder,'coordsend_x_y_z'), colClasses='numeric')),ncol=3)
segmid<-matrix(as.matrix(read.table(paste0(DataFolder,'coordsmid_x_y_z'), colClasses='numeric')),ncol=3)
elec<-matrix(as.matrix(read.table(paste0(DataFolder,'elcoord_x_y_z'), colClasses='numeric')),ncol=3)
seg.nb<-dim(segmid)[1]




SimTime<-as.matrix(read.table(paste0(DataFolder,'/time')))
TimeMax<-length(SimTime)

########!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#skCSD<-as.matrix(read.table(paste0(DataFolder,outname1,"/skCSD_plain")))#paste0(outname,outname1,"/skCSDall_M",256,"_R",55,"_SNR0_lambda0")))


SNR<-0

coordsEnd<-as.matrix(read.table(paste0(DataFolder,"/coordsend_x_y_z")))


seglength<-as.matrix(read.table(paste0(DataFolder,'/seglength')))
funaramvonal<-function(x) x/seglength

membCurr<-as.matrix(fread(paste0(DataFolder,'/membcurr')))
membcurr0<-apply(membCurr,2,funaramvonal)
somaPot<-as.matrix(read.table(paste0(DataFolder,'/somav.txt')))

smoothGT<-1
if(smoothGT==1){
  source("alprogik/KernSmoothDistance.R")
  memb.currents.smoothed<-KernSmoothDistance(32,DataFolder, DataFolder, outname1,"32")
  
  membcurr0<-memb.currents.smoothed
}




#Plotratio<-Xrange/max(segwidth)
segwidth1<-as.matrix(read.table(paste0(DataFolder,'segdiam_x_y_z'), colClasses='numeric'))
segwidth<-segwidth1/max(segwidth1)*30 # not realistic?
#membcurr<-membcurrO

#membcurrO<-membcurr
membcurr<-membcurr0[ToPlot,] 

Plotwidth<-10
############################################x
############
#########################
png(paste0('/home/csdori/ksCSD_2014/trunk/plots/summary/','GangDiffElMoreL1Predn.png'), width=12, height=14,units='in',res=300)



par(mfrow=c(3,3),oma=c(0,0,1,1),cex=1.13)
TimeInstant<-1981#1981
dimseg<-dim(segstart)[1]
Xrange<-range(segstart[,xA],elec[,xA])*1.05
Yrange<-range(segstart[,yA],elec[,yA])*1.05
XrangeEl<-range(segstart[,xAEl],elec[,xAEl])*1.05
YrangeEl<-range(segstart[,yAEl],elec[,yAEl])*1.05

plot(1,type='n',xlim=Xrange,ylim=Yrange,xlab='x (um)',ylab='y (um)',main="Ground Truth",asp=1)


rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray")
# col1<-ColoursDori(membcurr[,TimeInstant])[[1]]
# ExtVal<-ColoursDori(membcurr[,TimeInstant])[[2]]
col1<-ColoursDori(membcurr[,TimeInstant])[[1]]
ExtVal<-ColoursDori(membcurr[,TimeInstant])[[2]]

#highlight the interesting part of the plot
symbols(segmid[c(1:dimseg)[ToPlot],xA],segmid[c(1:dimseg)[ToPlot],yA], circles=c(diff(range(segstart[,xA]))/Plotwidth* abs(membcurr[,TimeInstant])/ExtVal), add=TRUE, bg=colors()[182], fg=colors()[182], inches=FALSE)
#points(elec[,xA],elec[,yA], pch=8) #diff(Xrange)

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


#####################################################
#####################################################
#TimeInstant<-TimeInstant+1



#multiple CSD
# locations<-c("gang_5x5_50","gang_5x5_100","gang_5x5_200", "gang_7x7_200","gang_9x9_200","gang_9x9_100", "gang_5x5_400") #,"gang_9x9_50"
# InterEl<-c(50,100,200,200,200,100,400)
locations<-c("gang_5x5_50","gang_5x5_100","gang_5x5_200","gang_5x5_400","gang_9x9_200","gang_9x9_100","Berd40_21x21" ) #,"gang_9x9_50"
InterEl<-c(50,100,200,400,200,100,40) #

MinErrors<-numeric()
MinErrorsCV<-numeric()
MinErrorsPred<-numeric()
for(datas in 1:7){
  
  DataFolder<-paste0("/home/csdori/ksCSD_2014/trunk/simulation/",locations[datas],"/")
  outname<-DataFolder
  #Error2Plot<-as.matrix(read.table(paste0(DataFolder,outname1,"/ErrorL1Smoothed_Noise_SNR0")))
  Error2Plot<-as.matrix(read.table( paste(outname,outname1,"/ErrorL1Smoothed_Noise","_SNR",SNR,"_",1,"_",1000,sep="")))
  if(file.exists(paste(DataFolder,outname1,"/ErrorCV_NoiseElec","_SNR",SNR,sep=""))==TRUE){
    CVTable<-as.matrix(read.table(paste(DataFolder,outname1,"/ErrorCV_NoiseElec","_SNR",SNR,sep="")))
    
    CurrMinCV<- min(CVTable)
    MinErrorsCV<-c(MinErrorsCV,CurrMinCV)
    MinIndexCV<-which(CVTable==CurrMinCV,arr.ind=TRUE)
    
    
    
  } else {MinErrorsCV<-c(MinErrorsCV,NA)}
  #32-es símításos
  CurrMin<- min(Error2Plot)
  MinErrors<-c(MinErrors,CurrMin)
  #MinIndexL1<-which(Error2Plot==CurrMin,arr.ind=TRUE)
  
  # Predict parameters from part of the data, use it for the other sets
  El1ErrorRawAllGang<-as.matrix(read.table( paste(outname,outname1,"/ErrorL1Smoothed_Noise","_SNR",SNR,"_",1001,"_",6800,sep="")))
  MinIndex<-as.matrix(read.table(paste(outname,outname1,"/L1Compare_BestParams",sep="")))
  MinErrorsPred<-c(MinErrorsPred, El1ErrorRawAllGang)
  #MinIndex<-which(Error2Plot==CurrMin,arr.ind=TRUE)
  

  MinIndexCV<-which(CVTable==CurrMinCV,arr.ind=TRUE)
  
  
  
  elec<-matrix(as.matrix(read.table(paste0(DataFolder,'elcoord_x_y_z'), colClasses='numeric')),ncol=3)
  
  #MinIndex<-data.matrix(read.table(paste(DataFolder,outname1,"/MinIndex_Noise","_SNR",SNR,sep="")))
  M<-512
  lambda.all<-c("1e-05","1e-04","0.001","0.01","0.1" )
  R.all<-2^(3:7)
  Rread<-R.all[MinIndex[1]]
  Lread<-lambda.all[MinIndex[2]]
  C.calc<-as.matrix(fread(paste0(DataFolder,outname1,"/ksCSD_Matr512_R",Rread,"lambda",Lread,"SNR0")))
  #} else {
  #  C.calc<-as.matrix(read.table(paste0(DataFolder,outname1,"/skCSDall_M512_R30lambda1e-04")))#30lambda1e-04")))#80lambda0.001")))
  #}
  SameplaceAll<-readLines(paste0(DataFolder,outname1,"/sameplace.txt"))
  skCSD.all<-array(0,c(seg.nb, dim(C.calc)[2]))
  for(i in 1: seg.nb){
    sameplace<-numeric() 
    sameplace<-as.integer(unlist(strsplit(SameplaceAll," ")[i]))
    #Plot the current related ti the different segments
    #  matplot(cbind(t(C.calc[c(sameplace),])),t='l', xlab="Time steps", ylab="ksCSD", main="Values of current at the same location")
    if(length(sameplace)==1) skCSD.all[i,]<-C.calc[sameplace,]
    if(length(sameplace)>1) skCSD.all[i,]<-as.matrix(colSums(C.calc[c(sameplace),]),nrow=1)
  }
  skCSD<-skCSD.all[ToPlot,]
  
  
  
  plot(1,type='n',xlim=Xrange,ylim=Yrange,xlab='x (um)',ylab='y (um)',main=paste("IED", InterEl[datas], "um"),asp=1)
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray")
  # col1<-ColoursDori(skCSD[,TimeInstant])[[1]]
  # ExtVal<-ColoursDori(skCSD[,TimeInstant])[[2]]
  col1<-ColoursDori(skCSD[,TimeInstant])[[1]]
  ExtVal<-ColoursDori(skCSD[,TimeInstant])[[2]]#ColoursDori(membcurr[,TimeInstant])[[2]] #ColoursDori(skCSD)[[2]]
  symbols(segmid[c(1:dimseg)[ToPlot],xA],segmid[c(1:dimseg)[ToPlot],yA], circles=c(diff(range(segstart[,xA]))/Plotwidth* abs(skCSD[,TimeInstant])/ExtVal), add=TRUE, bg=colors()[182], fg=colors()[182], inches=FALSE)
  
  
  #col.limits<-max(abs(range(skCSD, membrane.curr.line)))#,memb.currents.vonal[,min.indexes])))
  
  szinskala2<-color.scale(c(skCSD[,TimeInstant]),col=col1,zlim=c(-ExtVal*1.01,ExtVal*1.01))
  points(elec[,xA],elec[,yA], pch=8)
  segments(segstart[c(1:dimseg)[ToPlot],xA], segstart[c(1:dimseg)[ToPlot],yA], segend[c(1:dimseg)[ToPlot],xA], segend[c(1:dimseg)[ToPlot],yA], col= szinskala2,lwd=segwidth[c(1:dimseg)[ToPlot]],lend=1)
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
  #par(oma=c(0,0,2,0))
  title(paste("Time:", round(SimTime[TimeInstant],3), "ms") , outer=TRUE)
}



plot(1:length(MinErrors),MinErrors,pch=20, ylab="L1 Error",xaxt="n",xlab="Setup",main="Reconstruction Error",ylim=c(0.2,0.82))
axis(1,at=c(1:length(MinErrors)),c("B","C","D","E","F","G","H"))
points(MinErrorsPred,col="GREEN",pch=8)
MinErrorsCV<-MinErrorsCV[-7]
TransCVError<-(MinErrorsCV-min(MinErrorsCV))/diff(range(MinErrorsCV))*diff(range(MinErrors,MinErrorsPred))+min(MinErrors,MinErrorsPred)
#TransCVError<-(MinErrorsCV-min(MinErrorsCV))/diff(range(MinErrorsCV))*diff(c(0.2,0.8))+0.2

points(TransCVError,pch=4,col="RED")
axis(4,at=c(min(MinErrors,MinErrorsPred),max(MinErrors,MinErrorsPred)),labels=c(round(min(MinErrorsCV),1),round(max(MinErrorsCV),1)))
mtext("CV Error",4, line=2,col="RED")
legend("bottomleft",c("L1 - T","L1 - V", "CV"), col=c("BLACK", "GREEN", "RED"), pch=c(20,8,4))

dev.off()












##########################################################x
#Mindenfele Ganglion
##########################################################
TimeInstant<-1981
dimseg<-dim(segstart)[1]
Xrange<-range(segstart[,xA],elec[,xA])*1.05
Yrange<-range(segstart[,yA],elec[,yA])*1.05
XrangeEl<-range(segstart[,xAEl],elec[,xAEl])*1.05
YrangeEl<-range(segstart[,yAEl],elec[,yAEl])*1.05

#multiple CSD
locations<-c("gang_5x5_25","gang_5x5_50","gang_5x5_100","gang_5x5_200", "gang_7x7_200","gang_9x9_200","gang_9x9_100", "gang_5x5_400", "gang_5x5_200", "gang_10x10_50","gang_20x20_25","gang_20x20_50","gang_10x10_100","gang_9x9_50","gang_17x17_50","Berd40_21x21") #,"gang_9x9_50"


#multiple CSD
InterEl<-c(25,50,100,200,200,200,100,400,200,50,25,50,100,50,50,40)
MinErrors<-numeric()
for(datas in 17:length(locations)){
  png(paste0('/home/csdori/ksCSD_2014/trunk/plots/summary/',locations[datas],'.png'), width=5, height=5,units='in',res=300)
  DataFolder<-paste0("/home/csdori/ksCSD_2014/trunk/simulation/",locations[datas],"/")
  if(file.exists(paste0(DataFolder,outname1,"/ErrorL1Smoothed_Noise_SNR0"))==FALSE) next
  Error2Plot<-as.matrix(read.table(paste0(DataFolder,outname1,"/ErrorL1Smoothed_Noise_SNR0")))
  #try(CVTable<-as.matrix(read.table(paste(outname,outname1,"/ErrorCV_Noise","_SNR",SNR,sep=""))))
  #32-es símításos
  CurrMin<- min(Error2Plot)
  MinErrors<-c(MinErrors,CurrMin)
  MinIndex<-which(Error2Plot==CurrMin,arr.ind=TRUE)
  
  elec<-matrix(as.matrix(read.table(paste0(DataFolder,'elcoord_x_y_z'), colClasses='numeric')),ncol=3)
  
  #MinIndex<-data.matrix(read.table(paste(DataFolder,outname1,"/MinIndex_Noise","_SNR",SNR,sep="")))
  M<-512
  lambda.all<-c("1e-05","1e-04","0.001","0.01","0.1" )
  R.all<-2^(3:7)
  Rread<-R.all[MinIndex[1]]
  Lread<-lambda.all[MinIndex[2]]
  C.calc<-as.matrix(fread(paste0(DataFolder,outname1,"/ksCSD_Matr512_R",Rread,"lambda",Lread,"SNR0")))
  #} else {
  #  C.calc<-as.matrix(read.table(paste0(DataFolder,outname1,"/skCSDall_M512_R30lambda1e-04")))#30lambda1e-04")))#80lambda0.001")))
  #}
  SameplaceAll<-readLines(paste0(DataFolder,outname1,"/sameplace.txt"))
  skCSD.all<-array(0,c(seg.nb, dim(C.calc)[2]))
  for(i in 1: seg.nb){
    sameplace<-numeric() 
    sameplace<-as.integer(unlist(strsplit(SameplaceAll," ")[i]))
    #Plot the current related ti the different segments
    #  matplot(cbind(t(C.calc[c(sameplace),])),t='l', xlab="Time steps", ylab="ksCSD", main="Values of current at the same location")
    if(length(sameplace)==1) skCSD.all[i,]<-C.calc[sameplace,]
    if(length(sameplace)>1) skCSD.all[i,]<-as.matrix(colSums(C.calc[c(sameplace),]),nrow=1)
  }
  skCSD<-skCSD.all[ToPlot,]
  
  
  
  plot(1,type='n',xlim=Xrange,ylim=Yrange,xlab='x (um)',ylab='y (um)',main=paste("IED", InterEl[datas], "um"),asp=1)
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray")
  # col1<-ColoursDori(skCSD[,TimeInstant])[[1]]
  # ExtVal<-ColoursDori(skCSD[,TimeInstant])[[2]]
  col1<-ColoursDori(skCSD[,TimeInstant])[[1]]
  ExtVal<-ColoursDori(skCSD[,TimeInstant])[[2]]#ColoursDori(membcurr[,TimeInstant])[[2]] #ColoursDori(skCSD)[[2]]
  symbols(segmid[c(1:dimseg)[ToPlot],xA],segmid[c(1:dimseg)[ToPlot],yA], circles=c(diff(range(segstart[,xA]))/Plotwidth* abs(skCSD[,TimeInstant])/ExtVal), add=TRUE, bg=colors()[182], fg=colors()[182], inches=FALSE)
  
  
  #col.limits<-max(abs(range(skCSD, membrane.curr.line)))#,memb.currents.vonal[,min.indexes])))
  
  szinskala2<-color.scale(c(skCSD[,TimeInstant]),col=col1,zlim=c(-ExtVal*1.01,ExtVal*1.01))
  points(elec[,xA],elec[,yA], pch=8)
  segments(segstart[c(1:dimseg)[ToPlot],xA], segstart[c(1:dimseg)[ToPlot],yA], segend[c(1:dimseg)[ToPlot],xA], segend[c(1:dimseg)[ToPlot],yA], col= szinskala2,lwd=segwidth[c(1:dimseg)[ToPlot]],lend=1)
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
  
  
  dev.off()
  
}










